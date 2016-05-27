##Time-stamp: <2016-03-07 10:01:46 cwang68>

rm(list=ls());
options(error = recover, warn=2);
source("ps_toolkit.R");

ARGS.INX <- as.numeric(commandArgs(trailingOnly=TRUE));
if (0 == length(ARGS.INX)) {
    ARGS.INX <- 0;
    print(paste("No argument entered"));
}

SIMU.ONLY <- (0 == ARGS.INX);
##quit if argument exceeds number of nodes
if (N.NODES < ARGS.INX) {
    quit(save="no");
}

##set random seed;
set.seed(10*ARGS.INX + 1000);

##----------------------------------------------------------
##                 SIMULATION
##----------------------------------------------------------
##simulate data with two copies for binary and continuous analysis
if (SIMU.ONLY) {
    for (s in 1:length(SIMU.EXT.ALL)) {
        SIMU.EXT <- SIMU.EXT.ALL[s];

        ##simulation parameters
        make.global(set.simu.par(SIMU.EXT));
        save(SIMU.PAR.LST, SIMU.NOTE, SIMU.PS.COV,
             file=get.f.name(simu.ext=SIMU.EXT, lsts=c("simu_par")));

        ##simu data
        for (i in 1:N.REPS) {
            cur.f <- get.f.name(simu.ext=SIMU.EXT,
                                cur.rep=i, lsts=c("datarct"));

            if (!file.exists(cur.f) | IGNORE.EXISTING) {
                print(paste("Simu:", cur.f, sep=""));

                simu.data <- NULL;
                for (j in 1:length(SIMU.PAR.LST)) {
                    pt.study.j  <- GenDataMatrix(j, SIMU.PAR.LST);
                    simu.data   <- rbind(simu.data, pt.study.j);
                }

                write.table(simu.data, file=cur.f,
                            sep=" ", quote=F, row.names=F, col.names=T)
            }
        }
    }
} else {
    ##replications
    rep.start <- (ARGS.INX-1)*N.EACH+1;
    rep.end   <- ARGS.INX*N.EACH;

    for (rr in rep.start:rep.end) {
        s.inx    <- ceiling(rr / N.REPS);
        SIMU.EXT <- SIMU.EXT.ALL[s.inx];
        r        <- rr - (s.inx-1)*N.REPS;

        ##----------load simu parameters------------
        load(file=get.f.name(SIMU.EXT, lsts=c("simu_par")));

        ##print info
        print(paste( "---CUR SIM.EXT:--", SIMU.EXT, "---CUR REP:",r, "---", sep=""));

        ##------------read data----------------------
        rct.data <- read.table(get.f.name(simu.ext=SIMU.EXT, cur.rep=r, lsts=c("datarct")), header=T);
        STUDIES  <- sort(unique(rct.data$Study));
        ##HDP.COV  <- paste("score", 1:length(STUDIES), sep="");
        HDP.COV  <- "score";

        ##------------analysis-----------------------
        f.rst <- get.f.name(simu.ext=SIMU.EXT, cur.rep=r, lsts=c("result"));

        if (file.exists(f.rst) & IGNORE.EXISTING)
  	    next;

        ##---1. observed results
        rst.obs <- get.obs.effect(rct.data, study="Study", group="Z", y="Y");

        ##set work dir
        old.wd <- getwd();
        setwd(get.f.name(SIMU.EXT, cur.rep=r));

        rst.ps   <- list(NULL);
        rst.hdp  <- list(NULL);
        for (d in 1:length(DELTA.ALL)) {

            ##sensitivity parameter
            delta  <- DELTA.ALL[d];

            ##---2. propensity score analysis
            rst.ps[[d]]  <- get.ps.effect(rct.data,
                                          study="Study", group="Z",
                                          ps.cov=SIMU.PS.COV,
                                          y="Y",
                                          n.breaks=PS.BREAKS,
                                          delta=delta,
                                          n.boots=5);

            ##---3. hdp model analysis
            ##all ps scores
            hdp.ps   <- get.all.ps(rct.data, study="Study", group="Z",
                                   ps.cov=SIMU.PS.COV,
                                   delta=delta,
                                   take.logit=TRUE);
            pred.pts <- as.matrix(hdp.ps[, c("Study", HDP.COV), drop=FALSE]);

            ##options
            ops <- list(n.iter=N.ITER, n.discard=N.DISCARD,
                        n.batch=20,
                        mcmc.eps=0, eps=1,
                        m.prior=1,
                        B.prior=1, S.prior=1, alpha.prior=1,
                        alpha=50, q=10, cc=10);


            ##fit and predict
            for (grp in TRT.NUMBER) {
                ##go to grp folder
                cur.owd <- getwd();
                setwd(set.simu.dir(paste("grp", grp, sep="")));

                ##assign data
                cur.d <- subset(hdp.ps, grp == Z);
                ##sort by study no first
                cur.d <- cur.d[order(cur.d$Study),];

                ##fitting
                cur.z     <- as.matrix(cur.d[, c("Y", HDP.COV)]);
                cur.study <- cur.d$Study;
                do.call(hdpmn, c(list(Z=cur.z, study=cur.study), px=ncol(cur.z)-1, ops));

                ##predicting
                for (sy in STUDIES) {
                    cur.pts  <- pred.pts[which(pred.pts[,'Study'] == sy),];
                    cur.pred <- R.hdpmnPredict(j=sy,
                                               r=PRED.R,
                                               nsim=0,
                                               idx.x=1+1:length(HDP.COV),
                                               X=cur.pts[, HDP.COV,drop=FALSE]);

                    ## kk       <- 0;
                    ## while (kk < nrow(cur.pts)) {
                    ##     k.inxs   <- (kk+1):min(nrow(cur.pts), kk+PRED.BATCH);
                    ##     cur.pred <- R.hdpmnPredict(j=sy,
                    ##                                r=PRED.R,
                    ##                                nsim=0,
                    ##                                idx.x=1+1:length(HDP.COV),
                    ##                                X=cur.pts[k.inxs, HDP.COV,drop=FALSE]);

                    ##     ##adjust pt id
                    ##     cur.pred[,3] <- cur.pred[,3] + kk;
                    ##     pred.rst     <- rbind(pred.rst, cbind(grp, sy, cur.pred[,3:5]));
                    ##     kk           <- kk + PRED.BATCH;
                    ##     gc(verbose=TRUE);
                    ## }
                }

                ##go back
                setwd(cur.owd);
            }

            ##post analysis sum hdp
            for (s in STUDIES) {
                colnames(pred.rst) <- c("Grp", "Study", "ID", "Iter", "Pred");
                pred.rst <- read.table(pred.rst);
                cur.sy <- subset(pred.rst, Study == s);

                if (0 == nrow(cur.sy))
                    next;

                all.iter <- unique(cur.sy$Iter);
                sy.diff  <- NULL;
                for (k in all.iter) {
                    py <- NULL;
                    for (g in TRT.NUMBER) {
                        cur.p <- subset(cur.sy, Grp == g & Iter == k);
                        py    <- cbind(py, cur.p[order(cur.p$ID), 'Pred']);
                    }
                    sy.diff <- rbind(sy.diff, mean(py[,2] - py[,1]));
                }

                rst.hdp[[d]] <- rbind(rst.hdp[[d]],
                                      c(mean(sy.diff),
                                        quantile(sy.diff, c(0.025, 0.975))));
            }
        }

        ##return to previous dir
        setwd(old.wd);

        ##-----------save results------------------
        save(rst.obs, rst.ps, rst.hdp,
             DELTA.ALL, SIMU.PAR.LST, SIMU.NOTE, HDP.COV,
             file=f.rst);
    }
}

