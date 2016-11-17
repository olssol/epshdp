##Time-stamp: <2016-06-05 23:47:34 cwang68>

rm(list=ls());
options(error = recover, warn=0);
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

##----------------------------------------------------------
##                 SIMULATION
##----------------------------------------------------------
##simulate data with two copies for binary and continuous analysis
if (SIMU.ONLY) {
    for (s in 1:length(SIMU.EXT.ALL)) {
        SIMU.EXT <- SIMU.EXT.ALL[s];

        ##simulation parameters
        make.global(set.simu.par(SIMU.EXT));
        save(SIMU.PAR.LST, SIMU.NOTE,
             SIMU.PS.COV, SIMU.SCE,
             SIMU.HDP.COV, SIMU.OPS,
             file=get.f.name(simu.ext=SIMU.EXT, lsts=c("simu_par")));

        ##simu data
        set.seed(SIMU.SCE*1000);
        for (i in 1:N.REPS) {
            cur.f <- get.f.name(simu.ext=SIMU.EXT,
                                cur.rep=i, lsts=c("datarct"));

            if (!file.exists(cur.f) | !KEEP.EXISTING) {
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

    ##set random seed;
    set.seed(10*ARGS.INX + 1000);

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
        rct.data   <- read.table(get.f.name(simu.ext=SIMU.EXT, cur.rep=r, lsts=c("datarct")),
                                 header=T);
        STUDIES    <- sort(unique(rct.data$Study));
        rct.data$Y <- rct.data$Y;

        ##------------analysis-----------------------
        f.rst <- get.f.name(simu.ext=SIMU.EXT, cur.rep=r, lsts=c("result"));
        if (file.exists(f.rst) & KEEP.EXISTING)
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
                                          study="Study",
                                          group="Z",
                                          random="Randomized",
                                          ps.cov=SIMU.PS.COV,
                                          y="Y",
                                          n.class=PS.CLASS,
                                          delta=delta);

            ##---3. hdp model analysis

            ##all ps scores for observational studies
            hdp.ps <- get.all.ps(rct.data,
                                 study="Study",
                                 group="Z",
                                 random="Randomized",
                                 ps.cov=SIMU.PS.COV,
                                 delta=delta,
                                 take.logit=TRUE);

            HDP.COV  <- SIMU.HDP.COV[which(SIMU.HDP.COV %in% colnames(hdp.ps))];
            pred.pts <- as.matrix(hdp.ps[, c("Study", HDP.COV), drop=FALSE]);

            ##fit and predict
            for (grp in TRT.NUMBER) {
                ##go to grp folder
                cur.owd <- getwd();
                setwd(set.simu.dir(paste("grp", grp, sep="")));

                ##assign data
                cur.d <- subset(hdp.ps, grp == Z);
                cur.d[,"Y"] <- cur.d[,"Y"] * CONST.Y;
                ##sort by study no first
                cur.d <- cur.d[order(cur.d$Study),];

                ##fitting
                cur.z     <- as.matrix(cur.d[, c("Y", HDP.COV)]);
                cur.study <- cur.d$Study;
                do.call(hdpmn, c(list(Z=cur.z, study=cur.study), px=ncol(cur.z)-1, SIMU.OPS));

                ##predicting
                for (sy in STUDIES) {
                    cur.pts  <- pred.pts[which(pred.pts[,'Study'] == sy),];
                    cur.pred <- R.hdpmnPredict(j=sy,
                                               r=PRED.R,
                                               nsim=0,
                                               idx.x=1+1:length(HDP.COV),
                                               X=cur.pts[, HDP.COV, drop=FALSE]);
                }

                ##go back
                setwd(cur.owd);
            }

            ##post analysis sum hdp
            for (sy in STUDIES) {
                pred.g <- NULL;
                for (g in 1:length(TRT.NUMBER)) {
                    simu.dir        <- set.simu.dir(paste("grp", TRT.NUMBER[g], sep=""));
                    cur.g           <- read.table(paste(simu.dir, "/pred_", sy, "_", PRED.R, sep=""));
                    colnames(cur.g) <- c("ID", "Iter", "Pred");
                    pred.g[[g]]     <- cur.g;
                }

                all.iter <- unique(pred.g[[1]][,"Iter"]);
                sy.diff  <- NULL;
                for (k in all.iter) {
                    py <- NULL;
                    for (g in 1:length(TRT.NUMBER)) {
                        cur.p <- subset(pred.g[[g]], Iter == k);
                        py    <- cbind(py, cur.p[order(cur.p$ID), 'Pred']);
                    }

                    ##transfer back
                    py <- py / CONST.Y;

                    sy.diff <- c(sy.diff, mean(py[,2] - py[,1]));
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
             DELTA.ALL, SIMU.PAR.LST, SIMU.NOTE, HDP.COV, HDP.OPS,
             file=f.rst);
    }
}

