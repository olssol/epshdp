##
##  check hdp model
##
##

rm(list=ls());
options(error = recover, warn=2);
source("ps_toolkit.R");

ARGS.INX <- as.numeric(commandArgs(trailingOnly=TRUE));
if (0 == length(ARGS.INX)) {
    ARGS.INX <- 1;
    print(paste("No argument entered"));
}

SIMU.ONLY <- (0 == ARGS.INX);

##quit if argument exceeds number of nodes
if (N.NODES < ARGS.INX) {
    quit(save="no");
}

##set random seed;
set.seed(10*ARGS.INX + 1000);

SIMU.EXT <- 20;
r        <- 2;
delta    <- 0;

##----------load simu parameters------------
load(file=get.f.name(SIMU.EXT, lsts=c("simu_par")));

simu.data <- NULL;
for (j in 1:length(SIMU.PAR.LST)) {
    pt.study.j  <- GenDataMatrix(j, SIMU.PAR.LST);
    simu.data   <- rbind(simu.data, pt.study.j);
}


##----------load data-----------------------
rct.data <- read.table(get.f.name(simu.ext=SIMU.EXT,
                                  cur.rep=r,
                                  lsts=c("datarct")), header=T);
rct.data <- simu.data;

STUDIES  <- sort(unique(rct.data$Study));
##HDP.COV  <- paste("score", 1:length(STUDIES), sep="");
HDP.COV  <- c("score");

hdp.ps   <- get.all.ps(rct.data, study="Study", group="Z",
                       ps.cov=SIMU.PS.COV,
                       delta=delta,
                       take.logit=TRUE);
pred.pts <- as.matrix(hdp.ps[, c("Study", HDP.COV), drop=FALSE]);


##---1. observed results
rst.obs <- get.obs.effect(rct.data, study="Study", group="Z", y="Y");

rst.ps  <- get.ps.effect(rct.data,
                         study="Study", group="Z",
                         ps.cov=SIMU.PS.COV,
                         y="Y",
                         n.breaks=PS.BREAKS,
                         delta=delta,
                         n.boots=5);

##options
ops <- list(n.iter=10000, n.discard=4000,
            n.batch=20,
            mcmc.eps=0, eps=1,
            m.prior=1,
            B.prior=1, S.prior=1, alpha.prior=1,
            alpha=50, q=10, cc=10);

if (1) {
    ##fitting
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

        ##go back
        setwd(cur.owd);
    }
}

##predict
p.batch  <- 40;
pred.rst <- NULL;
for (grp in TRT.NUMBER) {
    ##go to grp folder
    cur.owd <- getwd();
    setwd(set.simu.dir(paste("grp", grp, sep="")));

    ##predicting
    tic();
    for (sy in 1:2) {
        cur.pts  <- pred.pts[which(pred.pts[,'Study'] == sy),];
        kk <- 0;
        while (kk < nrow(cur.pts)) {
            k.inxs   <- (kk+1):min(nrow(cur.pts), kk+p.batch);
            cur.pred <- R.hdpmnPredict(j=sy,
                                       r=PRED.R,
                                       nsim=0,
                                       idx.x=1+1:length(HDP.COV),
                                       X=cur.pts[k.inxs, HDP.COV,drop=FALSE]);
            cur.pred[,3] <- cur.pred[,3] + kk;
            pred.rst     <- rbind(pred.rst, cbind(grp, sy, cur.pred[,3:5]));
            kk           <- kk + p.batch;
            gc(verbose=TRUE);
        }
    }
    toc();
    ##go back
    setwd(cur.owd);
}

colnames(pred.rst) <- c("Grp", "Study", "ID", "Iter", "Pred");
pred.rst <- data.frame(pred.rst);
##save(pred.rst, file="tmp_chk.Rdata");

##sum
rst.hdp <- NULL;
for (s in STUDIES) {
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

    rst.hdp <- rbind(rst.hdp,
                     c(mean(sy.diff),
                       quantile(sy.diff, c(0.025, 0.975))));
}

print(rst.obs);
print(rst.ps);
print(rst.hdp);

##plot
plot(NULL, xlim=c(-5,5), ylim=c(-40, 20));
for (s in 1:2) {
    for (g in 0:1) {
        cur.d <- subset(hdp.ps, Study == s & Z == g);
        points(cur.d$score, cur.d$Y, col=c("black","red", "blue", "yellow")[(s-1)*2+g+1]);
    }
}

par(mfrow=c(1,2));
for (s in 1:2) {
    cur.d <- subset(hdp.ps, Study == s & Z == 0);
    plot(density(cur.d$Y), ylim=c(0,0.4), xlim=c(-60,40), main=paste("study", s));
    cur.d <- subset(hdp.ps, Study == s & Z == 1);
    lines(density(cur.d$Y), col="red");

    cur.p <- subset(pred.rst, Study == s & Grp == 0);
    lines(density(cur.p$Pred), lty=2);

    cur.p <- subset(pred.rst, Study == s & Grp == 1);
    lines(density(cur.p$Pred), lty=2, col="red");
}

