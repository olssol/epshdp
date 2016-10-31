

##
##  check hdp model
##
##

rm(list=ls());
options(error = recover, warn=2);
source("ps_toolkit.R");

make.global(set.simu.par(11));


## check rep 1
rct.data <- read.table(get.f.name(simu.ext=20,
                                  cur.rep=1,
                                  lsts=c("datarct")), header=T);

rst.obs <- get.obs.effect(rct.data, study="Study", group="Z", y="Y");
hdp.ps  <- get.all.ps(rct.data, study="Study", group="Z",
                      ps.cov=paste("X", 1:4, sep=""),
                      delta=0,
                      take.logit=TRUE);

chk.ps <- subset(hdp.ps, Study == 1);

rst.ps  <- get.ps.effect(rct.data,
                         study="Study", group="Z",
                         ps.cov=paste("X", 1:4, sep=""),
                         y="Y",
                         n.breaks=5,
                         delta=0,
                         n.boots=5);

grp0 <- read.table("simu_20/rep_1/simu_grp0/pred_1_0");
grp1 <- read.table("simu_20/rep_1/simu_grp1/pred_1_0");

g0.pred <- sqldf("select V3, avg(V5) as Pred
                  from grp0
                  group by V3");

##nls
ny <- nls(Y ~ score, data=chk.ps);

plot(chk.ps$score, chk.ps$Y);
points(chk.ps$score, g0.pred$Pred, col="red");

chk <- subset(grp0, 101 == V3);
plot(density(chk$V5));


mean(grp1$V5) - mean(grp0$V5)



## checking
cur.par <- set.simu.par.19(sizes = c(10000,10000,10000))$SIMU.PAR;
for (i in 1:length(cur.par)) {
    cur.d <- GenDataMatrix(i, cur.par);
    print(table(cur.d[,"Z"]));
    m.t0  <- mean(cur.d[which(0 == cur.d[,"Z"]),"Y"]);
    m.t1  <- mean(cur.d[which(1 == cur.d[,"Z"]),"Y"]);
    print(m.t1 - m.t0);
}

##---1. observed results


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

