rm(list=ls());
options(error = recover, warn=2);
source("ps_toolkit.R");

##----------------------------------------------------------------------------
##                FUNCTIONS
##----------------------------------------------------------------------------

get.zj.data <- function(fname = "zhenjing.csv") {
    all.data <- read.table(file=fname, header=T, sep=",");

    ##percentage wrong in raw data
    inx.spo2                <- which(all.data$SPO2 <=1);
    all.data$SPO2[inx.spo2] <- 100*all.data$SPO2[inx.spo2];
    all.data$spo2           <- as.integer(all.data$SPO2 >= 95);

    all.data$asa            <- as.integer(1 == all.data$ASA.grade);
    all.data$allergy        <- as.integer(0 != all.data$allergy);
    all.data$surgery        <- as.integer(0 != all.data$surgery);
    all.data$medicine       <- as.integer(0 != all.data$medicine);
    all.data$hypertension   <- as.integer(0 != all.data$history.of.hypertension);
    all.data$diabetes       <- as.integer(0 != all.data$diabetes);
    all.data$heart          <- as.integer(0 != all.data$heart.disease);
    all.data$ae             <- all.data$adverse.events;
    all.data$comfort        <- all.data$degree.of.comfort;

    ##exclud the following columns
    ## reasons: all non-ane sujects have reason 0
    ## examination: only 5 subject have 2
    all.data <- subset(all.data, 1 == examination & 0 == SAS
                       & 0 == COPD & 0 == asthma & 0 == cerebral.infarction);

    ##return
    rst <- all.data[, c("hospital",   "group",   "sex",     "age", "asa",
                        "diagnosis",  "BMI",  "spo2", "allergy", "surgery",
                        "medicine",   "hypertension", "diabetes", "heart",
                        "ae", "comfort")]

    inx.comp <- complete.cases(rst);
    rst[inx.comp,];

}


##----------------------------------------------------------------------------
##                CONSTANTS
##----------------------------------------------------------------------------

##analysis covariates for sedation data
PS.COV      <- c("sex",      "age",     "asa",     "diagnosis",
                 "BMI",      "spo2",    "allergy", "surgery",
                 "medicine", "hypertension",
                 "diabetes", "heart");


ARGS.INX <- as.numeric(commandArgs(trailingOnly=TRUE));
if (0 == length(ARGS.INX)) {
    ARGS.INX <- 1;
    print(paste("No argument entered"));
}
set.seed(as.numeric(Sys.time()) + ARGS.INX * 1e4);

N.NODES  <- 5000; ##no of nodes
SUFFIX   <- "_D";
PRED.R   <- 0;
DO.FIT   <- FALSE;
DO.PRED  <- TRUE;
DO.POST  <- FALSE;

##------------------------------------------------
##                  get data
##------------------------------------------------
##zhenjing study
all.data <- get.zj.data();

##tbl.baseline(all.data);
STUDIES <- sort(unique(all.data$hospital));


##------------------------------------------------
##            TYPICAL AE Analysis
##------------------------------------------------
##observed ae
obs.ae <- get.obs.effect(all.data, study="hospital", group="group", y="ae");

##ps adjusted
ps.ae  <- get.ps.effect(all.data,
                        study="hospital", group="group", y="ae",
                        n.breaks=PS.BREAKS);
for (i in 1:500) {
    print(i);
    cur.dta <- sample(1:nrow(all.data), nrow(all.data), replace=TRUE);
    cur.ae  <- get.ps.adj(all.data[cur.dta,], n.breaks=N.BREAKS);
    ps.ae   <- cbind(ps.ae, cur.ae);
}

rst.ps.ae <- rbind(ps.ae[,1],
                   apply(ps.ae, 1, quantile, c(0.025, 0.975), na.rm=TRUE));
rst.ps.ae <- t(rst.ps.ae);

##------------------------------------------------
##            HDP Analysis
##------------------------------------------------
##hdp analysis
hdp.ps <- get.all.ps(all.data);

if (DO.FIT) {
    ##options
    ops <- list(n.iter=200,
                n.discard=20,
                n.batch=10,
                mcmc.eps=1,
                eps=0.5,
                m.prior=1,
                B.prior=1,
                S.prior=1,
                alpha.prior=1,
                alpha=50,
                q=10,
                cc=10);

    ##directory
    old.wd <- getwd();
    setwd(set.simu.dir(paste("tmp", SUFFIX, sep="")));
    file.copy("../zhenjing.R", paste("script_zhenjing", SUFFIX, ".R",sep=""));

    for (grp in 0:1) {
        ##go to grp folder
        cur.owd <- getwd();
        setwd(set.simu.dir(paste("grp", grp, sep="")));

        ##assign data
        cur.d <- subset(hdp.ps, grp == group);
        ##sort by study no first
        cur.d <- cur.d[order(cur.d$hospital),];
        ##fitting
        cur.z     <- as.matrix(cur.d[, c("ae", paste("score", 1:length(STUDIES), sep=""))]);
        cur.study <- cur.d$hospital;
        do.call(hdpmn, c(list(Z=cur.z,study=cur.study), px=ncol(cur.z)-1, ops));

        ##go back
        setwd(cur.owd);
    }

    ##set back
    setwd(old.wd);
}

if (DO.PRED) {
    pts.all  <- as.matrix(hdp.ps[, paste("score", 1:length(STUDIES), sep="")]);

    ##parallel
    ntot  <- nrow(pts.all);
    neach <- floor(ntot/N.NODES);
    rep.1 <- (ARGS.INX-1)*ntot+1;
    rep.2 <- min(ntot, rep.1 + neach);

    pred.rst <- NULL;
    for (grp in 0:1) {
        for (s in 1:length(STUDIES)) {
            ## cur.pred <- hdpmnPredict(j=STUDIES[s],
            ##                            r=PRED.R,
            ##                            nsim=10,
            ##                            px=length(STUDIES),
            ##                            X=pts.all[rep.1:rep.2,,drop=FALSE],
            ##                            header=T,
            ##                            work.dir=paste(getwd(), "/tmp", SUFFIX, "/grp", grp, sep=""));
            cur.pred <- R.hdpmnPredict(j=STUDIES[s],
                                       r=PRED.R,
                                       nsim=10,
                                       idx.x=1+1:length(STUDIES),
                                       X=pts.all[rep.1:rep.2,,drop=FALSE],
                                       work.dir=paste(getwd(), "/tmp", SUFFIX, "/grp", grp, sep=""));
            pred.rst <- rbind(pred.rst, cbind(grp, s, cur.pred));
        }
    }

    ##save
    old.wd <- getwd();
    setwd(set.simu.dir(paste("tmp", SUFFIX, sep="")));
    save(pred.rst, file=paste("predrst_", ARGS.INX, ".Rdata", sep=""));
    setwd(old.wd);
}

if (DO.POST) {
    post.f <- paste("tmp", SUFFIX, "/predrst.Rdata", sep="");
    if (!file.exists(post.f)) {
        pred.rst <- post.combine(n=500,
                                 prefix = paste("tmp", SUFFIX, "/predrst/predrst_", sep=""));
        save(pred.rst, file=post.f);
    } else {
        load(post.f);
    }

    dtap               <- data.frame(pred.rst);
    dtap[,'cur.pred']  <- as.numeric(dtap[,'cur.pred'] > 0);
    colnames(dtap)     <- c("inx", "grp", "study", "j", "r", "pt", "v", "pred");


    dtap.avg <- sqldf('select inx, grp, study, pt, avg(pred) as pred
                       from dtap
                       group by inx, grp, study, pt');

    rst.hdp <- NULL;
    for (stu in 1:length(STUDIES)) {
        cur.s.0 <- subset(dtap.avg, grp == 0 & stu == study);
        cur.s.1 <- subset(dtap.avg, grp == 1 & stu == study);
        diff    <- cur.s.0$pred - cur.s.1$pred;
        rst.hdp <- rbind(rst.hdp,
                         c(mean(diff),
                           quantile(diff, c(0.025, 0.975))));
    }
    rst.hdp
}

