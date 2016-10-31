##Time-stamp: <2016-06-07 18:18:47 cwang68>

library(cgtool);
library(hdpmn);
library(sqldf);
library(mvtnorm);

##----------------------------------------------------------------------------
##                CONSTANTS
##----------------------------------------------------------------------------

##labels of treatment arms
TRT.NUMBER  <- c(0,1);

##analysis covariates for sedation data
PS.COV      <- c("sex",      "age",     "asa",     "diagnosis",
                 "BMI",      "spo2",    "allergy", "surgery",
                 "medicine", "hypertension",
                 "diabetes", "heart");

##simulation scenarios
SIMU.EXT.ALL    <- c(71:73);

##sensitivity
DELTA.ALL       <- c(0);

##ignore existing
IGNORE.EXISTING <- TRUE;

##propensity score breaks
PS.CLASS        <- 5;

##number bootstraps
N.BOOTS         <- 100;

##MCMC
HDP.OPS <- list(n.iter=6000,
                n.discard=2000,
                n.batch=20,
                mcmc.eps=1,
                eps=0.5,
                m.prior=1,
                B.prior=1, S.prior=1, alpha.prior=1,
                alpha=20,
                q=10,
                cc=10);

## Simulation set up
N.REPS          <- 500; ##number of total replications
N.NODES         <- 500; ##number of computation nodes
N.EACH          <- ceiling(N.REPS*length(SIMU.EXT.ALL) / N.NODES);

##HDP prediction
PRED.R          <- 0;

##ratio for Y to make Y comparable to scores
CONST.Y         <- 0.02;

##----------------------------------------------------------------------------
##                TOOL
##----------------------------------------------------------------------------
set.simu.dir <- function(simu.ext="tmp", rep=NULL, cur.wd=getwd()) {
    rst.dir <- paste(cur.wd, "/simu_", simu.ext, sep="");
    if (!is.null(rep)) {
        rst.dir <- paste(rst.dir, "/rep_", rep, sep="");
    }
    if (!file.exists(rst.dir)) {
        dir.create(rst.dir);
    }
    rst.dir
}

get.f.name <- function(simu.ext=NULL, cur.rep=NULL, lsts=NULL, suffix=".Rdata", cur.wd=getwd()) {
    fdir <- set.simu.dir(simu.ext, cur.rep, cur.wd);

    if (!is.null(lsts)) {
        fn <- paste(fdir, paste(lsts, collapse="_"), sep="/");

        if (!is.null(suffix)) {
            fn <- paste(fn, suffix, sep="");
        }
    } else {
        fn <- fdir
    }

    fn
}

##------------------------------------------------------
##
##            DATA MANIPULATION
##
##------------------------------------------------------

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

##print pvalue
prtpval <- function(pv, cut=0.0001, fmt="$%0.4f$") {
    if (pv < cut) {
        rst <- "$\\\\leq 0.0001$";
    } else {
        rst <- sprintf(fmt, pv);
    }

    rst
}

add.row <- function(dta, varname, fmt="%2.0f", nstudy=3) {
    rst <- NULL;
    for (j in 1:2) {
        for (i in 1:nstudy) {
            cur.hos  <- subset(dta, i == hospital);
            lst.d    <- rep(list(NULL),2);
            for (trt in 0:1) {
                cur.d <- subset(cur.hos, trt == group);

                if (is.null(varname)) {
                    if (2 == j) next;
                    rst <- c(rst, nrow(cur.d));
                } else {
                    cur.d          <- cur.d[[varname]];
                    lst.d[[trt+1]] <- cur.d;

                    if (1 == j) {
                        cur.r <- sprintf(fmt, mean(cur.d, na.rm=TRUE));
                    } else {
                        cur.r <- sprintf(fmt, range(cur.d, na.rm=TRUE));
                    }
                    rst <- c(rst, cur.r);

                    if (1 == j & 1 == trt) {
                        cur.test <- t.test(lst.d[[1]], lst.d[[2]])$p.value;
                        rst      <- c(rst, prtpval(cur.test));
                    }
                }
            }
        }
    }
    rst;
}

add.yesno <- function(dta, varname, les=0:1, fmt="%2.0f", nstudy=3) {
    rst <- NULL;
    for (j in 1:length(les)) {
        for (i in 1:nstudy) {
            cur.hos  <- subset(dta, i == hospital);
            for (trt in 0:1) {
                cur.d   <- subset(cur.hos, trt == group);
                cur.d   <- cur.d[[varname]];
                cur.n   <- length(which(les[j]==cur.d));
                cur.t   <- length(cur.d);
                rst     <- c(rst, cur.n, sprintf(fmt, cur.n/cur.t*100));

                if (1 == j & 1 == trt) {
                    cur.test <- chisq.test(cur.hos$group, cur.hos[[varname]])$p.value;
                    rst      <- c(rst, prtpval(cur.test));
                }
            }
        }
    }
    rst
}



##generate baseline covariate summary
tbl.baseline <- function(dta=all.data, out.f="tbl_base.tex", temp.f="tpl_base.tex") {
    ##N
    out.text <- add.row(dta, NULL, nrow);
    out.text <- c(out.text, add.yesno(dta, "sex"));
    out.text <- c(out.text, add.row(dta, "age"));
    out.text <- c(out.text, add.row(dta, "BMI"));
    out.text <- c(out.text, add.yesno(dta, "asa"));
    out.text <- c(out.text, add.yesno(dta, "spo2"));
    out.text <- c(out.text, add.yesno(dta, "diagnosis"));
    out.text <- c(out.text, add.yesno(dta, "medicine"));
    out.text <- c(out.text, add.yesno(dta, "hypertension"));
    out.text <- c(out.text, add.yesno(dta, "diabetes"));
    out.text <- c(out.text, add.yesno(dta, "heart"));

    ##tex table
    sub.txt.w.no(temp.f, out.text, out.f=out.f);
}


##------------------------------------------------------
##
##            TRADITIONAL ANALYSIS
##
##------------------------------------------------------

##observed treatment effect
get.obs.effect <- function(dta=all.data, study="Study", group="Z", y="Y") {

    ##binary or continuous
    is.bin   <- all(dta[, y] %in% c(0,1));
    d.study  <- sort(unique(dta[, study]));
    trts     <- sort(unique(dta[, group]));

    rst      <- NULL;
    for (i in 1:length(d.study)) {
        s.y   <- NULL;
        s.grp <- NULL;
        for (j in 1:length(trts)) {
            inx     <- which(d.study[i] == dta[, study] & trts[j] == dta[, group]);
            cur.n   <- length(inx);
            cur.d   <- dta[inx,];
            s.y     <- c(s.y, cur.n);

            if (is.bin) {
                c.y1   <- sum(cur.d[, y]);
                cur.b  <- binom.test(c(c.y1, cur.n-c.y1));
                s.grp  <- rbind(s.grp, c(c.y1, cur.n));
                s.y    <- c(s.y, cur.b$estimate, cur.b$conf.int);
            } else {
                cur.sd <- sd(cur.d[,y]);
                cur.m  <- mean(cur.d[,y]);
                s.grp  <- rbind(s.grp, cbind(grp=trts[j], y=cur.d[, y]));
                s.y    <- c(s.y, cur.m, cur.m-1.96*cur.sd, cur.m+1.96*cur.sd);
            }
        }

        ##difference
        s.y <- c(s.y, s.y[6] - s.y[2]);
        if (is.bin) {
            s.diff <- prop.test(s.grp[,1], s.grp[,2]);
            s.y    <- c(s.y, s.diff$conf.int, s.diff$p.value);
        } else {
            s.grp  <- data.frame(s.grp);
            s.diff <- t.test(y ~ grp, data = s.grp);
            s.y    <- c(s.y, s.diff$conf.int, s.diff$p.value);
        }

        rst <- rbind(rst, c(i, s.y));
    }

    colnames(rst) <- c("study",
                       "n0", "ey0", "lb0", "ub0",
                       "n1", "ey1", "lb1", "ub1",
                       "dif", "lb", "ub", "pval");
    rst
}

##get propensity score
get.ps <- function(dta, ps.cov=PS.COV, grp="group", delta=0) {
    ps.fml    <- as.formula(paste(grp, "~", paste(ps.cov, collapse="+"), sep=""));
    glm.fit   <- glm(ps.fml, family=binomial, data=dta);
    est.ps    <- glm.fit$fitted;

    ##exponential tilting
    e.delta      <- exp(delta);
    est.ps.delta <- sapply(est.ps, function(x) {x*e.delta/(x*e.delta+1-x)});

    ##for z=1 only
    inx.z1 <- which(1 == dta[, grp]);
    est.ps[inx.z1] <- est.ps.delta[inx.z1];

    list(score=est.ps, glm.fit=glm.fit);
}

##ps stratification
## random: whether the study is a randomized study
get.ps.strat.est <- function(dta,
                             grp="group",
                             random=0,
                             ps.cov=PS.COV,
                             y="ae",
                             n.class=5,
                             delta=0) {

    ##randomized study
    ## if (1 == random) {
    ##     which.0 <- which(0 == dta[,grp]);
    ##     trt.est <- mean(dta[-which.0, y]) - mean(dta[which.0, y]);
    ##     return(trt.est);
    ## }

    ##observational study
    ps.score <- get.ps(dta, ps.cov, grp, delta=delta)$score;
    cuts     <- quantile(ps.score, seq(0,1,length=n.class+1));
    cuts[1]  <- cuts[1]-0.01;

    a.stra <- NULL;
    for (i in 1:n.class) {
        cur.break <- which(ps.score <= cuts[i+1] & ps.score > cuts[i]);
        cur.d     <- dta[cur.break,];
        t.grps    <- sort(unique(cur.d[, grp]));

        if (1 == length(t.grps))
            next;

        inx.1 <- which(t.grps[1] == cur.d[, grp]);
        cur.e <- mean(cur.d[-inx.1, y]) - mean(cur.d[inx.1,y]);

        ##cur.t   <- table(dta[cur.inx, c(grp, y)]);
        ##if (nrow(cur.t) < 2)
        ##    next;
        ## cur.e <- apply(cur.t, 1, function(x) {
        ##     est <- x[2]/sum(x);
        ##     var <- est*(1-est)/sum(x);
        ##     c(est, var);
        ## });
        a.stra <- c(a.stra, cur.e);
    }

    trt.est <- mean(a.stra);
    trt.est
}


##typical PS stratification
get.ps.effect <- function(dta=all.data,
                          study="Study",
                          group="group",
                          random="Randomized",
                          ...,
                          n.boots=0) {

    d.study  <- unique(dta[, c(study, random)]);
    trts     <- sort(unique(dta[, group]));

    rst <- NULL;
    for (i in 1:nrow(d.study)) {
        cur.s   <- dta[which(d.study[i, study] == dta[, study]),];
        cur.rst <- get.ps.strat.est(cur.s, grp=group, random=d.study[i, random], ...);

        if (n.boots > 0) {
            d.grps <- list(NULL);
            for (g in 1:length(trts)) {
                d.grps[[g]] <- cur.s[which(trts[g] == cur.s[,group]),];
            }

            bs.est <- NULL;
            for (j in 1:n.boots) {
                cur.bs <- NULL;
                for (g in 1:length(trts)) {
                    inx <- sample(1:nrow(d.grps[[g]]),
                                  nrow(d.grps[[g]]),
                                  replace=TRUE);
                    cur.bs <- rbind(cur.bs, d.grps[[g]][inx,]);
                }

                cur.est  <- get.ps.strat.est(cur.bs, grp=group, ...);
                bs.est   <- rbind(bs.est, cur.est);
            }

            cur.rst <- c(cur.rst, sd(bs.est));
        } else {
            cur.rst <- c(cur.rst, 0);
        }

        rst <- rbind(rst, c(d.study[i, study], cur.rst));
    }

    colnames(rst) <- c("study", "est", "bssd");
    rst
}


##------------------------------------------------------
##
##            HDP ANALYSIS
##
##------------------------------------------------------

##get ps score for all subjects on all studies
get.all.ps <- function(dta=all.data,
                       study="hospital",
                       group="group",
                       random="randomized",
                       ps.cov=PS.COV,
                       prefix="score",
                       delta=0,
                       take.logit=TRUE) {

    d.study <- unique(dta[,c(study,random)]);
    osd     <- NULL;
    rst     <- NULL;
    for (i in 1:nrow(d.study)) {
        ##randomized study
        if (1 == d.study[i,random])
            next;

        osd     <- c(osd, i);
        cur.h   <- subset(dta, dta[, study] == d.study[i, study]);
        cur.glm <- get.ps(cur.h, ps.cov, group, delta=delta)$glm.fit;
        cur.ps  <- predict(cur.glm, dta, type="response");

        ##check
        ##tmp.1   <- cur.glm$fitted;
        ##tmp.2   <- predict(cur.glm, cur.h, type="response");
        ##stopifnot(all(tmp.1 == tmp.2));
        rst     <- cbind(rst, cur.ps);
    }

    colnames(rst) <- paste(prefix, osd, sep="");
    if (take.logit)
        rst <- log(rst/(1-rst));

    ##add score
    dta.score <- cbind(dta, rst);
    ps        <- apply(dta.score, 1,
                       function(x) {
                           if (x[study] %in% osd) {
                               rst <- x[paste(prefix, x[study], sep="")];
                           } else {
                               rst <- NA;
                           }
                           rst
                       });
    dta.score[, "ps"] <- ps;
    ##return
    dta.score
}


##------------------------------------------------------------------------------
##
##                   SIMULATION
##
##------------------------------------------------------------------------------


##portal for all simu parameters
set.simu.par <- function(simu.ext=SIMU.EXT) {
    eval(parse(text=paste("rst <- set.simu.par.", simu.ext, "()", sep="")));
    rst
}


##get true trt effects and assignment ratio simulation
get.simu.true <- function(simu.par.lst, sizes, nPat=100000) {
    for (j in 1:length(simu.par.lst)) {
        pt.j    <- GenDataMatrix(j, simu.par.lst, nPat=nPat, simu.y=FALSE);
        muCovar <- simu.par.lst[[j]]$muCovar;
        SdCovar <- simu.par.lst[[j]]$StDevCovar;
        mu.cov  <- c(1, muCovar, muCovar[length(muCovar)]^2 + SdCovar[length(muCovar)]^2);
        coeff   <- simu.par.lst[[j]]$RegressCoeffs;


        simu.par.lst[[j]]$nPat.tot      <- sizes[j];
        simu.par.lst[[j]]$true.pz0      <- mean(pt.j[,"Z"]);
        simu.par.lst[[j]]$true.effect   <- sum(mu.cov*coeff[2,]) - sum(mu.cov*coeff[1,]);
        simu.par.lst[[j]]$is.randomized <- identical(0, simu.par.lst[[j]]$ZBeta);
    }

    simu.par.lst
}

set.simu.par.71 <- function(sizes=c(500,500,500,200,200), simu.note="full analysis") {

    reg.coeff <- rbind(c(b0=0, u=0, v1=0, v2=0, v3=1.5, v4=1.5, v5=1.5, bv52=-2),
                       c(b0=3, u=0, v1=0, v2=0, v3=1.5, v4=1.5, v5=1.5, bv52=-2));

    study.1 <- list(muCovar=rep(5, 6), StDevCovar=rep(0.5, 6), corrCovar=0.1,
                    Ysig=1,
                    ZBeta=c(a0=-20, u=0, v1=1, v2=1, v3=1, v4=1, v5=0),
                    RegressCoeffs=reg.coeff);
    study.2 <- list(muCovar=rep(5, 6), StDevCovar=rep(0.5, 6), corrCovar=0.1,
                    Ysig=1,
                    ZBeta=c(a0=-40, u=0, v1=2, v2=2, v3=2, v4=2, v5=0),
                    RegressCoeffs=reg.coeff);
    study.3 <- list(muCovar=rep(5, 6), StDevCovar=rep(0.5, 6), corrCovar=0.1,
                    Ysig=1,
                    ZBeta=c(a0=-55, u=0, v1=3, v2=3, v3=3, v4=3, v5=0),
                    RegressCoeffs=reg.coeff);
    study.4 <- list(muCovar=rep(5, 6), StDevCovar=rep(0.1, 6), corrCovar=0.1,
                    Ysig=1, ZBeta=0, RegressCoeffs=reg.coeff);

    study.5 <- list(muCovar=rep(5, 6), StDevCovar=rep(0.1, 6), corrCovar=0.1,
                    Ysig=1, ZBeta=0, RegressCoeffs=reg.coeff);

    SIMU.PAR.LST <- get.simu.true(list(study.1, study.2, study.3, study.4),
                                  sizes);

    SIMU.PS.COV  <- paste("X", 1:4, sep="");
    SIMU.SCE     <- 71;
    SIMU.HDP.COV <- paste("score", 1:3, sep="");

    SIMU.OPS     <- HDP.OPS;
    SIMU.NOTE    <- simu.note;


    ##return
    rst <- make.list(environment());
    rst
}

set.simu.par.72 <- function() {
    rst <- set.simu.par.71(simu.note="no sharing");
    rst$SIMU.OPS$mcmc.eps <- 0;
    rst$SIMU.OPS$eps      <- 1;
    rst
}

set.simu.par.73 <- function() {
    rst <- set.simu.par.71(simu.note="no sharing, specific score");
    rst$SIMU.HDP.COV      <- "score1";
    rst$SIMU.OPS$mcmc.eps <- 0;
    rst$SIMU.OPS$eps      <- 1;
    rst
}

set.simu.par.12 <- function(sizes=c(500,500,500,200,200), simu.note="") {

    reg.coeff <- CONST.Y*rbind(c(b0=0, u=0, v1=0, v2=0, v3=1.5, v4=1.5, v5=1.5, bv52=-2),
                               c(b0=-15, u=0, v1=0, v2=0, v3=1, v4=1, v5=1, bv52=-1));

    study.1 <- list(muCovar=rep(5, 6), StDevCovar=rep(3, 6), corrCovar=0.1,
                    Ysig=1,
                    ZBeta=c(a0=-3, u=0, v1=0.2, v2=0.2, v3=0.1, v4=0.1, v5=0),
                    RegressCoeffs=reg.coeff);
    study.2 <- list(muCovar=rep(5, 6), StDevCovar=rep(3, 6), corrCovar=0.1,
                    Ysig=1,
                    ZBeta=c(a0=-4, u=0, v1=0.3, v2=0.3, v3=0.1, v4=0.1, v5=0),
                    RegressCoeffs=reg.coeff);
    study.3 <- list(muCovar=rep(5, 6), StDevCovar=rep(3, 6), corrCovar=0.1,
                    Ysig=1,
                    ZBeta=c(a0=-4, u=0, v1=0.2, v2=0.2, v3=0.3, v4=0.3, v5=0),
                    RegressCoeffs=reg.coeff);
    study.4 <- list(muCovar=rep(5, 6), StDevCovar=rep(1, 6), corrCovar=0.1,
                    Ysig=1, ZBeta=0, RegressCoeffs=reg.coeff);
    study.5 <- list(muCovar=rep(5, 6), StDevCovar=rep(1, 6), corrCovar=0.1,
                    Ysig=1, ZBeta=0, RegressCoeffs=reg.coeff);

    SIMU.PAR.LST <- get.simu.true(list(study.1, study.2, study.3, study.4, study.5), sizes);
    SIMU.PS.COV  <- paste("X", 1:4, sep="");
    SIMU.NOTE    <- simu.note;

    ##return
    rst <- make.list(environment());
    rst
}

set.simu.par.22 <- function() {
    rst <- set.simu.par.12(sizes = c(1000,1000,1000, 400, 400));
}


set.simu.par.13 <- function(sizes=c(500,500,500,200,200), simu.note="") {

    reg.coeff <- CONST.Y*rbind(c(b0=0,   u=0, v1=0, v2=0, v3=1.5, v4=1.5, v5=1.5, bv52=-2),
                               c(b0=-15, u=0, v1=0, v2=0, v3=1,   v4=1,   v5=1,   bv52=-1));

    study.1 <- list(muCovar=rep(0, 6), StDevCovar=rep(3, 6), corrCovar=0.1,
                    Ysig=1,
                    ZBeta=c(a0=0, u=0, v1=0.2, v2=0.2, v3=0.1, v4=0.1, v5=0),
                    RegressCoeffs=reg.coeff);
    study.2 <- list(muCovar=rep(5, 6), StDevCovar=rep(3, 6), corrCovar=0.1,
                    Ysig=1,
                    ZBeta=c(a0=-4, u=0, v1=0.3, v2=0.3, v3=0.1, v4=0.1, v5=0),
                    RegressCoeffs=reg.coeff);
    study.3 <- list(muCovar=rep(10, 6), StDevCovar=rep(3, 6), corrCovar=0.1,
                    Ysig=1,
                    ZBeta=c(a0=-10, u=0, v1=0.2, v2=0.2, v3=0.3, v4=0.3, v5=0),
                    RegressCoeffs=reg.coeff);
    study.4 <- list(muCovar=rep(5, 6), StDevCovar=rep(1, 6), corrCovar=0.1,
                    Ysig=1, ZBeta=0, RegressCoeffs=reg.coeff);
    study.5 <- list(muCovar=rep(5, 6), StDevCovar=rep(1, 6), corrCovar=0.1,
                    Ysig=1, ZBeta=0, RegressCoeffs=reg.coeff);

    SIMU.PAR.LST <- get.simu.true(list(study.1, study.2, study.3, study.4, study.5), sizes);
    SIMU.PS.COV  <- paste("X", 1:4, sep="");
    SIMU.NOTE    <- simu.note;

    ##return
    rst <- make.list(environment());
    rst
}

set.simu.par.23 <- function() {
    rst <- set.simu.par.13(sizes = c(1000,1000,1000, 400, 400));
}

set.simu.par.14 <- function(sizes=c(500,500,500,200,200), simu.note="") {

    reg.coeff <- CONST.Y*rbind(c(b0=0, u=1, v1=0, v2=0, v3=1.5, v4=1.5, v5=1.5, bv52=-2),
                               c(b0=-15, u=1, v1=0, v2=0, v3=1, v4=1, v5=1, bv52=-1));

    study.1 <- list(muCovar=rep(0, 6), StDevCovar=rep(3, 6), corrCovar=0.1,
                    Ysig=1,
                    ZBeta=c(a0=0, u=0.3, v1=0.2, v2=0.2, v3=0.1, v4=0.1, v5=0),
                    RegressCoeffs=reg.coeff);
    study.2 <- list(muCovar=rep(5, 6), StDevCovar=rep(3, 6), corrCovar=0.1,
                    Ysig=1,
                    ZBeta=c(a0=-5, u=0.3, v1=0.3, v2=0.3, v3=0.1, v4=0.1, v5=0),
                    RegressCoeffs=reg.coeff);
    study.3 <- list(muCovar=rep(10, 6), StDevCovar=rep(3, 6), corrCovar=0.1,
                    Ysig=1,
                    ZBeta=c(a0=-12, u=0.3, v1=0.2, v2=0.2, v3=0.3, v4=0.3, v5=0),
                    RegressCoeffs=reg.coeff);
    study.4 <- list(muCovar=rep(5, 6), StDevCovar=rep(1, 6), corrCovar=0.1,
                    Ysig=1, ZBeta=0, RegressCoeffs=reg.coeff);
    study.5 <- list(muCovar=rep(5, 6), StDevCovar=rep(1, 6), corrCovar=0.1,
                    Ysig=1, ZBeta=0, RegressCoeffs=reg.coeff);

    SIMU.PAR.LST <- get.simu.true(list(study.1, study.2, study.3, study.4, study.5), sizes);
    SIMU.PS.COV  <- paste("X", 1:4, sep="");
    SIMU.NOTE    <- simu.note;

    ##return
    rst <- make.list(environment());
    rst
}

set.simu.par.24 <- function() {
    rst <- set.simu.par.14(sizes = c(1000,1000,1000, 400, 400));
}

##assign treatment
set.trt.z <- function(covx, zbeta, trt.val=TRT.NUMBER) {

    n.x    <- ncol(covx) + 1;
    n.beta <- length(zbeta);

    if (n.beta < n.x) {
        zbeta <- c(zbeta, rep(0, n.x-n.beta));
    }

    covx.beta <- apply(covx, 1, function(x){sum(c(1,x)*zbeta)});
    pz        <- expit(covx.beta);
    simu.rst  <- rbinom(length(pz), 1, pz);
    rst       <- trt.val[simu.rst+1];
    rst
}


##set mean
set.xbeta <- function(covx, beta) {
    trt <- covx[1];
    ##add quaratic effect
    c.x <- covx[-1];
    c.x <- c(1, c.x, c.x[length(c.x)]^2);
    sum(c.x * beta[trt+1,]);
}

GenDataMatrix <- function(j, par.list, Study=NULL, cov.x=NULL, nPat=NULL, z.val=NULL, simu.y=TRUE) {

    ##expose parameters in par.list
    make.global(par.list[[j]], environment());

    if (!exists("is.randomized"))
        is.randomized <- -1;

    ##study
    if (is.null(Study)) Study <- j;

    ##sample size
    if (is.null(nPat)) nPat <- nPat.tot;

    ##Xbeta
    n.x <- length(StDevCovar);

    if (is.null(cov.x)) {
        Vars     <- StDevCovar*StDevCovar;
        CovarMat <- matrix(0, n.x, n.x);
        for (i in 1:n.x) {
            CovarMat[i,i] <- Vars[i];
            for (j in (i+1):n.x) {
                if ( j > n.x) break;
                CovarMat[i, j] <- corrCovar*StDevCovar[i]*StDevCovar[j];
                CovarMat[j, i] <- CovarMat[i, j];
            }
        }

        ##simulate continuous X
        COV.X <- rmvnorm(nPat, mean=muCovar, sigma=CovarMat);

        ##treatment
        if (!is.null(z.val)) {
            Z <- z.val;
        } else {
            Z <- set.trt.z(COV.X, ZBeta);
        }

        ##covariates X without intercept
        X <- cbind(Z, COV.X);
    } else {
        X <- NULL;
        for (i in 1:length(cov.x)) {
            X <- cbind(X, rep(cov.x[i], nPat));
        }
    }

    ##simulate Y
    if (simu.y) {
        ##beta0 + beta * x
        xbeta    <- apply(X, 1, function(x) {set.xbeta(x, RegressCoeffs)});
        Y        <- rnorm(nPat, xbeta, Ysig);
    } else {
        Y <- NA;
    }

    ##return
    Data           <- cbind(1:nPat, Study, is.randomized, Y, X);
    colnames(Data) <- c("pid", "Study", "Randomized", "Y", "Z",
                        paste("X", 1:(ncol(X)-1), sep=""));

    Data
}

##------------------------------------------------------------------------------
##
##              COMBINE SIMULATION RESULTS
##
##------------------------------------------------------------------------------
simu.post.combine <- function(simu.ext, cmb.reps=NULL, chk.exist=FALSE) {
    fsum <- function(mat, true.e) {
        rst <- apply(cbind(true.e, mat),
                     1,
                     function(x) {
            c((x[2]-x[1])^2,
               x[1]>=x[3] & x[1]<=x[4]);
                     });
        cbind(mat, t(rst));
    }

    if (is.null(cmb.reps)) {
        cmb.reps <- 1:N.REPS;
    }

    ##simu results files
    f.rst <- paste("simu_rst_", simu.ext, ".Rdata", sep="");
    if (file.exists(f.rst) & chk.exist) {
        print(paste(f.rst, " exists...", sep=""));
        return(NULL);
    }

    ##load simu setting
    load(get.f.name(simu.ext=simu.ext, lsts=c("simu_par")));
    true.e <- NULL;

    for (i in 1:length(SIMU.PAR.LST)) {
        true.e <- c(true.e, SIMU.PAR.LST[[i]]$true.effect);
    }

    ##combine
    all.rst <- NULL;
    for (cr in 1:length(cmb.reps)) {
        r     <- cmb.reps[cr];
        f.fit <- get.f.name(simu.ext=simu.ext, cur.r = cr, lsts=c("result"));
        if (!file.exists(f.fit)) {
            print(paste(f.fit, " does not exist...", sep=""));
            next;
        }

        load(f.fit);
        cur.obs <- rst.obs[, c('dif', 'lb', 'ub')];
        cur.obs <- fsum(cur.obs, true.e);
        cur.obs <- cbind("obs",
                         cr,
                         1:length(SIMU.PAR.LST),
                         0,
                         cur.obs);

        all.rst <- rbind(all.rst, cur.obs);

        for (j in 1:length(DELTA.ALL)) {
            ##ps
            cur.ps  <- rst.ps[[j]][,2:3];
            cur.ps  <- cbind(cur.ps[,1],
                             cur.ps[,1]-1.96*cur.ps[,2],
                             cur.ps[,1]+1.96*cur.ps[,2]);
            cur.ps  <- fsum(cur.ps, true.e);
            cur.ps  <- cbind("ps",
                             cr,
                             1:length(SIMU.PAR.LST),
                             DELTA.ALL[j],
                             cur.ps);
            all.rst <- rbind(all.rst, cur.ps);

            ##hdp
            cur.hdp <- fsum(rst.hdp[[j]], true.e);
            cur.hdp <- cbind("hdp",
                             cr,
                             1:length(SIMU.PAR.LST),
                             DELTA.ALL[j],
                             cur.hdp);
            all.rst <- rbind(all.rst, cur.hdp);
        }
    }

    ##summarize
    colnames(all.rst) <- c("method", "rep", "study", "delta", "est",
                           "lb",     "ub",  "se",    "inci");

    allrst  <- data.frame(all.rst);
    sum.rst <- sqldf("select method, study, delta,
                      avg(est) as est,
                      avg(se) as mse,
                      avg(inci) as inci
                      from allrst
                      group by method, study, delta");

    ##save results
    print(paste(f.rst, " generated...."));
    save(SIMU.PAR.LST, SIMU.NOTE, SIMU.PS.COV,
         cmb.reps, all.rst, sum.rst,
         file=f.rst);
}

