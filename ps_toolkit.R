library(hdpmn);
library(sqldf);
library(mvtnorm);

##library(parallel);
##require(mice);
##library(MatchIt);
##library(cgtool);

source("./ps_simu_sce.R");
source("./rwe_private.R");
source("./rwe_simu.R");

##----------------------------------------------------------------------------
##                CONSTANTS
##----------------------------------------------------------------------------
##labels of treatment arms
TRT.NUMBER  <- c(0,1);
##sensitivity
DELTA.ALL <- c(0);
##ignore existing
KEEP.EXISTING <- FALSE;
##propensity score breaks
PS.CLASS <- 5;
##number bootstraps
## N.BOOTS <- 2;
##HDP prediction
PRED.R  <- 0;
##ratio for Y to make Y comparable to scores
CONST.Y <- 1;

## batch set up
N.REPS       <- 2000; ##number of total replications
SIMU.EXT.ALL <- c(781:785); ##simulation scenarios

##----------------------------------------------------------------------------
##                TOOL
##----------------------------------------------------------------------------
##tic toc count time
tic <- function() {
    assign("PROCSTARTTIME", proc.time(), ".GlobalEnv");
}

toc <- function() {
    if (!exists("PROCSTARTTIME")) return;
    print(proc.time() - PROCSTARTTIME);
}

expit <- function(x) {
    rst <- exp(x)/(1+exp(x));
    rst[which(is.nan(rst))] <- 1;
    rst
}

##substitue text files with numbers, the sub.str is the string to be placed by nubmers
sub.txt.w.no <- function(template.f, numbers, out.f="rst.txt", sub.str="AA") {
    if (!file.exists(template.f)) {
        return;
    }
    ##read
    tpla <- readChar(template.f, file.info(template.f)$size);

    ##subsitute
    for (i in 1:length(numbers)) {
        tpla <- sub(sub.str, numbers[i], tpla);
    }

    ##write out
    write(tpla, file=out.f);
}

##make variables in a list global variables
make.global <- function(alist, dest.env='.GlobalEnv') {
    for (i in 1:length(alist)) {
        assign(names(alist[i]), alist[[i]], dest.env );
    }
}

##choose all captialized variables and make them a list
make.list <- function(par.env) {
    rst     <- list();
    objlst  <- objects(par.env);

    for (o in 1:length(objlst)) {
        curo <- objlst[o];
        if (toupper(curo) == curo) {
            rst[[curo]] <- get(curo, envir=par.env);
        }
    }
    rst
}


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

get.f.name <- function(simu.ext=NULL, cur.rep=NULL, lsts=NULL,
                       suffix=".Rdata", cur.wd=getwd()) {

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
    out.text <- c(out.text, add.row(dta,   "age"));
    out.text <- c(out.text, add.row(dta,   "BMI"));
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

##exponential tilting
get.expt <- function(ps, delta=0) {
    ##exponential tilting
    e.delta  <- exp(delta);
    ps.delta <- sapply(ps, function(x) {
        x*e.delta/(x*e.delta+1-x)
    });

    ps.delta;
}


##observed treatment effect
get.obs.effect <- function(dta=all.data, study="Study", group="Z", y="Y") {
    ##binary or continuous
    is.bin   <- all(dta[, y] %in% c(0,1));
    d.study  <- sort(unique(dta[, study]));
    trts     <- sort(unique(dta[, group]));

    rst <- NULL;
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
get.ps <- function(dta, ps.cov=PS.COV, grp="group", delta=0, ps.fml=NULL) {
    if (is.null(ps.fml)) {
        ps.fml <- as.formula(paste(grp, "~", paste(ps.cov, collapse="+"), sep=""));
    }

    glm.fit   <- glm(ps.fml, family=binomial, data=dta);
    est.ps    <- glm.fit$fitted;

    list(score   = get.expt(est.ps, delta),
         glm.fit = glm.fit);
}

##ps stratification
## random: whether the study is a randomized study
get.ps.strat.est <- function(dta,
                             grp="group",
                             y="ae",
                             ps.cov=PS.COV,
                             n.class=5,
                             delta=0,
                             ...) {

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

##ps match
get.ps.match.est <- function(dta,
                             grp="group",
                             random=0,
                             ps.cov=PS.COV,
                             y="ae",
                             n.class=5,
                             delta=0) {

    ps.fml    <- as.formula(paste(grp, "~", paste(ps.cov, collapse="+"), sep=""));
    mtd       <- matchit(ps.fml, data = dta, method="nearest");
    mtd.dta   <- match.data(mtd);
    inx.0     <- which(TRT.NUMBER[1] == mtd.dta[,grp]);
    trt.est   <- mean(dta[-inx.0, y]) - mean(dta[inx.0, y]);
    trt.est
}


##typical PS stratification
get.ps.effect <- function(dta     = all.data,
                          study   = "Study",
                          group   = "group",
                          ...,
                          n.boots=100) {

    studies  <- sort(unique(dta[, study]));
    trts     <- sort(unique(dta[, group]));

    rst <- NULL;
    for (s in studies) {
        cur.s   <- dta[which(s == dta[, study]),];
        cur.rst <- get.ps.strat.est(cur.s, grp=group, ...);

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

        rst <- rbind(rst, c(s, cur.rst));
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
                       ps.ignore.random=FALSE,
                       ps.cov=PS.COV,
                       prefix="score",
                       delta=0,
                       take.logit=TRUE) {

    d.study <- unique(dta[,c(study,random)]);

    osd     <- NULL;
    rst     <- NULL;
    for (i in 1:nrow(d.study)) {
        ##randomized study
        if (1 == d.study[i,random] &
            ps.ignore.random)
            next;

        osd       <- c(osd, d.study[i, study]);
        cur.h     <- subset(dta, dta[, study] == d.study[i, study]);
        cur.glm   <- get.ps(cur.h, ps.cov=ps.cov, grp=group)$glm.fit;
        cur.ps    <- predict(cur.glm, dta, type="response");
        cur.extps <- get.expt(cur.ps, delta);
        rst       <- cbind(rst, cur.extps);
        ##check
        ##tmp.1   <- cur.glm$fitted;
        ##tmp.2   <- predict(cur.glm, cur.h, type="response");
        ##stopifnot(all(tmp.1 == tmp.2));
    }

    colnames(rst) <- paste(prefix, osd, sep="");
    if (take.logit)
        rst <- log(rst/(1-rst));

    ##add score
    dta.score <- cbind(dta, rst);
    ps        <- apply(dta.score, 1, function(x) {
                           if (x[study] %in% osd) {
                               rst <- x[paste(prefix, x[study], sep="")];
                           } else {
                               rst <- NA;
                           }
                           rst
                       });
    dta.score[, "ps"] <- as.numeric(ps);

    ##return
    dta.score
}

## get treatment effect based on predicted valutes
get.trt.effect <- function(studies = STUDIES,
                           trt.number = TRT.NUMBER,
                           y.type  = "continuous",
                           y.sd    = 1,
                           pred.r  = PRED.R,
                           sub.inx = NULL,
                           wd      = "./") {

    old.wd <- setwd(wd);

    n.trt <- length(trt.number);
    rst   <- NULL;
    for (sy in studies) {
        pred.g <- NULL;
        for (g in 1:n.trt) {
            simu.dir        <- set.simu.dir(paste("grp", trt.number[g], sep=""));
            cur.g           <- read.table(paste(simu.dir, "/spred_", sy, "_", pred.r, sep=""));
            colnames(cur.g) <- c("ID", "Iter", "Pred");
            pred.g[[g]]     <- cur.g;
        }

        all.iter <- unique(pred.g[[1]][,"Iter"]);
        all.sy   <- NULL;
        for (k in all.iter) {
            py <- NULL;
            for (g in 1:length(trt.number)) {
                cur.p <- subset(pred.g[[g]], Iter == k);
                if (!is.null(sub.inx)) {
                    cur.p <- subset(cur.p, ID %in% sub.inx);
                }
                py <- cbind(py, cur.p[order(cur.p$ID), 'Pred']);
            }

            if ("binary" == y.type) {
                py <- py > 0;
            } else {
                ##scale back
                py <- y.sd[min(length(y.sd), sy)] * py;
            }

            ## use different subjects for prediction to reduce correlation
            cur.effect <- apply(py, 2, mean);
            all.sy <- rbind(all.sy, cur.effect);
        }

        py.diff <- all.sy[,2] - all.sy[,1];
        rst     <- rbind(rst,
                         c(mean(py.diff),
                           quantile(py.diff, c(0.025, 0.975))));
    }

    setwd(old.wd);
    rst
}

##estimate for DP models
get.dp.all <- function(hdp.ps, pred.pts, studies=STUDIES, simu.ops=SIMU.OPS,
                       pred.r=PRED.R, trt.number=TRT.NUMBER, center.y = TRUE) {

    ##single study
    simu.ops$mcmc.eps <- 0;
    simu.ops$eps      <- 1;

    all.sd <- NULL;
    for (sy in studies) {
        pred.g    <- NULL;
        cur.sdata <- subset(hdp.ps, Study == sy);

        y.sd <- 1;
        if (center.y) {
            y.sd            <- sd(cur.sdata[,"Y"]);
            cur.sdata[,"Y"] <- (cur.sdata[,"Y"] - mean(cur.sdata[,"Y"]))/y.sd;
        }

        all.sd <- c(all.sd, y.sd);

        for (grp in trt.number) {
            ##go to grp folder
            cur.owd <- getwd();
            setwd(set.simu.dir(paste("grp", grp, sep="")));

            ##assign data
            cur.d       <- subset(cur.sdata, grp == Z);
            cur.d$Study <- 1;

            ##fitting
            cur.z     <- as.matrix(cur.d[, c("Y", "ps")]);
            cur.study <- cur.d$Study;

            do.call(hdpmn::hdpmn, c(list(Z=cur.z, study=cur.study),
                             px=ncol(cur.z)-1,
                             simu.ops));

            ##predicting
            cur.pts  <- pred.pts[which(pred.pts[,'Study'] == sy),];
            cur.pred <- hdpmn::R.hdpmnPredict(j=1,
                                       r=pred.r,
                                       nsim=0,
                                       idx.x=2,
                                       X=cur.pts[, "ps", drop=FALSE]);

            ##rename files
            file.rename(from = paste("pred_1_", pred.r, sep=""),
                        to = paste("spred_", sy, "_", pred.r, sep=""));

            ##go back
            setwd(cur.owd);
        }
    }

    all.sd;
}

##estimate for hdp models
get.hdp.all <- function(hdp.ps, pred.pts,
                        studies=STUDIES, simu.ops=SIMU.OPS,
                        pred.r=PRED.R, hdp.cov=HDP.COV,
                        trt.number=TRT.NUMBER, center.y = TRUE) {

    n.trt <- length(trt.number);
    y.sd  <- 1;
    ##remove mean of y by study grp
    if (center.y) {
        ## oldhp   <- data.frame(hdp.ps);
        ## ymeans  <- sqldf("select Study, Z,
        ##                   avg(Y) as __M__,
        ##                   stdev(Y) as __SD__
        ##                   from oldhp group by Study, Z
        ##                   order by Study, Z");
        ## oldhp   <- sqldf("select a.*, b.__M__, b.__SD__ from oldhp a
        ##                   left join ymeans b on (a.Study = b.Study and a.Z = b.Z)");
        ## oldhp$Y <- (oldhp$Y - oldhp[["__M__"]]) / oldhp[["__SD__"]];
        y.sd      <- sd(hdp.ps$Y);
        hdp.ps$Y  <- (hdp.ps$Y - mean(hdp.ps$Y))/y.sd;
    }

    for (grp in trt.number) {
        ##go to grp folder
        cur.owd <- getwd();
        setwd(set.simu.dir(paste("grp", grp, sep="")));

        ##assign data
        cur.d  <- subset(hdp.ps, grp == Z);

        ##sort by study no first
        cur.d <- cur.d[order(cur.d$Study),];

        ##fitting
        cur.z     <- as.matrix(cur.d[, c("Y", hdp.cov)]);
        cur.study <- sapply(cur.d$Study, function(x) {which(x == studies)});

        do.call(hdpmn,
                c(list(Z=cur.z, study=cur.study),
                  px=ncol(cur.z)-1,
                  simu.ops));

        ##predicting
        for (sy in studies) {
            cur.pts  <- pred.pts[which(pred.pts[,'Study'] == sy),];
            cur.pred <- R.hdpmnPredict(j=which(sy == studies),
                                       r=pred.r,
                                       nsim=0,
                                       idx.x=1+1:length(hdp.cov),
                                       X=cur.pts[, hdp.cov, drop=FALSE]);
            ##rename files
            file.rename(from = paste("pred_", which(sy == studies), "_", pred.r, sep=""),
                        to   = paste("spred_", sy, "_", pred.r, sep=""));
        }

        ##go back
        setwd(cur.owd);
    }

    ## ##post analysis sum hdp
    ## rst <- NULL;
    ## for (sy in studies) {
    ##     pred.g <- NULL;
    ##     for (g in 1:n.trt) {
    ##         simu.dir        <- set.simu.dir(paste("grp", trt.number[g], sep=""));
    ##         cur.g           <- read.table(paste(simu.dir, "/pred_", sy, "_", pred.r, sep=""));
    ##         colnames(cur.g) <- c("ID", "Iter", "Pred");
    ##         pred.g[[g]]     <- cur.g;
    ##     }

    ##     all.iter <- unique(pred.g[[1]][,"Iter"]);
    ##     all.sy   <- NULL;
    ##     for (k in all.iter) {
    ##         py <- NULL;
    ##         for (g in 1:length(trt.number)) {
    ##             cur.p <- subset(pred.g[[g]], Iter == k);
    ##             py    <- cbind(py, cur.p[order(cur.p$ID), 'Pred']);
    ##         }

    ##         ##transfer back
    ##         all.sy <- rbind(all.sy, apply(py, 2, mean));
    ##         ##all.sy <- rbind(all.sy, py);
    ##     }

    ##     sy.diff <- y.sd*(sample(all.sy[,2],100000, replace = TRUE) -
    ##                      sample(all.sy[,1],100000, replace = TRUE));
    ##     rst     <- rbind(rst,
    ##                      c(mean(sy.diff),
    ##                        quantile(sy.diff, c(0.025, 0.975))));
    ## }

    ## rst
    y.sd
}

##------------------------------------------------------------------------------
##
##                   SIMULATION
##
##------------------------------------------------------------------------------
get.simu.data <- function(study, par.list) {
    cur.par        <- par.list[[study]];
    rst            <- do.call(rweSimuTwoArm, cur.par)$simu.data;
    rst$Study      <- study;
    rst$Randomized <- as.numeric(all(0 == cur.par$regCoeff.z));
    rst
}


##portal for all simu parameters
set.simu.par <- function(simu.ext=SIMU.EXT) {
    eval(parse(text=paste("rst <- set.simu.par.", simu.ext, "()", sep="")));
    rst
}


##get true trt effects, covariate imbalance and assignment ratio simulation
## get.simu.true <- function(simu.par.lst, nPat=100000, corCov = 0.1) {
##     for (j in 1:length(simu.par.lst)) {
##         print(j);
##         simu.par.lst[[j]]$corCov <- corCov;
##         cur.par      <- simu.par.lst[[j]];
##         cur.par$nPat <- nPat;
##         b0.ysig      <- do.call(rweSimuTwoArm,
##                                 c(cur.par, do.simu=FALSE, bin.mu = 0.5))$b0ysig;
##         simu.par.lst[[j]] <- c(simu.par.lst[[j]], b0.ysig[1], b0.ysig[2]);
##     }

##     simu.par.lst
## }

## get unbalance of the covariates
get.unbalance <- function(simu.par.lst, nPat=100000, seed=10000, ...) {
    set.seed(seed);
    rst <- NULL;
    for (j in 1:length(simu.par.lst)) {
        par.study <- simu.par.lst[[j]];
        par.study$nPat <- nPat;
        nx             <- length(par.study$muCov);
        pt.j           <- do.call(rweSimuTwoArm, par.study);

        ##unbalance
        unb <- rweUnbalance(nPat, pts = pt.j, ...);
        rst <- rbind(rst, cbind(Study = paste("Study", j, sep=" "),
                                unb));
    }
    rst
}

##assign treatment
set.trt.z <- function(covx, zbeta, trt.val=TRT.NUMBER) {

    ##n.x    <- ncol(covx) + 1;
    ##n.beta <- length(zbeta);
    ##if (n.beta < n.x) {
    ##    zbeta <- c(zbeta, rep(0, n.x-n.beta));
    ##}

    covx.beta <- apply(covx, 1, function(x){sum(c(1,x)*zbeta)});
    pz        <- expit(covx.beta);
    simu.rst  <- rbinom(length(pz), 1, pz);
    rst       <- trt.val[simu.rst+1];
    rst
}


##set mean
set.xbeta <- function(covx, effect, beta) {
    ##trt <- covx[1];
    ##add quaratic effect
    ##c.x <- covx[-1];
    ##c.x <- c(1, c.x, c.x[length(c.x)]^2);
    sum(covx * c(effect,beta));
}

##------------------------------------------------------------------------------
##
##                   SIMULATION BATCH
##
##------------------------------------------------------------------------------
get.simu.par <- function(simu.ext.all=SIMU.EXT.ALL) {
    for (s in 1:length(simu.ext.all)) {
        SIMU.EXT <- simu.ext.all[s];
        print(SIMU.EXT);
        ##simulation parameters
        make.global(set.simu.par(SIMU.EXT));
        save(SIMU.PAR.LST, SIMU.PS.COV, SIMU.EXT, SIMU.HDP.COV,
             file = get.f.name(simu.ext=SIMU.EXT, lsts=c("simu_par")));
    }
}

## simulate data for all scenarios
get.simu.all <- function(simu.ext.all=SIMU.EXT.ALL, n.reps = 1:N.REPS, seed = NULL) {

    if (!is.null(seed)) {
        set.seed(seed);
    }

    for (s in 1:length(simu.ext.all)) {
        SIMU.EXT <- simu.ext.all[s];
        par.file <- get.f.name(simu.ext=SIMU.EXT, lsts=c("simu_par"));

        if (!file.exists(par.file)) {
            get.simu.par(simu.ext.all = SIMU.EXT);
        }

        ## load simu parameters
        load(file = par.file);

        ##simu data
        for (i in n.reps) {
            print(i);
            cur.f <- get.f.name(simu.ext=SIMU.EXT,
                                cur.rep=i,
                                lsts=c("datarct"));

            print(paste("Simu:", cur.f, sep=""));
            simu.data <- NULL;
            for (j in 1:length(SIMU.PAR.LST)) {
                pt.study.j  <- get.simu.data(j, SIMU.PAR.LST);
                simu.data   <- rbind(simu.data, pt.study.j);
            }

            write.table(simu.data, file=cur.f,
                        sep=" ", quote=F, row.names=F, col.names=T);

        }
    }
}

##simulation analysis
batch.simu.ana <- function(rep, simu.ext, delta = 0) {

    ##MCMC
    SIMU.OPS <- list(n.iter=5000,
                     n.discard=2000,
                     n.batch=20,
                     mcmc.eps=1,
                     eps=0.5,
                     m.prior=1,
                     verbose=3,
                     B.prior=1, S.prior=1, alpha.prior=1,
                     alpha=20,
                     q=10,
                     cc=10);

    r <- rep;
    ##----------load simu parameters------------
    load(file=get.f.name(simu.ext, lsts=c("simu_par")));

    ##print info
    print(paste( "---CUR SIM.EXT:--", simu.ext, "---CUR REP:",r, "---", sep=""));

    ##------------read data----------------------
    rct.data   <- read.table(get.f.name(simu.ext=simu.ext, cur.rep=r, lsts=c("datarct")),
                             header=T);
    STUDIES    <- sort(unique(rct.data$Study));
    ##rct.data$Y <- rct.data$Y;

    ##------------analysis-----------------------
    f.rst <- get.f.name(simu.ext=simu.ext, cur.rep=r, lsts=c("result"));
    if (file.exists(f.rst) & KEEP.EXISTING)
        return();

    ##set work dir
    old.wd <- getwd();
    setwd(get.f.name(simu.ext, cur.rep=r));

    ##---1. observed results
    rst.obs <- get.obs.effect(rct.data, study="Study", group="Z", y="Y");

    rst.ps   <- list(NULL);
    rst.hdp  <- list(NULL);
    rst.dp   <- list(NULL);

    ##---2. propensity score analysis
    rst.ps  <- get.ps.effect(rct.data,
                             study="Study",
                             group="Z",
                             random="Randomized",
                             ps.cov=SIMU.PS.COV,
                             y="Y",
                             n.class=PS.CLASS,
                             delta=delta,
                             n.boots=0);

    ##---3. dp and hdp model analysis
    ##all ps scores for observational studies
    hdp.ps <- get.all.ps(rct.data,
                         study="Study",
                         group="Z",
                         random="Randomized",
                         ps.cov=SIMU.PS.COV,
                         delta=delta,
                         take.logit=TRUE);

    ## hdp fit and predict
    HDP.COV  <- SIMU.HDP.COV[which(SIMU.HDP.COV %in% colnames(hdp.ps))];
    pred.pts <- as.matrix(hdp.ps[, c("Study", HDP.COV), drop=FALSE]);
    sd.hdp   <- get.hdp.all(hdp.ps,
                            pred.pts,
                            studies = STUDIES,
                            hdp.cov = HDP.COV,
                            center.y = TRUE,
                            simu.ops= SIMU.OPS);
    rst.hdp  <- get.trt.effect(y.sd = sd.hdp,
                               studies = STUDIES);

    ##fit and predict
    pred.pts <- as.matrix(hdp.ps[, c("Study", "ps"), drop=FALSE]);
    sd.dp    <- get.dp.all(hdp.ps,
                           pred.pts,
                           studies  = STUDIES,
                           center.y = TRUE,
                           simu.ops = SIMU.OPS);
    rst.dp   <- get.trt.effect(studies = STUDIES,
                               y.sd    = sd.dp);


    ##return to previous dir
    setwd(old.wd);

    ##-----------save results------------------
    save(rst.obs, rst.ps, rst.hdp, rst.dp, delta,
         SIMU.PAR.LST, HDP.COV, file=f.rst);
}


##------------------------------------------------------------------------------
##
##              COMBINE SIMULATION RESULTS
##
##------------------------------------------------------------------------------

##plot unbalance in covariates
simu.plot.x <- function(simu.par.lst, main="") {

    f.x <- function(x) {
        x.min <- range.x[1] - (range.x[2]-range.x[1])*0.05;
        x.max <- range.x[2] + (range.x[2]-range.x[1])*0.05;
        cur.left + (x - x.min)/(x.max - x.min) * (cur.right - cur.left);
    }

    n.study  <- length(simu.par.lst);
    n.x      <- ncol(simu.par.lst[[1]]$unbalance);
    x.offset <- 0.05;
    y.offset <- 0.9;
    x.width  <- (1-x.offset)/n.study;
    y.width  <- y.offset/n.x;

    all.x <- NULL;
    for (j in 1:n.study) {
        all.x <- c(all.x, simu.par.lst[[j]]$unbalance);
    }
    range.x <- range(all.x);

    plot(NULL, xlim = c(0,1), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "", main=main);
    axis.x    <- NULL;
    axis.mark <- NULL;
    for (j in 1:n.study) {
        cur.x      <- simu.par.lst[[j]]$unbalance;
        cur.left   <- x.offset + (j-1)*x.width;
        cur.right  <- x.offset + j*x.width;
        cur.center <- (cur.left + cur.right)/2;
        text(cur.center, 1, paste("Study", j));
        lines(c(f.x(0), f.x(0)), c(0,y.offset), lty=2, col="gray");

        for (k in 1:n.x) {
            cur.y <- (k-1)*y.width + y.width/2;
            if (1 == j) {
                text(0, cur.y, colnames(cur.x)[k]);
            }

            lines(f.x(cur.x[2:3,k]), c(cur.y, cur.y), lty=1);
            points(f.x(cur.x[1,k]), cur.y, pch=15, cex=1.2);
        }

        ##
        cur.marks <- c(ceiling(range.x[1]), 0, floor(range.x[2]));
        axis.x    <- c(axis.x, f.x(cur.marks));
        axis.mark <- c(axis.mark, cur.marks)
    }
    axis(1, at=axis.x, labels=axis.mark);
}


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
        true.e <- c(true.e, SIMU.PAR.LST[[i]]$Effect);
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

        ##ps
        cur.ps  <- rst.ps[,2:3];
        cur.ps  <- cbind(cur.ps[,1],
                         cur.ps[,1]-1.96*cur.ps[,2],
                         cur.ps[,1]+1.96*cur.ps[,2]);
        cur.ps  <- fsum(cur.ps, true.e);
        cur.ps  <- cbind("ps",
                         cr,
                         1:length(SIMU.PAR.LST),
                         0,
                         cur.ps);
        all.rst <- rbind(all.rst, cur.ps);

        ##hdp
        cur.hdp <- fsum(rst.hdp, true.e);
        cur.hdp <- cbind("hdp",
                         cr,
                         1:length(SIMU.PAR.LST),
                         0,
                         cur.hdp);
        all.rst <- rbind(all.rst, cur.hdp);

        ##dp
        cur.dp <- fsum(rst.dp, true.e);
        cur.dp <- cbind("dp",
                        cr,
                        1:length(SIMU.PAR.LST),
                        0,
                        cur.dp);
        all.rst <- rbind(all.rst, cur.dp);
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
    save(SIMU.PAR.LST, SIMU.PS.COV,
         cmb.reps, all.rst, sum.rst,
         file=f.rst);
}

## print simulation parameters
simu.par.prt <- function(simu.par.lst, sce = "I",
                         template.f = "./Figures/template_simu.tex", out.f = "out.tex") {

    r.int <- function(x) {
        round(x)
    }

    r.f <- function(x) {
        round(x, digits = 2);
    }

    fill.all <- NULL;
    for (j in 1:length(simu.par.lst)) {
        simu.par <- simu.par.lst[[j]];
        cur.fill <- matrix("", 10, 15);

        cur.fill[1  ,1] <- sce;
        cur.fill[1  ,2] <- j;
        cur.fill[   ,3]  <- paste("$\\\\mu_{", 1:10, "}$", sep = "");
        cur.fill[   ,4]  <- r.f(simu.par$muCov);
        cur.fill[   ,5]  <- paste("$\\\\sigma_{", 1:10, "}$", sep = "");
        cur.fill[   ,6]  <- r.f(simu.par$sdCov);
        cur.fill[1:6,7]  <- paste("$\\\\alpha_", 1:6, "$", sep = "");
        cur.fill[1:6,8]  <- r.f(simu.par$regCoeff.z);
        cur.fill[3:10,9]   <- paste("$\\\\beta_{1,", 3:10, "}$", sep = "");
        cur.fill[3:10,10]  <- r.f(simu.par$regCoeff.y[2:9]);
        cur.fill[3:10,11]  <- paste("$\\\\beta_{2,", 3:10, "}$", sep = "");
        cur.fill[3:10,12]  <- r.f(simu.par$regCoeff.y[10:17]);

        cur.fill[1,13]  <- paste("$\\\\beta_{1,0}$");
        cur.fill[1,14]  <- r.f(simu.par$regCoeff.y[1]);
        cur.fill[2,13]  <- paste("$\\\\beta_{0,0}$");
        cur.fill[2,14]  <- r.f(0);
        cur.fill[3,13]  <- paste("$\\\\alpha_0$");
        cur.fill[3,14]  <- r.f(simu.par$b0);
        cur.fill[4,13]  <- paste("$\\\\sigma$");
        cur.fill[4,14]  <- r.f(simu.par$ysig);
        cur.fill[5,13]  <- paste("$\\\\rho$");
        cur.fill[5,14]  <- r.f(simu.par$corCov);

        cur.fill[1:9,15] <- "\\\\cline{3-14}"
        cur.fill[10,15] <- "\\\\hline";
        fill.all <- rbind(fill.all, cur.fill);
    }
    fill.all <- c(as.vector(t(fill.all)), sce);
    sub.txt.w.no(fill.all, template.f = template.f, out.f = out.f);
}

simu.post.prt <- function(simu.ext, sces = c("I","II","III","IV","V"), m3=c("obs", "ps", "dp", "hdp"),
                          template.f = "./Figures/template_simurst.tex",
                          out.f = "out.tex") {
    r.int <- function(x) {
        round(x)
    }

    r.f <- function(x) {
        formatC(x, format="f", digits = 2);
    }

    fill.all <- NULL;
    for (i in 1:length(simu.ext)) {
        load(paste("./SimuRst/simu_rst_", simu.ext[i], ".Rdata", sep=""));

        n.study                          <- length(SIMU.PAR.LST);
        cur.fill                         <- matrix("", n.study, 15);
        cur.fill[1,1]                    <- sces[i];
        cur.fill[n.study,ncol(cur.fill)] <- "\\\\hline";
        for (s in 1:n.study) {
            cur.fill[s,2:3] <- r.int(c(s, SIMU.PAR.LST[[s]]$nPat));
            cur.fill[s,4]   <- r.f(SIMU.PAR.LST[[s]]$ysig);
            c.row <- NULL;
            for (m in 1:length(m3)) {
                cur.v <- subset(sum.rst, method==m3[m] & study ==s, select=est);
                c.row <- c(c.row, cur.v[1,1] - SIMU.PAR.LST[[s]]$Effect);
            }
            for (m in 1:length(m3)) {
                cur.v <- subset(sum.rst, method==m3[m] & study ==s, select=mse);
                c.row <- c(c.row, cur.v[1,1]);
            }

            if (n.study == s) {
                ##randomized
                c.row <- c(c.row, c.row[7:8]/c.row[5]);
            } else {
                c.row <- c(c.row, c.row[7:8]/c.row[6]);
            }

            cur.fill[s, 5:14] <- r.f(c.row);
        }

        cur.fill[n.study, c(6,7,10,11,13)] <- "";
        cur.fill[n.study, 14] <- paste(cur.fill[n.study, 14], "*", sep="");

        fill.all <- rbind(fill.all, cur.fill);
    }

    fill.all <- as.vector(t(fill.all));
    sub.txt.w.no(fill.all, template.f = template.f, out.f = out.f);
}


##------------------------------------------------------------------------------
##
##                   SOLVD ANALYSIS
##
##------------------------------------------------------------------------------
get.solvd.data <- function(f.sdata   = "./SOLVD/all_trials.csv",
                           f.regdata = "./SOLVD/RBF_REG.csv") {
    ##SDATATMP <- read.csv("./SOLVD/SolvdTv8.csv");
    ##ANADATA$covhosply <- ANADATA$hosplastyr;
    ##ANADATA$covheight <- ANADATA$heightcm;

    ## NICK CODE
    ## form <- ttodthorchfhosp ~ deathorchfhosp + trt + study + age + gend + smoke + diabet + lvef + himi +
    ##     depedema + pulmedema + crackles + weightkg + anydiurbl + avsys + avdia + beat +
    ##     creatinine + sodium + copd + histk
    ## trials_droppedobs <- model.frame(form, all_trials)
    ## table(trials_droppedobs$study)

    SDATAALL <- read.csv(f.sdata);
    REGDATA  <- read.csv(f.regdata);
    SOLVD    <- subset(SDATAALL,
                       study %in% c("SOLVD Prevention", "SOLVD Treatment", "SOLVD Registry"));

    ANADATA  <- sqldf("select a.*, b.RBF4Z1 as chf, b.RBF6 as sex2
                   from SOLVD a
                   left join REGDATA b on (a.id = b.id_reg)")

    ANADATA[which(ANADATA$study == "SOLVD Prevention"), "chf"] <- "N";
    ANADATA[which(ANADATA$study == "SOLVD Treatment"),  "chf"] <- "Y";

    ANADATA$Study <- sapply(as.character(ANADATA$study), function(x) {
        switch(x,
               "SOLVD Treatment"  = 1,
               "SOLVD Prevention" = 2,
               "SOLVD Registry"   = 3,
               NA)
    });

    ANADATA$Randomized <- as.integer(3 != ANADATA$Study);
    ANADATA$Y  <- as.integer(ANADATA$ttodthorchfhosp < 365.25);
    ANADATA$Z  <- sapply(as.character(ANADATA$trtment), function(x) {
        switch(x,
               "Any ACE"      = 1,
               "Enalapril"    = 1,
               "Trandolapril" = 1,
               "No ACE"       = 0,
               "Placebo"      = 0,
               NA
               )
    })

    ANADATA$covchf    <- as.integer("Y" == ANADATA$chf);
    ANADATA$covage    <- ANADATA$age;
    ANADATA$covsex    <- ANADATA$gend;
    ANADATA$covsmoke  <- as.factor(ANADATA$smoke);
    ANADATA$covdiabet <- ANADATA$diabet;
    ANADATA$covlvef   <- ANADATA$lvef;
    ANADATA$covhimi   <- ANADATA$himi;
    ANADATA$covdepe   <- ANADATA$depedema;
    ##ANADATA$covpulm   <- ANADATA$pulmedema;
    ANADATA$covcrac   <- ANADATA$crackles;
    ANADATA$covweight <- ANADATA$weightkg;
    ANADATA$covdiurbl <- ANADATA$anydiurbl;
    ANADATA$covsys    <- ANADATA$avsys;
    ANADATA$covdia    <- ANADATA$avdia;
    ANADATA$covbeat   <- ANADATA$beat;
    ANADATA$covcreat  <- ANADATA$creatinine;
    ANADATA$covsodium <- ANADATA$sodium;
    ANADATA$covcopd   <- ANADATA$copd;
    ANADATA$covhisk   <- ANADATA$histk;
    ANADATA           <- ANADATA[complete.cases(ANADATA[, c("Study","Z","Y", "Randomized")]),];
    ps.cov            <- names(ANADATA)[grep("^cov", names(ANADATA))];
    rst               <- ANADATA[,c("Study", "Y", "Z", "Randomized", ps.cov)];
}

##draw bootstrap sample and conduct multiple imputation to fill in the missing covariates
get.solvd.sample <- function(solvd.data, bs.inx = 1, m = 5, seed = 10000) {

    pMatrix       <- 1 - diag(ncol(solvd.data));
    pMatrix[1:4,] <- 0;
    pMatrix[,1:4] <- 0;

    set.seed(seed);

    rst <- rep(list(NULL), 4);
    ##randomized study ignoring arms
    for (s in 1:2) {
        cat("study ", s, "\n");
        cur.raw <- subset(solvd.data, Study == s);
        if (bs.inx > 1) {
            cur.raw <- cur.raw[sample(1:nrow(cur.raw), replace = TRUE),];
        }

        rst[[s]] <- mice(cur.raw, m=m, method = "cart",
                         predictorMatrix = pMatrix);
    }

    s <- 3;
    for (grp in 0:1) {
        cat("study ", s, " group ", grp, "\n");
        cur.raw <- subset(solvd.data, Study == s & Z == grp);
        if (bs.inx > 1) {
            cur.raw <- cur.raw[sample(1:nrow(cur.raw), replace = TRUE),];
        }

        rst[[3+grp]] <- mice(cur.raw, m=m, method = "cart",
                             predictorMatrix = pMatrix);
    }

    fdir <- paste("./SOLVD/anasolvd_", bs.inx, sep="");
    dir.create(fdir);
    for (j in 1:m) {
        cfdir <- paste(fdir, "/", j, sep = "");
        dir.create(cfdir);

        cur.complete <- NULL;
        for (k in 1:length(rst)) {
            cur.complete <- rbind(cur.complete, complete(rst[[k]], j));
        }

        write.table(cur.complete,
                    file=paste(cfdir, "/datarct.txt", sep = ""),
                    sep=" ", quote=F, row.names=F, col.names=T);
    }
}

## create 1 complete dataset and draw bootstrap sample from the single
## complete dataset
get.solvd.sample.2 <- function(solvd.data, nbs = 1000, seed = 10000) {

    pMatrix       <- 1 - diag(ncol(solvd.data));
    pMatrix[1:4,] <- 0;
    pMatrix[,1:4] <- 0;

    set.seed(seed);

    rst <- rep(list(NULL), 4);
    for (s in 1:2) {
        cur.raw  <- subset(solvd.data, Study == s);
        rst[[s]] <- complete(mice(cur.raw, m=1, method = "cart",
                                  predictorMatrix = pMatrix),
                             1);
    }

    s <- 3;
    for (grp in 0:1) {
        cur.raw      <- subset(solvd.data, Study == s & Z == grp);
        rst[[3+grp]] <- complete(mice(cur.raw, m=1, method = "cart",
                                      predictorMatrix = pMatrix),
                                 1);
    }

    for (bs.inx in 1:nbs) {
        print(bs.inx);
        fdir <- paste("./SOLVD/anasolvd_", bs.inx, sep="");
        dir.create(fdir);

        j <- 1;
        cfdir <- paste(fdir, "/", j, sep = "");
        dir.create(cfdir);

        cur.complete <- NULL;
        for (k in 1:length(rst)) {
            cur.raw <- rst[[k]];
            if (bs.inx > 1) {
                cur.raw <- cur.raw[sample(1:nrow(cur.raw), replace = TRUE),];
            }
            cur.complete <- rbind(cur.complete, cur.raw);
        }

        write.table(cur.complete,
                    file=paste(cfdir, "/datarct.txt", sep = ""),
                    sep=" ", quote=F, row.names=F, col.names=T);
    }
}


## plot propensity score with different sensitivity parameters
plot.solvd.ps <- function(solvd.data, deltas=c(-2,0,2)) {
    reg.data <- subset(solvd.data, 3 == Study & 0 == covchf);
    inx.com  <- complete.cases(reg.data);
    reg.data <- reg.data[inx.com,];
    ps.cov   <- names(reg.data)[grep("^cov", names(reg.data))];
    gfit     <- get.ps(reg.data, ps.cov = ps.cov, grp = "Z")$glm.fit;
    ps.score <- predict(gfit, reg.data, type = "response");

    rst <- NULL;
    for (d in deltas) {
        cur.eps <- get.expt(ps.score, d);
        rst     <- rbind(rst, cbind(d, Z=reg.data$Z, ps=cur.eps));
    }
    rst   <- data.frame(rst, row.names = NULL);
    rst$Z <- c("Placebo", "Enalapril")[rst$Z+1];
    rst$d <- as.factor(rst$d);
    levels(rst$d) <- paste("Delta == ", deltas, sep = "");

    ggplot(rst, aes(x=ps)) +
        geom_density(alpha = 0.2, fill = "gray", aes(group=Z, linetype=Z)) +
        scale_x_continuous(breaks = c(0,0.5,1)) +
        scale_y_continuous(breaks = NULL) +
        labs(x = "Propensity Score", y="Density", title="") +
    ##scale_color_manual(values=c("black", "gray")) +
    theme_bw()+
        theme(strip.background = element_blank(),
              strip.text = element_text(size=12),
              panel.grid = element_blank(),
              legend.title=element_blank(),
              legend.position = c(0.9, 0.88),
              panel.spacing = unit(0.5, "lines"))+
        facet_grid(~d, labeller = label_parsed);
}

##print solvd baseline characteristics
tbl.solvd.baseline <- function(solvd.data) {

    sdata   <- solvd.data;
    sdata$S <- sdata$Study;
    sdata$S[2 == sdata$Study] <- 3;
    sdata$S[3 == sdata$Study & 1 == sdata$covchf] <- 2;
    sdata$S[3 == sdata$Study & 0 == sdata$covchf] <- 4;

    fmt <- function(x, ns=1) {
        format(x, digits = 2, nsmall = ns);
    }

    f.mis <- function(var) {
        rst <- "& Missing";
        for (s in 1:4) {
            for (t in 1:0) {
                cur.d <- subset(sdata, S == s & Z == t);
                rst <- paste(rst, "& ",
                             round(100*mean(is.na(cur.d[[var]])), 0.1))
            }
            if (s < 4)
                rst <- paste(rst, "&");
        }
        rst <- paste(rst, "\\", "hline       ", sep="")
        print(rst);
    }

    f.n <- function() {
        rst <- NULL;
        for (s in 1:4) {
            for (t in 1:0) {
                rst <- c(rst, "& ", nrow(subset(sdata, S == s & Z == t)))
            }
            if (s < 4)
                rst <- c(rst, " & ");
        }
        print(paste(rst, collapse = ""));
    }

    f.cont <- function(var, lbl, ns=0) {
        rst <- paste("{", "it ", lbl, "} & Mean", sep="");
        for (s in 1:4) {
            for (t in 1:0) {
                cur.d <- subset(sdata, S == s & Z == t);
                rst <- paste(rst, "& ",
                             fmt(mean(cur.d[[var]], na.rm = T), ns=ns),
                             "(", fmt(sd(cur.d[[var]], na.rm = T), ns=ns), ")", sep = "")
            }
            if (s <4)
                rst <- paste(rst, "&");
        }
        rst <- paste(rst, "\\     ");
        print(rst);
        f.mis(var);
    }

    f.cat <- function(var, lbl, levels=0:1) {
        for (l in 1:length(levels)) {
            if (1 == l) {
                rst <- paste("{", "it ", lbl, "} &", sep="");
            } else {
                rst <- "&";
            }

            rst <- paste(rst, levels[l], sep = " ");

            for (s in 1:4) {
                for (t in 1:0) {
                    cur.d <- subset(sdata, S == s & Z == t);
                    cur.v <- levels[l] == cur.d[[var]];
                    rst   <- paste(rst, "&", sum(cur.v, na.rm = T), "(",
                                   fmt(100*mean(cur.v,na.rm = T),ns=0), ")", sep="");
                }
                if (s <4)
                    rst <- paste(rst, "&");
            }
            rst <- paste(rst, "\\   ");
            print(rst);
        }

        f.mis(var);
    }

    f.n();
    f.cat("covchf", "CHF");
    f.cont("covage", "Age");
    f.cont("covweight", "Weight");
    f.cat("covsex", "Gender");
    f.cat("covsmoke", "Smoking", c("Current", "Former", "Never"));
    f.cont("covlvef", "LVEF");
    f.cont("covcreat", "Creatinine");
    f.cont("covsodium", "Sodium");
    f.cat("covdepe", "Dependent edema");
    f.cat("covdiabet", "Diabetes");
    f.cat("covhimi", "Myocardial");
    f.cont("covcrac", "Cardiothoracic");
    f.cat("covdiurbl", "Diuretic");
    f.cat("covcopd", "COPD");
    f.cont("covsys", "Blood");
    f.cont("covdia", "Diastolic");
    f.cont("covbeat", "Pulse");
    f.cat("covhisk", "Stroke");
}

##simulation analysis
batch.solvd.ana <- function(bs.inx, m.inx, d.inx = 0, d.start = 0, d.step = 0,  prefix = "./SOLVD/") {
    SOLVD.OPS <- list(n.iter=6000,
                      n.discard=4000,
                      n.batch=20,
                      mcmc.eps=1,
                      eps=0.5,
                      m.prior=1,
                      verbose=3,
                      B.prior=1, S.prior=1, alpha.prior=1,
                      alpha=1,
                      q=10,
                      cc=10,
                      y.type = "binary");

    delta <- d.start + d.step * d.inx;

    ##print info
    print(paste("---CUR BS:--", bs.inx,
                "---CUR M:", m.inx,
                "---CUR DELTA:", delta, sep=""));

    ##set work dir
    old.wd <- getwd();
    setwd(paste(prefix, "anasolvd_", bs.inx, "/", m.inx, sep=""));

    ##------------read data----------------------
    rct.data   <- read.table("datarct.txt",header=T);

    ##load(file = "SOLVD/anasolvd.Rdata");
    ##inx.com      <- complete.cases(solvd.data);
    ##rct.data     <- solvd.data[inx.com,];

    ## get into sub directory of delta
    delta.dir <- paste(getwd(), "/delta_", delta, sep="");
    if (!file.exists(delta.dir)) {
        dir.create(delta.dir);
    }

    old2.wd <- getwd();
    setwd(delta.dir);

    study1       <- subset(rct.data, Study == 1);
    study2       <- subset(rct.data, Study == 2);
    study3       <- subset(rct.data, Study == 3);
    study4       <- subset(rct.data, Study == 3 & 1 == covchf);
    study4$Study <- 4;
    study5       <- subset(rct.data, Study == 3 & 0 == covchf);
    study5$Study <- 5;
    ps.cov       <- names(rct.data)[grep("^cov", names(rct.data))];
    ps.cov.nochf <- ps.cov[ps.cov != "covchf"];

    ## observed
    rst.obs  <- get.obs.effect(rbind(study1, study2, study3, study4, study5),
                               study="Study", group="Z", y="Y");
    ## ps stratification
    ## rst.ps.1  <- get.ps.effect(study3,
    ##                            study="Study",
    ##                            group="Z",
    ##                            ps.cov=ps.cov,
    ##                            y="Y",
    ##                            delta=delta,
    ##                            n.boots=0);
    rst.ps  <- get.ps.effect(rbind(study4, study5),
                             study="Study",
                             group="Z",
                             ps.cov=ps.cov.nochf,
                             y="Y",
                             delta=delta,
                             n.boots=0);

    ## hdp.model
    ## hdp.ps.1 <- get.all.ps(rbind(study1, study4),
    ##                        study="Study",
    ##                        group="Z",
    ##                        random="Randomized",
    ##                        ps.cov=ps.cov.nochf,
    ##                        delta=delta,
    ##                        ps.ignore.random=TRUE,
    ##                        take.logit=TRUE);

    ## pred.pts <- as.matrix(hdp.ps.1[, c("Study", "score4"), drop=FALSE]);
    ## sd.hdp.1 <- get.hdp.all(hdp.ps.1,
    ##                         pred.pts,
    ##                         studies = c(1,4),
    ##                         center.y = FALSE,
    ##                         hdp.cov = "score4",
    ##                         simu.ops= SOLVD.OPS);
    ## rst.hdp.1 <- get.trt.effect(hdp.ps.1,
    ##                             studies = c(1,4),
    ##                             y.type = "binary");

    hdp.ps.2 <- get.all.ps(rbind(study1, study2, study3),
                           study="Study",
                           group="Z",
                           random="Randomized",
                           ps.cov=ps.cov,
                           delta=delta,
                           ps.ignore.random=TRUE,
                           take.logit=TRUE);

    pred.pts  <- as.matrix(hdp.ps.2[, c("Study", "score3"), drop=FALSE]);
    sd.hdp.2  <- get.hdp.all(hdp.ps.2,
                             pred.pts,
                             studies = c(1,2,3),
                             hdp.cov = "score3",
                             center.y = FALSE,
                             simu.ops= SOLVD.OPS);

    rst.hdp.2 <- get.trt.effect(studies = c(1,2,3),
                                y.type = "binary");

    study3.chf   <- which(1 == study3$covchf);
    study3.nochf <- which(0 == study3$covchf);
    rst.hdp.chf  <- get.trt.effect(studies = 3,
                                   sub.inx = study3.chf,
                                   y.type  = "binary");

    rst.hdp.nochf <- get.trt.effect(studies = 3,
                                    sub.inx = study3.nochf,
                                    y.type  = "binary");

    rst.hdp <- rbind(rst.hdp.2, rst.hdp.chf, rst.hdp.nochf);

    ## dp model for chf and non-chf two subgroups
    ## dp.ps <- get.all.ps(rbind(study4, study5),
    ##                     study="Study",
    ##                     group="Z",
    ##                     random="Randomized",
    ##                     ps.cov=ps.cov.nochf,
    ##                     delta=delta,
    ##                     ps.ignore.random=TRUE,
    ##                     take.logit=TRUE);

    ## pred.pts  <- as.matrix(dp.ps[,c("Study", "ps")]);
    ## dp.sd     <- get.dp.all(dp.ps,
    ##                         pred.pts,
    ##                         studies  = 4:5,
    ##                         center.y = FALSE,
    ##                         simu.ops = SOLVD.OPS);
    ## rst.dp   <- get.trt.effect(dp.ps,
    ##                            studies = 4:5,
    ##                            y.type = "binary");
    rst.dp <- NULL;

    ##-----------save results------------------
    setwd(old2.wd);
    save(rst.obs, rst.ps, rst.dp, rst.hdp,
         file=paste("solvdrst_", delta, ".Rdata", sep = ""));

    ##return to previous dir
    setwd(old.wd);
}

solvd.post.combine <- function(cmb.reps = 1:1000, m = 5, d.start = -2, d.step = 0.1, d.inx=0:60) {

    deltas <- d.start + d.step * d.inx;

    ##simu results files
    f.rst <- paste("solvd_rst.Rdata", sep="");
    ##combine
    all.rst <- NULL;
    for (cr in 1:length(cmb.reps)) {
        print(cr);
        r       <- cmb.reps[cr];
        cur.rep <- NULL;
        for (d in deltas) {
            cur.rst <- NULL;
            for (i in 1:m) {
                f.fit <- paste("SOLVD/anasolvd_", r, "/", i, "/solvdrst_", d, ".Rdata", sep = "");
                if (!file.exists(f.fit)) {
                    print(paste(f.fit, " does not exist...", sep=""));
                    cur.m <- NA;
                } else {
                    load(f.fit);
                    cur.m <- rst.obs[, c('dif')];
                    cur.m <- c(cur.m, rst.ps[,"est"]);
                    ##cur.m <- c(cur.m, rst.dp[,1]);
                    cur.m <- c(cur.m, rst.hdp[,1]);
                }
                cur.rst <- cbind(cur.rst, cur.m);
            }
            cur.rep <- c(cur.rep, apply(cur.rst, 1, mean));
        }
        all.rst <- cbind(all.rst, cur.rep);
    }

    all.rst <- cbind(rep(deltas, each = nrow(all.rst)/length(deltas)),  all.rst);

    ##save results
    print(paste(f.rst, " generated...."));
    save(all.rst, file=f.rst);
}

## for generating figure 4 in the paper
solvd.post.comparison <- function(prefix = "./SOLVD/anasolvd_1/1/") {
    rct.data   <- read.table(paste(prefix, "datarct.txt", sep = ""), header=T);
    study3     <- subset(rct.data, Study == 3);
    study3$pid <- 1:nrow(study3);
    ps.cov     <- names(rct.data)[grep("^cov", names(rct.data))];
    ps.cov.nochf <- ps.cov[ps.cov != "covchf"];

    rst      <- NULL;
    old.wd   <- setwd(paste(prefix, "delta_0", sep=""));
    for (chf in 0:1) {
        cur.study <- subset(study3, chf == covchf);
        cur.ps    <- get.ps(cur.study, ps.cov.nochf, grp = "Z")$score;
        cuts      <- quantile(cur.ps, seq(0,1,length=6));
        cuts[1]   <- cuts[1]-0.01;

        for (stra in 1:5) {
            cur.break <- which(cur.ps <= cuts[stra+1] & cur.ps > cuts[stra]);
            cur.study[cur.break,"strata"] <- stra;
        }

        eff.stra <- NULL;
        eff.hdp  <- NULL;
        for (stra in 1:5) {
            cur.grp0 <- subset(cur.study, 0 == Z & stra == strata);
            cur.grp1 <- subset(cur.study, 1 == Z & stra == strata);
            eff.stra <- c(eff.stra, mean(cur.grp1$Y) - mean(cur.grp0$Y));

            cur.ehdp <- get.trt.effect(studies = 3,
                                       sub.inx = c(cur.grp0$pid, cur.grp1$pid),
                                       y.type  = "binary");
            eff.hdp  <- c(eff.hdp, cur.ehdp[1]);
        }

        rst <- rbind(rst, data.frame(CHF    = c("Without CHF", "With CHF")[chf + 1],
                                     Strata = 1:5,
                                     Method = "PS Stratification",
                                     Effect = eff.stra
                                     ));

        rst <- rbind(rst, data.frame(CHF    = c("Without CHF", "With CHF")[chf + 1],
                                     Strata = 1:5,
                                     Method = "EPS-HDP",
                                     Effect = eff.hdp
                                     ));

    }
    setwd(old.wd);
    data.frame(rst)
}

## for testing if there are differences between stratum 5 vs 1-4
solvd.post.diagnosis <- function(prefix = "./SOLVD/anasolvd_1/1/") {
    rct.data     <- read.table(paste(prefix, "datarct.txt", sep = ""), header=T);
    study3       <- subset(rct.data, Study == 3);
    study3$pid   <- 1:nrow(study3);
    ps.cov       <- names(rct.data)[grep("^cov", names(rct.data))];
    ps.cov.nochf <- ps.cov[ps.cov != "covchf"];

    rst <- NULL;
    for (chf in 0:1) {
        cur.study <- subset(study3, chf == covchf);
        cur.ps    <- get.ps(cur.study, ps.cov.nochf, grp = "Z")$score;
        cuts      <- quantile(cur.ps, seq(0,1,length=6));
        cuts[1]   <- cuts[1]-0.01;

        for (stra in 1:5) {
            cur.break <- which(cur.ps <= cuts[stra+1] & cur.ps > cuts[stra]);
            cur.study[cur.break,"sgrp"] <- ifelse(stra == 5, 1, 0);
        }
        rst <- rbind(cur.study, rst);
        ## for (j in 1:length(ps.cov.nochf)) {
        ##     if (ps.cov.nochf[j] %in% c("covage", "covlvef", "covweight", "covsys", "covdia",
        ##                                "covbeat", "covcreat", "covsodium")) {
        ##         fmla     <- formula(paste(ps.cov.nochf[j], " ~ sgrp"));
        ##         cur.test <- t.test(fmla, data = cur.study);
        ##     } else {
        ##         cur.test <- fisher.test(cur.study[[ps.cov.nochf[j]]], cur.study$sgrp);
        ##     }
        ##     rst <- rbind(rst, c(chf, ps.cov.nochf[j], cur.test$p.value));
        ## }
    }

    rst
}
