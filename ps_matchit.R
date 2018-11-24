rm(list=ls());
options(error = recover, warn=0);
source("ps_toolkit.R");
require(parallel);
require(MatchIt);
require(randomForest);
require(hdpmn);

simu.data <- function(nsize, beta.y, beta.z) {
    xx   <- rnorm(nsize, 1, 1);
    bzx  <- beta.z[1] + beta.z[2] * xx;
    bzy  <- beta.y[1] + beta.y[2] * xx + beta.y[3] * xx^2;
    y    <- rnorm(nsize, bzy);
    z    <- rbinom(nsize, 1, exp(bzx)/(1+exp(bzx)));

    glm.fit   <- glm(z~xx, family=binomial);
    ps        <- glm.fit$fitted;

    data.frame(x=xx, ps=ps, z=z, y=y);
}

get.hdp <- function(dta) {
    d <- tempfile(tmpdir=".");
    dir.create(d);

    oldwd   <- setwd(d);
    HDP.OPS <- list(n.iter=3000,
                    n.discard=1000,
                    n.batch=20,
                    mcmc.eps=0,
                    eps=1,
                    m.prior=1,
                    B.prior=1, S.prior=1, alpha.prior=1,
                    alpha=10,
                    q=5,
                    verbose=0,
                    cc=5);

    dta$eps <- log(dta$ps/(1-dta$ps));

    rpy <- NULL;
    for (t in 0:1) {
        cur.d     <- subset(dta, z == t);
        cur.z     <- as.matrix(cur.d[, c("y", "eps")]);
        cur.study <- rep(1, nrow(cur.z));
        do.call(hdpmn, c(list(Z=cur.z, study=cur.study), px=1, HDP.OPS));
        cur.pred <- R.hdpmnPredict(j=1, r=0, nsim=0, idx.x=2,
                                   X= as.matrix(dta[, "eps", drop=FALSE]));
        cur.pred <- read.table("pred_1_0");
        rpy[t+1] <- mean(cur.pred[,3]);
    }

    setwd(oldwd);
    unlink(d, recursive=TRUE);

    rpy[2] - rpy[1];
}

get.rf <- function(dta) {
    ##random forest
    rpy <- NULL;
    rpx <- NULL;
    for (t in 0:1) {
        cur.d <- subset(dta, z == t);
        rfst  <- randomForest(y ~ ps, data = cur.d, samptype="swor", ntree=2000);
        cur.y <- predict(rfst, newdata = dta);
        rpy   <- c(rpy, mean(cur.y));

        ##rfst  <- randomForest(y ~ x, data = cur.d, samptype="swor", ntree=2000);
        ##cur.y <- predict(rfst, newdata = dta);
        ##rpx   <- c(rpx, mean(cur.y));

    }
    est.rf  <- rpy[2] - rpy[1];
    ##est.rfx <- rpx[2] - rpx[1];

    c(est.rf);
}

get.match <- function(dta) {
    ##matching
    mtd       <- matchit(z ~ x, data = dta, method="nearest");
    mtd.dta   <- match.data(mtd);
    inx.0     <- which(0 == mtd.dta[,"z"]);
    est.mt    <- mean(dta[-inx.0, "y"]) - mean(dta[inx.0, "y"]);

    est.mt
}

get.est <- function(dta) {

    ##hdp
    est.hdp <- get.hdp(dta);

    ##naive
    inx.0     <- which(0 == dta[,"z"]);
    est.naive <- mean(dta[-inx.0, "y"]) - mean(dta[inx.0, "y"]);


    c(est.naive, get.match(dta), get.rf(dta), est.hdp);
}

simu.all <- function(nreps, ..., n.cores=3) {

    rst <- parallel::mclapply(1:nreps,
                              function(x) {
                         print(x);
                         cur.dta <- simu.data(...);
                         rst     <- tryCatch({
                             get.est(cur.dta);
                         }, error = function(e) {
                             print(e);
                         })
                     }, mc.cores=n.cores);

    rst      <- matrix(unlist(rst),
                       ncol=length(rst[[1]]),
                       byrow = TRUE);

    rst.mean <- apply(rst, 2, mean);
    rst.mse  <- apply(rst, 2, function(x){mean(x^2)});

    rbind(rst.mean, rst.mse);
}

NREP     <- 1000;
rst      <- list();
rst[[1]] <- simu.all(NREP, 50,  c(0,10, 0), c(0, 0));
rst[[2]] <- simu.all(NREP, 200, c(0,10, 0), c(0, 0));
rst[[3]] <- simu.all(NREP, 50,  c(0,10,0), c(-5, 5));
rst[[4]] <- simu.all(NREP, 200, c(0,10,0), c(-5, 5));

##save(rst, NREP, file="rst_toyexp.Rdata");
