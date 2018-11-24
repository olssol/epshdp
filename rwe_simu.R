#' Simulate covariates following a mixture of multivariate normal distribution
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuCov <- function(nPat, muCov, sdCov, corCov, mix.phi = 1, seed = NULL, cov.breaks = NULL) {

    f.cur <- function(x, i) {
        if (is.array(x)) {
            rst <- x[min(i, nrow(x)),];
        } else {
            rst <- x;
        }

        rst
    }

    stopifnot(is.numeric(mix.phi) | any(mix.phi < 0));
    stopifnot(nPat > 0);

    if (!is.null(seed))
        set.seed(seed);

    n.pts <- rmultinom(1, nPat, mix.phi);
    cov.x <- NULL;
    for (i in 1:length(mix.phi)) {

        if (0 == n.pts[i])
            next;

        cur.mu  <- f.cur(muCov, i);
        cur.sd  <- f.cur(sdCov, i);
        cur.cor <- corCov[min(i, length(corCov))];
        cur.x   <- rmvnorm(n.pts[i],
                           mean=cur.mu,
                           sigma=get.covmat(cur.sd, cur.cor));
        cov.x   <- rbind(cov.x, cur.x);
    }

    colnames(cov.x) <- paste("V", 1:ncol(cov.x), sep="");
    cov.x           <- data.frame(cov.x);
    cov.x           <- get.cov.cat(cov.x, cov.breaks);
}

#' Simulate X multiplied by Beta
#'
#' @inheritParams simupara
#' @param ... Parameters for simulating covariates by function
#'     \code{\link{rweSimuCov}}
#'
#'
#' @export
#'
rweXBeta <- function(..., regCoeff, cov.x = NULL, fmla = NULL) {

    stopifnot(inherits(fmla,"formula") | is.null(fmla));

    if (is.null(cov.x))
        cov.x <- rweSimuCov(...);

    if (is.null(fmla)) {
        fmla <- formula(paste("~",
                              paste(colnames(cov.x), collapse = "+")));
    }

    d.matrix <- model.matrix(fmla, cov.x);
    xbeta    <- get.xbeta(d.matrix, regCoeff);
    xbeta
}


#' Compute standard error of the random error
#'
#' @inheritParams simupara
#' @inheritParams rweGetYSig
#'
#' @return mean of xbeta and standard error or the random error term
#'
#' @export
#'
rweGetYSig <- function(..., regCoeff, nPat=500000, xbeta = NULL, sig2Ratio = 1) {
    if (is.null(xbeta))
        xbeta   <- rweXBeta(nPat=nPat, regCoeff = regCoeff, ...);

    v.xbeta <- var(xbeta);
    ysig    <- sqrt(v.xbeta * sig2Ratio);

    c(mean(xbeta), ysig);
}

#' Get intercept for a binary outcome.
#'
#' The binary outcome may be an outcome or a treatment assignment.
#'
#' @inheritParams simupara
#'
#' @param ... Parameters for simulating covariates by function
#'     \code{\link{rweXBeta}} and \code{\link{rweSimuCov}}
#'
#' @return standard error or the random error term
#'
#' @export
#'
rweGetBinInt <- function(..., regCoeff, nPat=500000, xbeta = NULL, bin.mu = 0.5) {
    ## fill in 0 for intercept temporarily
    if (is.null(xbeta))
        ey <- rweXBeta(nPat, regCoeff = regCoeff, ...);

    fx <- function(b0) {
        logp <- (b0 + ey) - log(1 + exp(b0+ey));
        m    <- mean(exp(logp));
        abs(m - bin.mu);
    }

    mey <- max(abs(ey));
    rst <- optimize(fx, c(-100-max(ey),100+max(ey)));
    rst$minimum
}

#' Simulate random errors
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuError <- function(nPat,
                         error.type = c("normal", "skewed"),
                         ysig   = 1,
                         skew.n = NULL,skew.p = NULL, skew.noise = 0.0001,
                         ...) {

    type <- match.arg(error.type);
    rst <- switch(type,
                  normal = {rnorm(nPat, 0, ysig)},
                  skewed = {
        mu        <- skew.n * (1-skew.p) / skew.p;
        va        <- skew.n * (1-skew.p) / skew.p^2;
        noise.sig <- skew.noise;
        rst       <- rnbinom(nPat, skew.n, skew.p);
        rst       <- rst - mu + rnorm(nPat, 0, noise.sig);
        rst       <- rst/sqrt(va + noise.sig^2)*ysig;
    });

    rst
}

#' Simulate continuous outcomes for a two arm study
#'
#' @inheritParams simupara
#'
#' @export
#'
rweSimuTwoArm <- function(nPat, muCov, sdCov, corCov, trt.effect = 0,
                          regCoeff.y, regCoeff.z=0,
                          mix.phi = 1, cov.breaks = NULL,
                          fmla.y = NULL, fmla.z = NULL, ysig = NULL,
                          b0 = NULL, z1.p = 0.5,
                          sig2Ratio = 2,  ..., do.simu=TRUE) {
    ##treatment assignment
    if (is.null(b0)) {
        b0  <- rweGetBinInt(bin.mu   = z1.p,
                            muCov    = muCov,
                            sdCov    = sdCov,
                            corCov   = corCov,
                            regCoeff = regCoeff.z,
                            mix.phi  = mix.phi,
                            fmla     = fmla.z);
    }

    if (is.null(ysig)) {
        ysig <- rweGetYSig(muCov      = muCov,
                           sdCov      = sdCov,
                           corCov     = corCov,
                           mix.phi    = mix.phi,
                           cov.breaks = cov.breaks,
                           regCoeff   = regCoeff.y,
                           sig2Ratio  = sig2Ratio,
                           fmla       = fmla.y)[2];
    }

    simu.data <- NULL;
    if (do.simu) {
        ##covariates
        COV.X   <- rweSimuCov(nPat       = nPat,
                              muCov      = muCov,
                              sdCov      = sdCov,
                              corCov     = corCov,
                              mix.phi    = mix.phi,
                              cov.breaks = cov.breaks);

        if (identical(0, regCoeff.z)) {
            Z  <- rbinom(nPat, 1, z1.p);
        } else {
            xbeta.z <- rweXBeta(cov.x    = COV.X,
                                regCoeff = regCoeff.z,
                                fmla     = fmla.z);
            Z       <- rbinom(nPat, 1, expit(b0 + xbeta.z));
        }

        xbeta.y <- rweXBeta(cov.x    = COV.X,
                            regCoeff = regCoeff.y,
                            fmla     = fmla.y);

        epsilon.y <- rweSimuError(nPat, ysig = ysig, ...);
        Y         <- Z * trt.effect + xbeta.y + epsilon.y;
        simu.data <- cbind(pid=1:nPat, Y=Y, Z=Z, COV.X);
    }

    list(true.effect = trt.effect,
         simu.data   = simu.data,
         b0ysig      = c(b0 = b0, ysig = ysig));
}



#' Combining simulation results
#'
#' @param lst.rst List of simulation results. Each element represents a
#'     replication
#' @export
#'
rweSimuCombine <- function(lst.rst, fun = mean, ignore.error = TRUE, ...) {

    if (ignore.error) {
        err.inx <- NULL;
        for (i in 1:length(lst.rst)) {
            if ("try-error" == class(lst.rst[[i]]))
                err.inx <- c(err.inx, i);
        }

        if (!is.null(err.inx))
            lst.rst <- lst.rst[-err.inx];
    }

    nreps       <- length(lst.rst);
    rep1        <- lst.rst[[1]];
    lst.combine <- rep(list(NULL), length(rep1));
    for (i in 1:nreps) {
        for (j in 1:length(rep1)) {
            cur.value        <- lst.rst[[i]][[j]];
            lst.combine[[j]] <- rbind(lst.combine[[j]],
                                      as.vector(as.matrix(cur.value)));
        }
    }

    for (j in 1:length(lst.combine)) {
        cur.rst           <- apply(lst.combine[[j]], 2, fun, ...);
        dim(cur.rst)      <- dim(as.matrix(rep1[[j]]));
        dimnames(cur.rst) <- dimnames(as.matrix(rep1[[j]]));
        lst.combine[[j]]  <- cur.rst;
    }

    names(lst.combine) <- names(rep1);
    lst.combine
}
