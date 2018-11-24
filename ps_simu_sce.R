##------------------------------------------------------------------------------
##
##              SIMULATION SCENARIOS
##
##------------------------------------------------------------------------------
simu.formula.10 <- function() {
    formula.z <- as.formula(paste(c("~-1", paste("V", 1:6, sep="")),
                                  collapse = " + "));

    formula.y <- as.formula(paste(c("~-1",
                                    paste("V", 3:10, sep=""),
                                    paste("I(V", 3:10, "^2)", sep="")),
                                  collapse = " + "));

    c(fz = formula.z,
      fy = formula.y)
}

set.simu.par.base <- function(sizes      = c(1000,800,400,200),
                              mus        = c(2,2,2,2),
                              sds        = c(0.5,0.5,0.5,0.5),
                              zbeta      = list(c(rep(4,3), rep(-4,3)),
                                                c(rep(4,3), rep(-4,3)),
                                                c(rep(4,3), rep(-4,3)),
                                                0),
                              corCov     = 0.1,
                              nu         = 1,
                              nx         = 10,
                              trt.effect = 3,
                              regcoeff   = c(rep(0.2,4), rep(-0.2, 4),  rep(-0.5,4), rep(0.5, 4)),
                              fmls       = simu.formula.10(),
                              prop.z1    = 0.5,
                              ps.cov     = NULL,
                              simu.npat  = 100000,
                              simu.seed  = 100000) {

    studies <- list();
    for (i in 1:length(sizes)) {
        set.seed(simu.seed);
        studies[[i]] <- list(Effect     = trt.effect,
                             nPat       = sizes[i],
                             regCoeff.y = regcoeff,
                             fmla.y     = fmls$fy,
                             fmla.z     = fmls$fz,
                             sig2Ratio  = nu,
                             regCoeff.z = zbeta[[i]],
                             muCov      = rep(mus[i], nx),
                             sdCov      = rep(sds[i], nx),
                             corCov     = corCov);
        cur.par      <- studies[[i]];
        cur.par$nPat <- simu.npat;
        b0.ysig      <- do.call(rweSimuTwoArm,
                                c(cur.par, do.simu=FALSE, bin.mu = prop.z1))$b0ysig;
        studies[[i]] <- c(studies[[i]], b0.ysig["b0"], b0.ysig["ysig"]);
    }

    SIMU.PAR.LST <- studies;

    if (is.null(ps.cov))
        ps.cov <-  1:nx;
    SIMU.PS.COV  <- paste("V", ps.cov, sep="");
    SIMU.HDP.COV <- paste("score", seq_along(studies), sep="");

    ##return
    rst <- make.list(environment());
    rst
}


set.simu.par.671 <- function() {
    set.simu.par.base();
}

set.simu.par.672 <- function() {
    rst <- set.simu.par.base(sds = c(2,1,1,0.5));
}

set.simu.par.673 <- function() {
    rst <- set.simu.par.base(mus = c(-1,1,1.5,2));
}

set.simu.par.674 <- function() {
    rst <- set.simu.par.base(zbeta = list(c(rep(4,3), rep(-4,3)),
                                         c(rep(-4,3), rep(4,3)),
                                         c(4,-4,4,-4,4,-4),
                                         0));
}

set.simu.par.675 <- function() {
    rst <- set.simu.par.base(mus   = c(-1,1,1.5,2),
                            sds   = c(2,1,1,0.5),
                            zbeta = list(c(rep(4,3), rep(-4,3)),
                                         c(rep(-4,3), rep(4,3)),
                                         c(4,-4,4,-4,4,-4),
                                         0));
}

set.simu.par.691 <- function() {
    set.simu.par.base(ps.cov=1:6);
}

set.simu.par.692 <- function() {
    rst <- set.simu.par.base(sds = c(2,1,1,0.5),
                             ps.cov=1:6);
}

set.simu.par.693 <- function() {
    rst <- set.simu.par.base(mus = c(-1,1,1.5,2),
                             ps.cov=1:6);
}

set.simu.par.694 <- function() {
    rst <- set.simu.par.base(zbeta = list(c(rep(4,3), rep(-4,3)),
                                         c(rep(-4,3), rep(4,3)),
                                         c(4,-4,4,-4,4,-4),
                                         0),
                             ps.cov=1:6);
}

set.simu.par.695 <- function() {
    rst <- set.simu.par.base(mus   = c(-1,1,1.5,2),
                            sds   = c(2,1,1,0.5),
                            zbeta = list(c(rep(4,3), rep(-4,3)),
                                         c(rep(-4,3), rep(4,3)),
                                         c(4,-4,4,-4,4,-4),
                                         0),
                            ps.cov=1:6);
}


## 3 studies 771-775 791-795
set.simu.par.771 <- function(...) {
    rst <- set.simu.par.base(sizes      = c(1000,1000,200),
                             mus        = c(2,2,2),
                             sds        = c(0.5,0.5,0.5),
                             zbeta      = list(c(rep(4,3), rep(-4,3)),
                                               c(rep(4,3), rep(-4,3)),
                                               0),
                             ...);
}

set.simu.par.772 <- function(...) {
    rst <- set.simu.par.base(sizes      = c(1000,1000,200),
                             mus        = c(2,2,2),
                             sds        = c(1,0.75,0.5),
                             zbeta      = list(c(rep(4,3), rep(-4,3)),
                                               c(rep(4,3), rep(-4,3)),
                                               0),
                             ...);
}

set.simu.par.773 <- function(...) {
    rst <- set.simu.par.base(sizes      = c(1000,1000,200),
                             mus        = c(1.75,2,2.25),
                             sds        = c(0.5,0.5,0.5),
                             zbeta      = list(c(rep(4,3), rep(-4,3)),
                                               c(rep(4,3), rep(-4,3)),
                                               0),
                             ...);
}

set.simu.par.774 <- function(...) {
    rst <- set.simu.par.base(sizes      = c(1000,1000,200),
                             mus        = c(2,2,2),
                             sds        = c(0.5,0.5,0.5),
                             zbeta      = list(c(rep(4,3), rep(-4,3)),
                                               c(4,-4,4,-4,4,-4),
                                               0),
                             ...);
}

set.simu.par.775 <- function(...) {
    rst <- set.simu.par.base(sizes      = c(1000,1000,200),
                             mus        = c(1.75,2,2.25),
                             sds        = c(1,0.75,0.5),
                             zbeta      = list(c(rep(4,3), rep(-4,3)),
                                               c(4,-4,4,-4,4,-4),
                                               0),
                             ...);
}

set.simu.par.776 <- function(...) {
    rst <- set.simu.par.base(sizes      = c(1000,1000,200),
                             mus        = c(2,2,2),
                             sds        = c(1,1,0.5),
                             zbeta      = list(c(rep(4,3), rep(-4,3)),
                                               c(rep(4,3), rep(-4,3)),
                                               0),
                             ...);
}

## ps model is mis-specified 
set.simu.par.781 <- function() {
  set.simu.par.771(ps.cov = 1:5);
}

set.simu.par.782 <- function() {
  set.simu.par.772(ps.cov = 1:5);
}

set.simu.par.783 <- function() {
  set.simu.par.773(ps.cov = 1:5);
}

set.simu.par.784 <- function() {
  set.simu.par.774(ps.cov = 1:5);
}

set.simu.par.785 <- function() {
  set.simu.par.775(ps.cov = 1:5);

}

## v7-v10 are not included but ps model is correctly specified 
set.simu.par.791 <- function() {
    set.simu.par.771(ps.cov = 1:6);
}

set.simu.par.792 <- function() {
    set.simu.par.772(ps.cov = 1:6);
}

set.simu.par.793 <- function() {
    set.simu.par.773(ps.cov = 1:6);
}

set.simu.par.794 <- function() {
    set.simu.par.774(ps.cov = 1:6);
}

set.simu.par.795 <- function() {
    set.simu.par.775(ps.cov = 1:6);

}

set.simu.par.796 <- function() {
    set.simu.par.776(ps.cov = 1:6);

}

