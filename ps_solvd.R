rm(list=ls());
options(error = recover, warn=2);
source("ps_toolkit.R");

##--------------------------------------------------------------------------------
##                        BATCH
##--------------------------------------------------------------------------------
ARGS.INX <- as.numeric(commandArgs(trailingOnly=TRUE));
if (0 == length(ARGS.INX)) {
    ARGS.INX <- 1;
    print(paste("No argument entered"));
}

N.NODE   <- 40;
M        <- 1;
ANA.BOOT <- 1000;

if (1) {
    load(file = "SOLVD/solvd_ana.Rdata");
    ## n.each <- ANA.BOOT / N.NODE;
    ## reps   <- (ARGS.INX-1)*n.each+(1:n.each);
    ## parallel::mclapply(reps,
    ##                    function(i) {
    ##               get.solvd.sample(solvd.data,
    ##                                bs.inx = i,
    ##                                m      = M,
    ##                                seed   = ARGS.INX*100000+i*10)
    ##           },
    ##           mc.cores = 20);
    get.solvd.sample.2(solvd.data, nbs = ANA.BOOT, seed = 1000000);

} else {
    simu.all  <- expand.grid(1:ANA.BOOT, 1:M, 0);
    n.each    <- ceiling(nrow(simu.all)/N.NODE);
    reps      <- (ARGS.INX-1)*n.each+(1:n.each);
    print(reps);
    parallel::mclapply(reps,
                       function(i) {
                  batch.solvd.ana(bs.inx    = simu.all[i,1],
                                  m.inx     = simu.all[i,2],
                                  d.inx     = simu.all[i,3],
                                  d.start   = 0,
                                  d.step    = 0);
              },
              mc.cores = 20);
}

