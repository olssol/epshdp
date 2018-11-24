rm(list=ls());
options(error = recover, warn=2);
source("ps_toolkit.R");

ARGS.INX <- as.numeric(commandArgs(trailingOnly=TRUE));
if (0 == length(ARGS.INX)) {
    ARGS.INX <- 1;
    print(paste("No argument entered"));
}

SIMU.BATCH   <- expand.grid(1:N.REPS, SIMU.EXT.ALL);
N.NODES      <- 500; ##number of computation nodes
N.EACH       <- ceiling(nrow(SIMU.BATCH) / N.NODES);

##----------------------------------------------------------
##                 SIMULATION
##----------------------------------------------------------
##----------------------------------------------------------
##                 SIMULATION
##----------------------------------------------------------
TODO <- c("SIMUPARONLY", "SIMUDATA", "ANALYSIS", "COMBINE")[4];

## current reps
reps <- ((ARGS.INX-1)*N.EACH+1):min(ARGS.INX*N.EACH, nrow(SIMU.BATCH));
print(reps);

if ("SIMUPARONLY" == TODO) {
    get.simu.par();
} else if ("SIMUDATA" == TODO) {
    for (i in reps) {
      get.simu.all(simu.ext.all = SIMU.BATCH[i,2],
                   n.reps       = SIMU.BATCH[i,1],
                   seed         = i * 10);
    }
} else if ("ANALYSIS" == TODO) {
    ## parallel::mclapply(reps,
    ##                    function(i) {
    ##               print(i);
    ##               batch.simu.ana(rep      = SIMU.BATCH[i,1],
    ##                              simu.ext = SIMU.BATCH[i,2]);
    ##           },
    ##           mc.cores = min(N.EACH,parallel::detectCores())
    ##           );
    for (i in reps) {
        batch.simu.ana(rep      = SIMU.BATCH[i,1],
                       simu.ext = SIMU.BATCH[i,2]);
    }
} else if ("COMBINE" == TODO) {
    ## parallel::mclapply(c(671:675, 691:695),
    ##                    function(i) {
    ##               simu.post.combine(i, 1:2000);
    ##           },
    ##           mc.cores = 10
    ##           );
    for (i in c(781:785)) {
        simu.post.combine(i, 1:2000);
    }
}

