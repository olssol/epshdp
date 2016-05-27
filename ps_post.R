rm(list=ls());
options(error = recover, warn=2);
source("ps_toolkit.R");


##combine data
SIMU.EXT <- 18;
sum.rst <- simu.post.combine(SIMU.EXT);




##sips -s format pdf simu_rst_338_pt_2.png --out s338p2.pdf
ARGS.INX <- as.numeric(commandArgs(trailingOnly=TRUE));
if (0 == length(ARGS.INX)) {
    ARGS.INX <- 1;
    print(paste("No argument entered"));
}

for (i in 1:length(SIMU.EXT.ALL)) {
    SIMU.EXT <- SIMU.EXT.ALL[i];
    ##simu.post.combine(SIMU.EXT);

    load(file=get.f.name(SIMU.EXT, lsts=c("simu_par")));
    load(file=get.f.name(SIMU.EXT, 1, c("result")));
    rst.true <- get.true(SIMU.PAR.LST);

    print(paste("treatment:----model ", SIMU.EXT, "---------"));
    pdf(file=get.f.name(SIMU.EXT, 1, c("result"), ".pdf"), width=16, height=9);
    par(mfrow=c(1,2), mar=c(2,2,2,1), cex=1.5);
    for (z in 1:2) {
        plot.all(rst.true[[z]],
                 rst.hdp.sharing[[z]],
                 ##rst.hdp.no.sharing[[z]],
                 rst.reg[[z]],
                 main=paste("Y(",z-1,")",sep=""),
                 ##legends=c("True", "Register+Random", "Register", "Regression"));
                 legends=c("True", "HDP", "Regression"));
    }
    dev.off();

    print(paste("Difference:----model ", SIMU.EXT, "---------"));
    pdf(file=get.f.name(SIMU.EXT, 1, c("result_diff"), ".pdf"), width=16, height=9);
    par(cex=1.5);
    plot.all(get.diff(rst.true),
             get.diff(rst.hdp.sharing),
             ##get.diff(rst.hdp.no.sharing),
             get.diff(rst.reg),
             main="Y(1)-Y(0)",
             legends=c("True", "HDP", "Regression"));
    dev.off();

    print(paste("KL Distance:----model ", SIMU.EXT, "---------"));
    print(KL.dist(get.diff(rst.true), get.diff(rst.hdp.sharing), k=200));
}



if (0) {
    print(paste("Sensitivity:----model ", SIMU.EXT, "---------"));
    rst.d <- NULL;
    for (i in 1:length(DELTA.ALL)) {
        load(get.f.name(SIMU.EXT, 1, lsts=c("result", DELTA.ALL[i])));
        rst.d[[i]] <- rst.hdp.sharing;
    }
    pdf(file=get.f.name(SIMU.EXT, 1, c("result_delta"), ".pdf"), width=16, height=9);
    par(cex=1.5);
    plot.all(get.diff(rst.d[[1]]),
             get.diff(rst.d[[2]]),
             get.diff(rst.hdp.sharing),
             get.diff(rst.d[[3]]),
             get.diff(rst.d[[4]]),
             main="Y(1)-Y(0)",
             legends=c("delta=-2", "delta=-0.5",
             "delta=0", "delta=0.5", "delta=2"));
    dev.off();
}
