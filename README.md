# Source codes for "A Bayesian Non-Parametric Causal Inference Model for Synthesizing Randomized Clinical Trial and Real-World Evidence"

## ps_simu.R

The main R file for the simulation study. The code can run on parallel computing
systems. The following is an example:

```
#!/bin/bash

#PBS -N job
#PBS -q checkpt
#PBS -l nodes=1:ppn=20
#PBS -l walltime=32:00:00
#PBS -j oe

cd $PBS_O_WORKDIR;
R --vanilla --args $jid < ps_simu.R > ${PBS_JOBNAME}_$jid.out 
wait
```

## ps_simu_sce.R

Simulation scenarios. Scenarios 781-785 correspond to Scenarios I-V in the manuscript.

## ps_toolkit.R

Toolkit functions used by simulation and analysis. 

## ps_solvd.R

The code for analyzing the SOLVD data. SOLVD data is available from
https://biolincc.nhlbi.nih.gov/studies/solvd/

## R package: hdpmn_3.0.tar.gz

The R package that implement Hierarchical Dirichlet Process models. The package
was originally developed by Peter Mueller <pm@wotan.mdacc.tmc.edu>. Version 3.0
contains several modifications for this manuscript.
