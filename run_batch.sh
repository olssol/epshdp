#!/bin/bash

#R --vanilla --args 0 < ps_simu.R > log_0.txt

for i in `seq 1 2`; do
    R --vanilla --args $i < ps_simu.R > log_$i.txt &
done

for job in `jobs -p`; do
    wait $job
done

R --vanilla < ps_post.R > log_post.txt

##latexit temp
