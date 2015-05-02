#!/bin/sh
for i in 0 1 2 3 4 5
do
qsub -v runnum=$i,savefn=$1 pbs-numtest.sh
done
