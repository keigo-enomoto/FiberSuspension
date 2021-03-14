#!/bin/sh

# don't use OpenMP
#./../../src/fp_fiber3/mps input_cluster.txt

#use OpenMP
export OMP_NUM_THREADS=2
./../../src/fp_fiber3/mps_omp input_cluster.txt

