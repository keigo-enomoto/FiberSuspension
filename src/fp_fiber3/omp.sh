#!/bin/sh 

# make clean
rm *.o
make -f omp.mk mps_omp
