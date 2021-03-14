#!/bin/sh

# change here ==========
xx_arr=( 1 )
yy_arr=( 60 )
ap_arr=( 10 )
nf_arr=( 42 )
ss_arr=( 98 )

# ========================

for x_mult in "${xx_arr[@]}"; do
for yy in "${yy_arr[@]}"; do
for ap in "${ap_arr[@]}"; do
for nf in "${nf_arr[@]}"; do
for ss in "${ss_arr[@]}"; do

    xx=`expr ${yy} \* ${x_mult}`
    sed -e "s/xx/${xx}/g ; s/yy/${yy}/g ; s/ap/${ap}/g; s/nf/${nf}/g; s/ss/${ss}/g; " < input_tmp.txt > input.txt
    ./../../src/mk_fiber/2d input.txt > mk_log.txt
    dat=`ls *.dat`
 
    sed -i -e "s/@dat@/${dat}/g" input_cluster.txt

    ./../../src/fp_fiber3/init input_cluster.txt
    rm input.txt

done
done
done
done
done
