#!/bin/bash

#seq FIRST STEP LAST
vertice=(10 100 1000)
edge=(10 100 1000 10000)

for V in ${vertice[@]}
do
    for E in ${edge[@]}
    do
        if [ $E -lt $(($V*$V)) ]
        then
            echo "Start test. V = ${V}, E = ${E}" | tee -a "./cache_log.txt"
            3>>cache_log.txt perf stat -e L1-dcache-loads,L1-dcache-load-misses,L1-dcache-stores,L1-dcache-store-misses,L1-dcache-prefetches,L1-dcache-prefetch-misses,LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses,LLC-prefetch-misses,dTLB-loads,dTLB-load-misses,dTLB-stores,dTLB-store-misses,dTLB-prefetches,dTLB-prefetch-misses,iTLB-loads,iTLB-load-misses --log-fd 3 --append -d ./build/src/main $V $E | tee -a "./cache_log.txt"
            echo "Start test. V = ${V}, E = ${E}" | tee -a "./flops_log.txt"
            3>>flops_log.txt perf stat -e r5301c7,r5302c7 --log-fd 3 --append -d ./build/src/main $V $E | tee -a "./flops_log.txt"
            perf record -o "./perf_v${V}_e${E}.data"
            3>>"./perf_v${V}_e${E}.data" perf report -i "./perf_v${V}_e${E}.data"
        fi
    done
done

