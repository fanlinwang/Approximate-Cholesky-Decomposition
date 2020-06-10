#!/bin/bash

# fix deg = 10
vertice=(100 1000 5000 10000 50000)
edge=(1000 10000 50000 100000 500000)

for V in ${vertice[@]}
do
    for E in ${edge[@]}
    do
        if [[ $E -eq $(($V*10)) ]]
        then
            # to get exact cycles, run main without perf. 
            valgrind --tool=cachegrind --cache-sim=yes --branch-sim=yes ./build/test/runtime $V $E 1 | tee -a "./log/cache.txt"
         fi
    done
done

# cg_annotate --auto=yes --include=/home/fanlin/team008 cachegrind.out.xxx ./src/approxChol.cpp
# valgrind --tool=cachegrind --cache-sim=yes --branch-sim=yes ./build/test/runtime 50000 500000 1
