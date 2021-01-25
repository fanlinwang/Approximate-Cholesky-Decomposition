#!/bin/bash

vertice=(100 1000 5000 10000 50000 100000)
edge=(200 500 1000 2000 5000 10000 20000 25000 40000 50000 80000 100000 500000 1000000)

for V in ${vertice[@]}
do
    for E in ${edge[@]}
    do
        if [[ $E -gt $(($V)) ]] && [[ $E -lt $(($V*50)) ]]
        then
            # to get exact cycles, run main without perf. 
            ./build/src/main $V $E | tee -a "./log/exact_cycle.txt"
         fi
    done
done

./build/src/main 100000 1000000 | tee -a "./log/exact_cycle.txt"
