#!/bin/bash

vertice=(10 100 1000 5000 10000 20000)
edge=(100 500 1000 5000 10000 20000 40000)

for V in ${vertice[@]}
do
    for E in ${edge[@]}
    do
        if [[ $E -ge $(($V)) ]] && [[ $E -lt $(($V*10)) ]]
        then
            # to get exact cycles, run main without perf. 
            ./build/src/main $V $E | tee -a "./log/exact_cycle.txt"

         fi
    done
done

