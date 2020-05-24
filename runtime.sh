#!/bin/bash

vertice=(100 1000 5000 8000 10000)
edge=(200 500 1000 2000 5000 10000 16000 20000 25000 40000 50000 80000 100000 160000 500000)

for V in ${vertice[@]}
do
    for E in ${edge[@]}
    do
        if [[ $E -gt $(($V)) ]] && [[ $E -lt $(($V*50)) ]]
        then
            # to get exact cycles, run main without perf. 
            ./build/test/runtime $V $E 2 | tee -a "./log/runtime.txt"
         fi
    done
done

