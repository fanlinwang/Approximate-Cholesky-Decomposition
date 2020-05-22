#!/bin/bash

vertice=(100 1000 5000 10000 20000 50000)
edge=(500 1000 5000 10000 20000 50000 100000 500000)

for V in ${vertice[@]}
do
    for E in ${edge[@]}
    do
        if [[ $E -gt $(($V)) ]] && [[ $E -lt $(($V*50)) ]]
        then
            # to get general measures from perf
            echo "$V $E" | tee -a "./log/julia_runtime.txt"
            julia ./julia/runtime.jl $V $E | tee -a "./log/julia_runtime.txt"

        fi
    done
done

