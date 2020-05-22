#!/bin/bash

vertice=(10 100 1000 5000 10000 20000)
edge=(100 500 1000 5000 10000 20000 40000)

for V in ${vertice[@]}
do
    for E in ${edge[@]}
    do
        if [[ $E -ge $(($V)) ]] && [[ $E -lt $(($V*10)) ]]
        then
            # to get general measures from perf
            echo "Start test. V = ${V}, E = ${E}" | tee -a "./log/runtime.txt"
            2>>"./log/runtime.txt" sudo perf stat --log-fd 2 --append -d ./build/test/test_solver $V $E | tee -a "./log/runtime.txt"

        fi
    done
done

