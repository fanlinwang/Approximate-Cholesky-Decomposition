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
            echo "Start test. V = ${V}, E = ${E}" | tee -a "./log/general_metrics.txt"
            3>>"./log/general_metrics.txt" perf stat --log-fd 3 --append -d ./build/test/test_solver $V $E | tee -a "./log/general_metrics.txt"

            # to get cache loads from perf
            echo "Start test. V = ${V}, E = ${E}" | tee -a "./log/cache_log.txt"
            3>>"./log/cache_log.txt" perf stat -e LLC-loads,LLC-load-misses,LLC-stores,LLC-store-misses,dTLB-loads,dTLB-load-misses,dTLB-stores,dTLB-store-misses,iTLB-loads,iTLB-load-misses --log-fd 3 --append -d ./build/test/test_solver $V $E | tee -a "./log/cache_log.txt"

            # to get flops counts from perf
            echo "Start test. V = ${V}, E = ${E}" | tee -a "./log/flops_log.txt"
            3>>"./log/flops_log.txt" perf stat -e r5301c7,r5302c7 --log-fd 3 --append -d ./build/test/test_solver $V $E | tee -a "./log/flops_log.txt"

            if [[ $1 == 'p' ]] 
            then
                perf record -o "./log/bin/perf_v${V}_e${E}.data" ./build/src/main $V $E
                perf report --stdio --dsos=main -i "./log/bin/perf_v${V}_e${E}.data" | tee -a "./log/report_main/pf_v${V}_e${E}.txt"
                perf annotate --stdio --dsos=main -i "./log/bin/perf_v${V}_e${E}.data"  | tee -a "./log/annotate_main/p_v${V}_e${E}.txt"
                1>>"./log/annotate_approxChol/perf_annotate_v${V}_e${E}.txt" perf annotate --stdio --dsos=main --symbol=approxChol -i "./log/bin/perf_v${V}_e${E}.data"  
            fi
        fi
    done
done

