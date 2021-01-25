#!/bin/bash

vertice=(1000 50000)
edge=(5000 100000 250000)

for V in ${vertice[@]}
do
    for E in ${edge[@]}
    do
        if [[ $E -ge $(($V)) ]] && [[ $E -lt $(($V*50)) ]]
        then
            if [[ $1 == 'p' ]]
            then
                perf record -o "./log/bin/perf_v${V}_e${E}.data" ./build/test/runtime $V $E 1
                perf report --stdio --dsos=runtime -i "./log/bin/perf_v${V}_e${E}.data" | tee -a "./log/report_main/pf_v${V}_e${E}.txt"
                perf annotate --stdio --dsos=runtime -i "./log/bin/perf_v${V}_e${E}.data"  | tee -a "./log/annotate_main/p_v${V}_e${E}.txt"
                1>>"./log/annotate/baseline_v${V}_e${E}.txt" perf annotate --stdio --dsos=runtime --symbol=approxChol -i "./log/bin/perf_v${V}_e${E}.data"
                1>>"./log/annotate/VecMg_v${V}_e${E}.txt" perf annotate --stdio --dsos=runtime --symbol=approxChol_vector2_merge -i "./log/bin/perf_v${V}_e${E}.data"
                1>>"./log/annotate/VecMgSearch_v${V}_e${E}.txt" perf annotate --stdio --dsos=runtime --symbol=approxChol_vector2_merge_search -i "./log/bin/perf_v${V}_e${E}.data"
                1>>"./log/annotate/VecMgSearchSIMD_v${V}_e${E}.txt" perf annotate --stdio --dsos=runtime --symbol=approxChol_vector2_merge_search_opt2 -i "./log/bin/perf_v${V}_e${E}.data"
                1>>"./log/annotate/VecStructMgSearch_v${V}_e${E}.txt" perf annotate --stdio --dsos=runtime --symbol=approxChol_vector2_struct_merge_search -i "./log/bin/perf_v${V}_e${E}.data"
                1>>"./log/annotate/VecStructMgSearchSIMD_v${V}_e${E}.txt" perf annotate --stdio --dsos=runtime --symbol=approxChol_vector2_struct_merge_search_simd -i "./log/bin/perf_v${V}_e${E}.data"
            fi
        fi
    done
done
