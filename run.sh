#!/bin/bash

#seq FIRST STEP LAST
vertice=(10 100) 
edge=(10 200)
logfile="./"

for V in ${vertice[@]}
do
    for E in ${edge[@]} 
    do 
	if [ $E -lt $(($V*$V)) ]
	then
	    echo "Start test. V = ${V}, E = ${E}" | tee -a "./log.txt"
	    perf stat -d ./build/src/main $V $E | tee -a "./log.txt"
            perf stat -e r5302b1,r53010e,r5302c7,r5300c4,r5308c7,r5320c7,r5301c7,r5304c7,r5310c7,r5303c7,r533cc7,r532ac7,r5315c7 ./build/src/main $V $E | tee -a "./log.txt" 
            echo "========================" | tee -a "./log.txt" 
	fi    
    done
done


