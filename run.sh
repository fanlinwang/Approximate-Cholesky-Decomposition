#!/bin/bash

#seq FIRST STEP LAST
vertice=(10 100 1000) 
edge=(10 100)
logfile="./"

for V in ${vertice[@]}
do
    for E in ${edge[@]} 
    do 
        ./build/src/main $V $E | tee -a "./log.txt" 
        echo "========================" | tee -a "./log.txt" 
    done
done


