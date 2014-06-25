#!/bin/bash

cd ../plots
for f in ../input/*1.fq
do
    e=${f#../input/art}
    e=${e%.fa1.fq}
    awk ' NR % 4 == 1 { print $2; }' $f > ../input/$e.dat
    cat ../scripts/plot.gpt | sed "s/datafile/"$e"/g" | gnuplot
    #gnuplot -e "datafile='"$e"'" plot.gpt 
done
