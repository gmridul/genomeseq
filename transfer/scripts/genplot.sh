#!/bin/bash

cd ../plots
for f in ../dat/*
do
    e=${f#../dat/}
    e=${e%.dat}
    cat ../scripts/plot.gpt | sed "s/datafile/"$e"/g" | gnuplot
    #gnuplot -e "datafile='"$e"'" plot.gpt 
done
