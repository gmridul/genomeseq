#!/bin/bash

cd ../input
for f in *dat
do
    e=${f%.dat}
    echo $f
    awk 'BEGIN {i=0} {if($1<='$1') {i++}} END {print i/NR}' $f
done

