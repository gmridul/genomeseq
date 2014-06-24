#!/bin/bash

for f in AC*
do
    tr -d '\n' < $f > line$f
    perl -pe 's/CATG/ "CATG\n>seq " . ++$n . "\n" /ge' line$f > final$f 
    rm line$f
    mv final$f $f
    sed -i '1s/^/>seq 0\n/' $f
done
