#!/bin/bash

cd ../input
for f in AC*
do
    tr -d '\n' < $f > line$f
    perl -pe 's/('$1')('$2')/ "$1\n>" . ++$n . "\n$2" /ge' line$f > final$f 
    rm line$f
    mv final$f re$f\.fa
    sed -i '1s/^/>0\n/' re$f\.fa
done
