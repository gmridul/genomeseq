#!/bin/bash

cd ../../input
for f in $(ls --ignore=AC*)
do
    ../emulate/art*/art_illumina -i $f -o art$f -l 250 -amp -p -f $1
done
