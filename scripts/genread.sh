#!/bin/bash

cd ../input
for f in *fa
do
    ../emulate/art*/art_illumina -i $f -o art$f -l 250 -amp -p -f $1
done
