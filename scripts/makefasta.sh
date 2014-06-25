#!/bin/bash

cd ../input
# deleting last file in getbac.sh because that is not fetched properly
#sed '2q;d' ../scripts/getbac.sh | awk 'BEGIN { FS=" " } {print $NF}' | xargs rm

for f in AC*
do
    tr -d '\n' < $f > line$f
    perl -pe 's/('$1')('$2')/ "$1\n>" . ++$n . "\n$2" /ge' line$f > final$f 
    rm line$f
    mv final$f $1_$2_$f\.fa
    sed -i '1s/^/>0\n/' $1_$2_$f\.fa
done
