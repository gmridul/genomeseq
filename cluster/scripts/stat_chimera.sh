#!/bin/bash

../bowtie2-2.2.3/bowtie2-build ../genomeseq/input/$1 $1

../bowtie2-2.2.3/bowtie2 -x $1 -1 ../database/BAC-A-secondDatabaseFirstRE/$1*-1.fq -2 ../database/BAC-A-secondDatabaseFirstRE/$1*-2.fq --un-conc $1.unseq > $1.sam

../bowtie2-2.2.3/bowtie2 --local -x $1 -1 $1.1.unseq -2 $1.2.unseq > $1.local_unmapped.sam

python ../genomeseq/cluster/scripts/cutSlocal.py $1.local_unmapped.sam > $1.localScut.fq

../bowtie2-2.2.3/bowtie2 --local -x $1 -U $1.localScut.fq > $1.detect_chimera.sam

chim=$(awk '{print $3}' < $1.detect_chimera.sam | grep -c $1)
echo $chim > chimera_$1
num_line=$(wc -l $1.detect_chimera.sam | cut -d" " -f1)
echo $num_line-$chim | bc >> chimera_$1
rm $1.*
