../bowtie2-2.2.3/bowtie2-build ../genomeseq/oldinput/$1 $1

../bowtie2-2.2.3/bowtie2 -x $1 -1 ../database/BAC-A-secondDatabaseFirstRE/$1*-1.fq -2 ../database/BAC-A-secondDatabaseFirstRE/$1*-2.fq --un-conc $1.unseq > $1.sam

../bowtie2-2.2.3/bowtie2 --local -x $1 -1 $1.1.unseq -2 $1.2.unseq > $1.local_unmapped.sam

python ../genomeseq/cluster/cutSlocal.py $1.local_unmapped.sam > $1.localScut.fq

../bowtie2-2.2.3/bowtie2 --local -x $1 -U $1.localScut.fq > $1.detect_chimera.sam

awk '{print $3}' < detect_chimera.sam | grep -c 'gi'
