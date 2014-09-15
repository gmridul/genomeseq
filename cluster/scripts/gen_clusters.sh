#!/bin/bash
#export PATH=${PATH}:/home/mitra/Downloads/samtools-0.1.19/:../scripts
echo $PATH
samtools view -b -h -S ../../../bowtieindex/AC217266.sam | samtools sort - AC217266_sorted
samtools view -b -h -f 0x10 AC217266_sorted.bam > all_rev_aligned.bam
samtools view -b -h -F 20 AC217266_sorted.bam > all_for_aligned.bam
samtools view -b -h -f 4 AC217266_sorted.bam > all_not_aligned.bam
samtools index all_rev_aligned.bam
samtools index all_for_aligned.bam 
samtools index all_not_aligned.bam
samtools view all_for_aligned.bam | python gen_sam_clusters.py > orig_cluster_forward.txt
samtools view all_rev_aligned.bam | python gen_sam_clusters.py > orig_cluster_reverse.txt
#samtools view all_not_aligned.bam | cut -f 1 >> orig_cluster_forward.txt

