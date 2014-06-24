#!usr/bin/bash

#This pipeline is used to connected paired reads by COPE software
#And Clustering single end reads
for i in `seq 1 1 12`
do
	for j in NlaIII 
	do
	RE=$j;
#	Indexname=adaptTest;
# Find paired-end reads to COPE
#	cat /mnt/02/hungying/Ray-v2.1.0/BAC$i.$j.adaptTest/BAC.adaptTest.Orien.$j.nocut.p2.fasta > ./BAC$i.nocut.p2.fasta;
#	perl fasta_to_fastq.pl BAC$i.nocut.p2.fasta > BAC$i.nocut.p2.fastq;
#	cat /mnt/02/hungying/Ray-v2.1.0/BAC$i.$j.adaptTest/BAC.adaptTest.Orien.$j.nocut.p1.fasta > ./BAC$i.nocut.p1.fasta;

#Copy reads (Data from the cut.chimerica pipeline)
cp /mnt/02/hungying/Ray-v2.1.0/BAC$i.NlaIII.adaptTest/adaptTest.NlaIII.EcoliOut.fasta ./BAC$i.clean
python subsampler.py -fasta -n 12000 -in ./BAC$i.clean -out Nla.temp
cp /mnt/02/hungying/Ray-v2.1.0/BAC$i.Bsp1286Iset2.bsp.real/bsp.real.Bsp1286Iset2.EcoliOut.fasta ./BAC$i.clean
python subsampler.py -fasta -n 3000 -in ./BAC$i.clean -out Bsp.temp
cat Nla.temp > ./BAC$i.clean;
cat Bsp.temp >> ./BAC$i.clean;
bowtie2 --very-fast -p 12 -f -x R1.bless.750.ref -U BAC$i.clean -S BAC$i.sam;
perl BACseq_contig2wgs.check3.pl -i BAC$i.clean -sam BAC$i.sam > BAC$i.clean.wgscheck.fasta;
perl fasta_more.pl BAC$i.clean.wgscheck.fasta | perl fasta_more_1.pl > BAC$i.clean.wgscheck.more.fasta
perl relocate.read.pl -out BAC$i.clean2.ori BAC$i.clean.wgscheck.more.fasta
perl find_pair_singleton.pl -1 BAC.BAC$i.clean2.ori.1.fasta -2 BAC.BAC$i.clean2.ori.2.fasta -index BAC$i.clean.wgscheck.more.ori

#COPE
perl fasta_to_fastq.pl BAC$i.clean.wgscheck.more.ori.p2.fasta > BAC$i.R1.partial.2.fq
perl fasta_to_fastq.pl BAC$i.clean.wgscheck.more.ori.p1.fasta > BAC$i.R1.partial.1.fq
./cope -a BAC$i.R1.partial.1.fq -b BAC$i.R1.partial.2.fq -o BAC$i.connect.fq -2 BAC$i.nocut.p1.left.fastq -3 BAC$i.nocut.p2.left.fastq -l 30 -u 250 -d 0.95 -m 0;
perl fastq2fasta.pl BAC$i.connect.fq | perl cut_num.pl > BAC$i.connect.fasta;
perl fastq2fasta.pl BAC$i.nocut.p1.left.fastq | perl cut_num.pl > BAC$i.nocut.p1.left.fasta;
perl fastq2fasta.pl BAC$i.nocut.p2.left.fastq | perl cut_num.pl > BAC$i.nocut.p2.left.fasta;

#PS checking and clustering
bowtie2 --very-fast -p 12 -f -x R1.bless.750.ref -U BAC$i.connect.fasta -S BAC$i.sam;
perl BACseq_contig2wgs.check.pl -i BAC$i.connect.fasta -sam BAC$i.sam > BAC$i.connect.wgscheck.fasta;
bowtie2 --very-fast -p 12 -f -x arf.1.ref -U BAC$i.connect.wgscheck.fasta -S BAC$i.sam;
perl BACseq_contig2wgs.check.4.pl -i BAC$i.connect.wgscheck.fasta -sam BAC$i.sam > BAC$i.connect.wgscheck.arf.fasta;
./cd-hit-est -d 100 -c 0.98 -g 1 -i BAC$i.connect.wgscheck.arf.fasta -o BAC$i.connect.wgscheck.arf.est
perl cluster.select.pl -1 BAC$i.connect.wgscheck.arf.est.clstr -2 BAC$i.connect.wgscheck.arf.est -n 1 > BAC$i.connect.wgscheck.arf.est4

#single end 
cat BAC$i.clean.wgscheck.more.ori.single.fasta > BAC$i.single.all.fa;
cat BAC$i.nocut.p1.left.fasta >> BAC$i.single.all.fa;
cat BAC$i.nocut.p2.left.fasta >> BAC$i.single.all.fa;
./cd-hit-est -d 100 -c 0.99 -g 1 -i BAC$i.single.all.fa -o BAC$i.single.all.est;
#perl cluster.select.pl -1 BAC$i.single.all.est.clstr -2 BAC$i.single.all.est -n 0 > BAC$i.single.est8;
#In this version, I use cap3 to calculate consensus
perl cap3.consensus.pl -I BAC$i.single.all.fa -c BAC$i.single.all.est.clstr -ch 40 -co 50 -cp 98 -o contigs.temp;
perl change_fastq.tasr.pl contigs.temp > BAC$i.cap3.consen.fa;
bowtie2 --very-fast -p 12 -f -x R1.bless.750.ref -U BAC$i.cap3.consen.fa -S BAC$i.sam;
perl BACseq_contig2wgs.check2.pl -i BAC$i.cap3.consen.fa -sam BAC$i.sam > BAC$i.single.est8.wgscheck;
echo "/mnt/02/hungying/cope/BAC$i.clean.wgscheck.more.fasta" > tasr.fof
./TASR -f tasr.fof -s BAC$i.single.est8.wgscheck -b BAC$i.single.est8.tasr -m 100 -o 4 -r 0.8 -u 1 -k 35; 
#bowtie2 --very-fast -p 12 -f -x R2.bless.750.ref -U BAC$i.single.all.fa -S BAC$i.sam;

#Remove 
rm contigs.temp
rm BAC$i.clean
rm BAC$i.clean.wgscheck.fasta;
rm BAC$i.nocut.p1.left.*
rm BAC$i.nocut.p2.left.*
rm BAC$i.R1.partial.1.fq
rm BAC$i.R1.partial.2.fq
rm BAC$i.connect.wgscheck.fasta
rm BAC.BAC$i.clean2.ori.*
rm BAC$i*.log
rm BAC$i.sam
#Moving
mkdir BAC$i.BACseq.$j.sub1
mv BAC$i.clean.wgscheck* ./BAC$i.BACseq.$j.sub1
mv BAC$i.connect* ./BAC$i.BACseq.$j.sub1
mv BAC$i.nocut.p* ./BAC$i.BACseq.$j.sub1
mv BAC$i.single.* ./BAC$i.BACseq.$j.sub1
mv BAC$i.single.est10.tasr* ./BAC$i.BACseq.$j.sub1
done
done

