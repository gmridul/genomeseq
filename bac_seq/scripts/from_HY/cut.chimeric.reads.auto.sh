#!/bin/bash

for i in `seq 1 -1 1`
do 
	for j in BanII #NlaIII NspI Bsp1286Iset2
	do
	Input=BAC0$i.$j.1.cutadapter.2ends.fastq;
	RE=$j;
	Indexname=adaptTest;

#Trim barcode
#perl trim.bc.pl $Input > $Input.bc;

#Start point check
    perl $RE.start.qc.pl output/$Input > $Input.start;

#Cut reads
#perl readcut_$RE.pl $Input.start > $Indexname.$RE.cut.fasta;
	perl readcut_$RE.pl $Input.start > $Indexname.$RE.cut.fastq;

#Size selection
	perl fasta_more.pl $Indexname.$RE.cut.fastq | perl fasta_more_1.pl > $Indexname.$RE.cut.more20.fasta;
#E.coli Filter
	ncbi-blast-2.2.29+/bin/blastn -query $Indexname.$RE.cut.more20.fasta -db bacteria.all -evalue 1e-5 -perc_identity 98 -num_alignments 1 -num_threads 10 -outfmt 6 > testy
	perl bacteria.result.pl -b testy -c $Indexname.$RE.cut.more20.fasta > $Indexname.$RE.cut.EcoliOut.fasta;
#Musket error correction
	./musket $Indexname.$RE.cut.EcoliOut.fasta -omulti $Indexname.$RE.cut.EcoliOut.ec.fasta -inorder;
	mv $Indexname.$RE.cut.EcoliOut.ec.fasta.0 $Indexname.$RE.cut.EcoliOut.fasta;
#Finding no cutting size reads
	perl trim.cut.reads.$j.pl -in $Input -cut $Indexname.$RE.cut.EcoliOut.fasta -index $Indexname > $Indexname.$RE.nocut.fasta;
	perl relocate.read.pl -out $Indexname.$RE.nocut $Indexname.$RE.nocut.fasta;
	perl find_pair_singleton.pl -1 BAC.$Indexname.$RE.nocut.1.fasta -2 BAC.$Indexname.$RE.nocut.2.fasta -index BAC.$Indexname.Orien.$RE.nocut;


#Move Files
	mkdir BAC$i.$RE.$Indexname;
	mv *Orien* BAC$i.$RE.$Indexname;
#	mv $Input.bc $RE.$Indexname;
	mv $Indexname.$RE.cut.fasta BAC$i.$RE.$Indexname;
	mv $Indexname.$RE.nocut.fasta BAC$i.$RE.$Indexname;
	mv $Indexname.$RE.cut.more20.fasta BAC$i.$RE.$Indexname;
	mv $Indexname.$RE.cut.EcoliOut.fasta BAC$i.$RE.$Indexname;
	mv $Input.start BAC$i.$RE.$Indexname;
#Clear the files
	rm BAC.$Indexname.sinlge.fasta;
	rm BAC.$Indexname.$RE.nocut.2.fasta;
	rm BAC.$Indexname.$RE.nocut.1.fasta;
	done;
done;
