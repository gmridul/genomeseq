#!usr/bin/perl -w 

%barcode = (
	BAC1.BanII.cleaned.paired.1 => 'TGCGATAT',
	BAC1.BanII.cleaned.paired.2 => 'CCTAT',
	BAC2.BanII.cleaned.paired.1 => 'TTAGAGTGTG',
	BAC2.BanII.cleaned.paired.2 => 'CCTAT',
	BAC3.BanII.cleaned.paired.1 => 'ACGAGTG',
	BAC3.BanII.cleaned.paired.2 => 'AGAGAT',
	BAC4.BanII.cleaned.paired.1 => 'GGCTGTGTAT',
	BAC4.BanII.cleaned.paired.2 => 'AGAGAT',
	BAC5.BanII.cleaned.paired.1 => 'TATCTCGT',
	BAC5.BanII.cleaned.paired.2 => 'TTGCGAT',
	BAC6.BanII.cleaned.paired.1 => 'ACATACTC',
	BAC6.BanII.cleaned.paired.2 => 'TTGCGAT',
	BAC7.BanII.cleaned.paired.1 => 'CTGTGTATAC',
	BAC7.BanII.cleaned.paired.2 => 'AACGTCTC',
	BAC8.BanII.cleaned.paired.1 => 'GGTCTCAGT',
	BAC8.BanII.cleaned.paired.2 => 'AACGTCTC',
	BAC9.BanII.cleaned.paired.1 => 'CTACTCG',
	BAC9.BanII.cleaned.paired.2 => 'GTACACCTA',
	BAC10.BanII.cleaned.paired.1 => 'ACGATGCAT',
	BAC10.BanII.cleaned.paired.2 => 'GTACACCTA',
	BAC11.BanII.cleaned.paired.1 => 'CCTAGAC',
	BAC11.BanII.cleaned.paired.2 => 'TACTCGACTG',
	BAC12.BanII.cleaned.paired.1 => 'TGATGTGA',
	BAC12.BanII.cleaned.paired.2 => 'TACTCGACTG',
	BAC1.NlaIII.cleaned.paired.1 => 'GGATACTGT',
	BAC1.NlaIII.cleaned.paired.2 => 'CCTAT',
	BAC2.NlaIII.cleaned.paired.1 => 'CTCTACATGT',
	BAC2.NlaIII.cleaned.paired.2 => 'CCTAT',
	BAC3.NlaIII.cleaned.paired.1 => 'TTGTGTCGC',
	BAC3.NlaIII.cleaned.paired.2 => 'AGAGAT',
	BAC4.NlaIII.cleaned.paired.1 => 'TCTGCGT',
	BAC4.NlaIII.cleaned.paired.2 => 'AGAGAT',
	BAC5.NlaIII.cleaned.paired.1 => 'CATCATCAGC',
	BAC5.NlaIII.cleaned.paired.2 => 'TTGCGAT',
	BAC6.NlaIII.cleaned.paired.1 => 'CTACGAGT',
	BAC6.NlaIII.cleaned.paired.2 => 'TTGCGAT',
	BAC7.NlaIII.cleaned.paired.1 => 'TACTACGCT',
	BAC7.NlaIII.cleaned.paired.2 => 'AACGTCTC',
	BAC8.NlaIII.cleaned.paired.1 => 'ATCACGCAT',
	BAC8.NlaIII.cleaned.paired.2 => 'AACGTCTC',
	BAC9.NlaIII.cleaned.paired.1 => 'TCAGACTG',
	BAC9.NlaIII.cleaned.paired.2 => 'GTACACCTA',
	BAC10.NlaIII.cleaned.paired.1 => 'CGCACAT',
	BAC10.NlaIII.cleaned.paired.2 => 'GTACACCTA',
	BAC11.NlaIII.cleaned.paired.1 => 'TTGTGCAT',
	BAC11.NlaIII.cleaned.paired.2 => 'TACTCGACTG',
	BAC12.NlaIII.cleaned.paired.1 => 'GAGCATAGT',
	BAC12.NlaIII.cleaned.paired.2 => 'TACTCGACTG',
	BAC1.NspI.cleaned.paired.1 => 'TGCATGTCG',
	BAC1.NspI.cleaned.paired.2 => 'CCTAT',
	BAC2.NspI.cleaned.paired.1 => 'CGCATGTCAT',
	BAC2.NspI.cleaned.paired.2 => 'CCTAT',
	BAC3.NspI.cleaned.paired.1 => 'ATATGCTAGC',
	BAC3.NspI.cleaned.paired.2 => 'AGAGAT',
	BAC4.NspI.cleaned.paired.1 => 'TCTGTGAGC',
	BAC4.NspI.cleaned.paired.2 => 'AGAGAT',
	BAC5.NspI.cleaned.paired.1 => 'ATGTGCGATG',
	BAC5.NspI.cleaned.paired.2 => 'TTGCGAT',
	BAC6.NspI.cleaned.paired.1 => 'TCTGTCATGT',
	BAC6.NspI.cleaned.paired.2 => 'TTGCGAT',
	BAC7.NspI.cleaned.paired.1 => 'TCGTAGA',
	BAC7.NspI.cleaned.paired.2 => 'AACGTCTC',
	BAC8.NspI.cleaned.paired.1 => 'TCTAGAGA',
	BAC8.NspI.cleaned.paired.2 => 'AACGTCTC',
	BAC9.NspI.cleaned.paired.1 => 'TCACGCTGTA',
	BAC9.NspI.cleaned.paired.2 => 'GTACACCTA',
	BAC10.NspI.cleaned.paired.1 => 'TTAGTGCTC',
	BAC10.NspI.cleaned.paired.2 => 'GTACACCTA',
	BAC11.NspI.cleaned.paired.1 => 'CACGTCGTAT',
	BAC11.NspI.cleaned.paired.2 => 'TACTCGACTG',
	BAC12.NspI.cleaned.paired.1 => 'ACAGTCGT',
	BAC12.NspI.cleaned.paired.2 => 'TACTCGACTG',
	BAC1.Bsp1286Iset2.cleaned.paired.1 => 'CTAGTACAG',
	BAC1.Bsp1286Iset2.cleaned.paired.2 => 'TTGCGAT',
	BAC2.Bsp1286Iset2.cleaned.paired.1 => 'TCATAGACGC',
	BAC2.Bsp1286Iset2.cleaned.paired.2 => 'TTGCGAT',
	BAC3.Bsp1286Iset2.cleaned.paired.1 => 'TGTGATG',
	BAC3.Bsp1286Iset2.cleaned.paired.2 => 'TTGCGAT',
	BAC4.Bsp1286Iset2.cleaned.paired.1 => 'CGCTACTAGT',
	BAC4.Bsp1286Iset2.cleaned.paired.2 => 'TTGCGAT',
	BAC5.Bsp1286Iset2.cleaned.paired.1 => 'TGTGCATCAC',
	BAC5.Bsp1286Iset2.cleaned.paired.2 => 'TTGCGAT',
	BAC6.Bsp1286Iset2.cleaned.paired.1 => 'TCACTATC',
	BAC6.Bsp1286Iset2.cleaned.paired.2 => 'TTGCGAT',
	BAC7.Bsp1286Iset2.cleaned.paired.1 => 'CTAGTACAG',
	BAC7.Bsp1286Iset2.cleaned.paired.2 => 'AACGTCTC',
	BAC8.Bsp1286Iset2.cleaned.paired.2 => 'AACGTCTC',
	BAC9.Bsp1286Iset2.cleaned.paired.2 => 'AACGTCTC',
	BAC10.Bsp1286Iset2.cleaned.paired.2 => 'AACGTCTC',
	BAC11.Bsp1286Iset2.cleaned.paired.2 => 'AACGTCTC',
	BAC12.Bsp1286Iset2.cleaned.paired.2 => 'AACGTCTC',
	BAC8.Bsp1286Iset2.cleaned.paired.1 => 'TCATAGACGC',
	BAC9.Bsp1286Iset2.cleaned.paired.1 => 'TGTGATG',
	BAC10.Bsp1286Iset2.cleaned.paired.1 => 'CGCTACTAGT',
	BAC11.Bsp1286Iset2.cleaned.paired.1 => 'TGTGCATCAC',
	BAC12.Bsp1286Iset2.cleaned.paired.1 => 'TCACTATC',
	);
#print"$barcode{BAC7.BanII.cleaned.paired.1}\n"
# build RE array
@RE = qw(BanII Bsp1286Iset2 NspI NlaIII);
@BAC = qw(BAC1 BAC2 BAC3 BAC4 BAC5 BAC6 BAC7 BAC8 BAC9 BAC10 BAC11 BAC12);
foreach $re (@RE){
	foreach $bac (@BAC){
		$r_1 = reverse($barcode{$bac.$re.cleaned.paired.1});
		$r_2 = reverse($barcode{$bac.$re.cleaned.paired.2});
        $r_1 =~ tr/ACGTacgt/TGCAtgca/;
        $r_2 =~ tr/ACGTacgt/TGCAtgca/;
        
		system("cutadapt -f fasta -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -g ^$barcode{$bac.$re.cleaned.paired.1} -g ^$barcode{$bac.$re.cleaned.paired.2} -g ^$r_1 -g ^$r_2 -n 2 -e 0.15 $bac.$re.cleaned.paired.all.fasta > $bac.$re.cleaned.paired.cutadapter.fasta");
		system("perl head.trans.tail.pl $bac.$re.cleaned.paired.cutadapter.fasta > $bac.$re.cleaned.paired.cutadapter.R.fasta");
		system("cutadapt -f fasta -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -g ^$barcode{$bac.$re.cleaned.paired.1} -g ^$barcode{$bac.$re.cleaned.paired.2} -g ^$r_1 -g ^$r_2 -n 2 -e 0.15 $bac.$re.cleaned.paired.cutadapter.R.fasta > $bac.$re.cleaned.paired.cutadapter.2ends.fasta");
		system("perl head.trans.tail.pl $bac.$re.cleaned.paired.cutadapter.2ends.fasta > $bac.$re.cutadapter.2ends.fasta");
		#system("mv $bac.$re.cleaned.paired.cutadapter.2ends.fasta $bac.$re.cutadapter.2ends.fasta");
		system("rm $bac.$re.cleaned.paired.cutadapter.*");
		system("mv $bac.$re.cutadapter.2ends.fasta output");
		#error correction#
		#system("./musket -p 24 $bac.$re.cutadapter.2ends.fasta -omulti $bac.$re.cutadapter.2ends.ec.fasta -inorder");
		#system("mv $bac.$re.cutadapter.2ends.ec.fasta.0 $bac.$re.cutadapter.2ends.ec.fasta");
		}
	}

