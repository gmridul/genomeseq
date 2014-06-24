#!/usr/bin/perl -w

use POSIX;

%barcode = (
    "BAC01.BanII.1" => "TGCGATAT",
    "BAC01.BanII.2" => "CCTAT",
    "BAC02.BanII.1" => "TTAGAGTGTG",
    "BAC02.BanII.2" => "CCTAT",
    "BAC03.BanII.1" => "ACGAGTG",
    "BAC03.BanII.2" => "AGAGAT",
    "BAC04.BanII.1" => "GGCTGTGTAT",
    "BAC04.BanII.2" => "AGAGAT",
    "BAC05.BanII.1" => "TATCTCGT",
    "BAC05.BanII.2" => "TTGCGAT",
    "BAC06.BanII.1" => "ACATACTC",
    "BAC06.BanII.2" => "TTGCGAT",
    "BAC07.BanII.1" => "CTGTGTATAC",
    "BAC07.BanII.2" => "AACGTCTC",
    "BAC08.BanII.1" => "GGTCTCAGT",
    "BAC08.BanII.2" => "AACGTCTC",
    "BAC09.BanII.1" => "CTACTCG",
    "BAC09.BanII.2" => "GTACACCTA",
    "BAC10.BanII.1" => "ACGATGCAT",
    "BAC10.BanII.2" => "GTACACCTA",
    "BAC11.BanII.1" => "CCTAGAC",
    "BAC11.BanII.2" => "TACTCGACTG",
    "BAC12.BanII.1" => "TGATGTGA",
    "BAC12.BanII.2" => "TACTCGACTG",
    "BAC01.NlaIII.1" => "GGATACTGT",
    "BAC01.NlaIII.2" => "CCTAT",
    "BAC02.NlaIII.1" => "CTCTACATGT",
    "BAC02.NlaIII.2" => "CCTAT",
    "BAC03.NlaIII.1" => "TTGTGTCGC",
    "BAC03.NlaIII.2" => "AGAGAT",
    "BAC04.NlaIII.1" => "TCTGCGT",
    "BAC04.NlaIII.2" => "AGAGAT",
    "BAC05.NlaIII.1" => "CATCATCAGC",
    "BAC05.NlaIII.2" => "TTGCGAT",
    "BAC06.NlaIII.1" => "CTACGAGT",
    "BAC06.NlaIII.2" => "TTGCGAT",
    "BAC07.NlaIII.1" => "TACTACGCT",
    "BAC07.NlaIII.2" => "AACGTCTC",
    "BAC08.NlaIII.1" => "ATCACGCAT",
    "BAC08.NlaIII.2" => "AACGTCTC",
    "BAC09.NlaIII.1" => "TCAGACTG",
    "BAC09.NlaIII.2" => "GTACACCTA",
    "BAC10.NlaIII.1" => "CGCACAT",
    "BAC10.NlaIII.2" => "GTACACCTA",
    "BAC11.NlaIII.1" => "TTGTGCAT",
    "BAC11.NlaIII.2" => "TACTCGACTG",
    "BAC12.NlaIII.1" => "GAGCATAGT",
    "BAC12.NlaIII.2" => "TACTCGACTG",
    "BAC01.NspI.1" => "TGCATGTCG",
    "BAC01.NspI.2" => "CCTAT",
    "BAC02.NspI.1" => "CGCATGTCAT",
    "BAC02.NspI.2" => "CCTAT",
    "BAC03.NspI.1" => "ATATGCTAGC",
    "BAC03.NspI.2" => "AGAGAT",
    "BAC04.NspI.1" => "TCTGTGAGC",
    "BAC04.NspI.2" => "AGAGAT",
    "BAC05.NspI.1" => "ATGTGCGATG",
    "BAC05.NspI.2" => "TTGCGAT",
    "BAC06.NspI.1" => "TCTGTCATGT",
    "BAC06.NspI.2" => "TTGCGAT",
    "BAC07.NspI.1" => "TCGTAGA",
    "BAC07.NspI.2" => "AACGTCTC",
    "BAC08.NspI.1" => "TCTAGAGA",
    "BAC08.NspI.2" => "AACGTCTC",
    "BAC09.NspI.1" => "TCACGCTGTA",
    "BAC09.NspI.2" => "GTACACCTA",
    "BAC10.NspI.1" => "TTAGTGCTC",
    "BAC10.NspI.2" => "GTACACCTA",
    "BAC11.NspI.1" => "CACGTCGTAT",
    "BAC11.NspI.2" => "TACTCGACTG",
    "BAC12.NspI.1" => "ACAGTCGT",
    "BAC12.NspI.2" => "TACTCGACTG",
);

#@RE = qw(BanII NspI NlaIII);
#@BAC = qw(BAC01 BAC02 BAC03 BAC04 BAC05 BAC06 BAC07 BAC08 BAC09 BAC10 BAC11 BAC12);
@RE = qw(BanII);
@BAC = qw(BAC01);
foreach $bac (@BAC) {
    foreach $re (@RE) {
        $var1 = "${bac}.${re}.1";
        $var2 = "${bac}.${re}.2";
        print STDERR $barcode{$var1}, "-", $barcode{$var2}, "\n";
        #$r_1 = reverseComplement($barcode{$bac.$re.cleaned.paired.1});
        #$r_2 = reverseComplement($barcode{$bac.$re.cleaned.paired.2});
        #system("./cutadapt -f fasta -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -g ^$barcode{$bac.$re.cleaned.paired.1} -g ^$barcode{$bac.$re.cleaned.paired.2} -g ^$r_1 -g ^$r_2 -n 2 -e 0.15 $bac.$re.cleaned.paired.all.fasta > $bac.$re.cleaned.paired.cutadapter.fasta");
        #system("perl head.trans.tail.pl $bac.$re.cleaned.paired.cutadapter.fasta > $bac.$re.cleaned.paired.cutadapter.R.fasta");
        #system("./cutadapt -f fasta -b AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT -b AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -g ^$barcode{$bac.$re.cleaned.paired.1} -g ^$barcode{$bac.$re.cleaned.paired.2} -g ^$r_1 -g ^$r_2 -n 2 -e 0.15 $bac.$re.cleaned.paired.cutadapter.R.fasta > $bac.$re.cleaned.paired.cutadapter.2ends.fasta");
        #system("perl head.trans.tail.pl $bac.$re.cleaned.paired.cutadapter.2ends.fasta > $bac.$re.cutadapter.2ends.fasta");
        #system("rm $bac.$re.cleaned.paired.cutadapter.*");
        #system("mv $bac.$re.cutadapter.2ends.fasta /mnt/02/hungying/Ray-v2.1.0");
    }
}




























































