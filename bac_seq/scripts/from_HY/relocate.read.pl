#!usr/bin/perl

use Getopt::Long;
$test = &GetOptions("out=s{1}" => \$output);

while(<>){
    chomp;
    if ($_ =~ /^>/){
        $name = $_;
        $hash1{$name} = "";
        if($name =~ /.*-1/){
            push @name1, $name;
        }
        if($name =~ /.*-2/){
            push @name2, $name;
        }
    }
    if($_ =~ /^\w/){
        $hash1{$name} .= $_;
    }
}

foreach $qname1S (@name1){
    open FILE1, ">>BAC.$output.1.fasta";
    print FILE1 "$qname1S\n$hash1{$qname1S}\n";
} 
foreach $qname2S (@name2){
    open FILE2, ">>BAC.$output.2.fasta";
    print FILE2 "$qname2S\n$hash1{$qname2S}\n";
}

