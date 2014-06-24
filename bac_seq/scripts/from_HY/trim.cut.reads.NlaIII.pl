#!usr/bin/perl

use Getopt::Long;

$result = &GetOptions("in=s{1}" => \$infile1,
"cut=s{1}" => \$infile2,
                      "index|i=s{1}" => \$index);

              open FILE1, $infile1 or die $!;
              while(<FILE1>){
                  chomp;
                  if($_ =~ /^>/){
                      $read_id = $_;
                      $hash1{$read_id};
                  }
                  if($_ =~ /^\w/){
                      $hash1{$read_id} .= $_;
                  }
              }

              open FILE2, $infile2 or die $!;
              while(<FILE2>){
                  chomp;
                  if($_ =~ /^>/){
                      $name = $_;
                      if($name =~ /^(.*)-3$/){
#			print "$1\n";
                          push @chimeric, $1;
                      }
                  }
              }

              open FILE2, $infile2 or die $!;
              while(<FILE2>){
                  chomp;
                  if($_ =~ /^>/){
                      $name1 = $_;
                      $name1 =~ /^(.*)-(\d)$/;
                      $num = $2;
                      $hash2{$name1} = "";
                  }
                  if($_ =~ /^\w/){
                      $hash2{$name1} .= $_;
                      if ($num > 2){
                          open FILE3, ">>BAC.$index.sinlge.fasta";
                          print FILE3 "$name1\n$hash2{$name1}\n";
                      }
                  }
              }

              foreach $chimeric (@chimeric){
                  $hash1{$chimeric} = "";
              }

              foreach $key (sort keys %hash1){
                  if($hash1{$key}){
                      print"$key\n$hash1{$key}\n";
                  }
              }


