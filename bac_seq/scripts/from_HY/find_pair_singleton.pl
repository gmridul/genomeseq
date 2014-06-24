#!/usr/bin/perl
#File paired and single reads
use Getopt::Long;

$result = &GetOptions("1=s{1}" => \$infile1,
"2=s{1}" => \$infile2,
                      "index|i=s{1}" => \$index);
              if($result){
                  open FILE1, $infile1 or die $!;
                  while(<FILE1>){
                      $n++;
                      chomp;
                      if($_ =~ /^>/){
                          $name1 = $_;
                          $name1 =~ s/-1//;
                          push @name1, $name1;
                          $hash1{$name1};
                      }else{
                          $hash1{$name1} .= $_;
                      }
                  }
#foreach $k (keys %hash1){
#print "$k\n\n$hash1{$k}\n";
#	}
                  open FILE2, $infile2  or die $!;
                  while(<FILE2>){
                      $k++;
                      chomp;
                      if($_ =~ /^>/){
                          $name2 = $_;
                          $name2 =~ s/-2//;
                          push @name2, $name2;
                          $hash2{$name2};
                      }else{
                          $hash2{$name2} .= $_;
                      }
                  }
#foreach $k (keys %hash2){
#print "$k\n$hash2{$k}\n"
#	};
                  @name3 = ();
                  push @name3, @name1, @name2;
#print"@name3\n";
                  foreach $name3 (@name3){
                      $hash3{$name3}++;
                  }
#foreach $k (keys %hash3){
#	print"$k\n$hash3{$k}\n";
#	}
                  foreach $names (sort keys %hash3){
                      if ($hash3{$names}== "2"){
#		$id_1 = $names."-1";
#		$id_2 = $names."-2";
#		print"$id_1\n";
                          if($hash1{$names}){
                              open FILE3, ">>$index.p1.fasta";
                              print FILE3 "$names"."-1"."\n$hash1{$names}\n";
                          }
                          if($hash2{$names}){
                              open FILE4, ">>$index.p2.fasta";
                              print FILE4 "$names"."-2"."\n$hash2{$names}\n";
                          }
                      }
                      if ($hash3{$names}=="1"){
#		$id_s1 = $names."-1";
#		$id_s2 = $names."-2";
                          if ($hash1{$names}){
                              open FILE5, ">>$index.single.fasta";
                              print FILE5 "$names"."-1"."\n$hash1{$names}\n";
                          }
                          if ($hash2{$names}){
                              open FILE5, ">>$index.single.fasta";
                              print FILE5 "$names"."-2"."\n$hash2{$names}\n";
                          }
                      }
                  }
              }

