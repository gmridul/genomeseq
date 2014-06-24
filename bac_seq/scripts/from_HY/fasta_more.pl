#!/usr/bin/perl

my $id;
my @seq;
my $len;
while(<>){
          chomp;
          if(/^[\>]/){
                      if ($len > 100){
					  $n ++;
                      print"$id\n@seq\n";
					  };
                      @seq="";
                      $len = 0;
                      $id = $_;
                      } else {
                      push @seq, $_;
                      s/\s//g;
                      $len+=length($_);
                      }
          };

