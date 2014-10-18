#!/usr/bin/perl -w

#USAGE:  SeqExtractor.pl {FASTA file} 
#Output:  # Extracts FASTA sequence 

use strict;

if( exists $ARGV[0] ) {

    my $n_seq=0;
    my $seq=""; 
    
    open(Ffile,"$ARGV[0]");
    while(<Ffile>) {
         chomp;
         if(/>/) {
			$n_seq++;
			if($n_seq>1 ) {
 				print "$seq\n";
				$seq="";
			}
		   print ">$n_seq\n";
	       #  print "$_\n";
           next;
         }
         if($_ ne "") {
			   if(/>/) {
               	$seq = "";
			  } else {
			  	tr/BbDdEeFfHhIiJjKkLlMmNnOoPpQqRrSsUuVvYyZz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/;
                $seq = $seq . "$_";
			  }
         }
    } #end while
    if($n_seq>=1 ) {
 	print "$seq\n";
    }

}
else {
	print "#USAGE:   SeqExtractor.pl {FASTA file} \n";
}

