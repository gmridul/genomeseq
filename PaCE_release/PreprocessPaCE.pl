 #*****************************************************************************
 #  If you are a U.S. Governmental agency or a non-profit educational
 #  institution, you may use this Product, royalty free, if you agree to
 #  the below terms.  If you agree to these terms, contact Dr. Srinivas Aluru
 #  at aluru@iastate.edu with the following information:
 #  (i)   that you agree with the terms of this license;
 #  (ii)  identification of the U.S. Governmental agency or non-profit
 #  educational institution agreeing to this license; and
 #  (iii) the contact name, address, phone number, and email at your place of
 #  business.
 #  All other interested parties should contact licensing@iastate.edu.
 #  ________________________________________________________________________
 #                  Research License for PaCE Software ("Product")
 #  Copyright ¨ 2003 Iowa State University Research Foundation, Inc.
 #  All rights reserved
 #  This material is based on work supported by the National Science
 #  Foundation under Grant Nos. EIA-0130861 and ACI-0203782.
 #  READ THIS LICENSE AGREEMENT CAREFULLY BEFORE USING THIS PRODUCT. BY
 #  USING THIS PRODUCT YOU INDICATE YOUR ACCEPTANCE OF THE TERMS OF THIS
 #  LICENSE. THESE TERMS APPLY TO YOU AND ANY SUBSEQUENT USER-LICENSEE
 #  OF THIS PRODUCT.
 #  Iowa State University Research Foundation, Inc. ("Licensor") retains the
 #  ownership of this copy and any subsequent copies of the Product. Licensor
 #  grants to Licensee, a U.S. Governmental agency or a non-profit educational
 #  institution ("Licensee"), a non-exclusive, royalty free, non-transferable
 #  license to use the copy of the Product in accordance with the terms and
 #  conditions of this License Agreement.
 #  1. Permitted Uses.  Licensee may:
 #  a) use the Product solely for Licensee's own internal research purposes.
 #  b) alter, modify, or adapt the Product for Licensee's own internal
 #  research purposes.
 #
 #  2. Prohibited Uses.  Licensee may not:
 #  a) transfer, distribute, lease or sub-license the Product.
 #  b) use the Product for any commercial purpose.
 #  c) charge for access or viewing of the Product.
 #  d) alter, modify, or adapt the Product or documentation, or portions
 #  thereof except as provided herein.
 #  Limited Warranty and Liability:
 #  LICENSOR MAKES NO WARRANTY OR REPRESENTATION, EITHER EXPRESS OR IMPLIED,
 #  WITH RESPECT TO THE PRODUCT, INCLUDING ITS QUALITY, PERFORMANCE,
 #  MERCHANTABILITY, OR FITNESS FOR A PARTICULAR PURPOSE OR THAT THE USE OF
 #  THE PRODUCT OR DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS,
 #  TRADEMARKS, OR OTHER RIGHTS. THE PRODUCT PROVIDED HEREUNDER IS ON AN
 #  "AS IS" BASIS.  IN NO EVENT WILL LICENSOR BE LIABLE FOR DIRECT, INDIRECT,
 #  SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR
 #  INABILITY TO USE THE PRODUCT OR DOCUMENTATION, EVEN IF ADVISED OF THE
 #  POSSIBILITY OF SUCH DAMAGES.
 #  General:
 #  Licensor retains all rights not expressly granted herein. Nothing in this
 #  License Agreement constitutes a waiver of Licensor's rights under United
 #  States copyright law. This License and Licensee's right to use the Product
 #  automatically terminate without notice from Licensor if Licensee fails to
 #  comply with any provision of this License Agreement, or any terms and
 #  conditions associated with the transfer of this Product. Upon termination,
 #  you will remove and destroy all copies of the Product. This License
 #  Agreement is governed by the laws of the State of Iowa.
 #*****************************************************************************

#!/usr/bin/perl -w

#USAGE:  PreprocessPaCE.pl {FASTA file}  {1: prune PolyA/T, 0 otherwise}

use strict;

if(! exists $ARGV[1] ) {
	print "#USAGE:  PreprocessPaCE.pl {FASTA file}  {1: prune PolyA/T, 0 otherwise}\n";
	exit;
}

my $n_seq=0;
my $bSkip = 0;
my $gi = "";
my $seq = "";
my $prune = 0;
$prune=1 if($ARGV[1]==1);
my $As = "AAAAAAAAAA";
my $Ts = "TTTTTTTTTT";
my $cutoff_factor = 20; # As/Ts should begin or span the first/last len/20 portions
    
open(Ffile,"$ARGV[0]");
while(<Ffile>) {
         chomp;
         if(/>/) {
			$n_seq++;
			if($n_seq>1) {
				if((length $seq)>0 && $bSkip==0) {
					my $sequence = $seq;
					$sequence = stripTs($seq) if($prune==1);
					if((length $sequence)>0) {
						print "$gi\n";
						print "$sequence\n";
					}
				}
			}
			my @ar = split;
			$gi = $ar[0];
			$seq = "";
			$bSkip=0;
           	next;
         }
	    next if($bSkip==1);
	    my @ar = split;
	    my $l = @ar;
	    if($l >1) {
	    	#print "Skipping Line 1 : #$_#\n";
	    	$seq = "";
		$bSkip = 1;
		next;
	    }

		tr/BbDdEeFfHhIiJjKkLlMmNnOoPpQqRrSsUuVvWwXxYyZz/NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/;

	    my $temp = uc;
	    tr/ACGTN/     /;
	    #tr/ACGTNRYKMSWBVDHU/                /;
	    #tr/acgtnrykmswbvdhu/                /;
	    
	    my @temp_ar = split;
	    my $temp_l = @temp_ar;
	    if($temp_l>0) {
	    	#print "Skipping Line 2 : #$_#\n";
	    	$seq = "";
		$bSkip = 1;
		next;
	    }

	    $seq = $seq . "$temp";
} #end while

if((length $seq)>0 && $bSkip==0) {
	my $sequence = $seq;
	$sequence = stripTs($seq) if($prune==1);
	if((length $sequence)>0) {
		print "$gi\n";
		print "$sequence\n";
	}
}
exit;
# END OF PROGRAM

sub stripTs {
	my $s = shift;

	my $l =length $s;
	return $s if($l<=10);
	my $cutoff = $l/$cutoff_factor ;


	# leading A/Ts
	my $b_A = 1;
	my $leading_i = index $s,$As;
	if($leading_i<0) {
		$leading_i = index $s,$Ts;
		$b_A = 0;
	}
	my $leading_j = -1;
	if($leading_i>=0 && $leading_i<=$cutoff) {
		my $i=$leading_i+ (length $As);
		while($i<$l) {
			last if($b_A==1 && substr($s,$i,1) ne "A");
			last if($b_A==0 && substr($s,$i,1) ne "T");
			$i++;
		}
		$leading_j = $i-1;
	}
	else {
		$leading_i = -1;
		$leading_j = -1;
	}

	# trailing A/Ts
	$b_A = 1;
	my $trailing_i = rindex $s,$As;
	if($trailing_i<0) {
		$trailing_i = rindex $s,$Ts;
		$b_A = 0;
	}
	my $trailing_j=-1;

	if($trailing_i>=0) {
		$trailing_j = $trailing_i + (length $As) - 1;
	}

	if($trailing_i>=0 && $trailing_j>=($l-$cutoff) && $trailing_j<$l ) {
		my $i=$trailing_i - 1;
		while($i>=0 && $i>$leading_j) {
			last if($b_A==1 && substr($s,$i,1) ne "A");
			last if($b_A==0 && substr($s,$i,1) ne "T");
			$i--;
		}
		$trailing_i = $i + 1;
	}
	else {
		$trailing_i = -1;
		$trailing_j = -1;
	}

	#print "[$leading_i ... $leading_j] [$trailing_i ... $trailing_j]\n";

	my $temp_s=$s;
	if($leading_i>=0 && $leading_j<$l) {
		$temp_s = substr($s,$leading_j+1); # remaining string
	}
	my $temp_l = length $temp_s;
	my $diff = $l - $temp_l;

	$trailing_i = $trailing_i - $diff;
	$trailing_j = $trailing_j - $diff;

	$s = $temp_s;
	$l = length $s;

	if($trailing_i>=0 && $trailing_j<$l) {
		$temp_s = substr($s,0,$trailing_i);
	}
	$s = $temp_s;
	$l = length $s;

	return $s;
	
}# end stripTs
