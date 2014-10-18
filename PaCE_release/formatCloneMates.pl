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

#USAGE:  formatCloneMates.pl {Clone Mates file} {PaCE FASTA file}
#Output:  Converts CMF file containing clone mate information to PaCE-formatted
#         Clone Mates format.
#Input:
#         CMF format:   
#			>clone mate id #
#			  FASTA seq header # (in separate lines)
#			  ....
#Output:
#         PaCE_CMF format:   
#			>clone mate PaCE number#
#			  seq PaCE number # as present in the PaCE FASTA file (in separate lines)
#			  ....


if( exists $ARGV[0] && exists $ARGV[1]) {


    #1. Hash the sequences by their FASTA header and store their PaCE sequence number
 
    $firstword="";
    $thishead="";
    $PaCE_seq_num=0;
    open(FASTAFile,"$ARGV[1]");
    while(<FASTAFile>) {
	 chomp;
         if(/>/) {
		@header=split;
		($firstword,@rest)=@header;
                $thishead = substr($firstword,1);
		#print "$thishead\n";
		$hashEST{$thishead}=$PaCE_seq_num;
		$PaCE_seq_num=$PaCE_seq_num+1;
	 }
    }
    
    #2. Scan each clone mate cluster from the CDF file and output for FASTA sequence its corresponding PaCE seq number.
    
    open(CFILE,"$ARGV[0]");
    while(<CFILE>) {
	chomp;
         if(/>/) {
		print "$_"."\n";
         }
	else {
		@header=split;
		($thishead,@rest)=@header;
		if(defined($hashEST{$thishead})) {
		          print "$hashEST{$thishead}\n";
		}
	}
    } #endwhile
}
else {
	print "#USAGE:   formatCloneMates.pl {Clone Mates file} {PaCE FASTA file}\n";
}

