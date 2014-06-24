#!usr/bin/perl

while(<>){
	chomp;
	if($_ =~ /^>/){
			if($name){
		$seq =~ s/ //;
		print"$name\n$seq\n";
		}
		$name = $_;
		}else{
		$seq = $_;
#Bsp1286I 1 CGCTACTAGT 2 ACTAGTAGCG
		$seq =~ s/TGTGCATCAC$//;
		$seq =~ s/TTGCGAT$//;
		$seq =~ s/TGTGCATCAC(\w+)$//;
		$seq =~ s/TTGCGAT(\w+)$//;
#NspI 1 TTAGTGCTC 2 GAGCACTAA
		$seq =~ s/ATGTGCGATG$//;
		$seq =~ s/TTGCGAT$//;
		$seq =~ s/ATGTGCGATG(\w+)$//;
		$seq =~ s/TTGCGAT(\w+)$//;
#NlaIII 1 CGCACAT 2 ATGTGCG
		$seq =~ s/CATCATCAGC$//;
		$seq =~ s/TTGCGAT$//;
		$seq =~ s/CATCATCAGC(\w+)$//;
		$seq =~ s/TTGCGAT(\w+)$//;
#BanII 1 ACGATGCAT 2 ATGCATCGT		
		$seq =~ s/TATCTCGT(\w+)$//;
		$seq =~ s/TTGCGAT(\w+)$//;
		$seq =~ s/TATCTCGT$//;
		$seq =~ s/TTGCGAT$//;
		$seq =~ s/TAGGGGTAC$//;
		$seq =~ s/TAGGGGT$//;
		$seq =~ s/TAGGGTAC$//;
		$seq =~ s/GACATGCATCG$//;
		}
		}

