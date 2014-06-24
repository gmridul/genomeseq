#!usr/bin/perl

while(<>){
	chomp;
	if ($_ =~ /^>/){
		print"$_\n";
		}else{
		$seq_r = reverse($_);
        $seq_r =~ tr/ACGTacgt/TGCAtgca/;
		print"$seq_r\n";
		}
	}

