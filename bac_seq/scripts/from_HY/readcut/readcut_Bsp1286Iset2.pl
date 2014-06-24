#!usr/bin/perl/
while(<>){
    #print "$name poooooooooooooooooooooooooooooooooooooooooooooooop $seq\n";
    chomp;
    if($_ =~ /^>/){
        if($name){
            $fragment_len = length($seq);
#		if ($seq =~ /([A||G]CATG[T||C])/){
#			$cut_seq = $1;
#			}
            $seq =~ s/(G[A||G||T]GC[T||C||A])C/$1\nC/g;
            @seq = split/\n/, $seq;
            $r = 0;
            foreach(@seq){
                $r++;
                #print"$name-$r\n$_\n";
            }
            $seq = "";			
            @seq = "";
        }
        $name = $_;
    }else{
        $seq .= $_;
    }
}

