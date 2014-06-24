#!usr/bin/perl/

while(<>){
    chomp;
    if($_ =~ /^>/){
        print "$name\n";
        if($name){
            $fragment_len = length($seq);
            $seq =~ s/(G[A||G]GC[C||T])C/$1\nC/g;
            @seq = split/\n/, $seq;
            $r = 0;
            foreach (@seq){
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
