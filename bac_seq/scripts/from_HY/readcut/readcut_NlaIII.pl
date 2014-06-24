#!usr/bin/perl/

while(<>){
    chomp;
    if($_ =~ /^>/){
        if($name){
            $fragment_len = length($seq);
            $seq =~ s/CATG/CATG\n/g;
            @seq = split/\n/, $seq;
            $r = 0;
            foreach (@seq){
                $r++;
                print"$name-$r\n$_\n";
            }
            $seq = "";
            @seq = "";
        }
        $name = $_;
    }else{
        $seq .= $_;
    }
}

