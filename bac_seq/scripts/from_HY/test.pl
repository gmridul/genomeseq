#!usr/bin/perl
$intput="kabhi jo badal barse";
$intput=~ s/b[ad]dal/$&\n/g;
print($intput);
