#!/usr/bin/perl


use strict;
use warnings;


my $file= $ARGV[0];
open(FILE, "gunzip -c $file |") || die "cannot open pipe to $file";
while(my $line = <FILE>){
  chomp($line);
  print $line;
}
close(FILE);

