#!/usr/bin/perl


use strict;
use warnings;


open(FILE,$ARGV[0]) or die "cannot open ".$ARGV[0];
my $line = <FILE>;
print $line;
my $state=0; #0 before block, 1=in block
my $lenblock;
my $lastghost=0;
my $flushedlastghost=0;

while($line = <FILE>){
  chomp($line);


  if($state==0){
    if(length($line) == 0){
      print $line."\n";
    }else{
      if($line =~ /locus(\d+)\s(\d+)\s(\d+)$/){
	print "locus".$1." ".($2+1)." ".$3."\n";
	$state    =1;
	$lastghost=0;
	$lenblock =$3;
	$flushedlastghost=0;
      }
    }
  }else{
    if($state==1){
      if(length($line) == 0){#end of state
	print "ghost".($lastghost+1)."\t"."N"x$lenblock."\n";
	print $line."\n";
	$flushedlastghost=1;
	$state=0;
      }else{
	if($line =~ /ghost(\d+)\s/){
	  $lastghost=$1;
	  print $line."\n";
	}else{
	  print $line."\n";
	}
      }
    }else{
      die "internal error\n";
    }
  }

}
close(FILE);

if(!$flushedlastghost){
  print "ghost".($lastghost+1)."\t"."N"x$lenblock."\n";
}
