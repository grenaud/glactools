#!/usr/bin/perl


use strict;
use warnings;
use Data::Dumper;

my $mode=-1;
my $file= $ARGV[0];
my $numberOfPairs=0;

my @arrayDist1mid;
my @arrayDist1min;
my @arrayDist1max;

my @arrayDist2mid;
my @arrayDist2min;
my @arrayDist2max;

my @arrayDist3mid;
my @arrayDist3min;
my @arrayDist3max;

my @arrayDist4mid;
my @arrayDist4min;
my @arrayDist4max;

my @arrayDist5mid;
my @arrayDist5min;
my @arrayDist5max;

my @arrayDist6mid;
my @arrayDist6min;
my @arrayDist6max;

my %uniqpops;
my %uniqcomp;

open(FILE, "gunzip -c $file |") || die "cannot open pipe to $file";
while(my $line = <FILE>){
  chomp($line);
  #print $line;

  if($line =~ /^$/){
    next;
  }
  if($line =~ /^---------/){
    next;
  }

  my @array=split("\t",$line);

  #print $array[ 0]."\n";
  if($array[0] =~ /^(\S+)-(\S+)$/){
    my $pop1=$1;
    my $pop2=$2;
    warn $array[0]."\n";
    if($pop1 eq "root" || $pop2 eq "root" ){  warn "skipping line ".$array[0]."\n"; next; }
    if($pop1 eq "anc"  || $pop2 eq "anc"  ){  warn "skipping line ".$array[0]."\n"; next; }
    if($pop1 eq "ref"  || $pop2 eq "ref"  ){  warn "skipping line ".$array[0]."\n"; next; }

    if(exists $uniqcomp{$pop2."-".$pop1}){    warn "skipping dupline ".$array[0]."\n"; next; }
    $uniqcomp{$pop1."-".$pop2}=1;
    if(!exists $uniqpops{$pop1}){ $uniqpops{$pop1}=1; }
    if(!exists $uniqpops{$pop2}){ $uniqpops{$pop2}=1; }
  }


  if($#array == 126){
    $mode=2;
    push(@arrayDist1mid, $array[ 19*1+0]);
    push(@arrayDist1min, $array[ 19*1+1]);
    push(@arrayDist1max, $array[ 19*1+2]);

    push(@arrayDist2mid, $array[ 20*2+0]);
    push(@arrayDist2min, $array[ 20*2+1]);
    push(@arrayDist2max, $array[ 20*2+2]);

    push(@arrayDist3mid, $array[ 20*3+1]);
    push(@arrayDist3min, $array[ 20*3+2]);
    push(@arrayDist3max, $array[ 20*3+3]);

    push(@arrayDist4mid, $array[ 20*4+2]);
    push(@arrayDist4min, $array[ 20*4+3]);
    push(@arrayDist4max, $array[ 20*4+4]);

    push(@arrayDist5mid, $array[ 20*5+3]);
    push(@arrayDist5min, $array[ 20*5+4]);
    push(@arrayDist5max, $array[ 20*5+5]);

    #print "---\n";

    $numberOfPairs++;
  } else {
    if ($#array == 114) {
      $mode=1;

      push(@arrayDist1mid, $array[ 19*1]);
      push(@arrayDist2mid, $array[ 19*2]);
      push(@arrayDist3mid, $array[ 19*3]);
      push(@arrayDist4mid, $array[ 19*4]);
      push(@arrayDist5mid, $array[ 19*5]);
      push(@arrayDist6mid, $array[ 19*6]);

      $numberOfPairs++;
    } else {
      die "ERROR in line, the number of fields should be 127 or 115 line->".$line."<-, contains ".($#array+1)."\n";
    }
  }
}

close(FILE);


if($mode==1){
  my @pi=(0,0,0,0,0,0);
  #die (Dumper %uniqpops);
  my $freq=1/(scalar keys %uniqpops);
  $freq =   $freq * $freq;
  #essentially the mean for now
  for(my $i=0;$i<=$#arrayDist1mid;$i++){
    $pi[0]  += ($freq*$arrayDist1mid[$i]);
    $pi[1]  += ($freq*$arrayDist2mid[$i]);
    $pi[2]  += ($freq*$arrayDist3mid[$i]);
    $pi[3]  += ($freq*$arrayDist4mid[$i]);
    $pi[4]  += ($freq*$arrayDist5mid[$i]);
    $pi[5]  += ($freq*$arrayDist6mid[$i]);
  }

  print $pi[0]."\t".$pi[1]."\t".$pi[2]."\t".$pi[3]."\t".$pi[4]."\t".$pi[5]."\n";
}

if($mode==2){
  my @pi=(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
  #die (Dumper %uniqpops);
  my $freq=1/(scalar keys %uniqpops);
  $freq =   $freq * $freq;
  #essentially the mean for now
  for(my $i=0;$i<=$#arrayDist1mid;$i++){
    $pi[0]  += ($freq*$arrayDist1mid[$i]);
    $pi[1]  += ($freq*$arrayDist1min[$i]);
    $pi[2]  += ($freq*$arrayDist1max[$i]);

    $pi[3]  += ($freq*$arrayDist2mid[$i]);
    $pi[4]  += ($freq*$arrayDist2min[$i]);
    $pi[5]  += ($freq*$arrayDist2max[$i]);

    $pi[6]  += ($freq*$arrayDist3mid[$i]);
    $pi[7]  += ($freq*$arrayDist3min[$i]);
    $pi[8]  += ($freq*$arrayDist3max[$i]);

    $pi[9]  += ($freq*$arrayDist4mid[$i]);
    $pi[10] += ($freq*$arrayDist4min[$i]);
    $pi[11] += ($freq*$arrayDist4max[$i]);

    $pi[12] += ($freq*$arrayDist5mid[$i]);
    $pi[13] += ($freq*$arrayDist5min[$i]);
    $pi[14] += ($freq*$arrayDist5max[$i]);


  }

  print $pi[0]."\t(".$pi[1].",".$pi[2].")\t".$pi[3]."\t(".$pi[4].",".$pi[5].")\t".$pi[6]."\t(".$pi[7].",".$pi[8].")\t".$pi[9]."\t(".$pi[10].",".$pi[11].")\t".$pi[12]."\t(".$pi[13].",".$pi[14].")\n";
}


