#!/usr/bin/perl


use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Cwd 'abs_path';
use File::Temp qw(tempfile);

my $mock=0;

my @humanchr=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22);

sub fileExistsExec{
  my ($exeFile) = @_;

  if (!( -x $exeFile)) {
    die "Executable ".$exeFile." does not exist\n";
  }
}

sub fileExists{
  my ($exeFile) = @_;

  if (!( -e $exeFile)) {
    return 0;
  }else{
    return 1;
  }

}

sub detectCMD{
  my ($cmdtodetect) = @_;

  my $tmpf = new File::Temp( UNLINK => 1 );

  system('bash', '-i', '-c', "type $cmdtodetect > $tmpf");

  open(TMP, $tmpf )     or die "cannot open temp file $tmpf\n";
  while(my $line = <TMP>) {
    chomp $line;
    #print "line $line\n";
    my @array = split(" ",$line);
    my $lastfield = $array[ $#array ];
    #print $lastfield."\n";

    if($line =~ /aliased/){
      #my $cmd=$lastfield;
      my $cmd=join(" ",@array[ 4 ..$#array ]);
      #print Dumper(@array);
      #print "cmd= #".$cmd."#\n";
      if(substr($cmd,0,1) eq "`"){
	$cmd = substr($cmd,1);
      }
      if(substr($cmd,-1) eq "'"){
	$cmd = substr($cmd,0,-1);
      }
      return $cmd;
    }

    if($line =~ /hashed/){
      my $cmd=$lastfield;
      if(substr($cmd,0,1) eq "("){
	$cmd = substr($cmd,1);
      }
      if(substr($cmd,-1) eq ")"){
	$cmd = substr($cmd,0,-1);
      }
      return $cmd;
    }

    #print $lastfield."\n";
    if($#array>=2){
      return $lastfield;
    }

  }
  close(TMP);
  die "The command $cmdtodetect was not detected\n";

}

sub containsPops{
  my ($file,$pops) = @_;
  #todo
  my $cmd = "zcat $file |grep -m 1 \"^#chr\"";
  my $out = `$cmd`;

  chomp $out;
#  warn "#".$out."#\n";
  my @arrayp = split(",",$pops);

  foreach my $pop (@arrayp){
    if(index($out,$pop) == -1){
      die "Population $pop not found in file $file\n";
    }
  }
#  return 1;
}

sub ismstfile{
  my ($file) = @_;

  if(!fileExists($file)){
    #die "Mistar file $file does not exists\n";
    return 0;
  }

  #todo
  my $cmd = "zcat $file |head -1";
  my $out = `$cmd`;
  chomp $out;
  #warn "-".$out."-";

  if($out eq "#MISTAR"){
    return 1;
  }

  return 0
}


my $lengthmerge =      0;
my $minlength   =   1000;
my $outputf     = "none";

my $help;

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "\n\n
 This script is creates a Makefile to run intersect
 on many mst files by spliting on chrs.

\n\n usage:\t".$0." <options> [mst file 1] [populations for file 1] [mst file 2] [populations for file 2] ... [out prefix mst file]


For example:
\t".$0."  pop1.mst.gz pop1,pop2 pop2.mst.gz pop3,pop4 pop3.mst.gz all out.mst.gz
\n\n".

" Options:\n".

# "\t--mock\t\t\t\t\tDo nothing, just print the commands that will be run\n".
#  "\t-l [X bp]\t\t\t\tConsider overlap between regions if they are Xbp away (default ".$lengthmerge." bp)\n".
#  "\t-m [X bp]\t\t\t\tOnly print regions longer than Xbp away (default ".$minlength." bp)\n".
#  "\t-o [output]\t\t\t\tOutput file (default: ".$outputf.")\n".

    "\n\n";
  die;
}
#use \"all\" without the quotes to say use all pops from the mst file



 usage() if ( @ARGV < 1 or
             ! GetOptions('help|?' => \$help)
          or defined $help );


my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);

my $mistarintersect        = $pathdir."/mistarintersect";
fileExistsExec($mistarintersect);
warn "Found $mistarintersect\n";
my $mistarcat        = $pathdir."/mistarcat";
fileExistsExec($mistarcat);
warn "Found $mistarcat\n";

my $mistarfilter        = $pathdir."/mistarfilter";
fileExistsExec($mistarfilter);
warn "Found $mistarfilter\n";


my $tabixcmd = detectCMD("tabix");
warn "Found $tabixcmd\n";

my $bgzipcmd = detectCMD("bgzip");
warn "Found $bgzipcmd\n";

#print detectCMD("dtu");


#print "#".(-e "/home/gabriel/dataeos/khoisan/mst/15304_15304_all.inter.mst.gz")."#\n";
#print "#".fileExists("/home/gabriel/dataeos/khoisan/mst/15304_15304_all.inter.mst.gz")."#\n";
#print $#ARGV;
my @filesToAnalyze;

#for(my $i=0;$i<$#ARGV;$i++){
#  print $i."\t".$ARGV[$i]."\n";
#}
for(my $i=0;$i<$#ARGV;){
  #print $i."\t".$ARGV[$i]."\t".$#ARGV."\n";
  if(ismstfile($ARGV[$i])){

    if(!fileExists($ARGV[$i].".tbi")){
      die "Your file has to be tabix indexed, index not found ".$ARGV[$i].".tbi\n";
    }

    if($i<($#ARGV-1)){

      if(ismstfile($ARGV[$i+1])){#all
	warn "Using file ".$ARGV[$i]." with all pops\n";
	push(@filesToAnalyze,[$ARGV[$i],""]);
	$i++;
	next;
      }else{
	warn "Using file ".$ARGV[$i]." with all the following pops ".$ARGV[$i+1]."\n";
	containsPops($ARGV[$i],$ARGV[$i+1]);
	push(@filesToAnalyze,[ $ARGV[$i],$ARGV[$i+1] ]);
	$i+=2;
	next;
      }
    }else{#all
      	warn "Using file ".$ARGV[$i]." with all pops\n";
      push(@filesToAnalyze,[$ARGV[$i],""]);
      $i++;
      next;
    }

  }else{
    die "ERROR: argument ".$ARGV[$i]." is not a mistar file\n";
  }
}



my $outprefix = $ARGV[ $#ARGV ] ;
my @arrayTargets;
my @arrayAllChrCombined;
my @arrayAllChrUNDEFCombined;

my @arrayDeleteIntermediate;


my $stringPrintlater="";


foreach my $chr (@humanchr){
  #push(@arrayTargets,        $outprefix."_".$chr.".mst.gz.tbi");
  push(@arrayAllChrCombined,     $outprefix."_".$chr.".mst.gz.tbi");
  push(@arrayDeleteIntermediate, $outprefix."_".$chr.".mst.gz");
  push(@arrayDeleteIntermediate, $outprefix."_".$chr.".mst.gz.tbi");


  my $cmd;
  if($#filesToAnalyze == 0){
    $cmd = "cat ";
  }else{
    $cmd = $mistarintersect." ";
  }


  foreach my $rec (@filesToAnalyze){
    if($rec->[1] eq ""){
      $cmd = $cmd. " <(".$tabixcmd." -h ".$rec->[0]." ".$chr." ) ";
    }else{
      $cmd = $cmd. " <(".$tabixcmd." -h ".$rec->[0]." ".$chr." | ".$mistarfilter." popsub /dev/stdin ".$rec->[1]."  ) ";
    }
  }

  $cmd = $cmd. " | ".$bgzipcmd." -c > ".$outprefix."_".$chr.".mst.gz";

  $stringPrintlater=$stringPrintlater."\n".$outprefix."_".$chr.".mst.gz.tbi:\n\t".$cmd."\n\t".$tabixcmd."  -s 1 -b 2 -e 2 ".$outprefix."_".$chr.".mst.gz\n\n";

  #add no undef
  push(@arrayAllChrUNDEFCombined, $outprefix."_".$chr.".noundef.mst.gz.tbi");

  push(@arrayDeleteIntermediate,  $outprefix."_".$chr.".noundef.mst.gz");
  push(@arrayDeleteIntermediate,  $outprefix."_".$chr.".noundef.mst.gz.tbi");


  $stringPrintlater=$stringPrintlater."\n".$outprefix."_".$chr.".noundef.mst.gz.tbi: ".$outprefix."_".$chr.".mst.gz.tbi\n\t".$mistarfilter." noundef  ".$outprefix."_".$chr.".mst.gz | ".$bgzipcmd." -c > ".$outprefix."_".$chr.".noundef.mst.gz\n\t".$tabixcmd."  -s 1 -b 2 -e 2 ".$outprefix."_".$chr.".noundef.mst.gz\n\n";

}



my $cmdcat = $mistarcat." ";
foreach my $chr (@humanchr){
  $cmdcat = $cmdcat. " ".$outprefix."_".$chr.".mst.gz ";
}
$cmdcat = $cmdcat. " | ".$bgzipcmd." -c > ".$outprefix."_all.mst.gz";

$stringPrintlater=$stringPrintlater."\n".$outprefix."_all.mst.gz.tbi: ".join(" ",@arrayAllChrCombined)."\n\t".$cmdcat."\n\t".$tabixcmd."  -s 1 -b 2 -e 2 ".$outprefix."_all.mst.gz\n\n";


push(@arrayTargets,        $outprefix."_all.mst.gz.tbi");
#push(@arrayDeleteCombined, $outprefix."_all.mst.gz.tbi");





#noundef
my $cmdcat = $mistarcat." ";
foreach my $chr (@humanchr){
  $cmdcat = $cmdcat. " ".$outprefix."_".$chr.".noundef.mst.gz ";
}
$cmdcat = $cmdcat. " | ".$bgzipcmd." -c > ".$outprefix."_all.noundef.mst.gz";

$stringPrintlater=$stringPrintlater."\n".$outprefix."_all.noundef.mst.gz.tbi: ".join(" ",@arrayAllChrUNDEFCombined)."\n\t".$cmdcat."\n\t".$tabixcmd."  -s 1 -b 2 -e 2 ".$outprefix."_all.noundef.mst.gz\n\n";

push(@arrayTargets,        $outprefix."_all.noundef.mst.gz.tbi");
#push(@arrayDeleteCombined, $outprefix."_all.noundef.mst.gz.tbi");


print "SHELL := /bin/bash\n\nall:\t".join(" ",@arrayTargets)."\n\ncleanall:\n\trm -vf ".join(" ",@arrayTargets)."\n\ncleanintermediate:\n\trm -vf ".join(" ",@arrayDeleteIntermediate)."\n\n".$stringPrintlater."\n";


#print Dumper(@filesToAnalyze);
