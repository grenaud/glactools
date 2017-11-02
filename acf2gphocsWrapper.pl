#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';

my $mock=0;
#my $noo=0;

sub runcmd{
  my ($cmdtorun) = @_;

  warn  "running cmd ". $cmdtorun."\n";
  if($mock != 1){
    my @argstorun = ( "bash", "-c", $cmdtorun );

    if(system(@argstorun) != 0){
      die "system  cmd $cmdtorun failed: $?"
    }else{
    }
  }

}



sub fileExistsExec{
  my ($exeFile) = @_;

  if (!( -x $exeFile)) {
    die "Executable ".$exeFile." does not exist\n";
  }
}

sub fileExists{
  my ($exeFile) = @_;

  if (!( -e $exeFile)) {
    die "Executable ".$exeFile." does not exist\n";
  }
}


my $lengthmerge =0;
my $minlength   =1000;
my $outputf     ="none";
my $allowCpg    = 0;

my $help;

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "\n\n
 This script is a wrapper to convert ACF files into the seq format
 of g-phocs. Normally, you have to run acf2gphocs using a bed file
 of regions to consider. This script finds these continuous regions.

\n\n usage:\t".$0." <options> [ACF file] \n\n".

" Options:\n".

 "\t--mock\t\t\t\t\tDo nothing, just print the commands that will be run\n".
  "\t-l [X bp]\t\t\t\tConsider overlap between regions if they are Xbp away (default ".$lengthmerge." bp)\n".
  "\t-m [X bp]\t\t\t\tOnly print regions longer than Xbp away (default ".$minlength." bp)\n".
  "\t-o [output]\t\t\t\tOutput file (default: ".$outputf.")\n".
  "\t--allowCpg [output]\t\t\t\tDo not mask CpGs (default: ".$outputf.")\n".

    "\n\n";
  die;
}




 usage() if ( @ARGV < 1 or
             ! GetOptions('help|?' => \$help, 'mock' => \$mock, 'allowCpg' => \$allowCpg, 'l=i' => \$lengthmerge, 'm=i' => \$minlength,'o=s' => \$outputf)
          or defined $help );

if($outputf eq "none"){
  die "ERROR: output file must be specified";
}

my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);

my $acf2bed        = $pathdir."/glactools acf2bed";
my $acf2gphocs     = $pathdir."/glactools acf2gphocs";
my $inputACFfile      = $ARGV[ $#ARGV  ];

#fileExistsExec($acf2bed);
#fileExistsExec($acf2gphocs);
fileExists($inputACFfile);

my $cmd;

#call acf2bed





$cmd=$acf2bed." ".$inputACFfile." ";

#detect mergeBed
my $cmdmergebeddetect="which mergeBed";
my $mergebebpath=`$cmdmergebeddetect`;
chomp($mergebebpath);
if(!$mergebebpath){
  die "ERROR: cannot detect command mergeBed, please install it and make it accessible in your path\n";
}


if($lengthmerge>0){
  $cmd = $cmd ." | ".$mergebebpath." -d ".$lengthmerge." -i /dev/stdin "
}
#print $lengthmerge;

#filter by length
if($minlength>0){
  $cmd = $cmd ." | awk '{if( (\$3-\$2)>=".$minlength." ){print \$0}}' ";
}

my $outputftempBED=$outputf."bed_";


#call acf2gphocs
$cmd = $cmd ."  |gzip  > ".$outputftempBED;
runcmd($cmd);


#running acf2phocs
my $outputftemp=$outputf."_";

$cmd = $acf2gphocs;
if($allowCpg){
  $cmd = $cmd." --allowCpg ";
}

$cmd = $cmd."  ".$inputACFfile." <(zcat $outputftempBED) > ".$outputftemp;
runcmd($cmd);


my $numberofLoci=0;
open(FILE,$outputftemp) or die "cannot open ".$outputftemp."\n";
while(my $line = <FILE>){
  chomp($line);
  if($line =~ /^locus(\d+)/){
    $numberofLoci++;
  }
}
close(FILE);



open(FILEO,">".$outputf) or die "cannot write to ".$outputf."\n";
print FILEO $numberofLoci."\n\n";

open(FILE,$outputftemp) or die "cannot open ".$outputftemp."\n";

while(my $line = <FILE>){
  print FILEO $line;
}
close(FILE);

close(FILEO);

unlink($outputftemp) or die "cannot remove temp file ".$outputftemp."\n";


warn "acf2gphocsWrapper.pl finished successfully, wrote data to ".$outputf."\n";
