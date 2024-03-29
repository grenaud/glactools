
  glactools: command-line toolset for the management of Genotype Likelihoods and Allele Counts
==========================================================

QUESTIONS :
   gabriel [dot] reno [ at sign ] gmail [dot] com


About
----------------------

glactools is a set of command-line tools for the management of Genotype Likelihood (GL) and Allele Counts (AC).


Description
----------------------

glactools is a suite of utilities to:
* convert various file formats (VCF,BAM,23andme) into **genotype likelihood (GLF)** or **allele count (ACF)** files. 
  * GLF files contains genotype likelihoods for a single individual
  * ACF files contains allele counts for either a single individual or a group of individuals (population)
* GLF/ACF contain both variant and invariant sites
* GLF/ACF are binary compressed files that borrow heavily from the BAM format.
* filter, combined GLF/ACF from various individuals, merge individuals into populations
* create subsets (only a single population or retain transversions)
* index for rapid retrieval
* compute summary statistics on those matrices. 
* export genotype likelihood (GLF) or allele count (ACF) to various formats for population genetics applications (treemix,fasta,EIGENSTRAT,G-PhoCS,PLINK).

![diagram](https://github.com/grenaud/glactools/raw/master/doc/diagram.png "Diagram")


Examples of uses
----------------------

glactools aims at allowing users to convert genetic data into an intermediate format (either GLF for genotype likelihoods or ACF for allele counts), perform operations and export to different format. This set of tools enables users to perform various tasks without knowledge of scripting, here are some examples:

* Compute D-statistics but use the Gorilla from UCSC as the ancestral allele.
* Retain sites in the human genome where a certain population has the ancestral allele and a different population has the derived allele.
* Produce Treemix input using a mixture of BAM and VCF files using transversions only
* Combine a BAM file from a single individual and 1000Genomes data to G-PhoCS or ADMIXTURE input



Downloading:
----------------------

Go to https://github.com/grenaud/glactools and either:

1) Download the ZIP 

or

2) Do a "git clone --depth 1 https://github.com/grenaud/glactools.git"

Installation
----------------------

1) make sure you have "git" installed, check for it by typing " git --version" and "cmake --version".  

For Ubuntu:

     sudo apt-get install git

For MacOS, if you have Homebrew (https://brew.sh/) installed: 

     brew install git

2) make sure you have gcc that supports -std=c++11, gcc version 4.7 or above. Type "gcc -v" and check the version. Have the standard zlib, bzlib.h and lzma. Under Ubuntu type:

     sudo apt install g++
     sudo apt install zlib1g-dev
     sudo apt install libbz2-dev
     sudo apt install liblzma-dev
     sudo apt install libcurl4-openssl-dev
     
For MacOS, if you have Homebrew (https://brew.sh/) installed:

     brew install gcc
     brew install zlib
     brew install zv
     
    

3) As the makefile uses "git clone" to download subpackages, please make sure that the computer on which you are installing glactools has access to the internet. Once done, simply type :

     cd glactools
     make

For MacOS, if you get the problem: fatal error: 'lzma.h' file not found, this is a problem building htslib with homebrew, please refer to the following htslib page: https://github.com/samtools/htslib/issues/493


4) (optional) Either put the executable in the overall path or add the path to your $PATH environment or add an alias to be able to call "glactools" from any directory.

Quick start
-----------------
For the impatients, you can download some ACF and GLF data:

     wget -O 2arch3modern.acf.gz https://www.dropbox.com/s/n4su20ghlb3jqni/2arch3modern.acf.gz?dl=0
     wget -O YorubaB.glf.gz https://www.dropbox.com/s/w9af2n6mr9nafw8/YorubaB.glf.gz?dl=0

The first file contains chromosome 21 from 2 archaic homins and 3 modern humans. The second are the genotype likelihoods for a Yoruba individual.

You can now view the first lines:
 
     glactools view 2arch3modern.acf.gz |head -n 20 
     glactools view YorubaB.glf.gz       |head -n 20 

You can index them:

     glactools index 2arch3modern.acf.gz
     glactools index YorubaB.glf.gz       

You can view a chunk:

     glactools view  2arch3modern.acf.gz       21:9675190-9675199 
     glactools view  YorubaB.glf.gz            21:9560830-9560840

You can view the defline:
 
     glactools view -h 2arch3modern.acf.gz |head -n 20 
   
or just view which populations are defined:
 
     glactools view -p 2arch3modern.acf.gz 

or view the full header:
 
     glactools view -P 2arch3modern.acf.gz 

The header is particularly useful as it defines which operations were used to produce this file.

Documentation
-----------------

The documentation is found here:

     doc/reference.pdf


Example of usage
-----------------

We will download 5 different, single individual VCF files as testData (those were mapped to hg19, see "how do I specify the root and ancestral population?" for info about hg38):

    mkdir -p testData/
    cd testData/
    wget http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/AltaiNea.hg19_1000g.21.mod.vcf.gz
    wget http://cdna.eva.mpg.de/neandertal/altai/Denisovan/DenisovaPinky.hg19_1000g.21.mod.vcf.gz
    wget http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004468.hg19_1000g.21.mod.vcf.gz
    wget http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004475.hg19_1000g.21.mod.vcf.gz
    wget http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004477.hg19_1000g.21.mod.vcf.gz
    wget ftp://ftp.healthtech.dtu.dk/public/EPO/all_hg19.epo.gz
    wget ftp://ftp.healthtech.dtu.dk/public/EPO/all_hg19.epo.gz.tbi
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
    cd ..

- Convert the VCF files to ACF files:

      glactools vcf2acf --fai testData/human_g1k_v37.fasta.fai --epo testData/all_hg19.epo.gz testData/AltaiNea.hg19_1000g.21.mod.vcf.gz      AltaiNean    > testData/AltaiNean.acf.gz
      glactools vcf2acf --fai testData/human_g1k_v37.fasta.fai --epo testData/all_hg19.epo.gz testData/DenisovaPinky.hg19_1000g.21.mod.vcf.gz Denisova     > testData/Denisova.acf.gz
      glactools vcf2acf --fai testData/human_g1k_v37.fasta.fai --epo testData/all_hg19.epo.gz testData/SS6004468.hg19_1000g.21.mod.vcf.gz     FrenchB      > testData/FrenchB.acf.gz
      glactools vcf2acf --fai testData/human_g1k_v37.fasta.fai --epo testData/all_hg19.epo.gz testData/SS6004475.hg19_1000g.21.mod.vcf.gz     YorubaB      > testData/YorubaB.acf.gz
      glactools vcf2acf --fai testData/human_g1k_v37.fasta.fai --epo testData/all_hg19.epo.gz testData/SS6004477.hg19_1000g.21.mod.vcf.gz     AustralianB  > testData/AustralianB.acf.gz

- glactools index them:

      glactools index testData/AltaiNean.acf.gz
      glactools index testData/Denisova.acf.gz
      glactools index testData/FrenchB.acf.gz
      glactools index testData/YorubaB.acf.gz
      glactools index testData/AustralianB.acf.gz

- Create an intersection:

      glactools intersect testData/AltaiNean.acf.gz testData/Denisova.acf.gz testData/AustralianB.acf.gz testData/FrenchB.acf.gz testData/YorubaB.acf.gz > testData/2arch3modern.acf.gz

These commands are found in testData/Makefile

- Visualize the intersection:

      glactools view    testData/2arch3modern.acf.gz |less -S # view data
      glactools view -h testData/2arch3modern.acf.gz |less -S # view data+defline
      glactools view -H testData/2arch3modern.acf.gz |less -S # view data+full header


- Merge the modern humans and archaic as one population:

      glactools meld -u testData/2arch3modern.acf.gz   "AltaiNean,Denisova" "Archaics"  |./glactools meld /dev/stdin   "AustralianB,FrenchB,YorubaB" "Modern"  > testData/all.merged.acf.gz

- Get basic statistics:
     
      glactools stats testData/2arch3modern.acf.gz

- Index the file
     
      glactools index testData/2arch3modern.acf.gz

- View a genomic region:

      glactools view testData/2arch3modern.acf.gz 21:25098220-25098230

- View an entire chromosome:

      glactools view testData/2arch3modern.acf.gz 21

- Visualize sites where the archaics hominin and modern ones differ:

      glactools snosharing -u testData/all.merged.acf.gz "Archaics"  "Modern" |./glactools view - 

- Visualize sites where the archaics hominin and modern ones differ and the archaic is ancestral and the modern humans are derived:

      glactools  snosharing -u testData/all.merged.acf.gz "Archaics"  "Modern"  | ./glactools sharing -u /dev/stdin  "root" "Archaics"|./glactools view -

- Visualize sites where the archaics hominin and modern differ and the archaic is derived and the modern humans are ancestral:

      glactools  snosharing -u testData/all.merged.acf.gz "Archaics"  "Modern"  | ./glactools sharing -u /dev/stdin  "root" "Modern"|./glactools view - 

- Export to treemix:

      glactools acf2treemix testData/2arch3modern.acf.gz    |gzip > testData/all.treemix.gz


Problems/feature request
----------------------

If you have a Github account, I recommend that you create an issue. That way other users can see what you wrote, comment on it and I can keep track of it more easily. I more than welcome pull requests!

Otherwise, send me a mail gabriel [dot] reno [at sign here] gmail [dot] com

Tips
----------------------

* when working with VCF files called from bcftools call, make sure that the -v option is not used because this will only print variable sites. If you use GATK, make sure you output every site using --output_mode EMIT_ALL_SITES.
* Do not store ACF/GLF in raw text, it is a waste of disk space.


FAQ
----------------------

### Why do I have a bunch of garbled characters printed on the terminal when I use glactools?

glactools **ALWAYS** prints compressed binary. The only thing to modify is the ability to print as uncompressed binary (-u option), this is recommended when using UNIX pipes. However if you wish to view an ACF/GLF as a text file, simply use "glactools view"

### why do you have data import from single VCF and multi VCF at the same time?

A single VCF usually carries extra information for the single individual such as depth of coverage and additional information in the INFO fields. Ideally you should have a consistent set of filters that does not generate any reference/alternative allele/heterozygous site bias and generate your GLF or ACF files.

### I got the following: Warning: No EOF marker, likely due to an I/O error  [E::bgzf_read] Read block operation failed with error -1 after 0 of 2 bytes

This is likely an input/output error and the file was not written properly to begin with. Try to regenerate it.

### what is the difference between "union" and "intersection"? 

"union" will allow sites to be undefined in a specific population or individual.  "intersection" will require all sides to be defined in every population or individual.

### what is the -u option and what does it do?

how many program there is a -u option which allows users to get a uncompressed glactools output.  This is useful when UNIX piping from one program to another.  If not, one program would compress whereas the second one would decompress,  This is wasteful in terms of CPU.  Therefore when piping into another glactools program, we recommend using the -u.

### what is the difference between root and ancestral (anc)?

The root is an individual or population that is an outgroup to all other individuals/populations in the file. The ancestor is the most recent common ancestor to the root population and all other individuals/populations in a file.

### how do I specify the root and ancestral population

<<<<<<< HEAD
if you're dealing with hominin samples,  we recommend using the -epo option which uses  EPO alignments from Ensembl which are alignments to different primate species. otherwise simply transform  a VCF file  from the roof population  using the program "usepopsrootanc". The EPO information can be found here for hg19:

    ftp://ftp.healthtech.dtu.dk/public/EPO/all_hg19.epo.gz
    ftp://ftp.healthtech.dtu.dk/public/EPO/all_hg19.epo.gz.tbi

and here for hg38:

    ftp://ftp.healthtech.dtu.dk/public/EPO/all_hg38.epo.gz
    ftp://ftp.healthtech.dtu.dk/public/EPO/all_hg38.epo.gz.tbi
=======
if you're dealing with hominin samples,  we recommend using the -epo option which uses  EPO alignments from Ensembl which are alignments to different primate species. otherwise simply transform  a VCF file  from the roof population  using the program "usepopsrootanc". Here is an example:

Suppose you have a file with the populations you want to use as outgroup and ancestor. In the case of hominin, it would be the chimp as outgroup and the Chimp/human ancestor as the ancestor. For example we have the following:

file1.acf.gz:

    #ACF
    #PG:union AltaiNean.acf.gz Denisova.acf.gz FrenchB.acf.gz YorubaB.acf.gz AustralianB.acf.gz 
    #GITVERSION: 6e48e0efbfb131a6723830fa98f655a36da8c574
    #DATE: 2019-12-02
    #SQ	SN:1	LN:249250621
    #SQ	SN:2	LN:243199373
    #chr	coord	REF,ALT	root	anc	chimp	ChimpHumanAncestor
    2	100022	A,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100023	G,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100024	A,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100025	C,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100026	A,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100027	G,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100028	A,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100029	A,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100030	C,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100031	T,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100032	T,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100033	C,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100034	A,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100035	T,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100036	T,N	0,0:0	0,0:0	1,0:0	1,0:0
    2	100037	G,N	0,0:0	0,0:0	1,0:0	1,0:0

And we would have a second file without root and ancestral information:

file2.acf.gz:

    #ACF
    #PG:union AltaiNean.acf.gz Denisova.acf.gz FrenchB.acf.gz YorubaB.acf.gz AustralianB.acf.gz 
    #GITVERSION: 6e48e0efbfb131a6723830fa98f655a36da8c574
    #DATE: 2019-12-02
    #SQ	SN:1	LN:249250621
    #SQ	SN:2	LN:243199373
    #chr	coord	REF,ALT	root	anc	AltaiNean	Denisova	FrenchB	YorubaB	AustralianB
    2	100022	A,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100023	G,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100024	A,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100025	C,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100026	A,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100027	G,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100028	A,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100029	A,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100030	C,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100031	T,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100032	T,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100033	C,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100034	A,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100035	T,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100036	T,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100037	G,N	0,0:0	0,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    
First we will go into file1 and it take our two populations as our root and ancestor and put it into file3.acf.gz:

`glactools usepopsrootanc file1.acf.gz chimp ChimpHumanAncestor  > file3.acf.gz
`

file3.acf.gz will look like this:

    #ACF
    #PG:usepopsrootanc /tmp/temp.acf.gz chimp ChimpHumanAncestor 
    #GITVERSION: 9d3e4107ea445a16737fb841e2181dabf31acac9
    #DATE: 2020-04-02
    #USEPOPASROOTANC: chimp chimp
    #USEPOPASROOTANC#1
    #	#ACF
    #	#PG:union AltaiNean.acf.gz Denisova.acf.gz FrenchB.acf.gz YorubaB.acf.gz AustralianB.acf.gz 
    #	#GITVERSION: 6e48e0efbfb131a6723830fa98f655a36da8c574
    #	#DATE: 2019-12-02
    #	#chr	coord	REF,ALT	root	anc	chimp	ChimpHumanAncestor
    #SQ	SN:1	LN:249250621
    #SQ	SN:2	LN:243199373
    #chr	coord	REF,ALT	root	anc
    2	100022	A,N	1,0:0	1,0:0
    2	100023	G,N	1,0:0	1,0:0
    2	100024	A,N	1,0:0	1,0:0
    2	100025	C,N	1,0:0	1,0:0
    2	100026	A,N	1,0:0	1,0:0
    2	100027	G,N	1,0:0	1,0:0
    2	100028	A,N	1,0:0	1,0:0
    2	100029	A,N	1,0:0	1,0:0
    2	100030	C,N	1,0:0	1,0:0
    2	100031	T,N	1,0:0	1,0:0
    2	100032	T,N	1,0:0	1,0:0
    2	100033	C,N	1,0:0	1,0:0
    2	100034	A,N	1,0:0	1,0:0
    2	100035	T,N	1,0:0	1,0:0
    2	100036	T,N	1,0:0	1,0:0
    2	100037	G,N	1,0:0	1,0:0

Then you can use "replaceanc" to put the ancestral information into file2:

`glactools replaceanc  file2.acf.gz  file3.acf.gz  > file4.acf.gz
`

file4 should contain the same information as file1 except with the ancestral information of file1:

    #ACF
    #PG:replaceanc /tmp/temp2.acf.gz /tmp/file3.acf.gz 
    #GITVERSION: 9d3e4107ea445a16737fb841e2181dabf31acac9
    #DATE: 2020-04-02
    #REPLACEANC:
    #REPLACEANC#1
    #	#ACF
    #	#PG:union AltaiNean.acf.gz Denisova.acf.gz FrenchB.acf.gz YorubaB.acf.gz AustralianB.acf.gz 
    #	#GITVERSION: 6e48e0efbfb131a6723830fa98f655a36da8c574
    #	#DATE: 2019-12-02
    #	#chr	coord	REF,ALT	root	anc	AltaiNean	Denisova	FrenchB	YorubaB	AustralianB
    #REPLACEANC#2
    #	#ACF
    #	#PG:usepopsrootanc /tmp/temp.acf.gz chimp ChimpHumanAncestor 
    #	#GITVERSION: 9d3e4107ea445a16737fb841e2181dabf31acac9
    #	#DATE: 2020-04-02
    #	#USEPOPASROOTANC: chimp chimp
    #	#USEPOPASROOTANC#1
    #	#	#ACF
    #	#	#PG:union AltaiNean.acf.gz Denisova.acf.gz FrenchB.acf.gz YorubaB.acf.gz AustralianB.acf.gz 
    #	#	#GITVERSION: 6e48e0efbfb131a6723830fa98f655a36da8c574
    #	#	#DATE: 2019-12-02
    #	#	#chr	coord	REF,ALT	root	anc	chimp	ChimpHumanAncestor
    #	#chr	coord	REF,ALT	root	anc
    #SQ	SN:1	LN:249250621
    #SQ	SN:2	LN:243199373
    #chr	coord	REF,ALT	root	anc	AltaiNean	Denisova	FrenchB	YorubaB	AustralianB
    2	100022	A,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100023	G,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100024	A,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100025	C,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100026	A,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100027	G,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100028	A,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100029	A,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100030	C,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100031	T,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100032	T,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100033	C,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100034	A,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100035	T,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100036	T,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    2	100037	G,N	1,0:0	1,0:0	2,0:0	2,0:0	2,0:0	2,0:0	2,0:0
    

In this case I have created 1 superfluous files, depending on the level of comfort with Unix file descriptors, you can do the same process without creating intermediate files:

`glactools replaceanc  file2.acf.gz  <(glactools usepopsrootanc -u file1.acf.gz chimp ChimpHumanAncestor )  > file4.acf.gz
`

>>>>>>> 6d59fd67972381a17ace6dafa267e8e49a57f81c


### Can glactools handle data coming from simulations?

Yes, but I recommend using msprime as it can produce directly VCF output.  The one issue is that it does not,  as of this writing, combine individuals together and merely reports haploid data. So one would need to modify the GT field accordingly.

### Can glactools handle BCF?

Yes, simply use "bcftools view" and pipe into glactools as such:

     glactools vcf2acf <(bcftools view in.bcf) put_name_of_sample_here 

### Can glactools handle CRAM?

Yes, simply use "samtools view -b in.cram" and pipe into glactools:

     glactools vcf2acf put_reference_here.fa <(samtools view -b in.cram) put_name_of_sample_here 


### Got the following error message while building the package:  CMake Error at CMakeLists.txt:9 (cmake_minimum_required): CMake 3.0 or higher is required.  You are running version 2.X.XXX

Update cmake, in Ubuntu, the following seems to work:

     sudo -E add-apt-repository -y ppa:george-edison55/cmake-3.x
     sudo -E apt-get update
     sudo apt-get install cmake


### When I import data from EIGENSTRAT, all sites are flagged as non-CpGs, why?

We cannot know for sure whether a C or G had a G or a C after it. 

### When I run "compute", I get 6 types of results, what are the differences?

The summary stats in glactools compute are reported as such:
* all: all sites
* no Cpg: sites marked as CpG are ignored        
* only Cpg: only sites marked as CpG are considered        
* transitions: only transitions are considered        
* transversions: only transversions are considered        
* no damage: sites that are potentially damaged are ignored. The difference between transversions is that ancestral C -> derived T is ignored but not ancestral T -> derived C.
 

### I get: "Error, missing data in the EPO file" when adding the ancestral alleles, why?

The likely reason is a discrepancy between chromosome names (e.g. "1" and "chr1"). 

 

