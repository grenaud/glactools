==========================================================
  glactools: a suite of utilities for the management of Genotype likelihoods and allele counts
==========================================================

QUESTIONS :
   gabriel [dot] reno [ at sign ] gmail [dot] com


About
----------------------

glactools is a set of command-line tools for the management of Genotype likelihood (GL) and allele counts (AC).


Description
----------------------

glactools is a suite of utilities to:
* convert various file formats (VCF,BAM,23andme) into genotype likelihood (GLF) or allele count (ACF) matrices. 
** GLF files contains genotype likelihoods for a single individual
** ACF files contains allele counts for either a single individual or a group of individuals (population)
* GLF/ACF contain both variant and invariant sites
* filter, combined GLF/ACF from various individuals, merge individuals into populations
* create subsets (only a single population or retain transversions)
* index for rapid retrieval
* compute summary statistics on those matrices. 
* export genotype likelihood (GLF) or allele count (ACF) to various formats for population genetics applications (treemix,fasta,EIGENSTRAT,G-PhoCS,PLINK).


Downloading:
----------------------

Go to https://github.com/grenaud/glactools and either:

1) Download the ZIP 

or

2) Do a "git clone https://github.com/grenaud/glactoolstools.git"

Installation
----------------------

1) make sure you have "cmake" and "git" installed, check for it by typing " git --version" and "cmake --version"

2) Make sure you are connected to the internet and type :
   cd glactools
   make



Documentation
-----------------

The documentation is found here:

     doc/reference.pdf


Example of usage
-----------------

We will download 5 different, single individual VCF files as testData:

    mkdir -p testData/
    cd testData/
    wget http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/AltaiNea.hg19_1000g.21.mod.vcf.gz
    wget http://cdna.eva.mpg.de/neandertal/altai/Denisovan/DenisovaPinky.hg19_1000g.21.mod.vcf.gz
    wget http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004468.hg19_1000g.21.mod.vcf.gz
    wget http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004475.hg19_1000g.21.mod.vcf.gz
    wget http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004477.hg19_1000g.21.mod.vcf.gz
    wget http://dna.ku.dk/~gabriel/epo/all.epo.gz
    wget http://dna.ku.dk/~gabriel/epo/all.epo.gz.tbi
    wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai
    cd ..

- Convert the VCF files to ACF files:

    ./glactools vcf2acf --fai testData/human_g1k_v37.fasta.fai --epo testData/all.epo.gz testData/AltaiNea.hg19_1000g.21.mod.vcf.gz      AltaiNean    > testData/AltaiNean.acf.gz
    ./glactools vcf2acf --fai testData/human_g1k_v37.fasta.fai --epo testData/all.epo.gz testData/DenisovaPinky.hg19_1000g.21.mod.vcf.gz Denisova     > testData/Denisova.acf.gz
    ./glactools vcf2acf --fai testData/human_g1k_v37.fasta.fai --epo testData/all.epo.gz testData/SS6004468.hg19_1000g.21.mod.vcf.gz     FrenchB      > testData/FrenchB.acf.gz
    ./glactools vcf2acf --fai testData/human_g1k_v37.fasta.fai --epo testData/all.epo.gz testData/SS6004475.hg19_1000g.21.mod.vcf.gz     YorubaB      > testData/YorubaB.acf.gz
    ./glactools vcf2acf --fai testData/human_g1k_v37.fasta.fai --epo testData/all.epo.gz testData/SS6004477.hg19_1000g.21.mod.vcf.gz     AustralianB  > testData/AustralianB.acf.gz


- glactools index them:
    ./glactools index testData/AltaiNean.acf.gz
    ./glactools index testData/Denisova.acf.gz
    ./glactools index testData/FrenchB.acf.gz
    ./glactools index testData/YorubaB.acf.gz
    ./glactools index testData/AustralianB.acf.gz

- Create an intersection:

    ./glactools intersect testData/AltaiNean.acf.gz testData/Denisova.acf.gz testData/AustralianB.acf.gz testData/FrenchB.acf.gz testData/YorubaB.acf.gz > testData/2arch3modern.acf.gz

These commands are found in testData/Makefile

- Visualize the intersection:

    ./glactools view    testData/2arch3modern.acf.gz |less -S # view data
    ./glactools view -h testData/2arch3modern.acf.gz |less -S # view data+defline
    ./glactools view -H testData/2arch3modern.acf.gz |less -S # view data+full header


- Merge the modern humans and archaic as one population:

     ./glactools meld -u testData/2arch3modern.acf.gz   "AltaiNean,Denisova" "Archaics"  |./glactools meld /dev/stdin   "AustralianB,FrenchB,YorubaB" "Modern"  > testData/all.merged.acf.gz


- Visualize sites where the archaics and modern differ:

     ./glactools snosharing -u testData/all.merged.acf.gz "Archaics"  "Modern" |./glactools view - 

- Visualize sites where the archaics and modern differ and the archaic is ancestral and the modern humans are derived:

     ./glactools  snosharing -u testData/all.merged.acf.gz "Archaics"  "Modern"  | ./glactools sharing -u /dev/stdin  "root" "Archaics"|./glactools view -


- Visualize sites where the archaics and modern differ and the archaic is derived and the modern humans are ancestral:

     ./glactools  snosharing -u testData/all.merged.acf.gz "Archaics"  "Modern"  | ./glactools sharing -u /dev/stdin  "root" "Modern"|./glactools view - 


- Export to treemix:

     ./glactools acf2treemix testData/2arch3modern.acf.gz    |gzip > testData/all.treemix.gz

Problems/feature request
----------------------

If you have a Github account, I recommend that you create an issue. That way other users can see what you wrote comma comment on it and I can keep track of it more easily. I welcome pull requests!

Otherwise, send me a mail gabriel [dot] reno [at sign here] gmail [dot] com


FAQ
----------------------

why do you have data import from single VCF and multi VCF at the same time?

A single VCF usually carries extra information for the single individual such as depth of coverage and additional information in the INFO fields. Ideally you should have a consistent set of filters that does not generate any reference/alternative allele/heterozygous site bias and generate your GLF or ACF files.


