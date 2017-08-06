==========================================================
  glactools: a suite of utilities for the management of Genotype likelihoods and allele counts
==========================================================

QUESTIONS :
   gabriel [dot] reno [ at sign ] gmail.com


About
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

We will download 5 VCF files as testData:

mkdir -p testData/
cd testData/
wget http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/AltaiNea.hg19_1000g.21.mod.vcf.gz
wget http://cdna.eva.mpg.de/neandertal/altai/Denisovan/DenisovaPinky.hg19_1000g.21.mod.vcf.gz
wget http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004468.hg19_1000g.21.mod.vcf.gz
wget http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004475.hg19_1000g.21.mod.vcf.gz
wget http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004477.hg19_1000g.21.mod.vcf.gz

TODO
For the Altai Neandertal, the Denisova, a French individual and a Yoruba individual, all high coverage genomes. 

- Convert the VCF files to MISTARtools files:

./vcf2mistar testData/Altai.vcf.gz    AltaiNeandertal testData/chr21.epo.gz |bgzip -c > testData/Altai.mst.gz 
./vcf2mistar testData/Denisova.vcf.gz Denisova        testData/chr21.epo.gz |bgzip -c > testData/Denisova.mst.gz 
./vcf2mistar testData/French.vcf.gz   French          testData/chr21.epo.gz |bgzip -c > testData/French.mst.gz 
./vcf2mistar testData/Yoruba.vcf.gz   Yoruba          testData/chr21.epo.gz |bgzip -c > testData/Yoruba.mst.gz 


- Tabix index them:

tabix -s 1 -b 2 -e 2 testData/Altai.mst.gz 
tabix -s 1 -b 2 -e 2 testData/Denisova.mst.gz 
tabix -s 1 -b 2 -e 2 testData/French.mst.gz 
tabix -s 1 -b 2 -e 2 testData/Yoruba.mst.gz 

- Create the intersection:

./mistarintersect  testData/{Altai,Denisova,French,Yoruba}.mst.gz |bgzip -c > testData/all.mst.gz


- Merge the modern humans and archaic as one population:

./mistarmeld testData/all.mst.gz   "AltaiNeandertal,Denisova" "Archaics"  |./mistarmeld /dev/stdin   "French,Yoruba" "Modern"|bgzip -c > testData/all.merged.mst.gz


- Visualize sites where the archaics and modern differ:

./mistarfilter  znosharing testData/all.merged.mst.gz "Archaics"  "Modern"


- Visualize sites where the archaics and modern differ and the archaic is ancestral and the modern humans are derived:

./mistarfilter  znosharing testData/all.merged.mst.gz "Archaics"  "Modern"  | ./mistarfilter sharing /dev/stdin  "root" "Archaics"


- Visualize sites where the archaics and modern differ and the archaic is derived and the modern humans are ancestral:

./mistarfilter  znosharing testData/all.merged.mst.gz "Archaics"  "Modern"  | ./mistarfilter sharing /dev/stdin  "root" "Modern"


- Build a neighbor-joining tree:

./mistarcompute nj  --model none testData/all.mst.gz  > testData/all.nw 


- Export to treemix:

./mistar2treemix  testData/all.mst.gz  |gzip > testData/all.treemix.gz


