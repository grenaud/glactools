

all : 2arch3modern.acf.gz


2arch3modern.acf.gz: AltaiNean.acf.gz Denisova.acf.gz AustralianB.acf.gz FrenchB.acf.gz YorubaB.acf.gz
	../glactools intersect AltaiNean.acf.gz Denisova.acf.gz AustralianB.acf.gz FrenchB.acf.gz YorubaB.acf.gz > 2arch3modern.acf.gz

AltaiNean.acf.gz: human_g1k_v37.fasta.fai 
	wget -O AltaiNea.hg19_1000g.21.mod.vcf.gz http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/AltaiNea.hg19_1000g.21.mod.vcf.gz
	../glactools vcf2acf --fai human_g1k_v37.fasta.fai --epo all.epo.gz AltaiNea.hg19_1000g.21.mod.vcf.gz AltaiNean > AltaiNean.acf.gz

Denisova.acf.gz: human_g1k_v37.fasta.fai 
	wget -O DenisovaPinky.hg19_1000g.21.mod.vcf.gz http://cdna.eva.mpg.de/neandertal/altai/Denisovan/DenisovaPinky.hg19_1000g.21.mod.vcf.gz
	../glactools vcf2acf --fai human_g1k_v37.fasta.fai --epo all.epo.gz DenisovaPinky.hg19_1000g.21.mod.vcf.gz Denisova > Denisova.acf.gz

FrenchB.acf.gz: human_g1k_v37.fasta.fai 
	wget -O SS6004468.hg19_1000g.21.mod.vcf.gz http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004468.hg19_1000g.21.mod.vcf.gz
	../glactools vcf2acf --fai human_g1k_v37.fasta.fai --epo all.epo.gz SS6004468.hg19_1000g.21.mod.vcf.gz FrenchB  > FrenchB.acf.gz

YorubaB.acf.gz: human_g1k_v37.fasta.fai 
	wget -O SS6004475.hg19_1000g.21.mod.vcf.gz http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004475.hg19_1000g.21.mod.vcf.gz
	../glactools vcf2acf --fai human_g1k_v37.fasta.fai --epo all.epo.gz SS6004475.hg19_1000g.21.mod.vcf.gz YorubaB  > YorubaB.acf.gz

AustralianB.acf.gz: human_g1k_v37.fasta.fai 
	wget -O SS6004477.hg19_1000g.21.mod.vcf.gz http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004477.hg19_1000g.21.mod.vcf.gz
	../glactools vcf2acf --fai human_g1k_v37.fasta.fai --epo all.epo.gz SS6004477.hg19_1000g.21.mod.vcf.gz AustralianB  > AustralianB.acf.gz

2arch3modern.glf.gz: AltaiNean.glf.gz Denisova.glf.gz AustralianB.glf.gz FrenchB.glf.gz YorubaB.glf.gz
	../glactools intersect AltaiNean.glf.gz Denisova.glf.gz AustralianB.glf.gz FrenchB.glf.gz YorubaB.glf.gz > 2arch3modern.glf.gz

AltaiNean.glf.gz: human_g1k_v37.fasta.fai 
	wget -O AltaiNea.hg19_1000g.21.mod.vcf.gz http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/AltaiNea.hg19_1000g.21.mod.vcf.gz
	../glactools vcf2glf --fai human_g1k_v37.fasta.fai --epo all.epo.gz AltaiNea.hg19_1000g.21.mod.vcf.gz AltaiNean > AltaiNean.glf.gz

Denisova.glf.gz: human_g1k_v37.fasta.fai 
	wget -O DenisovaPinky.hg19_1000g.21.mod.vcf.gz http://cdna.eva.mpg.de/neandertal/altai/Denisovan/DenisovaPinky.hg19_1000g.21.mod.vcf.gz
	../glactools vcf2glf --fai human_g1k_v37.fasta.fai --epo all.epo.gz DenisovaPinky.hg19_1000g.21.mod.vcf.gz Denisova > Denisova.glf.gz

FrenchB.glf.gz: human_g1k_v37.fasta.fai 
	wget -O SS6004468.hg19_1000g.21.mod.vcf.gz http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004468.hg19_1000g.21.mod.vcf.gz
	../glactools vcf2glf --fai human_g1k_v37.fasta.fai --epo all.epo.gz SS6004468.hg19_1000g.21.mod.vcf.gz FrenchB  > FrenchB.glf.gz

YorubaB.glf.gz: human_g1k_v37.fasta.fai 
	wget -O SS6004475.hg19_1000g.21.mod.vcf.gz http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004475.hg19_1000g.21.mod.vcf.gz
	../glactools vcf2glf --fai human_g1k_v37.fasta.fai --epo all.epo.gz SS6004475.hg19_1000g.21.mod.vcf.gz YorubaB  > YorubaB.glf.gz

AustralianB.glf.gz: human_g1k_v37.fasta.fai 
	wget -O SS6004477.hg19_1000g.21.mod.vcf.gz http://cdna.eva.mpg.de/neandertal/altai/ModernHumans/vcf/SS6004477.hg19_1000g.21.mod.vcf.gz
	../glactools vcf2glf --fai human_g1k_v37.fasta.fai --epo all.epo.gz SS6004477.hg19_1000g.21.mod.vcf.gz AustralianB  > AustralianB.glf.gz

human_g1k_v37.fasta.fai:
	wget -O all.epo.gz http://dna.ku.dk/~gabriel/epo/all.epo.gz
	wget -O all.epo.gz.tbi http://dna.ku.dk/~gabriel/epo/all.epo.gz.tbi
	wget -O human_g1k_v37.fasta.fai ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai

