#!/bin/bash

RVAL=0
set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

#TODO delete previous data
#
echo "Retrieving data to test:"

echo -n "Retrieving EPO data:"
if wget -O epochr2.gz  -o /dev/null https://www.dropbox.com/s/4gtfqzgawtquq7t/epochr2.gz?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi

echo -n "Retrieving EPO index data:"
if wget -O epochr2.gz.tbi   -o /dev/null https://www.dropbox.com/s/gt6z2m5azwdqgmn/epochr2.gz.tbi?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi

echo -n "Retrieving genome index data:"
if wget -O human_g1k_v37.fasta.fai  -o /dev/null  https://www.dropbox.com/s/nw0at1hq8in42nl/human_g1k_v37.fasta.fai?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi


echo -n "Retrieving test data file 1 of 16:"
if wget -O AltaiNea.hg19_1000g.2.mod.vcf.gz  -o /dev/null https://www.dropbox.com/s/sdxoxdd2t2agl8p/AltaiNea.hg19_1000g.2.mod.vcf.gz?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi


echo -n "Retrieving test data file 2 of 16:"
if wget -O DenisovaPinky.hg19_1000g.2.mod.vcf.gz  -o /dev/null https://www.dropbox.com/s/gjsougijkp1tcdf/DenisovaPinky.hg19_1000g.2.mod.vcf.gz?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi


echo -n "Retrieving test data file 3 of 16:"
if wget -O SS6004468.hg19_1000g.2.mod.vcf.gz  -o /dev/null https://www.dropbox.com/s/lgjsv9vfzk7pqog/SS6004468.hg19_1000g.2.mod.vcf.gz?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi


echo -n "Retrieving test data file 4 of 16:"
if wget -O SS6004475.hg19_1000g.2.mod.vcf.gz  -o /dev/null https://www.dropbox.com/s/f96lfi6gh2a1cd9/SS6004475.hg19_1000g.2.mod.vcf.gz?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi


echo -n "Retrieving test data file 5 of 16:"
if wget -O SS6004477.hg19_1000g.2.mod.vcf.gz  -o /dev/null  https://www.dropbox.com/s/nl1ctxqid9b80mr/SS6004477.hg19_1000g.2.mod.vcf.gz?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi


echo -n "Retrieving test data file 6 of 16:"
if wget -O ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz  -o /dev/null https://www.dropbox.com/s/rn20i08jk22o9ea/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi


echo -n "Retrieving test data file 7 of 16:"
if wget -O human_MT.fa  -o /dev/null https://www.dropbox.com/s/jzy6bzhufuqk85e/human_MT.fa?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi


echo -n "Retrieving test data file 8 of 16:"
if wget -O human_MT.fa.fai   -o /dev/null https://www.dropbox.com/s/7c3r6pvve5w9s1v/human_MT.fa.fai?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi


echo -n "Retrieving test data file 9 of 16:"
if wget -O tiny.bam  -o /dev/null https://www.dropbox.com/s/6wdma5kfgye3fq4/tiny.bam?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi


echo -n "Retrieving test data file 10 of 16:"
if wget -O tiny.bam  -o /dev/null https://www.dropbox.com/s/zvdaws3n8q2fww0/tiny.bed.gz?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi




echo -n "Retrieving test data file 11 of 16:"
if wget -O tiny.bam.bai  -o /dev/null https://www.dropbox.com/s/ektzxgj7gpczhm3/tiny.bam.bai?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi


echo -n "Retrieving test data file 12 of 16:"
if wget -O chrM.hg19.panTro5.net.axt.gz  -o /dev/null https://www.dropbox.com/s/02rygzvagccox94/chrM.hg19.panTro5.net.axt.gz?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi


echo -n "Retrieving test data file 13 of 16:"
if wget -O chrM.fa.fai  -o /dev/null https://www.dropbox.com/s/ve8ualt1b3q9m2f/chrM.fa.fai?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi



echo -n "Retrieving test data file 14 of 16:"
if wget -O smallPublic23andMeData.gz  -o /dev/null https://www.dropbox.com/s/j0vyeepec9uhnbj/smallPublic23andMeData.gz?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi

echo -n "Retrieving test data file 15 of 16:"
if wget -O epochr1.gz  -o /dev/null https://www.dropbox.com/s/zn831wkv3ljb24i/epochr1.gz?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi


echo -n "Retrieving test data file 16 of 16:"
if wget -O epochr1.gz.tbi  -o /dev/null https://www.dropbox.com/s/9m5rl0y5eetyqw0/epochr1.gz.tbi?dl=0; then
    echo -e " ${GREEN}ok${NC}";
else
    echo -e " ${RED}failed${NC}";
    exit 1;
fi




echo "done retrieving files"



echo "Testing vcf2acf"

./vcf2acf/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with vcf2acf exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing vcf2glf"

./vcf2glf/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with vcf2glf exit code: $?${NC}"    
    RVAL=1
fi


echo "Testing vcfm2acf"
./vcfm2acf/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with vcfm2acf exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing bam2acf"
./bam2acf/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with bam2acf exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing axt2acf"
./axt2acf/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with axt2acf exit code: $?${NC}"    
    RVAL=1
fi


echo "Testing 23andme2acf"
./23andme2acf/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with 23andme exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing intersect"
./intersect/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with intersect exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing union"
./union/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with union exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing noundef"
./noundef/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with noundef exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing bedfilter"
./bedfilter/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with bedfilter exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing segsite"
./segsite/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with segsite exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing sharing"
./sharing/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with sharing exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing nosharing"
./nosharing/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with nosharing exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing snosharing"
./snosharing/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with snosharing exit code: $?${NC}"    
    RVAL=1
fi



echo "Testing freqspec"

./freqspec/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with freqspec exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing closest"
./closest/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with closest exit code: $?${NC}"    
    RVAL=1
fi


echo "Testing stats"

./stats/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with stats exit code: $?${NC}"    
    RVAL=1
fi


echo "Testing reheader"

./reheader/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with reheader exit code: $?${NC}"    
    RVAL=1
fi


echo "Testing meld"
./meld/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with meld exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing popsub"
./popsub/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with popsub exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing removepop"
./removepop/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with removepop exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing usepopsrootanc"
./usepopsrootanc/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with usepopsrootanc exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing replaceanc"
./replaceanc/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with replaceanc exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing rename"
./rename/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with rename exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing down"
./down/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with down exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing glac2bed"
./glac2bed/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with glac2bed exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing glac2vcf"
./glac2vcf/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with glac2vcf exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing acf2bplink"
./acf2bplink/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with acf2bplink exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing acf2fasta"
./acf2fasta/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with acf2fasta exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing acf2gphocs"
./acf2gphocs/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with acf2gphocs exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing acf2nexus"
./acf2nexus/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with acf2nexus exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing acf2treemix"
./acf2treemix/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with acf2treemix exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing acf2betascan"
./acf2betascan/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with acf2betascan exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing acf2gross"
./acf2gross/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with acf2gross exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing acf2eigenstrat"
./acf2eigenstrat/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with acf2eigenstrat exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing eigen2acf"
./eigen2acf/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with eigen2acf exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing beagle2glf"
./beagle2glf/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with beagle2glf exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing bplink2acf"
./bplink2acf/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with bplink2acf exit code: $?${NC}"    
    RVAL=1
fi


echo "Testing glf2acf"
./glf2acf/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with glf2acf exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing index"
./index/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with index exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing idxstats"
./idxstats/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with idxstats exit code: $?${NC}"    
    RVAL=1
fi

echo "Testing view"
./view/test.sh
if [ ! $? -eq 0 ] ;then
    echo -e "${RED}Problem with view exit code: $?${NC}"    
    RVAL=1
fi

echo -n "cleaning up:"
#rm -f *md5sum *gz 2*reheader 2*closest 2*freqspec plink.{bed,bim,fam} eigenstrat.{geno,ind,snp}
echo -e "${GREEN}ok${NC}"    

echo "Test completed";
exit ${RVAL}
