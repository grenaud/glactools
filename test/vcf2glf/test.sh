#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color


MD5=md5sum

echo -n  "Converting the Altai Neanderthal to GLF:"
../glactools vcf2glf --fai human_g1k_v37.fasta.fai --epo epochr2.gz AltaiNea.hg19_1000g.2.mod.vcf.gz AltaiNean 2> /dev/null >  AltaiNean.glf.gz
echo -e " ${GREEN}ok${NC}"
echo -n  "Converting the Denisovan to GLF:"
../glactools vcf2glf --fai human_g1k_v37.fasta.fai --epo epochr2.gz DenisovaPinky.hg19_1000g.2.mod.vcf.gz Denisova 2> /dev/null >  Denisova.glf.gz
echo -e " ${GREEN}ok${NC}"
echo -n  "Converting the B Team French to GLF:"
../glactools vcf2glf --fai human_g1k_v37.fasta.fai --epo epochr2.gz SS6004468.hg19_1000g.2.mod.vcf.gz FrenchB  2> /dev/null >  FrenchB.glf.gz
echo -e " ${GREEN}ok${NC}"
echo -n  "Converting the B Team Yoruba to GLF:"
../glactools vcf2glf --fai human_g1k_v37.fasta.fai --epo epochr2.gz SS6004475.hg19_1000g.2.mod.vcf.gz YorubaB  2> /dev/null >  YorubaB.glf.gz
echo -e " ${GREEN}ok${NC}"
echo -n  "Converting the B Team Australian to GLF:"
../glactools vcf2glf --fai human_g1k_v37.fasta.fai --epo epochr2.gz SS6004477.hg19_1000g.2.mod.vcf.gz AustralianB  2> /dev/null >  AustralianB.glf.gz
echo -e " ${GREEN}ok${NC}"

echo -n  "testing md5sum:"
../glactools view AltaiNean.glf.gz |md5sum 2> /dev/null >  AltaiNean_glf.md5sum

if diff AltaiNean_glf.md5sum AltaiNean_glf.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED} test failed${NC}"
    exit 1
fi
