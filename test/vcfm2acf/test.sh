#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Converting the 1000 Genomes to ACF:"
../glactools vcfm2acf --fai human_g1k_v37.fasta.fai --epo epochr2.gz ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 2> /dev/null > 1000G.acf.gz
echo -e " ${GREEN}ok${NC}"

echo -n  "testing md5sum:"
../glactools view 1000G.acf.gz |md5sum > 1000G_acf.md5sum

if diff 1000G_acf.md5sum 1000G_acf.md5sum_  > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
