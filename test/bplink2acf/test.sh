#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Testing bplink2acf:"
sed -i "s/root/root1/g" plink.fam
sed -i "s/anc/anc1/g"   plink.fam
../glactools bplink2acf  --fai human_g1k_v37.fasta.fai --epo epochr2.gz plink  2> /dev/null > plink.acf.gz
echo -e " ${GREEN}ok${NC}"

echo -ne "testing md5sum:"
../glactools view plink.acf.gz  |md5sum > plink.acf.md5sum

if diff plink.acf.md5sum plink.acf.md5sum_  > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
