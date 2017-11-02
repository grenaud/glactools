#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Testing eigen2acf:"
../glactools eigen2acf  --fai human_g1k_v37.fasta.fai --epo epochr2.gz eigenstrat  2> /dev/null > eigenstrat.acf.gz
echo -e " ${GREEN}ok${NC}"

echo -ne "testing md5sum:"
../glactools view eigenstrat.acf.gz  |md5sum > eigenstrat.acf.md5sum

if diff eigenstrat.acf.md5sum eigenstrat.acf.md5sum_  > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
