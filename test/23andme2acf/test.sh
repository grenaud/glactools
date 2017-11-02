#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Converting 23andme to ACF:"
../glactools 23andme2acf --epo epochr1.gz --fai human_g1k_v37.fasta.fai smallPublic23andMeData.gz anon  2> /dev/null > 23andme.acf.gz
echo -e " ${GREEN}ok${NC}"

echo -n "testing md5sum:"
../glactools view 23andme.acf.gz |md5sum > 23andme_acf.md5sum

if diff 23andme_acf.md5sum 23andme_acf.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"    
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
