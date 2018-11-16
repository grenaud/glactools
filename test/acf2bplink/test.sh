#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Testing acf2bplink for ACF files:"
../glactools acf2bplink --homo 2arch3modern_inter.acf.gz plink 2> /dev/null
echo -e " ${GREEN}ok${NC}"

echo -ne "testing md5sum:"
cat plink.fam |md5sum > plink.md5sum

if diff plink.md5sum plink.md5sum_  > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
