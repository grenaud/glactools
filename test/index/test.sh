#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Creating index of an ACF file:"
../glactools index 2arch3modern_inter.acf.gz 2> /dev/null
echo -e " ${GREEN}ok${NC}"

echo -n "Retrieving region using the index:"
../glactools view 2arch3modern_inter.acf.gz 2:100030-100040 2> /dev/null > 2arch3modern_inter.index.md5sum
echo -e " ${GREEN}ok${NC}"

echo -n "testing md5sum:"
if diff 2arch3modern_inter.index.md5sum  2arch3modern_inter.index.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
