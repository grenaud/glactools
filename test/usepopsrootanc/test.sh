#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Testing usepopsrootanc for ACF files:"
../glactools usepopsrootanc 2arch3modern_union.acf.gz AltaiNean Denisova  2> /dev/null > 2arch3modern_union.usepopsrootanc.acf.gz
echo -e " ${GREEN}ok${NC}"

echo -ne "testing md5sum:"
../glactools view 2arch3modern_union.usepopsrootanc.acf.gz |md5sum > 2arch3modern_union.usepopsrootanc.acf.md5sum

if diff 2arch3modern_union.usepopsrootanc.acf.md5sum 2arch3modern_union.usepopsrootanc.acf.md5sum_  > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
