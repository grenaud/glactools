#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Testing glac2bed for ACF files:"
../glactools glac2bed 2arch3modern_inter.acf.gz 2> /dev/null | gzip -c > glac2bed.gz
echo -e " ${GREEN}ok${NC}"

echo -ne "testing md5sum:"
gzip -c -d glac2bed.gz |md5sum > glac2bed.md5sum

if diff glac2bed.md5sum glac2bed.md5sum_  > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
