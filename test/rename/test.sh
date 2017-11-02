#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Testing rename of ACF files:"
../glactools rename 2arch3modern_union.acf.gz "FrenchB,YorubaB,AustralianB" "FrenchBTeam,YorubaBTeam,AustralianBTeam"  2> /dev/null > 2arch3modern_union.rename.acf.gz
echo -e " ${GREEN}ok${NC}"

echo -ne "testing md5sum:"
../glactools view 2arch3modern_union.rename.acf.gz |md5sum > 2arch3modern_union.rename.acf.md5sum

if diff 2arch3modern_union.rename.acf.md5sum 2arch3modern_union.rename.acf.md5sum_  > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
