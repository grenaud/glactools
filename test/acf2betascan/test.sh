#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Testing acf2betascan for ACF files:"
../glactools acf2betascan --fold 2arch3modern_inter.acf.gz  2> /dev/null |gzip > 2arch3modern_inter.betascan.gz

echo -e " ${GREEN}ok${NC}"

echo -ne "testing md5sum:"
gzip -c -d 2arch3modern_inter.betascan.gz |md5sum > 2arch3modern_inter.betascan.md5sum

if diff 2arch3modern_inter.betascan.md5sum 2arch3modern_inter.betascan.md5sum_  > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
