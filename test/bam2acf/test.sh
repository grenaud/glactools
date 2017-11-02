#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Converting BAM to ACF:"
 ../glactools  bam2acf human_MT.fa tiny.bam  test 2> /dev/null > bam.acf.gz
echo -e " ${GREEN}ok${NC}"

echo -n "testing md5sum:"
../glactools view bam.acf.gz |md5sum > bam_acf.md5sum

if diff bam_acf.md5sum bam_acf.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"    
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
