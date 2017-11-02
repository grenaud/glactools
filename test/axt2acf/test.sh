#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Converting AXT to ACF:"
../glactools axt2acf --fai chrM.fa.fai chrM Pantro5 chrM.hg19.panTro5.net.axt.gz  2> /dev/null > axt.acf.gz
echo -e " ${GREEN}ok${NC}"

echo -n "testing md5sum:"
../glactools view axt.acf.gz |md5sum > axt_acf.md5sum

if diff axt_acf.md5sum axt_acf.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"    
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
