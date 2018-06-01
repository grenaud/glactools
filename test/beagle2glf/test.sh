#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Testing beagle2glf:"
../glactools beagle2glf  --fai human_MT.fa.fai test_haplo1  2> /dev/null > beagle2glf.glf.gz
echo -e " ${GREEN}ok${NC}"

echo -ne "testing md5sum:"
../glactools view beagle2glf.glf.gz  |md5sum > beagle2glf.glf.md5sum

if diff beagle2glf.glf.md5sum beagle2glf.glf.md5sum_  > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
