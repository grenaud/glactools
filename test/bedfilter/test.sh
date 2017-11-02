#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Filtering bedfilter sites:"
../glactools bedfilter AltaiNean.acf.gz filter.bed 2> /dev/null  > AltaiNean_filter.acf.gz
echo -e " ${GREEN}ok${NC}"


../glactools view AltaiNean_filter.acf.gz |md5sum > AltaiNean_filter.acf.md5sum

echo -ne "testing md5sum:"

if diff AltaiNean_filter.acf.md5sum AltaiNean_filter.acf.md5sum_  > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
