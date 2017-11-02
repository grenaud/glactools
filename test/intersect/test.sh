#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Creating the intersection of ACF files:"
../glactools intersect AltaiNean.acf.gz  Denisova.acf.gz FrenchB.acf.gz YorubaB.acf.gz AustralianB.acf.gz 2> /dev/null > 2arch3modern_inter.acf.gz
echo -e " ${GREEN}ok${NC}"

echo -n "Creating the intersection of GLF files:"
../glactools intersect AltaiNean.glf.gz  Denisova.glf.gz FrenchB.glf.gz YorubaB.glf.gz AustralianB.glf.gz 2> /dev/null > 2arch3modern_inter.glf.gz
echo -e " ${GREEN}ok${NC}"

echo -n "testing md5sum:"
../glactools view 2arch3modern_inter.acf.gz |md5sum > 2arch3modern_inter_acf.md5sum

if diff 2arch3modern_inter_acf.md5sum 2arch3modern_inter_acf.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "testing md5sum:"
../glactools view 2arch3modern_inter.glf.gz |md5sum > 2arch3modern_inter_glf.md5sum

if diff 2arch3modern_inter_glf.md5sum 2arch3modern_inter_glf.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
