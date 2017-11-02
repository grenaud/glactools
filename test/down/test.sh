#!/bin/bash

set -e

RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Testing down of ACF files:"
../glactools down 2arch3modern_union.acf.gz "FrenchB,YorubaB,AustralianB"  2> /dev/null > 2arch3modern_union.down.acf.gz
echo -e " ${GREEN}ok${NC}"

echo -e "testing md5sum: ${YELLOW}not applicable due to use of random()${NC}"
