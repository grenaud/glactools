#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Testing acf2fasta for ACF files:"
../glactools acf2fasta 2arch3modern_inter.acf.gz  2> /dev/null|gzip > acf2fasta.fasta.gz
echo -e " ${GREEN}ok${NC}"

echo -e "testing md5sum: ${YELLOW}not applicable due to use of random()${NC}"
