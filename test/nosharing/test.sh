#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Filtering nosharing:"
../glactools nosharing 2arch3modern_union.acf.gz  AltaiNean Denisova 2> /dev/null  > 2arch3modern_union.nosharing.acf.gz
echo -e " ${GREEN}ok${NC}"



echo -e "testing md5sum: ${YELLOW}not applicable due to use of random()${NC}"
