#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Checking basic stats of an ACF file:"
../glactools stats AltaiNean.acf.gz > AltaiNean.stats
md5sum AltaiNean.stats > AltaiNean.stats.md5sum
echo -e " ${GREEN}ok${NC}"

echo -n "testing md5sum:"
if diff AltaiNean.stats.md5sum AltaiNean.stats.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
