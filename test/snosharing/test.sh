#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

MD5=md5sum

echo -n "Filtering snosharing:"
../glactools snosharing 2arch3modern_union.acf.gz  AltaiNean Denisova 2> /dev/null  > 2arch3modern_union.snosharing.acf.gz
echo -e " ${GREEN}ok${NC}"



../glactools view 2arch3modern_union.snosharing.acf.gz |md5sum > 2arch3modern_union.snosharing.acf.md5sum

echo -ne "testing md5sum:"
if diff 2arch3modern_union.snosharing.acf.md5sum 2arch3modern_union.snosharing.acf.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
