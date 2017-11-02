#!/bin/bash

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

MD5=md5sum


echo -n "Testing view:"
../glactools view 2arch3modern_inter.acf.gz  2> /dev/null > 2arch3modern_inter.view.md5sum
echo -e " ${GREEN}ok${NC}"

echo -n "Testing view -p:"
../glactools view 2arch3modern_inter.acf.gz  2> /dev/null > 2arch3modern_inter.viewp.md5sum
echo -e " ${GREEN}ok${NC}"

echo -n "Testing view -P:"
../glactools view 2arch3modern_inter.acf.gz  2> /dev/null > 2arch3modern_inter.viewP.md5sum
echo -e " ${GREEN}ok${NC}"

echo -n "Testing view -h:"
../glactools view 2arch3modern_inter.acf.gz  2> /dev/null > 2arch3modern_inter.viewh.md5sum
echo -e " ${GREEN}ok${NC}"

echo -n "Testing view -H:"
../glactools view 2arch3modern_inter.acf.gz  2> /dev/null > 2arch3modern_inter.viewH.md5sum
echo -e " ${GREEN}ok${NC}"



echo -n "testing md5sum for view:"
if diff 2arch3modern_inter.view.md5sum 2arch3modern_inter.view.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "testing md5sum for view -p:"
if diff 2arch3modern_inter.viewp.md5sum 2arch3modern_inter.viewp.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "testing md5sum for view -P:"
if diff 2arch3modern_inter.viewP.md5sum 2arch3modern_inter.viewP.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "testing md5sum for view -h:"
if diff 2arch3modern_inter.viewh.md5sum 2arch3modern_inter.viewh.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi

echo -n "testing md5sum for view -H:"
if diff 2arch3modern_inter.viewH.md5sum 2arch3modern_inter.viewH.md5sum_ > /dev/null
then
    echo -e " ${GREEN}test passed${NC}"
else
    echo -e " ${RED}test failed${NC}"
    exit 1
fi
