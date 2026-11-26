#!/bin/bash
# Offline tests of "glactools compute" on small made-up data, with expected values that can be checked by hand.
#   test/testCompute.sh [path/to/glactools]
set -euo pipefail
G=$(realpath "${1:-$(dirname "$0")/../glactools}")
W=$(mktemp -d); trap 'rm -rf "$W"' EXIT; cd "$W"

fail() { echo "FAILED: $1" >&2; exit 1; }
header() { printf '#ACF\n#PG:test\n#GITVERSION: x\n#DATE: x\n#SQ\tSN:1\tLN:1000000\n#chr\tcoord\tREF,ALT\troot\tanc'; printf '\t%s' "$@"; printf '\n'; }
binary() { "$G" view -b "$1.txt" > "$1.acf.gz" 2>/dev/null; }
compute() { timeout 120 "$G" compute "$@" 2>/dev/null; }
# field NAME N: the Nth non-empty value of the row called NAME (read from stdin)
field() { awk -F'\t' -v name="$1" -v n="$2" '$1==name { c=0; for(i=2;i<=NF;i++) if($i!="" && ++c==n) { print $i; exit } }'; }
equals() { [ -n "$2" ] && awk -v a="$2" -v b="$3" 'BEGIN{ d=a-b; if(d<0) d=-d; exit !(d<1e-5) }' || fail "$1: got '$2', expected $3"; }

# f2 and f3: one site, 4 alleles each, reference allele frequency P1=0.75, P2=0.25, P3=1
{ header P1 P2 P3; printf '1\t100\tA,C\t1,0:0\t1,0:0\t3,1:0\t1,3:0\t4,0:0\n'; } > one.txt; binary one
compute -p f2 --noboot one.acf.gz > f2.out || fail "f2 does not run"
equals "f2(P1,P2) sum"       "$(field 'P1 - P2' 1 < f2.out)" 0.25
equals "f2(P1,P2) corrected" "$(field 'P1 - P2' 8 < f2.out)" 0.15625     # 0.25 - 0.1875/4 - 0.1875/4
equals "f2(P2,P3) corrected" "$(field 'P2 - P3' 8 < f2.out)" 0.515625    # 0.5625 - 0.1875/4
compute -p f3 --noboot one.acf.gz > f3.out || fail "f3 does not run"
equals "f3(P1;P2,P3) sum"       "$(field 'P1 ; P2 , P3' 1 < f3.out)" -0.125
equals "f3(P1;P2,P3) corrected" "$(field 'P1 ; P2 , P3' 8 < f3.out)" -0.1875    # -0.125 - 0.1875/3

# -u: a second site where P3 is undefined counts for P1-P2 but not for P1-P3
{ cat one.txt; printf '1\t200\tA,C\t1,0:0\t1,0:0\t2,2:0\t1,1:0\t0,0:0\n'; } > undef.txt; binary undef
compute -u -p f2 --noboot undef.acf.gz > f2u.out
equals "-u f2(P1,P2) sites"     "$(field 'P1 - P2' 2 < f2u.out)" 2
equals "-u f2(P1,P2) corrected" "$(field 'P1 - P2' 8 < f2u.out)" -0.00260417  # (0.25 - 0.4375/4 - 0.4375/3)/2
equals "-u f2(P1,P3) sites"     "$(field 'P1 - P3' 2 < f2u.out)" 1
equals "-u f2(P1,P3) corrected" "$(field 'P1 - P3' 8 < f2u.out)" 0.015625

# dstat noDamage: site 100 is ancestral C -> derived T (to be skipped), sites 200 and 300 are A -> C
{ header X Y Z
  printf '1\t100\tC,T\t2,0:0\t2,0:0\t0,2:0\t2,0:0\t0,2:0\n'
  printf '1\t%s\tA,C\t2,0:0\t2,0:0\t0,2:0\t0,2:0\t2,0:0\n' 200 300; } > dam.txt; binary dam
compute -p dstat --noboot dam.acf.gz > dam.out
# the 6th group of 5 values is noDamage: AA AD DA DD D
got=$(for n in 26 27 28 29; do field 'Y-Z@X' $n < dam.out; done | tr '\n' ' ')
[ "$got" = "0 0 2 0 " ] || fail "dstat noDamage counts: got '$got', expected '0 0 2 0'"

# Jackknife: three identical blocks. Leaving one out changes nothing, so every interval must be
# the estimate itself, and every combined count must be three times the count of one block.
block() {  # every combination of substitution, CpG flag, ancestral allele and derived state of X, Y and Z,
           # each repeated 1 to 5 times so that the categories have different values
    awk 'BEGIN{ OFS="\t"; split("A,C C,T G,A A,G C,A T,C",S," "); p=0
      for(b=0;b<3;b++) for(s=1;s<=6;s++) for(c=0;c<2;c++) for(a=0;a<2;a++) for(m=0;m<8;m++){
        A=(a?"0,2:":"2,0:") c; D=(a?"2,0:":"0,2:") c
        for(r=0;r<=(m*s+m*m*c+a*s+int(m/3))%5;r++)
          print 1, 100+(p++)*10, S[s], A, A, (int(m/4)%2?D:A), (int(m/2)%2?D:A), (m%2?D:A) } }'
}
{ header X Y Z; block; } > blocks.txt; binary blocks
sites=$(( $(block | wc -l) / 3 ))
# check STAT VALUES_PER_CATEGORY COUNTS "column:interval_column ..." (columns are 1-based within a category)
check() {
    compute -p "$1" -s "$sites" blocks.acf.gz > "$1.out" || fail "$1 does not run"
    awk -F'\t' -v stat="$1" -v w="$2" -v ncount="$3" -v pairs="$4" '
      function abs(x){ return x<0?-x:x }
      /^---/ { all = /ALL/; nb += !all; next }
      NF>1 && $1 !~ /ref/ {
        n=0; for(i=2;i<=NF;i++) if($i!="") v[++n]=$i
        if(!all){ if(nb==1) for(i=1;i<=n;i++) first[$1,i]=v[i]; next }
        np=split(pairs,P," "); width=w+2*np; rows++
        for(cat=0;cat<6;cat++){
          for(i=1;i<=ncount;i++) if(v[cat*width+i]+0 != 3*first[$1,cat*w+i]){ print stat, $1, "category", cat+1, "count", i, "is", v[cat*width+i], "not 3 x", first[$1,cat*w+i]; bad=1 }
          for(k=1;k<=np;k++){ split(P[k],q,":"); est=v[cat*width+q[1]]; lo=v[cat*width+q[2]]; hi=v[cat*width+q[2]+1]
            if(abs(lo-est)>1e-4 || abs(hi-est)>1e-4){ print stat, $1, "category", cat+1, "estimate", est, "interval", lo, hi; bad=1 } } } }
      END { if(nb!=3 || rows==0){ print stat, "unexpected output: blocks", nb, "rows", rows; bad=1 } exit bad }' "$1.out" >&2 || fail "$1 jackknife"
}
check dstat   5 4  "5:6"
check paircoa 10 4 "5:11 8:13"
check dist    19 18 "19:20"

# fst is not implemented and must say so
if compute -p fst blocks.acf.gz > /dev/null; then fail "fst should be rejected"; fi
msg=$("$G" compute -p fst blocks.acf.gz 2>&1 || true)
grep -q "not implemented" <<<"$msg" || fail "fst message"

echo "compute tests passed"
