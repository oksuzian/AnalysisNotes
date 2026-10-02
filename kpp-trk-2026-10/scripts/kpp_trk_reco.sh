#!/bin/bash
# Combined reco (fcl/kpp_trk_pass1.fcl) and EventNtuple with CRV pulses for a list of
# raw kpp_trk files, one job per file, all in parallel, then merged.
#
#   kpp_trk_reco.sh <muse-workdir> <filelist.txt> <outdir>
#
# Writes <outdir>/fNN/{rec.art,nts.root,crvdqm.root,*.log}, then <outdir>/nts.root and
# <outdir>/crvdqm.root. Exits non-zero if a file is not on disk or any job fails.
[ $# -ge 3 ] || { echo "usage: kpp_trk_reco.sh <muse-workdir> <filelist.txt> <outdir>"; exit 2; }
LIST=$(realpath "$2"); OUT=$(realpath -m "$3")
source "$(dirname "$(realpath "$0")")/env.sh" "$1" "$OUT/cfg"
mkdir -p "$OUT"

while read -r RAW; do
  n=$(basename "$RAW")
  case $(cat "$(dirname "$RAW")/.(get)($n)(locality)" 2>&1) in
    ONLINE*) ;;
    *) echo "$n is not on disk (prestage it first)"; exit 1;;
  esac
done < "$LIST"

i=0
while read -r RAW; do
  d=$OUT/f$(printf %02d $i); mkdir -p "$d"
  ( cd "$d"
    echo "$RAW" > input.txt
    mu2e -c kpp_trk_pass1.fcl -s "$RAW" -o rec.art -T crvdqm.root > reco.log 2>&1 &&
    mu2e -c nts_combined_crvpulses.fcl -s rec.art -T nts.root > nts.log 2>&1
    echo "EXIT $?" >> nts.log ) &
  i=$((i + 1))
done < "$LIST"
wait

bad=$(grep -L '^EXIT 0' "$OUT"/f*/nts.log)
[ -z "$bad" ] || { echo "failed jobs: $bad"; exit 1; }
hadd -f "$OUT/nts.root" "$OUT"/f*/nts.root > "$OUT/hadd.log" 2>&1 &&
hadd -f "$OUT/crvdqm.root" "$OUT"/f*/crvdqm.root >> "$OUT/hadd.log" 2>&1 || { echo "hadd failed"; exit 1; }
awk '/TrigReport Events total/ {n += $5} /TrigReport +[0-9]+ +[0-9]+ .*KLFilter$/ && !s[FILENAME]++ {k += $4}
     END {print "events", n, "  passing KLFilter (a line track)", k}' "$OUT"/f*/reco.log
echo "done: $OUT/nts.root"
