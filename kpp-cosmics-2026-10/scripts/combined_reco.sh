#!/bin/bash
# Combined tracker+calo+CRV reco and EventNtuple for one raw KPP file, split into
# parallel event-range chunks, then merged.
#
#   combined_reco.sh <muse-workdir> <raw.art> <outdir> [nchunks=12]
#
# Writes <outdir>/cNN/{rec.art,nts.root,crvdqm.root,*.log}, then
# <outdir>/nts.root and <outdir>/crvdqm.root. Exits non-zero if any job fails.
# no "set -u": setupmu2e-art.sh and muse reference unset variables
[ $# -ge 3 ] || { echo "usage: combined_reco.sh <muse-workdir> <raw.art> <outdir> [nchunks=12]"; exit 2; }
WORK=$(realpath "$1"); RAW=$2; OUT=$(realpath -m "$3"); K=${4:-12}
NOTE=$(dirname "$(realpath "$0")")/..
NTSFCL=$(realpath "$NOTE/fcl/nts_combined_crvpulses.fcl")

source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh > /dev/null 2>&1
cd "$WORK" && muse setup > /dev/null 2>&1 || { echo "muse setup failed in $WORK"; exit 1; }
mkdir -p "$OUT" && cd "$OUT"
NEV=$(mu2e -c Offline/Print/fcl/count.fcl -s "$RAW" 2>&1 | awk '/Events total/ {print $5}')
[ -n "$NEV" ] || { echo "could not count events in $RAW"; exit 1; }
PER=$(( (NEV + K - 1) / K ))
echo "$RAW: $NEV events, $K chunks of $PER"

for i in $(seq 0 $((K - 1))); do
  d=$OUT/c$(printf %02d $i); mkdir -p "$d"
  ( cd "$d"
    mu2e -c PassN/Combined_Pass1.fcl -s "$RAW" --nskip $((i * PER)) -n $PER \
         -o rec.art -T crvdqm.root > reco.log 2>&1 &&
    mu2e -c "$NTSFCL" -s rec.art -T nts.root > nts.log 2>&1
    echo "EXIT $?" >> nts.log ) &
done
wait

bad=$(grep -L '^EXIT 0' "$OUT"/c*/nts.log)
[ -z "$bad" ] || { echo "failed chunks: $bad"; exit 1; }
hadd -f "$OUT/nts.root" "$OUT"/c*/nts.root > "$OUT/hadd.log" 2>&1 &&
hadd -f "$OUT/crvdqm.root" "$OUT"/c*/crvdqm.root >> "$OUT/hadd.log" 2>&1 || { echo "hadd failed"; exit 1; }
echo "done: $OUT/nts.root"
