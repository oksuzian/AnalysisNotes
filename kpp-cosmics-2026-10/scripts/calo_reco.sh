#!/bin/bash
# Calo-only reco and EventNtuple (PassN Calo/ fcl, P. Girotti) for a list of raw files,
# one job per file, then merged.
#
#   calo_reco.sh <muse-workdir> <filelist.txt> <outdir>
#
# Writes <outdir>/fNN/{caloReco.art,nts.root,*.log} and <outdir>/nts.root.
# Exits non-zero if any job fails.
# no "set -u": setupmu2e-art.sh and muse reference unset variables
[ $# -ge 3 ] || { echo "usage: calo_reco.sh <muse-workdir> <filelist.txt> <outdir>"; exit 2; }
WORK=$(realpath "$1"); LIST=$(realpath "$2"); OUT=$(realpath -m "$3")

source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh > /dev/null 2>&1
cd "$WORK" && muse setup > /dev/null 2>&1 || { echo "muse setup failed in $WORK"; exit 1; }
mkdir -p "$OUT"

i=0
while read -r RAW; do
  d=$OUT/f$(printf %02d $i); mkdir -p "$d"
  ( cd "$d"
    mu2e -c PassN/Calo/Calo_CRR_pass1.fcl -s "$RAW" -o caloReco.art > reco.log 2>&1 &&
    mu2e -c PassN/Calo/Calo_ntuple.fcl -s caloReco.art -T nts.root > nts.log 2>&1
    echo "EXIT $?" >> nts.log ) &
  i=$((i + 1))
done < "$LIST"
wait

bad=$(grep -L '^EXIT 0' "$OUT"/f*/nts.log)
[ -z "$bad" ] || { echo "failed jobs: $bad"; exit 1; }
hadd -f "$OUT/nts.root" "$OUT"/f*/nts.root > "$OUT/hadd.log" 2>&1 || { echo "hadd failed"; exit 1; }
echo "done: $OUT/nts.root"
