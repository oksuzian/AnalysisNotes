#!/bin/bash
# EventNtuple for the production extracted-cosmic MC, every Nth mcs file.
#
#   nts_mc.sh <muse-workdir> <outdir> [njobs=6] [every=25]
#
# every=25 takes 100 of the 2500 files (about 2.9M events). Writes
# <outdir>/files.txt and <outdir>/fNNN/nts.root.
[ $# -ge 2 ] || { echo "usage: nts_mc.sh <muse-workdir> <outdir> [njobs=6] [every=25]"; exit 2; }
WORK=$(realpath "$1"); OUT=$(realpath -m "$2"); NJOBS=${3:-6}; EVERY=${4:-25}
HERE=$(dirname "$(realpath "$0")")
mkdir -p "$OUT"

source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh > /dev/null 2>&1
setup sam_web_client > /dev/null 2>&1
samweb list-files "dh.dataset mcs.mu2e.CosmicCRYExtracted.MDC2025au_best_v1_5.art" | sort |
  python3 "$HERE/filelist.py" --every "$EVERY" - > "$OUT/files.txt" || exit 1
echo "$(wc -l < "$OUT/files.txt") files"

"$HERE/run_files.sh" "$WORK" "$OUT/files.txt" "$OUT" "$NJOBS" nts_mc.fcl
