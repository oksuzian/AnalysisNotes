#!/bin/bash
# Combined reco and EventNtuple for all of run 124155 (225 raw files, 10.6M events).
#
#   reco_124155.sh <muse-workdir> <outdir> [njobs=20] [every=1]
#
# every=N processes every Nth file only (a quick look). Writes <outdir>/files.txt,
# <outdir>/cfg/ (relabelled tracker tables) and <outdir>/fNNN/nts.root.
[ $# -ge 2 ] || { echo "usage: reco_124155.sh <muse-workdir> <outdir> [njobs=20] [every=1]"; exit 2; }
WORK=$(realpath "$1"); OUT=$(realpath -m "$2"); NJOBS=${3:-20}; EVERY=${4:-1}
HERE=$(dirname "$(realpath "$0")")
mkdir -p "$OUT/cfg"

source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh > /dev/null 2>&1
setup sam_web_client > /dev/null 2>&1
samweb list-files "run_number 124155 and data_tier raw" | sort |
  python3 "$HERE/filelist.py" --every "$EVERY" - > "$OUT/files.txt" || exit 1
echo "$(wc -l < "$OUT/files.txt") files"

# PassN#21 tracker tables: relabel the run-123680 tables to cover 124000-129999
TABLES=$WORK/PassN/Tracker/trackercalibrations.txt
[ -f "$TABLES" ] || { echo "missing $TABLES"; exit 1; }
sed -E 's/^(TABLE [A-Za-z]+) 123680-123680/\1 124000-129999/' "$TABLES" > "$OUT/cfg/trackercalibrations_124xxx.txt"
grep -q '124000-129999' "$OUT/cfg/trackercalibrations_124xxx.txt" || { echo "relabel failed"; exit 1; }

EXTRA_SEARCH=$OUT/cfg "$HERE/run_files.sh" "$WORK" "$OUT/files.txt" "$OUT" "$NJOBS" reco_124155.fcl nts_data.fcl
