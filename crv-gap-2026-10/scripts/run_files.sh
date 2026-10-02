#!/bin/bash
# Run one art job chain per input file, NJOBS at a time, at low priority.
#
#   run_files.sh <muse-workdir> <filelist.txt> <outdir> <njobs> <fcl> [<ntuple fcl>]
#
# With two fcl files: reco (<fcl>, writes rec.art) then EventNtuple (<ntuple fcl>);
# rec.art is deleted once the ntuple is written. With one: <fcl> writes nts.root.
# fcl names are looked up in this note's fcl/ directory first. Set EXTRA_SEARCH to
# append a directory to MU2E_SEARCH_PATH (DbService text files).
#
# Writes <outdir>/fNNN/{nts.root,*.log,job.log}. A file whose job.log already says
# "EXIT 0" is skipped, so a rerun only redoes failures. Exits non-zero if any job failed.
# no "set -u": setupmu2e-art.sh and muse reference unset variables
[ $# -ge 5 ] || { echo "usage: run_files.sh <muse-workdir> <filelist> <outdir> <njobs> <fcl> [<ntuple fcl>]"; exit 2; }
WORK=$(realpath "$1"); LIST=$(realpath "$2"); OUT=$(realpath -m "$3"); NJOBS=$4
FCL=$5; FCL2=$6
NOTE=$(dirname "$(realpath "$0")")/..

source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh > /dev/null 2>&1
cd "$WORK" && muse setup > /dev/null 2>&1 || { echo "muse setup failed in $WORK"; exit 1; }
export FHICL_FILE_PATH=$(realpath "$NOTE/fcl"):$FHICL_FILE_PATH
[ -n "$EXTRA_SEARCH" ] && export MU2E_SEARCH_PATH=$MU2E_SEARCH_PATH:$EXTRA_SEARCH
mkdir -p "$OUT"

job() {
  local d=$OUT/f$(printf %03d "$1")
  grep -qs '^EXIT 0' "$d/job.log" && return 0
  rm -rf "$d"; mkdir -p "$d"; cd "$d" || return 1
  if [ -n "$FCL2" ]; then
    nice -n 10 mu2e -c "$FCL" -s "$2" -o rec.art > reco.log 2>&1 &&
    nice -n 10 mu2e -c "$FCL2" -s rec.art -T nts.root > nts.log 2>&1 &&
    rm -f rec.art
  else
    nice -n 10 mu2e -c "$FCL" -s "$2" -T nts.root > nts.log 2>&1
  fi
  echo "EXIT $? $2" > job.log
}
export -f job; export OUT FCL FCL2

grep -v '^\s*$' "$LIST" | nl -v0 -w1 -s' ' | xargs -P "$NJOBS" -n 2 bash -c 'job "$@"' _

n=$(grep -cv '^\s*$' "$LIST")
bad=$(for i in $(seq 0 $((n - 1))); do f=$OUT/f$(printf %03d $i)/job.log; grep -qs '^EXIT 0' "$f" || echo "$f"; done)
[ -z "$bad" ] || { echo "failed jobs:"; echo "$bad"; exit 1; }
echo "done: $n files in $OUT"
