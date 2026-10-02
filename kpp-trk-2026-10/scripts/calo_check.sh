#!/bin/bash
# Does a raw file carry calorimeter data? Runs the calo decoder alone with its diagnostic
# on, and counts the events for which it finds no calorimeter subsystem data.
#
#   calo_check.sh <muse-workdir> <raw.art> <outdir> [nevents=100]
#
# Writes <outdir>/calo_check_<run>.json.
[ $# -ge 3 ] || { echo "usage: calo_check.sh <muse-workdir> <raw.art> <outdir> [nevents=100]"; exit 2; }
RAW=$2; OUT=$(realpath -m "$3"); N=${4:-100}
source "$(dirname "$(realpath "$0")")/env.sh" "$1" "$OUT/cfg"
mkdir -p "$OUT" && cd "$OUT"
run=$(basename "$RAW" | sed -E 's/^raw\.mu2e\.[^.]+\.[^.]+\.([0-9]+)_.*/\1/')
mu2e -c calo_check.fcl -s "$RAW" -n $N > calo_check_$run.log 2>&1 || { echo "job failed, see calo_check_$run.log"; exit 1; }
nev=$(awk '/TrigReport Events total/ {print $5}' calo_check_$run.log)
nocalo=$(grep -c 'found no Calorimeter decoders' calo_check_$run.log)
printf '{"file": "%s", "events": %s, "events_without_calo_data": %s}\n' "$(basename "$RAW")" "$nev" "$nocalo" | tee calo_check_$run.json
