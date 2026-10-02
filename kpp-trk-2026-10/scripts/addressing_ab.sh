#!/bin/bash
# Tracker decoder A/B on one raw file: addressing byLink (PassN default) vs byMnid.
# Decodes every event twice, hashes every StrawDigi (straw ID, TDCs, TOTs, PMP) per event,
# and compares. A control run uses byMnid with two panel MnIDs swapped in the map; its
# hashes must differ, which shows the comparison is sensitive to the panel assignment.
#
#   addressing_ab.sh <muse-workdir> <raw.art> <outdir> [nchunks=12]
#
# Writes <outdir>/addressing_ab.json. Exits non-zero if any job fails.
[ $# -ge 3 ] || { echo "usage: addressing_ab.sh <muse-workdir> <raw.art> <outdir> [nchunks=12]"; exit 2; }
RAW=$2; OUT=$(realpath -m "$3"); K=${4:-12}
HERE=$(dirname "$(realpath "$0")")
source "$HERE/env.sh" "$1" "$OUT/cfg"
mkdir -p "$OUT" && cd "$OUT"

# control map: swap the MnIDs of the first two rows of the last IoV block
MAP=$WORK/Offline/TrackerConditions/data/TrkPanelMap.txt
last=$(grep -n '^TABLE TrkPanelMap' "$MAP" | tail -1 | cut -d: -f1)
r1=$(awk -v s=$last 'NR>s && !/^#/ && NF>=7 {print NR; exit}' "$MAP")
r2=$(awk -v s=$r1 'NR>s && !/^#/ && NF>=7 {print NR; exit}' "$MAP")
awk -F, -v OFS=, -v a=$r1 -v b=$r2 'NR==FNR {if (FNR==a) ma=$1; if (FNR==b) mb=$1; next}
     {if (FNR==a) $1=mb; else if (FNR==b) $1=ma; print}' "$MAP" "$MAP" > "$CFG/kpp-trk-2026-10/TrkPanelMap_swapped.txt"

NEV=$(mu2e -c Offline/Print/fcl/count.fcl -s "$RAW" 2>&1 | awk '/Events total/ {print $5}')
[ -n "$NEV" ] || { echo "could not count events in $RAW"; exit 1; }
PER=$(( (NEV + K - 1) / K ))
echo "$RAW: $NEV events, $K chunks of $PER per mode"

run() {  # run <fcl> <dir> <nskip> <nevents>
  mkdir -p "$2" && cd "$2"
  mu2e -c "$1" -s "$RAW" --nskip $3 -n $4 2> err.log |
    tee >(grep -aE '^(ERROR|WARNING)' > decoder_msgs.txt) |
    python3 "$HERE/digi_digest.py" digest events.txt > digest.log
  echo "EXIT ${PIPESTATUS[0]}" >> digest.log
}
for m in bylink bymnid; do
  for i in $(seq 0 $((K - 1))); do
    ( run decode_$m.fcl "$OUT/$m/c$(printf %02d $i)" $((i * PER)) $PER ) &
  done
done
( run decode_bymnid_swapped.fcl "$OUT/control" 0 500 ) &
wait

bad=$(grep -L '^EXIT 0' "$OUT"/*/digest.log "$OUT"/*/c*/digest.log)
[ -z "$bad" ] || { echo "failed jobs: $bad"; exit 1; }
python3 "$HERE/digi_digest.py" compare "$OUT" "$RAW" > "$OUT/addressing_ab.json"
cat "$OUT/addressing_ab.json"
