# Sourced by the other scripts:  source env.sh <muse-workdir> <cfgdir>
# Sets up muse in the workdir, writes the stopgap tracker calibration file into
# <cfgdir>/kpp-trk-2026-10/, and puts the note's fcl/ and <cfgdir> on the search paths.
# no "set -u": setupmu2e-art.sh and muse reference unset variables
WORK=$(realpath "$1"); CFG=$(realpath -m "$2")
NOTE=$(realpath "$(dirname "${BASH_SOURCE[0]}")/..")
source /cvmfs/mu2e.opensciencegrid.org/setupmu2e-art.sh > /dev/null 2>&1
cd "$WORK" && muse setup > /dev/null 2>&1 || { echo "muse setup failed in $WORK"; exit 1; }
[ -f "$WORK/PassN/Combined_Pass1.fcl" ] || { echo "no PassN/Combined_Pass1.fcl in $WORK (clone PassN#21 there)"; exit 1; }
echo "PassN at $(git -C "$WORK/PassN" log -1 --format=%h), Offline at $(git -C "$WORK/Offline" log -1 --format=%h)"

mkdir -p "$CFG/kpp-trk-2026-10"
CAL=$CFG/kpp-trk-2026-10/trackercalibrations_kpp_trk.txt
sed 's/^\(TABLE [A-Za-z]*\) 123680-123680$/\1 124000-129999/' \
    "$WORK/PassN/Tracker/trackercalibrations.txt" > "$CAL"
[ "$(grep -c '^TABLE .* 124000-129999$' "$CAL")" -eq 3 ] ||
  { echo "expected 3 run-123680 tables to relabel in PassN/Tracker/trackercalibrations.txt"; exit 1; }

export FHICL_FILE_PATH=$NOTE/fcl:$FHICL_FILE_PATH
export MU2E_SEARCH_PATH=$MU2E_SEARCH_PATH:$CFG
