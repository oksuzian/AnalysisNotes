#!/usr/bin/env python3
"""Turn Mu2e file names into /pnfs paths, and check that every file is on disk.

    samweb list-files "run_number 124155 and data_tier raw" | sort | python3 filelist.py - > files.txt
    python3 filelist.py --every 25 - < names.txt          # every 25th name only

Files are tape-backed and read in place; they are never copied. The directory is
the first two byte pairs of sha256(file name). Exits non-zero, listing the files,
if any file is only on tape (reading it would stall on a tape recall).
"""
import argparse
import hashlib
import sys

AREA = {"raw": "phy-raw", "dig": "phy-sim", "mcs": "phy-sim", "nts": "phy-nts"}


def pnfs_path(name):
    tier, owner, desc, conf, _, ext = name.split(".")
    h = hashlib.sha256(name.encode()).hexdigest()
    return f"/pnfs/mu2e/tape/{AREA[tier]}/{tier}/{owner}/{desc}/{conf}/{ext}/{h[:2]}/{h[2:4]}/{name}"


def locality(path):
    d, name = path.rsplit("/", 1)
    with open(f"{d}/.(get)({name})(locality)") as f:
        return f.read().strip()


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("names", nargs="+", help="file names, or - to read them from stdin")
    ap.add_argument("--every", type=int, default=1, help="keep every Nth name")
    a = ap.parse_args()
    names = sys.stdin.read().split() if a.names == ["-"] else a.names
    paths = [pnfs_path(n) for n in names[:: a.every]]
    offline = [p for p in paths if not locality(p).startswith("ONLINE")]
    if offline:
        sys.exit("not on disk (tape only):\n" + "\n".join(offline))
    print("\n".join(paths))
