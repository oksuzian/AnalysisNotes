#!/usr/bin/env python3
"""Print the /pnfs path of Mu2e raw files from their names (tape-backed, read in place).

    python3 rawpath.py $(grep 124984 ../data/files.txt) > files_124984.txt
    python3 rawpath.py --check <name> ...     # also print the dCache locality

The directory is the first two byte pairs of sha256(file name). A file must be on disk
(ONLINE or ONLINE_AND_NEARLINE), not only on tape, before it is read.
"""
import hashlib
import os
import sys


def pnfs_path(name):
    tier, owner, desc, conf = name.split(".")[:4]
    h = hashlib.sha256(name.encode()).hexdigest()
    return f"/pnfs/mu2e/tape/phy-raw/{tier}/{owner}/{desc}/{conf}/art/{h[:2]}/{h[2:4]}/{name}"


if __name__ == "__main__":
    check = sys.argv[1:2] == ["--check"]
    for n in sys.argv[2 if check else 1:]:
        p = pnfs_path(n)
        if check:
            with open(os.path.join(os.path.dirname(p), f".(get)({n})(locality)")) as f:
                print(f.read().strip(), p)
        else:
            print(p)
