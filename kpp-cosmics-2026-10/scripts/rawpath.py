#!/usr/bin/env python3
"""Print the /pnfs path of Mu2e raw files from their names (tape-backed, read in place).

    python3 rawpath.py raw.mu2e.cosmics.kpp.123680_000001.art
    python3 rawpath.py --kpp-calo 123681 1 201 10 > calo_123681.txt   # subruns 1,11,...,201

The directory is the first two byte pairs of sha256(file name). Check that a file
is on disk, not only on tape, before reading it:
    cat "<dir>/.(get)(<name>)(locality)"     # want ONLINE or ONLINE_AND_NEARLINE
"""
import hashlib
import sys


def pnfs_path(name):
    tier, owner, desc, conf = name.split(".")[:4]
    h = hashlib.sha256(name.encode()).hexdigest()
    return f"/pnfs/mu2e/tape/phy-raw/{tier}/{owner}/{desc}/{conf}/art/{h[:2]}/{h[2:4]}/{name}"


if __name__ == "__main__":
    if sys.argv[1] == "--kpp-calo":
        run, first, last, step = map(int, sys.argv[2:6])
        for s in range(first, last + 1, step):
            print(pnfs_path(f"raw.mu2e.cosmics.kpp_calo.{run}_{s:06d}-01.art"))
    else:
        for n in sys.argv[1:]:
            print(pnfs_path(n))
