#!/usr/bin/env python3
"""Per-event StrawDigi digests for scripts/addressing_ab.sh.

    ... | digi_digest.py digest events.txt     # PrintModule output on stdin
    digi_digest.py compare <outdir> <raw.art>  # JSON summary of an addressing_ab.sh run

digest: reads PrintModule output (strawDigiPrinter, verbose 1) and writes one line per
event, "run subrun event ndigi md5", where the md5 is over the digis in decoder order
(straw ID, TDC0, TDC1, TOT0, TOT1, PMP).
"""
import glob
import hashlib
import json
import re
import sys
from collections import Counter

EVENT = re.compile(rb"PrintModule Run/Subrun/Event\s+(\d+)\s+(\d+)\s+(\d+)")
DIGI = re.compile(rb"^\s*\d+\s+(\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s+(-?\d+)\s*$")


def digest(path):
    key, h, n = None, None, 0
    with open(path, "w") as out:
        def flush():
            if key is not None:
                out.write(f"{key} {n} {h.hexdigest()}\n")
        for line in sys.stdin.buffer:
            m = EVENT.search(line)
            if m:
                flush()
                key, h, n = " ".join(x.decode() for x in m.groups()), hashlib.md5(), 0
                continue
            d = DIGI.match(line)
            if d and key is not None:
                h.update(b" ".join(d.groups()) + b"\n")
                n += 1
        flush()


def read(files):
    ev = {}
    for f in sorted(files):
        for line in open(f):
            r, s, e, n, md5 = line.split()
            ev[(int(r), int(s), int(e))] = (int(n), md5)
    return ev


def compare(outdir, raw):
    a = read(glob.glob(f"{outdir}/bylink/c*/events.txt"))
    b = read(glob.glob(f"{outdir}/bymnid/c*/events.txt"))
    c = read(glob.glob(f"{outdir}/control/events.txt"))
    msgs = {m: sum(1 for f in glob.glob(f"{outdir}/{m}/c*/decoder_msgs.txt") for _ in open(f))
            for m in ("bylink", "bymnid")}
    same = [k for k in a if k in b and a[k] == b[k]]
    ctl_diff = [k for k in c if k in a and c[k] != a[k]]
    S = {
        "file": raw.split("/")[-1],
        "events_bylink": len(a), "events_bymnid": len(b),
        "digis_bylink": sum(v[0] for v in a.values()),
        "digis_bymnid": sum(v[0] for v in b.values()),
        "events_identical": len(same),
        "events_differing": len(set(a) | set(b)) - len(same),
        "decoder_messages_bylink": msgs["bylink"], "decoder_messages_bymnid": msgs["bymnid"],
        "control_events": len(c), "control_events_differing_from_bylink": len(ctl_diff),
    }
    print(json.dumps(S, indent=1))


if __name__ == "__main__":
    if sys.argv[1] == "digest":
        digest(sys.argv[2])
    elif sys.argv[1] == "compare":
        compare(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)
