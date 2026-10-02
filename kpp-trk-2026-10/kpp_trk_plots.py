#!/usr/bin/env python3
"""Tracker and CRV results for the kpp_trk cosmic runs (124984, 124986, 124989) from
EventNtuple files made by scripts/kpp_trk_reco.sh.

    source /cvmfs/mu2e.opensciencegrid.org/bin/pyenv.sh ana
    python kpp_trk_plots.py $OUT/run124984/nts.root -o plots/run124984

Writes to the output directory:
  summary.json            event counts, events with entries per branch, CRV per-sector
                          timing and residuals, T(CRV)-T(trk)
  rates_vs_subrun.png     fraction of events with a good track / a CRV coincidence, per subrun
  slide04_timing.png      T(CRV)-T(trk) for in-time track-CRV pairs, all sectors
  crv_timing_by_sector.png  T(CRV)-T(trk) per sector, two-Gaussian fits
  slide10_crv.png         track-CRV residuals across the bars

The track and CRV analysis is the one of ../kpp-cosmics-2026-10/kpp_cosmic_plots.py
(DocDB 58468, slides 4 and 10), imported from there so both notes use the same code.
These runs have no calorimeter data, so its calo sections do not apply.
"""
import argparse
import json
import os
import sys

import awkward as ak
import numpy as np
import uproot

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "kpp-cosmics-2026-10"))
import kpp_cosmic_plots as K  # noqa: E402

plt = K.plt
COUNTED = ["trk", "crvcoincs", "crvpulses", "crvdigis", "calohits", "caloclusters"]


def counts(files):
    """Events with at least one entry, and total entries, per top-level branch."""
    t0 = uproot.open(f"{files[0]}:{K.TREE}")
    leaf = {b: next((k for k in t0.keys() if k.startswith(b + "/")), None) for b in COUNTED}
    br = [v for v in leaf.values() if v] + ["evtinfo/run", "evtinfo/subrun", "trk/trk.nactive"]
    a = uproot.concatenate([f"{f}:{K.TREE}" for f in files], sorted(set(br)), library="ak")
    out = {}
    for b, k in leaf.items():
        if k is None:
            out[b] = None
            continue
        n = ak.num(a[k], axis=1)
        out[b] = {"events_with_entries": int(ak.sum(n > 0)), "entries": int(ak.sum(n))}
    good = ak.sum(a["trk/trk.nactive"] >= K.MIN_ACTIVE, axis=1)
    out["trk_nactive_ge_%d" % K.MIN_ACTIVE] = {"events_with_entries": int(ak.sum(good > 0)),
                                                "entries": int(ak.sum(good))}
    return a, out, good


def rates_vs_subrun(a, good, S, path, tag):
    run = K.flat(a["evtinfo/run"]); sub = K.flat(a["evtinfo/subrun"])
    crv = ak.num(a[next(k for k in a.fields if k.startswith("crvcoincs/"))], axis=1)
    has_trk, has_crv = np.asarray(good > 0), np.asarray(crv > 0)
    fig, ax = plt.subplots(figsize=(8, 4))
    S["rates_vs_subrun"] = {}
    for r in np.unique(run):
        q = run == r
        subs = np.unique(sub[q])
        ft = np.array([has_trk[q & (sub == s)].mean() for s in subs])
        fc = np.array([has_crv[q & (sub == s)].mean() for s in subs])
        ax.plot(subs, ft, ".", label=f"run {r}: good track")
        ax.plot(subs, fc, ".", label=f"run {r}: CRV coincidence")
        low = subs[fc < 0.5 * np.median(fc)]  # subruns where the CRV data stop
        S["rates_vs_subrun"][int(r)] = {"subruns": [int(subs.min()), int(subs.max())], "n_subruns": int(len(subs)),
                                         "track_fraction_min_max": [float(ft.min()), float(ft.max())],
                                         "crv_fraction_median": float(np.median(fc)),
                                         "crv_fraction_min_max": [float(fc.min()), float(fc.max())],
                                         "subruns_with_crv_below_half_median": [int(x) for x in low]}
    ax.set(xlabel="subrun", ylabel="fraction of events", ylim=(0, 1), title=f"{tag}  per-subrun rates")
    ax.legend(fontsize=8)
    K.save(fig, path)


def two_gauss(x, n1, m1, s1, n2, d, s2, c):
    """Two Gaussians plus a constant; the second peak sits d ns after the first."""
    return (n1 * np.exp(-0.5 * ((x - m1) / s1) ** 2) + n2 * np.exp(-0.5 * ((x - m1 - d) / s2) ** 2) + c)


def fit_two(ax, v, rng, bins, title, xlabel):
    """Fit two narrow peaks (sigma 1.5-10 ns) 12-35 ns apart, starting from the tallest bin taken
    once as the early and once as the late peak; keep the better chi2."""
    h, e = np.histogram(v, bins=bins, range=rng)
    xc = 0.5 * (e[1:] + e[:-1])
    w = np.sqrt(np.maximum(h, 1))
    mode = xc[h.argmax()]
    lo = [0, rng[0], 1.5, 0, 12, 1.5, 0]
    hi = [np.inf, rng[1], 10, np.inf, 35, 10, np.inf]
    best = None
    for m1 in (mode, mode - 22):
        p0 = [h.max(), min(max(m1, rng[0] + 1), rng[1] - 1), 4, h.max() / 2, 22, 4, 0]
        try:
            p, cov = K.curve_fit(two_gauss, xc, h, p0=p0, sigma=w, bounds=(lo, hi), maxfev=20000)
        except RuntimeError:
            continue
        chi2 = float(np.sum(((two_gauss(xc, *p) - h) / w) ** 2))
        if best is None or chi2 < best[2]:
            best = (p, cov, chi2)
    p, cov, chi2 = best
    err = np.sqrt(np.diag(cov))
    a1, a2 = p[0] * p[2], p[3] * p[5]
    r = {"n": int(len(v)), "peak1_ns": float(p[1]), "peak1_err": float(err[1]), "sigma1_ns": float(p[2]),
         "peak2_ns": float(p[1] + p[4]), "separation_ns": float(p[4]), "separation_err": float(err[4]),
         "sigma2_ns": float(p[5]), "fraction_peak2": float(a2 / (a1 + a2)), "chi2": chi2, "ndf": int(len(h) - 7)}
    ax.hist(v, bins=bins, range=rng, histtype="stepfilled", alpha=0.3)
    ax.hist(v, bins=bins, range=rng, histtype="step", color="C0")
    xf = np.linspace(*rng, 400)
    ax.plot(xf, two_gauss(xf, *p), "r-")
    ax.text(0.97, 0.95, f"N = {r['n']}\npeaks {r['peak1_ns']:.1f}, {r['peak2_ns']:.1f} ns\n"
            f"sigmas {r['sigma1_ns']:.1f}, {r['sigma2_ns']:.1f} ns\nlate fraction {r['fraction_peak2']:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=8,
            bbox=dict(boxstyle="square", fc="white", ec="0.7"))
    ax.set(title=title, xlabel=xlabel)
    return r


def timing_by_sector(V, vmatch, S, path, tag):
    """T(CRV)-T(trk) per sector, fitted with two Gaussians: in these runs every sector
    shows two peaks about 22 ns apart. (The sum over sectors mixes four peak pairs.)"""
    secs = sorted(int(x) for x in np.unique(V["sector"][vmatch]))
    fig, ax = plt.subplots(1, len(secs), figsize=(4 * len(secs), 3.6), squeeze=False)
    S["crv_timing_two_peaks"] = {}
    for j, sec in enumerate(secs):
        v = V["dt"][vmatch & (V["sector"] == sec)]
        c0 = S["crv"]["sectors"][sec]["dt_peak_ns"]
        c0 = float(np.median(v[np.abs(v - c0) < 50]))
        S["crv_timing_two_peaks"][sec] = fit_two(ax[0, j], v, (c0 - 50, c0 + 50), 50,
                                                 f"{tag}  sector {sec}", "T(CRV) - T(trk), ns")
    K.save(fig, path)


def late_fraction_by_subrun(V, vmatch, subrun, S):
    """Fraction of in-time pairs in the later of the two peaks, per subrun: is the split a
    property of the run, or does it drift?"""
    late = np.full(len(V["dt"]), -1)
    for sec, p in S["crv_timing_two_peaks"].items():
        q = vmatch & (V["sector"] == int(sec))
        late[q] = (V["dt"][q] > p["peak1_ns"] + 0.5 * p["separation_ns"]).astype(int)
    m = late >= 0
    sb, L = subrun[V["evt"][m]], late[m]
    subs = np.unique(sb)
    fr = np.array([L[sb == x].mean() for x in subs])
    n = np.array([(sb == x).sum() for x in subs])
    S["late_fraction_by_subrun"] = {"subruns": int(len(subs)), "mean": float(fr.mean()),
                                    "rms": float(fr.std()),
                                    "rms_expected_from_statistics": float(np.sqrt(fr.mean() * (1 - fr.mean()) / n).mean()),
                                    "min_max": [float(fr.min()), float(fr.max())]}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("-o", "--out", default="plots")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    a_cnt, cnt, good = counts(args.files)
    runs = sorted(set(K.flat(a_cnt["evtinfo/run"]).tolist()))
    tag = "run " + ",".join(map(str, runs))
    S = {"files": [os.path.basename(f) for f in args.files], "runs": runs, "events": len(a_cnt),
         "counts": cnt, "skipped": []}
    rates_vs_subrun(a_cnt, good, S, f"{args.out}/rates_vs_subrun.png", tag)

    a = K.load(args.files, ["trk", "crv"])
    trk = K.tracks(a)
    S["good_tracks"] = int(ak.sum(ak.num(trk)))
    cv = ak.zip({"sector": a["crvcoincs.sectorType"], "t": a["crvcoincs.time"],
                 **{c.lower(): a[f"crvcoincs.pos.fCoordinates.f{c}"] for c in "XYZ"}})
    V, vmatch = K.crv_sectors(trk, cv, S, args.out, tag)
    K.timing(None, None, V, vmatch, S, args.out, tag)
    timing_by_sector(V, vmatch, S, f"{args.out}/crv_timing_by_sector.png", tag)
    late_fraction_by_subrun(V, vmatch, K.flat(a_cnt["evtinfo/subrun"]), S)
    S["skipped"].append("calorimeter sections of kpp_cosmic_plots.py: no calorimeter data in kpp_trk")

    with open(f"{args.out}/summary.json", "w") as f:
        json.dump(S, f, indent=1, default=float)
    print(json.dumps(S, indent=1, default=float))


if __name__ == "__main__":
    main()
