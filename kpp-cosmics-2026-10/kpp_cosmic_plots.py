#!/usr/bin/env python3
"""Reproduce the KPP cosmic-run plots of P. Murat, DocDB 58468, from EventNtuple files.

    source /cvmfs/mu2e.opensciencegrid.org/bin/pyenv.sh ana
    python kpp_cosmic_plots.py nts.combined.123680.root -o plots_123680   # tracker+calo+CRV
    python kpp_cosmic_plots.py calo_123681/f*/nts.root -o plots_123681    # calo only

Slides reproduced, each when the ntuple has the branches it needs:
     4  timing: T(cal)-T(trk), T(cal)-T(CRV) (triple coincidences), T(CRV)-T(trk)  [trk+calo, trk+crv]
   6-8  calorimeter disk position: track-cluster residuals vs track slope           [trk+calo]
    10  CRV position: track-CRV coincidence residuals                               [trk+crv]
 13-14  Michel decays: time between two hits in one crystal                         [calo]
    18  energy of the cluster holding the second (Michel) hit                       [calo]

Positions are in the tracker frame. The tracks are KinematicLine fits, and the
straight line is taken from the first trkseg (position, direction, time).
EventNtuple stores caloclusters.cog_ in the disk front-face frame. The offset
of each disk to the tracker frame is measured from single-crystal clusters,
because calohits.crystalPos_ is in the tracker frame.

Nothing is assumed about alignment or timing; all of it is measured:
  - time offsets: the peak of the raw time difference, per CRV sector;
  - disk Z: the plane shift that maximises track-cluster matches, refined by
    zeroing the residual-vs-slope (Pasha's method, slides 6-7).

Calorimeter data times need Offline#2022 (CaloDigi t0 in 5 ns ticks). Without it,
calo times come out 5x too small: positions still match tracks, but the time
pre-selection that the disk fit starts from finds no calo-track peak.
"""
import argparse
import json
import os

import awkward as ak
import matplotlib
import numpy as np
import uproot
from scipy.optimize import curve_fit

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

C_MM_PER_NS = 299.792458
TREE = "EventNtuple/ntuple"
MIN_ACTIVE = 10          # slide 4: tracks with N(active hits) >= 10
DT_WINDOW = 30.0         # slide 4: |dT| < 30 ns
CALO_DT_WINDOW = 1000.0  # calo-track pre-selection window around the peak, ns
MATCH_MM = 150.0         # track-cluster spatial match after the Z correction
MICHEL_EMIN = 10.0       # slide 13: hit EDep > 10 MeV
MICHEL_NSIPM = 2         # slide 13: 2-SiPM hits
OVERFLOW_NS = 100.0      # slides 14-15: dT below this is mostly overflow waveforms

GROUPS = {
    "trk": ["trk.nactive", "trksegs"],
    "calo": ["caloclusters.diskID_", "caloclusters.time_", "caloclusters.energyDep_",
             "caloclusters.cog_.fCoordinates.fX", "caloclusters.cog_.fCoordinates.fY",
             "caloclusters.cog_.fCoordinates.fZ", "caloclusters.size_", "caloclusters.hits_",
             "calohits.crystalId_", "calohits.nSiPMs_", "calohits.time_", "calohits.eDep_",
             "calohits.clusterIdx_", "calohits.crystalPos_.fCoordinates.fX",
             "calohits.crystalPos_.fCoordinates.fY", "calohits.crystalPos_.fCoordinates.fZ"],
    "crv": ["crvcoincs.sectorType", "crvcoincs.time", "crvcoincs.pos.fCoordinates.fX",
            "crvcoincs.pos.fCoordinates.fY", "crvcoincs.pos.fCoordinates.fZ"],
}


def available_groups(path):
    t = uproot.open(f"{path}:{TREE}")
    names = {k.split("/")[-1] for k in t.keys()}
    return [g for g, br in GROUPS.items() if all(b in names for b in br)]


def load(files, groups):
    # evtinfo leaves are stored as "evtinfo/run", which an expression cannot name directly
    br = ["evtrun"] + [b for g in groups for b in GROUPS[g]]
    return uproot.concatenate([f"{f}:{TREE}" for f in files], br,
                              aliases={"evtrun": "evtinfo/run"}, library="ak")


def flat(x):
    return np.asarray(ak.to_numpy(ak.flatten(x, axis=None)))


def tracks(a):
    """One straight line per good track: point, unit direction, time at the point."""
    segs = a["trksegs"]
    good = (ak.num(segs, axis=2) > 0) & (a["trk.nactive"] >= MIN_ACTIVE)
    s0 = ak.firsts(segs[good], axis=2)
    p, m = s0["pos"]["fCoordinates"], s0["mom"]["fCoordinates"]
    pm = np.sqrt(m.fX**2 + m.fY**2 + m.fZ**2)
    return ak.zip({"x": p.fX, "y": p.fY, "z": p.fZ, "t": s0["time"],
                   "dx": m.fX / pm, "dy": m.fY / pm, "dz": m.fZ / pm})


def clusters(a):
    """Clusters in the tracker frame. The front-face-to-tracker offset per disk comes
    from single-crystal clusters, whose cog is the crystal position."""
    one = a["caloclusters.size_"] == 1
    hidx = ak.firsts(a["caloclusters.hits_"][one], axis=2)
    disk1 = flat(a["caloclusters.diskID_"][one])
    off = {0: [0.0] * 3, 1: [0.0] * 3}
    for i, c in enumerate("XYZ"):
        d = flat(a[f"calohits.crystalPos_.fCoordinates.f{c}"][hidx]
                 - a[f"caloclusters.cog_.fCoordinates.f{c}"][one])
        for k in (0, 1):
            off[k][i] = float(np.median(d[disk1 == k]))
            spread = float(np.percentile(np.abs(d[disk1 == k] - off[k][i]), 99))
            if spread > 1.0:
                raise RuntimeError(f"disk {k} axis {c}: offset is not a pure translation "
                                   f"(99% spread {spread:.2f} mm)")
    d = a["caloclusters.diskID_"]
    cl = ak.zip({"disk": d, "t": a["caloclusters.time_"], "e": a["caloclusters.energyDep_"],
                 **{c.lower(): a[f"caloclusters.cog_.fCoordinates.f{c}"] + ak.where(d == 0, off[0][i], off[1][i])
                    for i, c in enumerate("XYZ")}})
    return cl, off


def to_plane(trk, axis, value):
    """Intersect the track line with the plane axis=value; returns x, y, z, t there."""
    s = (value - trk[axis]) / trk["d" + axis]
    return trk.x + s * trk.dx, trk.y + s * trk.dy, trk.z + s * trk.dz, trk.t + s / C_MM_PER_NS


def peak(dt, lo=-1.2e5, hi=1.2e5, width=50.0):
    """Peak of a raw time difference (mode, refined with the median within +-2 bins)."""
    h, e = np.histogram(dt, bins=max(1, int((hi - lo) / width)), range=(lo, hi))
    j = h.argmax()
    mode = 0.5 * (e[j] + e[j + 1])
    significance = h[j] / (np.median(h[h > 0]) if (h > 0).any() else 1.0)
    return float(np.median(dt[np.abs(dt - mode) < 2 * width])), float(significance)


def gauss(x, n, mu, sig, c):
    return n * np.exp(-0.5 * ((x - mu) / sig) ** 2) + c


def fit_gauss(ax, v, rng, bins, title, xlabel):
    h, e = np.histogram(v, bins=bins, range=rng)
    xc = 0.5 * (e[1:] + e[:-1])
    core = v[(v > rng[0]) & (v < rng[1])]
    med = np.median(core)
    p0 = [h.max(), med, 1.4826 * np.median(np.abs(core - med)) + 1e-3, 0.0]
    try:
        p, cov = curve_fit(gauss, xc, h, p0=p0, sigma=np.sqrt(np.maximum(h, 1)), maxfev=20000)
        err = np.sqrt(np.diag(cov))
    except RuntimeError:
        p, err = np.array(p0), np.full(4, np.nan)
    ax.stairs(h, e, color="tab:blue", fill=True, alpha=0.3)
    ax.stairs(h, e, color="tab:blue")
    xf = np.linspace(*rng, 400)
    ax.plot(xf, gauss(xf, *p), "r-", lw=1.2)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(xlabel)
    ax.text(0.97, 0.95, f"N = {len(core)}\nmean = {p[1]:.1f} ± {err[1]:.1f}\n"
            f"sigma = {abs(p[2]):.1f} ± {err[2]:.1f}", transform=ax.transAxes,
            ha="right", va="top", fontsize=8, bbox=dict(fc="white", ec="0.7"))
    return {"n": int(len(core)), "mean": float(p[1]), "sigma": float(abs(p[2])),
            "mean_err": float(err[1]), "sigma_err": float(err[2])}


def fit_slope(slope, resid, cut=200.0):
    """Clipped straight-line fit resid = a + b*slope."""
    m = np.abs(resid) < cut
    for _ in range(4):
        b, a0 = np.polyfit(slope[m], resid[m], 1)
        r = resid - (a0 + b * slope)
        m = np.abs(r) < max(3 * 1.4826 * np.median(np.abs(r[m])), 20.0)
    return float(a0), float(b)


def save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def calo_disks(trk, cl, S, out, tag):
    """Slides 6-8: disk Z from a plane scan, then zero the residual-vs-slope."""
    tc = ak.cartesian({"t": trk, "c": cl}, axis=1)
    X, Y, _, Tt = to_plane(tc.t, "z", tc.c.z)
    P = {"evt": flat(ak.local_index(tc, axis=0) + 0 * tc.c.disk), "disk": flat(tc.c.disk),
         "t_cal": flat(tc.c.t), "dt": flat(tc.c.t - Tt), "dX": flat(X - tc.c.x), "dY": flat(Y - tc.c.y),
         "sx": flat(tc.t.dx / tc.t.dz), "sy": flat(tc.t.dy / tc.t.dz)}
    off_tc, sig_tc = peak(P["dt"])
    S["calo_minus_trk_peak_ns"] = {"value": off_tc, "peak_over_median": sig_tc}
    inwin = np.abs(P["dt"] - off_tc) < CALO_DT_WINDOW
    S["disks"] = {}
    matched = np.zeros(len(P["dt"]), bool)
    for disk in (0, 1):
        q = inwin & (P["disk"] == disk)
        scan = [(int(((np.abs(P["dX"][q] + z * P["sx"][q]) < 60) &
                      (np.abs(P["dY"][q] + z * P["sy"][q]) < 60)).sum()), z) for z in range(-3000, 3001, 25)]
        n_best, z0 = max(scan)
        _, bx = fit_slope(P["sx"][q], P["dX"][q] + z0 * P["sx"][q])
        _, by = fit_slope(P["sy"][q], P["dY"][q] + z0 * P["sy"][q])
        dz = z0 - 0.5 * (bx + by)
        rx, ry = P["dX"] + dz * P["sx"], P["dY"] + dz * P["sy"]
        m = q & (np.abs(rx) < MATCH_MM) & (np.abs(ry) < MATCH_MM)
        matched |= m
        sel = q & (np.abs(rx) < 1000) & (np.abs(ry) < 1000)
        fig, ax = plt.subplots(2, 2, figsize=(10, 8))
        for j, (r, s, nm) in enumerate([(rx, P["sx"], "x"), (ry, P["sy"], "y")]):
            ax[0, j].scatter(s[sel], r[sel], s=2, alpha=0.4)
            ax[0, j].set(xlim=(-1, 1), ylim=(-1000, 1000), xlabel=f"track d{nm}/dz",
                         ylabel=f"Δ{nm.upper()} = {nm.upper()}(trk) - {nm.upper()}(cluster), mm",
                         title=f"{tag}  disk{disk}: d{nm}_calc vs d{nm}dz, Z shifted {dz:+.0f} mm")
        res = {c: fit_gauss(ax[1, j], r[sel], (-1000, 1000), 100, f"disk{disk}: Δ{c.upper()}", "residual, mm")
               for j, (c, r) in enumerate([("x", rx), ("y", ry)])}
        save(fig, f"{out}/" + ("slide06_07_disk0.png" if disk == 0 else "slide08_disk1.png"))
        S["disks"][disk] = {"z_shift_mm": dz, "scan_best_mm": z0, "scan_matches": n_best,
                            "scan_matches_at_nominal": next(n for n, z in scan if z == 0),
                            "slope_x_after_scan_mm": bx, "slope_y_after_scan_mm": by,
                            "matched_pairs": int(m.sum()), "resid_x": res["x"], "resid_y": res["y"]}
    return P, matched


def crv_sectors(trk, cv, S, out, tag):
    """Slide 10: track-CRV residuals, across the bars of the sectors that see tracks."""
    tv = ak.cartesian({"t": trk, "v": cv}, axis=1)
    VX, _, VZ, Tv = to_plane(tv.t, "y", tv.v.y)
    V = {"evt": flat(ak.local_index(tv, axis=0) + 0 * tv.v.sector), "sector": flat(tv.v.sector),
         "t_crv": flat(tv.v.t), "dt": flat(tv.v.t - Tv), "dX": flat(VX - tv.v.x), "dZ": flat(VZ - tv.v.z)}
    vmatch = np.zeros(len(V["dt"]), bool)
    S["crv"] = {"sectors": {}}
    for sec in np.unique(V["sector"]):
        q = V["sector"] == sec
        off_s, sig_s = peak(V["dt"][q])
        m = q & (np.abs(V["dt"] - off_s) < DT_WINDOW)
        info = {"pairs": int(q.sum()), "dt_peak_ns": off_s, "peak_over_median": sig_s, "in_time": int(m.sum())}
        if sig_s >= 10 and m.sum() >= 20:
            vmatch |= m
            wx = 1.4826 * np.median(np.abs(V["dX"][m] - np.median(V["dX"][m])))
            wz = 1.4826 * np.median(np.abs(V["dZ"][m] - np.median(V["dZ"][m])))
            info.update({"precise": "X" if wx < wz else "Z", "robust_width_x_mm": float(wx),
                         "robust_width_z_mm": float(wz)})
        else:
            info["precise"] = None  # no tracks reach this sector in time
        S["crv"]["sectors"][int(sec)] = info
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for j, c in enumerate("XZ"):
        secs = [s for s, i in S["crv"]["sectors"].items() if i["precise"] == c]
        m = vmatch & np.isin(V["sector"], secs)
        if not m.any():
            ax[j].text(0.5, 0.5, f"no sector with in-time tracks\nmeasures {c} across its bars",
                       ha="center", va="center", transform=ax[j].transAxes)
            ax[j].set_title(f"{tag}  d{c.lower()}_crvc")
            S["crv"]["fit_" + c] = None
            continue
        r = V["d" + c][m]
        h, e = np.histogram(r, bins=120, range=(-3000, 3000))
        mode = 0.5 * (e[h.argmax()] + e[h.argmax() + 1])
        S["crv"]["fit_" + c] = fit_gauss(ax[j], r, (mode - 1000, mode + 1000), 100,
                                         f"{tag}  d{c.lower()}_crvc, sectors {secs}",
                                         f"Δ{c} = {c}(trk) - {c}(CRV), mm")
    save(fig, f"{out}/slide10_crv.png")
    return V, vmatch


def timing(P, matched, V, vmatch, S, out, tag):
    """Slide 4. T(cal)-T(trk) for spatially matched pairs; T(cal)-T(CRV) for triple
    coincidences (one entry per event: first matched pair); T(CRV)-T(trk) as a bonus."""
    r = {}
    if P is not None:
        r["T(cal)-T(trk)"] = P["dt"][matched]
    if V is not None:
        r["T(CRV)-T(trk)"] = V["dt"][vmatch]
    if P is not None and V is not None:
        cal, crv = {}, {}
        for e_, t_ in zip(P["evt"][matched], P["t_cal"][matched]):
            cal.setdefault(e_, t_)
        for e_, t_ in zip(V["evt"][vmatch], V["t_crv"][vmatch]):
            crv.setdefault(e_, t_)
        triple = sorted(set(cal) & set(crv))
        S["triple_coincidence_events"] = len(triple)
        r["T(cal)-T(CRV)"] = np.array([cal[e_] - crv[e_] for e_ in triple])
    fig, ax = plt.subplots(1, len(r), figsize=(5 * len(r), 4), squeeze=False)
    S["timing"] = {}
    for k, (nm, v) in enumerate(r.items()):
        if len(v) < 10:
            ax[0, k].set_title(f"{tag}  {nm}: {len(v)} entries")
            S["timing"][nm] = {"n": int(len(v))}
            continue
        c0, _ = peak(v, v.min() - 1, v.max() + 1, 10.0)
        S["timing"][nm] = fit_gauss(ax[0, k], v, (c0 - 100, c0 + 100), 100, f"{tag}  {nm}, tallest peak", "ΔT, ns")
        S["timing"][nm]["p5_p95_ns"] = [float(np.percentile(v, 5)), float(np.percentile(v, 95))]
    save(fig, f"{out}/slide04_timing.png")
    if P is not None and matched.any():
        c = S["calo_minus_trk_peak_ns"]["value"]
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.hist(r["T(cal)-T(trk)"], bins=200, range=(c - 1000, c + 1000), histtype="step")
        ax.set(xlabel="T(cal) - T(trk), ns", title=f"{tag}  T(cal)-T(trk), matched pairs, full range")
        save(fig, f"{out}/slide04_timing_calo_wide.png")


def michel(a, S, out, tag):
    """Slides 13-14, 18: two hits in one crystal; dT between the 1st and 2nd, and the
    energy of the cluster holding the 2nd."""
    h = ak.zip({"cry": a["calohits.crystalId_"], "n": a["calohits.nSiPMs_"], "t": a["calohits.time_"],
                "e": a["calohits.eDep_"], "ci": a["calohits.clusterIdx_"]})
    h = h[(h.n == MICHEL_NSIPM) & (h.e > MICHEL_EMIN)]
    h = h[ak.argsort(h.t, axis=1)]
    h = h[ak.argsort(h.cry, axis=1, stable=True)]
    first = ak.concatenate([ak.ones_like(h.cry[:, :1], dtype=bool), h.cry[:, 1:] != h.cry[:, :-1]], axis=1)
    pair = (h.cry[:, 1:] == h.cry[:, :-1]) & first[:, :-1]
    second = h[:, 1:][pair]
    dt = flat(second.t - h[:, :-1][pair].t)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    hh, e = np.histogram(dt, bins=120, range=(0, 30000))
    xc = 0.5 * (e[1:] + e[:-1])
    fr = (xc > 1000) & (xc < 15000) & (hh > 0)
    (n0, tau), cov = curve_fit(lambda x, n, t: n * np.exp(-x / t), xc[fr], hh[fr],
                               p0=[hh[fr][0], 2200.0], sigma=np.sqrt(hh[fr]), maxfev=20000)
    ax[0].stairs(hh, e)
    ax[0].set_yscale("log")
    xf = np.linspace(1000, 15000, 200)
    ax[0].plot(xf, n0 * np.exp(-xf / tau), "r-")
    ax[0].set(xlabel="dT = T2 - T1, ns", title=f"{tag}: deltaT(pulse-pulse), same crystal")
    ax[0].text(0.95, 0.95, f"N = {len(dt)}\ntau = {tau:.0f} ± {np.sqrt(cov[1, 1]):.0f} ns",
               transform=ax[0].transAxes, ha="right", va="top", bbox=dict(fc="white", ec="0.7"))
    ax[1].hist(dt, bins=100, range=(0, 1000), histtype="step")
    ax[1].set_yscale("log")
    ax[1].set(xlabel="dT, ns", title=f"{tag}: deltaT, 0-1000 ns")
    save(fig, f"{out}/slide13_14_michel_dt.png")
    S["michel"] = {"pairs": int(len(dt)), "tau_ns": float(tau), "tau_err_ns": float(np.sqrt(cov[1, 1])),
                   "fit_range_ns": [1000, 15000], "pairs_dt_lt_100ns": int((dt < OVERFLOW_NS).sum())}

    ci = second.ci
    has = flat(ci >= 0)
    safe = ak.where(ci >= 0, ci, 0)
    eclu, dclu = flat(a["caloclusters.energyDep_"][safe])[has], flat(a["caloclusters.diskID_"][safe])[has]
    keep = dt[has] > OVERFLOW_NS
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for disk, col in ((0, "tab:blue"), (1, "tab:red")):
        ax[0].hist(eclu[keep & (dclu == disk)], bins=40, range=(0, 100), histtype="step", color=col, label=f"disk {disk}")
    ax[0].legend()
    ax[0].set(xlabel="E, MeV", title=f"{tag}: cluster energy, per disk (dT > {OVERFLOW_NS:.0f} ns)")
    ax[1].hist(eclu[keep], bins=40, range=(0, 100), histtype="step")
    ax[1].set(xlabel="E, MeV", title=f"{tag}: cluster energy, both disks")
    save(fig, f"{out}/slide18_michel_energy.png")
    S["michel"]["clusters_plotted"] = int(keep.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+")
    ap.add_argument("-o", "--out", default="plots")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    groups = available_groups(args.files[0])
    a = load(args.files, groups)
    runs = sorted(set(flat(a["evtrun"]).tolist()))
    tag = "run " + ",".join(map(str, runs))
    S = {"files": args.files, "runs": runs, "events": len(a), "branch_groups": groups, "skipped": []}

    trk = cl = cv = P = V = None
    if "trk" in groups:
        trk = tracks(a)
        S["good_tracks"] = int(ak.sum(ak.num(trk)))
    if "calo" in groups:
        cl, off = clusters(a)
        S["disk_frame_offset_mm"] = off
    if "crv" in groups:
        cv = ak.zip({"sector": a["crvcoincs.sectorType"], "t": a["crvcoincs.time"],
                     **{c.lower(): a[f"crvcoincs.pos.fCoordinates.f{c}"] for c in "XYZ"}})

    if trk is not None and cl is not None:
        P, matched = calo_disks(trk, cl, S, args.out, tag)
    else:
        S["skipped"].append("slides 6-8 and calo timing: need trk and calo branches")
    if trk is not None and cv is not None:
        V, vmatch = crv_sectors(trk, cv, S, args.out, tag)
    else:
        S["skipped"].append("slide 10 and CRV timing: need trk and crv branches")
    if P is not None or V is not None:
        timing(P, matched if P is not None else None, V, vmatch if V is not None else None, S, args.out, tag)
    else:
        S["skipped"].append("slide 4: needs trk plus calo or crv branches")
    if cl is not None:
        michel(a, S, args.out, tag)
    else:
        S["skipped"].append("slides 13-14, 18: need calo branches")

    with open(f"{args.out}/summary.json", "w") as f:
        json.dump(S, f, indent=1, default=float)
    print(json.dumps(S, indent=1, default=float))
    for s in S["skipped"]:
        print("SKIPPED:", s)


if __name__ == "__main__":
    main()
