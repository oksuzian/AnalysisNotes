#!/usr/bin/env python3
"""CRV module-gap inefficiency from tracks extrapolated to the EX sector, data vs MC.

Reproduces R. Mina, "Estimating CRV module gap inefficiency using extrapolated cosmic
tracks" (Sep 3 2026), from EventNtuple files:

    source /cvmfs/mu2e.opensciencegrid.org/bin/pyenv.sh ana
    python crv_gap.py --data 'run124155/f*/nts.root' --mc 'mc_au/f*/nts.root' -o plots

Method, in the order the code runs it:
  1. Geometry. The EX counter positions come from the ntuple itself (crvpulses.pos is
     the counter centre): 4 layers of 64 counters, stacked along z. Counter width,
     thickness and the nominal module gap are those of crv_counters_extracted_v04.txt.
  2. Tracks. KinematicLine fits; the straight line is the first trkseg (point,
     direction, time). Selections as on the slides: data status>0, nactive>=10,
     chi2/ndof<3, nplanes>=3; MC the same with nactive>=20 and fitcon>0.01.
  3. Offsets, measured on the first files of each sample: the CRV-track time offset
     (coincidences and pulses) and, for data only, the tracker z shift (peak of the EX
     coincidence z minus the extrapolated track z). Data tracks are moved by that shift.
  4. Per track and EX layer: the crossing point at the layer centre; its counter, bin
     type (counter, inter-dicounter gap, inter-module gap) and signed distance to the
     nearest module gap. The layer is "hit" if an in-time pulse sits in the counter the
     track points at or in either neighbour (as on the slides, used for the efficiency
     maps), or within 6 counters (the "wide" hit, used for the gap width). Multiple
     scattering between tracker and CRV moves the crossing by about 50 mm, so the narrow
     hit loses about 16% of tracks that the wide one keeps.
  5. MC truth: the muon's crossing of the EX inner face (crvcoincsmcplane_ex),
     carried straight to each layer, with the same hit rules, for every MC event.
  6. Effective gap width w(theta_z) = D/eps0, D = (1/rho) sum_d [eps0 N_den(d) - N_pass(d)]
     over |d| < GAP_WINDOW, eps0 the wide-hit efficiency in the plateau next to it, rho
     the mean denominator density. The integral does not depend on how the extrapolation
     smears the dip, as long as the window holds the smeared dip. The model
     max(0, W - t tan theta + 2 l sin theta) (a pulse needs a path of at least l in a
     counter) is fitted to MC truth for W and l, then to each track sample for W with l
     fixed. R = W(data) / W(MC track).

Writes the plots and summary.json to -o. Every number in the note comes from summary.json.
"""
import argparse
import glob
import json
import os
import re
from multiprocessing import Pool

import awkward as ak
import matplotlib
import numpy as np
import uproot

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

TREE = "EventNtuple/ntuple"
C_MM_PER_NS = 299.792458
EX_SECTOR, EX_TYPE = 0, 1            # crs.sectorNames index and crs.sectorTypeEX
COUNTER_WIDTH = 51.3                 # crs.scintillatorBarWidth, mm
LAYER_THICKNESS = 19.8               # crs.scintillatorBarThickness, mm
NOMINAL_GAP = 4.433                  # crs.gapBetweenModules, mm
EX_HALF_LENGTH = 3000.0              # crs.scintillatorBarLengthEX / 2, mm
X_MARGIN = 100.0                     # stay this far from the counter ends
EVENT_WINDOW_NS = 1.0e5              # data live time per event (track t0 spans 0-100 us)
IN_TIME_NS = 30.0                    # |dt - peak| for an in-time coincidence
PULSE_WINDOW = (-50.0, 100.0)        # pulse dt - peak; pulses read at the far bar end come ~50 ns late
COINC_DZ_MM = 250.0                  # coincidence-track |dz| for slide 7
NARROW, WIDE = 1, 6                  # hit = pulse within this many counters of the crossing
# Multiple scattering between the tracker and the CRV moves the extrapolated crossing by
# about 50 mm (MC track minus truth; 12 mm above 20 GeV, 230 mm below 1 GeV). So the gap
# width uses the wide hit (+-6 counters, about +-310 mm), integrates the deficit over
# +-300 mm and takes the plateau from 300-400 mm. Module gaps are 826 mm apart, so no
# crossing is more than 413 mm from its nearest gap.
GAP_WINDOW, PLATEAU = 300.0, (300.0, 400.0)   # mm, for the gap width
FLAT_Z = (-2700.0, -700.0)           # the flat region of slides 16-17
SEL = {"data": dict(nactive=10, fitcon=None), "mc": dict(nactive=20, fitcon=0.01)}
TYPES = ("counter", "dicounter gap", "module gap")
COLORS = {"data": "navy", "mc": "crimson", "truth": "darkorange", "mcw": "seagreen"}
NOFFSET_FILES = 10

LEAVES = ["trk.status", "trk.nactive", "trk.chisq", "trk.ndof", "trk.nplanes", "trk.fitcon",
          "trksegs", "crvpulses.barId", "crvpulses.sectorId", "crvpulses.time", "crvpulses.PEs",
          "crvpulses.pos.fCoordinates.fX", "crvpulses.pos.fCoordinates.fY",
          "crvpulses.pos.fCoordinates.fZ", "crvcoincs.sectorType", "crvcoincs.time",
          "crvcoincs.pos.fCoordinates.fX", "crvcoincs.pos.fCoordinates.fY",
          "crvcoincs.pos.fCoordinates.fZ"]
TRUTH = ["crvcoincsmcplane_ex.pdgId", "crvcoincsmcplane_ex.time", "crvcoincsmcplane_ex.kineticEnergy"] + [
    f"crvcoincsmcplane_ex.{v}.fCoordinates.f{c}" for v in ("pos", "dir") for c in "XYZ"]


# ---------------------------------------------------------------------------- geometry
def ex_geometry(files):
    """EX counters per layer from pulse positions: centres, edges, gaps, bar lookups."""
    seen = {}
    for f in files:
        a = uproot.open(f)[TREE].arrays([l for l in LEAVES if l.startswith("crvpulses.")], library="np")
        sec = np.concatenate(a["crvpulses.sectorId"]) if len(a["crvpulses.sectorId"]) else np.array([])
        if not len(sec):
            continue
        bar = np.concatenate(a["crvpulses.barId"])
        xyz = [np.concatenate(a[f"crvpulses.pos.fCoordinates.f{c}"]) for c in "XYZ"]
        for b, s, x, y, z in zip(bar, sec, *xyz):
            if s == EX_SECTOR:
                seen.setdefault(int(b), (float(x), float(y), float(z)))
        if len(seen) == 256:
            break
    if len(seen) != 256:
        raise RuntimeError(f"found {len(seen)} of 256 EX counters in the pulses")
    bars = np.array(sorted(seen))
    pos = np.array([seen[b] for b in bars])
    ys = np.unique(np.round(pos[:, 1], 2))
    if len(ys) != 4 or np.any(np.abs(pos[:, 0]) > 1e-3):
        raise RuntimeError(f"unexpected EX layout: layer y {ys}, |x| max {np.abs(pos[:, 0]).max()}")
    nbar = int(bars.max()) + 1
    g = {"y": ys, "layer_of_bar": np.full(nbar, -1), "counter_of_bar": np.full(nbar, -1), "layers": []}
    for L, yl in enumerate(ys):
        m = np.round(pos[:, 1], 2) == yl
        order = np.argsort(pos[m, 2])
        zc = pos[m, 2][order]
        if len(zc) != 64:
            raise RuntimeError(f"layer {L}: {len(zc)} counters")
        g["layer_of_bar"][bars[m][order]] = L
        g["counter_of_bar"][bars[m][order]] = np.arange(64)
        lo, hi = zc - COUNTER_WIDTH / 2, zc + COUNTER_WIDTH / 2
        width = lo[1:] - hi[:-1]                     # gap after counter k
        mod = np.where(width > 1.0)[0]               # module gaps: 4.433 mm
        dic = np.where((width > 0.05) & (width <= 1.0))[0]   # inter-dicounter gaps: 0.2 mm
        if len(mod) != 3 or not np.allclose(width[mod], NOMINAL_GAP, atol=0.01):
            raise RuntimeError(f"layer {L}: module gaps {width[mod]}")
        g["layers"].append({"zc": zc, "lo": lo, "hi": hi, "width": width, "mod": mod, "dic": dic,
                            "gap_centre": 0.5 * (hi[mod] + lo[mod + 1])})
    # the box all four layers cover, for the edge distance of slide 7
    g["zbox"] = (max(l["lo"][0] for l in g["layers"]), min(l["hi"][-1] for l in g["layers"]))
    return g


def classify(g, L, z):
    """Counter index the track points at, bin type (0,1,2 or -1 outside), bin index,
    signed distance to the nearest module gap centre."""
    lay = g["layers"][L]
    k = np.clip(np.searchsorted(lay["zc"], z), 1, 63)
    k = np.where(np.abs(z - lay["zc"][k - 1]) < np.abs(z - lay["zc"][k]), k - 1, k)
    inside = (z >= lay["lo"][k]) & (z <= lay["hi"][k])
    # if not inside the nearest counter, z sits in the gap on its left or right
    left = np.clip(np.where(z < lay["lo"][k], k - 1, k), 0, 62)
    typ = np.where(inside, 0, np.where(lay["width"][left] > 1.0, 2, 1))
    out = (z < lay["lo"][0]) | (z > lay["hi"][-1])
    typ = np.where(out, -1, typ)
    idx = np.where(typ == 0, k // 2, left)          # dicounter number, or gap after counter
    d = z[:, None] - lay["gap_centre"][None, :]
    j = np.abs(d).argmin(axis=1)
    return k, typ, idx, d[np.arange(len(z)), j], j


# ------------------------------------------------------------------------ per-file pass
def line_at(tr, y):
    """Straight line at height y: x, z and time there."""
    s = (y - tr["y"]) / tr["dy"]
    return tr["x"] + s * tr["dx"], tr["z"] + s * tr["dz"], tr["t"] + s / C_MM_PER_NS


def angles(dx, dy, dz):
    return {"zenith": np.degrees(np.arccos(np.abs(dy))),
            "thz": np.degrees(np.arctan2(np.abs(dz), np.abs(dy))),
            "phi": np.degrees(np.arctan2(dz, dx))}


def process(job):
    """One ntuple file -> flat numpy arrays. job = (file, kind, geometry, offsets)."""
    f, kind, g, off = job
    t = uproot.open(f)[TREE]
    a = t.arrays(LEAVES, library="ak")
    out = {"nevents": len(a)}
    # every track (for the census), then the analysis selection
    st, na = a["trk.status"], a["trk.nactive"]
    c2 = a["trk.chisq"] / ak.where(a["trk.ndof"] > 0, a["trk.ndof"], 1)
    s0 = ak.firsts(a["trksegs"], axis=2)
    ok = ~ak.is_none(s0, axis=1)
    sel = SEL[kind]
    qual = (st > 0) & (na >= sel["nactive"]) & (c2 < 3)
    if sel["fitcon"] is not None:
        qual = qual & (a["trk.fitcon"] > sel["fitcon"])
    ana = qual & (a["trk.nplanes"] >= 3) & ok
    m = s0["mom"]["fCoordinates"]
    pm = np.sqrt(m.fX ** 2 + m.fY ** 2 + m.fZ ** 2)
    p = s0["pos"]["fCoordinates"]
    trk = {"x": p.fX, "y": p.fY, "z": p.fZ + off["zshift"], "t": s0["time"],
           "dx": m.fX / pm, "dy": m.fY / pm, "dz": m.fZ / pm}
    fl = lambda v: np.asarray(ak.to_numpy(ak.flatten(ak.fill_none(v, np.nan), axis=None)), dtype=float)
    out.update({"census_" + k: fl(v[ok]) for k, v in
                {"qual": qual, "ana": ana, "dy": trk["dy"], "dx": trk["dx"], "dz": trk["dz"],
                 "t0": trk["t"]}.items()})
    out["ntracks"] = int(ak.sum(ak.num(st)))
    tr = {k: v[ana] for k, v in trk.items()}
    out.update({k: fl(v) for k, v in angles(tr["dx"], tr["dy"], tr["dz"]).items()})
    ntr = ak.num(tr["x"])
    evt = np.repeat(np.arange(len(a)), ak.to_numpy(ntr))
    # ---- coincidences (slides 5-7)
    ex = a["crvcoincs.sectorType"] == EX_TYPE
    cc = {c: a[f"crvcoincs.pos.fCoordinates.f{c.upper()}"][ex] for c in "xyz"}
    cc["t"] = a["crvcoincs.time"][ex]
    pairs = ak.cartesian({"t": ak.zip(tr), "c": ak.zip(cc)}, axis=1, nested=True)
    xt, zt, tt = line_at(pairs.t, pairs.c.y)
    dt = pairs.c.t - tt
    out["res_dt"] = fl(dt)
    intime = np.abs(dt - off["dt_coinc"]) < IN_TIME_NS
    out["res_dx"], out["res_dz"] = fl((pairs.c.x - xt)[intime]), fl((pairs.c.z - zt)[intime])
    near = intime & (np.abs(pairs.c.z - zt) < COINC_DZ_MM)
    out["coinc_hit"] = fl(ak.any(near, axis=2))
    ymid = float(np.mean(g["y"]))
    xm, zm, _ = line_at(tr, ymid)
    xm, zm = fl(xm), fl(zm)
    out["edge"] = np.minimum.reduce([EX_HALF_LENGTH - np.abs(xm), zm - g["zbox"][0], g["zbox"][1] - zm])
    # ---- pulses, per layer (slides 11-18, 24)
    pex = a["crvpulses.sectorId"] == EX_SECTOR
    pb = a["crvpulses.barId"][pex]
    lookup = lambda arr: ak.unflatten(arr[ak.to_numpy(ak.flatten(pb))], ak.num(pb))
    pul = {"t": a["crvpulses.time"][pex], "L": lookup(g["layer_of_bar"]), "k": lookup(g["counter_of_bar"])}
    out.update(layer_hits(g, tr, ntr, pul, off["dt_pulse"], fl, ""))
    # ---- MC truth: the muon crossing the EX inner face, every event
    if kind == "mc":
        b = t.arrays(TRUTH, library="ak")
        mu = np.abs(b["crvcoincsmcplane_ex.pdgId"]) == 13
        first = lambda v: ak.to_numpy(ak.fill_none(ak.firsts(v[mu]), np.nan))
        cols = {"x": "pos.fCoordinates.fX", "y": "pos.fCoordinates.fY", "z": "pos.fCoordinates.fZ",
                "t": "time", "dx": "dir.fCoordinates.fX", "dy": "dir.fCoordinates.fY", "dz": "dir.fCoordinates.fZ"}
        tru = {k: first(b[f"crvcoincsmcplane_ex.{v}"]) for k, v in cols.items()}
        has = np.isfinite(tru["x"])
        # multiple scattering: analysis track minus truth at the EX inner face, vs muon energy
        T = {k: fl(v) for k, v in tr.items()}
        yt = tru["y"][evt]
        s_ = (yt - T["y"]) / T["dy"]
        out["ms_dz"] = T["z"] + s_ * T["dz"] - tru["z"][evt]
        out["ms_ke"] = first(b["crvcoincsmcplane_ex.kineticEnergy"])[evt]
        ntu = ak.Array(np.ones(int(has.sum()), dtype=np.int64))
        tru = {k: ak.unflatten(v[has], ntu) for k, v in tru.items()}
        pul_t = {k: v[has] for k, v in pul.items()}
        out["truth_thz"] = fl(angles(tru["dx"], tru["dy"], tru["dz"])["thz"])
        out.update(layer_hits(g, tru, ntu, pul_t, off["dt_truth"], fl, "truth_"))
    return out


def layer_hits(g, tr, ntr, pul, dt0, fl, pre):
    """Crossing of each EX layer, its bin, and whether an in-time pulse is there."""
    res = {}
    n = int(ak.sum(ntr))
    gidx = ak.unflatten(np.arange(n), ak.to_numpy(ntr))
    k_all = np.full((n, 4), -1)
    for L, yl in enumerate(g["y"]):
        x, z, _ = line_at(tr, yl)
        x, z = fl(x), fl(z)
        k, typ, idx, d, j = classify(g, L, z)
        typ = np.where(np.abs(x) < EX_HALF_LENGTH - X_MARGIN, typ, -1)
        k_all[:, L] = np.where(typ >= 0, k, -1)
        res.update({f"{pre}z{L}": z, f"{pre}typ{L}": typ, f"{pre}idx{L}": idx, f"{pre}d{L}": d, f"{pre}gap{L}": j})
    # pulses: track time at the pulse's layer, and its counter
    pr = ak.cartesian({"t": ak.zip({**tr, "g": gidx}), "p": ak.zip(pul)}, axis=1)
    yl = np.asarray(g["y"])[ak.to_numpy(ak.flatten(pr.p.L))] if len(ak.flatten(pr.p.L)) else np.array([])
    tp = ak.to_numpy(ak.flatten(pr.p.t))
    gg = ak.to_numpy(ak.flatten(pr.t.g))
    Ls = ak.to_numpy(ak.flatten(pr.p.L))
    ks = ak.to_numpy(ak.flatten(pr.p.k))
    cols = {c: ak.to_numpy(ak.flatten(pr.t[c])) for c in ("y", "dy", "t")}
    ttrk = cols["t"] + (yl - cols["y"]) / cols["dy"] / C_MM_PER_NS
    dtp = tp - ttrk
    res[f"{pre}pulse_dt"] = dtp[np.abs(dtp) < 2000][:200000]
    intime = (dtp - dt0 > PULSE_WINDOW[0]) & (dtp - dt0 < PULSE_WINDOW[1]) & (k_all[gg, Ls] >= 0)
    for name, nk in (("hit", NARROW), ("hitw", WIDE)):
        hit = np.zeros((n, 4), bool)
        good = intime & (np.abs(ks - k_all[gg, Ls]) <= nk)
        hit[gg[good], Ls[good]] = True
        for L in range(4):
            res[f"{pre}{name}{L}"] = hit[:, L]
    return res


def run(files, kind, g, off, nproc):
    with Pool(nproc) as pool:
        parts = pool.map(process, [(f, kind, g, off) for f in files])
    out = {}
    for k in parts[0]:
        v = [p[k] for p in parts]
        out[k] = sum(v) if np.isscalar(v[0]) else np.concatenate(v)
    return out


# --------------------------------------------------------------------------- numbers
def peak(v, lo, hi, width):
    h, e = np.histogram(v, bins=int((hi - lo) / width), range=(lo, hi))
    j = h.argmax()
    mode = 0.5 * (e[j] + e[j + 1])
    return float(np.median(v[np.abs(v - mode) < 2 * width]))


def core(v, centre, half):
    """Peak and core sigma: iterate mean/rms within +-1.5 sigma, starting from +-half."""
    mu, sig = centre, half / 1.5
    for _ in range(20):
        w = v[np.abs(v - mu) < 1.5 * sig]
        mu, sig = float(np.mean(w)), float(np.std(w)) / 0.7737   # rms of a gaussian cut at 1.5 sigma
    return mu, sig


def eff(k, n):
    n = np.asarray(n, float)
    e = np.divide(k, n, out=np.full_like(n, np.nan), where=n > 0)
    return e, np.sqrt(np.divide(e * (1 - e), n, out=np.full_like(n, np.nan), where=n > 0))


def stack(s, pre):
    """Concatenate the four layers: crossing z, bin type and index, distance to the nearest
    module gap and which gap, narrow and wide hit, layer, theta_z."""
    ths = s[pre + "thz"] if pre else s["thz"]
    cols = {c: np.concatenate([s[f"{pre}{c}{L}"] for L in range(4)]) for c in ("z", "typ", "idx", "d", "hit", "hitw", "gap")}
    cols["L"] = np.repeat(np.arange(4), len(ths))
    cols["thz"] = np.tile(ths, 4)
    m = cols["typ"] >= 0
    return {k: v[m] for k, v in cols.items()}


def gap_width(d, hit, wts=None):
    """w = D / eps0 (see the module docstring), with its statistical error."""
    wts = np.ones_like(d) if wts is None else wts
    pl = (np.abs(d) > PLATEAU[0]) & (np.abs(d) < PLATEAU[1])
    inn = np.abs(d) < GAP_WINDOW
    n0, k0 = wts[pl].sum(), (wts * hit)[pl].sum()
    if n0 < 10 or inn.sum() < 10:
        return np.nan, np.nan, int(inn.sum())
    e0 = k0 / n0
    rho = wts[pl].sum() / (2 * (PLATEAU[1] - PLATEAU[0]))      # denominators per mm
    nin, kin = wts[inn].sum(), (wts * hit)[inn].sum()
    D = (e0 * nin - kin) / rho
    w = D / e0
    # binomial errors on the inner pass count and on eps0
    ein = kin / nin
    var_k = nin * ein * (1 - ein)
    var_e0 = e0 * (1 - e0) / n0
    dw = np.sqrt(var_k / (rho * e0) ** 2 + (kin / (rho * e0 ** 2)) ** 2 * var_e0)
    return float(w), float(dw), int(inn.sum())


def width_model(W, ell, tan):
    """Effective width for one track: a layer gives no pulse if the track's path through
    either counter is shorter than ell. Across a gap of width W in a layer of thickness t
    that is a band of width W - t tan(theta) + 2 ell sin(theta), or nothing."""
    return np.maximum(0.0, W - LAYER_THICKNESS * tan + 2 * ell * tan / np.sqrt(1 + tan ** 2))


def fit_width(w, dw, tan_by_bin, ell=None):
    """Fit W (and ell, if not given) to w(theta), the model averaged over each bin's tracks."""
    ok = np.isfinite(w) & np.isfinite(dw) & (dw > 0)

    def chi2(W, L):
        model = np.array([np.mean(width_model(W, L, tb)) if len(tb) else np.nan for tb in tan_by_bin])
        return float(np.sum(((w - model)[ok] / dw[ok]) ** 2))
    if ell is None:
        Ws, Ls = np.linspace(0, 10, 201), np.linspace(0, 10, 201)
        c = np.array([[chi2(W, L) for L in Ls] for W in Ws])
        i, j = np.unravel_index(c.argmin(), c.shape)
        inside = c < c[i, j] + 1
        dW = 0.5 * (Ws[inside.any(1)].max() - Ws[inside.any(1)].min())
        dL = 0.5 * (Ls[inside.any(0)].max() - Ls[inside.any(0)].min())
        return {"W": float(Ws[i]), "dW": float(dW), "ell": float(Ls[j]), "dell": float(dL),
                "chi2": float(c[i, j]), "nbins": int(ok.sum())}
    Ws = np.linspace(0, 15, 3001)
    c = np.array([chi2(W, ell) for W in Ws])
    i = c.argmin()
    inside = Ws[c < c[i] + 1]
    return {"W": float(Ws[i]), "dW": float(0.5 * (inside.max() - inside.min())), "ell": ell, "dell": None,
            "chi2": float(c[i]), "nbins": int(ok.sum())}


# ------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", required=True, help="glob of data ntuples")
    ap.add_argument("--mc", required=True, help="glob of MC ntuples")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("-j", "--nproc", type=int, default=8)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    files = {"data": sorted(glob.glob(a.data)), "mc": sorted(glob.glob(a.mc))}
    for k, v in files.items():
        if not v:
            raise SystemExit(f"no {k} files match")
    g = ex_geometry(files["mc"])
    gd = ex_geometry(files["data"])
    for L in range(4):
        if not np.allclose(g["layers"][L]["zc"], gd["layers"][L]["zc"]):
            raise SystemExit(f"data and MC EX geometry differ in layer {L}")
    S = {"files": {k: len(v) for k, v in files.items()},
         "geometry": {"layer_y": g["y"].tolist(), "module_gap_centres": [l["gap_centre"].tolist() for l in g["layers"]],
                      "zbox": list(g["zbox"])}}

    # pass 1: offsets from the first files
    off = {}
    for kind in ("data", "mc"):
        o0 = {"zshift": 0.0, "dt_coinc": 0.0, "dt_pulse": 0.0, "dt_truth": 0.0}
        p = run(files[kind][:NOFFSET_FILES], kind, g, o0, a.nproc)
        dtc = peak(p["res_dt"], -2000, 2000, 5)
        o0["dt_coinc"] = dtc
        p = run(files[kind][:NOFFSET_FILES], kind, g, o0, a.nproc)
        zs = peak(p["res_dz"], -3000, 3000, 10) if kind == "data" else 0.0   # coinc z - track z
        o = {"zshift": zs, "dt_coinc": dtc, "dt_pulse": peak(p["pulse_dt"], -2000, 2000, 2),
             "dt_truth": peak(p["truth_pulse_dt"], -2000, 2000, 2) if kind == "mc" else 0.0}
        off[kind] = o
    S["offsets"] = off

    # pass 2: everything, with the offsets applied
    R = {k: run(files[k], k, g, off[k], a.nproc) for k in ("data", "mc")}
    R0 = {k: run(files[k][:NOFFSET_FILES], k, g, dict(off[k], zshift=0.0), a.nproc) for k in ("data", "mc")}
    S["livetime_s"] = livetime(files)
    S["census"] = census(R, S)
    S["residuals"] = residuals(R0, R, off, a.out)
    S["tracker_shift_mm"] = off["data"]["zshift"]
    plot_census(R, S, a.out)
    S["coinc_eff"] = coinc_edge(R, a.out)
    S["by_type"], S["thz"] = by_type(R, g, a.out)
    S["gap_width"] = widths(R, a.out)
    plot_vs_z(R, g, a.out)
    plot_folded(R, a.out)
    S["narrow_vs_wide"] = narrow_vs_wide(R, a.out)
    S["scattering"] = scattering(R, a.out)
    with open(os.path.join(a.out, "summary.json"), "w") as f:
        json.dump(S, f, indent=1, default=float)
    print(json.dumps({k: S[k] for k in ("offsets", "census", "coinc_eff", "gap_width")}, indent=1, default=float))


def livetime(files):
    """Data: events processed x the 100 us event window (events from each job's reco.log).
    MC: the summed cosmic live time of every subrun."""
    nev = 0
    for f in files["data"]:
        log = os.path.join(os.path.dirname(f), "reco.log")
        m = re.search(r"Events total = (\d+)", open(log).read())
        if not m:
            raise RuntimeError(f"no event count in {log}")
        nev += int(m.group(1))
    lt_mc, nev_mc = 0.0, 0
    for f in files["mc"]:
        F = uproot.open(f)
        lt_mc += float(F["EventNtuple/subrunNtuple"]["cosmicLivetime"].array(library="np").sum())
        nev_mc += F[TREE].num_entries
    return {"data": nev * EVENT_WINDOW_NS * 1e-9, "data_events": nev, "mc": lt_mc, "mc_events": nev_mc}


def census(R, S):
    out = {}
    for k in ("data", "mc"):
        r, lt = R[k], S["livetime_s"][k]
        n = len(r["census_qual"])
        out[k] = {"events": S["livetime_s"][k + "_events"], "livetime_s": lt, "tracks": n,
                  "tracks_per_s": n / lt, "quality_frac": float(np.mean(r["census_qual"])),
                  "nplanes_frac_of_quality": float(np.sum(r["census_ana"]) / max(1, np.sum(r["census_qual"]))),
                  "analysis_tracks": int(np.sum(r["census_ana"])), "analysis_tracks_per_s": float(np.sum(r["census_ana"]) / lt)}
    return out


def plot_census(R, S, out):
    """Slides 4 and 22: analysis-track directions and t0, area-normalised and per live time."""
    fig, axs = plt.subplots(2, 4, figsize=(16, 7))
    spec = [("zenith", (0, 90), "zenith angle from -y [deg]"), ("thz", (0, 90), r"$\theta_z$ [deg]"),
            ("phi", (-180, 180), r"azimuth $\phi$ [deg]"), ("t0", (0, 1.0e5), "track t at the tracker [ns]")]
    for k in ("data", "mc"):
        r = R[k]
        a = r["census_ana"] > 0
        ang = angles(r["census_dx"][a], r["census_dy"][a], r["census_dz"][a])
        ang["t0"] = r["census_t0"][a]
        lab = "run 124155" if k == "data" else "MC MDC2025au"
        lt = S["livetime_s"][k]
        for i, (v, rng, xl) in enumerate(spec):
            axs[0, i].hist(ang[v], bins=60, range=rng, density=True, histtype="step", color=COLORS[k], label=lab)
            axs[1, i].hist(ang[v], bins=60, range=rng, weights=np.full(len(ang[v]), 1 / lt), histtype="step", color=COLORS[k], label=lab)
            axs[1, i].set_xlabel(xl)
    axs[0, 0].legend(fontsize=8)
    axs[0, 0].set_ylabel("tracks (area normalised)")
    axs[1, 0].set_ylabel("tracks / s / bin")
    c = S["census"]
    fig.suptitle(f"analysis tracks: data {c['data']['analysis_tracks_per_s']:.1f}/s, MC {c['mc']['analysis_tracks_per_s']:.1f}/s")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "slide04_census.png"), dpi=110)
    plt.close(fig)


def residuals(R0, R, off, out):
    """Slide 5: EX coincidence minus extrapolated track; data before the z shift."""
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.6))
    res = {}
    for k in ("data", "mc"):
        r0 = R0[k]
        lab = "run 124155" if k == "data" else "MC MDC2025au"
        for i, (v, rng, b, xl, half) in enumerate([
                ("res_dt", (-1000, 1000), 400, r"$\Delta t$ = t(coinc) - t(track) [ns]", 30),
                ("res_dx", (-6000, 6000), 120, r"$\Delta x$ along the bars [mm]", 1500),
                ("res_dz", (-2000, 2000), 200, r"$\Delta z$ stacking [mm]", 150)]):
            x = r0[v]
            c = off[k]["dt_coinc"] if v == "res_dt" else peak(x, *rng, (rng[1] - rng[0]) / b)
            mu, sig = core(x, c, half)
            res.setdefault(k, {})[v] = {"peak": mu, "core_sigma": sig, "n": int(len(x))}
            ax[i].hist(x, bins=b, range=rng, density=True, histtype="step", color=COLORS[k],
                       label=f"{lab}: {mu:+.0f} $\\pm$ {sig:.0f}")
            ax[i].set_xlabel(xl)
            ax[i].legend(fontsize=8)
    # after the shift (data, all files): should sit at zero
    mu, sig = core(R["data"]["res_dz"], 0.0, 150)
    res["data"]["res_dz_after_shift"] = {"peak": mu, "core_sigma": sig, "n": int(len(R["data"]["res_dz"]))}
    fig.suptitle(f"EX coincidence minus extrapolated track (first {NOFFSET_FILES} files of each sample, data tracks not shifted)")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "slide05_residuals.png"), dpi=110)
    plt.close(fig)
    return res


def coinc_edge(R, out):
    bins = np.array([0, 25, 50, 75, 100, 150, 200, 250, 300, 400, 500, 700, 1000, 1300, 1700])
    fig, ax = plt.subplots(figsize=(7, 4))
    res = {}
    for k in ("data", "mc"):
        r = R[k]
        e = r["edge"]
        n, _ = np.histogram(e, bins)
        h, _ = np.histogram(e[r["coinc_hit"] > 0], bins)
        y, dy = eff(h, n)
        xc = 0.5 * (bins[1:] + bins[:-1])
        ax.errorbar(xc, y, dy, fmt="o", ms=4, color=COLORS[k], label="run 124155" if k == "data" else "MC MDC2025au")
        far = e >= 800
        ep, dep = eff(np.sum(r["coinc_hit"][far] > 0), np.sum(far))
        res[k] = {"plateau_eff": float(ep), "plateau_err": float(dep), "n_plateau": int(far.sum()),
                  "bins": bins.tolist(), "eff": y.tolist()}
    ax.set_xlabel("distance of the track to the nearest EX edge, at mid-sector [mm]")
    ax.set_ylabel("EX coincidence efficiency")
    ax.set_ylim(0, 1.05)
    ax.legend()
    ax.set_title(f"plateau (d >= 800 mm): MC {100 * res['mc']['plateau_eff']:.1f}%, data {100 * res['data']['plateau_eff']:.1f}%")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "slide07_coinc_eff_edge.png"), dpi=110)
    plt.close(fig)
    return res


def thz_weights(R):
    """Per-crossing weights that give MC reco the data theta_z distribution (flat region)."""
    bins = np.arange(0, 72, 2.0)
    sd, sm = stack(R["data"], ""), stack(R["mc"], "")
    fd = (sd["z"] > FLAT_Z[0]) & (sd["z"] < FLAT_Z[1])
    fm = (sm["z"] > FLAT_Z[0]) & (sm["z"] < FLAT_Z[1])
    hd, _ = np.histogram(sd["thz"][fd], bins, density=True)
    hm, _ = np.histogram(sm["thz"][fm], bins, density=True)
    w = np.divide(hd, hm, out=np.zeros_like(hd), where=hm > 0)
    return sm, w[np.clip(np.digitize(sm["thz"], bins) - 1, 0, len(w) - 1)], bins, hd, hm


def by_type(R, g, out):
    """Slides 14, 16, 17: efficiency per bin type and layer in the flat region."""
    sm, wm, tb, hd, hm = thz_weights(R)
    samples = {"truth": stack(R["mc"], "truth_"), "mc": sm, "data": stack(R["data"], ""), "mcw": sm}
    res = {}
    fig, ax = plt.subplots(1, 4, figsize=(17, 3.8), sharey=True)
    for i, (k, s) in enumerate(samples.items()):
        w = wm if k == "mcw" else np.ones(len(s["z"]))
        flat_ = (s["z"] > FLAT_Z[0]) & (s["z"] < FLAT_Z[1])
        res[k] = {}
        for t, name in enumerate(TYPES):
            ys, es = [], []
            for L in range(4):
                m = flat_ & (s["typ"] == t) & (s["L"] == L)
                nn, kk = w[m].sum(), (w * s["hit"])[m].sum()
                neff = nn ** 2 / max((w[m] ** 2).sum(), 1e-9)       # effective entries for the error
                e = kk / nn if nn > 0 else np.nan
                ys.append(e)
                es.append(np.sqrt(e * (1 - e) / neff) if nn > 0 else np.nan)
            res[k][name] = {"eff": ys, "err": es}
            ax[i].errorbar(np.arange(4) + 0.05 * (t - 1), 100 * np.array(ys), 100 * np.array(es), fmt="o-", ms=4, label=name)
        res[k]["counter_minus_module_gap_pp"] = [100 * (res[k]["counter"]["eff"][L] - res[k]["module gap"]["eff"][L]) for L in range(4)]
        ax[i].set_title({"truth": "MC truth", "mc": "MC reco", "data": "data (run 124155)", "mcw": r"MC reco, data $\theta_z$"}[k])
        ax[i].set_xticks(range(4), [f"L{L}" for L in range(4)])
        ax[i].grid(alpha=0.3)
    ax[0].set_ylabel("pulse efficiency [%]")
    ax[0].legend(fontsize=8)
    fig.suptitle(f"folded over {FLAT_Z[0]:.0f} < z < {FLAT_Z[1]:.0f} mm")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "slide14_17_eff_by_type.png"), dpi=110)
    plt.close(fig)
    # slide 16: the denominators' theta_z
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.stairs(hd * 2, tb, color=COLORS["data"], label="data")
    ax.stairs(hm * 2, tb, color=COLORS["mc"], label="MC")
    crit = np.degrees(np.arctan(NOMINAL_GAP / LAYER_THICKNESS))
    ax.axvline(crit, ls="--", color="k", lw=1)
    ax.set_xlabel(r"$\theta_z$ = atan(|p$_z$|/|p$_y$|) [deg]")
    ax.set_ylabel("fraction / 2 deg")
    ax.set_title(f"layer crossings, {FLAT_Z[0]:.0f} < z < {FLAT_Z[1]:.0f} mm; dashed: full-miss angle {crit:.1f} deg")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out, "slide16_thetaz.png"), dpi=110)
    plt.close(fig)
    # slide 24: efficiency vs theta_z, counters and module gaps
    fig, ax = plt.subplots(1, 2, figsize=(12, 4), sharey=True)
    thz_res = {}
    for k in ("mc", "data"):
        s = samples[k]
        flat_ = (s["z"] > FLAT_Z[0]) & (s["z"] < FLAT_Z[1])
        for i, t in enumerate((0, 2)):
            m = flat_ & (s["typ"] == t)
            b = np.arange(0, 52, 2.0)
            n, _ = np.histogram(s["thz"][m], b)
            h, _ = np.histogram(s["thz"][m & (s["hit"] > 0)], b)
            y, dy = eff(h, n)
            ok = n >= 20
            ax[i].errorbar((0.5 * (b[1:] + b[:-1]))[ok], 100 * y[ok], 100 * dy[ok], fmt="o", ms=4, color=COLORS[k],
                           label="data" if k == "data" else "MC")
            ax[i].set_title(TYPES[t])
            thz_res.setdefault(k, {})[TYPES[t]] = {"bins": b.tolist(), "eff": y.tolist(), "n": n.tolist()}
    for i in range(2):
        ax[i].axvline(crit, ls="--", color="k", lw=1)
        ax[i].set_xlabel(r"$\theta_z$ [deg]")
        ax[i].grid(alpha=0.3)
    ax[0].set_ylabel("pulse efficiency [%]")
    ax[0].legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out, "slide24_eff_vs_thetaz.png"), dpi=110)
    plt.close(fig)
    return res, thz_res


def widths(R, out):
    """Slide 18: effective gap width per theta_z bin, and the fitted gap width."""
    samples = {"truth": (stack(R["mc"], "truth_"), 2.0), "mc": (stack(R["mc"], ""), 2.0), "data": (stack(R["data"], ""), 4.0)}
    res = {}
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for k, (s, step) in samples.items():
        b = np.arange(0, 30 + step, step)
        near = np.abs(s["d"]) < PLATEAU[1]
        W, dW, tans, N = [], [], [], 0
        for lo, hi in zip(b[:-1], b[1:]):
            m = near & (s["thz"] >= lo) & (s["thz"] < hi)
            w, dw, n = gap_width(s["d"][m], s["hitw"][m].astype(float))
            W.append(w)
            dW.append(dw)
            tans.append(np.tan(np.radians(s["thz"][m & (np.abs(s["d"]) < GAP_WINDOW)])))
            N += n
        W, dW = np.array(W), np.array(dW)
        fit = fit_width(W, dW, tans, None if k == "truth" else res["truth"]["fit"]["ell"])
        res[k] = {"bins": b.tolist(), "w": W.tolist(), "dw": dW.tolist(), "fit": fit, "n_inner": N}
        # all theta_z together, and each gap on its own
        res[k]["all_thz"] = gap_width(s["d"][near], s["hitw"][near].astype(float))[:2]
        res[k]["per_gap"] = [gap_width(s["d"][near & (s["gap"] == j)], s["hitw"][near & (s["gap"] == j)].astype(float))[:2]
                             for j in range(3)]
        xc = 0.5 * (b[1:] + b[:-1])
        lab = {"data": "data, track", "mc": "MC, track", "truth": "MC, truth"}[k]
        ax.errorbar(xc, W, dW, fmt="o", ms=4, color=COLORS[k], label=f"{lab}: W = {fit['W']:.2f} $\\pm$ {fit['dW']:.2f} mm")
        th = np.linspace(0, 30, 200)
        ax.plot(th, width_model(fit["W"], fit["ell"], np.tan(np.radians(th))), color=COLORS[k], lw=1, alpha=0.7)
    ell = res["truth"]["fit"]["ell"]
    ax.plot([], [], "k-", lw=1, label=f"model, minimum path {ell:.1f} mm (fitted to truth)")
    ax.axhline(0, color="grey", lw=0.5)
    ax.set_xlabel(r"$\theta_z$ [deg]")
    ax.set_ylabel(r"$\hat w$ = D / $\varepsilon_0$ [mm]")
    ax.legend(fontsize=8)
    ax.set_title("effective module-gap width, three gaps and four layers combined")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "slide18_gap_width.png"), dpi=110)
    plt.close(fig)
    fd, fm = res["data"]["fit"], res["mc"]["fit"]
    if fm["W"] > 0:
        r = fd["W"] / fm["W"]
        dr = np.hypot(fd["dW"] / fm["W"], r * fm["dW"] / fm["W"])
    else:
        r = dr = float("nan")
    res["ratio_data_mc"] = {"R": float(r), "dR": float(dr)}
    template_fits(R, samples, res, out)
    return res


def template_fits(R, samples, res, out):
    """Cross-check of the integral: fit the folded wide-hit profile (|d| < 400 mm) with a
    quadratic background minus a dip shaped like the MC smearing (track minus truth at
    EX), so slow trends in efficiency cannot pose as a gap. Also the width each track
    sample should show if its gaps were those of MC truth: the truth model averaged over
    the sample's own theta_z."""
    b = np.arange(-400, 401, 10.0)
    xc = 0.5 * (b[1:] + b[:-1])
    dz = R["mc"]["ms_dz"]
    dz = dz[np.isfinite(dz)]
    K, _ = np.histogram(np.concatenate([dz, -dz]), b, density=True)
    tf = res["truth"]["fit"]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    for i, k in enumerate(("mc", "data")):
        s = samples[k][0]
        n, _ = np.histogram(s["d"], b)
        h, _ = np.histogram(s["d"][s["hitw"] > 0], b)
        e = h / n
        se = np.sqrt(e * (1 - e) / n)
        X = np.c_[np.ones_like(xc), xc, xc ** 2, -K]
        Wt = 1 / se ** 2
        cov = np.linalg.inv(X.T @ (X * Wt[:, None]))
        p = cov @ (X.T @ (Wt * e))
        chi2 = float(np.sum(((X @ p - e) / se) ** 2))
        inner = np.abs(s["d"]) < GAP_WINDOW
        expect = float(np.mean(width_model(tf["W"], tf["ell"], np.tan(np.radians(s["thz"][inner])))))
        res[k]["template"] = {"w": float(p[3] / p[0]), "dw": float(np.sqrt(cov[3, 3]) / p[0]), "chi2": chi2,
                              "ndf": len(e) - 4, "expected_from_truth_model": expect}
        ax[i].errorbar(xc, 100 * e, 100 * se, fmt="o", ms=2, color=COLORS[k], label="wide-hit efficiency")
        ax[i].plot(xc, 100 * (X[:, :3] @ p[:3]), "k:", lw=1, label="background")
        ax[i].plot(xc, 100 * (X @ p), "k-", lw=1, label=f"fit: w = {p[3] / p[0]:.2f} $\\pm$ {np.sqrt(cov[3, 3]) / p[0]:.2f} mm")
        ax[i].plot(xc, 100 * (X[:, :3] @ p[:3] - p[0] * expect * K), color="darkorange", lw=1,
                   label=f"expected from MC truth: w = {expect:.2f} mm")
        ax[i].set_title("MC, track" if k == "mc" else "data, track")
        ax[i].set_xlabel("signed distance to the nearest module gap [mm]")
        ax[i].legend(fontsize=7)
        ax[i].grid(alpha=0.3)
    ax[0].set_ylabel("pulse efficiency, +-6 counters [%]")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "extra_template_fit.png"), dpi=110)
    plt.close(fig)


def plot_vs_z(R, g, out):
    """Slides 11, 13, 15: efficiency vs z per layer: truth (fine), and per counter/gap bin."""
    fig, ax = plt.subplots(4, 3, figsize=(17, 9), sharex=True, sharey=True)
    for c, (k, pre) in enumerate([("truth", "truth_"), ("mc", ""), ("data", "")]):
        s = stack(R["mc" if k != "data" else "data"], pre)
        for L in range(4):
            a = ax[L, c]
            lay = g["layers"][L]
            for t, col in enumerate(("navy", "crimson", "darkorange")):
                m = (s["L"] == L) & (s["typ"] == t)
                nb = {0: 32, 1: 63, 2: 63}[t]
                n = np.bincount(s["idx"][m], minlength=nb)
                h = np.bincount(s["idx"][m], weights=s["hit"][m].astype(float), minlength=nb)
                y, dy = eff(h, n)
                pos = 0.5 * (lay["zc"][0::2] + lay["zc"][1::2]) if t == 0 else 0.5 * (lay["hi"][:-1] + lay["lo"][1:])
                ok = n >= 20
                a.errorbar(pos[ok], 100 * y[ok], 100 * dy[ok], fmt="o", ms=3, color=col,
                           label=TYPES[t] if (L == 0 and c == 0) else None)
            for gc in lay["gap_centre"]:
                a.axvline(gc, color="darkorange", alpha=0.3, lw=3)
            a.set_ylim(30, 105)
            a.grid(alpha=0.3)
            if c == 0:
                a.set_ylabel(f"L{L} eff [%]")
        ax[0, c].set_title({"truth": "MC truth", "mc": "MC reco", "data": "data (run 124155), tracker shifted"}[k])
        ax[3, c].set_xlabel("stacking coordinate z [mm]")
    ax[0, 0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "slide13_15_eff_vs_z.png"), dpi=110)
    plt.close(fig)


def plot_folded(R, out):
    """Slide 12: efficiency vs signed distance to the nearest module gap. Truth in 1 mm
    bins; the track extrapolation smears the dip over about +-50 mm, so 10 mm bins."""
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.8))
    for i, (k, pre, b) in enumerate([("truth", "truth_", np.arange(-60, 61, 1.0)),
                                     ("mc", "", np.arange(-400, 401, 10.0)),
                                     ("data", "", np.arange(-400, 401, 10.0))]):
        s = stack(R["mc" if k != "data" else "data"], pre)
        n, _ = np.histogram(s["d"], b)
        h, _ = np.histogram(s["d"][s["hitw"] > 0], b)
        y, dy = eff(h, n)
        xc = 0.5 * (b[1:] + b[:-1])
        ax[i].errorbar(xc, 100 * y, 100 * dy, fmt="o", ms=2, color=COLORS[k])
        ax[i].axvspan(-NOMINAL_GAP / 2, NOMINAL_GAP / 2, color="gold", alpha=0.4)
        if k != "truth":
            for e in (-GAP_WINDOW, GAP_WINDOW):
                ax[i].axvline(e, color="grey", ls=":", lw=1)
        ax[i].set_title({"truth": "MC truth (1 mm bins)", "mc": "MC, track extrapolation (10 mm bins)",
                         "data": "data, track extrapolation (10 mm bins)"}[k])
        ax[i].set_xlabel("signed distance to the nearest module gap [mm]")
        ax[i].grid(alpha=0.3)
    ax[0].set_ylabel("pulse efficiency [%]")
    fig.suptitle("folded over all layers and module gaps; dotted: the window the gap width integrates")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "slide12_folded_gap.png"), dpi=110)
    plt.close(fig)


def narrow_vs_wide(R, out):
    """Is a drop in efficiency along z the CRV, or the pointing? The wide hit tolerates
    +-6 counters of extrapolation error; the narrow hit only +-1."""
    b = np.arange(-2800, 801, 200.0)
    fig, ax = plt.subplots(figsize=(8, 4))
    res = {}
    for k in ("mc", "data"):
        s = stack(R[k], "")
        n, _ = np.histogram(s["z"], b)
        res[k] = {"bins": b.tolist()}
        for hk, ls in (("hit", "-"), ("hitw", "--")):
            h, _ = np.histogram(s["z"][s[hk] > 0], b)
            y, dy = eff(h, n)
            res[k][hk] = y.tolist()
            ax.errorbar(0.5 * (b[1:] + b[:-1]), 100 * y, 100 * dy, fmt="o" + ls, ms=3, color=COLORS[k],
                        label=f"{'data' if k == 'data' else 'MC'}, {'+-1 counter' if hk == 'hit' else '+-6 counters'}")
    ax.set_xlabel("stacking coordinate z [mm] (data tracks shifted)")
    ax.set_ylabel("pulse efficiency, all layers [%]")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(out, "extra_narrow_vs_wide_z.png"), dpi=110)
    plt.close(fig)
    return res


def scattering(R, out):
    """MC: extrapolated track minus true crossing at the EX inner face, by muon energy."""
    dz, ke = R["mc"]["ms_dz"], R["mc"]["ms_ke"]
    ok = np.isfinite(dz) & np.isfinite(ke)
    edges = [0, 1000, 2000, 4000, 8000, 20000, 1e9]
    res = []
    fig, ax = plt.subplots(figsize=(7, 4))
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = ok & (ke >= lo) & (ke < hi)
        if m.sum() < 50:
            continue
        q = np.percentile(dz[m], [16, 84])
        res.append({"ke_lo": lo, "ke_hi": hi, "n": int(m.sum()), "half_width_68": float(0.5 * (q[1] - q[0]))})
        ax.hist(np.clip(dz[m], -400, 400), bins=80, range=(-400, 400), density=True, histtype="step",
                label=f"{lo / 1000:g}-{hi / 1000:g} GeV: {res[-1]['half_width_68']:.0f} mm")
    q = np.percentile(dz[ok], [16, 84])
    ax.set_xlabel("extrapolated track z minus true z at the EX inner face [mm]")
    ax.set_ylabel("tracks (area normalised)")
    ax.set_title(f"MC, analysis tracks: 68% half-width {0.5 * (q[1] - q[0]):.0f} mm")
    ax.legend(fontsize=8, title="muon kinetic energy at EX")
    fig.tight_layout()
    fig.savefig(os.path.join(out, "extra_ms_resolution.png"), dpi=110)
    plt.close(fig)
    return {"half_width_68_all": float(0.5 * (q[1] - q[0])), "by_energy": res}


if __name__ == "__main__":
    main()
