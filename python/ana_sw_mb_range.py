#!/usr/bin/env python3
"""
ana_sw_mb_range.py
==================
Finds the shockwave (SW) and muzzle blast (MB) in each recording by band-pass
filtering, measures the SW->MB delay Dt, and converts Dt to shooter range with
a drag-aware ballistic model. If the folder has a metadata.md (like
sound/Cipatat/metadata.md), the computed range is compared with the reference.

Detection (blind -- does NOT use the reference distance)
--------------------------------------------------------
SW/MB detection uses the pipeline's own classifier (gs_det_events.py; see its
docstring): every event is labelled shockwave or muzzle blast by its frequency
content and time/duration. A SW and an MB are each optional and the MB does not
have to follow a SW. Dt and range are computed only when both are found and the
MB comes later.

Thresholds were tuned on the Cipatat recordings (51.2 kHz, 4 channels).

Range model
-----------
Shooter at distance D, straight bullet path passing the array at miss distance
b (DCPA, taken from metadata.md -- a real detector would estimate it), speed
v(s) = v0 exp(-s / lam). The SW arrival is the minimum over emission points of
t_bullet(s) + |array - p(s)| / c (Fermat), the MB arrives at D / c, and D is
solved from Dt = D/c - t_SW. v0 and lam come from the BALLISTIC table below:
approximate public values chosen independently of the generator's library.

Usage
-----
  python ana_sw_mb_range.py ..\\sound\\Cipatat                  # every *.wav + metadata.md
  python ana_sw_mb_range.py ..\\sound\\Cipatat -c 349           # speed of sound (default 347 m/s)
  python ana_sw_mb_range.py shot.wav -cal 7.62x54R -b 10       # single file, no metadata
  python ana_sw_mb_range.py ..\\sound\\Cipatat -plot            # also save per-file plots
  python ana_sw_mb_range.py ..\\sound\\Cipatat -v               # list every classified event

Output: a table on screen, plus sw_mb_range_results.md / .csv and
sw_mb_range_overview.png in the folder of the first wav (or -o DIR).
"""

import argparse
import csv
import glob
import os
import re

import numpy as np
import soundfile as sf
from scipy.optimize import brentq

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- ballistics
# Detector-side table: v0 = muzzle velocity (m/s), lam = speed decay length (m) for
# v(s) = v0 exp(-s/lam); lam None = constant speed. Approximate public data; the
# cannon entries are rough because the exact round types are not known.
BALLISTIC = {
    "5.56x45":  dict(v0=940,  lam=950,  note="M855/SS109 from M16"),
    "7.62x39":  dict(v0=718,  lam=860,  note="M43 ball"),
    "7.62x54R": dict(v0=830,  lam=1150, note="7N1/LPS"),
    "7.62x51":  dict(v0=850,  lam=1050, note="M80 ball"),
    "8.6x70":   dict(v0=900,  lam=2000, note=".338 Lapua 250 gr"),
    "12.7x99":  dict(v0=890,  lam=2500, note="M33 ball"),
    "14.5x114": dict(v0=1000, lam=3200, note="B-32 (KPV)"),
    "20mm":     dict(v0=1030, lam=2500, note="20 mm HEI/TP, type unknown"),
    "25x137":   dict(v0=1100, lam=3000, note="TP-T/HEI-T, type unknown"),
    "30x173":   dict(v0=1080, lam=3500, note="TP-T"),
    "T40":      dict(v0=1500, lam=6000, note="CTA 40 APFSDS assumed"),
    "105HEP":   dict(v0=730,  lam=4000, note="105 mm HEP (M393-like)"),
    "120BSCC":  dict(v0=1700, lam=3000, note="120 mm sub-calibre training round assumed"),
    "RPG7":     dict(v0=295,  lam=None, note="PG-7V, subsonic: no shockwave"),
}


def ballistic_key(caliber_text):
    """Map a caliber string from metadata.md to a BALLISTIC key."""
    t = caliber_text.lower().replace(" ", "").replace("×", "x").replace(",", ".")
    for pat, key in (("5.56x45", "5.56x45"), ("7.62x39", "7.62x39"), ("7.62x54", "7.62x54R"),
                     ("7.62x51", "7.62x51"), ("8.6x70", "8.6x70"), ("12.7x99", "12.7x99"),
                     ("14.5x114", "14.5x114"), ("25x137", "25x137"), ("30x173", "30x173"),
                     ("t40", "T40"), ("105", "105HEP"), ("120", "120BSCC"), ("pg-7", "RPG7"),
                     ("73mm", "RPG7"), ("20mm", "20mm")):
        if pat in t:
            return key
    return None


def t_bullet(s, v0, lam):
    return s / v0 if lam is None else lam * np.expm1(s / lam) / v0


def v_bullet(s, v0, lam):
    return v0 if lam is None else v0 * np.exp(-s / lam)


def t_shockwave(D, b, v0, lam, c):
    """(first SW arrival time after the shot, local Mach at emission) or (None, None)."""
    L = np.sqrt(max(D * D - b * b, 0.0))
    s = np.linspace(0.0, L + 3 * b, 20001)
    tt = t_bullet(s, v0, lam) + np.sqrt((L - s) ** 2 + b * b) / c
    k = int(np.argmin(tt))
    m = v_bullet(s[k], v0, lam) / c
    return (float(tt[k]), float(m)) if m > 1.0 else (None, None)


def dt_model(D, b, v0, lam, c):
    ts, m = t_shockwave(D, b, v0, lam, c)
    return (None if ts is None else D / c - ts), m


def range_from_dt(dt, b, v0, lam, c, d_max=5000.0):
    """Shooter distance D that reproduces dt (None if no solution)."""
    def f(D):
        v = dt_model(D, b, v0, lam, c)[0]
        return np.nan if v is None else v - dt
    grid = np.linspace(max(1.01 * b, 1.0), d_max, 500)
    vals = np.array([f(D) for D in grid])
    for i in range(len(grid) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] <= 0 < vals[i + 1]:
            return brentq(f, grid[i], grid[i + 1], xtol=0.01)
    return None


def v0_to_match(dt, D, b, lam, c):
    """Muzzle velocity (same lam) that makes the model give exactly D for dt."""
    def f(v0):
        v = dt_model(D, b, v0, lam, c)[0]
        return np.nan if v is None else v - dt
    vs = np.linspace(c * 1.05, 2500, 300)
    vals = np.array([f(v) for v in vs])
    for i in range(len(vs) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and np.sign(vals[i]) != np.sign(vals[i + 1]):
            return brentq(f, vs[i], vs[i + 1], xtol=0.5)
    return None


# --------------------------------------------------------------------------- signal
# The SW/MB classifier lives in the detection pipeline (gs_det_events.py, used by
# gs_det_signal_prepare.py with [detection] method = classify); this script uses the same code.
from gs_det_events import bandpass, envelope, detect_sw_mb   # noqa: E402


def detect(x, fs, verbose=False, mb_in_tail=False):
    ev = detect_sw_mb(x, fs, dict(mb_in_tail=mb_in_tail), verbose=verbose)
    return dict(t_sw=ev["t_sw"], t_mb=ev["t_mb"], band=ev["mb_band"], n_sw=ev["n_sw"],
                sw_snr=ev["sw_features"]["snr_hf"] if ev["sw_features"] else None,
                mb_snr=ev["mb_features"]["snr_lf"] if ev["mb_features"] else None,
                t_settle=ev["t_sw_tail_end"], event_class=ev["event_class"])


# --------------------------------------------------------------------------- metadata
def read_metadata(path):
    """{wav_name: dict(weapon, caliber, D, b)} from the Reference values table of metadata.md."""
    meta = {}
    if not os.path.isfile(path):
        return meta
    for line in open(path, encoding="utf-8"):
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) < 10 or not cells[1].lower().startswith("shot-"):
            continue
        wav = cells[1].split()[0]

        def num(s):
            s = s.replace("−", "-").replace(" ", "")
            return float(s) if re.fullmatch(r"-?\d+(\.\d+)?", s) else None
        meta[wav] = dict(weapon=cells[2], caliber=cells[3], rounds=num(cells[5]),
                         D=num(cells[8]), b=num(cells[9]))
    return meta


# --------------------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description="SW/MB delay and range for each recording.", add_help=False)
    p.add_argument("-h", "-help", "--help", action="help", help="show this help message and exit")
    p.add_argument("inputs", nargs="+", help="wav files and/or folders")
    p.add_argument("-meta", default=None, help="metadata.md (default: metadata.md next to the wavs)")
    p.add_argument("-c", type=float, default=347.0, help="speed of sound m/s (default 347, ~27 C)")
    p.add_argument("-cal", default=None, help="ballistic key for files without metadata, e.g. 7.62x54R")
    p.add_argument("-b", type=float, default=None, help="miss distance (m) for files without metadata")
    p.add_argument("-o", default=None, metavar="DIR", help="folder for the results")
    p.add_argument("-plot", action="store_true", help="also save a plot per file")
    p.add_argument("-v", action="store_true", help="list every classified event and its features")
    p.add_argument("-mb_in_tail", action="store_true",
                   help="if no MB is found after the SW echo tail, accept one inside it (close shots)")
    args = p.parse_args()

    wavs = []
    for item in args.inputs:
        wavs += sorted(glob.glob(os.path.join(item, "*.wav"))) if os.path.isdir(item) else [item]
    wavs = [w for w in wavs if not re.search(r"_(sw|mb)\.wav$", w)]
    if not wavs:
        raise SystemExit("No wav files found.")
    folder = os.path.dirname(os.path.abspath(wavs[0]))
    out_dir = args.o or folder
    os.makedirs(out_dir, exist_ok=True)
    meta = read_metadata(args.meta or os.path.join(folder, "metadata.md"))
    c = args.c

    rows = []
    overview = []
    for w in wavs:
        name = os.path.basename(w)
        x, fs = sf.read(w, dtype="float64", always_2d=True)
        if args.v:
            print(name)
        det = detect(x, fs, verbose=args.v, mb_in_tail=args.mb_in_tail)
        m = meta.get(name, {})
        dup_of = None
        if not m:                                            # e.g. shot-06a -> shot-06
            alt = re.sub(r"^(shot-\d+)[a-z]_", r"\1_", name)
            if alt != name and alt in meta:
                m, dup_of = meta[alt], alt
        key = ballistic_key(m.get("caliber", "")) if m else args.cal
        b = m.get("b") if m else args.b
        D_ref = m.get("D")
        row = dict(file=name, weapon=m.get("weapon", ""), key=key or "", D_ref=D_ref, b=b,
                   t_sw=det["t_sw"], t_mb=det["t_mb"], band=det["band"],
                   sw_snr=det["sw_snr"], mb_snr=det.get("mb_snr"))
        row["found"] = ("SW+MB" if det["t_sw"] is not None and det["t_mb"] is not None else
                        "SW only" if det["t_sw"] is not None else "MB only" if det["t_mb"] is not None else "none")
        row["dt"] = (det["t_mb"] - det["t_sw"]) if row["found"] == "SW+MB" else None
        if row["dt"] is not None and row["dt"] <= 0:
            row["dt"] = None
            row["found"] = "MB before SW"
        bl = BALLISTIC.get(key) if key else None
        row["D_calc"] = row["dt_exp"] = row["err_pct"] = row["v0_fit"] = row["note"] = None
        if bl is None:
            row["note"] = "no ballistic data"
        elif bl["lam"] is None:
            row["note"] = bl["note"]
        elif b is None:
            row["note"] = "no DCPA"
        elif row["dt"] is None:
            row["note"] = f"{row['found']}: no range from Dt"
        else:
            row["D_calc"] = range_from_dt(row["dt"], b, bl["v0"], bl["lam"], c)
            if D_ref:
                row["dt_exp"] = dt_model(D_ref, b, bl["v0"], bl["lam"], c)[0]
                row["v0_fit"] = v0_to_match(row["dt"], D_ref, b, bl["lam"], c)
                if row["D_calc"]:
                    row["err_pct"] = 100 * (row["D_calc"] - D_ref) / D_ref
        if dup_of:
            row["note"] = (row["note"] + "; " if row["note"] else "") + f"no own entry, uses {dup_of.split('_')[0]}"
        rows.append(row)
        overview.append((name, x, fs, det, row))
        f = lambda v, fmt: "--" if v is None else fmt % v
        print(f"{name:28s} {row['found']:8s} SW {f(row['t_sw'], '%.4f')}  MB {f(row['t_mb'], '%.4f')}  "
              f"Dt {f(row['dt'], '%.4f')}  "
              f"D {f(row['D_calc'], '%.0f')} m (ref {f(D_ref, '%.0f')}, err {f(row['err_pct'], '%+.1f')}%)"
              + (f"  [{row['note']}]" if row["note"] else ""))

        if args.plot:
            plot_one(os.path.join(out_dir, os.path.splitext(name)[0] + "_range.png"), name, x, fs, det, row, c)

    write_results(out_dir, rows, c)
    plot_overview(os.path.join(out_dir, "sw_mb_range_overview.png"), overview, c)
    print(f"\nSaved: {os.path.join(out_dir, 'sw_mb_range_results.md')} (+ .csv, overview png)")


def write_results(out_dir, rows, c):
    def f(v, fmt):
        return "–" if v is None else fmt % v
    with open(os.path.join(out_dir, "sw_mb_range_results.md"), "w", encoding="utf-8") as fh:
        fh.write("# SW/MB delay and range check\n\n"
                 f"Generated by `python/ana_sw_mb_range.py` (c = {c:g} m/s). Dt = muzzle-blast onset − "
                 "shockwave onset; D = range solved from Dt with the drag-aware model, using the reference "
                 "DCPA as miss distance. *v0 to match* = muzzle velocity (same decay length) that would "
                 "make D equal the reference.\n\n"
                 "| File | Weapon | Ballistic model | Found | SW (s) | MB (s) | Dt (s) | Dt expected (s) "
                 "| D computed (m) | D reference (m) | Error | v0 used (m/s) | v0 to match (m/s) | Note |\n"
                 "|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n")
        for r in rows:
            bl = BALLISTIC.get(r["key"]) if r["key"] else None
            fh.write(f"| {r['file'].split('_')[0]} | {r['weapon']} | {r['key'] or '–'} | {r['found']} | "
                     f"{f(r['t_sw'], '%.4f')} | {f(r['t_mb'], '%.4f')} | {f(r['dt'], '%.4f')} | {f(r['dt_exp'], '%.4f')} | "
                     f"{f(r['D_calc'], '%.0f')} | {f(r['D_ref'], '%.0f')} | {f(r['err_pct'], '%+.1f %%')} | "
                     f"{bl['v0'] if bl else '–'} | {f(r['v0_fit'], '%.0f')} | {r['note'] or ''} |\n")
    with open(os.path.join(out_dir, "sw_mb_range_results.csv"), "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["file", "weapon", "ballistic", "D_ref_m", "DCPA_m", "found", "t_sw_s", "t_mb_s", "mb_band_hz",
                     "dt_s", "dt_expected_s", "D_calc_m", "err_pct", "v0_to_match_mps", "note"])
        for r in rows:
            wr.writerow([r["file"], r["weapon"], r["key"], r["D_ref"], r["b"], r["found"],
                         "" if r["t_sw"] is None else f"{r['t_sw']:.4f}",
                         "" if r["t_mb"] is None else f"{r['t_mb']:.4f}",
                         "" if r["band"] is None else f"{r['band'][0]}-{r['band'][1]}",
                         "" if r["dt"] is None else f"{r['dt']:.4f}",
                         "" if r["dt_exp"] is None else f"{r['dt_exp']:.4f}",
                         "" if r["D_calc"] is None else f"{r['D_calc']:.1f}",
                         "" if r["err_pct"] is None else f"{r['err_pct']:.1f}",
                         "" if r["v0_fit"] is None else f"{r['v0_fit']:.0f}", r["note"] or ""])


def plot_one(path, name, x, fs, det, row, c):
    t = np.arange(len(x)) / fs
    fig, ax = plt.subplots(2, 1, figsize=(12, 5), sharex=True)
    lo, hi = det["band"] or (100, 500)
    for a, (blo, bhi, sm, lab) in zip(ax, ((1000, 4000, 0.0005, "SW band 1-4 kHz"),
                                          (lo, hi, 0.005, f"MB band {lo}-{hi} Hz"))):
        e = envelope(bandpass(x, fs, blo, bhi), fs, sm)
        a.semilogy(t, e + 1e-12, lw=0.6, color="k")
        a.set_ylabel(lab, fontsize=8)
        if det["t_sw"] is not None:
            a.axvline(det["t_sw"], color="red", ls="--", lw=1, label="SW")
        if det["t_mb"] is not None:
            a.axvline(det["t_mb"], color="green", ls="--", lw=1, label="MB")
        if row["dt_exp"] is not None and det["t_sw"] is not None:
            a.axvline(det["t_sw"] + row["dt_exp"], color="orange", ls=":", lw=1.5, label="MB expected (ref)")
        a.grid(alpha=0.3)
    ax[0].legend(fontsize=7, loc="upper right")
    ax[1].set_xlabel("Time (s)")
    fig.suptitle(f"{name}: Dt {row['dt'] if row['dt'] is None else round(row['dt'], 4)} s, "
                 f"D {row['D_calc'] if row['D_calc'] is None else round(row['D_calc'])} m "
                 f"(ref {row['D_ref']} m)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def plot_overview(path, items, c):
    n = len(items)
    fig, axes = plt.subplots(n, 1, figsize=(14, 1.35 * n + 1), squeeze=False)
    for ax, (name, x, fs, det, row) in zip(axes[:, 0], items):
        lo, hi = det["band"] or (100, 500)
        ref_t = det["t_sw"] if det["t_sw"] is not None else (det["t_mb"] or 0.0)
        t0 = ref_t - 0.2
        t1 = max(det["t_mb"] or 0, ref_t + (row["dt_exp"] or 0)) + 0.4
        i0, i1 = max(0, int(t0 * fs)), min(len(x), int(t1 * fs))
        t = np.arange(i0, i1) / fs
        e_sw = envelope(bandpass(x[i0:i1], fs, 1000, 4000), fs, 0.0005)
        e_mb = envelope(bandpass(x[i0:i1], fs, lo, hi), fs, 0.005)
        ax.semilogy(t, e_sw / np.median(e_sw), lw=0.5, color="#F44336", label="SW band")
        ax.semilogy(t, e_mb / np.median(e_mb), lw=0.8, color="#2196F3", label="MB band")
        if det["t_sw"] is not None:
            ax.axvline(det["t_sw"], color="red", ls="--", lw=1)
        if det["t_mb"] is not None:
            ax.axvline(det["t_mb"], color="green", ls="--", lw=1.2)
        if row["dt_exp"] is not None and det["t_sw"] is not None:
            ax.axvline(det["t_sw"] + row["dt_exp"], color="orange", ls=":", lw=1.5)
        err = "" if row["err_pct"] is None else f", err {row['err_pct']:+.1f}%"
        dcalc = "--" if row["D_calc"] is None else f"{row['D_calc']:.0f}"
        ax.text(0.005, 0.85, f"{name.split('_')[0]}  {row['weapon']}  [{row['found']}]  "
                f"D {dcalc} m (ref {row['D_ref']:.0f}{err})" if row["D_ref"] else name,
                transform=ax.transAxes, fontsize=8, va="top",
                bbox=dict(fc="white", ec="none", alpha=0.8))
        ax.set_ylim(0.5, None)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=7, loc="upper right", ncol=2)
    axes[-1, 0].set_xlabel("Time (s)")
    fig.suptitle(f"SW (red dashed), MB detected (green dashed), MB expected from reference (orange dotted); "
                 f"c = {c:g} m/s", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


if __name__ == "__main__":
    main()
