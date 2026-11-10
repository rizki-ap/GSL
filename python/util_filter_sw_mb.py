#!/usr/bin/env python3
"""
util_filter_sw_mb.py
====================
Separates the two parts of a supersonic gunshot by frequency, to find them
in a recording:
  - shockwave (SW, the bullet's crack): short, high-frequency  -> 1-4 kHz band
  - muzzle blast (MB, the gun's bang):  longer, low-frequency  -> 100-500 Hz band

Each band is filtered with a zero-phase Butterworth band-pass (sosfiltfilt,
so arrival times are NOT shifted), and its envelope is compared against the
band's own noise floor. The SW is the strongest SW-band event in the file;
the MB is the strongest MB-band event AFTER the SW (from -guard to -dt_max
seconds after it, so the SW's own tail is skipped).
If no SW is found above -snr, the strongest MB-band event in the whole file
is reported instead.

Output
------
  <name>_swmb.png     raw signal + both bands with envelopes, full length and
                      zoomed on SW and MB (opened in a viewer unless --no-show)
  <name>_sw.wav, <name>_mb.wav
                      filtered signals, float32, only with -save_wav

Usage
-----
  python util_filter_sw_mb.py shot.wav
  python util_filter_sw_mb.py shot.wav -sw 1000 4000 -mb 100 500     # the defaults
  python util_filter_sw_mb.py shot.wav -chan 0 1 3 5                 # channels to use
  python util_filter_sw_mb.py shot.wav -dt_max 2 -snr 10
  python util_filter_sw_mb.py shot.wav -save_wav --no-show
"""

import argparse
import os
import subprocess
import sys

import numpy as np
import soundfile as sf
from scipy.signal import butter, hilbert, sosfiltfilt

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS = ['#2196F3', '#F44336', '#4CAF50', '#FF9800', '#9C27B0', '#00BCD4',
          '#795548', '#607D8B']


def open_image(path):
    try:
        if sys.platform.startswith("win"):
            os.startfile(os.path.abspath(path))
        elif sys.platform == "darwin":
            subprocess.Popen(["open", path])
        else:
            subprocess.Popen(["xdg-open", path])
    except Exception as e:
        print(f"Could not open {path}: {e}")


def bandpass(x, fs, lo, hi, order=4):
    sos = butter(order, [lo, hi], btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, x, axis=0)


def envelope(x, fs, smooth_s):
    """Hilbert envelope per channel, smoothed with a moving average."""
    env = np.abs(hilbert(x, axis=0))
    w = max(1, int(smooth_s * fs))
    k = np.ones(w) / w
    return np.stack([np.convolve(env[:, i], k, mode="same") for i in range(env.shape[1])], axis=1)


def noise_floor(env):
    """Robust noise level of an envelope (median), so short events don't bias it."""
    return float(np.median(env)) + 1e-20


def find_event(env_sum, fs, noise, i0=0, i1=None, snr_db=10.0):
    """Strongest envelope peak in [i0, i1); onset = first sample before it that
    rises above snr_db over the noise. Returns dict or None."""
    i1 = len(env_sum) if i1 is None else min(i1, len(env_sum))
    if i1 - i0 < 2:
        return None
    seg = env_sum[i0:i1]
    p = i0 + int(np.argmax(seg))
    snr = 20 * np.log10(env_sum[p] / noise)
    if snr < snr_db:
        return dict(peak_idx=p, onset_idx=None, snr_db=snr, ok=False)
    thr = noise * 10 ** (snr_db / 20)
    below = np.flatnonzero(env_sum[i0:p] < thr)
    onset = i0 + int(below[-1]) + 1 if below.size else i0
    return dict(peak_idx=p, onset_idx=onset, snr_db=snr, ok=True)


def per_channel_onsets(env, noise_ch, center, fs, snr_db, half_win_s):
    """Onset per channel near `center` (first rise above snr_db over that channel's noise)."""
    out = []
    a = max(0, center - int(half_win_s * fs))
    b = min(len(env), center + int(half_win_s * fs))
    for i in range(env.shape[1]):
        thr = noise_ch[i] * 10 ** (snr_db / 20)
        hits = np.flatnonzero(env[a:b, i] >= thr)
        out.append((a + int(hits[0])) / fs if hits.size else None)
    return out


def main():
    p = argparse.ArgumentParser(description="Locate shockwave and muzzle blast by band-pass filtering.",
                                add_help=False)
    p.add_argument("-h", "-help", "--help", action="help", help="show this help message and exit")
    p.add_argument("wav_path")
    p.add_argument("-sw", type=float, nargs=2, default=[1000.0, 4000.0], metavar=("LO", "HI"),
                   help="shockwave band in Hz (default: 1000 4000)")
    p.add_argument("-mb", type=float, nargs=2, default=[100.0, 500.0], metavar=("LO", "HI"),
                   help="muzzle-blast band in Hz (default: 100 500)")
    p.add_argument("-chan", type=int, nargs="+", default=None, metavar="CH",
                   help="channels to use (0-based, default: all)")
    p.add_argument("-dt_max", type=float, default=3.0, metavar="SEC",
                   help="search the MB up to this many seconds after the SW (default: 3)")
    p.add_argument("-guard", type=float, default=0.05, metavar="SEC",
                   help="skip this long after the SW peak before searching the MB, so the SW's "
                        "own tail/reflections aren't taken as the MB (default: 0.05)")
    p.add_argument("-snr", type=float, default=10.0, metavar="DB",
                   help="minimum envelope level above the band's noise floor to count "
                        "as an event (default: 10 dB)")
    p.add_argument("-zoom", type=float, default=0.05, metavar="SEC",
                   help="half-width of the zoom plots (default: 0.05)")
    p.add_argument("-o", default=None, metavar="PNG", help="output image (default: <name>_swmb.png)")
    p.add_argument("-save_wav", action="store_true", help="also save the filtered signals as wav")
    p.add_argument("--no-show", action="store_true", help="don't open the image")
    args = p.parse_args()

    x, fs = sf.read(args.wav_path, dtype="float64", always_2d=True)
    n_ch = x.shape[1]
    chans = list(range(n_ch)) if args.chan is None else list(dict.fromkeys(args.chan))
    bad = [c for c in chans if not 0 <= c < n_ch]
    if bad:
        p.error(f"-chan {bad} out of range: file has channels 0..{n_ch - 1}")
    for name, (lo, hi) in (("-sw", args.sw), ("-mb", args.mb)):
        if not 0 < lo < hi < fs / 2:
            p.error(f"{name} band must satisfy 0 < LO < HI < {fs / 2:g} Hz")
    x = x[:, chans]
    t = np.arange(len(x)) / fs
    print(f"{args.wav_path}: {n_ch} ch @ {fs} Hz, {len(x) / fs:.3f} s, using ch {chans}")

    sw = bandpass(x, fs, *args.sw)
    mb = bandpass(x, fs, *args.mb)
    env_sw = envelope(sw, fs, 0.0005)      # SW is ~0.2-0.5 ms long: light smoothing
    env_mb = envelope(mb, fs, 0.005)       # MB lasts several ms: heavier smoothing
    sum_sw, sum_mb = env_sw.mean(axis=1), env_mb.mean(axis=1)
    n_sw, n_mb = noise_floor(sum_sw), noise_floor(sum_mb)
    nch_sw = [noise_floor(env_sw[:, i]) for i in range(len(chans))]
    nch_mb = [noise_floor(env_mb[:, i]) for i in range(len(chans))]

    ev_sw = find_event(sum_sw, fs, n_sw, snr_db=args.snr)
    if ev_sw and ev_sw["ok"]:
        i0 = ev_sw["peak_idx"] + int(args.guard * fs)   # skip the SW and its tail
        ev_mb = find_event(sum_mb, fs, n_mb, i0, ev_sw["peak_idx"] + int(args.dt_max * fs), args.snr)
    else:
        ev_mb = find_event(sum_mb, fs, n_mb, snr_db=args.snr)

    def report(label, ev, env, noise_ch, half_win):
        if ev is None or not ev["ok"]:
            snr = f"{ev['snr_db']:.1f} dB" if ev else "n/a"
            print(f"{label}: NOT FOUND (strongest candidate only {snr} above noise)")
            return None
        on = per_channel_onsets(env, noise_ch, ev["peak_idx"], fs, args.snr, half_win)
        on_txt = ", ".join(f"ch{c}={v:.5f}" if v is not None else f"ch{c}=--" for c, v in zip(chans, on))
        print(f"{label}: onset {ev['onset_idx'] / fs:.5f} s, peak {ev['peak_idx'] / fs:.5f} s, "
              f"{ev['snr_db']:.1f} dB above noise")
        print(f"    per-channel onset (s): {on_txt}")
        return ev

    print()
    ev_sw = report(f"SW  ({args.sw[0]:g}-{args.sw[1]:g} Hz)", ev_sw, env_sw, nch_sw, 0.005)
    ev_mb = report(f"MB  ({args.mb[0]:g}-{args.mb[1]:g} Hz)", ev_mb, env_mb, nch_mb, 0.02)
    if ev_sw and ev_mb:
        dt = (ev_mb["onset_idx"] - ev_sw["onset_idx"]) / fs
        print(f"\nSW -> MB delay: {dt:.4f} s")

    # ---------------- plot ----------------
    rows = [("Raw", x, None, None), (f"SW band {args.sw[0]:g}-{args.sw[1]:g} Hz", sw, sum_sw, ev_sw),
            (f"MB band {args.mb[0]:g}-{args.mb[1]:g} Hz", mb, sum_mb, ev_mb)]
    fig, axes = plt.subplots(3, 3, figsize=(18, 9), squeeze=False,
                             gridspec_kw=dict(width_ratios=[2.2, 1, 1]))
    zooms = [("SW zoom", ev_sw), ("MB zoom", ev_mb)]
    for r, (title, sig, env_sum, ev) in enumerate(rows):
        for c in range(3):
            ax = axes[r, c]
            if c == 0:
                mask = slice(None)
            else:
                zev = zooms[c - 1][1]
                if zev is None:
                    ax.text(0.5, 0.5, f"{zooms[c - 1][0]}: not found", ha="center", va="center",
                            transform=ax.transAxes)
                    ax.set_xticks([]); ax.set_yticks([])
                    continue
                ctr = zev["peak_idx"]
                h = int(args.zoom * fs)
                mask = slice(max(0, ctr - h), min(len(t), ctr + h))
            for i, ch in enumerate(chans):
                ax.plot(t[mask], sig[mask, i], lw=0.5 if c == 0 else 0.8,
                        color=COLORS[ch % len(COLORS)], label=f"Ch {ch}" if (r == 0 and c == 0) else None)
            if env_sum is not None:
                ax.plot(t[mask], env_sum[mask], color="k", lw=1.0, alpha=0.8,
                        label="envelope (mean)" if c == 0 else None)
            for lab, e, col in (("SW", ev_sw, "red"), ("MB", ev_mb, "green")):
                if e is not None and t[mask][0] <= e["onset_idx"] / fs <= t[mask][-1]:
                    ax.axvline(e["onset_idx"] / fs, color=col, ls="--", lw=1)
                    if c == 0:
                        ax.text(e["onset_idx"] / fs, 0.95, f" {lab}", color=col, fontsize=9,
                                transform=ax.get_xaxis_transform(), va="top", fontweight="bold")
            ax.grid(alpha=0.3)
            if c == 0:
                ax.set_ylabel(title, fontsize=9, fontweight="bold")
            if r == 0:
                ax.set_title("Full length" if c == 0 else zooms[c - 1][0], fontsize=10, fontweight="bold")
            if r == 2:
                ax.set_xlabel("Time (s)")
            if c == 0 and r in (0, 1):
                ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.0, 1), borderaxespad=0.2)
    dt_txt = (f"  |  SW->MB {(ev_mb['onset_idx'] - ev_sw['onset_idx']) / fs:.3f} s"
              if ev_sw and ev_mb else "")
    fig.suptitle(f"{args.wav_path}  |  {fs} Hz, ch {chans}{dt_txt}", fontsize=11, fontweight="bold")
    fig.tight_layout()
    out_png = args.o or os.path.splitext(args.wav_path)[0] + "_swmb.png"
    try:
        fig.savefig(out_png, dpi=130, bbox_inches="tight")
    except OSError:
        out_png = os.path.splitext(out_png)[0] + "_1.png"
        fig.savefig(out_png, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved: {out_png}")

    if args.save_wav:
        stem = os.path.splitext(args.wav_path)[0]
        for tag, sig in (("sw", sw), ("mb", mb)):
            path = f"{stem}_{tag}.wav"
            sf.write(path, sig.astype(np.float32), fs, subtype="FLOAT")
            print(f"Saved: {path}")

    if not args.no_show:
        open_image(out_png)


if __name__ == "__main__":
    main()
