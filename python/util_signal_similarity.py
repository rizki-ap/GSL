#!/usr/bin/env python3
"""
util_signal_similarity.py
=========================
Similarity between MEASURED and SYNTHESIZED gunshot recordings, following the
evaluation of publication/Paper-Generator.md, Section 4.3:

  4.3.1 waveform morphology   Pearson r, SNR_synth, peak error, duration error
                              (T_N of the N-wave / positive-phase T_mb of the blast,
                              MEASURED from each signal)
  4.3.2 frequency domain      Tukey-window event PSD (4096-point zero padding), averaged
                              over channels and recordings; LSD; blast peak frequency
                              (Friedlander f_pk = 1 / (2 pi T_mb)); N-wave first null
                              (f_null ~ 1.43 / T_N); ensemble (pooled) coherence and its
                              > 0.8 bandwidth -- reported ONLY when enough segments are
                              pooled (default K >= 8), otherwise it is trivially ~1
  4.3.3 timing / TDOA         inter-channel TDOAs of each event (GCC-PHAT or plain
                              cross-correlation, sinc-upsampled peak) of the synthesized
                              signal vs the measured one (or vs reference TDOAs):
                              MAE, RMSE, Bias over all pairs i < j, plus the N-1
                              independent pairs (reference channel 0)
  4.3.4 shock-to-blast delay  per-channel delay t_mb - t_sw, synthesized vs measured

Both signals are high-passed (2nd-order, zero-phase) at the lower edge of each event's band
before windowing -- the generator does not model wind, recoil or drift (--no-highpass).
Event windows are found automatically with the detector's own classifier
(gs_det_events.detect_sw_mb) or given by hand. Amplitudes are never normalized
(the paper evaluates amplitude fidelity); give both files in the same units.

Usage
-----
  python util_signal_similarity.py measured.wav synthetic.wav
  python util_signal_similarity.py meas.wav synth.wav --event nwave --channels 1 2 3 4
  python util_signal_similarity.py --pair m1.wav s1.wav --pair m2.wav s2.wav --json table.json
  python util_signal_similarity.py meas.wav synth.wav --meas-sw-ms 1699.3 --synth-sw-ms 50.2
  python util_signal_similarity.py meas.wav synth.wav --plot out_dir

Options (main ones; see --help)
  --event nwave|blast|both   which events to evaluate (default both)
  --channels 1 2 ...         1-based channels (default: all channels common to both files)
  --nwave-window-ms 1.5      N-wave window length (starts 0.25 ms before the onset)
  --blast-window-ms 20       blast window length (starts 2 ms before the onset)
  --nwave-band 200 F         N-wave band (default 200 Hz to min(20 kHz, 0.45 fs): the null
                             1.43/T_N lies near 5-7 kHz for T_N = 0.2-0.3 ms)
  --blast-band 50 4000       blast band (paper Table 2)
  --min-segments 8           coherence is reported only with at least this many segments
  --tdoa gcc|xcorr           TDOA estimator (default xcorr: GCC-PHAT fails on recordings
                             with common-mode noise, e.g. Cipatat)
  --synth-json GEN.json      generator JSON (repeat per pair): synthetic event times per
                             channel and model TDOA labels; TDOAs are then reported as
                             MEASURED vs model labels (paper 4.3.3)
  --meas-scale / --synth-scale  multiply the samples (e.g. --meas-scale 1520: PILAR-V FS to Pa;
                             the generator's *_clean.wav is already in Pa)
  --ref-tdoa FILE.json       reference TDOAs {"nwave": {"0,1": s, ...}, "blast": {...}}

As a library:
  from util_signal_similarity import load_wav, evaluate
  per_pair, summary = evaluate([(meas_x, meas_fs, synth_x, synth_fs)], event="both")
"""

import argparse
import json
import os
import sys
import warnings
from math import gcd

import numpy as np
from scipy.signal import (butter, correlate, find_peaks, resample_poly, savgol_filter,
                          sosfiltfilt, welch)
from scipy.signal.windows import tukey

_trapz = getattr(np, "trapezoid", None) or np.trapz          # numpy 2.x / 1.x


# ===========================================================================
# I/O
# ===========================================================================
def load_wav(path, channels=None):
    """(fs, x[n, ch]) as float64 in the file's own scale (integer PCM -> +-1 full scale,
    float WAV unchanged). No normalization. channels: 0-based list or None (all)."""
    try:
        import soundfile as sf
        x, fs = sf.read(path, dtype="float64", always_2d=True)
    except ImportError:
        from scipy.io import wavfile
        fs, x = wavfile.read(path)
        if np.issubdtype(x.dtype, np.integer):
            x = x.astype(np.float64) / float(np.iinfo(x.dtype).max)
        x = np.atleast_2d(np.asarray(x, dtype=np.float64))
        if x.shape[0] < x.shape[1]:
            x = x.T
    if channels is not None:
        x = x[:, [c for c in channels if c < x.shape[1]]]
    return int(fs), x


def to_rate(x, fs, fs_target):
    if int(fs) == int(fs_target):
        return x
    g = gcd(int(fs), int(fs_target))
    return resample_poly(x, int(fs_target) // g, int(fs) // g, axis=0)


def bandpass(x, fs, lo, hi, order=4):
    hi = min(hi, 0.45 * fs)
    sos = butter(order, [lo, hi], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, x, axis=0)


# ===========================================================================
# Events: automatic windows
# ===========================================================================
def detect_events(x, fs):
    """Array-level (t_sw, t_mb) in s from gs_det_events (None when absent)."""
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import gs_det_events as ge
    except ImportError:
        warnings.warn("gs_det_events not found -- give the event times by hand")
        return None, None
    r = ge.detect_sw_mb(x - x.mean(axis=0), fs)
    return r.get("t_sw"), r.get("t_mb")


def channel_onset(x1, fs, t0, kind, search_s=0.0015):
    """Onset (s) of one channel near t0: first sample above 30 % of the local peak of the
    band-passed signal (N-wave 1-8 kHz, blast 150-2000 Hz) within t0 +- search_s."""
    band = (1000.0, 8000.0) if kind == "nwave" else (150.0, 2000.0)
    y = np.abs(bandpass(x1, fs, *band))
    i0, i1 = max(0, int((t0 - search_s) * fs)), min(len(y), int((t0 + search_s) * fs))
    if i1 - i0 < 3:
        return t0
    seg = y[i0:i1]
    k = int(np.argmax(seg > 0.3 * seg.max()))
    return (i0 + k) / fs


def cut(x1, fs, t_start, dur):
    i0 = int(round(t_start * fs)); n = int(round(dur * fs))
    if i0 < 0 or i0 + n > len(x1) or n < 4:
        return None
    return np.asarray(x1[i0:i0 + n], dtype=np.float64).copy()


# ===========================================================================
# Pulse parameters measured from the signal
# ===========================================================================
def polarity(w):
    """+1 if the first major excursion is positive, else -1."""
    a = np.abs(w)
    if a.max() <= 0:
        return 1
    k = int(np.argmax(a > 0.5 * a.max()))
    return 1 if w[k] >= 0 else -1


def nwave_duration(w, fs):
    """T_N (s): front shock (steepest rise) to rear shock (steepest rise after the
    minimum), on the polarity-corrected window."""
    w = polarity(w) * np.asarray(w)
    d = np.diff(w)
    if len(d) < 5:
        return np.nan
    dd = d[: max(5, int(0.6 * len(d)))]
    i_front = int(np.argmax(dd > 0.5 * dd.max()))          # EARLIEST steep rise (the rear shock is as steep)
    i_front += int(np.argmax(d[i_front:i_front + 4]))
    # the minimum must belong to THIS N-wave (small arms: T_N < ~0.8 ms), not to a later
    # ground reflection or echo, which on real recordings follows ~0.5-3 ms after the SW
    i_end = min(len(w), i_front + max(4, int(0.0008 * fs)))
    i_min = i_front + int(np.argmin(w[i_front:i_end]))
    if i_min >= len(d) - 1:
        return np.nan
    i_rear = i_min + int(np.argmax(d[i_min:i_min + max(3, int(0.00015 * fs))]))
    return (i_rear - i_front) / fs if i_rear > i_front else np.nan


def blast_duration(w, fs):
    """Positive-phase duration T_mb (s): onset (10 % of the peak, before it) to the first
    zero crossing after the peak, on the polarity-corrected window (Friedlander T)."""
    w = polarity(w) * np.asarray(w)
    ip = int(np.argmax(w))
    if w[ip] <= 0:
        return np.nan
    pre = np.where(w[:ip + 1] < 0.1 * w[ip])[0]
    i0 = int(pre[-1]) + 1 if len(pre) else 0
    zc = np.where(w[ip:] <= 0)[0]
    return (ip + int(zc[0]) - i0) / fs if len(zc) else np.nan


# ===========================================================================
# 4.3.1 waveform morphology
# ===========================================================================
def best_lag(x, y, fs, max_lag_s):
    """(lag_s, signed peak correlation): lag > 0 means y arrives after x."""
    n = min(len(x), len(y))
    a, b = x[:n] - x[:n].mean(), y[:n] - y[:n].mean()
    c = correlate(a, b, mode="full")
    lags = np.arange(-n + 1, n)
    keep = np.abs(lags) <= int(round(max_lag_s * fs))
    c, lags = c[keep], lags[keep]
    k = int(np.argmax(np.abs(c)))
    return float(-lags[k] / fs), float(np.sign(c[k]))


def align(x, y, fs, lag_s):
    k = int(round(lag_s * fs))
    a, b = (x[:len(x) - k], y[k:]) if k > 0 else ((x[-k:], y[:len(y) + k]) if k < 0 else (x, y))
    n = min(len(a), len(b))
    return np.asarray(a[:n]), np.asarray(b[:n])


def pearson_r(x, y):
    x0, y0 = x - x.mean(), y - y.mean()
    den = np.sqrt(np.sum(x0 * x0) * np.sum(y0 * y0))
    return float(np.sum(x0 * y0) / den) if den > 1e-30 else np.nan


def snr_synth_db(x, y):
    num, den = np.sum(x * x), np.sum((x - y) ** 2)
    return float(10 * np.log10(num / den)) if num > 1e-30 and den > 1e-30 else np.nan


def pct_err(meas, synth):
    return float(abs(synth - meas) / abs(meas) * 100) if np.isfinite(meas) and meas != 0 and np.isfinite(synth) else np.nan


# ===========================================================================
# 4.3.2 frequency domain
# ===========================================================================
def event_psd(w, fs, nfft=4096, tukey_alpha=0.25):
    w = (w - w.mean()) * tukey(len(w), alpha=tukey_alpha)
    return welch(w, fs=fs, window="boxcar", nperseg=len(w), noverlap=0,
                 nfft=max(int(nfft), len(w)), detrend=False, scaling="density")


def mean_psd(windows, fs, nfft, tukey_alpha):
    sp = [event_psd(w, fs, nfft, tukey_alpha) for w in windows]
    return sp[0][0], np.mean([s[1] for s in sp], axis=0)


def lsd_db(f, Pm, Ps, band):
    k = (f >= band[0]) & (f <= band[1]) & (Pm > 0) & (Ps > 0)
    if k.sum() < 3:
        return np.nan
    d = 10 * np.log10(Ps[k] / Pm[k])
    return float(np.sqrt(_trapz(d * d, f[k]) / (f[k][-1] - f[k][0])))


def smooth_log(P, window=11):
    w = min(window | 1, len(P) - (1 - len(P) % 2))
    return 10 ** savgol_filter(np.log10(np.maximum(P, 1e-30)), w, 2) if w >= 5 else P


def peak_freq(f, P, band, window=11):
    k = (f >= band[0]) & (f <= band[1])
    return float(f[k][np.argmax(smooth_log(P[k], window))]) if k.sum() >= 3 else np.nan


def null_freq(f, P, T_N, band, window=11):
    """First non-trivial N-wave null: minimum of the smoothed spectrum within 0.6-1.6 x
    the expected 1.43 / T_N (T_N measured from the same signal)."""
    if not np.isfinite(T_N) or T_N <= 0:
        return np.nan
    fe = 1.43 / T_N
    k = (f >= max(band[0], 0.6 * fe)) & (f <= min(band[1], 1.6 * fe))
    if k.sum() < 3:
        return np.nan
    return float(f[k][np.argmin(smooth_log(P[k], window))])


def pooled_coherence(pairs, fs, nperseg, nfft):
    """Ensemble MSC |sum Sxy|^2 / (sum Sxx sum Syy) over all segments of all pairs (paper
    eq. of 4.3.2). Returns (f, C, K)."""
    sxy = sx = sy = None
    K = 0
    win = tukey(nperseg, alpha=0.25)
    for x, y in pairs:
        n = min(len(x), len(y))
        for i in range(0, n - nperseg + 1, max(1, nperseg // 2)):
            a = (x[i:i + nperseg] - x[i:i + nperseg].mean()) * win
            b = (y[i:i + nperseg] - y[i:i + nperseg].mean()) * win
            X, Y = np.fft.rfft(a, nfft), np.fft.rfft(b, nfft)
            if sxy is None:
                sxy = np.zeros(len(X), complex); sx = np.zeros(len(X)); sy = np.zeros(len(X))
            sxy += X * np.conj(Y); sx += np.abs(X) ** 2; sy += np.abs(Y) ** 2
            K += 1
    if K == 0:
        return None, None, 0
    f = np.fft.rfftfreq(nfft, 1 / fs)
    C = np.where(sx * sy > 1e-30, np.abs(sxy) ** 2 / np.maximum(sx * sy, 1e-30), 0.0)
    return f, np.clip(C, 0, 1), K


def coherence_bandwidth(f, C, band, thr=0.8):
    k = (f >= band[0]) & (f <= band[1])
    ff, cc = f[k], C[k]
    return float(np.sum(np.diff(ff)[(cc[:-1] > thr) & (cc[1:] > thr)])) if len(ff) > 1 else np.nan


# ===========================================================================
# 4.3.3 TDOA
# ===========================================================================
def tdoa(a, b, fs, max_lag_s, method="xcorr", up=16):
    """t_a - t_b (s) by plain cross-correlation or GCC-PHAT, sinc-upsampled peak."""
    n = 1 << int(np.ceil(np.log2(2 * len(a))))
    A, B = np.fft.rfft(a - a.mean(), n), np.fft.rfft(b - b.mean(), n)
    G = A * np.conj(B)
    if method == "gcc":
        G = G / np.maximum(np.abs(G), 1e-20)
    cc = np.fft.irfft(G, n * up)                       # sinc interpolation by zero padding
    m = int(round(max_lag_s * fs * up))
    cc = np.concatenate([cc[-m:], cc[:m + 1]])
    return (int(np.argmax(cc)) - m) / (fs * up)


def tdoa_table(windows, fs, max_lag_s, method):
    out = {}
    for i in range(len(windows)):
        for j in range(i + 1, len(windows)):
            out[f"{i},{j}"] = tdoa(windows[i], windows[j], fs, max_lag_s, method)
    return out


def tdoa_errors(t_synth, t_ref):
    keys = [k for k in t_synth if k in t_ref]
    e = np.array([t_synth[k] - t_ref[k] for k in keys])
    ind = np.array([t_synth[k] - t_ref[k] for k in keys if k.startswith("0,")])
    stat = lambda v: dict(mae_us=float(np.mean(np.abs(v)) * 1e6), rmse_us=float(np.sqrt(np.mean(v ** 2)) * 1e6),
                          bias_us=float(np.mean(v) * 1e6), n=int(len(v))) if len(v) else None
    return dict(all_pairs=stat(e), independent_ref0=stat(ind), per_pair_us={k: float(v * 1e6) for k, v in zip(keys, e)})


# ===========================================================================
# Generator ground truth (gs_gen_* JSON)
# ===========================================================================
def read_generator_json(path):
    """Event times inside the synthetic WAV and model TDOA labels from the generator's
    JSON (<name>.json of gs_gen_apply_adc.py / <name>_clean.json). Returns
    (times dict, reference TDOAs dict) -- times in s from the start of the WAV."""
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    g = d.get("clean_signal_ground_truth", d)
    t0 = g.get("t_master_start_s", 0.0)
    per = g.get("per_mic", [])
    sw = [m["shockwave"]["t_arrive_s"] - t0 if m.get("shockwave", {}).get("generated") else None for m in per]
    mb = [m["muzzle_blast"]["t_arrive_s"] - t0 if m.get("muzzle_blast", {}).get("generated") else None for m in per]
    times = {}
    if per and all(v is not None for v in sw):
        times.update(synth_sw=min(sw), synth_sw_ch=sw)
    if per and all(v is not None for v in mb):
        times.update(synth_mb=min(mb), synth_mb_ch=mb)
    ref = {}
    for kind, key in (("nwave", "tdoa_matrix_sw_us"), ("blast", "tdoa_matrix_mb_us")):
        Mx = g.get(key)
        if Mx:
            ref[kind] = {f"{i},{j}": Mx[i][j] * 1e-6 for i in range(len(Mx)) for j in range(i + 1, len(Mx))}
    return times, ref


# ===========================================================================
# Evaluation
# ===========================================================================
DEFAULTS = dict(nwave_window_ms=1.5, nwave_pre_ms=0.25, blast_window_ms=20.0, blast_pre_ms=2.0,
                nwave_band=None, blast_band=(50.0, 4000.0), nfft=4096, tukey_alpha=0.25,
                max_lag_ms=1.0, min_segments=8, coherence_threshold=0.8, smooth_window=11,
                tdoa_method="xcorr", tdoa_max_lag_ms=1.2, highpass=True)


def _events_for(x, fs, given_sw, given_mb):
    t_sw, t_mb = detect_events(x, fs) if (given_sw is None or given_mb is None) else (None, None)
    return (given_sw if given_sw is not None else t_sw), (given_mb if given_mb is not None else t_mb)


def evaluate(pairs, event="both", times=None, ref_tdoa=None, opts=None, verbose=True):
    """pairs: list of (x_meas[n, ch], fs_meas, x_synth[n, ch], fs_synth); times: optional list of
    dicts {meas_sw, meas_mb, synth_sw, synth_mb} in s. Returns (per_pair list, summary dict)."""
    o = dict(DEFAULTS, **(opts or {}))
    kinds = ["nwave", "blast"] if event == "both" else [event]
    pool = {k: dict(meas=[], synth=[], coh=[], T_meas=[], T_synth=[], rows=[], tdoa=[], fs=None) for k in kinds}
    delay_err = []
    per_pair = []

    for ip, (xm, fm, xs, fsn) in enumerate(pairs):
        fs = max(fm, fsn)
        xm, xs = to_rate(xm, fm, fs), to_rate(xs, fsn, fs)
        nch = min(xm.shape[1], xs.shape[1])
        xm, xs = xm[:, :nch], xs[:, :nch]
        t = (times[ip] if times else {}) or {}
        m_sw, m_mb = _events_for(xm, fs, t.get("meas_sw"), t.get("meas_mb"))
        s_sw, s_mb = _events_for(xs, fs, t.get("synth_sw"), t.get("synth_mb"))
        rec = dict(pair=ip, fs=fs, channels=nch, meas_sw_s=m_sw, meas_mb_s=m_mb, synth_sw_s=s_sw, synth_mb_s=s_mb)
        onsets = {}
        for kind in kinds:
            tm, ts = (m_sw, s_sw) if kind == "nwave" else (m_mb, s_mb)
            if tm is None or ts is None:
                rec[kind] = "not found in " + ("both" if tm is None and ts is None else ("measured" if tm is None else "synthetic"))
                continue
            L = o[f"{kind}_window_ms"] / 1e3; pre = o[f"{kind}_pre_ms"] / 1e3
            P = pool[kind]; P["fs"] = fs
            band = tuple(o["nwave_band"] or (200.0, min(20000.0, 0.45 * fs))) if kind == "nwave" else tuple(o["blast_band"])
            P["band"] = band
            wins_m, wins_s, rows = [], [], []
            if o["highpass"]:      # remove what the generator does not model (wind, recoil, drift)
                sos = butter(2, band[0], btype="high", fs=fs, output="sos")
                xmk, xsk = sosfiltfilt(sos, xm, axis=0), sosfiltfilt(sos, xs, axis=0)
            else:
                xmk, xsk = xm, xs
            for c in range(nch):
                om = channel_onset(xm[:, c], fs, tm, kind)
                given = t.get(f"synth_{'sw' if kind == 'nwave' else 'mb'}_ch")
                os_ = given[c] if given is not None and c < len(given) else channel_onset(xs[:, c], fs, ts, kind)
                onsets[(kind, c)] = (om, os_)
                wm, ws = cut(xmk[:, c], fs, om - pre, L), cut(xsk[:, c], fs, os_ - pre, L)
                if wm is None or ws is None:
                    continue
                wins_m.append(cut(xmk[:, c], fs, tm - pre - 0.0012, L + 0.0024))   # common window for TDOA
                wins_s.append(cut(xsk[:, c], fs, ts - pre - 0.0012, L + 0.0024))
                lag, sgn = best_lag(wm, ws, fs, o["max_lag_ms"] / 1e3)
                am, as_ = align(wm, ws, fs, lag)
                Tm = nwave_duration(wm, fs) if kind == "nwave" else blast_duration(wm, fs)
                Ts = nwave_duration(ws, fs) if kind == "nwave" else blast_duration(ws, fs)
                row = dict(channel=c + 1, lag_ms=lag * 1e3, pearson_r=pearson_r(am, as_), snr_synth_db=snr_synth_db(am, as_),
                           peak_meas=float(np.max(np.abs(am))), peak_synth=float(np.max(np.abs(as_))),
                           T_meas_ms=Tm * 1e3, T_synth_ms=Ts * 1e3,
                           polarity_meas=polarity(wm), polarity_synth=polarity(ws))
                row["peak_error_pct"] = pct_err(row["peak_meas"], row["peak_synth"])
                row["T_error_pct"] = pct_err(Tm, Ts)
                row["polarity_mismatch"] = row["polarity_meas"] != row["polarity_synth"] or sgn < 0
                rows.append(row)
                P["meas"].append(am); P["synth"].append(as_); P["coh"].append((am, as_))
                P["T_meas"].append(Tm); P["T_synth"].append(Ts)
            rec[kind] = rows
            P["rows"] += rows
            if len(wins_m) >= 2 and all(w is not None for w in wins_m + wins_s):
                lagmax = o["tdoa_max_lag_ms"] / 1e3
                t_m = tdoa_table(wins_m, fs, lagmax, o["tdoa_method"])
                t_s = tdoa_table(wins_s, fs, lagmax, o["tdoa_method"])
                t_lab = (t.get("ref_tdoa") or ref_tdoa or {}).get(kind)
                if t_lab:     # paper 4.3.3: measured TDOAs against the model labels
                    rec[f"{kind}_tdoa"] = dict(measured_s=t_m, synthetic_s=t_s, model_label_s=t_lab,
                                               errors_measured_vs_model=tdoa_errors(t_m, t_lab),
                                               errors_synthetic_vs_model=tdoa_errors(t_s, t_lab))
                    P["tdoa"].append((t_m, t_lab))
                else:         # no labels: synthetic TDOAs against the measured ones
                    rec[f"{kind}_tdoa"] = dict(measured_s=t_m, synthetic_s=t_s, errors=tdoa_errors(t_s, t_m))
                    P["tdoa"].append((t_s, t_m))
        if "nwave" in kinds and "blast" in kinds:
            for c in range(nch):
                if ("nwave", c) in onsets and ("blast", c) in onsets:
                    dm = onsets[("blast", c)][0] - onsets[("nwave", c)][0]
                    ds = onsets[("blast", c)][1] - onsets[("nwave", c)][1]
                    delay_err.append(ds - dm)
        per_pair.append(rec)

    summary = {}
    for kind in kinds:
        P = pool[kind]
        if not P["rows"]:
            summary[kind] = "no events"
            continue
        fs, band = P["fs"], P["band"]
        mean = lambda key: float(np.nanmean([r[key] for r in P["rows"]]))
        s = dict(n_windows=len(P["rows"]), band_hz=list(band),
                 pearson_r=mean("pearson_r"), snr_synth_db=mean("snr_synth_db"),
                 peak_error_pct=mean("peak_error_pct"), T_meas_ms=mean("T_meas_ms"), T_synth_ms=mean("T_synth_ms"),
                 T_error_pct=mean("T_error_pct"),
                 polarity_mismatches=int(sum(r["polarity_mismatch"] for r in P["rows"])))
        f, Pm = mean_psd(P["meas"], fs, o["nfft"], o["tukey_alpha"])
        _, Ps = mean_psd(P["synth"], fs, o["nfft"], o["tukey_alpha"])
        s["lsd_db"] = lsd_db(f, Pm, Ps, band)
        if kind == "blast":
            s["f_pk_meas_hz"], s["f_pk_synth_hz"] = peak_freq(f, Pm, band, o["smooth_window"]), peak_freq(f, Ps, band, o["smooth_window"])
            s["f_pk_error_pct"] = pct_err(s["f_pk_meas_hz"], s["f_pk_synth_hz"])
            s["f_pk_expected_from_T_hz"] = dict(meas=1 / (2 * np.pi * s["T_meas_ms"] / 1e3), synth=1 / (2 * np.pi * s["T_synth_ms"] / 1e3))
        else:
            Tm, Ts = np.nanmedian(P["T_meas"]), np.nanmedian(P["T_synth"])
            s["f_null_meas_hz"], s["f_null_synth_hz"] = null_freq(f, Pm, Tm, band, o["smooth_window"]), null_freq(f, Ps, Ts, band, o["smooth_window"])
            s["f_null_error_pct"] = pct_err(s["f_null_meas_hz"], s["f_null_synth_hz"])
            s["f_null_expected_hz"] = dict(meas=1.43 / Tm if Tm > 0 else np.nan, synth=1.43 / Ts if Ts > 0 else np.nan)
        L = min(len(w) for w in P["meas"])
        nseg = max(16, min(256, L // 4))
        fc, C, K = pooled_coherence(P["coh"], fs, nseg, o["nfft"])
        s["coherence_segments_K"], s["coherence_segment_len"] = K, nseg
        if K >= o["min_segments"]:
            k = (fc >= band[0]) & (fc <= band[1])
            s["coherence_mean"] = float(np.mean(C[k]))
            s["coherence_bw_above_%.1f_hz" % o["coherence_threshold"]] = coherence_bandwidth(fc, C, band, o["coherence_threshold"])
        else:
            s["coherence_mean"] = None
            s["coherence_note"] = f"only {K} segments (< {o['min_segments']}): coherence would be biased towards 1; pool more channels/recordings"
        if P["tdoa"]:
            e = [v for t_s, t_r in P["tdoa"] for v in [t_s[k] - t_r[k] for k in t_s if k in t_r]]
            ei = [v for t_s, t_r in P["tdoa"] for v in [t_s[k] - t_r[k] for k in t_s if k in t_r and k.startswith("0,")]]
            e, ei = np.array(e), np.array(ei)
            s["tdoa"] = dict(method=o["tdoa_method"],
                             compared="measured vs model labels" if (ref_tdoa or any((tt or {}).get("ref_tdoa") for tt in (times or []))) else "synthetic vs measured",
                             mae_us=float(np.mean(np.abs(e)) * 1e6), rmse_us=float(np.sqrt(np.mean(e ** 2)) * 1e6),
                             bias_us=float(np.mean(e) * 1e6), n_pairs=int(len(e)),
                             independent_mae_us=float(np.mean(np.abs(ei)) * 1e6) if len(ei) else None)
        summary[kind] = s
    if delay_err:
        d = np.array(delay_err)
        summary["sw_to_mb_delay"] = dict(mae_us=float(np.mean(np.abs(d)) * 1e6), bias_us=float(np.mean(d) * 1e6), n=int(len(d)))
    if verbose:
        print_summary(summary)
    return per_pair, summary


# ===========================================================================
# Output
# ===========================================================================
def _fmt(v):
    if v is None:
        return "-"
    if isinstance(v, float):
        return f"{v:.4g}" if np.isfinite(v) else "nan"
    if isinstance(v, (list, tuple)):
        return "-".join(_fmt(a) for a in v)
    return str(v)


def print_summary(summary):
    for kind in ("nwave", "blast"):
        s = summary.get(kind)
        if s is None:
            continue
        title = "N-wave (shockwave)" if kind == "nwave" else "Muzzle blast"
        print(f"\n=== {title} ===")
        if isinstance(s, str):
            print("  " + s); continue
        print("  4.3.1  windows {n_windows} | Pearson r {r} | SNR_synth {snr} dB | peak error {pe} % | "
              "T meas/synth {tm}/{ts} ms (error {te} %) | polarity mismatches {pm}".format(
                  n_windows=s["n_windows"], r=_fmt(s["pearson_r"]), snr=_fmt(s["snr_synth_db"]), pe=_fmt(s["peak_error_pct"]),
                  tm=_fmt(s["T_meas_ms"]), ts=_fmt(s["T_synth_ms"]), te=_fmt(s["T_error_pct"]), pm=s["polarity_mismatches"]))
        line = f"  4.3.2  band {_fmt(s['band_hz'])} Hz | LSD {_fmt(s['lsd_db'])} dB"
        if kind == "blast":
            line += f" | f_pk meas/synth {_fmt(s['f_pk_meas_hz'])}/{_fmt(s['f_pk_synth_hz'])} Hz (error {_fmt(s['f_pk_error_pct'])} %)"
        else:
            line += f" | f_null meas/synth {_fmt(s['f_null_meas_hz'])}/{_fmt(s['f_null_synth_hz'])} Hz (error {_fmt(s['f_null_error_pct'])} %)"
        print(line)
        if s.get("coherence_mean") is None:
            print(f"         coherence: {s.get('coherence_note')}")
        else:
            bw = [v for k, v in s.items() if k.startswith("coherence_bw")][0]
            print(f"         coherence: K {s['coherence_segments_K']} segments, mean {_fmt(s['coherence_mean'])}, bandwidth > 0.8: {_fmt(bw)} Hz")
        if "tdoa" in s:
            t = s["tdoa"]
            print(f"  4.3.3  TDOA ({t['method']}, {t['compared']}) | MAE {_fmt(t['mae_us'])} us | RMSE {_fmt(t['rmse_us'])} us | "
                  f"Bias {_fmt(t['bias_us'])} us | pairs {t['n_pairs']} | independent (ref ch 1) MAE {_fmt(t['independent_mae_us'])} us")
    if "sw_to_mb_delay" in summary:
        d = summary["sw_to_mb_delay"]
        print(f"\n=== 4.3.4 shock-to-blast delay === MAE {_fmt(d['mae_us'])} us | Bias {_fmt(d['bias_us'])} us | channels {d['n']}")


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def plot_pairs(per_pair, pairs, summary, out_dir, opts):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(out_dir, exist_ok=True)
    o = dict(DEFAULTS, **(opts or {}))
    for rec, (xm, fm, xs, fsn) in zip(per_pair, pairs):
        fs = rec["fs"]; xm0, xs0 = to_rate(xm, fm, fs), to_rate(xs, fsn, fs)
        fig, axs = plt.subplots(2, 2, figsize=(14, 7))
        for col, kind in enumerate(("nwave", "blast")):
            tm, ts = (rec["meas_sw_s"], rec["synth_sw_s"]) if kind == "nwave" else (rec["meas_mb_s"], rec["synth_mb_s"])
            if tm is None or ts is None:
                axs[0, col].set_title(f"{kind}: not found"); continue
            sk = summary.get(kind) if isinstance(summary.get(kind), dict) else {}
            lo = (sk.get("band_hz") or [50.0])[0]
            if o["highpass"]:
                sos = butter(2, lo, btype="high", fs=fs, output="sos")
                xm, xs = sosfiltfilt(sos, xm0, axis=0), sosfiltfilt(sos, xs0, axis=0)
            else:
                xm, xs = xm0, xs0
            L, pre = o[f"{kind}_window_ms"] / 1e3, o[f"{kind}_pre_ms"] / 1e3
            wm, ws = cut(xm[:, 0], fs, channel_onset(xm[:, 0], fs, tm, kind) - pre, L), cut(xs[:, 0], fs, channel_onset(xs[:, 0], fs, ts, kind) - pre, L)
            if wm is None or ws is None:
                continue
            lag, _ = best_lag(wm, ws, fs, o["max_lag_ms"] / 1e3); am, as_ = align(wm, ws, fs, lag)
            t = np.arange(len(am)) / fs * 1e3
            axs[0, col].plot(t, am, label="measured"); axs[0, col].plot(t, as_, label="synthetic", alpha=.8)
            axs[0, col].set_title(f"{kind}, channel 1 (aligned, lag {lag*1e3:+.3f} ms)"); axs[0, col].set_xlabel("ms"); axs[0, col].legend(); axs[0, col].grid(alpha=.3)
            f, Pm = event_psd(am, fs, o["nfft"]); _, Ps = event_psd(as_, fs, o["nfft"])
            axs[1, col].semilogx(f[1:], 10 * np.log10(Pm[1:] + 1e-30), label="measured")
            axs[1, col].semilogx(f[1:], 10 * np.log10(Ps[1:] + 1e-30), label="synthetic")
            b = (summary.get(kind) or {}).get("band_hz") if isinstance(summary.get(kind), dict) else None
            if b:
                axs[1, col].axvspan(b[0], b[1], color="0.9", zorder=0)
            axs[1, col].set_xlabel("Hz"); axs[1, col].set_ylabel("dB"); axs[1, col].legend(); axs[1, col].grid(alpha=.3, which="both")
        fig.tight_layout(); p = os.path.join(out_dir, f"similarity_pair{rec['pair']:02d}.png"); fig.savefig(p, dpi=110); plt.close(fig)
        print(f"Saved {p}")


def main():
    ap = argparse.ArgumentParser(description="Measured vs synthesized gunshot similarity (Paper-Generator Section 4.3).")
    ap.add_argument("measured", nargs="?"); ap.add_argument("synthetic", nargs="?")
    ap.add_argument("--pair", nargs=2, action="append", metavar=("MEASURED", "SYNTHETIC"), help="repeat to pool recordings")
    ap.add_argument("--event", choices=["nwave", "blast", "both"], default="both")
    ap.add_argument("--channels", type=int, nargs="*", help="1-based channels (default: all)")
    for k in ("meas_sw", "meas_mb", "synth_sw", "synth_mb"):
        ap.add_argument(f"--{k.replace('_', '-')}-ms", type=float, help="event time by hand (first pair only)")
    ap.add_argument("--nwave-window-ms", type=float, default=DEFAULTS["nwave_window_ms"])
    ap.add_argument("--blast-window-ms", type=float, default=DEFAULTS["blast_window_ms"])
    ap.add_argument("--nwave-band", type=float, nargs=2)
    ap.add_argument("--blast-band", type=float, nargs=2, default=list(DEFAULTS["blast_band"]))
    ap.add_argument("--max-lag-ms", type=float, default=DEFAULTS["max_lag_ms"])
    ap.add_argument("--min-segments", type=int, default=DEFAULTS["min_segments"])
    ap.add_argument("--tdoa", choices=["xcorr", "gcc"], default=DEFAULTS["tdoa_method"])
    ap.add_argument("--tdoa-max-lag-ms", type=float, default=DEFAULTS["tdoa_max_lag_ms"])
    ap.add_argument("--ref-tdoa", help="JSON with reference TDOAs per event (seconds, keys 'i,j' 0-based)")
    ap.add_argument("--synth-json", action="append", help="generator JSON per pair (repeat in pair order): "
                    "synthetic event times + model TDOA labels")
    ap.add_argument("--meas-scale", type=float, default=1.0, help="multiply measured samples (e.g. 1520 = Pa for PILAR-V)")
    ap.add_argument("--synth-scale", type=float, default=1.0, help="multiply synthetic samples")
    ap.add_argument("--no-highpass", action="store_true", help="do not high-pass both signals at the band's lower edge")
    ap.add_argument("--json", help="write per-pair results and the summary to this file")
    ap.add_argument("--plot", metavar="DIR", help="save overlay/spectrum plots (channel 1) per pair")
    a = ap.parse_args()

    files = list(a.pair or [])
    if a.measured and a.synthetic:
        files.insert(0, (a.measured, a.synthetic))
    if not files:
        ap.error("give measured and synthetic files (or --pair ...)")
    ch = [c - 1 for c in a.channels] if a.channels else None
    pairs = []
    for m, s in files:
        fm, xm = load_wav(m, ch); fsn, xs = load_wav(s, ch)
        xm, xs = xm * a.meas_scale, xs * a.synth_scale
        print(f"{m}: {fm} Hz, {xm.shape[1]} ch, {len(xm)/fm:.3f} s   |   {s}: {fsn} Hz, {xs.shape[1]} ch, {len(xs)/fsn:.3f} s")
        pairs.append((xm, fm, xs, fsn))
    times = [dict(meas_sw=a.meas_sw_ms and a.meas_sw_ms / 1e3, meas_mb=a.meas_mb_ms and a.meas_mb_ms / 1e3,
                  synth_sw=a.synth_sw_ms and a.synth_sw_ms / 1e3, synth_mb=a.synth_mb_ms and a.synth_mb_ms / 1e3)] + [{}] * (len(pairs) - 1)
    ref = json.load(open(a.ref_tdoa)) if a.ref_tdoa else None
    for k, gj in enumerate(a.synth_json or []):
        if k < len(times):
            tg, rg = read_generator_json(gj)
            times[k] = {**{kk: vv for kk, vv in times[k].items() if vv is not None}, **tg, "ref_tdoa": rg}
    opts = dict(nwave_window_ms=a.nwave_window_ms, blast_window_ms=a.blast_window_ms,
                nwave_band=tuple(a.nwave_band) if a.nwave_band else None, blast_band=tuple(a.blast_band),
                max_lag_ms=a.max_lag_ms, min_segments=a.min_segments, tdoa_method=a.tdoa, tdoa_max_lag_ms=a.tdoa_max_lag_ms,
                highpass=not a.no_highpass)
    per_pair, summary = evaluate(pairs, a.event, times, ref, opts)
    if a.json:
        with open(a.json, "w", encoding="utf-8") as f:
            json.dump(_clean(dict(files=files, options=opts, summary=summary, per_pair=per_pair)), f, indent=2)
        print(f"\nJSON: {a.json}")
    if a.plot:
        plot_pairs(per_pair, pairs, summary, a.plot, opts)


if __name__ == "__main__":
    main()
