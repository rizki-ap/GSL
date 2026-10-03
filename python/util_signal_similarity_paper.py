#!/usr/bin/env python3
"""
Paper-specific signal evaluation for the Gunshot Audio Generator paper.

Implements MD Sections 4.3.1-4.3.2:
  4.3.1: Pearson r, SNR_synth, peak-amplitude error, duration error
  4.3.2: pooled magnitude-squared coherence, >0.8 bandwidth,
          Tukey-window event PSD, averaged PSD, LSD, spectral peak/null errors.

This file intentionally preserves signal amplitude. Do not independently
peak-normalize measured and synthesized WAV files before calling these
functions, because the paper evaluates amplitude fidelity.

The MD leaves numerical event-window, smoothing, Tukey-alpha, and coherence
segment parameters unspecified. They are therefore explicit arguments here
and must be fixed before producing final paper results.
"""

import argparse
import json
from pathlib import Path
import numpy as np
from scipy.io import wavfile
from scipy.signal import correlate, find_peaks, resample_poly, savgol_filter, welch
from scipy.signal.windows import tukey


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def load_wav(path):
    """Load mono WAV without independent peak normalization."""
    fs, x = wavfile.read(path)
    if x.ndim > 1:
        x = x[:, 0]
    if np.issubdtype(x.dtype, np.integer):
        info = np.iinfo(x.dtype)
        x = x.astype(np.float64) / max(abs(info.min), info.max)
    else:
        x = x.astype(np.float64)
    return int(fs), x


def match_sample_rates(x1, fs1, x2, fs2, target_fs=None):
    """Resample both signals to a common rate."""
    if target_fs is None:
        target_fs = fs1 if fs1 == fs2 else min(fs1, fs2)
    target_fs = int(target_fs)

    def rs(x, fs):
        if int(fs) == target_fs:
            return np.asarray(x, dtype=np.float64)
        from math import gcd
        g = gcd(int(fs), target_fs)
        return resample_poly(x, target_fs // g, int(fs) // g)

    return rs(x1, fs1), rs(x2, fs2), target_fs


def extract_event(x, fs, start_s, duration_s):
    """Extract a fixed event window."""
    i = int(round(start_s * fs))
    n = int(round(duration_s * fs))
    if i < 0 or n <= 0 or i + n > len(x):
        raise ValueError("Event window is outside the WAV duration.")
    return np.asarray(x[i:i+n], dtype=np.float64).copy()


# ---------------------------------------------------------------------------
# Alignment
# ---------------------------------------------------------------------------

def best_fit_lag(x, y, fs, max_lag_s=0.005):
    """Return lag in seconds; positive means y arrives after x."""
    n = min(len(x), len(y))
    x = np.asarray(x[:n]) - np.mean(x[:n])
    y = np.asarray(y[:n]) - np.mean(y[:n])
    c = correlate(x, y, mode="full")
    lags = np.arange(-len(y) + 1, len(x))
    keep = np.abs(lags) <= int(round(max_lag_s * fs))
    lags, c = lags[keep], c[keep]
    return float(-lags[np.argmax(np.abs(c))] / fs)


def align_by_lag(x, y, fs, lag_s):
    """Crop overlapping portions after applying the supplied lag."""
    k = int(round(lag_s * fs))
    if k > 0:
        a, b = x[:max(0, len(x)-k)], y[k:]
    elif k < 0:
        k = -k
        a, b = x[k:], y[:max(0, len(y)-k)]
    else:
        a, b = x, y
    n = min(len(a), len(b))
    if n < 2:
        raise ValueError("No usable overlap after alignment.")
    return np.asarray(a[:n]), np.asarray(b[:n])


# ---------------------------------------------------------------------------
# Section 4.3.1 -- waveform morphology
# ---------------------------------------------------------------------------

def pearson_r(x, y):
    """Pearson correlation coefficient from MD Section 4.3.1."""
    n = min(len(x), len(y))
    x, y = np.asarray(x[:n]), np.asarray(y[:n])
    x0, y0 = x-x.mean(), y-y.mean()
    den = np.sqrt(np.sum(x0*x0) * np.sum(y0*y0))
    return np.nan if den <= 1e-30 else float(np.sum(x0*y0)/den)


def snr_synth_db(measured, synthesized):
    """MD SNR_synth = 10 log10(sum(x^2) / sum((x-xhat)^2))."""
    n = min(len(measured), len(synthesized))
    x, y = measured[:n], synthesized[:n]
    num = np.sum(x*x)
    den = np.sum((x-y)**2)
    if num <= 1e-30:
        return np.nan
    if den <= 1e-30:
        return np.inf
    return float(10*np.log10(num/den))


def peak_overpressure(x):
    """Delta-p estimate from the supplied event window."""
    return float(np.max(np.abs(x))) if len(x) else np.nan


def peak_amplitude_error_percent(measured, synthesized):
    pm, ps = peak_overpressure(measured), peak_overpressure(synthesized)
    if not np.isfinite(pm) or pm <= 1e-30:
        return np.nan
    return float(abs(ps-pm)/abs(pm)*100.0)


def duration_error_percent(T_meas_s, T_synth_s):
    if not np.isfinite(T_meas_s) or T_meas_s <= 0:
        return np.nan
    return float(abs(T_synth_s-T_meas_s)/T_meas_s*100.0)


# ---------------------------------------------------------------------------
# Section 4.3.2 -- event PSD and LSD
# ---------------------------------------------------------------------------

def event_psd(x, fs, nfft=4096, tukey_alpha=0.25):
    """
    PSD using a Tukey-windowed event and 4096-point zero padding.

    nfft is interpolation only; it does not increase physical resolution.
    """
    x = np.asarray(x, dtype=np.float64)
    w = tukey(len(x), alpha=tukey_alpha)
    xw = (x-x.mean()) * w
    nfft = max(int(nfft), len(x))
    return welch(
        xw, fs=fs, window="boxcar", nperseg=len(x), noverlap=0,
        nfft=nfft, detrend=False, scaling="density"
    )


def average_psd(events, fs, f_lo=50.0, f_hi=8000.0,
                nfft=4096, tukey_alpha=0.25):
    """Average PSD over supplied event windows."""
    spectra = [event_psd(e, fs, nfft, tukey_alpha) for e in events]
    if not spectra:
        raise ValueError("No events supplied.")
    f = spectra[0][0]
    P = np.mean(np.vstack([s[1] for s in spectra]), axis=0)
    keep = (f >= f_lo) & (f <= f_hi)
    if keep.sum() < 3:
        raise ValueError("Too few PSD bins in requested frequency band.")
    return f[keep], P[keep]


def log_spectral_distance_from_psd(f, P_meas, P_synth,
                                   f_lo=50.0, f_hi=8000.0):
    """MD LSD from the ratio of averaged PSDs."""
    f = np.asarray(f)
    Pm, Ps = np.asarray(P_meas), np.asarray(P_synth)
    keep = ((f >= f_lo) & (f <= f_hi) &
            np.isfinite(Pm) & np.isfinite(Ps) & (Pm > 0) & (Ps > 0))
    if keep.sum() < 3:
        return np.nan
    ff = f[keep]
    d = 10*np.log10(Ps[keep]/Pm[keep])
    return float(np.sqrt(np.trapezoid(d*d, ff)/(ff[-1]-ff[0])))


def _smooth_psd(P, window=11, polyorder=2):
    """Log-domain Savitzky-Golay smoothing for spectral feature picking."""
    P = np.asarray(P)
    if len(P) < 5:
        return P
    w = int(window)
    if w % 2 == 0:
        w += 1
    w = min(w, len(P) if len(P) % 2 else len(P)-1)
    if w < 5:
        return P
    p = min(polyorder, w-2)
    return 10**savgol_filter(np.log10(np.maximum(P, 1e-30)), w, p)


def spectral_peak_frequency(f, P, f_lo=50.0, f_hi=8000.0,
                           smooth_window=11):
    keep = (f >= f_lo) & (f <= f_hi) & np.isfinite(P) & (P > 0)
    if keep.sum() < 3:
        return np.nan
    ff, PP = f[keep], _smooth_psd(P[keep], smooth_window)
    return float(ff[np.argmax(PP)])


def spectral_first_nontrivial_null(f, P, expected_hz=None,
                                   f_lo=50.0, f_hi=8000.0,
                                   smooth_window=11):
    """
    Pick the N-wave null from the smoothed spectrum.

    If expected_hz is given, search 0.5--2.0 times that value. This follows
    the MD's stated f_null ~= 1.43/T_N while making the actual measured
    minimum explicit and reproducible.
    """
    keep = (f >= f_lo) & (f <= f_hi) & np.isfinite(P) & (P > 0)
    if keep.sum() < 5:
        return np.nan
    ff, PP = f[keep], _smooth_psd(P[keep], smooth_window)

    if expected_hz is not None and np.isfinite(expected_hz):
        region = (ff >= max(f_lo, .5*expected_hz)) & \
                 (ff <= min(f_hi, 2.0*expected_hz))
        if region.sum() >= 3:
            return float(ff[region][np.argmin(PP[region])])

    peaks, _ = find_peaks(PP)
    if len(peaks) == 0:
        return np.nan
    dominant = peaks[np.argmax(PP[peaks])]
    minima, _ = find_peaks(-PP)
    minima = minima[minima > dominant]
    return float(ff[minima[0]]) if len(minima) else np.nan


def frequency_error_percent(f_meas, f_synth):
    if not np.isfinite(f_meas) or f_meas <= 0 or not np.isfinite(f_synth):
        return np.nan
    return float(abs(f_synth-f_meas)/f_meas*100.0)


# ---------------------------------------------------------------------------
# Ensemble magnitude-squared coherence
# ---------------------------------------------------------------------------

def _coherence_segments(x, y, fs, nperseg=256, noverlap=None, nfft=4096):
    """Return per-segment Sxy/Sxx/Syy for pooled ensemble coherence."""
    n = min(len(x), len(y))
    x, y = np.asarray(x[:n]), np.asarray(y[:n])
    if n < nperseg:
        nperseg, noverlap = n, 0
    elif noverlap is None:
        noverlap = nperseg//2
    step = nperseg-noverlap
    w = tukey(nperseg, alpha=.25)
    f = np.fft.rfftfreq(nfft, 1/fs)
    out = []
    for i in range(0, n-nperseg+1, step):
        a = (x[i:i+nperseg]-np.mean(x[i:i+nperseg]))*w
        b = (y[i:i+nperseg]-np.mean(y[i:i+nperseg]))*w
        X, Y = np.fft.rfft(a, nfft), np.fft.rfft(b, nfft)
        out.append((X*np.conj(Y), abs(X)**2, abs(Y)**2))
    return f, out


def ensemble_magnitude_squared_coherence(event_pairs, fs, nperseg=256,
                                         noverlap=None, nfft=4096):
    """
    MD coherence:
      |sum_k Sxy_k|^2 / [(sum_k Sxx_k)(sum_k Syy_k)].

    event_pairs may contain events from multiple sensors and recordings.
    K therefore represents all pooled segments supplied by the caller.
    """
    sx = sy = sxy = None
    f = None
    K = 0
    for x, y in event_pairs:
        f0, segs = _coherence_segments(
            x, y, fs, nperseg, noverlap, nfft
        )
        for a, b, c in segs:
            if sxy is None:
                f = f0
                sxy = np.zeros_like(a, dtype=complex)
                sx = np.zeros_like(b, dtype=float)
                sy = np.zeros_like(c, dtype=float)
            sxy += a
            sx += b
            sy += c
            K += 1
    if K == 0:
        raise ValueError("No coherence segments.")
    C = np.zeros_like(sx)
    good = sx*sy > 1e-30
    C[good] = abs(sxy[good])**2/(sx[good]*sy[good])
    return f, np.clip(C, 0, 1), K


def coherence_bandwidth(f, C, threshold=.8, f_lo=50.0, f_hi=8000.0):
    """Total frequency width for which Cbar > threshold."""
    keep = (f >= f_lo) & (f <= f_hi)
    ff, cc = f[keep], C[keep]
    if len(ff) < 2:
        return np.nan
    return float(np.sum(np.diff(ff)[(cc[:-1] > threshold) &
                                     (cc[1:] > threshold)]))


# ---------------------------------------------------------------------------
# Main paper evaluation
# ---------------------------------------------------------------------------

def paper_evaluate_pair(measured, synthesized, fs,
                        duration_meas_s, duration_synth_s,
                        event_type="nwave", max_lag_s=.005,
                        f_lo=50.0, f_hi=8000.0, nfft=4096,
                        tukey_alpha=.25, coherence_nperseg=256,
                        coherence_threshold=.8, smooth_window=11):
    """
    Evaluate one event pair.

    Inputs must already be isolated to the event window. Best-fit lag is
    removed before morphology/spectral comparison, as required by the MD.
    """
    lag = best_fit_lag(measured, synthesized, fs, max_lag_s)
    xm, xs = align_by_lag(measured, synthesized, fs, lag)

    out = {
        "best_fit_lag_ms": lag*1000,
        "pearson_r": pearson_r(xm, xs),
        "snr_synth_db": snr_synth_db(xm, xs),
        "peak_measured": peak_overpressure(xm),
        "peak_synthesized": peak_overpressure(xs),
        "peak_amplitude_error_percent":
            peak_amplitude_error_percent(xm, xs),
        "duration_measured_s": duration_meas_s,
        "duration_synthesized_s": duration_synth_s,
        "duration_error_percent":
            duration_error_percent(duration_meas_s, duration_synth_s),
    }

    fm, Pm = average_psd([xm], fs, f_lo, f_hi, nfft, tukey_alpha)
    fs_, Ps = average_psd([xs], fs, f_lo, f_hi, nfft, tukey_alpha)
    Ps = np.interp(fm, fs_, Ps)

    out["lsd_db"] = log_spectral_distance_from_psd(
        fm, Pm, Ps, f_lo, f_hi
    )

    fpm = spectral_peak_frequency(fm, Pm, f_lo, f_hi, smooth_window)
    fps = spectral_peak_frequency(fm, Ps, f_lo, f_hi, smooth_window)
    out["f_peak_measured_hz"] = fpm
    out["f_peak_synthesized_hz"] = fps
    out["f_peak_error_percent"] = frequency_error_percent(fpm, fps)

    if event_type.lower() == "nwave":
        en_m = 1.43/duration_meas_s if duration_meas_s > 0 else np.nan
        en_s = 1.43/duration_synth_s if duration_synth_s > 0 else np.nan
        fnm = spectral_first_nontrivial_null(
            fm, Pm, en_m, f_lo, f_hi, smooth_window
        )
        fns = spectral_first_nontrivial_null(
            fm, Ps, en_s, f_lo, f_hi, smooth_window
        )
        out["f_null_measured_hz"] = fnm
        out["f_null_synthesized_hz"] = fns
        out["f_null_error_percent"] = frequency_error_percent(fnm, fns)
    else:
        out["f_null_measured_hz"] = np.nan
        out["f_null_synthesized_hz"] = np.nan
        out["f_null_error_percent"] = np.nan

    fc, C, K = ensemble_magnitude_squared_coherence(
        [(xm, xs)], fs, coherence_nperseg, None, nfft
    )
    out["coherence_segments_K"] = K
    keep = (fc >= f_lo) & (fc <= f_hi)
    out["coherence_mean"] = float(np.mean(C[keep])) if keep.any() else np.nan
    out["coherence_bandwidth_hz"] = coherence_bandwidth(
        fc, C, coherence_threshold, f_lo, f_hi
    )
    return out


def paper_evaluate_many(measured_events, synthesized_events, fs,
                        durations_meas_s, durations_synth_s,
                        event_type="nwave", **kwargs):
    """
    Preferred API for publication results.

    PSDs are averaged across all supplied events. Coherence pools segments
    across all supplied event pairs, so event_pairs can represent sensors
    and recordings.
    """
    if not (len(measured_events) == len(synthesized_events) ==
            len(durations_meas_s) == len(durations_synth_s)):
        raise ValueError("All event lists must have equal length.")

    aligned = []
    per_event = []
    for xm, xs, tm, ts in zip(
        measured_events, synthesized_events,
        durations_meas_s, durations_synth_s
    ):
        lag = best_fit_lag(
            xm, xs, fs, kwargs.get("max_lag_s", .005)
        )
        am, ass = align_by_lag(xm, xs, fs, lag)
        aligned.append((am, ass))
        per_event.append(paper_evaluate_pair(
            am, ass, fs, tm, ts, event_type=event_type, **kwargs
        ))
        per_event[-1]["best_fit_lag_ms"] = lag*1000

    f1, P1 = average_psd(
        [p[0] for p in aligned], fs,
        kwargs.get("f_lo", 50), kwargs.get("f_hi", 8000),
        kwargs.get("nfft", 4096), kwargs.get("tukey_alpha", .25)
    )
    f2, P2 = average_psd(
        [p[1] for p in aligned], fs,
        kwargs.get("f_lo", 50), kwargs.get("f_hi", 8000),
        kwargs.get("nfft", 4096), kwargs.get("tukey_alpha", .25)
    )
    P2 = np.interp(f1, f2, P2)

    agg = {
        "n_events": len(aligned),
        "pearson_r_mean": float(np.nanmean([r["pearson_r"] for r in per_event])),
        "snr_synth_db_mean": float(np.nanmean([r["snr_synth_db"] for r in per_event])),
        "peak_amplitude_error_percent_mean":
            float(np.nanmean([r["peak_amplitude_error_percent"] for r in per_event])),
        "duration_error_percent_mean":
            float(np.nanmean([r["duration_error_percent"] for r in per_event])),
        "lsd_db": log_spectral_distance_from_psd(
            f1, P1, P2, kwargs.get("f_lo", 50), kwargs.get("f_hi", 8000)
        ),
    }

    fpm = spectral_peak_frequency(
        f1, P1, kwargs.get("f_lo", 50), kwargs.get("f_hi", 8000),
        kwargs.get("smooth_window", 11)
    )
    fps = spectral_peak_frequency(
        f1, P2, kwargs.get("f_lo", 50), kwargs.get("f_hi", 8000),
        kwargs.get("smooth_window", 11)
    )
    agg["f_peak_measured_hz"] = fpm
    agg["f_peak_synthesized_hz"] = fps
    agg["f_peak_error_percent"] = frequency_error_percent(fpm, fps)

    if event_type.lower() == "nwave":
        tm = np.mean(durations_meas_s)
        ts = np.mean(durations_synth_s)
        fnm = spectral_first_nontrivial_null(
            f1, P1, 1.43/tm, kwargs.get("f_lo", 50),
            kwargs.get("f_hi", 8000), kwargs.get("smooth_window", 11)
        )
        fns = spectral_first_nontrivial_null(
            f1, P2, 1.43/ts, kwargs.get("f_lo", 50),
            kwargs.get("f_hi", 8000), kwargs.get("smooth_window", 11)
        )
        agg["f_null_measured_hz"] = fnm
        agg["f_null_synthesized_hz"] = fns
        agg["f_null_error_percent"] = frequency_error_percent(fnm, fns)
    else:
        agg["f_null_measured_hz"] = np.nan
        agg["f_null_synthesized_hz"] = np.nan
        agg["f_null_error_percent"] = np.nan

    fc, C, K = ensemble_magnitude_squared_coherence(
        aligned, fs, kwargs.get("coherence_nperseg", 256),
        None, kwargs.get("nfft", 4096)
    )
    agg["coherence_segments_K"] = K
    keep = (fc >= kwargs.get("f_lo", 50)) & (fc <= kwargs.get("f_hi", 8000))
    agg["coherence_mean"] = float(np.mean(C[keep])) if keep.any() else np.nan
    agg["coherence_bandwidth_hz"] = coherence_bandwidth(
        fc, C, kwargs.get("coherence_threshold", .8),
        kwargs.get("f_lo", 50), kwargs.get("f_hi", 8000)
    )
    return per_event, agg


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def clean_json(d):
    out = {}
    for k, v in d.items():
        if isinstance(v, np.generic):
            v = v.item()
        if isinstance(v, float) and not np.isfinite(v):
            v = None
        out[k] = v
    return out


def main():
    p = argparse.ArgumentParser(
        description="Evaluate Gunshot Generator metrics from MD §§4.3.1-4.3.2."
    )
    p.add_argument("measured")
    p.add_argument("synthesized")
    p.add_argument("--event-type", choices=["nwave", "blast"], default="nwave")
    p.add_argument("--event-start-ms", type=float, required=True)
    p.add_argument("--event-duration-ms", type=float, required=True)
    p.add_argument("--measured-start-ms", type=float)
    p.add_argument("--synth-start-ms", type=float)
    p.add_argument("--measured-duration-ms", type=float)
    p.add_argument("--synth-duration-ms", type=float)
    p.add_argument("--max-lag-ms", type=float, default=5.0)
    p.add_argument("--f-lo", type=float, default=50.0)
    p.add_argument("--f-hi", type=float, default=8000.0)
    p.add_argument("--nfft", type=int, default=4096)
    p.add_argument("--tukey-alpha", type=float, default=.25)
    p.add_argument("--coherence-nperseg", type=int, default=256)
    p.add_argument("--coherence-threshold", type=float, default=.8)
    p.add_argument("--smooth-window", type=int, default=11)
    p.add_argument("--json", type=str)
    a = p.parse_args()

    fm, xm = load_wav(a.measured)
    fs, xs = load_wav(a.synthesized)
    xm, xs, fs = match_sample_rates(xm, fm, xs, fs)

    sm = (a.measured_start_ms if a.measured_start_ms is not None
          else a.event_start_ms)/1000
    ss = (a.synth_start_ms if a.synth_start_ms is not None
          else a.event_start_ms)/1000
    dm = (a.measured_duration_ms if a.measured_duration_ms is not None
          else a.event_duration_ms)/1000
    ds = (a.synth_duration_ms if a.synth_duration_ms is not None
          else a.event_duration_ms)/1000

    em = extract_event(xm, fs, sm, dm)
    es = extract_event(xs, fs, ss, ds)

    r = paper_evaluate_pair(
        em, es, fs, dm, ds, a.event_type,
        a.max_lag_ms/1000, a.f_lo, a.f_hi, a.nfft,
        a.tukey_alpha, a.coherence_nperseg,
        a.coherence_threshold, a.smooth_window
    )

    print("\n=== Section 4.3.1: Waveform morphology ===")
    for k in ("best_fit_lag_ms", "pearson_r", "snr_synth_db",
              "peak_measured", "peak_synthesized",
              "peak_amplitude_error_percent",
              "duration_measured_s", "duration_synthesized_s",
              "duration_error_percent"):
        print(f"{k:40s}: {r[k]}")

    print("\n=== Section 4.3.2: Frequency-domain similarity ===")
    for k in ("lsd_db", "f_peak_measured_hz", "f_peak_synthesized_hz",
              "f_peak_error_percent", "f_null_measured_hz",
              "f_null_synthesized_hz", "f_null_error_percent",
              "coherence_segments_K", "coherence_mean",
              "coherence_bandwidth_hz"):
        print(f"{k:40s}: {r[k]}")

    if a.json:
        Path(a.json).write_text(
            json.dumps(clean_json(r), indent=2), encoding="utf-8"
        )
        print(f"\nJSON report: {a.json}")


if __name__ == "__main__":
    main()
