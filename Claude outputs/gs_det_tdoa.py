#!/usr/bin/env python3
"""
gs_det_tdoa.py
================
DETECTION STAGE 2. Reads Stage 1's prepared signal + onset detections and
computes the time difference of arrival (TDOA) for all C(N,2) microphone pairs
(6 for the 4-mic arrays: PILAR-V, Tetrahedral-0,3 and UMA8) -- separately for
the shockwave and muzzle blast events. Either event may be missing (Stage 1's
classify method reports SW only / MB only): only the TDOAs of the events that
were found are computed, e.g. MB only -> bearing without range downstream.

TDOA method ([tdoa] method):
  xcorr     (default, Nov 2026) plain cross-correlation of the band-passed event
            (SW xcorr_sw_band_hz, MB xcorr_mb_band_hz), short Tukey-tapered windows,
            sinc-upsampled x interp_factor, positive peak only.
  gcc_phat  GCC-PHAT on the raw windows (the original method). It fails on real
            recordings whose background noise is the same on all channels (Cipatat:
            correlation 0.8-0.9): PHAT whitens the spectrum and locks onto that zero-delay
            noise, giving ~0 us for every pair. Fine for the generator's independent noise.
  A deconvolved signal (Stage 1 deconvolution) always uses peak timing.

Critical design point: windows are extracted around a SHARED reference
time (one reference channel's own detected onset), NOT each channel's own
independently-detected onset. Centering each channel's window on its own
onset silently cancels out the very TDOA this step exists to measure --
GCC-PHAT would then always report tau=0 regardless of true bearing. The
reference channel is chosen as whichever channel had the most reliable
detection in Stage 1.

Usage
-----
  python gs_det_tdoa.py det_result_prepared.json                # det_config.ini
  python gs_det_tdoa.py det_result_prepared.json my_det_config.ini
  python gs_det_tdoa.py det_result_prepared.json -mic_config=UMA8   # normally not needed

-mic_config defaults to the preset recorded by Stage 1 in *_prepared.json.

Several shots (Stage 1 [detection] multi_shot): the TDOAs are computed for every shot and
written to 'shots' in *_tdoa.json; the top-level fields are those of the first shot.
  python gs_det_tdoa.py det_result_prepared.json --method gcc_phat   # the original method
"""

import argparse
import json
from datetime import datetime, timezone
from itertools import combinations

import numpy as np
from scipy.io import wavfile
from scipy.signal import butter, resample, sosfiltfilt
from scipy.signal.windows import tukey

from gs_det_signal_prepare import load_det_config, add_mic_config_argument, upstream_mic_config


def gcc_phat(x1, x2, fs, interp_factor=16, max_lag_s=0.001):
    n = len(x1)
    X1 = np.fft.rfft(x1, n=n)
    X2 = np.fft.rfft(x2, n=n)
    G = X2 * np.conj(X1)
    G /= (np.abs(G) + 1e-12)
    nf = n * interp_factor
    cc = np.fft.fftshift(np.fft.irfft(G, n=nf))
    lags = np.arange(-nf // 2, nf // 2) / (fs * interp_factor)
    mask = np.abs(lags) <= max_lag_s
    idx = np.argmax(cc[mask])
    return float(lags[mask][idx])


def bandpass(x, fs, band, order=4):
    lo, hi = band
    hi = min(hi, 0.45 * fs)
    return sosfiltfilt(butter(order, [lo, hi], btype="band", fs=fs, output="sos"), x)


def xcorr_tdoa(x1, x2, fs, interp_factor=16, max_lag_s=0.001):
    """TDOA (s) of x2 relative to x1 (positive: x2 later) from the positive peak of the plain
    cross-correlation of two equal-length, already band-passed windows. Both windows are
    Tukey-tapered and sinc-upsampled by interp_factor before correlating."""
    n = min(len(x1), len(x2))
    w = tukey(n, 0.2)
    a = resample((x1[:n] - np.mean(x1[:n])) * w, n * interp_factor)
    b = resample((x2[:n] - np.mean(x2[:n])) * w, n * interp_factor)
    nf = 1 << int(np.ceil(np.log2(2 * len(a))))
    cc = np.fft.irfft(np.fft.rfft(b, nf) * np.conj(np.fft.rfft(a, nf)), nf)
    m = int(round(max_lag_s * fs * interp_factor))
    cc = np.concatenate([cc[-m:], cc[:m + 1]])          # lags -m .. +m
    return float((int(np.argmax(cc)) - m) / (fs * interp_factor))


def peak_timing_tdoa(signal_ref, signal_ch, t_axis, t_onset_ref, t_onset_ch,
                      fs, search_start_s=0.0, search_end_s=0.0005):
    """Direct impulse-peak TDOA estimator -- used instead of GCC-PHAT when
    deconvolution has been applied.

    Searches for the peak in [t_onset + search_start_s, t_onset + search_end_s]
    per channel. For the shockwave (impulse-like N-wave), use start≈0, end≈0.5ms.
    For the muzzle blast (Friedlander), the onset detection finds the rising edge,
    but the actual pressure peak arrives ~t_pos later -- use start≈0.5ms, end≈3ms."""
    def find_peak_time(sig, t_center):
        lo = np.searchsorted(t_axis, t_center + search_start_s)
        hi = np.searchsorted(t_axis, t_center + search_end_s)
        hi = min(hi, len(sig))
        if lo >= hi:
            return t_center + search_start_s
        return float(t_axis[lo + np.argmax(np.abs(sig[lo:hi]))])

    t_peak_ref = find_peak_time(signal_ref, t_onset_ref)
    t_peak_ch  = find_peak_time(signal_ch,  t_onset_ch)
    return t_peak_ch - t_peak_ref


def extract_window(signal, t_axis, t_center, t_pre, t_post):
    mask = (t_axis >= t_center - t_pre) & (t_axis <= t_center + t_post)
    return signal[mask]


def pick_reference_channel(per_channel):
    """Prefer a channel with both SW and reliable MB detected; fall back to
    SW-only, then to MB-only (blast-only recordings)."""
    candidates = [c for c in per_channel if c.get("sw_detected") and c.get("mb_detected") and c.get("mb_reliable")]
    if not candidates:
        candidates = [c for c in per_channel if c.get("sw_detected")]
    if not candidates:
        candidates = [c for c in per_channel if c.get("mb_detected") and c.get("mb_reliable")]
    if not candidates:
        return None
    return min(candidates, key=lambda c: c["channel"])


def _median_key(per_channel, key):
    """Median of |per_channel[key]| over the channels with a SW and a value, or None."""
    v = [abs(c[key]) for c in per_channel if c.get("sw_detected") and c.get(key)]
    return float(np.median(v)) if v else None


def _median_tn(per_channel):
    """Median N-wave duration over the channels that measured one (s), or None."""
    v = [c["T_N"] for c in per_channel if c.get("sw_detected") and c.get("T_N")]
    return float(np.median(v)) if v else None


def tdoa_one_shot(per_channel, ctx, verbose=True):
    """TDOAs of one shot from its Stage-1 per-channel detections. Returns a dict (reference
    channel, tdoa_sw_s, tdoa_mb_s, dt_sw_to_mb_s, amplitudes, T_N, consistency) or None."""
    signals, t_axis, fs, tw = ctx["signals"], ctx["t_axis"], ctx["fs"], ctx["tw"]
    method, use_peak_timing, max_lag = ctx["method"], ctx["use_peak_timing"], ctx["max_lag"]
    n_channels = len(signals)
    ref = pick_reference_channel(per_channel)
    if ref is None:
        return None
    ref_ch = ref["channel"]
    sw_used = bool(ref.get("sw_detected"))
    mb_used = bool(ref.get("mb_detected") and ref.get("mb_reliable"))
    if verbose:
        t_sw_txt = f"{ref['t_sw']:.5f}s" if sw_used else "none"
        print(f"[Det Stage 2] Reference channel: {ref_ch} (t_sw={t_sw_txt}, mb_reliable={ref.get('mb_reliable')})")

    def band_signals(key, band):
        if key not in ctx["cache"]:
            ctx["cache"][key] = [bandpass(s, fs, band) for s in signals]
        return ctx["cache"][key]

    pairs = list(combinations(range(n_channels), 2))

    def one(kind, t_ref, post=None):
        out = {}
        xpost = tw[f"xcorr_{kind}_window_post_s"] if post is None else post
        if method == "xcorr" and not use_peak_timing:
            sig = band_signals(kind, tw[f"xcorr_{kind}_band_hz"])
        for (a, b) in pairs:
            if use_peak_timing:
                ta = per_channel[a].get(f"t_{kind}") if per_channel[a].get(f"{kind}_detected") else t_ref
                tb = per_channel[b].get(f"t_{kind}") if per_channel[b].get(f"{kind}_detected") else t_ref
                lo, hi = (-0.0001, 0.0005) if kind == "sw" else (0.0005, 0.003)
                tau = peak_timing_tdoa(signals[a], signals[b], t_axis, ta, tb, fs,
                                       search_start_s=lo, search_end_s=hi)
            elif method == "xcorr":
                wa = extract_window(sig[a], t_axis, t_ref, tw[f"xcorr_{kind}_window_pre_s"], xpost)
                wb = extract_window(sig[b], t_axis, t_ref, tw[f"xcorr_{kind}_window_pre_s"], xpost)
                tau = xcorr_tdoa(wa, wb, fs, interp_factor=tw["interp_factor"], max_lag_s=max_lag)
            else:
                wa = extract_window(signals[a], t_axis, t_ref, tw[f"{kind}_window_pre_s"], tw[f"{kind}_window_post_s"])
                wb = extract_window(signals[b], t_axis, t_ref, tw[f"{kind}_window_pre_s"], tw[f"{kind}_window_post_s"])
                tau = gcc_phat(wa, wb, fs, interp_factor=tw["interp_factor"], max_lag_s=max_lag)
            out[f"{a},{b}"] = tau
        return out

    tdoa_sw = one("sw", ref["t_sw"]) if sw_used else {}
    def closure(td):
        return {f"{a},{b}": float(td.get(f"{a},{b}", 0) - (td.get(f"0,{b}", 0) - td.get(f"0,{a}", 0)))
                for (a, b) in (pairs if td else []) if a != 0}

    def worst(td):
        return max((abs(v) for v in closure(td).values()), default=0.0)

    tdoa_mb = one("mb", ref["t_mb"]) if mb_used else {}
    mb_post_used = tw.get("xcorr_mb_window_post_s")
    # The MB window is long (it must hold the blast's negative phase and ground reflection), so
    # another event can fall into it: in a burst, the SW or an echo of a later round. Then the
    # pair TDOAs stop adding up. Retry with shorter windows and keep the most consistent set.
    tol = tw.get("mb_closure_tol_s", 0.0)
    if tdoa_mb and method == "xcorr" and not use_peak_timing and tol > 0 and worst(tdoa_mb) > tol:
        full = tw["xcorr_mb_window_post_s"]
        best = (worst(tdoa_mb), full, tdoa_mb)
        for frac in (0.5, 0.25):
            td = one("mb", ref["t_mb"], post=full * frac)
            if worst(td) < best[0]:
                best = (worst(td), full * frac, td)
        if verbose:
            print(f"[Det Stage 2] MB pair TDOAs inconsistent ({worst(tdoa_mb)*1e6:.0f} us); "
                  f"MB window {best[1]*1e3:.1f} ms used ({best[0]*1e6:.0f} us)")
        tdoa_mb, mb_post_used = best[2], best[1]
    if verbose:
        if not mb_used:
            print("[Det Stage 2] WARNING: reference channel's MB detection wasn't reliable -- "
                  "skipping MB TDOA (would be meaningless windowed around a bad onset guess).")
        if not sw_used:
            print("[Det Stage 2] No shockwave detected -- MB TDOA only (bearing, no range).")
        for name, td in (("SW", tdoa_sw), ("MB", tdoa_mb)):
            if td:
                print(f"{name} TDOA (us):")
                for k, v in td.items():
                    print(f"  ({k}): {v*1e6:+.3f}")
    # Consistency check: tau(1,2) should ~= tau(0,2)-tau(0,1), etc. Reported, not corrected.
    consistency_sw, consistency_mb = closure(tdoa_sw), closure(tdoa_mb)
    return dict(reference_channel=ref_ch, sw_available=bool(tdoa_sw), mb_available=bool(tdoa_mb),
                dt_sw_to_mb_s=(ref["t_mb"] - ref["t_sw"]) if (mb_used and sw_used) else None,
                dp_sw_ref_channel=ref.get("dP_sw"), t_n_sw_ref_channel=ref.get("T_N"),
                t_n_sw_median=_median_tn(per_channel),
                sw_clipped=any(bool(c.get("sw_clipped")) for c in per_channel if c.get("sw_detected")),
                dp_sw_unclipped_est=_median_key(per_channel, "dP_sw_unclipped_est"),
                dp_mb_ref_channel=ref.get("dP_mb"), tdoa_sw_s=tdoa_sw, tdoa_mb_s=tdoa_mb,
                sw_consistency_residual_s=consistency_sw, mb_consistency_residual_s=consistency_mb,
                mb_window_post_s=mb_post_used)


def main():
    parser = argparse.ArgumentParser(
        description="Detection Stage 2: TDOA for all C(N,2) mic pairs (SW and MB).")
    parser.add_argument("prepared_json", help="Path to Stage 1's *_prepared.json")
    parser.add_argument("config", nargs="?", default="det_config.ini",
                         help="Path to detection .ini config (default: det_config.ini)")
    parser.add_argument("-o", "--output", default=None,
                         help="Output basename (overrides output.basename in config)")
    parser.add_argument("--method", choices=["xcorr", "gcc_phat"], default=None,
                        help="TDOA method (default: [tdoa] method in the .ini, else xcorr)")
    add_mic_config_argument(parser)
    args = parser.parse_args()

    with open(args.prepared_json) as f:
        prep = json.load(f)

    cfg = load_det_config(args.config, args.mic_config or upstream_mic_config(prep))
    if args.output:
        cfg["output"]["basename"] = args.output

    fs, data = wavfile.read(prep["prepared_wav"])
    n_channels = data.shape[1]
    if n_channels != cfg["geometry"]["n_mics"]:
        raise SystemExit(f"[Det Stage 2] ERROR: {prep['prepared_wav']} has {n_channels} channels but "
                         f"mic_config {cfg['geometry']['mic_config']} has {cfg['geometry']['n_mics']} "
                         f"mics -- re-run Stage 1 with the matching -mic_config.")
    print(f"[Det Stage 2] mic_config = {cfg['geometry']['mic_config']} "
          f"({n_channels} mics, {n_channels * (n_channels - 1) // 2} pairs)")
    signals = [data[:, i].astype(np.float64) for i in range(n_channels)]
    t_axis = np.arange(len(signals[0])) / fs

    tw = cfg["tdoa"]
    # the lag search must cover the array's largest possible delay (PILAR-V: 0.356 m -> ~1.04 ms,
    # more than the 1 ms default max_lag_s)
    max_lag = max(tw["max_lag_s"], 1.1 * cfg["geometry"]["aperture_m"] / cfg["physics"]["speed_of_sound"])
    if max_lag > tw["max_lag_s"]:
        print(f"[Det Stage 2] max_lag_s raised to {max_lag * 1e3:.2f} ms to cover the array aperture")
    use_peak_timing = prep.get("deconvolution_applied", False)
    method = (args.method or tw.get("method", "xcorr")).replace("-", "_")
    if method not in ("xcorr", "gcc_phat"):
        raise SystemExit(f"[Det Stage 2] unknown [tdoa] method '{method}' (xcorr or gcc_phat)")
    method_label = ("peak-timing (deconvolved signal)" if use_peak_timing else
                    "cross-correlation of the band-passed event" if method == "xcorr" else "GCC-PHAT")
    print(f"[Det Stage 2] TDOA method: {method_label}")
    ctx = dict(signals=signals, t_axis=t_axis, fs=fs, tw=tw, method=method,
               use_peak_timing=use_peak_timing, max_lag=max_lag, cache={})

    main_shot = tdoa_one_shot(prep["per_channel"], ctx, verbose=True)
    if main_shot is None:
        raise SystemExit("[Det Stage 2] No shot was detected in Stage 1 "
                         f"(event_class = {prep.get('event_class')}"
                         + (f"; {len(prep.get('blast_only_events') or [])} blast(s) without a shockwave "
                            "not reported: supersonic_only" if prep.get("blast_only_events") else "")
                         + ") -- nothing to compute TDOA on.")
    ref_ch = main_shot["reference_channel"]
    tdoa_sw, tdoa_mb = main_shot["tdoa_sw_s"], main_shot["tdoa_mb_s"]

    # several shots (Stage 1 [detection] multi_shot): the same for every shot
    shots_out = []
    for k, sh in enumerate(prep.get("shots") or []):
        r = tdoa_one_shot(sh["per_channel"], ctx, verbose=False)
        if r is None:
            continue
        r.update(index=sh.get("index", k + 1), event_class=sh.get("event_class"),
                 t_sw_s=sh.get("t_sw_s"), t_mb_s=sh.get("t_mb_s"))
        shots_out.append(r)
    if len(shots_out) > 1:
        print(f"[Det Stage 2] {len(shots_out)} shots: TDOA computed for each (see 'shots' in the .json)")

    out_json = f"{cfg['output']['basename']}_tdoa.json"
    meta = dict(
        generated_at_utc=datetime.now(timezone.utc).isoformat(),
        stage="tdoa",
        source_prepared_json=args.prepared_json,
        reference_channel=ref_ch,
        tdoa_method=("peak_timing" if use_peak_timing else method),
        sample_rate_hz=fs,
        n_channels=n_channels,
        sw_available=main_shot["sw_available"],
        mb_available=main_shot["mb_available"],
        event_class=prep.get("event_class"),
        dt_sw_to_mb_s=main_shot["dt_sw_to_mb_s"],
        dp_sw_ref_channel=main_shot["dp_sw_ref_channel"],
        # Stage 1 peak-normalizes the signal; this factor turns dp_* back into full-scale units
        normalization_peak=prep.get("normalization_peak"),
        t_n_sw_ref_channel=main_shot["t_n_sw_ref_channel"],
        t_n_sw_median=main_shot["t_n_sw_median"],
        sw_clipped=main_shot["sw_clipped"],
        dp_sw_unclipped_est=main_shot["dp_sw_unclipped_est"],
        dp_mb_ref_channel=main_shot["dp_mb_ref_channel"],
        tdoa_sw_s=tdoa_sw,
        tdoa_mb_s=tdoa_mb,
        sw_consistency_residual_s=main_shot["sw_consistency_residual_s"],
        mb_consistency_residual_s=main_shot["mb_consistency_residual_s"],
        n_shots=max(1, len(shots_out)),
        shots=shots_out,
        config_used=cfg,
    )
    with open(out_json, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"[Det Stage 2] Wrote {out_json}")


if __name__ == "__main__":
    main()
