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

    ref = pick_reference_channel(prep["per_channel"])
    if ref is None:
        raise SystemExit("[Det Stage 2] No shockwave or muzzle blast was detected in Stage 1 "
                         f"(event_class = {prep.get('event_class')}) -- nothing to compute TDOA on.")
    ref_ch = ref["channel"]
    sw_used = bool(ref.get("sw_detected"))
    t_sw_txt = f"{ref['t_sw']:.5f}s" if sw_used else "none"
    print(f"[Det Stage 2] Reference channel: {ref_ch} "
          f"(t_sw={t_sw_txt}, mb_reliable={ref.get('mb_reliable')}, "
          f"events: {prep.get('event_class', 'n/a')})")

    tw = cfg["tdoa"]
    # the lag search must cover the array's largest possible delay (PILAR-V: 0.356 m -> ~1.04 ms,
    # more than the 1 ms default max_lag_s)
    max_lag = max(tw["max_lag_s"], 1.1 * cfg["geometry"]["aperture_m"] / cfg["physics"]["speed_of_sound"])
    if max_lag > tw["max_lag_s"]:
        print(f"[Det Stage 2] max_lag_s raised to {max_lag * 1e3:.2f} ms to cover the array aperture")
    pairs = list(combinations(range(n_channels), 2))
    use_peak_timing = prep.get("deconvolution_applied", False)
    method = (args.method or tw.get("method", "xcorr")).replace("-", "_")
    if method not in ("xcorr", "gcc_phat"):
        raise SystemExit(f"[Det Stage 2] unknown [tdoa] method '{method}' (xcorr or gcc_phat)")
    method_label = ("peak-timing (deconvolved signal)" if use_peak_timing else
                    "cross-correlation of the band-passed event" if method == "xcorr" else "GCC-PHAT")
    print(f"[Det Stage 2] TDOA method: {method_label}")
    if method == "xcorr" and not use_peak_timing:
        sig_sw = [bandpass(s, fs, tw["xcorr_sw_band_hz"]) for s in signals] if sw_used else None
        mb_ok = ref.get("mb_detected") and ref.get("mb_reliable")
        sig_mb = [bandpass(s, fs, tw["xcorr_mb_band_hz"]) for s in signals] if mb_ok else None

    tdoa_sw = {}
    for (a, b) in (pairs if sw_used else []):
        if use_peak_timing:
            ch_a = prep["per_channel"][a]
            ch_b = prep["per_channel"][b]
            t_sw_a = ch_a.get("t_sw") if ch_a.get("sw_detected") else ref["t_sw"]
            t_sw_b = ch_b.get("t_sw") if ch_b.get("sw_detected") else ref["t_sw"]
            tau = peak_timing_tdoa(signals[a], signals[b], t_axis,
                                    t_sw_a, t_sw_b, fs,
                                    search_start_s=-0.0001, search_end_s=0.0005)
        elif method == "xcorr":
            wa = extract_window(sig_sw[a], t_axis, ref["t_sw"], tw["xcorr_sw_window_pre_s"], tw["xcorr_sw_window_post_s"])
            wb = extract_window(sig_sw[b], t_axis, ref["t_sw"], tw["xcorr_sw_window_pre_s"], tw["xcorr_sw_window_post_s"])
            tau = xcorr_tdoa(wa, wb, fs, interp_factor=tw["interp_factor"], max_lag_s=max_lag)
        else:
            sw_wins_a = extract_window(signals[a], t_axis, ref["t_sw"],
                                        tw["sw_window_pre_s"], tw["sw_window_post_s"])
            sw_wins_b = extract_window(signals[b], t_axis, ref["t_sw"],
                                        tw["sw_window_pre_s"], tw["sw_window_post_s"])
            tau = gcc_phat(sw_wins_a, sw_wins_b, fs,
                            interp_factor=tw["interp_factor"], max_lag_s=max_lag)
        tdoa_sw[f"{a},{b}"] = tau

    tdoa_mb = {}
    mb_used = ref.get("mb_detected") and ref.get("mb_reliable")
    if mb_used:
        for (a, b) in pairs:
            if use_peak_timing:
                ch_a = prep["per_channel"][a]
                ch_b = prep["per_channel"][b]
                t_mb_a = ch_a.get("t_mb") if ch_a.get("mb_detected") else ref["t_mb"]
                t_mb_b = ch_b.get("t_mb") if ch_b.get("mb_detected") else ref["t_mb"]
                # MB: search from 0.5ms to 3ms AFTER onset -- the Friedlander peak
                # arrives after the onset threshold crossing, not right at it.
                tau = peak_timing_tdoa(signals[a], signals[b], t_axis,
                                        t_mb_a, t_mb_b, fs,
                                        search_start_s=0.0005, search_end_s=0.003)
            elif method == "xcorr":
                wa = extract_window(sig_mb[a], t_axis, ref["t_mb"], tw["xcorr_mb_window_pre_s"], tw["xcorr_mb_window_post_s"])
                wb = extract_window(sig_mb[b], t_axis, ref["t_mb"], tw["xcorr_mb_window_pre_s"], tw["xcorr_mb_window_post_s"])
                tau = xcorr_tdoa(wa, wb, fs, interp_factor=tw["interp_factor"], max_lag_s=max_lag)
            else:
                mb_wins_a = extract_window(signals[a], t_axis, ref["t_mb"],
                                            tw["mb_window_pre_s"], tw["mb_window_post_s"])
                mb_wins_b = extract_window(signals[b], t_axis, ref["t_mb"],
                                            tw["mb_window_pre_s"], tw["mb_window_post_s"])
                tau = gcc_phat(mb_wins_a, mb_wins_b, fs,
                                interp_factor=tw["interp_factor"], max_lag_s=max_lag)
            tdoa_mb[f"{a},{b}"] = tau
    else:
        print("[Det Stage 2] WARNING: reference channel's MB detection wasn't reliable -- "
              "skipping MB TDOA (would be meaningless windowed around a bad onset guess).")
    if not sw_used:
        print("[Det Stage 2] No shockwave detected -- MB TDOA only (bearing, no range).")
    if tdoa_sw:
        print("SW TDOA (us):")
    for k, v in tdoa_sw.items():
        print(f"  ({k}): {v*1e6:+.3f}")
    if tdoa_mb:
        print("MB TDOA (us):")
        for k, v in tdoa_mb.items():
            print(f"  ({k}): {v*1e6:+.3f}")

    # Consistency check: for a truly consistent set, tau(1,2) should ~= tau(0,2)-tau(0,1), etc.
    # Reported for transparency, not used to "correct" anything.
    consistency_sw = {}
    for (a, b) in (pairs if tdoa_sw else []):
        if a == 0:
            continue
        expected = tdoa_sw.get(f"0,{b}", 0) - tdoa_sw.get(f"0,{a}", 0)
        actual = tdoa_sw.get(f"{a},{b}", 0)
        consistency_sw[f"{a},{b}"] = float(actual - expected)

    out_json = f"{cfg['output']['basename']}_tdoa.json"
    dt_sw_to_mb = (ref["t_mb"] - ref["t_sw"]) if (mb_used and sw_used) else None
    meta = dict(
        generated_at_utc=datetime.now(timezone.utc).isoformat(),
        stage="tdoa",
        source_prepared_json=args.prepared_json,
        reference_channel=ref_ch,
        tdoa_method=("peak_timing" if use_peak_timing else method),
        sample_rate_hz=fs,
        n_channels=n_channels,
        sw_available=bool(tdoa_sw),
        mb_available=bool(tdoa_mb),
        event_class=prep.get("event_class"),
        dt_sw_to_mb_s=dt_sw_to_mb,
        dp_sw_ref_channel=ref.get("dP_sw"),
        # Stage 1 peak-normalizes the signal; this factor turns dp_* back into full-scale units
        normalization_peak=prep.get("normalization_peak"),
        t_n_sw_ref_channel=ref.get("T_N"),
        dp_mb_ref_channel=ref.get("dP_mb"),
        tdoa_sw_s=tdoa_sw,
        tdoa_mb_s=tdoa_mb,
        sw_consistency_residual_s=consistency_sw,
        config_used=cfg,
    )
    with open(out_json, "w") as f:
        json.dump(meta, f, indent=2)
    print(f"[Det Stage 2] Wrote {out_json}")


if __name__ == "__main__":
    main()
