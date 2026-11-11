#!/usr/bin/env python3
"""
gs_det_signal_prepare.py
===========================
DETECTION STAGE 1. Reads a multichannel .wav file (e.g. the output of
gs_gen_apply_adc.py, or a real recording) and BLINDLY detects the
shockwave (SW) and muzzle blast (MB) onset per channel -- no access to any
generation-block ground truth, only the waveform itself and det_config.ini.

Two detection methods ([detection] method in det_config.ini):

  classify (default) -- array-level event classifier (gs_det_events.py): every
     event is labelled SW or MB by its frequency content (1-4 kHz vs
     100-1000 Hz) and time/duration (impulsiveness, rise time, width). SW and
     MB are each optional and the MB does not have to follow a SW, so blast-only
     (subsonic, slowed bullet, shot fired away) and SW-only recordings are
     handled. Each channel's onset is then refined around the array-level time.
     Tuned on real recordings (Cipatat test range).

  threshold -- the original per-channel two-stage detector:
  1. SW: first sample where |signal| exceeds sw_threshold_mult x the noise
     RMS measured from the first noise_window_s of the recording.
  2. MB: energy-based search within [dt_min_s, dt_max_s] AFTER the SW
     onset, looking for the first sustained rise above
     mb_energy_threshold_mult x noise_RMS^2.

Optional Stage-1 deconvolution ([deconvolution] enabled=true): if a
calibration recording (the array's own response to a known sweep, played
once at install time) is available, its per-channel room response is
estimated via regularized deconvolution and used to deconvolve the live
signal BEFORE onset detection. This is the calibration-sweep approach that
was validated (in earlier work this pipeline is derived from) to recover
sub-10us TDOA even in the RT60-vs-event-gap overlap case that otherwise
breaks GCC-PHAT entirely -- legitimate for a fixed-position sensor (the
room doesn't change shot-to-shot), and does NOT require or use any
generation-block ground truth: only a real (or, for testing, simulated)
recording of the array's own response to a known test signal.

Usage
-----
  python gs_det_signal_prepare.py input.wav                     # det_config.ini
  python gs_det_signal_prepare.py input.wav my_det_config.ini
  python gs_det_signal_prepare.py input.wav my_det_config.ini -o my_result
  python gs_det_signal_prepare.py input.wav -mic_config=UMA8         # miniDSP UMA-8 v2 (4 of its mics)

-mic_config selects the installed microphone array ([mic_config.<NAME>] section of
det_config.ini; default [geometry] mic_config = Tetrahedral-0,3; use PILAR-V for the real
PILAR-V sensor, e.g. the Cipatat recordings). A recording with exactly N channels
(N = number of mics, e.g. the generator's output) is used as is. A recording with more channels
is reduced with the preset's channel_map: those .wav channels (0-based) are taken as mics
0..N-1 (e.g. channels 0,1,3,5 of an 8-channel UMA-8 RAW recording); without a channel_map the
first N channels are used. The selected preset is recorded in the
output .json, so Stages 2-4 pick it up without repeating the option.

Output
------
  <output>_prepared.wav   multichannel, float32, peak-normalized copy of the
                           input (downstream stages load this, not the
                           original file, so normalization is consistent) --
                           deconvolved first, if [deconvolution] is enabled
  <output>_prepared.json  per-channel SW/MB onset times + amplitudes, and
                           whether each channel's MB detection was reliable
"""

import argparse
import configparser
import json
import sys
from datetime import datetime, timezone

import numpy as np
from scipy.io import wavfile
from scipy.signal import chirp as make_chirp

import gs_det_events


def make_calibration_sweep(fs, f0, f1, duration):
    """The KNOWN test signal played during calibration -- must match
    whatever was actually played when the calibration recording was made
    (same f0/f1/duration/fs, all in det_config.ini, not derived from
    anything generation-side)."""
    t = np.arange(0, duration, 1.0 / fs)
    return make_chirp(t, f0=f0, t1=duration, f1=f1, method="logarithmic")


def estimate_rir_from_calibration(sweep, recorded, rir_len, reg=1e-2):
    """Regularized frequency-domain deconvolution: given the known sweep and
    the array's recorded response to it, estimate that channel's room
    impulse response. Standard swept-sine RIR measurement technique."""
    n_fft = len(recorded)
    Rec = np.fft.rfft(recorded, n=n_fft)
    Sw = np.fft.rfft(sweep, n=n_fft)
    G = Rec * np.conj(Sw) / (np.abs(Sw) ** 2 + reg * np.max(np.abs(Sw) ** 2))
    rir_full = np.fft.irfft(G, n=n_fft)
    return rir_full[:rir_len]


def wiener_deconvolve(y, rir, reg=1e-2):
    """Regularized inverse filter: remove the estimated room response from a
    live recording. Same technique used to estimate the RIR above, applied
    in reverse."""
    n = len(y)
    n_fft = n + len(rir) - 1
    Y = np.fft.rfft(y, n=n_fft)
    H = np.fft.rfft(rir, n=n_fft)
    H_mag2 = np.abs(H) ** 2
    G = np.conj(H) / (H_mag2 + reg * np.max(H_mag2))
    x_hat = np.fft.irfft(Y * G, n=n_fft)
    return x_hat[:n]


def apply_calibration_deconvolution(signals, fs, decon_cfg, channel_map=None):
    """Load the calibration recording, estimate each channel's RIR from it,
    and deconvolve the corresponding live channel. Returns the deconvolved
    signals list, or the original signals unchanged (with a warning) if the
    calibration file can't be used. channel_map: the recording channels used
    as mics (same mapping as the live recording), applied when the
    calibration recording has more channels than there are mics."""
    try:
        cal_fs, cal_data = wavfile.read(decon_cfg["calibration_wav"])
    except FileNotFoundError:
        print(f"WARNING: calibration_wav '{decon_cfg['calibration_wav']}' not found -- "
              f"skipping deconvolution, using raw signal.", file=sys.stderr)
        return signals
    if cal_data.ndim == 1:
        cal_data = cal_data[:, None]
    if cal_fs != fs:
        print(f"WARNING: calibration recording fs ({cal_fs}) != input fs ({fs}) -- "
              f"skipping deconvolution.", file=sys.stderr)
        return signals
    if (cal_data.shape[1] > len(signals) and channel_map is not None
            and max(channel_map) < cal_data.shape[1]):
        cal_data = cal_data[:, channel_map]
    if cal_data.shape[1] != len(signals):
        print(f"WARNING: calibration recording has {cal_data.shape[1]} channels, "
              f"input has {len(signals)} -- skipping deconvolution.", file=sys.stderr)
        return signals

    is_int = np.issubdtype(cal_data.dtype, np.integer)
    full_scale = float(np.iinfo(cal_data.dtype).max) if is_int else 1.0
    cal_signals = [cal_data[:, i].astype(np.float64) / full_scale for i in range(cal_data.shape[1])]

    sweep = make_calibration_sweep(fs, decon_cfg["sweep_f0_hz"], decon_cfg["sweep_f1_hz"],
                                    decon_cfg["sweep_duration_s"])
    rir_len = int(decon_cfg["rir_length_s"] * fs)

    print(f"[Det Stage 1] Deconvolving with calibration recording "
          f"'{decon_cfg['calibration_wav']}' ({len(signals)} channels)...")
    deconvolved = []
    for i in range(len(signals)):
        rir_est = estimate_rir_from_calibration(sweep, cal_signals[i], rir_len,
                                                  reg=decon_cfg["regularization"])
        deconvolved.append(wiener_deconvolve(signals[i], rir_est, reg=decon_cfg["regularization"]))
    return deconvolved


# Microphone-array presets of the INSTALLED array (the detector's own hardware -- legitimately
# known, see det_config.ini). Selected with -mic_config=<NAME> on the command line or
# [geometry] mic_config = <NAME> in det_config.ini. Each preset is a [mic_config.<NAME>] section;
# any key it does not set falls back to [geometry]. These built-in copies make the presets work
# with an .ini that predates them; a [mic_config.<NAME>] section in the .ini overrides them key
# by key. Deliberately an independent copy of the generator's presets (fairness boundary: no
# imports from gs_gen_*).
DEFAULT_MIC_CONFIG = "Tetrahedral-0,3"
MIC_CONFIG_PRESETS = {
    # Regular 4-mic tetrahedron, 0.30 m edge (formerly "PilarV"; simulation default)
    "Tetrahedral-0,3": {"array_type": "tetrahedron", "l_array": "0.30"},
    # Real PILAR-V sensor ("Trackball"; Cipatat recordings): ring of 3 mics, radius
    # 0.2055 m at 0/120/240 deg + 1 mic 0.2412 m above. Frame x = north, y = west, z = up.
    "PILAR-V": {"array_type": "custom",
                "mic_positions": "0.2055,0,0; -0.1028,0.178,0; -0.1028,-0.178,0; 0,0,0.2412"},
    # miniDSP UMA-8 v2, 4 of its 7 mics: the centre mic (#0, origin) + ring mics #1, #3, #5
    # (42.5 mm radius, 120 deg apart: an equilateral triangle, 73.6 mm side), all at z = 0.
    # RAW firmware gives 8 channels: channel_map picks recording channels 0, 1, 3, 5.
    "UMA8": {"array_type": "custom",
             "mic_positions": "0,0,0; 0,0.0425,0; 0.036806,-0.02125,0; -0.036806,-0.02125,0",
             "channel_map": "0,1,3,5",
             "planar_source_side": "above"},
}
MIC_CONFIG_KEYS = ("array_type", "l_array", "mic_positions", "mic_positions_file",
                   "channel_map", "planar_source_side")


def add_mic_config_argument(parser):
    """Add the -mic_config option (also --mic_config / --mic-config) to a stage's argparse parser."""
    parser.add_argument("-mic_config", "--mic_config", "--mic-config", dest="mic_config", default=None,
                        metavar="NAME",
                        help="installed microphone array preset: PILAR-V (the real sensor), "
                             "Tetrahedral-0,3 (regular tetrahedron, 0.30 m edge) or UMA8 "
                             "(miniDSP UMA-8 v2: centre mic + ring mics #1, #3, #5) -- any [mic_config.NAME] section of the .ini. "
                             "Default: the preset recorded in the upstream stage's .json (Stages 2-4), "
                             "else [geometry] mic_config in the .ini (Tetrahedral-0,3).")


def upstream_mic_config(meta):
    """mic_config recorded in an upstream DETECTION stage's .json (None if absent)."""
    try:
        return meta["config_used"]["geometry"].get("mic_config")
    except (KeyError, TypeError, AttributeError):
        return None


def _parse_xyz_rows(text):
    rows = []
    for row in text.split(";"):
        if not row.strip():
            continue
        v = [float(p) for p in row.split(",") if p.strip() != ""]
        if len(v) != 3:
            raise ValueError(f"mic_positions: expected 'x,y,z' per mic, got '{row.strip()}'")
        rows.append(v)
    return rows


LEGACY_MIC_CONFIG_NAMES = {"pilarv": "Tetrahedral-0,3"}   # renamed Nov 2026


def _resolve_mic_config(cp, name):
    """Geometry dict of the [mic_config.<name>] preset (case-insensitive), falling back to [geometry]."""
    sections = {s.split(".", 1)[1].lower(): s for s in cp.sections() if s.startswith("mic_config.")}
    key = str(name).strip().lower()
    if key in LEGACY_MIC_CONFIG_NAMES and key not in sections:
        name = LEGACY_MIC_CONFIG_NAMES[key]
        print(f"NOTE: mic_config 'PilarV' was renamed to '{name}' (the real sensor is now 'PILAR-V').",
              file=sys.stderr)
        key = name.lower()
    if key not in sections:
        valid = ", ".join(s.split(".", 1)[1] for s in sections.values())
        raise ValueError(f"unknown mic_config '{name}' -- valid: {valid} "
                         f"(add a [mic_config.{name}] section to the .ini to define it)")
    section = sections[key]
    raw = {k: (cp.get(section, k) if cp.has_option(section, k) else cp.get("geometry", k))
           for k in MIC_CONFIG_KEYS}

    array_type = raw["array_type"].strip().lower()
    if array_type not in ("tetrahedron", "custom"):
        raise ValueError(f"mic_config {name}: array_type must be 'tetrahedron' or 'custom'")
    mic_positions = None
    if array_type == "custom":
        if raw["mic_positions"].strip():
            mic_positions = _parse_xyz_rows(raw["mic_positions"])
        elif raw["mic_positions_file"].strip():
            mic_positions = np.loadtxt(raw["mic_positions_file"].strip(), delimiter=",",
                                       comments="#", ndmin=2).tolist()
        if not mic_positions or len(mic_positions) < 2:
            raise ValueError(f"mic_config {name}: array_type = custom needs mic_positions "
                             f"or mic_positions_file with >= 2 mics")
    n_mics = len(mic_positions) if mic_positions else 4
    channel_map = [int(c) for c in raw["channel_map"].split(",") if c.strip() != ""] or None
    if channel_map is not None and len(channel_map) != n_mics:
        raise ValueError(f"mic_config {name}: channel_map has {len(channel_map)} entries "
                         f"but the array has {n_mics} mics")
    side = raw["planar_source_side"].strip().lower()
    if side not in ("above", "below"):
        raise ValueError(f"mic_config {name}: planar_source_side must be 'above' or 'below'")

    geometry = dict(mic_config=section.split(".", 1)[1], array_type=array_type,
                    l_array=float(raw["l_array"]), mic_positions=mic_positions,
                    channel_map=channel_map, planar_source_side=side, n_mics=n_mics)
    pos = build_mic_positions(dict(geometry=geometry))
    geometry["mic_positions_used"] = pos.tolist()
    geometry["aperture_m"] = array_aperture(pos, geometry)
    return geometry


def _parse_band(text):
    """'1000,4000' -> (1000.0, 4000.0)"""
    lo, hi = (float(v) for v in text.split(","))
    if not 0 < lo < hi:
        raise ValueError(f"band '{text}': need 0 < low < high (Hz)")
    return (lo, hi)


DETECTION_METHODS = ("classify", "threshold")


def load_det_config(path, mic_config=None):
    """mic_config: installed-array preset name ([mic_config.<NAME>] section); None uses
    [geometry] mic_config."""
    cp = configparser.ConfigParser()
    cp.read_dict({
        "geometry": {"mic_config": DEFAULT_MIC_CONFIG, "l_array": "0.30", "array_type": "tetrahedron",
                     "mic_positions": "", "mic_positions_file": "", "channel_map": "",
                     "planar_source_side": "above"},
        **{f"mic_config.{k}": v for k, v in MIC_CONFIG_PRESETS.items()},
        "physics": {"speed_of_sound": "343.0", "atmospheric_pressure_pa": "101325.0"},
        "detection": {"method": "classify",
                      "sw_threshold_mult": "10.0", "mb_energy_threshold_mult": "20.0",
                      "dt_min_s": "0.01", "dt_max_s": "5.0", "noise_window_s": "0.005",
                      "sw_band_hz": "1000,4000", "mb_bands_hz": "150,1000",
                      "sw_min_hf_snr_db": "25.0", "sw_min_crest_db": "14.0",
                      "sw_max_rise_ms": "0.08", "sw_max_lf_width_ms": "0.0",
                      "mb_min_lf_snr_db": "11.0", "mb_max_hf_excess_db": "10.0", "mb_min_prominence_db": "6.0",
                      "mb_in_tail": "false"},
        "tdoa": {"interp_factor": "16", "sw_window_pre_s": "0.003", "sw_window_post_s": "0.015",
                 "mb_window_pre_s": "0.006", "mb_window_post_s": "0.080", "max_lag_s": "0.001"},
        "bullet_reference_library": {"calibers": "5.56_NATO,7.62_NATO",
                                      "5.56_NATO_d": "0.0057", "5.56_NATO_L": "0.023", "5.56_NATO_v0": "940",
                                      "7.62_NATO_d": "0.00782", "7.62_NATO_L": "0.0285", "7.62_NATO_v0": "850"},
        "calibration": {"assumed_peak_pa_at_full_scale": "1520.0"},
        "deconvolution": {"enabled": "false", "calibration_wav": "",
                          "sweep_f0_hz": "40.0", "sweep_f1_hz": "20000.0",
                          "sweep_duration_s": "1.0", "rir_length_s": "2.0",
                          "regularization": "0.01"},
        "output": {"basename": "det_result"},
    })
    if path is not None:
        found = cp.read(path)
        if not found:
            print(f"WARNING: config file '{path}' not found -- using built-in defaults.", file=sys.stderr)

    calibers = [c.strip() for c in cp.get("bullet_reference_library", "calibers").split(",")]
    def _opt(cal, key):
        k = f"{cal}_{key}"
        return cp.getfloat("bullet_reference_library", k) if cp.has_option("bullet_reference_library", k) else None
    bullet_library = {
        cal: dict(d=cp.getfloat("bullet_reference_library", f"{cal}_d"),
                  L=cp.getfloat("bullet_reference_library", f"{cal}_L"),
                  v0=_opt(cal, "v0"), lam=_opt(cal, "lam"))
        for cal in calibers
    }

    geometry = _resolve_mic_config(cp, mic_config or cp.get("geometry", "mic_config"))
    # Full-scale calibration: a per-array override (assumed_peak_pa_at_full_scale.<mic_config>)
    # wins over the generic value.
    fs_key = f"assumed_peak_pa_at_full_scale.{geometry.get('mic_config', '')}"
    fs_pa = (cp.getfloat("calibration", fs_key) if cp.has_option("calibration", fs_key)
             else cp.getfloat("calibration", "assumed_peak_pa_at_full_scale"))

    return {
        "geometry": geometry,
        "physics": dict(speed_of_sound=cp.getfloat("physics", "speed_of_sound"),
                        atmospheric_pressure_pa=cp.getfloat("physics", "atmospheric_pressure_pa")),
        "detection": dict(method=cp.get("detection", "method").strip().lower(),
                           sw_band_hz=_parse_band(cp.get("detection", "sw_band_hz")),
                           mb_bands_hz=[_parse_band(b) for b in cp.get("detection", "mb_bands_hz").split(";")
                                        if b.strip()],
                           sw_min_hf_snr_db=cp.getfloat("detection", "sw_min_hf_snr_db"),
                           sw_min_crest_db=cp.getfloat("detection", "sw_min_crest_db"),
                           sw_max_rise_ms=cp.getfloat("detection", "sw_max_rise_ms"),
                           sw_max_lf_width_ms=cp.getfloat("detection", "sw_max_lf_width_ms"),
                           mb_min_lf_snr_db=cp.getfloat("detection", "mb_min_lf_snr_db"),
                           mb_max_hf_excess_db=cp.getfloat("detection", "mb_max_hf_excess_db"),
                           mb_min_prominence_db=cp.getfloat("detection", "mb_min_prominence_db"),
                           mb_in_tail=cp.getboolean("detection", "mb_in_tail"),
                           sw_threshold_mult=cp.getfloat("detection", "sw_threshold_mult"),
                           mb_energy_threshold_mult=cp.getfloat("detection", "mb_energy_threshold_mult"),
                           dt_min_s=cp.getfloat("detection", "dt_min_s"),
                           dt_max_s=cp.getfloat("detection", "dt_max_s"),
                           noise_window_s=cp.getfloat("detection", "noise_window_s")),
        "tdoa": dict(interp_factor=cp.getint("tdoa", "interp_factor"),
                     sw_window_pre_s=cp.getfloat("tdoa", "sw_window_pre_s"),
                     sw_window_post_s=cp.getfloat("tdoa", "sw_window_post_s"),
                     mb_window_pre_s=cp.getfloat("tdoa", "mb_window_pre_s"),
                     mb_window_post_s=cp.getfloat("tdoa", "mb_window_post_s"),
                     max_lag_s=cp.getfloat("tdoa", "max_lag_s")),
        "bullet_reference_library": bullet_library,
        "calibration": dict(assumed_peak_pa_at_full_scale=fs_pa),
        "deconvolution": dict(enabled=cp.getboolean("deconvolution", "enabled"),
                               calibration_wav=cp.get("deconvolution", "calibration_wav"),
                               sweep_f0_hz=cp.getfloat("deconvolution", "sweep_f0_hz"),
                               sweep_f1_hz=cp.getfloat("deconvolution", "sweep_f1_hz"),
                               sweep_duration_s=cp.getfloat("deconvolution", "sweep_duration_s"),
                               rir_length_s=cp.getfloat("deconvolution", "rir_length_s"),
                               regularization=cp.getfloat("deconvolution", "regularization")),
        "output": dict(basename=cp.get("output", "basename")),
    }


def make_tetrahedron(edge_length):
    raw = np.array([[0.000, 0.000, 1.000], [0.000, 0.943, -0.333],
                     [-0.816, -0.471, -0.333], [0.816, -0.471, -0.333]], dtype=float)
    scale = edge_length / np.linalg.norm(raw[0] - raw[1])
    return raw * scale


def build_mic_positions(cfg):
    """(N, 3) mic positions (m) of the installed array selected by mic_config."""
    g = cfg["geometry"]
    if g["array_type"] == "custom":
        pos = np.asarray(g["mic_positions"], dtype=float)
        if pos.ndim != 2 or pos.shape[1] != 3 or len(pos) < 2:
            raise ValueError("mic positions must be an (N, 3) array with N >= 2")
        return pos
    return make_tetrahedron(g["l_array"])


def array_aperture(mic_pos, geometry=None):
    """Largest mic-to-mic distance (m) -- aperture / c is the largest possible TDOA.
    For the tetrahedron this is l_array (its nominal edge)."""
    if geometry is not None and geometry.get("array_type") == "tetrahedron":
        return float(geometry["l_array"])
    p = np.asarray(mic_pos, dtype=float)
    return float(np.max(np.linalg.norm(p[:, None, :] - p[None, :, :], axis=-1)))


def find_settled_point(signal, sw_idx, fs, noise_rms, settle_mult=3.0,
                        guard_s=0.005, dwell_s=0.005, max_search_s=2.0):
    """Find where the post-SW energy settles back near the noise floor,
    WITHOUT assuming any particular reverb time -- a fixed dt_min gap (the
    original, simpler design) silently breaks whenever the true RT60 is
    long enough to still be ringing at that fixed gap, and a detector isn't
    allowed to just know the true RT60 (that's generation-block ground
    truth). This scans forward from a guard interval after the SW onset,
    looking for the first window whose RMS energy drops and STAYS below
    settle_mult x noise_rms for at least dwell_s -- adaptive, no assumed
    reverb time required."""
    win_step = max(1, int(0.001 * fs))
    dwell_samples = max(1, int(dwell_s * fs))
    start_idx = sw_idx + int(guard_s * fs)
    end_idx = min(len(signal), sw_idx + int(max_search_s * fs))
    threshold = settle_mult * noise_rms

    k = start_idx
    while k < end_idx - dwell_samples:
        window = signal[k:k + dwell_samples]
        if np.sqrt(np.mean(window ** 2)) < threshold:
            return k
        k += win_step
    return start_idx   # couldn't find a settled point -- fall back to the guard time alone


def detect_events(signal, t_axis, fs, sw_thresh_mult, mb_energy_thresh_mult,
                   dt_min, dt_max, noise_window_s):
    """Blind two-stage detector. Returns dict with t_sw, t_mb, dP_sw, T_N,
    dP_mb, noise_rms, mb_reliable (False if no settled point / no MB rise
    could be found -- signals a genuinely hard case, not a silent failure).

    dt_min is used as a GUARD floor (minimum gap regardless of settling),
    not the sole gate -- the actual search start additionally waits for
    find_settled_point() above, so this works whether the true reverb time
    is short or long, without needing to know which."""
    n_noise = int(noise_window_s * fs)
    noise_rms = float(np.sqrt(np.mean(signal[:n_noise] ** 2)))

    thresh_sw = sw_thresh_mult * noise_rms
    above = np.abs(signal) > thresh_sw
    crossings_sw = np.where(np.diff(above.astype(int)) > 0)[0]
    if len(crossings_sw) == 0:
        return None
    sw_idx = crossings_sw[0] + 1
    t_sw_onset = t_axis[sw_idx]

    win_end = min(len(signal) - 1, sw_idx + int(0.010 * fs))
    sw_win = signal[sw_idx:win_end]
    t_sw_win = t_axis[sw_idx:win_end]
    # Refine: the raw threshold-crossing sample can land on a noise blip a
    # few samples before the true N-wave rise (its leading edge is the
    # steepest, most reliable feature) -- snap to the local peak within a
    # short window instead of trusting the crossing sample directly.
    refine_win = min(len(sw_win), int(0.0003 * fs))   # ~0.3ms, wide enough for typical T_N
    if refine_win > 1:
        peak_offset = int(np.argmax(np.abs(sw_win[:refine_win])))
        if peak_offset > 0:
            sw_idx = sw_idx + peak_offset
            t_sw_onset = t_axis[sw_idx]
            win_end = min(len(signal) - 1, sw_idx + int(0.010 * fs))
            sw_win = signal[sw_idx:win_end]
            t_sw_win = t_axis[sw_idx:win_end]
    dP_measured = float(sw_win[0]) if len(sw_win) else float("nan")
    zc = np.where(np.diff(np.sign(sw_win)))[0]
    T_N_meas = float(2.0 * (t_sw_win[zc[0]] - t_sw_win[0])) if len(zc) >= 1 else float("nan")

    settled_idx = find_settled_point(signal, sw_idx, fs, noise_rms, max_search_s=dt_max)
    search_start_idx = max(settled_idx, sw_idx + int(dt_min * fs))
    t_search_start = t_axis[search_start_idx]
    t_search_end = t_sw_onset + dt_max
    mask_search = (t_axis >= t_search_start) & (t_axis <= t_search_end)
    sig_search = signal[mask_search]
    t_search = t_axis[mask_search]
    if len(sig_search) < 10:
        return dict(t_sw=t_sw_onset, t_mb=None, dP_sw=dP_measured, T_N=T_N_meas,
                    dP_mb=None, noise_rms=noise_rms, mb_reliable=False)

    win_e = max(1, int(0.001 * fs))
    energy = np.array([np.sum(sig_search[k:k + win_e] ** 2) / win_e
                        for k in range(len(sig_search) - win_e)])
    above_mb = energy > mb_energy_thresh_mult * noise_rms ** 2
    mb_reliable = True
    if len(above_mb) == 0:
        return dict(t_sw=t_sw_onset, t_mb=None, dP_sw=dP_measured, T_N=T_N_meas,
                    dP_mb=None, noise_rms=noise_rms, mb_reliable=False)
    if above_mb[0]:
        mb_idx = 0
        mb_reliable = False   # still above threshold right at the settled point -- genuinely hard case
    else:
        crossings_mb = np.where(np.diff(above_mb.astype(int)) > 0)[0]
        if len(crossings_mb) == 0:
            return dict(t_sw=t_sw_onset, t_mb=None, dP_sw=dP_measured, T_N=T_N_meas,
                        dP_mb=None, noise_rms=noise_rms, mb_reliable=False)
        mb_idx = crossings_mb[0] + 1
    t_mb_onset = t_search[mb_idx]

    mb_win_idx = np.searchsorted(t_axis, t_mb_onset)
    mb_win = signal[mb_win_idx:mb_win_idx + int(0.005 * fs)]
    dP_mb_est = float(mb_win.max()) if len(mb_win) else float("nan")

    return dict(t_sw=float(t_sw_onset), t_mb=float(t_mb_onset), dP_sw=dP_measured,
                T_N=T_N_meas, dP_mb=dP_mb_est, noise_rms=noise_rms, mb_reliable=mb_reliable)


def _jsonable(d):
    """Copy of an event-feature dict with plain Python numbers (for the .json)."""
    if d is None:
        return None
    out = {}
    for k, v in d.items():
        if isinstance(v, (np.floating, np.integer)):
            v = v.item()
        out[k] = v
    return out


def classify_channels(normalized, t_axis, fs, cfg, noise_window_s):
    """[detection] method = classify: array-level SW/MB classification (gs_det_events),
    then each channel's onset/amplitude refined around the array-level event.
    Returns (per_channel list in the same format as the threshold method, array summary)."""
    det = cfg["detection"]
    params = {k: det[k] for k in gs_det_events.DEFAULT_PARAMS}
    x = np.stack(normalized, axis=1)
    ev = gs_det_events.detect_sw_mb(x, fs, params)
    max_tdoa = cfg["geometry"]["aperture_m"] / cfg["physics"]["speed_of_sound"]

    per_channel = []
    for i, sig in enumerate(normalized):
        n_noise = max(1, int(noise_window_s * fs))
        noise_rms = float(np.sqrt(np.mean(sig[:n_noise] ** 2)))
        ch = dict(channel=i, sw_detected=False, mb_detected=False, mb_reliable=False,
                  t_sw=None, t_mb=None, dP_sw=None, T_N=None, dP_mb=None, noise_rms=noise_rms)
        if ev["t_sw"] is not None:
            k = gs_det_events.channel_onset(sig, fs, ev["sw_idx"], max_tdoa + 0.0005, 1e-4,
                                            det["sw_band_hz"])
            # same refinement as the threshold method: snap to the local peak (~0.3 ms)
            r = max(1, int(0.0003 * fs))
            k = k + int(np.argmax(np.abs(sig[k:k + r])))
            win_end = min(len(sig) - 1, k + int(0.010 * fs))
            sw_win = sig[k:win_end]
            zc = np.where(np.diff(np.sign(sw_win)))[0]
            ch.update(sw_detected=True, t_sw=float(t_axis[k]),
                      dP_sw=float(sw_win[0]) if len(sw_win) else None,
                      T_N=float(2.0 * (t_axis[k + zc[0]] - t_axis[k])) if len(zc) else None)
        if ev["t_mb"] is not None:
            # search back as far as the array-level onset lies before its peak (+ max TDOA)
            back = (ev["mb_idx"] / fs - ev["t_mb"]) + max_tdoa + 0.002
            k = gs_det_events.channel_onset(sig, fs, ev["mb_idx"], max_tdoa + 0.005, 1e-3,
                                            ev["mb_band"], back_s=back, causal=True)
            mb_win = sig[k:k + int(0.005 * fs)]
            ch.update(mb_detected=True, mb_reliable=True, t_mb=float(t_axis[k]),
                      dP_mb=float(mb_win.max()) if len(mb_win) else None)
        per_channel.append(ch)

    summary = dict(
        event_class=ev["event_class"], t_sw_s=ev["t_sw"], t_mb_s=ev["t_mb"],
        dt_sw_to_mb_s=(ev["t_mb"] - ev["t_sw"]) if ev["event_class"] == "SW+MB" else None,
        n_sw_in_burst=ev["n_sw"], sw_tail_end_s=ev["t_sw_tail_end"],
        mb_band_hz=list(ev["mb_band"]) if ev["mb_band"] else None,
        mb_may_be_hidden_in_sw_tail=bool(ev["mb_hidden_possible"]),
        sw_features=_jsonable(ev["sw_features"]), mb_features=_jsonable(ev["mb_features"]),
        events=[_jsonable({k: e[k] for k in ("t", "cls", "snr_hf", "snr_lf", "crest", "rise_ms",
                                              "lf_width_ms", "hf_excess", "prominence")})
                for e in ev["events"] if e["cls"] not in ("-", "SW-tail")],
    )
    return per_channel, summary


def main():
    parser = argparse.ArgumentParser(
        description="Detection Stage 1: blind SW/MB onset detection and normalization.")
    parser.add_argument("input_wav", help="Path to the multichannel recording (.wav)")
    parser.add_argument("config", nargs="?", default="det_config.ini",
                         help="Path to detection .ini config (default: det_config.ini)")
    parser.add_argument("-o", "--output", default=None,
                         help="Output basename (overrides output.basename in config)")
    parser.add_argument("--calibration-wav", default=None,
                         help="Path to a calibration recording (overrides det_config.ini's "
                              "deconvolution.calibration_wav). A calibration recording is the "
                              "array's own response to a known log-sweep played at install time "
                              "(same sweep parameters as det_config.ini [deconvolution]). "
                              "Enables room-response deconvolution before onset detection -- "
                              "the single most effective fix for reverb corrupting TDOA.")
    add_mic_config_argument(parser)
    args = parser.parse_args()

    cfg = load_det_config(args.config, args.mic_config)
    if args.output:
        cfg["output"]["basename"] = args.output
    if args.calibration_wav:
        cfg["deconvolution"]["enabled"] = True
        cfg["deconvolution"]["calibration_wav"] = args.calibration_wav

    fs, data = wavfile.read(args.input_wav)
    if data.ndim == 1:
        data = data[:, None]

    # --- Map recording channels -> mics of the selected mic_config ---------------
    geom = cfg["geometry"]
    n_mics = geom["n_mics"]
    n_in = data.shape[1]
    # A recording with exactly n_mics channels (e.g. the generator's output) is used as is;
    # channel_map only applies to a recording with MORE channels (e.g. 8-channel UMA-8 RAW).
    channel_map = geom["channel_map"]
    if n_in < n_mics:
        raise SystemExit(f"[Det Stage 1] ERROR: {args.input_wav} has {n_in} channels but "
                         f"mic_config {geom['mic_config']} has {n_mics} mics.")
    if n_in == n_mics:
        channel_map = list(range(n_mics))
    elif channel_map is None:
        print(f"[Det Stage 1] NOTE: {n_in} channels in the recording, using the first {n_mics} "
              f"for mic_config {geom['mic_config']} (set channel_map to choose others).")
        channel_map = list(range(n_mics))
    elif max(channel_map) >= n_in or min(channel_map) < 0:
        raise SystemExit(f"[Det Stage 1] ERROR: channel_map {channel_map} of mic_config "
                         f"{geom['mic_config']} does not fit the {n_in}-channel {args.input_wav}.")
    data = data[:, channel_map]
    n_channels = data.shape[1]
    print(f"[Det Stage 1] mic_config = {geom['mic_config']} ({n_mics} mics, "
          f"aperture {geom['aperture_m'] * 1e3:.1f} mm), recording channels {channel_map}")
    is_int = np.issubdtype(data.dtype, np.integer)
    full_scale = float(np.iinfo(data.dtype).max) if is_int else 1.0
    signals = [data[:, i].astype(np.float64) / full_scale for i in range(n_channels)]

    # --- Optional deconvolution (before normalization and onset detection) ----
    decon_applied = False
    if cfg["deconvolution"]["enabled"]:
        if not cfg["deconvolution"]["calibration_wav"]:
            print("[Det Stage 1] WARNING: deconvolution.enabled=true but no "
                  "calibration_wav set -- skipping deconvolution.", file=sys.stderr)
        else:
            signals = apply_calibration_deconvolution(signals, fs, cfg["deconvolution"], channel_map)
            decon_applied = True

    peak = max(np.abs(s).max() for s in signals)
    normalized = [s / peak if peak > 0 else s for s in signals]

    t_axis = np.arange(len(normalized[0])) / fs
    det = cfg["detection"]

    decon_note = " (with room deconvolution)" if decon_applied else ""
    print(f"[Det Stage 1] Detecting SW/MB onsets in {args.input_wav} "
          f"({n_channels} channels, {fs} Hz){decon_note} -- BLIND")

    # --- SW consistency check: do all channels agree on direction? -----------
    # A strongly inconsistent set (cross-channel TDOA spread >> expected max
    # from array geometry) is a sign the SW detection caught reverb or noise
    # rather than the real wavefront -- flagged here so downstream stages can
    # weight accordingly, rather than silently computing a wrong bearing.
    def sw_consistency_flag(t_sw_list, fs_rate, aperture, c):
        valid = [t for t in t_sw_list if t is not None]
        if len(valid) < 2:
            return False, float("nan")
        spread_us = (max(valid) - min(valid)) * 1e6
        max_expected_us = aperture / c * 1e6   # aperture / c = max possible TDOA
        return spread_us > 3 * max_expected_us, spread_us

    if det["method"] not in DETECTION_METHODS:
        raise SystemExit(f"[Det Stage 1] ERROR: [detection] method must be one of {DETECTION_METHODS}")
    print(f"[Det Stage 1] Detection method: {det['method']}")
    array_summary = None
    per_channel = []
    t_sw_all = []
    if det["method"] == "classify":
        per_channel, array_summary = classify_channels(normalized, t_axis, fs, cfg, det["noise_window_s"])
        a = array_summary
        f = lambda v: "--" if v is None else f"{v:.5f}s"
        print(f"[Det Stage 1] Array events: {a['event_class']}  (SW {f(a['t_sw_s'])}, MB {f(a['t_mb_s'])}"
              + (f", MB band {a['mb_band_hz'][0]:g}-{a['mb_band_hz'][1]:g} Hz" if a["mb_band_hz"] else "")
              + (f", burst of {a['n_sw_in_burst']} SWs" if a["n_sw_in_burst"] > 1 else "") + ")")
        if a["event_class"] == "SW only" and a["mb_may_be_hidden_in_sw_tail"]:
            print("  NOTE: no MB after the SW's echo tail; an MB inside the tail can't be told apart from "
                  "echoes (set mb_in_tail = true for very close shots).")
        for ch in per_channel:
            print(f"  Ch{ch['channel']}: t_sw={f(ch['t_sw'])}  t_mb={f(ch['t_mb'])}")
            t_sw_all.append(ch["t_sw"])
    for i in (range(n_channels) if det["method"] == "threshold" else []):
        result = detect_events(normalized[i], t_axis, fs,
                                det["sw_threshold_mult"], det["mb_energy_threshold_mult"],
                                det["dt_min_s"], det["dt_max_s"], det["noise_window_s"])
        if result is None:
            print(f"  Ch{i}: SW NOT DETECTED", file=sys.stderr)
            per_channel.append(dict(channel=i, sw_detected=False))
            t_sw_all.append(None)
            continue
        flag = "" if result["mb_reliable"] else "  <-- MB UNRELIABLE"
        print(f"  Ch{i}: t_sw={result['t_sw']:.5f}s  t_mb={result['t_mb']}s{flag}")
        per_channel.append(dict(channel=i, sw_detected=True, mb_detected=result["t_mb"] is not None,
                                 **result))
        t_sw_all.append(result["t_sw"])

    sw_inconsistent, sw_spread_us = sw_consistency_flag(
        t_sw_all, fs, cfg["geometry"]["aperture_m"], cfg["physics"]["speed_of_sound"])
    if sw_inconsistent:
        print(f"  WARNING: SW spread across channels = {sw_spread_us:.1f} us, "
              f"exceeds 3x the array's maximum possible TDOA "
              f"({cfg['geometry']['aperture_m']/cfg['physics']['speed_of_sound']*1e6:.1f} us) -- "
              f"SW detection may have caught reverb or noise, not the true wavefront. "
              f"Consider enabling [deconvolution] in det_config.ini.")

    wav_path = f"{cfg['output']['basename']}_prepared.wav"
    json_path = f"{cfg['output']['basename']}_prepared.json"

    out_data = np.stack(normalized, axis=1).astype(np.float32)
    wavfile.write(wav_path, fs, out_data)

    meta = dict(
        generated_at_utc=datetime.now(timezone.utc).isoformat(),
        stage="signal_prepared",
        source_wav=args.input_wav,
        mic_config=geom["mic_config"],
        source_channels_used=channel_map,
        prepared_wav=wav_path,
        deconvolution_applied=decon_applied,
        sample_rate_hz=fs,
        n_channels=n_channels,
        n_samples=len(normalized[0]),
        normalization_peak=float(peak),
        sw_spread_us=float(sw_spread_us),
        sw_consistency_flag=bool(sw_inconsistent),
        detection_method=det["method"],
        event_class=(array_summary["event_class"] if array_summary else
                     ("SW+MB" if any(c.get("mb_detected") and c.get("mb_reliable") for c in per_channel)
                      else "SW only") if any(c.get("sw_detected") for c in per_channel) else "none"),
        array_events=array_summary,
        config_used=cfg,
        per_channel=per_channel,
    )
    with open(json_path, "w") as f:
        json.dump(meta, f, indent=2)

    print(f"[Det Stage 1] Wrote {wav_path}")
    print(f"[Det Stage 1] Wrote {json_path}")


if __name__ == "__main__":
    main()
