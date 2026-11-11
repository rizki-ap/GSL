#!/usr/bin/env python3
"""
util_gen_cipatat_clean.py
=========================
Clean (noise-free, reverberation-free) twins of the Cipatat development recordings.

For every development shot (01-06, 06a, 15, 16, 17) the scene of the reference table
(sound/Cipatat/metadata.md: calibre, azimuth, elevation, distance, miss distance, air
temperature) is rebuilt with the generator's Stage 1 (gs_gen_clean_signal.py: Whitham N-wave
with exponential drag, Friedlander muzzle blast, ISO 9613-1 absorption, PILAR-V geometry) and
fitted to the real recording:

  1. Bullet side: the trajectory passes the array on the left or the right at the miss
     distance. The side whose shockwave TDOAs match the measured ones is used.
  2. Levels and durations: shockwave peak and N-wave duration, muzzle-blast peak and
     positive-phase duration are scaled to the values measured on the recording
     (shockwave / muzzle_blast amp_scale and dur_scale). Each round of a burst is fitted on
     its own.
  3. Timing: each shockwave and each muzzle blast is placed at its measured arrival time
     (mean over the 4 channels). The TDOAs between channels are the model's.

Nothing else is added: no noise, no reverberation, no ground reflection, no echoes. The
signal is synthesised at 10 x 51 200 Hz and decimated (band-limited, like the real ADC).

Output (one per recording, same length and format as the real cut):
  <out>/shot-NN_audio_4ch_clean.wav   4 ch, 51 200 Hz, PCM_24, Cipatat full-scale units (Pa / 1520)
  <out>/shot-NN_audio_4ch_clean.json  scene, fitted scales, measured and model values, ground truth
  <out>/clean_summary.csv             one line per event (measured, model and fitted values)
  <out>/clean_vs_real_overview.png    real (high-passed) and clean signal, channel 1

Usage
  python util_gen_cipatat_clean.py                       # ../sound/Cipatat -> ../sound/Cipatat_clean
  python util_gen_cipatat_clean.py --shots 01 04 17
  python util_gen_cipatat_clean.py --no-fit              # pure model levels (no amplitude/duration fit)
  python util_gen_cipatat_clean.py --fs-pa 1520 --in DIR --out DIR

Needs gs_gen_physic.py, gs_gen_clean_signal.py, util_signal_similarity.py, gs_det_events.py
and gen_config.ini in this folder.
"""
import argparse
import copy
import csv
import json
import os
import sys
import warnings
from datetime import datetime, timezone

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, resample_poly, sosfiltfilt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import gs_gen_physic as gp                      # noqa: E402
from gs_gen_clean_signal import generate_clean  # noqa: E402
import util_signal_similarity as us             # noqa: E402

FS = 51200
UP = 10                       # synthesis at UP x FS, then decimated
FS_PA = 1520.0                # Pa per full scale (estimate, sound/Cipatat/metadata.md)

# Reference scene per recording (sound/Cipatat/metadata.md) and approximate event times in
# the cut file (shockwave, muzzle blast; refined on the signal). Azimuth clockwise from
# north, array frame x = north, y = west, z = up.
SCENES = {
    "01":  dict(weapon="M16", cal="5.56x45", az=269, el=0, D=306, b=20, T=26,
                ev=[(1.6993, 2.1443)]),
    "02":  dict(weapon="AK47", cal="7.62x39", az=85, el=-1, D=306, b=20, T=9,
                ev=[(1.6755, 2.0207)]),
    "03":  dict(weapon="Dragunov", cal="7.62x54R", az=355, el=6, D=500, b=10, T=20,
                ev=[(2.4329, 3.0419)]),
    "04":  dict(weapon="MINIMI (burst)", cal="7.62x51", az=1, el=1, D=304, b=30, T=14,
                ev=[(2.5435, 2.9822), (2.6196, 3.0583), (2.6885, 3.1332), (2.7703, 3.2081),
                    (2.8417, 3.2809)]),
    "05":  dict(weapon="AX338", cal="8.6x70", az=0, el=7, D=517, b=35, T=27,
                ev=[(2.0421, 2.7470)]),
    "06":  dict(weapon="STEYR", cal="12.7x99", az=11, el=2, D=764, b=30, T=35,
                ev=[(2.0427, 3.1400)]),
    "06a": dict(weapon="STEYR", cal="12.7x99", az=11, el=2, D=764, b=30, T=35,
                ev=[(10.0427, None)], mb_from="06"),
    "15":  dict(weapon="AK47", cal="7.62x39", az=0, el=-10, D=300, b=20, T=33,
                ev=[(1.6398, 1.9602)]),
    "16":  dict(weapon="M16", cal="5.56x45", az=175, el=0, D=304, b=24, T=24,
                ev=[(1.3699, 1.7982)]),
    "17":  dict(weapon="AK47 (5 single shots)", cal="7.62x39", az=28, el=0, D=88, b=4, T=26,
                ev=[(1.9229, 2.0457), (8.9673, 9.0910), (16.0529, 16.1774), (23.3789, 23.5051),
                    (30.9191, 31.0446)]),
}
SKIPPED = "07-14 (14.5 mm to 120 mm cannons and RPG-7): outside the development rule and the generator's bullet library"


# ---------------------------------------------------------------------------------------------
# measurement (the same estimators are applied to the real and the clean signal)
# ---------------------------------------------------------------------------------------------
def hpf(x, fs, fc, order=4):
    return sosfiltfilt(butter(order, fc, btype="high", fs=fs, output="sos"), x, axis=0)


def envelope(x, fs, lo, hi, smooth_s):
    y = us.bandpass(x, fs, lo, hi)
    return uniform_filter1d((y ** 2).sum(axis=1), max(1, int(smooth_s * fs)))


def find_event(x, fs, t_guess, kind, search_s):
    """Array-level time of the event's energy peak near t_guess."""
    lo, hi, sm = (1000.0, 8000.0, 0.0005) if kind == "nwave" else (150.0, 1000.0, 0.003)
    i0, i1 = max(0, int((t_guess - search_s) * fs)), min(len(x), int((t_guess + search_s) * fs))
    e = envelope(x[i0:i1], fs, lo, hi, sm)
    return (i0 + int(np.argmax(e))) / fs


def measure_event(x, fs, t_peak, kind):
    """Per channel: onset (s), peak (signal units) and duration (s).
    N-wave: 100 Hz high-pass, peak at the front shock, T_N front-to-rear shock.
    Blast: 50 Hz high-pass, positive peak, positive-phase duration."""
    nch = x.shape[1]
    out = dict(onset_s=[], peak=[], duration_s=[], polarity=[])
    xf = hpf(x, fs, 100.0 if kind == "nwave" else 50.0)
    for c in range(nch):
        on = us.channel_onset(x[:, c], fs, t_peak if kind == "nwave" else t_peak - 0.001,
                              kind, 0.0015 if kind == "nwave" else 0.0025)
        if kind == "nwave":
            w = us.cut(xf[:, c], fs, on - 0.00025, 0.0015)
        else:
            w = us.cut(xf[:, c], fs, on - 0.002, 0.020)
        if w is None:
            for k in out:
                out[k].append(np.nan)
            continue
        wu = resample_poly(w, 8, 1)
        fsu = 8 * fs
        pol = us.polarity(wu[: int((0.00025 + 0.0003) * fsu)] if kind == "nwave" else wu)
        wu = pol * wu
        if kind == "nwave":
            pk = float(np.max(wu[int(0.00015 * fsu): int(0.0004 * fsu)]))
            dur = us.nwave_duration(wu, fsu)
        else:
            pk = float(np.max(wu[int(0.0015 * fsu): int(0.008 * fsu)]))
            dur = us.blast_duration(wu, fsu)
        out["onset_s"].append(float(on))
        out["peak"].append(pk)
        out["duration_s"].append(float(dur))
        out["polarity"].append(int(pol))
    return out


# ---------------------------------------------------------------------------------------------
# synthesis
# ---------------------------------------------------------------------------------------------
def scene_vectors(sc, side):
    """Shooter position and bullet direction: the bullet is aimed at the array and turned
    horizontally by asin(b / D) to the given side (+1 / -1)."""
    az, el = np.radians(sc["az"]), np.radians(sc["el"])
    u = np.array([np.cos(el) * np.cos(az), -np.cos(el) * np.sin(az), np.sin(el)])
    shooter = sc["D"] * u
    d = np.arcsin(sc["b"] / sc["D"]) * side
    v = -u
    rot = np.array([[np.cos(d), -np.sin(d), 0.0], [np.sin(d), np.cos(d), 0.0], [0.0, 0.0, 1.0]])
    return shooter, rot @ v


def make_cfg(base, sc, side, comp, amp=1.0, dur=1.0):
    cfg = copy.deepcopy(base)
    shooter, v = scene_vectors(sc, side)
    cfg["signal"]["components"] = comp
    cfg["adc"]["fs"] = FS * UP
    cfg["bullet"].update(caliber=sc["cal"], mach=0.0, muzzle_velocity_mps=0.0, diameter_mm=0.0,
                         length_mm=0.0, drag_model="exponential", decay_length_m=0.0)
    cfg["trajectory"].update(mode="explicit", shooter_pos=shooter.tolist(), direction=v.tolist())
    cfg["atmosphere"].update(enabled=True, temp_c=float(sc["T"]))
    cfg["propagation"]["wind_mps"] = [0.0, 0.0, 0.0]
    if "ground" in cfg:                     # clean twins: no ground reflection
        cfg["ground"]["enabled"] = False
    key = "shockwave" if comp == "shockwave" else "muzzle_blast"
    cfg[key]["amp_scale"], cfg[key]["dur_scale"] = float(amp), float(dur)
    cfg["timing"].update(pre_roll=0.010, post_roll=0.150 if comp == "muzzle_blast" else 0.010)
    return cfg


def synth(cfg):
    """(signal in Pa at FS x UP [n, ch], time of its first sample in s, ground truth)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        sig, fs, gt = generate_clean(cfg)
    n = min(len(s) for s in sig)
    return np.stack([s[:n] for s in sig], axis=1), gt["t_master_start_s"], gt


def decimate_block(y_hi, t0_hi, shift):
    """Band-limit and decimate a high-rate block that is to start at t0_hi + shift.
    Returns (block at FS, index of its first sample at FS)."""
    start_hi = int(round((t0_hi + shift) * FS * UP))
    pad = start_hi % UP
    y = np.vstack([np.zeros((pad + 20 * UP, y_hi.shape[1])), y_hi, np.zeros((20 * UP, y_hi.shape[1]))])
    yd = resample_poly(y, 1, UP, axis=0)
    return yd, (start_hi - pad) // UP - 20


def measure_block(yd, i0, kind, t_peak):
    """Measure a decimated block with the same estimators as the recording (times in file time)."""
    pre = int(0.05 * FS)
    buf = np.vstack([np.zeros((pre, yd.shape[1])), yd, np.zeros((pre, yd.shape[1]))])
    m = measure_event(buf, FS, t_peak - (i0 - pre) / FS, kind)
    m["onset_s"] = [o + (i0 - pre) / FS for o in m["onset_s"]]
    return m


def fit_component(base, sc, side, comp, meas, n_iter=3, fit=True):
    """Scale amplitude and duration to the measured means and find the array time shift that
    puts the mean onset at the measured mean onset.
    Returns (high-rate block in FS units, its start time, shift, amp, dur, ground truth, model time)."""
    kind = "nwave" if comp == "shockwave" else "blast"
    key = "shockwave" if comp == "shockwave" else "muzzle_blast"
    amp = dur = 1.0
    t_meas = np.nanmean(meas["onset_s"])
    for it in range(n_iter + 1):
        y, t0, gt = synth(make_cfg(base, sc, side, comp, amp, dur))
        t_model = float(np.nanmean([pm[key]["t_arrive_s"] for pm in gt["per_mic"]]))
        yd, i0 = decimate_block(y / FS_PA, t0, 0.0)
        g = measure_block(yd, i0, kind, t_model + (0.0 if kind == "nwave" else 0.002))
        shift = float(t_meas - np.nanmean(g["onset_s"]))
        if not fit or it == n_iter:
            break
        amp *= np.nanmean(meas["peak"]) / np.nanmean(g["peak"])
        dur *= np.nanmedian(meas["duration_s"]) / np.nanmedian(g["duration_s"])
    return y / FS_PA, t0, shift, float(amp), float(dur), gt, t_model


def channel_lags(x, y_hi, t0, shift, t_ev, kind, max_lag_s, common=False):
    """Per-channel residual shift (s) that maximizes the cross-correlation between the recording
    and the array-placed clean event (high-passed, upsampled x 8, lag limited to +-max_lag_s)."""
    a, b, fc = (0.0003, 0.0006, 100.0) if kind == "nwave" else (0.003, 0.012, 50.0)
    i0, i1 = max(0, int((t_ev - 0.06) * FS)), min(len(x), int((t_ev + 0.06) * FS))
    seg = np.zeros((i1 - i0, x.shape[1]))
    yd, j0 = decimate_block(y_hi, t0, shift)
    add(seg, yd, j0 - i0)
    xr, xc = hpf(x[i0:i1], FS, fc), hpf(seg, FS, fc)
    k0, k1 = int((t_ev - a - 0.0015) * FS) - i0, int((t_ev + b + 0.0015) * FS) - i0   # same window for every lag range
    up, m = 8, int(max_lag_s * FS * 8)
    k0, k1 = max(0, k0), min(len(xr), k1)
    ws = []
    for c in range(x.shape[1]):
        r = resample_poly(xr[k0:k1, c], up, 1)
        q = resample_poly(xc[k0:k1, c], up, 1)
        cc = np.correlate(r, q, mode="full")
        mid = len(q) - 1
        ws.append(cc[mid - m: mid + m + 1] / (np.linalg.norm(r) * np.linalg.norm(q) + 1e-30))
    if common:                                   # one lag for all channels
        return (int(np.argmax(np.sum(ws, axis=0))) - m) / (FS * up)
    # a channel whose best lag is at the edge of the range has no clear optimum: keep the common shift
    return [0.0 if int(np.argmax(w)) in (0, 2 * m) else (int(np.argmax(w)) - m) / (FS * up) for w in ws]


def place(buf, y_hi, t0, shifts):
    """Add the high-rate block to buf, channel c shifted by shifts[c] (s)."""
    for c, sh in enumerate(shifts):
        yd, i0 = decimate_block(y_hi[:, c:c + 1], t0, sh)
        add(buf[:, c:c + 1], yd, i0)


def tdoa_us(onsets):
    o = np.asarray(onsets, dtype=float)
    return ((o[1:] - o[0]) * 1e6).tolist()


def corr(a, b):
    a, b = a - a.mean(), b - b.mean()
    d = np.sqrt((a * a).sum() * (b * b).sum())
    return float((a * b).sum() / d) if d > 0 else float("nan")


# ---------------------------------------------------------------------------------------------
# per recording
# ---------------------------------------------------------------------------------------------
def process(shot, base, in_dir, out_dir, fit=True, cache=None, align="channel"):
    sc = SCENES[shot]
    wav_in = os.path.join(in_dir, f"shot-{shot}_audio_4ch_cut.wav")
    fs, x = us.load_wav(wav_in)
    if fs != FS:
        raise ValueError(f"{wav_in}: {fs} Hz, expected {FS}")
    n, nch = x.shape
    clean = np.zeros((n, nch))
    events, side, side_rms = [], None, {}

    def component(comp, meas, t_peak, amp_dur=None):
        kind = "nwave" if comp == "shockwave" else "blast"
        key = comp
        if amp_dur is None:
            y, t0, sh, amp, dur, gt, t_mod = fit_component(base, sc, side, comp, meas, fit=fit)
            t_ev = float(np.nanmean(meas["onset_s"]))
            common = channel_lags(x, y, t0, sh, t_ev, kind, 0.0003 if kind == "nwave" else 0.0015, common=True)
            sh += common
            lags = (channel_lags(x, y, t0, sh, t_ev, kind, 0.0002) if align == "channel"
                    else [0.0] * nch)
        else:                                   # no measurement: model timing relative to the SW
            amp, dur, sh = amp_dur
            y, t0, gt = synth(make_cfg(base, sc, side, comp, amp, dur))
            y = y / FS_PA
            t_mod = float(np.mean([pm[key]["t_arrive_s"] for pm in gt["per_mic"]]))
            lags = [0.0] * nch
        place(clean, y, t0, [sh + l for l in lags])
        t_model = [pm[key]["t_arrive_s"] for pm in gt["per_mic"]]
        rec = dict(amp_scale=amp, dur_scale=dur, time_shift_s=sh,
                   channel_shift_us=[1e6 * l for l in lags],
                   t_arrive_s=[t + sh + l for t, l in zip(t_model, lags)],
                   tdoa_model_us=tdoa_us(t_model),
                   tdoa_clean_us=tdoa_us([t + l for t, l in zip(t_model, lags)]))
        if comp == "shockwave":
            rec["model"] = dict(  # generator values without the fitted scales
                peak_pa=[pm[key]["peak_dP_pa"] / amp for pm in gt["per_mic"]],
                T_N_s=[pm[key]["duration_T_N_s"] / dur for pm in gt["per_mic"]],
                miss_distance_b_m=[pm[key]["perp_dist_b_m"] for pm in gt["per_mic"]],
                mach_emit=[pm[key].get("mach_emit") for pm in gt["per_mic"]])
        else:
            rec["model"] = dict(peak_pa=[pm[key]["peak_dP_pa"] / amp for pm in gt["per_mic"]],
                                t_pos_s=[pm[key]["positive_phase_t_pos_s"] / dur for pm in gt["per_mic"]],
                                range_m=[pm[key]["range_m"] for pm in gt["per_mic"]])
        if meas is not None:
            rec["measured"] = dict(t_s=float(np.nanmean(meas["onset_s"])), onset_s=meas["onset_s"],
                                   peak_pa=[p * FS_PA for p in meas["peak"]],
                                   **({"T_N_s": meas["duration_s"]} if kind == "nwave" else {"t_pos_s": meas["duration_s"]}),
                                   polarity=meas["polarity"], tdoa_us=tdoa_us(meas["onset_s"]))
            rec["_t_peak"] = t_peak
        else:
            rec["measured"] = None
        return rec, gt, t_mod

    for k, (tsw_g, tmb_g) in enumerate(sc["ev"]):
        ev = dict(index=k + 1)
        tsw = find_event(x, fs, tsw_g, "nwave", 0.01)
        msw = measure_event(x, fs, tsw, "nwave")
        if side is None:                                   # choose the bullet side on event 1
            for s in (+1, -1):
                _, _, gt = synth(make_cfg(base, sc, s, "shockwave"))
                t_mod = [pm["shockwave"]["t_arrive_s"] for pm in gt["per_mic"]]
                e = np.array(tdoa_us(t_mod)) - np.array(tdoa_us(msw["onset_s"]))
                side_rms[s] = float(np.sqrt(np.nanmean(e ** 2)))
            side = min(side_rms, key=side_rms.get)
        ev["shockwave"], gt_sw, tmod_sw = component("shockwave", msw, tsw)
        if tmb_g is not None:
            tmb = find_event(x, fs, tmb_g, "blast", 0.015)
            mmb = measure_event(x, fs, tmb, "blast")
            ev["muzzle_blast"], _, tmod_mb = component("muzzle_blast", mmb, tmb)
            ev["sw_to_mb_delay_s"] = dict(
                measured=ev["muzzle_blast"]["measured"]["t_s"] - ev["shockwave"]["measured"]["t_s"],
                model=tmod_mb - tmod_sw)
        else:   # no measurable blast (the cut ends before it): model delay, levels of the reference shot
            ref = (cache or {}).get(sc.get("mb_from"))
            a_mb, d_mb = (ref["amp"], ref["dur"]) if ref else (1.0, 1.0)
            mb, _, tmod_mb = component("muzzle_blast", None, None, (a_mb, d_mb, ev["shockwave"]["time_shift_s"]))
            mb["levels_from"] = f"shot-{sc['mb_from']}" if ref else "model"
            mb["t_s"] = float(np.mean(mb["t_arrive_s"]))
            mb["inside_file"] = bool(mb["t_s"] < n / FS)
            ev["muzzle_blast"] = mb
        events.append(ev)

    # ---- the clean signal measured like the recording, and similarity (aligned windows)
    xf_sw, cf_sw = hpf(x, FS, 100.0), hpf(clean, FS, 100.0)
    xf_mb, cf_mb = hpf(x, FS, 50.0), hpf(clean, FS, 50.0)
    for ev in events:
        for comp, kind, xf, cf, (a_, w_) in (("shockwave", "nwave", xf_sw, cf_sw, (0.00025, 0.0010)),
                                             ("muzzle_blast", "blast", xf_mb, cf_mb, (0.002, 0.012))):
            r = ev[comp]
            if not r.get("measured"):
                continue
            g = measure_event(clean, FS, r.pop("_t_peak"), kind)
            r["clean"] = dict(onset_s=g["onset_s"], peak_pa=[p * FS_PA for p in g["peak"]],
                              **({"T_N_s": g["duration_s"]} if kind == "nwave" else {"t_pos_s": g["duration_s"]}))
            r["pearson_r"] = [corr(us.cut(xf[:, c], FS, o - a_, w_), us.cut(cf[:, c], FS, o - a_, w_))
                              for c, o in enumerate(r["measured"]["onset_s"])]
            if kind == "nwave":   # onset picking is not reliable enough on the real muzzle blasts
                r["tdoa_err_us"] = [float(v) for v in (np.array(tdoa_us(g["onset_s"]))
                                                       - np.array(tdoa_us(r["measured"]["onset_s"])))]

    # ---- write
    peak = float(np.abs(clean).max())
    if peak >= 1.0:
        raise ValueError(f"shot-{shot}: clean signal clips ({peak:.3f} FS)")
    import soundfile as sf
    base_name = f"shot-{shot}_audio_4ch_clean"
    sf.write(os.path.join(out_dir, base_name + ".wav"), clean.astype(np.float32), FS, subtype="PCM_24")
    shooter, v = scene_vectors(sc, side)
    meta = dict(
        generated_at_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        generator="util_gen_cipatat_clean.py (gs_gen_clean_signal.generate_clean)",
        source_recording=os.path.basename(wav_in), wav_file=base_name + ".wav",
        format=dict(channels=nch, sample_rate_hz=FS, subtype="PCM_24", samples=n,
                    value_scale=f"full scale = {FS_PA:g} Pa (estimate)", synthesis_rate_hz=FS * UP),
        contents="shockwave + muzzle blast only: no noise, no reverberation, no ground reflection, no echoes",
        fitted_to_recording=bool(fit),
        alignment=("per channel: each event at the recording's arrival on every channel (TDOAs of the recording)"
                   if align == "channel" else "array: each event at the recording's mean arrival, TDOAs of the model scene"),
        scene=dict(weapon=sc["weapon"], calibre=sc["cal"], azimuth_deg=sc["az"], elevation_deg=sc["el"],
                   distance_m=sc["D"], miss_distance_m=sc["b"], air_temp_c=sc["T"],
                   speed_of_sound_mps=gp.speed_of_sound(sc["T"]),
                   bullet=dict(gp.BULLET_LIBRARY[gp.BULLET_ALIASES.get(sc["cal"], sc["cal"])]),
                   drag_model="exponential", humidity_pct=base["atmosphere"]["humidity_pct"],
                   bullet_side=("+" if side > 0 else "-"), side_tdoa_rms_us={("+" if s > 0 else "-"): v_ for s, v_ in side_rms.items()},
                   shooter_pos_m=shooter.tolist(), trajectory_direction=v.tolist(),
                   mic_positions_m=gp.build_mic_positions(make_cfg(base, sc, side, "shockwave")).tolist(),
                   frame="x = north, y = west, z = up; azimuth clockwise from north"),
        events=events,
        ground_truth_shockwave_model=dict((k_, gt_sw[k_]) for k_ in
                                          ("bullet_speed_mps", "mach", "decay_length_m", "tdoa_matrix_sw_us", "per_mic")),
        peak_fs=peak)
    with open(os.path.join(out_dir, base_name + ".json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2, default=_jd)
    return meta, x, clean


def add(buf, yd, i0):
    i1 = i0 + len(yd)
    a, b = max(0, i0), min(len(buf), i1)
    if b > a:
        buf[a:b] += yd[a - i0:b - i0]


def _jd(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, np.bool_):
        return bool(o)
    raise TypeError(type(o))


# ---------------------------------------------------------------------------------------------
def summary_rows(shot, meta):
    rows = []
    for ev in meta["events"]:
        sw, mb = ev["shockwave"], ev["muzzle_blast"]
        m = mb.get("measured")
        rows.append(dict(
            shot=shot, event=ev["index"], calibre=meta["scene"]["calibre"], side=meta["scene"]["bullet_side"],
            t_sw_s=round(sw["measured"]["t_s"], 4),
            sw_peak_pa=round(float(np.nanmean(sw["measured"]["peak_pa"])), 1),
            sw_model_peak_pa=round(float(np.nanmean([p for p in sw["model"]["peak_pa"] if p])), 1),
            sw_T_N_us=round(1e6 * float(np.nanmedian(sw["measured"]["T_N_s"]))),
            sw_model_T_N_us=round(1e6 * float(np.nanmean([t for t in sw["model"]["T_N_s"] if t]))),
            sw_amp_scale=round(sw["amp_scale"], 3), sw_dur_scale=round(sw["dur_scale"], 3),
            sw_r=round(float(np.nanmean(sw["pearson_r"])), 2),
            sw_tdoa_rms_us=round(float(np.sqrt(np.nanmean(np.square(sw["tdoa_err_us"])))), 1),
            t_mb_s=round(m["t_s"], 4) if m else round(mb["t_s"], 4),
            mb_peak_pa=round(float(np.nanmean(m["peak_pa"])), 2) if m else "",
            mb_model_peak_pa=round(float(np.nanmean(mb["model"]["peak_pa"])), 1) if m else "",
            mb_t_pos_ms=round(1e3 * float(np.nanmedian(m["t_pos_s"])), 2) if m else "",
            mb_model_t_pos_ms=round(1e3 * float(np.nanmean(mb["model"]["t_pos_s"])), 2) if m else "",
            mb_amp_scale=round(mb["amp_scale"], 4), mb_dur_scale=round(mb["dur_scale"], 3),
            mb_r=round(float(np.nanmean(mb["pearson_r"])), 2) if m else "",
            delay_meas_ms=round(1e3 * ev["sw_to_mb_delay_s"]["measured"], 1) if m else "",
            delay_model_ms=round(1e3 * ev["sw_to_mb_delay_s"]["model"], 1) if m else "",
            note="" if m else ("MB after the end of the file" if not mb["inside_file"] else "MB at model delay")))
    return rows


def plot_overview(results, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nr = len(results)
    fig, ax = plt.subplots(nr, 2, figsize=(12, 1.9 * nr), squeeze=False)
    for r, (shot, meta, x, clean) in enumerate(results):
        ev = meta["events"][0]
        t_sw = ev["shockwave"]["measured"]["onset_s"][0]
        xs, cs = hpf(x, FS, 100.0), hpf(clean, FS, 100.0)
        i0 = int((t_sw - 0.0005) * FS); i1 = int((t_sw + 0.0035) * FS)
        t = (np.arange(i0, i1) / FS - t_sw) * 1e3
        a = ax[r, 0]
        a.plot(t, xs[i0:i1, 0] * FS_PA, color="0.6", lw=1, label="recording (HP 100 Hz)")
        a.plot(t, cs[i0:i1, 0] * FS_PA, color="C0", lw=1.2, label="clean")
        a.set_ylabel(f"shot-{shot}\nPa", fontsize=8)
        a.tick_params(labelsize=7)
        mb = ev["muzzle_blast"]
        a = ax[r, 1]
        if mb.get("measured"):
            t_mb = mb["measured"]["onset_s"][0]
            xm, cm = hpf(x, FS, 50.0), hpf(clean, FS, 50.0)
            i0 = int((t_mb - 0.005) * FS); i1 = int((t_mb + 0.025) * FS)
            t = (np.arange(i0, i1) / FS - t_mb) * 1e3
            a.plot(t, xm[i0:i1, 0] * FS_PA, color="0.6", lw=1, label="recording (HP 50 Hz)")
            a.plot(t, cm[i0:i1, 0] * FS_PA, color="C3", lw=1.2, label="clean")
        else:
            a.text(0.5, 0.5, "no muzzle blast in this file", ha="center", va="center", transform=a.transAxes, fontsize=8)
        a.tick_params(labelsize=7)
        if r == 0:
            ax[r, 0].set_title("Shockwave, channel 1 (ms from onset)", fontsize=9)
            ax[r, 1].set_title("Muzzle blast, channel 1 (ms from onset)", fontsize=9)
            ax[r, 0].legend(fontsize=7, loc="upper right")
            ax[r, 1].legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def main():
    global FS_PA
    ap = argparse.ArgumentParser(description="Clean twins of the Cipatat development recordings.")
    ap.add_argument("--in", dest="in_dir", default=os.path.join(HERE, "..", "sound", "Cipatat"))
    ap.add_argument("--out", dest="out_dir", default=os.path.join(HERE, "..", "sound", "Cipatat_clean"))
    ap.add_argument("--config", default=os.path.join(HERE, "gen_config.ini"))
    ap.add_argument("--shots", nargs="*", default=list(SCENES))
    ap.add_argument("--no-fit", action="store_true", help="keep the model's levels and durations")
    ap.add_argument("--fs-pa", type=float, default=FS_PA, help="Pa per full scale (default 1520)")
    ap.add_argument("--align", choices=["channel", "array"], default="channel",
                    help="channel (default): every channel at the recording's arrival time; "
                         "array: mean arrival only, model TDOAs")
    ap.add_argument("--no-plot", action="store_true")
    a = ap.parse_args()
    FS_PA = a.fs_pa
    os.makedirs(a.out_dir, exist_ok=True)
    base = gp.load_config(a.config, "PILAR-V", "both")
    shots = [s for s in SCENES if s in a.shots]
    if "06a" in shots and "06" not in shots:
        shots.insert(shots.index("06a"), "06")
    cache, results, rows = {}, [], []
    for s in shots:
        meta, x, clean = process(s, base, a.in_dir, a.out_dir, fit=not a.no_fit, cache=cache, align=a.align)
        mb = meta["events"][0]["muzzle_blast"]
        cache[s] = dict(amp=mb["amp_scale"], dur=mb["dur_scale"])
        results.append((s, meta, x, clean))
        rows += summary_rows(s, meta)
        r0 = rows[-len(meta["events"])]
        print(f"shot-{s}: side {meta['scene']['bullet_side']}, SW r {r0['sw_r']}, TDOA rms {r0['sw_tdoa_rms_us']} us, "
              f"MB r {r0['mb_r']}, peak {meta['peak_fs']:.3f} FS")
    with open(os.path.join(a.out_dir, "clean_summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    if not a.no_plot:
        plot_overview(results, os.path.join(a.out_dir, "clean_vs_real_overview.png"))
    print(f"Skipped: {SKIPPED}")
    print(f"Output folder: {os.path.abspath(a.out_dir)}")


if __name__ == "__main__":
    main()
