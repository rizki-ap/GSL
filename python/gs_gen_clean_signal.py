#!/usr/bin/env python3
"""
gs_gen_clean_signal.py
======================
STAGE 1. Noiseless multichannel gunshot pressure signal: ballistic shockwave (N-wave) plus muzzle
blast (Friedlander pulse), physics only. Writes <o>_clean.wav (float32, Pa) and <o>_clean.json
(configuration and exact ground truth: arrival times, TDOAs, geometry, checks).

Usage
-----
  python gs_gen_clean_signal.py gen_config.ini
  python gs_gen_clean_signal.py gen_config.ini -o my_shot

Next stages: gs_gen_add_noise.py, then gs_gen_apply_adc.py (both assume 4 channels; for arrays
with another number of sensors use the clean output directly).

Model summary (see gs_gen_physic.py): Mach-cone shock arrival in closed form, Whitham N-wave,
Friedlander blast with Hopkinson-Cranz/Kinney-Graham scaling, speed of sound from temperature,
optional first-order wind, causal ISO 9613-1 absorption per event, no shock below
shockwave.mach_threshold and at sensors the Mach cone does not reach.
"""

import argparse
import json
from datetime import datetime, timezone

import numpy as np
from scipy.io import wavfile

import gs_gen_physic as gp


def _none_if_nan(x):
    return None if (x is None or (isinstance(x, float) and np.isnan(x))) else float(x)


def _matrix_us(times):
    """N x N TDOA matrix in microseconds (row i minus column j); None where an event is absent."""
    t = np.asarray(times, dtype=float)
    m = (t[:, None] - t[None, :]) * 1e6
    return [[_none_if_nan(v) for v in row] for row in m]


def generate_clean(cfg):
    """Returns (list of per-sensor signals in Pa, sample rate, ground-truth dict)."""
    prop, shk, mbc = cfg["propagation"], cfg["shockwave"], cfg["muzzle_blast"]
    tim, atm_cfg = cfg["timing"], cfg["atmosphere"]
    FS = cfg["adc"]["fs"]

    atmosphere = None
    if atm_cfg["enabled"]:
        atmosphere = dict(temp_c=atm_cfg["temp_c"], humidity_pct=atm_cfg["humidity_pct"],
                          pressure_kpa=atm_cfg["pressure_kpa"])
    p_atm = atm_cfg["pressure_kpa"] * 1e3

    c = gp.speed_of_sound(atm_cfg["temp_c"])
    wind = np.asarray(prop["wind_mps"], dtype=float)
    wind = wind if np.linalg.norm(wind) > 0 else None

    mic_pos = gp.build_mic_positions(cfg)
    n_mics = len(mic_pos)
    bullet, M = gp.resolve_bullet(cfg, c)
    origin, v_hat = gp.build_trajectory(cfg)
    shooter = origin.copy()
    shock_on = M > max(shk["mach_threshold"], 1.0)

    # ---- arrival times ------------------------------------------------------------------
    geo = [gp.shock_geometry(p, origin, v_hat, M, c, wind) if shock_on else None for p in mic_pos]
    t_mb = np.zeros(n_mics)
    r_mb = np.zeros(n_mics)
    for i, p in enumerate(mic_pos):
        t_mb[i], r_mb[i], _ = gp.muzzle_arrival(p, shooter, c, wind)
    t_sw = np.array([g["t_arrive"] if (g is not None and g["in_cone"]) else np.nan for g in geo]) \
        if shock_on else np.full(n_mics, np.nan)

    # ---- master time axis -------------------------------------------------------------------
    t0 = np.nanmin(np.concatenate([t_sw, t_mb])) - tim["pre_roll"]
    t1 = t_mb.max() + tim["post_roll"]
    t_master = np.arange(t0, t1, 1.0 / FS)

    # ---- waveforms -----------------------------------------------------------------------------
    w_kg = mbc["w_g_tnt"] * 1e-3
    signals = []
    T_N = np.full(n_mics, np.nan)
    dP_sw = np.full(n_mics, np.nan)
    t_pos = np.zeros(n_mics)
    dP_mb = np.zeros(n_mics)
    att_sw = np.zeros(n_mics)
    att_mb = np.zeros(n_mics)

    for i in range(n_mics):
        sw = np.zeros(len(t_master))
        if shock_on and geo[i]["in_cone"]:
            g = geo[i]
            dP, T = gp.nwave_params(g["b"], M, c, bullet, shk["amp_scale"], shk["dur_scale"], p_atm)
            sw = gp.nwave_waveform(t_master, g["t_arrive"], dP, T)
            if atmosphere is not None:
                sw = gp.spectral_absorption(sw, FS, g["R"], atmosphere)
                att_sw[i] = float(gp.atmospheric_attenuation_db(1.0 / T, g["R"], **atmosphere))
            dP_sw[i], T_N[i] = dP, T

        dP, tp = gp.blast_params(r_mb[i], w_kg, mbc["amp_scale"], mbc["dur_scale"], p_atm)
        mb = gp.friedlander_waveform(t_master, t_mb[i], dP, tp)
        if atmosphere is not None:
            mb = gp.spectral_absorption(mb, FS, r_mb[i], atmosphere)
            att_mb[i] = float(gp.atmospheric_attenuation_db(1.0 / tp, r_mb[i], **atmosphere))
        dP_mb[i], t_pos[i] = dP, tp
        signals.append(sw + mb)

    # ---- ground truth ---------------------------------------------------------------------------
    per_mic = []
    for i in range(n_mics):
        g = geo[i]
        reason = None
        if not shock_on:
            reason = "below_mach_threshold"
        elif not g["in_cone"]:
            reason = "outside_mach_cone"
        sw_rec = dict(active=reason is None, inactive_reason=reason)
        if g is not None:
            sw_rec.update(along_track_a_m=g["a"], perp_dist_b_m=g["b"], emission_coord_s_m=g["s"])
        if reason is None:
            sw_rec.update(t_arrive_s=float(t_sw[i]), t_emit_s=float(g["t_emit"]),
                          acoustic_path_R_m=float(g["R"]), duration_T_N_s=float(T_N[i]),
                          peak_dP_pa=float(dP_sw[i]), atmospheric_attenuation_db=float(att_sw[i]))
        per_mic.append(dict(
            mic_index=i, shockwave=sw_rec,
            muzzle_blast=dict(t_arrive_s=float(t_mb[i]), range_m=float(r_mb[i]),
                              positive_phase_t_pos_s=float(t_pos[i]), peak_dP_pa=float(dP_mb[i]),
                              atmospheric_attenuation_db=float(att_mb[i]))))

    active = ~np.isnan(t_sw)
    sep = (t_mb - t_sw)[active] if active.any() else np.array([])
    checks = dict(
        shock_before_blast_all_active=bool(np.all(sep > 0)) if sep.size else None,
        min_blast_minus_shock_s=float(sep.min()) if sep.size else None,
        events_overlap_any=bool(np.any(sep < T_N[active])) if sep.size else None,
        n_sensors_with_shock=int(active.sum()), n_sensors=int(n_mics))

    ground_truth = dict(
        mic_pos=mic_pos.tolist(), shooter_pos=shooter.tolist(), bullet_origin=origin.tolist(),
        trajectory_direction=v_hat.tolist(), bullet_speed_mps=float(M * c), mach=float(M),
        speed_of_sound_mps=float(c), wind_mps=(wind.tolist() if wind is not None else [0.0, 0.0, 0.0]),
        sample_rate_hz=FS, atmosphere_applied=atmosphere,
        models=dict(shockwave="whitham", muzzle_blast="hopkinson_cranz_kinney_graham",
                    absorption="iso_9613_1_minimum_phase", shock_branch_enabled=bool(shock_on),
                    mach_threshold=shk["mach_threshold"]),
        t_master_start_s=float(t0), t_master_end_s=float(t1), n_samples=len(t_master),
        per_mic=per_mic,
        true_tdoa_sw_us=[None if (np.isnan(t_sw[i]) or np.isnan(t_sw[0])) else float((t_sw[i] - t_sw[0]) * 1e6)
                         for i in range(1, n_mics)],
        true_tdoa_mb_us=[float((t_mb[i] - t_mb[0]) * 1e6) for i in range(1, n_mics)],
        tdoa_matrix_sw_us=_matrix_us(t_sw), tdoa_matrix_mb_us=_matrix_us(t_mb),
        dt_mb_minus_sw_centroid_s=(float(np.nanmean(t_mb[active]) - np.nanmean(t_sw[active]))
                                   if active.any() else None),
        checks=checks)
    return signals, FS, ground_truth


def write_clean_wav(path, signals, fs):
    n = min(len(s) for s in signals)
    wavfile.write(path, fs, np.stack([s[:n] for s in signals], axis=1).astype(np.float32))


def _json_default(o):
    if isinstance(o, np.floating):
        return float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not JSON serialisable: {type(o)}")


def main():
    ap = argparse.ArgumentParser(description="Stage 1: noiseless multichannel gunshot signal.")
    ap.add_argument("config", nargs="?", default="gen_config.ini")
    ap.add_argument("-o", "--output", default=None, help="output basename (overrides output.basename)")
    args = ap.parse_args()

    cfg = gp.load_config(args.config)
    if args.output:
        cfg["output"]["basename"] = args.output
    signals, fs, gt = generate_clean(cfg)

    base = cfg["output"]["basename"]
    wav_path, json_path = f"{base}_clean.wav", f"{base}_clean.json"
    write_clean_wav(wav_path, signals, fs)
    with open(json_path, "w") as f:
        json.dump(dict(generated_at_utc=datetime.now(timezone.utc).isoformat(), stage="clean_signal",
                       generator="gs_gen_clean_signal.py", wav_file=wav_path,
                       wav_format="float32 PCM, unnormalized Pa units", sample_rate_hz=fs,
                       n_channels=len(signals), n_samples=len(signals[0]),
                       duration_s=len(signals[0]) / fs, config_used=cfg, ground_truth=gt),
                  f, indent=2, default=_json_default)

    ck = gt["checks"]
    print(f"[Stage 1] Wrote {wav_path} ({len(signals)} ch, {fs} Hz, {len(signals[0]) / fs:.3f} s) and {json_path}")
    print(f"[Stage 1] shock at {ck['n_sensors_with_shock']}/{ck['n_sensors']} sensors; "
          f"min(blast - shock) = {ck['min_blast_minus_shock_s']} s; events overlap = {ck['events_overlap_any']}")
    if len(signals) != 4:
        print("[Stage 1] NOTE: the noise and ADC stages assume 4 channels.")
    print(f"[Stage 1] Next: python gs_gen_add_noise.py {wav_path} {args.config}")


if __name__ == "__main__":
    main()
