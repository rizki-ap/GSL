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
  python gs_gen_clean_signal.py gen_config.ini -mic_config=UMA8     # miniDSP UMA-8 v2 (4 of its mics)
  python gs_gen_clean_signal.py gen_config.ini -components=shockwave     # shockwave only
  python gs_gen_clean_signal.py gen_config.ini -components=muzzle_blast  # muzzle blast only

-mic_config selects the microphone array ([mic_config.<NAME>] section of the .ini; default
[geometry] mic_config = Tetrahedral-0,3, the regular 4-mic tetrahedron with 0.30 m edge;
PILAR-V = the real sensor geometry). One channel per mic.

-components selects the acoustic events to synthesise: both (default, [signal] components in
the .ini), shockwave, or muzzle_blast (aliases sw, mb, blast). With one event the time axis
spans only that event's arrivals (pre_roll before the first, post_roll after the last). The
ground truth still describes both events; per_mic[i].shockwave.generated and
per_mic[i].muzzle_blast.generated say which ones are in the signal, and 'components' records
the choice. Asking for the shockwave when no microphone receives one (subsonic bullet, or
outside the Mach cone) is an error rather than a silent file.

Next stages: gs_gen_add_noise.py, then gs_gen_apply_adc.py (any number of channels; they reuse
the mic_config recorded in this stage's .json).

Model summary (see gs_gen_physic.py): Mach-cone shock arrival in closed form, Whitham N-wave,
Friedlander blast with Hopkinson-Cranz/Kinney-Graham scaling, speed of sound from temperature,
optional first-order wind, causal ISO 9613-1 absorption per event, no shock below
shockwave.mach_threshold and at sensors the Mach cone does not reach.

Ground reflection ([ground] enabled = true): flat ground at z = -array_height_m in the array
frame. Each sensor also receives the shockwave and the blast at its mirror image below the
ground (image-source method), scaled by reflection_coeff_sw / reflection_coeff_mb. The
reflected N-wave has its own arrival time, miss distance, amplitude and duration; it reaches
the higher sensors later. Its ground truth is in per_mic[i].shockwave_reflected and
per_mic[i].muzzle_blast_reflected. It follows -components (a reflected shockwave only with the
shockwave). A reflection is omitted where the emission point (shock) or the shooter (blast) is
not above the ground plane.
"""

import argparse
import sys
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


def _reflection_point(x_emit, img_pos, z_g):
    """Point where the straight ray from the emission point to the image sensor crosses the ground."""
    x_emit, img_pos = np.asarray(x_emit, float), np.asarray(img_pos, float)
    f = (x_emit[2] - z_g) / (x_emit[2] - img_pos[2])
    return (x_emit + f * (img_pos - x_emit)).tolist()


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
    drag = bullet.get("drag", "none") == "exponential"
    comps = cfg["signal"]["components"]
    gen_sw = comps in ("both", "shockwave")
    gen_mb = comps in ("both", "muzzle_blast")

    # ---- arrival times ------------------------------------------------------------------
    if drag:
        geo = [gp.shock_geometry_drag(p, origin, v_hat, M * c, bullet["lam"], c, wind,
                                      shk["mach_threshold"]) if shock_on else None for p in mic_pos]
    else:
        geo = [gp.shock_geometry(p, origin, v_hat, M, c, wind) if shock_on else None for p in mic_pos]
    t_mb = np.zeros(n_mics)
    r_mb = np.zeros(n_mics)
    for i, p in enumerate(mic_pos):
        t_mb[i], r_mb[i], _ = gp.muzzle_arrival(p, shooter, c, wind)
    t_sw = np.array([g["t_arrive"] if (g is not None and g["in_cone"]) else np.nan for g in geo]) \
        if shock_on else np.full(n_mics, np.nan)

    # ---- ground reflection: the same solutions at the image sensors ---------------------------
    gnd = cfg.get("ground", {"enabled": False})
    ground_on = bool(gnd["enabled"])
    t_swr = np.full(n_mics, np.nan)
    t_mbr = np.full(n_mics, np.nan)
    r_mbr = np.full(n_mics, np.nan)
    geo_r = [None] * n_mics
    swr_reason = ["ground_disabled"] * n_mics
    mbr_reason = ["ground_disabled"] * n_mics
    if ground_on:
        z_g = -gnd["array_height_m"]
        if np.min(mic_pos[:, 2]) <= z_g:
            raise ValueError(f"ground.array_height_m = {gnd['array_height_m']} m puts the ground above a "
                             f"microphone (lowest mic at z = {np.min(mic_pos[:, 2]):.3f} m)")
        img = gp.image_positions(mic_pos, z_g)
        for i in range(n_mics):
            if shock_on:
                if drag:
                    g = gp.shock_geometry_drag(img[i], origin, v_hat, M * c, bullet["lam"], c, wind,
                                               shk["mach_threshold"])
                else:
                    g = gp.shock_geometry(img[i], origin, v_hat, M, c, wind)
                geo_r[i] = g
                if not g["in_cone"]:
                    swr_reason[i] = "subsonic_at_emission" if g.get("subsonic") else "outside_mach_cone"
                elif g["x_emit"][2] <= z_g:
                    swr_reason[i] = "emission_below_ground"
                else:
                    swr_reason[i] = None
                    t_swr[i] = g["t_arrive"]
            else:
                swr_reason[i] = "below_mach_threshold"
            if shooter[2] > z_g:
                t_mbr[i], r_mbr[i], _ = gp.muzzle_arrival(img[i], shooter, c, wind)
                mbr_reason[i] = None
            else:
                mbr_reason[i] = "shooter_below_ground"

    # ---- master time axis: spans the events being generated --------------------------------
    if comps == "shockwave" and np.all(np.isnan(t_sw)):
        if not shock_on:
            why = f"the bullet is below the Mach threshold (Mach {M:.2f} <= {max(shk['mach_threshold'], 1.0):.2f})"
        elif drag and any(g is not None and g["subsonic"] for g in geo):
            why = "the bullet has slowed below the Mach threshold before its shock reaches the array"
        else:
            why = "no microphone is inside the Mach cone"
        raise ValueError(f"components = shockwave, but no microphone receives a shockwave: {why}")
    def _last(*arrs):
        v = np.concatenate([np.asarray(a, dtype=float) for a in arrs])
        return np.nanmax(v)
    if comps == "shockwave":
        t0 = np.nanmin(t_sw) - tim["pre_roll"]
        t1 = _last(t_sw, t_swr) + tim["post_roll"]
    elif comps == "muzzle_blast":
        t0 = t_mb.min() - tim["pre_roll"]
        t1 = _last(t_mb, t_mbr) + tim["post_roll"]
    else:
        t0 = np.nanmin(np.concatenate([t_sw, t_mb])) - tim["pre_roll"]
        t1 = _last(t_mb, t_mbr, t_swr) + tim["post_roll"]
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
    T_Nr = np.full(n_mics, np.nan)
    dP_swr = np.full(n_mics, np.nan)
    t_posr = np.full(n_mics, np.nan)
    dP_mbr = np.full(n_mics, np.nan)

    for i in range(n_mics):
        sw = np.zeros(len(t_master))
        if shock_on and geo[i]["in_cone"]:
            g = geo[i]
            M_e = g["M_emit"] if drag else M          # local Mach at the emission point
            dP, T = gp.nwave_params(g["b"], M_e, c, bullet, shk["amp_scale"], shk["dur_scale"], p_atm)
            if gen_sw:
                sw = gp.nwave_waveform(t_master, g["t_arrive"], dP, T)
                if atmosphere is not None:
                    sw = gp.spectral_absorption(sw, FS, g["R"], atmosphere)
            if atmosphere is not None:
                att_sw[i] = float(gp.atmospheric_attenuation_db(1.0 / T, g["R"], **atmosphere))
            dP_sw[i], T_N[i] = dP, T

        dP, tp = gp.blast_params(r_mb[i], w_kg, mbc["amp_scale"], mbc["dur_scale"], p_atm)
        mb = np.zeros(len(t_master))
        if gen_mb:
            mb = gp.friedlander_waveform(t_master, t_mb[i], dP, tp)
            if atmosphere is not None:
                mb = gp.spectral_absorption(mb, FS, r_mb[i], atmosphere)
        if atmosphere is not None:
            att_mb[i] = float(gp.atmospheric_attenuation_db(1.0 / tp, r_mb[i], **atmosphere))
        dP_mb[i], t_pos[i] = dP, tp

        refl = np.zeros(len(t_master))
        if ground_on and swr_reason[i] is None:
            g = geo_r[i]
            M_e = g["M_emit"] if drag else M
            dP, T = gp.nwave_params(g["b"], M_e, c, bullet, shk["amp_scale"], shk["dur_scale"], p_atm)
            dP *= gnd["reflection_coeff_sw"]
            if gen_sw:
                w = gp.nwave_waveform(t_master, g["t_arrive"], dP, T)
                if atmosphere is not None:
                    w = gp.spectral_absorption(w, FS, g["R"], atmosphere)
                refl += w
            dP_swr[i], T_Nr[i] = dP, T
        if ground_on and mbr_reason[i] is None:
            dP, tp = gp.blast_params(r_mbr[i], w_kg, mbc["amp_scale"], mbc["dur_scale"], p_atm)
            dP *= gnd["reflection_coeff_mb"]
            if gen_mb:
                w = gp.friedlander_waveform(t_master, t_mbr[i], dP, tp)
                if atmosphere is not None:
                    w = gp.spectral_absorption(w, FS, r_mbr[i], atmosphere)
                refl += w
            dP_mbr[i], t_posr[i] = dP, tp
        signals.append(sw + mb + refl)

    # ---- ground truth ---------------------------------------------------------------------------
    per_mic = []
    for i in range(n_mics):
        g = geo[i]
        reason = None
        if not shock_on:
            reason = "below_mach_threshold"
        elif drag and g["subsonic"]:
            reason = "subsonic_at_emission"
        elif not g["in_cone"]:
            reason = "outside_mach_cone"
        sw_rec = dict(active=reason is None, inactive_reason=reason, generated=bool(gen_sw and reason is None))
        if g is not None:
            sw_rec.update(along_track_a_m=g["a"], perp_dist_b_m=g["b"], emission_coord_s_m=g["s"])
        if drag and g is not None and g["M_emit"] is not None:
            sw_rec.update(mach_emit=float(g["M_emit"]), bullet_speed_emit_mps=float(g["v_emit"]))
        if reason is None:
            sw_rec.update(t_arrive_s=float(t_sw[i]), t_emit_s=float(g["t_emit"]),
                          acoustic_path_R_m=float(g["R"]), duration_T_N_s=float(T_N[i]),
                          peak_dP_pa=float(dP_sw[i]), atmospheric_attenuation_db=float(att_sw[i]))
        rec = dict(
            mic_index=i, shockwave=sw_rec,
            muzzle_blast=dict(generated=bool(gen_mb), t_arrive_s=float(t_mb[i]), range_m=float(r_mb[i]),
                              positive_phase_t_pos_s=float(t_pos[i]), peak_dP_pa=float(dP_mb[i]),
                              atmospheric_attenuation_db=float(att_mb[i])))
        if ground_on:
            r = swr_reason[i]
            swr = dict(active=r is None, inactive_reason=r, generated=bool(gen_sw and r is None))
            if r is None:
                g = geo_r[i]
                swr.update(t_arrive_s=float(t_swr[i]),
                           delay_after_direct_s=(float(t_swr[i] - t_sw[i]) if np.isfinite(t_sw[i]) else None),
                           perp_dist_b_m=float(g["b"]), acoustic_path_R_m=float(g["R"]),
                           reflection_point_m=_reflection_point(g["x_emit"], img[i], z_g),
                           duration_T_N_s=float(T_Nr[i]), peak_dP_pa=float(dP_swr[i]),
                           reflection_coeff=float(gnd["reflection_coeff_sw"]))
            r = mbr_reason[i]
            mbr = dict(active=r is None, inactive_reason=r, generated=bool(gen_mb and r is None))
            if r is None:
                mbr.update(t_arrive_s=float(t_mbr[i]), delay_after_direct_s=float(t_mbr[i] - t_mb[i]),
                           range_m=float(r_mbr[i]), positive_phase_t_pos_s=float(t_posr[i]),
                           peak_dP_pa=float(dP_mbr[i]), reflection_coeff=float(gnd["reflection_coeff_mb"]))
            rec.update(shockwave_reflected=swr, muzzle_blast_reflected=mbr)
        per_mic.append(rec)

    active = ~np.isnan(t_sw)
    sep = (t_mb - t_sw)[active] if active.any() else np.array([])
    checks = dict(
        shock_before_blast_all_active=bool(np.all(sep > 0)) if sep.size else None,
        min_blast_minus_shock_s=float(sep.min()) if sep.size else None,
        events_overlap_any=bool(np.any(sep < T_N[active])) if sep.size else None,
        n_sensors_with_shock=int(active.sum()), n_sensors=int(n_mics))

    ground_truth = dict(
        components=comps,
        mic_pos=mic_pos.tolist(), shooter_pos=shooter.tolist(), bullet_origin=origin.tolist(),
        trajectory_direction=v_hat.tolist(), bullet_speed_mps=float(M * c), mach=float(M),
        speed_of_sound_mps=float(c), wind_mps=(wind.tolist() if wind is not None else [0.0, 0.0, 0.0]),
        sample_rate_hz=FS, atmosphere_applied=atmosphere,
        drag_model=bullet.get("drag", "none"),
        decay_length_m=(float(bullet["lam"]) if drag else None),
        subsonic_distance_m=(gp.subsonic_distance(M * c, bullet["lam"], max(shk["mach_threshold"], 1.0) * c)
                             if drag else None),
        models=dict(shockwave="whitham", muzzle_blast="hopkinson_cranz_kinney_graham",
                    absorption="iso_9613_1_minimum_phase", shock_branch_enabled=bool(shock_on),
                    mach_threshold=shk["mach_threshold"],
                    ground_reflection=("image_source_flat_ground_constant_coefficient" if ground_on else None)),
        ground=dict(enabled=ground_on, **({"array_height_m": gnd["array_height_m"],
                                          "ground_plane_z_m": -gnd["array_height_m"],
                                          "reflection_coeff_sw": gnd["reflection_coeff_sw"],
                                          "reflection_coeff_mb": gnd["reflection_coeff_mb"]} if ground_on else {})),
        t_master_start_s=float(t0), t_master_end_s=float(t1), n_samples=len(t_master),
        per_mic=per_mic,
        true_tdoa_sw_us=[None if (np.isnan(t_sw[i]) or np.isnan(t_sw[0])) else float((t_sw[i] - t_sw[0]) * 1e6)
                         for i in range(1, n_mics)],
        true_tdoa_mb_us=[float((t_mb[i] - t_mb[0]) * 1e6) for i in range(1, n_mics)],
        tdoa_matrix_sw_us=_matrix_us(t_sw), tdoa_matrix_mb_us=_matrix_us(t_mb),
        **({"tdoa_matrix_sw_reflected_us": _matrix_us(t_swr),
            "tdoa_matrix_mb_reflected_us": _matrix_us(t_mbr)} if ground_on else {}),
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
    gp.add_mic_config_argument(ap)
    gp.add_components_argument(ap)
    args = ap.parse_args()

    try:
        cfg = gp.load_config(args.config, args.mic_config, args.components)
    except ValueError as e:
        sys.exit(f"[Stage 1] ERROR: {e}")
    if args.output:
        cfg["output"]["basename"] = args.output
    print(f"[Stage 1] mic_config = {cfg['geometry']['mic_config']}, components = {cfg['signal']['components']}")
    try:
        signals, fs, gt = generate_clean(cfg)
    except ValueError as e:
        sys.exit(f"[Stage 1] ERROR: {e}")

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
    if gt["ground"]["enabled"]:
        d = [pm["shockwave_reflected"].get("delay_after_direct_s") for pm in gt["per_mic"]]
        d = [v * 1e3 for v in d if v is not None]
        print(f"[Stage 1] ground reflection on (array {gt['ground']['array_height_m']} m above ground): "
              + (f"reflected shock {min(d):.2f}-{max(d):.2f} ms after the direct one" if d else "no reflected shock"))
    print(f"[Stage 1] Next: python gs_gen_add_noise.py {wav_path} {args.config}")


if __name__ == "__main__":
    main()
