#!/usr/bin/env python3
"""
gs_det_shooter_locator.py
============================
DETECTION STAGE 3. Reads Stage 2's TDOA output and computes the shooter's
bearing (from muzzle-blast TDOA) and, if muzzle blast TDOA is available,
range (via the PILAR-V formula) and full 3D position.

Self-contained: range needs a Mach-number estimate, which gs_det_
classify_bullet.py also produces (more thoroughly, via 3 independent
methods) -- but this script doesn't depend on that output, so it can run
right after gs_det_tdoa.py. It gets its own quick Mach estimate via the
angle-only relation (k_hat_sw . n_hat_shooter ~ 1/M), which needs nothing
but the TDOA vectors already computed in Stage 2 -- no waveform amplitude
measurements, no bullet reference library required. Run gs_det_classify_
bullet.py afterward for a more accurate, amplitude-informed Mach/caliber
estimate if you want a refined range.

Usage
-----
  python gs_det_shooter_locator.py det_result_tdoa.json                # det_config.ini
  python gs_det_shooter_locator.py det_result_tdoa.json my_det_config.ini
  python gs_det_shooter_locator.py det_result_tdoa.json -mic_config=UMA8   # normally not needed

-mic_config defaults to the preset recorded in *_tdoa.json (i.e. the one Stage 1 used).

Planar arrays (e.g. UMA8, all mics at z = 0) only measure the in-plane part of the
direction; the out-of-plane part is |k_perp| = sqrt(1 - |k_par|^2) with its sign taken from
the preset's planar_source_side ('above' = sources on the +z side). Azimuth is fully
determined; elevation is a magnitude-only estimate and is poorly conditioned near the horizon.
"""

import argparse
import json
from datetime import datetime, timezone

import numpy as np
from scipy.optimize import brentq

from gs_det_signal_prepare import (load_det_config, build_mic_positions,
                                   add_mic_config_argument, upstream_mic_config)


def solve_doa(mic_pos, tdoa_dict, c, ref_ch=0, planar_source_side="above"):
    """Solve the plane-wave linear system delta_p . k_hat = c*tau for the
    propagation direction k_hat, using the N-1 pairs relative to ref_ch
    (least squares; for the 4-mic tetrahedron this is the exact 3x3 solve).

    Planar array (delta_p of rank 2, e.g. UMA8): only the in-plane component
    k_par is observable. The out-of-plane component is sqrt(1 - |k_par|^2)
    along the array normal (oriented towards +z), signed by planar_source_side:
    'above' = source on the +normal side, so the wave travels towards -normal.

    Returns (k_hat, |k|). |k| should be ~1.0 for consistent TDOA; for a planar
    array it is |k_par| when that exceeds 1 (inconsistent TDOA), else 1.0."""
    mic_pos = np.asarray(mic_pos, dtype=float)
    other = [i for i in range(len(mic_pos)) if i != ref_ch]
    delta_p = np.array([mic_pos[i] - mic_pos[ref_ch] for i in other])
    d = []
    for i in other:
        a, b = min(ref_ch, i), max(ref_ch, i)
        key = f"{a},{b}"
        sign = 1.0 if a == ref_ch else -1.0
        d.append(sign * tdoa_dict[key] * c)
    d = np.array(d)

    _, s, vt = np.linalg.svd(delta_p)
    rank = int(np.sum(s > s[0] * 1e-6))
    if rank < 2:
        raise ValueError("mic positions are collinear -- a direction cannot be resolved")
    k_vec = np.linalg.lstsq(delta_p, d, rcond=None)[0]   # minimum-norm: in-plane if rank 2
    if rank == 3:
        k_hat = k_vec / np.linalg.norm(k_vec)
        return k_hat, float(np.linalg.norm(k_vec))

    normal = vt[2]
    if normal[np.argmax(np.abs(normal))] < 0:             # orient towards +z (or the dominant axis)
        normal = -normal
    k_par = k_vec - np.dot(k_vec, normal) * normal
    n_par = float(np.linalg.norm(k_par))
    if n_par >= 1.0:
        return k_par / n_par, n_par
    k_perp = np.sqrt(1.0 - n_par ** 2) * (-1.0 if planar_source_side == "above" else 1.0)
    k_hat = k_par + k_perp * normal
    return k_hat / np.linalg.norm(k_hat), 1.0


def pilar_v_exact(dt_val, b, M, c, r_max=5000.0):
    beta = np.sqrt(M ** 2 - 1.0)

    def f(r):
        return np.sqrt(r ** 2 + b ** 2) / c - (r + b * beta) / (M * c) - dt_val

    return brentq(f, 0.1, r_max)


def main():
    parser = argparse.ArgumentParser(
        description="Detection Stage 3: DOA + PILAR-V range -> 3D shooter location.")
    parser.add_argument("tdoa_json", help="Path to Stage 2's *_tdoa.json")
    parser.add_argument("config", nargs="?", default="det_config.ini",
                         help="Path to detection .ini config (default: det_config.ini)")
    parser.add_argument("-o", "--output", default=None,
                         help="Output basename (overrides output.basename in config)")
    parser.add_argument("--assumed-miss-distance", type=float, default=None,
                         help="Fallback miss distance b (m) for the range solve, ONLY used if "
                              "gs_det_classify_bullet.py hasn't already provided a measured one. "
                              "This is a genuine detector-design assumption (a typical/expected "
                              "engagement envelope), not generation ground truth.")
    add_mic_config_argument(parser)
    args = parser.parse_args()

    with open(args.tdoa_json) as f:
        tdoa_data = json.load(f)

    cfg = load_det_config(args.config, args.mic_config or upstream_mic_config(tdoa_data))
    if args.output:
        cfg["output"]["basename"] = args.output

    geom = cfg["geometry"]
    mic_pos = build_mic_positions(cfg)
    side = geom["planar_source_side"]
    if tdoa_data.get("n_channels", len(mic_pos)) != len(mic_pos):
        raise SystemExit(f"[Det Stage 3] ERROR: TDOA file has {tdoa_data['n_channels']} channels but "
                         f"mic_config {geom['mic_config']} has {len(mic_pos)} mics.")
    print(f"[Det Stage 3] mic_config = {geom['mic_config']} ({len(mic_pos)} mics)")
    c = cfg["physics"]["speed_of_sound"]
    ref_ch = tdoa_data["reference_channel"]

    result = dict(
        generated_at_utc=datetime.now(timezone.utc).isoformat(),
        stage="shooter_location",
        source_tdoa_json=args.tdoa_json,
        mic_config=geom["mic_config"],
    )

    k_hat_sw, k_norm_sw = solve_doa(mic_pos, tdoa_data["tdoa_sw_s"], c, ref_ch, side)
    result["k_hat_sw"] = k_hat_sw.tolist()
    result["k_hat_sw_norm_check"] = k_norm_sw  # should be ~1.0 for consistent TDOA
    print(f"[Det Stage 3] k_hat_sw = {k_hat_sw}  (|k| = {k_norm_sw:.4f}, should be ~1.0)")

    if not tdoa_data["mb_available"]:
        print("[Det Stage 3] No reliable muzzle-blast TDOA available -- bearing/range/position "
              "cannot be computed, only the bullet trajectory direction (k_hat_sw) above.")
        result["bearing_available"] = False
    else:
        n_hat_shooter, k_norm_mb = solve_doa(mic_pos, tdoa_data["tdoa_mb_s"], c, ref_ch, side)
        n_hat_shooter = -n_hat_shooter  # propagation direction -> direction TO the source
        az = float(np.degrees(np.arctan2(n_hat_shooter[1], n_hat_shooter[0])))
        el = float(np.degrees(np.arcsin(np.clip(n_hat_shooter[2], -1, 1))))
        result["bearing_available"] = True
        result["n_hat_shooter"] = n_hat_shooter.tolist()
        result["azimuth_deg"] = az
        result["elevation_deg"] = el
        print(f"[Det Stage 3] Shooter bearing: azimuth={az:.2f} deg, elevation={el:.2f} deg "
              f"(|k| = {k_norm_mb:.4f}, should be ~1.0)")

        # Self-contained rough Mach + range estimate (Method A: angle-only for
        # Mach, needs nothing but the TDOA-derived direction vectors above;
        # range then needs Dt and a miss-distance b, using this detector's own
        # reference bullet library with a DEFAULT assumed caliber -- NOT a
        # measurement of the true caliber, just a working assumption so this
        # stage can run standalone. gs_det_classify_bullet.py refines both the
        # Mach estimate (3 independent methods) and the caliber identification;
        # re-run this stage with --assumed-miss-distance from its output for a
        # better range if you want that refinement.)
        v_hat_approx = -n_hat_shooter
        dot_A = float(np.dot(k_hat_sw, v_hat_approx))
        M_A = 1.0 / dot_A if dot_A > 0 else None
        result["mach_estimate_angle_only"] = M_A

        dt_sw_to_mb = tdoa_data.get("dt_sw_to_mb_s")
        dp_sw_norm = tdoa_data.get("dp_sw_ref_channel")
        # dp_sw_norm is in NORMALIZED units (Stage 1 normalizes the wav);
        # the bullet library's dP0_sw is in Pa -- convert using this
        # detector's own assumed hardware calibration before mixing them.
        dp_sw = (abs(dp_sw_norm) * cfg["calibration"]["assumed_peak_pa_at_full_scale"]
                 if dp_sw_norm is not None else None)
        default_caliber = next(iter(cfg["bullet_reference_library"]))
        bl = cfg["bullet_reference_library"][default_caliber]

        if M_A is not None and M_A > 1.01 and dt_sw_to_mb is not None:
            if args.assumed_miss_distance is not None:
                b_est = args.assumed_miss_distance
                b_source = "user-provided --assumed-miss-distance"
            elif dp_sw is not None and dp_sw > 0:
                b_est = bl["b0_sw"] * (bl["dP0_sw"] / dp_sw) ** 2
                b_source = f"amplitude-derived, assuming default caliber {default_caliber}"
            else:
                b_est = 10.0
                b_source = "no amplitude measurement available -- arbitrary 10m fallback"

            try:
                r_est = pilar_v_exact(dt_sw_to_mb, b_est, M_A, c)
                r_3d = float(np.sqrt(r_est ** 2 + b_est ** 2))
                shooter_pos_est = (r_3d * n_hat_shooter).tolist()
                result["range_estimate_m"] = r_est
                result["miss_distance_assumed_m"] = b_est
                result["miss_distance_source"] = b_source
                result["shooter_position_estimate_m"] = shooter_pos_est
                print(f"[Det Stage 3] Mach (angle-only) = {M_A:.3f}, b assumed = {b_est:.2f} m "
                      f"({b_source})")
                print(f"[Det Stage 3] PILAR-V range = {r_est:.2f} m, "
                      f"3D position estimate = {shooter_pos_est}")
            except (ValueError, RuntimeError) as e:
                print(f"[Det Stage 3] Range solve failed ({e}) -- bearing-only result reported.")
                result["range_estimate_m"] = None
        else:
            reason = "M<=1 (not physical)" if (M_A is None or M_A <= 1.01) else "no Dt available"
            print(f"[Det Stage 3] Skipping range solve ({reason}). "
                  "Run gs_det_classify_bullet.py for a more robust Mach/range estimate.")
            result["range_estimate_m"] = None

    out_json = f"{cfg['output']['basename']}_location.json"
    with open(out_json, "w") as f:
        json.dump(result, f, indent=2)
    print(f"[Det Stage 3] Wrote {out_json}")


if __name__ == "__main__":
    main()
