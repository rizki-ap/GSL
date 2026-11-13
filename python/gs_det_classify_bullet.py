#!/usr/bin/env python3
"""
gs_det_classify_bullet.py
============================
DETECTION STAGE 4. Reads Stage 2's TDOA output and estimates the bullet's
Mach number and caliber via three independent methods (adapted from the
validated Section 9 velocity-estimation approach):

  Method A -- angle only: k_hat_sw . v_hat ~ 1/M, using v_hat ~ -n_hat_shooter
    as a cheap approximation of bullet direction. No amplitude calibration
    needed, but least accurate (biased by miss distance).
  Method B -- duration + amplitude: fits Whitham's N-wave equations
        dP  = K_P p_atm (M^2-1)^(1/8) d / (L^(1/4) b^(3/4))
        T_N = K_T M d (b/L)^(1/4) / (c (M^2-1)^(3/8))
    for M and the miss distance b, for EACH caliber in the reference
    library (d = diameter, L = length). Needs this detector's own amplitude
    calibration (det_config.ini [calibration]) to convert normalized dP_sw
    back to Pa. Eliminating b leaves T_N ~ M (M^2-1)^(-1/3) d^(4/3) /
    (L^(1/3) dP^(1/3) c): the Mach factor varies only ~+-5 % for M = 1.4-2.6,
    so T_N and dP fix the miss distance and the CALIBER (via T_N dP^(1/3))
    well, but the Mach number only weakly -- so Method A's Mach is used (capped
    at the caliber's muzzle Mach v0/c) and Method B returns the miss distance
    and a per-caliber duration misfit. Checked on the Cipatat recordings
    (Nov 2026): Whitham's duration law matched within ~7 % (the previous
    formula was ~7x too long); with the true Mach, the misfit picked the right
    size class (5.56 / 7.62-.338 / 12.7) in 7 of 9 shots, but cannot separate
    calibers of similar d^(4/3)/L^(1/3) (7.62x39, 7.62x51, 7.62x54R, .338).
  Method C -- geometric correction: uses Method B's miss-distance estimate
    and a self-contained PILAR-V range solve to correct the approximate
    bullet direction, rather than assuming v_hat ~ -n_hat_shooter directly.

Caliber classification: Method B is run against EVERY caliber in the
reference library (det_config.ini: the six small-arms calibers of the
Cipatat development set); the caliber with the smallest |T_N misfit| is
reported as the best guess. It is reliable at the size-class level at best
-- treat it as a demonstration of the approach, not a validated classifier.

Usage
-----
  python gs_det_classify_bullet.py det_result_tdoa.json                # det_config.ini
  python gs_det_classify_bullet.py det_result_tdoa.json my_det_config.ini
  python gs_det_classify_bullet.py det_result_tdoa.json -mic_config=UMA8   # normally not needed

-mic_config defaults to the preset recorded in *_tdoa.json (i.e. the one Stage 1 used).
"""

import argparse
import json
from datetime import datetime, timezone

import numpy as np

from gs_det_signal_prepare import (load_det_config, build_mic_positions,
                                   add_mic_config_argument, upstream_mic_config)
from gs_det_shooter_locator import solve_doa, pilar_v_exact


def method_a_angle_only(k_hat_sw, n_hat_shooter):
    v_hat_approx = -n_hat_shooter
    dot = float(np.dot(k_hat_sw, v_hat_approx))
    return 1.0 / dot if dot > 0 else None


WHITHAM_K_P = 0.53     # peak-overpressure coefficient (same published constant as the generator,
WHITHAM_K_T = 1.82     # duration coefficient            kept as an independent copy here)
M_SPLIT = np.sqrt(3.0)  # T_N(M) at fixed dP has its minimum here
M_MAX_SEARCH = 6.0


def whitham_b_from_dp(dP, M, bl, p_atm):
    """Miss distance (m) that gives peak overpressure dP (Pa) at Mach M (Whitham)."""
    return (WHITHAM_K_P * p_atm * (M * M - 1) ** 0.125 * bl["d"] / (bl["L"] ** 0.25 * dP)) ** (4.0 / 3.0)


def whitham_tn(M, b, c, bl):
    """N-wave duration (s) at Mach M and miss distance b (Whitham)."""
    return WHITHAM_K_T * M * bl["d"] * (b / bl["L"]) ** 0.25 / (c * (M * M - 1) ** 0.375)


def method_b_duration_amplitude(T_N, dP, c, bl, p_atm=101325.0, M_hint=None, M_max=None):
    """Whitham fit of the measured N-wave (T_N, dP) for one candidate caliber (needs bl["d"], bl["L"]).

    Eliminating b gives T_N ~ M (M^2-1)^(-1/3) * d^(4/3) / (L^(1/3) dP^(1/3) c). The Mach factor
    changes only about +-5 % over M = 1.4-2.6 (minimum at sqrt(3)), so (T_N, dP) pin the Mach
    number poorly, while T_N * dP^(1/3) is a caliber signature. Hence:
      * M: Method A's Mach (M_hint) when available, kept within (1, M_max]; otherwise the Mach
        that best fits T_N (the faster root when there are two);
      * b: the miss distance that gives dP at that M;
      * misfit_db = 20 log10(T_N model / T_N measured) at (M, b) -- the caliber score
        (0 dB = this caliber explains the N-wave).
    Also returns M_fit, the Mach that alone best fits T_N, and its root(s), for information.
    Returns None for invalid input."""
    if (T_N is None or not np.isfinite(T_N) or T_N <= 0 or dP is None or dP <= 0
            or not bl.get("d") or not bl.get("L")):
        return None
    m_hi = min(M_MAX_SEARCH, M_max) if M_max else M_MAX_SEARCH
    if m_hi <= 1.01:
        return None

    def log_err(M):
        return np.log(whitham_tn(M, whitham_b_from_dp(dP, M, bl, p_atm), c, bl) / T_N)

    Ms = np.linspace(1.01, m_hi, 600)
    err = np.array([log_err(M) for M in Ms])
    roots = [float(Ms[i]) for i in range(1, len(Ms)) if err[i - 1] * err[i] <= 0]
    M_fit = max(roots) if roots else float(Ms[int(np.argmin(np.abs(err)))])
    if M_hint is not None and np.isfinite(M_hint) and M_hint > 1:
        M, mode = float(min(max(M_hint, 1.01), m_hi)), "method_a_mach"
    else:
        M, mode = M_fit, "duration_fit"
    return dict(M=M, mach_source=mode, b=float(whitham_b_from_dp(dP, M, bl, p_atm)),
                misfit_db=float(20 * np.log10(np.exp(log_err(M)))),
                M_fit=float(M_fit), roots=[(m, float(whitham_b_from_dp(dP, m, bl, p_atm))) for m in roots],
                ambiguous=len(roots) > 1)


def method_c_geometric_correction(k_hat_sw, n_hat_shooter, b_corr, M_B, dt_sw_to_mb, c):
    """Refine v_hat using a self-contained PILAR-V range solve (independent
    of gs_det_shooter_locator.py's own range estimate -- this is its own
    internal solve, matching how the notebook's Method C was self-contained
    too)."""
    if dt_sw_to_mb is None or M_B is None or M_B <= 1.0:
        return None
    try:
        r_est = pilar_v_exact(dt_sw_to_mb, b_corr, M_B, c)
    except (ValueError, RuntimeError):
        return None
    sin_alpha = b_corr / np.sqrt(r_est ** 2 + b_corr ** 2)
    cos_alpha = np.sqrt(max(0.0, 1 - sin_alpha ** 2))
    k_sw_along = np.dot(k_hat_sw, n_hat_shooter) * n_hat_shooter
    e_lat = -(k_hat_sw - k_sw_along)
    norm = np.linalg.norm(e_lat)
    if norm < 1e-9:
        return None
    e_lat = e_lat / norm
    v_hat_corr = cos_alpha * (-n_hat_shooter) + sin_alpha * e_lat
    v_hat_corr /= np.linalg.norm(v_hat_corr)
    dot_C = float(np.dot(k_hat_sw, v_hat_corr))
    return 1.0 / dot_C if dot_C > 0 else None


def main():
    parser = argparse.ArgumentParser(
        description="Detection Stage 4: bullet Mach/caliber classification (3 methods).")
    parser.add_argument("tdoa_json", help="Path to Stage 2's *_tdoa.json")
    parser.add_argument("config", nargs="?", default="det_config.ini",
                         help="Path to detection .ini config (default: det_config.ini)")
    parser.add_argument("-o", "--output", default=None,
                         help="Output basename (overrides output.basename in config)")
    add_mic_config_argument(parser)
    args = parser.parse_args()

    with open(args.tdoa_json) as f:
        tdoa_data = json.load(f)

    cfg = load_det_config(args.config, args.mic_config or upstream_mic_config(tdoa_data))
    if args.output:
        cfg["output"]["basename"] = args.output

    mic_pos = build_mic_positions(cfg)
    side = cfg["geometry"]["planar_source_side"]
    if tdoa_data.get("n_channels", len(mic_pos)) != len(mic_pos):
        raise SystemExit(f"[Det Stage 4] ERROR: TDOA file has {tdoa_data['n_channels']} channels but "
                         f"mic_config {cfg['geometry']['mic_config']} has {len(mic_pos)} mics.")
    c = cfg["physics"]["speed_of_sound"]
    ref_ch = tdoa_data["reference_channel"]

    result = dict(
        generated_at_utc=datetime.now(timezone.utc).isoformat(),
        stage="bullet_classification",
        source_tdoa_json=args.tdoa_json,
    )

    if not tdoa_data.get("sw_available", bool(tdoa_data.get("tdoa_sw_s"))):
        print("[Det Stage 4] No shockwave TDOA available (blast-only event) -- classification needs "
              "the shockwave. Aborting.")
        result["classification_available"] = False
        result["best_guess_caliber_category"] = None
        result["best_guess_caliber_category_reason"] = "no shockwave: muzzle blast only"
        with open(f"{cfg['output']['basename']}_classification.json", "w") as f:
            json.dump(result, f, indent=2)
        return

    k_hat_sw, _ = solve_doa(mic_pos, tdoa_data["tdoa_sw_s"], c, ref_ch, side)

    if not tdoa_data["mb_available"]:
        print("[Det Stage 4] No muzzle-blast TDOA available -- classification needs bearing, "
              "which needs MB TDOA. Aborting.")
        result["classification_available"] = False
        result["best_guess_caliber_category"] = None
        result["best_guess_caliber_category_reason"] = "no muzzle blast: no bearing"
        with open(f"{cfg['output']['basename']}_classification.json", "w") as f:
            json.dump(result, f, indent=2)
        return

    n_hat_shooter, _ = solve_doa(mic_pos, tdoa_data["tdoa_mb_s"], c, ref_ch, side)
    n_hat_shooter = -n_hat_shooter

    T_N = tdoa_data.get("t_n_sw_ref_channel")
    dp_sw_norm = tdoa_data.get("dp_sw_ref_channel")
    norm_peak = tdoa_data.get("normalization_peak")
    if norm_peak is None:
        norm_peak = 1.0
        print("[Det Stage 4] WARNING: no normalization_peak in the TDOA file (older Stage 2) -- "
              "dP is taken as full-scale units; re-run Stage 2 for a correct amplitude.")
    dp_sw_pa = (abs(dp_sw_norm) * norm_peak * cfg["calibration"]["assumed_peak_pa_at_full_scale"]
                if dp_sw_norm is not None else None)
    dt_sw_to_mb = tdoa_data.get("dt_sw_to_mb_s")

    M_A = method_a_angle_only(k_hat_sw, n_hat_shooter)
    result["method_a_mach"] = M_A
    print(f"[Det Stage 4] Method A (angle only): Mach = {M_A}")

    print(f"[Det Stage 4] Testing {len(cfg['bullet_reference_library'])} candidate caliber(s) "
          f"against Method B (T_N={T_N}, dP={dp_sw_pa} Pa assumed)...")
    caliber_results = {}
    p_atm = cfg["physics"].get("atmospheric_pressure_pa", 101325.0)
    for cal, bl in cfg["bullet_reference_library"].items():
        M_max = 1.02 * bl["v0"] / c if bl.get("v0") else None
        mb_result = method_b_duration_amplitude(T_N, dp_sw_pa, c, bl, p_atm, M_hint=M_A, M_max=M_max)
        if mb_result is None:
            caliber_results[cal] = dict(mach=None, miss_distance_m=None, valid=False)
            print(f"    {cal}: no valid supersonic solution")
            continue
        M_B, b_corr = mb_result["M"], mb_result["b"]
        M_C = method_c_geometric_correction(k_hat_sw, n_hat_shooter, b_corr, M_B, dt_sw_to_mb, c)
        caliber_results[cal] = dict(mach_method_b=M_B, mach_method_c=M_C,
                                     miss_distance_m=b_corr, valid=True,
                                     duration_misfit_db=mb_result["misfit_db"],
                                     mach_source=mb_result["mach_source"], mach_duration_fit=mb_result["M_fit"],
                                     method_b_roots=mb_result["roots"], ambiguous=mb_result["ambiguous"])
        print(f"    {cal}: Method B Mach = {M_B:.3f} ({mb_result['mach_source']}), Method C Mach = {M_C}, "
              f"b = {b_corr:.2f} m, T_N misfit = {mb_result['misfit_db']:+.2f} dB")

    result["caliber_candidates"] = caliber_results

    valid_calibers = {k: v for k, v in caliber_results.items() if v["valid"]}
    if valid_calibers:
        # Best guess: the caliber whose Whitham model fits the measured N-wave duration best
        # (smallest |T_N misfit|, i.e. the best match of the T_N * dP^(1/3) signature).
        # Calibers of similar d^(4/3) / L^(1/3) (the 7.62 family, .338) are hard to tell
        # apart, so don't oversell this as a validated classifier (see module docstring).
        best_cal = min(valid_calibers, key=lambda k: abs(valid_calibers[k]["duration_misfit_db"]))
        # Reported calibre (project rule): category mandatory, class optional
        # ([output] report_caliber_class), variant not reported (only in caliber_candidates).
        best_cls = cfg["bullet_reference_library"][best_cal].get("cls", best_cal)
        best_cat = cfg["bullet_reference_library"][best_cal].get("cat", "other")
        result["best_guess_caliber_category"] = best_cat
        if cfg["output"]["report_caliber_class"]:
            result["best_guess_caliber_class"] = best_cls
        print(f"[Det Stage 4] Best-guess calibre category: {best_cat} (class {best_cls})")
    else:
        result["best_guess_caliber_category"] = None
        result["best_guess_caliber_category_reason"] = "no calibre produced a valid solution"
        print("[Det Stage 4] No caliber produced a valid solution.")

    out_json = f"{cfg['output']['basename']}_classification.json"
    result["classification_available"] = True
    with open(out_json, "w") as f:
        json.dump(result, f, indent=2)
    print(f"[Det Stage 4] Wrote {out_json}")


if __name__ == "__main__":
    main()
