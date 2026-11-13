#!/usr/bin/env python3
"""
gs_det_shooter_locator.py
============================
DETECTION STAGE 3. Reads Stage 2's TDOA output and computes the shooter's
bearing (from the muzzle-blast TDOA) and, when a muzzle blast is available,
the range, the bullet trajectory and the full 3D shooter position.

Range models ([range] model in det_config.ini, or --range-model):

  drag (default)
    Drag-aware fit, v(s) = v0 exp(-s / lambda), with v0 and lambda taken from
    the detector's OWN ballistic table ([bullet_reference_library] in det_config.ini;
    never from the generator). Measurements used:
      n_hat  direction to the shooter          (muzzle-blast TDOA)
      k_hat  shockwave propagation direction   (shockwave TDOA)
      Dt     SW -> MB delay on the reference channel
      T_N    N-wave duration (optional; ignored if physically implausible)
    The shooter, the trajectory and the array lie in one plane. For a trial
    shooter distance D and trajectory angle alpha (miss distance b = D sin alpha)
    the forward model finds the shock emission point from the LOCAL Mach-angle
    condition cos(phi) = c / v(s) (Fermat: first-arriving shock) and predicts
    theta (the angle between n_hat and -k_hat), Dt and T_N (Whitham). (D, alpha)
    are fitted per caliber by minimising
      J = (dtheta/sigma_theta)^2 + (dDt/sigma_dt)^2 + (ln T_N ratio/sigma_T)^2 + (b/sigma_b)^2
    (grid search, then Nelder-Mead). The caliber with the lowest J is reported
    with its range, miss distance, trajectory and bullet speed at emission.
    J <~ 1: consistent; J > 4: warning (wrong table, caliber or measurement).
    Without a usable T_N, the angle and Dt fit every fast-enough caliber
    equally well, so [range] default_caliber is used and the choice is flagged
    as ambiguous; the range itself depends only weakly on the caliber.
    validated against gs_gen_physic.shock_geometry_drag: test_det_range_drag.py.

  constant
    The original constant-speed PILAR-V solve (Mach from the angle-only
    relation, miss distance from SW amplitude). Kept for comparison: it
    overestimates range for real (slowing) bullets -- see drag_model.md.

Speed of sound: [physics] temperature_c (or --temperature-c) sets
c = 331.3 + 0.606 T; otherwise [physics] speed_of_sound is used.

Usage
-----
  python gs_det_shooter_locator.py det_result_tdoa.json                       # det_config.ini
  python gs_det_shooter_locator.py det_result_tdoa.json my_det_config.ini
  python gs_det_shooter_locator.py det_result_tdoa.json --temperature-c 26
  python gs_det_shooter_locator.py det_result_tdoa.json --caliber 7.62_NATO
  python gs_det_shooter_locator.py det_result_tdoa.json --range-model constant   # old solve
  python gs_det_shooter_locator.py det_result_tdoa.json -mic_config=UMA8      # normally not needed

-mic_config defaults to the preset recorded in *_tdoa.json (i.e. the one Stage 1 used).

Several shots (Stage 1 [detection] multi_shot, Stage 2 'shots'): every shot is located and
written to 'shots' in *_location.json; the top-level result is the first shot. A shot far
from the median bearing of the file is flagged 'bearing_outlier'.

Project limits ([range]): shots are expected at 50-2000 m. A shot located closer than
friendly_range_m (50 m) is flagged 'friendly', one beyond max_report_range_m (2000 m)
'beyond_max_range'. The range search itself runs from min_range_m to max_range_m.

Planar arrays (e.g. UMA8, all mics at z = 0) only measure the in-plane part of the
direction; the out-of-plane part is |k_perp| = sqrt(1 - |k_par|^2) with its sign taken from
the preset's planar_source_side ('above' = sources on the +z side). Azimuth is fully
determined; elevation is a magnitude-only estimate and is poorly conditioned near the horizon.
"""

import argparse
import configparser
import json
import math
from datetime import datetime, timezone

import numpy as np
from scipy.optimize import brentq, minimize

from gs_det_signal_prepare import (load_det_config, build_mic_positions,
                                   add_mic_config_argument, upstream_mic_config)

# Whitham N-wave duration coefficient (public constant, Whitham 1952):
#   T_N = K_T * M * d * (b / L)^(1/4) / (c * (M^2 - 1)^(3/8))
WHITHAM_K_T = 1.82
# Whitham peak-overpressure coefficient (legacy constant-speed solve only):
#   dP = K_P * p_atm * (M^2 - 1)^(1/8) * d / (L^(1/4) * b^(3/4))
WHITHAM_K_P = 0.53

# Detector-side defaults, used when det_config.ini has no [range] section. The ballistic
# table is the detector's [bullet_reference_library] (see det_config.ini), maintained
# independently of the generator's BULLET_LIBRARY (fairness rule).
P_ATM = 101325.0          # Pa, ambient pressure in the Whitham amplitude law
RANGE_DEFAULTS = {"model": "drag", "default_caliber": "5.56_NATO", "mach_min": "1.1",
                  "min_range_m": "5.0", "max_range_m": "2200.0", "n_grid": "240",
                  "friendly_range_m": "50.0", "max_report_range_m": "2000.0",
                  "lam_temperature_scaling": "true", "t_n_offset_s": "0.0",
                  "dp_sigma_log": "0.0",
                  "direction_sigma_deg": "1.5", "dt_sigma_ms": "0.2", "t_n_sigma_log": "0.3",
                  "miss_sigma_m": "50"}
# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

def load_range_config(path, bullet_library):
    """[range] and [physics] temperature_c from the detection .ini (read here so that
    gs_det_signal_prepare.load_det_config stays unchanged). The ballistic table is built
    from the detector's [bullet_reference_library] (bullet_library, as returned by
    load_det_config): every calibre with a muzzle velocity v0; lam missing or 0 = constant
    speed."""
    cp = configparser.ConfigParser()
    cp.read_dict({"range": RANGE_DEFAULTS, "physics": {"temperature_c": ""}})
    if path is not None:
        cp.read(path)
    if cp.has_section("ballistics"):
        print("[Det Stage 3] NOTE: [ballistics] is no longer read -- the drag fit uses "
              "[bullet_reference_library] (<cal>_v0, <cal>_lam, <cal>_d, <cal>_L)")

    table = {}
    for cal, bl in bullet_library.items():
        if not bl.get("v0"):
            print(f"[Det Stage 3] NOTE: {cal} has no v0 in [bullet_reference_library] -- left out of the drag fit")
            continue
        lam = bl.get("lam")
        table[cal] = dict(v0=float(bl["v0"]), lam=np.inf if not lam else float(lam),
                          d=float(bl["d"]), L=float(bl["L"]), cls=bl.get("cls", cal),
                          cat=bl.get("cat", "other"))
    t_raw = cp.get("physics", "temperature_c").strip()
    return dict(model=cp.get("range", "model").strip().lower(),
                default_caliber=cp.get("range", "default_caliber").strip(),
                mach_min=cp.getfloat("range", "mach_min"),
                min_range_m=cp.getfloat("range", "min_range_m"),
                max_range_m=cp.getfloat("range", "max_range_m"),
                friendly_range_m=cp.getfloat("range", "friendly_range_m"),
                max_report_range_m=cp.getfloat("range", "max_report_range_m"),
                lam_temperature_scaling=cp.getboolean("range", "lam_temperature_scaling"),
                t_n_offset_s=cp.getfloat("range", "t_n_offset_s"),
                dp_sigma_log=cp.getfloat("range", "dp_sigma_log"),
                n_grid=cp.getint("range", "n_grid"),
                direction_sigma_deg=cp.getfloat("range", "direction_sigma_deg"),
                dt_sigma_ms=cp.getfloat("range", "dt_sigma_ms"),
                t_n_sigma_log=cp.getfloat("range", "t_n_sigma_log"),
                miss_sigma_m=cp.getfloat("range", "miss_sigma_m"),
                temperature_c=float(t_raw) if t_raw else None,
                ballistics=table)


def speed_of_sound_from_temperature(t_c):
    """Speed of sound in dry air (m/s), linear approximation."""
    return 331.3 + 0.606 * t_c


# ---------------------------------------------------------------------------
# Directions
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Range: constant speed (legacy; also imported by gs_det_classify_bullet.py)
# ---------------------------------------------------------------------------

def pilar_v_exact(dt_val, b, M, c, r_max=5000.0):
    """Constant-speed SW->MB delay solved for r = along-track distance from the
    shooter to the closest point of approach (miss distance b, Mach M)."""
    beta = np.sqrt(M ** 2 - 1.0)

    def f(r):
        return np.sqrt(r ** 2 + b ** 2) / c - (r + b * beta) / (M * c) - dt_val

    return brentq(f, 0.1, r_max)


# ---------------------------------------------------------------------------
# Range: drag-aware
# ---------------------------------------------------------------------------

def bullet_speed(s, v0, lam):
    """Speed (m/s) after s metres, v0 exp(-s / lam); constant if lam is inf."""
    return v0 if not np.isfinite(lam) else v0 * np.exp(-s / lam)


def bullet_time(s, v0, lam):
    """Time of flight (s) to s metres, lam (exp(s / lam) - 1) / v0; s / v0 if lam is inf."""
    return s / v0 if not np.isfinite(lam) else lam * np.expm1(s / lam) / v0


def forward_model(D, alpha, c, v0, lam, mach_min=1.0, n_iter=60):
    """Forward model in the plane of the shooter, the trajectory and the array.

    Shooter at distance D from the array; the trajectory leaves the shooter at
    angle alpha from the line of sight, so it passes the array at miss distance
    b = D sin(alpha). The first shock reaching the array is emitted at the path
    length s* that minimises t(s) + |x - p(s)|/c (Fermat), i.e. where
    1/v(s) = cos(angle between path and the ray to the array)/c  -- the LOCAL
    Mach-angle condition, solved by bisection (unique root on [0, a]).

    Works element-wise on numpy arrays. Returns a dict of arrays:
      theta  angle at the array between the shooter and the apparent SW source
      dt     muzzle-blast minus shockwave arrival time
      s, v, rho (emission path length, bullet speed there, emission-to-array distance),
      miss, valid (shock reaches the array and is emitted while supersonic)."""
    D, alpha = np.broadcast_arrays(np.asarray(D, float), np.asarray(alpha, float))
    a = D * np.cos(alpha)
    b = D * np.sin(alpha)
    inside = (v0 * np.cos(alpha) > c) & (a > 0)
    lo, hi = np.zeros_like(a), np.maximum(a, 0.0)
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        fp = 1.0 / bullet_speed(mid, v0, lam) - (a - mid) / (c * np.hypot(a - mid, b) + 1e-300)
        neg = fp < 0
        lo = np.where(neg, mid, lo)
        hi = np.where(neg, hi, mid)
    s = 0.5 * (lo + hi)
    v = bullet_speed(s, v0, lam) * np.ones_like(s)
    xe, ye = D - s * np.cos(alpha), s * np.sin(alpha)          # emission point (shooter on +x)
    rho = np.hypot(xe, ye)
    theta = np.arctan2(ye, xe)
    dt = D / c - bullet_time(s, v0, lam) - rho / c
    return dict(theta=theta, dt=dt, s=s, v=v, rho=rho, miss=b,
                valid=inside & (v > mach_min * c))


def _forward_scalar(D, alpha, c, v0, lam, mach_min=1.0, n_iter=50):
    """Scalar (pure-Python, fast) version of forward_model; None if not valid."""
    a, b = D * math.cos(alpha), D * math.sin(alpha)
    if a <= 0 or v0 * math.cos(alpha) <= c:
        return None
    inv_lam = 0.0 if not np.isfinite(lam) else 1.0 / lam
    lo, hi = 0.0, a
    for _ in range(n_iter):
        mid = 0.5 * (lo + hi)
        if math.exp(mid * inv_lam) / v0 - (a - mid) / (c * math.hypot(a - mid, b) + 1e-300) < 0:
            lo = mid
        else:
            hi = mid
    s = 0.5 * (lo + hi)
    v = v0 * math.exp(-s * inv_lam)
    if v <= mach_min * c:
        return None
    xe, ye = D - s * math.cos(alpha), s * math.sin(alpha)
    rho = math.hypot(xe, ye)
    t = s / v0 if inv_lam == 0.0 else lam * math.expm1(s * inv_lam) / v0
    return dict(theta=math.atan2(ye, xe), dt=D / c - t - rho / c, s=s, v=v, rho=rho, miss=b)


def fit_range_drag(dt_meas, theta_meas, c, bl, rcfg, t_n_meas=None, max_starts=4, dp_meas=None):
    """Fit (D, alpha) for one caliber (bl: v0, lam, d, L) so the forward model
    reproduces the measurements. Cost (sum of squared normalised residuals):
      SW angle   (theta_pred - theta) / direction_sigma_deg
      delay      (Dt_pred - Dt) / dt_sigma_ms
      N-wave     ln(T_N_pred / T_N) / t_n_sigma_log      (only if T_N is usable)
      amplitude  ln(dP_pred / dP) / dp_sigma_log         (only with a SW peak in Pa)
      miss prior miss distance / miss_sigma_m             (0 disables)
    The angle and the delay alone fix (D, alpha) exactly for any caliber whose
    table is fast enough; the N-wave term mainly pins the miss distance and is
    what separates calibers. J <~ 1 means the caliber explains everything; a
    few or more means it does not. A coarse grid finds the basins, Nelder-Mead
    refines the best few. Returns the local minima found, best first."""
    v0, lam = bl["v0"], bl["lam"]
    sig_th = math.radians(rcfg["direction_sigma_deg"])
    sig_dt = rcfg["dt_sigma_ms"] * 1e-3
    sig_t = rcfg["t_n_sigma_log"]
    sig_b = rcfg["miss_sigma_m"]
    use_t = (t_n_meas is not None and np.isfinite(t_n_meas) and t_n_meas > 0 and sig_t > 0)
    sig_p = rcfg.get("dp_sigma_log", 0.0)
    use_p = (dp_meas is not None and np.isfinite(dp_meas) and dp_meas > 0 and sig_p > 0)
    mm = rcfg["mach_min"]

    def t_term(M, b):
        with np.errstate(divide="ignore", invalid="ignore"):
            T = WHITHAM_K_T * M * bl["d"] * (np.maximum(b, 1e-6) / bl["L"]) ** 0.25 / (
                c * np.sqrt(np.maximum(M ** 2 - 1.0, 1e-12)) ** 0.75)
            return (np.log(T / t_n_meas) / sig_t) ** 2, T

    def p_term(M, b):
        # Whitham peak overpressure; with T_N it fixes the miss distance AND the bullet size
        with np.errstate(divide="ignore", invalid="ignore"):
            P = WHITHAM_K_P * P_ATM * np.maximum(M ** 2 - 1.0, 1e-12) ** 0.125 * bl["d"] / (
                np.maximum(b, 1e-6) ** 0.75 * bl["L"] ** 0.25)
            return (np.log(P / dp_meas) / sig_p) ** 2, P

    Ds = np.geomspace(rcfg["min_range_m"], rcfg["max_range_m"], rcfg["n_grid"])
    als = np.radians(np.geomspace(0.02, 80.0, 160))
    f = forward_model(Ds[:, None], als[None, :], c, v0, lam, mm)
    J = ((f["theta"] - theta_meas) / sig_th) ** 2 + ((f["dt"] - dt_meas) / sig_dt) ** 2
    if use_t:
        J = J + t_term(f["v"] / c, f["miss"])[0]
    if use_p:
        J = J + p_term(f["v"] / c, f["miss"])[0]
    if sig_b > 0:
        J = J + (f["miss"] / sig_b) ** 2
    J = np.where(f["valid"], J, np.inf)
    prof = J.min(axis=1)
    best_a = J.argmin(axis=1)
    if not np.isfinite(prof).any():
        return []

    # local minima of the range profile, best few -> refine in (ln D, ln alpha)
    mins = [i for i in range(len(Ds)) if np.isfinite(prof[i])
            and (i == 0 or prof[i] <= prof[i - 1]) and (i == len(Ds) - 1 or prof[i] <= prof[i + 1])]
    mins = sorted(mins, key=lambda i: prof[i])[:max_starts]

    def terms(x):
        r = _forward_scalar(math.exp(x[0]), math.exp(x[1]), c, v0, lam, mm)
        if r is None:
            return None, None
        j = [((r["theta"] - theta_meas) / sig_th) ** 2, ((r["dt"] - dt_meas) / sig_dt) ** 2, 0.0,
             (r["miss"] / sig_b) ** 2 if sig_b > 0 else 0.0, 0.0]
        if use_t:
            j[2] = float(t_term(r["v"] / c, r["miss"])[0])
        if use_p:
            j[4] = float(p_term(r["v"] / c, r["miss"])[0])
        return j, r

    def cost(x):
        j, _ = terms(x)
        return 1e12 if j is None else sum(j)

    sols = []
    for i in mins:
        res = minimize(cost, [math.log(Ds[i]), math.log(als[best_a[i]])], method="Nelder-Mead",
                       options=dict(xatol=1e-7, fatol=1e-10, maxiter=800))
        j, r = terms(res.x)
        D, al = math.exp(res.x[0]), math.exp(res.x[1])
        if j is None or any(abs(D - s["D"]) < 1e-3 * D for s in sols):
            continue
        sols.append(dict(D=D, alpha=al, J=float(sum(j)), J_terms=[float(x) for x in j],
                         miss=r["miss"], s_emit=r["s"], rho=r["rho"], v_emit=r["v"], mach_emit=r["v"] / c,
                         theta_residual_deg=math.degrees(r["theta"] - theta_meas),
                         dt_residual_ms=(r["dt"] - dt_meas) * 1e3))
    # Best first. Among solutions that all fit (J < 1), prefer the smaller miss
    # distance (shots aimed at or near the array are the common case).
    sols.sort(key=lambda s: (s["J"] >= 1.0, s["miss"] if s["J"] < 1.0 else s["J"]))
    return sols


def whitham_duration(M, b, d, L, c):
    """Predicted N-wave duration (s)."""
    if M <= 1.0 or b <= 0:
        return np.nan
    return WHITHAM_K_T * M * d * (b / L) ** 0.25 / (c * (M ** 2 - 1.0) ** 0.375)


def t_n_plausible_range(ballistics, c, miss_range=(0.5, 100.0), mach_range=(1.05, 3.0), margin=1.4):
    """Widest band of N-wave durations any caliber in the table can produce for
    miss distances and local Mach numbers in the given ranges, widened by margin.
    A measured T_N outside it is a Stage 1 measurement failure, not physics."""
    Ms = np.linspace(mach_range[0], mach_range[1], 50)
    Ts = [whitham_duration(M, b, bl["d"], bl["L"], c)
          for bl in ballistics.values() for b in miss_range for M in Ms]
    return min(Ts) / margin, max(Ts) * margin


def trajectory_vectors(sol, n_hat, k_hat):
    """3D shooter position, emission point, bullet direction and closest point of
    approach, in the plane spanned by the shooter direction and the SW direction."""
    src = -k_hat                                            # apparent SW source direction
    e = src - np.dot(src, n_hat) * n_hat
    if np.linalg.norm(e) < 1e-9:
        return {}
    e /= np.linalg.norm(e)
    P = sol["D"] * n_hat
    u = -np.cos(sol["alpha"]) * n_hat + np.sin(sol["alpha"]) * e
    return dict(shooter_position_m=P.tolist(),
                emission_point_m=(P + sol["s_emit"] * u).tolist(),
                bullet_direction=u.tolist(),
                closest_point_of_approach_m=(P + sol["D"] * np.cos(sol["alpha"]) * u).tolist())


def locate_drag(rcfg, dt_meas, theta, c, n_hat, k_hat, t_n_meas, forced_caliber=None, dp_meas=None):
    """Fit every caliber in the ballistic table, pick one, return a result dict."""
    t_ok = (t_n_meas is not None and np.isfinite(t_n_meas) and t_n_meas > 0
            and rcfg["t_n_sigma_log"] > 0)
    t_note = None
    if t_ok:
        lo, hi = t_n_plausible_range(rcfg["ballistics"], c)
        if not (lo <= t_n_meas <= hi):
            t_ok = False
            t_note = (f"measured T_N {t_n_meas * 1e6:.0f} us is outside the physically possible "
                      f"{lo * 1e6:.0f}-{hi * 1e6:.0f} us for this table -- ignored (check Stage 1)")
            print(f"[Det Stage 3] WARNING: {t_note}")
    candidates = {}
    for cal, bl in rcfg["ballistics"].items():
        sols = fit_range_drag(dt_meas, theta, c, bl, rcfg, t_n_meas if t_ok else None, dp_meas=dp_meas)
        if not sols:
            candidates[cal] = dict(valid=False)
            print(f"    {cal}: no supersonic geometry fits")
            continue
        s0 = sols[0]
        t_pred = whitham_duration(s0["mach_emit"], s0["miss"], bl["d"], bl["L"], c)
        p_pred = (WHITHAM_K_P * P_ATM * (s0["mach_emit"] ** 2 - 1.0) ** 0.125 * bl["d"]
                  / (max(s0["miss"], 1e-6) ** 0.75 * bl["L"] ** 0.25)) if s0["mach_emit"] > 1 else float("nan")
        candidates[cal] = dict(valid=True, range_m=s0["D"], miss_distance_m=s0["miss"],
                               mach_at_emission=s0["mach_emit"], speed_at_emission_mps=s0["v_emit"],
                               fit_cost=s0["J"],
                               fit_cost_terms=dict(sw_angle=s0["J_terms"][0], dt=s0["J_terms"][1],
                                                   n_wave_duration=s0["J_terms"][2],
                                                   miss_prior=s0["J_terms"][3],
                                                   sw_amplitude=s0["J_terms"][4]),
                               sw_angle_residual_deg=s0["theta_residual_deg"],
                               dt_residual_ms=s0["dt_residual_ms"],
                               t_n_predicted_s=t_pred, t_n_measured_s=t_n_meas if t_ok else None,
                               dp_predicted_pa=p_pred, dp_measured_pa=dp_meas,
                               muzzle_velocity_assumed_mps=bl["v0"], decay_length_assumed_m=bl["lam"],
                               alternative_solutions=[dict(range_m=s["D"], miss_distance_m=s["miss"],
                                                           fit_cost=s["J"]) for s in sols[1:4]],
                               _sol=s0)
        print(f"    {cal:10s} range {s0['D']:7.1f} m  miss {s0['miss']:5.1f} m  Mach {s0['mach_emit']:.2f}  "
              f"cost {s0['J']:8.2f}  (angle {s0['theta_residual_deg']:+.2f} deg, Dt {s0['dt_residual_ms']:+.2f} ms, "
              f"T_N {t_pred * 1e6:.0f}" + (f"/{t_n_meas * 1e6:.0f} us" if t_ok else " us pred")
              + (f", dP {p_pred:.0f}/{dp_meas:.0f} Pa)" if dp_meas else ")")
              + (f"  +{len(sols) - 1} alt." if len(sols) > 1 else ""))

    valid = {k: v for k, v in candidates.items() if v["valid"]}
    if not valid:
        return dict(range_model="drag", range_estimate_m=None, caliber_candidates=candidates)

    fitting = [k for k in valid if valid[k]["fit_cost"] < 1.0]
    if forced_caliber:
        if forced_caliber not in valid:
            raise SystemExit(f"[Det Stage 3] --caliber {forced_caliber}: no valid solution "
                             f"(valid: {', '.join(valid)})")
        chosen, how = forced_caliber, "forced with --caliber"
    elif not t_ok and rcfg["default_caliber"] in fitting:
        chosen, how = rcfg["default_caliber"], "default_caliber (fits; no N-wave duration used)"
    else:
        chosen = min(valid, key=lambda k: valid[k]["fit_cost"])
        how = ("lowest cost: SW angle + delay" + (" + N-wave duration" if t_ok else "")
               + (" + SW amplitude" if dp_meas and rcfg.get("dp_sigma_log", 0) > 0 else "")
               + ("" if t_ok else " (no N-wave duration used)"))
        if not t_ok and len(fitting) > 1:
            how += f"; ambiguous, also fits: {', '.join(k for k in fitting if k != chosen)}"

    sol = valid[chosen]["_sol"]
    for v in candidates.values():
        v.pop("_sol", None)
    if sol["J"] > 4.0:
        print(f"[Det Stage 3] WARNING: best fit cost {sol['J']:.1f} > 4 -- the measurements are not "
              f"consistent with this caliber's ballistics; treat the range as rough.")
    if t_note:
        how += "; " + t_note
    # Reported calibre (project rule): category (small, medium, big, RPG, other; mandatory) and
    # class (bore number: 5.56, 7.62, 8.6, 12.7, ...; optional, [output] report_caliber_class).
    # The variant (best-fitting entry, e.g. 7.62x54R) is not reported as a result; it is only
    # named as range_ballistics, the ballistics the range was computed with. Per class and
    # category: the lowest fit cost of its variants, and per class / category the lowest fit cost of its
    # members. Cartridges of one bore (7.62x39 / x51 / x54R) are hard to tell apart acoustically.
    cls_of = {k: rcfg["ballistics"][k].get("cls", k) for k in candidates}
    cat_of = {k: rcfg["ballistics"][k].get("cat", "other") for k in candidates}
    class_costs, cat_costs = {}, {}
    for k, v in valid.items():
        class_costs[cls_of[k]] = min(class_costs.get(cls_of[k], np.inf), v["fit_cost"])
        cat_costs[cat_of[k]] = min(cat_costs.get(cat_of[k], np.inf), v["fit_cost"])
    for k, v in candidates.items():
        v["caliber_class"], v["caliber_category"] = cls_of[k], cat_of[k]
    order = lambda d: ", ".join(f"{c} {j:.2f}" for c, j in sorted(d.items(), key=lambda kv: kv[1]))
    print(f"[Det Stage 3] Calibre category {cat_of[chosen]}, class {cls_of[chosen]} "
          f"(range computed with the {chosen} ballistics)")
    print(f"[Det Stage 3]   cost per class: {order(class_costs)}")
    print(f"[Det Stage 3]   cost per category: {order(cat_costs)}")
    out = dict(range_model="drag", caliber_category=cat_of[chosen], caliber_class=cls_of[chosen],
               range_ballistics=chosen,
               caliber_class_costs={c: float(j) for c, j in class_costs.items()},
               caliber_category_costs={c: float(j) for c, j in cat_costs.items()},
               caliber_selection=how,
               range_estimate_m=sol["D"], miss_distance_m=sol["miss"],
               mach_at_emission=sol["mach_emit"], speed_at_emission_mps=sol["v_emit"],
               along_track_to_cpa_m=sol["D"] * np.cos(sol["alpha"]),
               fit_cost=sol["J"], fit_consistent=bool(sol["J"] <= 4.0),
               caliber_candidates=candidates)
    out.update(trajectory_vectors(sol, n_hat, k_hat))
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

# Several shots: a shot is flagged 'bearing_outlier' when its bearing is more than this many
# degrees from the median bearing of all shots in the file (needs at least 3 shots).
BEARING_OUTLIER_DEG = 5.0


def main():
    parser = argparse.ArgumentParser(
        description="Detection Stage 3: DOA + drag-aware range -> 3D shooter location.")
    parser.add_argument("tdoa_json", help="Path to Stage 2's *_tdoa.json")
    parser.add_argument("config", nargs="?", default="det_config.ini",
                        help="Path to detection .ini config (default: det_config.ini)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output basename (overrides output.basename in config)")
    parser.add_argument("--range-model", choices=["drag", "constant"], default=None,
                        help="Range model (default: [range] model in the .ini, else drag)")
    parser.add_argument("--temperature-c", type=float, default=None,
                        help="Air temperature (deg C) for c = 331.3 + 0.606 T "
                             "(overrides [physics] temperature_c / speed_of_sound)")
    parser.add_argument("--caliber", default=None,
                        help="Drag model: use this [bullet_reference_library] caliber instead of choosing by T_N")
    parser.add_argument("--assumed-miss-distance", type=float, default=None,
                        help="Constant model only: miss distance b (m) instead of the "
                             "amplitude-derived estimate.")
    add_mic_config_argument(parser)
    args = parser.parse_args()

    with open(args.tdoa_json) as f:
        tdoa_data = json.load(f)

    cfg = load_det_config(args.config, args.mic_config or upstream_mic_config(tdoa_data))
    rcfg = load_range_config(args.config, cfg["bullet_reference_library"])
    if args.output:
        cfg["output"]["basename"] = args.output
    model = args.range_model or rcfg["model"]
    if model not in ("drag", "constant"):
        raise SystemExit(f"[Det Stage 3] unknown range model '{model}' (drag or constant)")

    geom = cfg["geometry"]
    mic_pos = build_mic_positions(cfg)
    side = geom["planar_source_side"]
    if tdoa_data.get("n_channels", len(mic_pos)) != len(mic_pos):
        raise SystemExit(f"[Det Stage 3] ERROR: TDOA file has {tdoa_data['n_channels']} channels but "
                         f"mic_config {geom['mic_config']} has {len(mic_pos)} mics.")
    print(f"[Det Stage 3] mic_config = {geom['mic_config']} ({len(mic_pos)} mics)")

    t_c = args.temperature_c if args.temperature_c is not None else rcfg["temperature_c"]
    if t_c is not None:
        c = speed_of_sound_from_temperature(t_c)
        c_source = f"temperature {t_c:.1f} C"
    else:
        c = cfg["physics"]["speed_of_sound"]
        c_source = "[physics] speed_of_sound"
    print(f"[Det Stage 3] c = {c:.1f} m/s ({c_source})")
    # The library's lam is for 15 C (ICAO). Drag scales with air density, which at constant
    # pressure goes as 1/T: lam_eff = lam x (T + 273.15) / 288.15.
    if t_c is not None and rcfg["lam_temperature_scaling"]:
        f_lam = (t_c + 273.15) / 288.15
        for bl in rcfg["ballistics"].values():
            bl["lam"] = bl["lam"] * f_lam
        print(f"[Det Stage 3] bullet decay lengths scaled x{f_lam:.3f} for the air density at {t_c:.1f} C")
    result = dict(
        generated_at_utc=datetime.now(timezone.utc).isoformat(),
        stage="shooter_location",
        source_tdoa_json=args.tdoa_json,
        mic_config=geom["mic_config"],
        speed_of_sound_mps=c,
        speed_of_sound_source=c_source,
    )

    def locate_one(tdoa_data, result):
        """Bearing, SW direction and range of one shot (fills result)."""
        ref_ch = tdoa_data["reference_channel"]

        sw_available = tdoa_data.get("sw_available", bool(tdoa_data.get("tdoa_sw_s")))
        result["event_class"] = tdoa_data.get("event_class")
        k_hat_sw = None
        if sw_available:
            k_hat_sw, k_norm_sw = solve_doa(mic_pos, tdoa_data["tdoa_sw_s"], c, ref_ch, side)
            result["k_hat_sw"] = k_hat_sw.tolist()
            result["k_hat_sw_norm_check"] = k_norm_sw  # should be ~1.0 for consistent TDOA
            print(f"[Det Stage 3] k_hat_sw = {k_hat_sw}  (|k| = {k_norm_sw:.4f}, should be ~1.0)")
        else:
            print("[Det Stage 3] No shockwave: blast-only event -- bearing from the muzzle blast, "
                  "no bullet direction, Mach or range.")

        if not tdoa_data["mb_available"]:
            print("[Det Stage 3] No reliable muzzle-blast TDOA available -- bearing/range/position "
                  "cannot be computed, only the shockwave direction (k_hat_sw) above.")
            result["bearing_available"] = False
            if k_hat_sw is not None:
                # A close impulsive TEST sound (clap, balloon) also looks like a shockwave. For such a
                # source -- not a bullet -- the direction to it is simply -k_hat_sw. Informational only.
                src = -k_hat_sw
                az_i = float(np.degrees(np.arctan2(src[1], src[0])))
                el_i = float(np.degrees(np.arcsin(np.clip(src[2], -1, 1))))
                result["impulse_source_direction_if_not_a_bullet"] = dict(azimuth_deg=az_i, elevation_deg=el_i)
                print(f"[Det Stage 3] If this was a test impulse (clap, balloon), not a bullet: its source "
                      f"is at azimuth={az_i:.2f} deg, elevation={el_i:.2f} deg.")
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

            # theta: angle between the shooter direction and the apparent shockwave
            # source direction. 1/cos(theta) is the angle-only Mach (Method A), which
            # is biased low by the miss angle; the drag solve corrects for that.
            dt_sw_to_mb = tdoa_data.get("dt_sw_to_mb_s")
            result["dt_sw_to_mb_s"] = dt_sw_to_mb
            theta, M_A = None, None
            if k_hat_sw is not None:
                cos_theta = float(np.clip(np.dot(k_hat_sw, -n_hat_shooter), -1.0, 1.0))
                theta = float(np.arccos(cos_theta))
                M_A = 1.0 / cos_theta if cos_theta > 0 else None
                result["sw_angle_from_shooter_deg"] = float(np.degrees(theta))
            result["mach_estimate_angle_only"] = M_A

            if k_hat_sw is None or dt_sw_to_mb is None or M_A is None or M_A <= 1.01:
                reason = ("no shockwave (blast only)" if k_hat_sw is None
                          else "no Dt available" if dt_sw_to_mb is None
                          else "SW direction not consistent with a supersonic bullet")
                print(f"[Det Stage 3] Skipping range solve ({reason}).")
                result["range_estimate_m"] = None
            elif model == "drag":
                print(f"[Det Stage 3] Drag-aware range solve: Dt = {dt_sw_to_mb * 1e3:.2f} ms, "
                      f"SW {np.degrees(theta):.2f} deg off the shooter direction "
                      f"(angle-only Mach {M_A:.3f})")
                # median over the channels (more robust than one channel), else the reference's
                t_n = tdoa_data.get("t_n_sw_median") or tdoa_data.get("t_n_sw_ref_channel")
                # measured T_N includes the sensor's response (rounded shocks): remove it before
                # comparing with Whitham's T_N
                if t_n and rcfg["t_n_offset_s"]:
                    print(f"[Det Stage 3] T_N measured {t_n * 1e6:.0f} us, minus sensor offset "
                          f"{rcfg['t_n_offset_s'] * 1e6:.0f} us = {(t_n - rcfg['t_n_offset_s']) * 1e6:.0f} us")
                    t_n = t_n - rcfg["t_n_offset_s"]
                # SW peak in Pa (same calibration as Stage 4: assumed_peak_pa_at_full_scale)
                dp_pa = None
                if tdoa_data.get("dp_sw_ref_channel") and tdoa_data.get("normalization_peak"):
                    dp_pa = (abs(tdoa_data["dp_sw_ref_channel"]) * tdoa_data["normalization_peak"]
                             * cfg["calibration"]["assumed_peak_pa_at_full_scale"])
                    print(f"[Det Stage 3] SW peak {dp_pa:.0f} Pa (assumed full scale "
                          f"{cfg['calibration']['assumed_peak_pa_at_full_scale']:.0f} Pa)")
                out = locate_drag(rcfg, dt_sw_to_mb, theta, c, n_hat_shooter, k_hat_sw, t_n, args.caliber,
                                  dp_meas=dp_pa)
                result.update(out)
                if out.get("range_estimate_m") is None:
                    print("[Det Stage 3] No caliber in [bullet_reference_library] fits -- bearing-only result reported.")
                else:
                    print(f"[Det Stage 3] Calibre category {out['caliber_category']}, class "
                          f"{out['caliber_class']} ({out['caliber_selection']}): "
                          f"range = {out['range_estimate_m']:.1f} m, miss distance = "
                          f"{out['miss_distance_m']:.1f} m, bullet speed at the array = "
                          f"{out['speed_at_emission_mps']:.0f} m/s (Mach {out['mach_at_emission']:.2f})")
            else:
                # Legacy constant-speed solve. Miss distance from the SW amplitude by inverting
                # Whitham's dP law for the first [bullet_reference_library] caliber at the
                # angle-only Mach (the library no longer has dP0_sw / b0_sw reference points).
                dp_sw_norm = tdoa_data.get("dp_sw_ref_channel")
                norm_peak = tdoa_data.get("normalization_peak")
                dp_sw = (abs(dp_sw_norm) * (norm_peak if norm_peak else 1.0)
                         * cfg["calibration"]["assumed_peak_pa_at_full_scale"]
                         if dp_sw_norm is not None else None)
                if dp_sw is not None and not norm_peak:
                    print("[Det Stage 3] WARNING: no normalization_peak in the TDOA file -- SW amplitude "
                          "taken as normalized units x full scale")
                default_caliber = next(iter(cfg["bullet_reference_library"]))
                bl = cfg["bullet_reference_library"][default_caliber]
                p_atm = cfg["physics"].get("atmospheric_pressure_pa", 101325.0)
                if args.assumed_miss_distance is not None:
                    b_est, b_source = args.assumed_miss_distance, "user-provided --assumed-miss-distance"
                elif dp_sw is not None and dp_sw > 0:
                    b_est = (WHITHAM_K_P * p_atm * (M_A ** 2 - 1.0) ** 0.125 * bl["d"]
                             / (bl["L"] ** 0.25 * dp_sw)) ** (4.0 / 3.0)
                    b_source = f"amplitude-derived (Whitham), assuming default caliber {default_caliber}"
                else:
                    b_est, b_source = 10.0, "no amplitude measurement available -- arbitrary 10m fallback"
                try:
                    r_est = pilar_v_exact(dt_sw_to_mb, b_est, M_A, c)
                    r_3d = float(np.sqrt(r_est ** 2 + b_est ** 2))
                    result.update(range_model="constant", range_estimate_m=r_3d,
                                  along_track_to_cpa_m=r_est, miss_distance_m=b_est,
                                  miss_distance_source=b_source,
                                  shooter_position_m=(r_3d * n_hat_shooter).tolist())
                    print(f"[Det Stage 3] Constant-speed range = {r_3d:.2f} m "
                          f"(Mach {M_A:.3f} angle-only, b = {b_est:.2f} m, {b_source})")
                except (ValueError, RuntimeError) as e:
                    print(f"[Det Stage 3] Range solve failed ({e}) -- bearing-only result reported.")
                    result["range_estimate_m"] = None

        # Calibre (project rule): category always present, class only if configured
        if result.get("caliber_category") is None:
            result["caliber_category"] = None
            result["caliber_category_reason"] = (
                "no shockwave: muzzle blast only, the calibre cannot be determined"
                if not sw_available else
                "no muzzle blast: no SW -> MB delay, the calibre cannot be determined"
                if not result.get("bearing_available") or not tdoa_data.get("dt_sw_to_mb_s") else
                "no calibre fits the shockwave and the delay"
                if model == "drag" else "range model 'constant' does not determine the calibre")
        if not cfg["output"]["report_caliber_class"]:
            result.pop("caliber_class", None)
            result.pop("caliber_class_costs", None)

        # Project limits: 50-2000 m. Closer than friendly_range_m = own position (not a
        # non-friendly shot); beyond max_report_range_m = outside the project's coverage.
        rng = result.get("range_estimate_m")
        result["friendly"] = bool(rng is not None and rng < rcfg["friendly_range_m"])
        result["beyond_max_range"] = bool(rng is not None and rng > rcfg["max_report_range_m"])
        if result["friendly"]:
            print(f"[Det Stage 3] Range {rng:.0f} m < {rcfg['friendly_range_m']:.0f} m: flagged 'friendly' "
                  f"(not reported as a non-friendly shot)")
        if result["beyond_max_range"]:
            print(f"[Det Stage 3] Range {rng:.0f} m > {rcfg['max_report_range_m']:.0f} m: flagged "
                  f"'beyond_max_range' (outside the project's coverage)")
        return result

    locate_one(tdoa_data, result)
    shots = tdoa_data.get("shots") or []
    if len(shots) > 1:
        result["n_shots"] = len(shots)
        result["shots"] = []
        for sh in shots:
            print(f"\n[Det Stage 3] ---- shot {sh.get('index')} ({sh.get('event_class')}) ----")
            td = dict(tdoa_data)
            td.update({k: v for k, v in sh.items()})
            r = locate_one(td, dict(index=sh.get("index"), t_sw_s=sh.get("t_sw_s"), t_mb_s=sh.get("t_mb_s")))
            result["shots"].append(r)
        # The rounds of one burst come from one shooter. A shot whose bearing is far from the
        # median bearing of the burst was probably timed on an overlapping event (for example
        # the SW or an echo of another round inside its MB window): flagged, not changed.
        dirs = [np.asarray(r["n_hat_shooter"], float) for r in result["shots"]
                if r.get("bearing_available") and r.get("n_hat_shooter") is not None]
        if len(dirs) >= 3:
            med = np.median(np.array(dirs), axis=0)
            med /= np.linalg.norm(med)
            for r in result["shots"]:
                if r.get("bearing_available") and r.get("n_hat_shooter") is not None:
                    ang = float(np.degrees(np.arccos(np.clip(np.dot(med, r["n_hat_shooter"]), -1, 1))))
                    r["angle_to_burst_bearing_deg"] = ang
                    r["bearing_outlier"] = ang > BEARING_OUTLIER_DEG
            n_out = sum(bool(r.get("bearing_outlier")) for r in result["shots"])
            if n_out:
                print(f"[Det Stage 3] WARNING: {n_out} shot(s) more than {BEARING_OUTLIER_DEG:.0f} deg "
                      f"from the burst's median bearing (flag 'bearing_outlier')")
            for k in ("angle_to_burst_bearing_deg", "bearing_outlier"):
                if k in result["shots"][0]:
                    result[k] = result["shots"][0][k]
        print(f"\n[Det Stage 3] {len(shots)} shots located (top-level result = shot 1; all in 'shots')")

    out_json = f"{cfg['output']['basename']}_location.json"
    with open(out_json, "w") as f:
        json.dump(result, f, indent=2)
    print(f"[Det Stage 3] Wrote {out_json}")


if __name__ == "__main__":
    main()
