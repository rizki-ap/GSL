#!/usr/bin/env python3
"""Validate the drag-aware range solver in gs_det_shooter_locator.py against the
generator's own shock model (gs_gen_physic.shock_geometry_drag), at the level of
the three measured quantities: n_hat (to shooter), k_hat (SW propagation), Dt.
Test-only use of the generator -- the detector itself never imports it."""
import numpy as np
from gs_gen_physic import shock_geometry_drag, BULLET_LIBRARY
from gs_det_shooter_locator import fit_range_drag, pilar_v_exact, RANGE_DEFAULTS
RCFG = {k: (float(v) if k not in ('model', 'default_caliber') else v) for k, v in RANGE_DEFAULTS.items()}
RCFG['n_grid'] = int(RCFG['n_grid'])

C = 347.0
rng = np.random.default_rng(1)

def unit(v): return v / np.linalg.norm(v)

def rot(v, ang):   # rotate v by ang (rad) about a random axis perpendicular to it
    a = unit(np.cross(v, rng.standard_normal(3)))
    return unit(v * np.cos(ang) + np.cross(a, v) * np.sin(ang))

def scenario(cal, D, b, az_deg, el_deg=1.0):
    g = BULLET_LIBRARY[cal]
    az, el = np.radians([az_deg, el_deg])
    n = np.array([np.cos(el) * np.cos(az), -np.cos(el) * np.sin(az), np.sin(el)])
    P = D * n
    side = unit(np.cross(n, [0, 0, 1]))
    v = unit(b * side - P)                      # aim so the path passes b metres beside the array
    rec = shock_geometry_drag(np.zeros(3), P, v, g["v0"], g["lam"], C)
    if not rec["in_cone"]:
        return None
    k = unit(-rec["x_emit"])                    # SW propagation direction at the array
    dt = D / C - rec["t_arrive"]
    true_b = np.linalg.norm(np.cross(-P, v))
    return dict(n=n, k=k, dt=dt, v0=g["v0"], lam=g["lam"], b=true_b, M=rec["M_emit"])

def solve(sc, v0, lam, n_err_deg=0.0, k_err_deg=0.0):
    n = rot(sc["n"], np.radians(n_err_deg)) if n_err_deg else sc["n"]
    k = rot(sc["k"], np.radians(k_err_deg)) if k_err_deg else sc["k"]
    theta = np.arccos(np.clip(np.dot(k, -n), -1, 1))
    sols = fit_range_drag(sc["dt"], theta, C, dict(v0=v0, lam=lam, d=0.0057, L=0.023), RCFG)
    return sols[0] if sols else None, theta

def const_range(sc):
    theta = np.arccos(np.clip(np.dot(sc["k"], -sc["n"]), -1, 1))
    MA = 1 / np.cos(theta)
    r = pilar_v_exact(sc["dt"], sc["b"], MA, C)       # given the TRUE miss distance
    return np.hypot(r, sc["b"])

print(f"{'cal':10s} {'D':>5s} {'b':>4s} | {'M_emit':>6s} | {'drag, matched':>13s} {'miss est':>8s} | "
      f"{'lam -20%':>8s} {'lam +20%':>8s} {'v0 -3%':>7s} | {'1deg noise':>10s} | {'constant*':>9s}")
for cal in ["5.56_NATO", "7.62_NATO", "12.7_NATO"]:
    for D in [100, 200, 400, 600, 800]:
        for b in [5, 20]:
            sc = scenario(cal, D, b, 30)
            if sc is None:
                print(f"{cal:10s} {D:5d} {b:4d} | no shockwave at the array"); continue
            def err(s): return f"{100 * (s['D'] - D) / D:+7.1f}%" if s else "   n/a "
            m, _ = solve(sc, sc["v0"], sc["lam"])
            lo, _ = solve(sc, sc["v0"], sc["lam"] * 0.8)
            hi, _ = solve(sc, sc["v0"], sc["lam"] * 1.2)
            vs, _ = solve(sc, sc["v0"] * 0.97, sc["lam"])
            noisy = [solve(sc, sc["v0"], sc["lam"], 1.0 / np.sqrt(2), 1.0 / np.sqrt(2))[0] for _ in range(20)]
            ne = [100 * (s["D"] - D) / D for s in noisy if s]
            noise = f"{np.sqrt(np.mean(np.square(ne))):5.1f}% rms" if ne else "n/a"
            cr = const_range(sc)
            print(f"{cal:10s} {D:5d} {b:4d} | {sc['M']:6.2f} | {err(m):>13s} {m['miss'] if m else float('nan'):7.1f}m | "
                  f"{err(lo):>8s} {err(hi):>8s} {err(vs):>7s} | {noise:>10s} | {100 * (cr - D) / D:+8.1f}%")
print("\n* constant-speed PILAR-V solve, angle-only Mach, given the TRUE miss distance (old Stage 3, best case)")
