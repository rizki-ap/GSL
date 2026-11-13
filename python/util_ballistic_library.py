#!/usr/bin/env python3
"""Bullet library from ammunition data: muzzle speed v0 and speed-decay length lam for the
exponential drag model v(s) = v0 exp(-s / lam) used by the generator (gs_gen_physic.py
BULLET_LIBRARY) and the detector (det_config.ini [bullet_reference_library]).

Input: ammo_data.csv (one row per calibre): bullet mass, muzzle velocity, ballistic coefficient
with its drag model (G1 or G7, lb/in^2), diameter, length, and where the numbers come from.

Method:
  1. Flat-fire point-mass solution of the speed along the path with the standard G1 / G7 drag
     tables (embedded below) in the ICAO standard atmosphere (15 C, 101.325 kPa):
        dv/ds = -rho pi Cd_G(v / c) v / (8 BC),   BC in kg/m^2 = BC[lb/in^2] x 703.07
     Gravity is ignored (it changes the speed along the path very little).
  2. lam fitted to the TIME OF FLIGHT, t(s) = lam (exp(s / lam) - 1) / v0, over the part of
     the path where the bullet is still above Mach mach_end (1.1 = the generator's shockwave
     threshold), up to the project's 2000 m. Time of flight is what the SW -> MB delay depends on.
  lam is for the standard atmosphere: drag scales with air density, so at temperature T and
  pressure p, lam_eff = lam x (T + 273.15) / 288.15 x 101.325 / p (applied by both blocks).

Usage:
  python util_ballistic_library.py                     # table of results
  python util_ballistic_library.py --ini               # lines for det_config.ini
  python util_ballistic_library.py --temp-c 27         # same, another temperature (checks)
  python util_ballistic_library.py --csv my_ammo.csv   # your own ammunition data
"""
import argparse
import csv
import math
import os

import numpy as np
from scipy.optimize import minimize_scalar

# Standard drag tables (Mach, Cd), as published for the G1 and G7 reference projectiles.
TABLES = {}
TABLES["G1"] = [
    (0, 0.2629), (0.05, 0.2558), (0.1, 0.2487), (0.15, 0.2413), (0.2, 0.2344), (0.25, 0.2278),
    (0.3, 0.2214), (0.35, 0.2155), (0.4, 0.2104), (0.45, 0.2061), (0.5, 0.2032), (0.55, 0.202),
    (0.6, 0.2034), (0.7, 0.2165), (0.725, 0.223), (0.75, 0.2313), (0.775, 0.2417), (0.8, 0.2546),
    (0.825, 0.2706), (0.85, 0.2901), (0.875, 0.3136), (0.9, 0.3415), (0.925, 0.3734), (0.95, 0.4084),
    (0.975, 0.4448), (1, 0.4805), (1.025, 0.5136), (1.05, 0.5427), (1.075, 0.5677), (1.1, 0.5883),
    (1.125, 0.6053), (1.15, 0.6191), (1.2, 0.6393), (1.25, 0.6518), (1.3, 0.6589), (1.35, 0.6621),
    (1.4, 0.6625), (1.45, 0.6607), (1.5, 0.6573), (1.55, 0.6528), (1.6, 0.6474), (1.65, 0.6413),
    (1.7, 0.6347), (1.75, 0.628), (1.8, 0.621), (1.85, 0.6141), (1.9, 0.6072), (1.95, 0.6003),
    (2, 0.5934), (2.05, 0.5867), (2.1, 0.5804), (2.15, 0.5743), (2.2, 0.5685), (2.25, 0.563),
    (2.3, 0.5577), (2.35, 0.5527), (2.4, 0.5481), (2.45, 0.5438), (2.5, 0.5397), (2.6, 0.5325),
    (2.7, 0.5264), (2.8, 0.5211), (2.9, 0.5168), (3, 0.5133), (3.1, 0.5105), (3.2, 0.5084),
    (3.3, 0.5067), (3.4, 0.5054), (3.5, 0.504), (3.6, 0.503), (3.7, 0.5022), (3.8, 0.5016),
    (3.9, 0.501), (4, 0.5006), (4.2, 0.4998), (4.4, 0.4995), (4.6, 0.4992), (4.8, 0.499),
    (5, 0.4988),
]
TABLES["G7"] = [
    (0, 0.1198), (0.05, 0.1197), (0.1, 0.1196), (0.15, 0.1194), (0.2, 0.1193), (0.25, 0.1194),
    (0.3, 0.1194), (0.35, 0.1194), (0.4, 0.1193), (0.45, 0.1193), (0.5, 0.1194), (0.55, 0.1193),
    (0.6, 0.1194), (0.65, 0.1197), (0.7, 0.1202), (0.725, 0.1207), (0.75, 0.1215), (0.775, 0.1226),
    (0.8, 0.1242), (0.825, 0.1266), (0.85, 0.1306), (0.875, 0.1368), (0.9, 0.1464), (0.925, 0.166),
    (0.95, 0.2054), (0.975, 0.2993), (1, 0.3803), (1.025, 0.4015), (1.05, 0.4043), (1.075, 0.4034),
    (1.1, 0.4014), (1.125, 0.3987), (1.15, 0.3955), (1.2, 0.3884), (1.25, 0.381), (1.3, 0.3732),
    (1.35, 0.3657), (1.4, 0.358), (1.5, 0.344), (1.55, 0.3376), (1.6, 0.3315), (1.65, 0.326),
    (1.7, 0.3209), (1.75, 0.316), (1.8, 0.3117), (1.85, 0.3078), (1.9, 0.3042), (1.95, 0.301),
    (2, 0.298), (2.05, 0.2951), (2.1, 0.2922), (2.15, 0.2892), (2.2, 0.2864), (2.25, 0.2835),
    (2.3, 0.2807), (2.35, 0.2779), (2.4, 0.2752), (2.45, 0.2725), (2.5, 0.2697), (2.55, 0.267),
    (2.6, 0.2643), (2.65, 0.2615), (2.7, 0.2588), (2.75, 0.2561), (2.8, 0.2533), (2.85, 0.2506),
    (2.9, 0.2479), (2.95, 0.2451), (3, 0.2424), (3.1, 0.2368), (3.2, 0.2313), (3.3, 0.2258),
    (3.4, 0.2205), (3.5, 0.2154), (3.6, 0.2106), (3.7, 0.206), (3.8, 0.2017), (3.9, 0.1975),
    (4, 0.1935), (4.2, 0.1861), (4.4, 0.1793), (4.6, 0.173), (4.8, 0.1672), (5, 0.1618),
]

P_STD_KPA, T_STD_C = 101.325, 15.0
BC_SI = 703.0696          # (kg/m^2) per (lb/in^2)


def air(t_c, p_kpa):
    """(density kg/m^3, speed of sound m/s) of dry air."""
    return p_kpa * 1e3 / (287.05 * (t_c + 273.15)), 331.3 * math.sqrt((t_c + 273.15) / 273.15)


def speed_curve(v0, bc, model, rho, c, s_max=2000.0, ds=0.5):
    """Speed v(s) and time of flight t(s) along the path (RK4 in distance)."""
    tab = np.array(TABLES[model])
    k = rho * math.pi / (8.0 * bc * BC_SI)
    f = lambda v: -k * np.interp(v / c, tab[:, 0], tab[:, 1]) * v
    s = np.arange(0.0, s_max + ds, ds)
    v = np.empty_like(s)
    v[0] = v0
    for i in range(1, len(s)):
        a = v[i - 1]
        k1 = f(a); k2 = f(a + 0.5 * ds * k1); k3 = f(a + 0.5 * ds * k2); k4 = f(a + ds * k3)
        v[i] = a + ds * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
    t = np.concatenate([[0.0], np.cumsum(ds * 0.5 * (1.0 / v[1:] + 1.0 / v[:-1]))])
    return s, v, t


def expo(s, v0, lam):
    return v0 * np.exp(-s / lam), lam * np.expm1(s / lam) / v0


def fit_row(r, t_c=T_STD_C, p_kpa=P_STD_KPA, mach_end=1.1, s_max=2000.0):
    rho, c = air(t_c, p_kpa)
    v0 = float(r["v0_mps"])
    s, v, t = speed_curve(v0, float(r["bc_lb_in2"]), r["drag_model"].strip().upper(), rho, c, s_max)
    sup = v / c >= mach_end
    n = len(s) if sup.all() else max(int(np.argmin(sup)), 2)
    if n <= 2:                                   # (barely) not supersonic: fit over the first 50 m
        n = int(50 / (s[1] - s[0]))
    S, V, T = s[:n], v[:n], t[:n]
    lam = minimize_scalar(lambda L: np.mean((expo(S, v0, L)[1] - T) ** 2),
                          bounds=(50.0, 50000.0), method="bounded").x
    ve, te = expo(S, v0, lam)
    return dict(calibre=r["calibre"], load=r["load"], v0=v0, lam=lam, mach0=v0 / c,
                s_supersonic=float(S[-1]) if sup[0] else 0.0,
                v300=float(np.interp(300, s, v)), v1000=float(np.interp(1000, s, v)),
                dv_max_pct=float(np.max(np.abs(ve - V) / V) * 100),
                dt_max_ms=float(np.max(np.abs(te - T)) * 1e3),
                d_mm=float(r["d_mm"]), L_mm=float(r["L_mm"]), confidence=r["confidence"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "ammo_data.csv"))
    ap.add_argument("--temp-c", type=float, default=T_STD_C)
    ap.add_argument("--pressure-kpa", type=float, default=P_STD_KPA)
    ap.add_argument("--mach-end", type=float, default=1.1)
    ap.add_argument("--ini", action="store_true", help="print det_config.ini lines")
    a = ap.parse_args()
    with open(a.csv, newline="", encoding="utf-8") as f:
        rows = [fit_row(r, a.temp_c, a.pressure_kpa, a.mach_end) for r in csv.DictReader(f)]
    if a.ini:
        for o in rows:
            print(f"{o['calibre']}_d = {o['d_mm'] / 1000:g}\n{o['calibre']}_L = {o['L_mm'] / 1000:g}\n"
                  f"{o['calibre']}_v0 = {o['v0']:.0f}\n{o['calibre']}_lam = {o['lam']:.0f}")
        return
    rho, c = air(a.temp_c, a.pressure_kpa)
    print(f"{a.temp_c:.1f} C, {a.pressure_kpa:.3f} kPa: rho {rho:.3f} kg/m^3, c {c:.1f} m/s; "
          f"lam fitted to time of flight while Mach >= {a.mach_end:g} (max 2000 m)")
    print(f"{'calibre':15s} {'v0':>5s} {'Mach0':>5s} {'lam':>6s} {'supersonic to':>13s} {'v(300m)':>7s} "
          f"{'v(1000m)':>8s} {'max dv':>6s} {'max dt':>7s}  confidence")
    for o in rows:
        print(f"{o['calibre']:15s} {o['v0']:5.0f} {o['mach0']:5.2f} {o['lam']:6.0f} {o['s_supersonic']:11.0f} m "
              f"{o['v300']:7.0f} {o['v1000']:8.0f} {o['dv_max_pct']:5.1f}% {o['dt_max_ms']:5.1f}ms  {o['confidence']}")


if __name__ == "__main__":
    main()
