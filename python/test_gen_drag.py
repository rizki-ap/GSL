#!/usr/bin/env python3
"""
test_gen_drag.py -- checks for the bullet drag model in gs_gen_physic.py

    python3 test_gen_drag.py          (or: pytest test_gen_drag.py)

  1. limit       : with a huge decay length, shock_geometry_drag equals the constant-speed
                   closed form of shock_geometry (arrival time, emission point, cone test)
  2. local Mach  : at the emission point, cos(phi) = c / v(s*)  (local Mach-angle condition)
  3. Fermat      : the arrival time is the minimum of t(s) + |x - p(s)|/c over the path
  4. subsonic    : a sensor whose emission point is below the Mach threshold gets no shock
  5. kinematics  : bullet_speed(subsonic_distance) = threshold speed; t(s) = integral ds / v
  6. end to end  : generate_clean with drag -- shock later and weaker-Mach than without,
                   muzzle blast unchanged, ground truth carries the new fields
"""
import numpy as np

import gs_gen_physic as gp
import gs_gen_clean_signal as gc

C = 343.0
RNG = np.random.default_rng(7)
U = np.array([1.0, 0.0, 0.0])
O = np.zeros(3)


def _sensors(n=200):
    """Random sensors 50..1500 m down range, 2..80 m off the path."""
    x = RNG.uniform(50, 1500, n)
    r = RNG.uniform(2, 80, n)
    ang = RNG.uniform(0, 2 * np.pi, n)
    return np.stack([x, r * np.cos(ang), r * np.sin(ang)], axis=1)


def test_limit_constant_speed():
    M = 2.5
    worst_t, worst_s = 0.0, 0.0
    for p in _sensors():
        g0 = gp.shock_geometry(p, O, U, M, C)
        g1 = gp.shock_geometry_drag(p, O, U, M * C, 1e15, C)
        assert g0["in_cone"] == g1["in_cone"]
        if g0["in_cone"]:
            worst_t = max(worst_t, abs(g0["t_arrive"] - g1["t_arrive"]))
            worst_s = max(worst_s, abs(g0["s"] - g1["s"]))
    assert worst_t < 1e-9 and worst_s < 1e-6, (worst_t, worst_s)
    return f"max |dt| {worst_t:.1e} s, max |ds| {worst_s:.1e} m"


def test_local_mach_angle():
    worst = 0.0
    for p in _sensors():
        g = gp.shock_geometry_drag(p, O, U, 940.0, 1000.0, C)
        if g["in_cone"]:
            cosphi = (g["a"] - g["s"]) / g["R"]
            worst = max(worst, abs(cosphi - C / g["v_emit"]))
    assert worst < 1e-9, worst
    return f"max |cos(phi) - c/v| {worst:.1e}"


def test_fermat_minimum():
    worst = 0.0
    for p in _sensors(50):
        g = gp.shock_geometry_drag(p, O, U, 940.0, 1000.0, C)
        if not g["in_cone"]:
            continue
        s = np.linspace(0.0, g["a"], 200001)
        f = gp.bullet_time(s, 940.0, 1000.0) + np.hypot(g["a"] - s, g["b"]) / C
        assert g["t_arrive"] <= f.min() + 1e-12
        worst = max(worst, f.min() - g["t_arrive"])
    assert worst < 1e-7, worst
    return f"grid minimum within {worst:.1e} s of the solver"


def test_subsonic_cutoff():
    v0, lam, thr = 940.0, 1000.0, 1.1
    s_sub = gp.subsonic_distance(v0, lam, thr * C)
    near = gp.shock_geometry_drag([s_sub - 300.0, 5.0, 0.0], O, U, v0, lam, C, mach_min=thr)
    far = gp.shock_geometry_drag([s_sub + 200.0, 5.0, 0.0], O, U, v0, lam, C, mach_min=thr)
    assert near["in_cone"] and near["M_emit"] > thr
    assert (not far["in_cone"]) and far["subsonic"] and far["M_emit"] <= thr
    return f"5.56 (lam 1000 m): Mach {thr} reached at {s_sub:.0f} m"


def test_kinematics():
    v0, lam = 940.0, 1000.0
    s_sub = gp.subsonic_distance(v0, lam, 1.1 * C)
    assert abs(gp.bullet_speed(s_sub, v0, lam) - 1.1 * C) < 1e-9
    s = np.linspace(0, 800, 800001)
    integrate = getattr(np, 'trapezoid', None) or getattr(np, 'trapz')   # numpy 2 / 1
    num = integrate(1.0 / gp.bullet_speed(s, v0, lam), s)
    assert abs(num - gp.bullet_time(800.0, v0, lam)) < 1e-9
    return "speed and time of flight consistent"


def test_end_to_end():
    def run(drag):
        cfg = gp.load_config('gen_config.ini')
        cfg['bullet'].update(caliber='5.56_NATO', mach=0.0, drag_model=drag)
        cfg['trajectory'].update(range=600.0, y_miss=25.0)
        return gc.generate_clean(cfg)[2]
    g0, g1 = run('none'), run('exponential')
    t0 = min(p['shockwave']['t_arrive_s'] for p in g0['per_mic'])
    t1 = min(p['shockwave']['t_arrive_s'] for p in g1['per_mic'])
    m0 = [p['muzzle_blast']['t_arrive_s'] for p in g0['per_mic']]
    m1 = [p['muzzle_blast']['t_arrive_s'] for p in g1['per_mic']]
    assert t1 > t0, "a decelerating bullet must arrive later"
    assert m0 == m1, "muzzle blast must not depend on drag"
    me = g1['per_mic'][0]['shockwave']['mach_emit']
    assert 1.0 < me < g1['mach']
    assert g1['drag_model'] == 'exponential' and g1['subsonic_distance_m'] > 0
    return (f"5.56 at 600 m: SW {1e3 * (t1 - t0):.0f} ms later, Mach at emission {me:.2f} "
            f"(muzzle {g1['mach']:.2f}), subsonic at {g1['subsonic_distance_m']:.0f} m")


if __name__ == '__main__':
    import warnings
    warnings.simplefilter('ignore')
    ok = True
    for name, fn in list(globals().items()):
        if name.startswith('test_') and callable(fn):
            try:
                print(f"PASS  {name}: {fn()}")
            except AssertionError as e:
                ok = False
                print(f"FAIL  {name}: {e}")
    print("ALL DRAG TESTS PASSED" if ok else "SOME DRAG TESTS FAILED")
