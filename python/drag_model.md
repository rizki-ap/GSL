# Bullet drag model (`gs_gen_physic.py`)

The generator previously flew every bullet in a straight line at **constant speed**. Its own docstring noted that this "shifts arrival times and TDOAs noticeably beyond a few tens of metres". This change adds an optional speed loss.

- **Default is unchanged:** `drag_model = none` reproduces the previous output **bit-for-bit** (checked on 6 configurations: default, 600 m, 7.62, 9 mm, explicit shooter, wind).
- **Version:** 2026-10-09

---

## Configuration (`gen_config.ini`, `[bullet]`)

| Key | Values | Meaning |
|---|---|---|
| `drag_model` | `none` (default), `exponential` | Constant speed, or v(s) = v₀·e^(−s/λ) |
| `decay_length_m` | 0 = library value | λ, the distance over which speed falls by a factor e |
| `caliber` | adds `12.7_NATO` | d = 12.95 mm, L = 58 mm, v₀ = 890 m/s, λ = 2500 m |

Library decay lengths, as **approximate** fits to typical published velocity tables:

| Calibre | v₀ (m/s) | λ (m) | Reaches Mach 1.1 at |
|---|---|---|---|
| 5.56_NATO | 940 | 1000 | ≈ 910 m |
| 7.62_NATO | 850 | 1200 | ≈ 975 m |
| 12.7_NATO | 890 | 2500 | ≈ 2140 m |
| 9mm_Parabellum | 370 | 400 | subsonic at the muzzle (no shockwave, as before) |

**Calibrating λ from your own ballistic table:** pick two distances s₁ < s₂ with speeds v₁, v₂, then λ = (s₂ − s₁) / ln(v₁ / v₂). The exponential model assumes a constant drag coefficient. Real drag rises near Mach 1, so the model is least accurate in the transonic region; a G7 drag table would be the next refinement.

## What the model does

- **Kinematics:** speed v(s) = v₀·e^(−s/λ), time of flight t(s) = λ(e^(s/λ) − 1)/v₀ (`bullet_speed`, `bullet_time`).
- **Shockwave arrival** (`shock_geometry_drag`): the first shock to reach a sensor is emitted at the path point s* that minimises t(s) + |x − p(s)|/c (Fermat's principle). Its condition is the **local** Mach-angle relation cos φ = c / v(s*). The constant-speed closed form in `shock_geometry` is the λ → ∞ limit.
- **No shockwave** if the local Mach at s* is ≤ `shockwave.mach_threshold`. The ground truth gives `inactive_reason = subsonic_at_emission`.
- **N-wave amplitude and duration** use the local Mach at s*, not the muzzle Mach.
- **Unchanged:** the muzzle blast, and gravity, which is still ignored (straight path).

New ground-truth fields:
- **top level:** `drag_model`, `decay_length_m`, `subsonic_distance_m`
- **per mic, shockwave:** `mach_emit`, `bullet_speed_emit_mps`

## Validation (`python3 test_gen_drag.py`)

| Test | Result |
|---|---|
| λ → ∞ equals the constant-speed closed form | max 1.2 ps, 56 pm |
| Local Mach angle cos φ = c / v at emission | max error 1.5e-14 |
| Arrival = minimum of t(s) + r(s)/c on a 200k-point grid | within 67 ps |
| Sensor beyond the subsonic point gets no shock | ✓ (5.56: Mach 1.1 at 913 m) |
| Speed / time-of-flight consistency | ✓ |
| End to end: 5.56 at 600 m | SW 223 ms later than constant speed; Mach at emission 1.54 vs 2.74 at the muzzle; muzzle blast identical |

---

## Finding: the constant-speed range formula is not usable for real bullets

PILAR's d = c·Δt / (1 − m·s) is exact for a constant-speed bullet. With drag, using **perfect** inputs (true Δt, true shooter direction m, true SW direction s):

| Calibre | 200 m | 400 m | 600 m | 800 m | 1000 m | 1500 m |
|---|---|---|---|---|---|---|
| constant speed (any) | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% | 0.0% |
| 5.56 with drag | +8% | +21% | +48% | +118% | no SW | no SW |
| 7.62 with drag | +8% | +19% | +41% | +91% | +270% | no SW |
| 12.7 with drag | +3% | +7% | +11% | +16% | +24% | +57% |

Spec: ±10% ≤ 200 m, ±20% ≤ 600 m, ±30% ≤ 1200 m.

A slowing bullet makes the shockwave arrive later, which shrinks Δt, and the formula reads that as a farther shooter. The exact percentages depend on the drag model, but the conclusion doesn't: **range needs a drag-aware solver**, in `gs_det_shooter_locator.py` (Milestone C software).

**What the detector can and cannot know.** The array observes only the **local** Mach at the emission point (cone geometry: n·u = 1/M_local), the N-wave duration and amplitude, and the directions and Δt. Muzzle velocity and drag are never observed. So the solver needs a ballistic table (v₀, λ per calibre class) that belongs to the **detector**: its own `det_config.ini` section with values from independent published data. It must **never** import this generator's `BULLET_LIBRARY`.

The solver works as a 1-D search: for each candidate calibre class, find the shooter distance along the MB direction whose forward model reproduces the measured Δt and SW direction; then pick the class most consistent with the measured N-wave duration.

**Testing it:** run the generator with drag values deliberately different from the detector's table (for example `decay_length_m` ±20%). The resulting range error is the real effect of ballistic uncertainty. Testing with matching values would only measure the solver's numerics.

## Finding: blast-only events are normal at long range

With drag, 5.56 has no shockwave beyond about 900 m and 7.62 beyond about 1000 m. The FPGA already handles this: Milestone B scenario `d556_1000` gives 0 SW and 1 MB event, with bearing 0.14° / 0.58°. The processor needs a **blast-only path**: bearing from the MB TDOAs, with range unknown. Possible estimates are MB amplitude, a calibre assumption, or another sensor.
