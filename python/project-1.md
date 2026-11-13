# Gunshot-Locator — Project Memory

## Project rules (Nov 2026)

Design limits for the generator, the detector and every test. Code and configs follow them.

| Rule | Value | Where it is applied |
|---|---|---|
| Maximum range, any calibre | 2000 m | generator refuses a shooter beyond 2000 m; detector range search to 2200 m, flag `beyond_max_range` above 2000 m |
| Minimum range for non-friendly detection | 50 m | detector flags a shot located closer than 50 m `friendly` (`[range] friendly_range_m`); the generator allows < 50 m for tests, with a note |
| Maximum burst rate | 600 rpm (rounds ≥ 100 ms apart) | detector `[detection] min_shot_interval_s = 0.09` for blast-only shots; burst tests at ≤ 600 rpm |
| Maximum calibre | 12.7 mm | generator refuses a bullet diameter above 13.0 mm (the 12.7×99 bullet measures 12.95 mm); both bullet libraries stop at 12.7×99 |
| Environment | outdoor only | generator refuses `noise.model = realistic` and the `room` layer; outdoor reflections use the `echoes` layer |
| Microphone array height | 2.0 m above ground (default) | `gen_config.ini [ground] array_height_m` |
| Supersonic bullets only | Mach ≥ 1.1 where the shock is emitted | detector `[detection] supersonic_only = true`: a muzzle blast without a shockwave is not reported (listed in `blast_only_events`); Stage 3 `[range] mach_min = 1.1`; generator `[shockwave] mach_threshold = 1.1` (no shock below it) |
| Gun height | 1.5 m above ground (default) | `gen_config.ini [trajectory] shooter_height_m`; in `range_miss` mode the shooter is at z = 1.5 − 2.0 = −0.5 m |

Recorded datasets may lie outside these limits (Cipatat shot-04, a MINIMI burst, fires at about 800 rpm, and its array stood 1.2–1.3 m above ground). They are still used for testing, and the detector still finds every round of shot-04: supersonic rounds are separated by pairing each shockwave with its blast, which needs no interval between rounds.

## Purpose

An acoustic gunshot detection and shooter-localization system modeled on the
PILAR-V (product reference: Boomerang). A 4-microphone tetrahedral
array detects supersonic gunshots and determines shooter range and bearing
by exploiting the two-event physics of a supersonic shot:

- The bullet's **shockwave** (N-wave, Whitham weak-shock model) arrives
  first and encodes the bullet's **trajectory direction** via Mach-cone
  geometry.
- The **muzzle blast** (Friedlander pulse) arrives second and encodes the
  **shooter's bearing**.
- The time gap between them (`Dt`) encodes **range**. Stage 3 solves it with
  a drag-aware bullet model (`[range] model = drag`, Nov 2026); the original
  constant-speed PILAR-V formula remains available as `constant`.

**Main publication reference:**
https://pub.dega-akustik.de/DAGA_1999-2008/data/articles/001903.pdf
(additional references in `/reference`)

**Gunshot recording datasets:**
- https://zenodo.org/record/7004819 — Publication: https://doi.org/10.1016/j.dib.2023.109091
- https://cadreforensics.com/audio/
- Local copies in `/sound` (Nov 2026), each folder with its own `metadata.md`:
  - `Cipatat/` — PILAR-V test-range recordings (4 ch, 51.2 kHz, SW + MB, 88–1470 m); the main real dataset
  - `Cipatat_clean/` — clean synthetic twins of the Cipatat development shots (`util_gen_cipatat_clean.py`)
  - `NIJ/` — AR-15 and .308 at 500 kHz, 12 mics on a 3 m semicircle around the gun (muzzle blast + near-field N-wave)
  - `Zenodo/` — 40 mono UMA-8 recordings at 9–27 m, 229 labelled shots, muzzle blast only
  - `Cadre/` — sample of the Cadre forensics dataset (M16)

## Repository Layout

```
/python          -- all simulation, generation, and detection code (this file lives here)
/sound           -- gunshot recordings and synthetic twins, one folder per dataset (see above)
/devices         -- sensor documentation (PILAR-V, UMA-8), processing-chain diagram
/publication     -- paper drafts (generator, detector) and figures
/reference       -- supporting publications
/verilog         -- FPGA/RTL prototyping (TDOA estimation, in progress)
/matlab          -- early GCC-PHAT experiments
```

## Microphone-array presets (`mic_config`, Nov 2026)

| Name | Geometry | Use |
|---|---|---|
| `PILAR-V` | Real PILAR-V sensor ("Trackball"): ring of 3 mics (radius 0.2055 m, 0/120/240°) + 1 mic 0.2412 m above; max TDOA ≈ 1.04 ms. See `devices/mic_pilar-v.md`. | Real recordings (Cipatat) |
| `Tetrahedral-0,3` | Regular tetrahedron, 0.30 m **edge** (circumradius 0.184 m); formerly named `PilarV` (old .json files still resolve). | Default in `gen_config.ini` and `det_config.ini` (simulations, FPGA milestones) |
| `UMA8` | miniDSP UMA-8 v2, 4 of its 7 mics | UMA-8 recordings |

## Pipeline Architecture

Two independent pipelines, plus shared utilities. **Hard rule: no file in
the detection pipeline imports from or reads ground truth produced by the
generation pipeline.** This is enforced architecturally (no shared imports,
duplicated reference constants where needed), not just by convention, so
that detection accuracy numbers are honest — the detector only ever sees a
`.wav` file and its own `det_config.ini`, exactly like a real system would.

### Generation block (`gs_gen_*.py`, `gen_config.ini`)

Simulates the acoustic signal a gunshot produces at each of the 4 mics,
in three sequential stages so environment and hardware can be varied
independently of the underlying physics.

| File | Role |
|---|---|
| `gs_gen_physic.py` | Shared physics/config library. Whitham N-wave, Friedlander blast, ISO 9613-1 atmospheric absorption, flat-ground image positions, analog gain-chain math, noise/reverb helpers, config schema. Imported by all other `gs_gen_*` files. |
| `gs_gen_clean_signal.py` | Stage 1. Pure physics, no noise: geometry → bullet property → trajectory → N-wave + Friedlander waveforms, plus the **ground reflection** (`[ground]`, image-source method, Nov 2026: the reflected shock is the direct-shock solution at each mic's mirror image, so it has its own delay, TDOAs, miss distance, amplitude and duration; constant reflection coefficient per event; calibrated on Cipatat). Outputs `*_clean.wav` (float32, **actual Pa units**, not normalized) + `*_clean.json` (full per-mic ground truth: arrival times, T_N, dP, etc). |
| `gs_gen_add_noise.py` | Stage 2. Adds the ENVIRONMENT only. `model = layers` (Nov 2026) combines any of: `colored` (piecewise power law, part common to all mics), `wind` (gusty, incoherent), `tones` (machinery, from set directions), `chirps` (insects), `diffuse` (noise from all directions, correct inter-mic coherence), `real` (noise cut from real recordings), `impulses` (interfering clicks/blasts/N-waves, times in the JSON), `echoes` (discrete outdoor echoes, each with its own direction) and `room` (diffuse tail with set direct-to-reverberant ratio); levels absolute or set by a target SNR (`-snr`). The older `simple` and `realistic` models remain. No mic/ADC effects. Outputs `*_noisy.wav` (still float32 Pa units) + `*_noisy.json`. |
| `gs_gen_apply_adc.py` | Stage 3. Mic bandpass → sensitivity calibration → ADC bit-depth quantization. Two calibration modes: `peak_match` (empirical, rescale to a measured real-recording peak) or `gain_chain` (absolute, computed from mic sensitivity/preamp gain/ADC Vref — can genuinely clip, which is the point). Supports per-channel mic/preamp variation and shared-ADC effects (gain/offset mismatch, and critically, **multiplexed-ADC timing skew**). Outputs the final 4-channel 16-bit PCM `.wav` + `.json`. |
| `gen_config.ini` | All generation parameters: geometry, bullet (caliber/Mach), trajectory (range/miss distance), atmosphere, noise model, ADC/analog-chain hardware. The bullet library (`BULLET_LIBRARY` in `gs_gen_physic.py`) covers the six small-arms calibres of the Cipatat development set (5.56×45, 7.62×51, 7.62×39, 7.62×54R, 8.6×70, 12.7×99) plus 9 mm. |

Run order: `gs_gen_clean_signal.py` → `gs_gen_add_noise.py` → `gs_gen_apply_adc.py`.

### Detection block (`gs_det_*.py`, `det_config.ini`)

Blind pipeline: reads only a `.wav` file (any 4-channel recording, real or
simulated) and its own config.

| File | Role |
|---|---|
| `gs_det_signal_prepare.py` | Stage 1. Blind SW/MB onset detection, two methods (`[detection] method`): **`classify`** (default, Nov 2026) labels every event SW or MB by frequency content and time/duration via `gs_det_events.py` — SW and MB are each optional and the MB need not follow a SW (blast-only and SW-only events supported) — then refines each channel's onset; **`threshold`** is the original per-channel detector (amplitude threshold for SW, adaptive energy-settling search for MB). Normalizes signal. Outputs `*_prepared.wav` + `*_prepared.json` (per-channel onset times/amplitudes, MB reliability flag, `event_class`, array-level events; with `multi_shot`, a `shots` list, top-level fields = shot 1). N-wave duration T_N per channel (Nov 2026): front shock to rear shock (steepest rises), 100 Hz high-pass, ×8 upsampled (`gs_det_events.nwave_duration`); the old 2 × (peak → zero crossing) in whole samples was 20–40 % short. Real-time timing (Nov 2026): MB waiting deadline `[detection] max_dt_s = auto` = the longest possible SW → MB delay of any library calibre (bullet at Mach 1.1, coldest air `deadline_min_temp_c`, + `deadline_margin`): 1.84 s with the current library; each shot gets an `early_report` at the SW (supersonic, SW direction, calibre category from the SW alone, MB deadline) and `t_final_report_s` (MB + 30 ms, or the deadline + 30 ms for SW only). |
| `gs_det_events.py` | Library: array-level SW/MB event classifier used by Stage 1 `classify` (and by `ana_sw_mb_range.py`). Thresholds tuned on real recordings (Cipatat); `mb_max_hf_excess_db = 12` (Nov 2026; close blasts are sharp). Stages 2–4 handle `SW only` (bullet direction, no bearing) and `MB only` (bearing, no range/Mach). **`detect_shots`** (Nov 2026, `[detection] multi_shot`): every shot in the file, for bursts and repeated shots. Supersonic: each SW is paired with the MB at the **common SW→MB delay** (the earliest delay matched by the most SWs), so echoes and the bursts' own MB trains are not paired. Subsonic, with an array: candidates more than `direction_gate_deg` from the level-weighted shooter direction are rejected as echoes. Mono without SW: all candidates reported, `echo_risk` set on those that follow a ≥ 3 dB stronger event within 300 ms. No minimum gap between shots other than `min_dt_s` (30 ms), so 600–1000 rpm bursts are kept. |
| `gs_det_tdoa.py` | Stage 2. TDOA across all 6 mic pairs (SW and MB separately). `[tdoa] method` (Nov 2026): **`xcorr`** (default) — plain cross-correlation of the band-passed event (SW 1–8 kHz, MB 150–1000 Hz) in short windows, sinc-upsampled; **`gcc_phat`** — the original method, which returns ~0 µs on recordings whose noise is common to all channels (Cipatat). Both use **shared-reference windowing** (all channels windowed around ONE reference channel's onset — critical, see Known Limitations). Outputs `*_tdoa.json` (6 TDOA pairs × 2 events, plus `Dt`; one set per shot in `shots`). MB pair TDOAs must add up (`tau(a,b) = tau(0,b) − tau(0,a)`); above `mb_closure_tol_s` (50 µs) the MB window is shortened to 1/2 and 1/4 and the most consistent set kept, because in a burst the SW or an echo of a later round can fall into the 30 ms MB window. |
| `gs_det_shooter_locator.py` | Stage 3. Solves DOA (bearing) from MB TDOA; solves the SW direction `k_hat_sw` from SW TDOA. Range (`[range] model`, Nov 2026): **`drag`** (default) fits shooter distance and trajectory angle per caliber of the detector's own `[bullet_reference_library]` (every calibre with `v0`; the patch's separate `[ballistics]` table was merged into it) (v(s) = v0·e^(−s/λ), local Mach-angle shock emission) to the SW angle, the SW→MB delay and the N-wave duration, and reports the best caliber, range, miss distance, bullet speed and 3D trajectory; **`constant`** is the old constant-speed PILAR-V solve. Speed of sound from `[physics] temperature_c` / `--temperature-c`; the same temperature scales each calibre's `lam` with the air density (`lam × (T + 273.15) / 288.15`, `[range] lam_temperature_scaling`). Doesn't depend on Stage 4. Validated against the generator's shock model by `test_det_range_drag.py`. Calibre choice (Nov 2026): reported as category (`caliber_category`: small, medium, big, RPG, other; always present) and optionally class (`caliber_class`, bore number); the variant is not reported (`range_ballistics` names the one used for the range); the fit also uses the SW peak in Pa (Whitham dP, `[range] dp_sigma_log`, calibration `assumed_peak_pa_at_full_scale`), which with T_N fixes miss distance and bullet size; T_N is the channel median minus the sensor's response offset (`t_n_offset_s`, 35 µs for PILAR-V, from Cipatat). Several shots: each is located into `shots`; a shot more than 5° from the median bearing of the file is flagged `bearing_outlier`. |
| `gs_det_classify_bullet.py` | Stage 4. Three independent Mach-estimation methods (A: angle-only, calibration-independent; B: duration+amplitude, needs the detector's own assumed hardware calibration; C: geometric correction using B's result). Top-level shot only (shot 1 with several shots). Tests all calibers in the reference library (`det_config.ini` `[bullet_reference_library]`: the same six calibres, kept as an independent copy). Method B fits Whitham's N-wave equations (the generator's physics; checked on Cipatat, Nov 2026). Amplitude = Stage-1 normalized dP × `normalization_peak` × `assumed_peak_pa_at_full_scale` (1520 Pa, the PILAR-V estimate; per-array override `assumed_peak_pa_at_full_scale.<mic_config>`). |
| `det_config.ini` | Array geometry (its own installed hardware, legitimate to know), detection thresholds, TDOA settings (`xcorr` / `gcc_phat`), an **independent copy** of the bullet reference library, and the detector's own assumed sensitivity calibration (`assumed_peak_pa_at_full_scale`). |

Run order: `gs_det_signal_prepare.py` → `gs_det_tdoa.py` → `gs_det_shooter_locator.py` → `gs_det_classify_bullet.py` (last two can run independently once Stage 2 is done).

### Utilities (`util_*.py`)

| File | Role |
|---|---|
| `util_plot_signal.py` | Plots every channel of a multichannel `.wav`: full-duration + auto-zoomed view around the loudest event, overlay mode, handles any dtype (int PCM normalized, float kept in raw units — important since generation-block float WAVs hold real Pa values). |
| `util_get_accoustic_param.py` | Estimates noise/reverb/calibration parameters (RT60 via Schroeder backward integration with multi-window fit-quality reporting, noise floor, PSD slope, approximate mic passband) directly from a **real recording** — this is what should have produced `gen_config.ini`'s `[calibration]` constants in the first place. Also has `--check-channels` to flag bad channels (inverted polarity, low array-consensus correlation) across a multi-mic dataset before using any of them as a reference. Self-contained (no external project dependency). |
| `util_signal_similarity.py` | Measured vs synthesized gunshots, per Paper-Generator Section 4.3: Pearson r, SNR_synth, peak and duration errors (T_N / T_mb measured from the signals), event-PSD LSD, blast f_pk and N-wave f_null, pooled coherence (only with ≥ 8 segments), multichannel TDOA MAE/RMSE/Bias (vs generator labels with `--synth-json`), shock-to-blast delay error. Pools several recordings. See `util_signal_similarity.md`. |
| `ana_ground_reflection.py` | Measures the ground-reflected shockwave in the Cipatat recordings (residual after removing the clean twin's direct wave, template match): delay, reflected/direct ratio, and the array height from the ring/top delay ratio, independent of the trajectory. |
| `util_gen_cipatat_clean.py` | Clean twins of the Cipatat development recordings (01–06, 06a, 15, 16, 17) in `sound/Cipatat_clean/`: the reference scene rebuilt with Stage 1 (PILAR-V, exponential drag, per-shot temperature), bullet side chosen from the SW TDOAs, SW and MB levels and durations fitted to each event, events placed at the recording's per-channel arrivals. No noise, reverberation or ground reflection. `--align array` keeps the model TDOAs, `--no-fit` the model levels. See `sound/Cipatat_clean/metadata.md`. |
| `util_ballistic_library.py` | Computes the bullet libraries' `v0` and decay length `lam` from ammunition data (`ammo_data.csv`: muzzle velocity, G1/G7 ballistic coefficient, mass, d, L, source and confidence per calibre). Flat-fire point-mass solution with the standard G1/G7 drag tables (embedded), ICAO atmosphere (15 °C, 101.325 kPa); `lam` fitted to the time of flight while the bullet is above Mach 1.1 (max 2000 m). `--ini` prints the `det_config.ini` lines. Nov 2026. |

## Key Physics & Algorithms Implemented

- **Whitham weak-shock N-wave model** (Nov 2026):
  `dP = 0.53·p0·(M²−1)^(1/8)·d / (L^(1/4)·b^(3/4))`,
  `T_N = 1.82·M·d·(b/L)^(1/4) / (c·(M²−1)^(3/8))`, with bullet diameter d and
  length L from the bullet library (six Cipatat calibres + 9 mm). Checked on
  Cipatat: N-wave durations within ~10 %. Optional exponential drag
  `v(s) = v0·e^(−s/λ)` with the shock emitted where the local Mach angle fits.
- **Ground reflection** (optional, `[ground]`): image-source method, constant
  reflection coefficient per event; calibrated on Cipatat (array 1.2–1.3 m
  above ground, reflected/direct shockwave ≈ 0.27).
- **Friedlander muzzle blast** with Hopkinson-Cranz cube-root scaling.
- **Bullet libraries** (generator `BULLET_LIBRARY`, detector
  `[bullet_reference_library]`, kept as separate copies): `v0` and `lam` computed
  from ammunition data by `util_ballistic_library.py` (Nov 2026; before: rough
  hand-fitted values, time-of-flight errors up to 17 % for 12.7 mm). The exponential
  model then matches the G1/G7 solution within 2–25 ms of flight time (≤ 1.2 %).
  `lam` is for 15 °C and is scaled with air density by both blocks.
- **Range**: drag-aware fit of shooter distance and trajectory angle per
  calibre to the SW angle, `Dt` and T_N (Stage 3 `drag`, validated against
  the generator by `test_det_range_drag.py`); the constant-speed PILAR-V
  formula (`brentq` root-finding) is kept as `constant`.
- **Cross-correlation** of the band-passed event for TDOA (default since Nov 2026; GCC-PHAT optional), frequency-domain zero-padded for sub-sample
  resolution (validated to <10µs on clean signals).
- **ISO 9613-1 atmospheric absorption** — validated against published
  reference values (109.8 dB/km @ 4kHz/10%RH vs. reference 109; 23.1 vs.
  23 @ 70%RH, both <1% error). Applied to each event's waveform as a causal
  (minimum-phase) frequency-domain filter over its own path length.
- **Analog gain chain**: `PA_TO_NORM = mic_sensitivity_V/Pa · 10^(gain_dB/20) / V_ref`.
  Supports per-channel mic/preamp variation and shared-ADC effects
  including fractional-sample timing skew for multiplexed ADCs (FFT-based
  phase-shift delay, validated to <0.001 sample error).

## Known Limitations & Important Findings

These are validated, reproducible findings from extensive testing this
project's own tools — not guesses:

1. **The `realistic` noise model (RT60=1.25s) badly breaks GCC-PHAT
   at short range.** At the default 200m/Mach 2.5 config, `Dt≈324ms <
   RT60=1.25s`, so the shockwave's reverb tail is still ringing when the
   muzzle blast arrives — the two events' reverberant energy overlaps.
   Confirmed via exact-RIR deconvolution (recovers sub-10µs TDOA, proving
   it's a reverb problem, not a fundamental one) and blind WPE
   dereverberation (does NOT help — there isn't enough reverberant
   excitation in one brief impulsive event for WPE's frame-prediction
   approach to work; it needs continuous excitation like speech). The
   practical fix for a real fixed-position sensor: measure the room
   response once with a calibration sweep, reuse it for deconvolution on
   every subsequent live detection. (Since Nov 2026 `gen_config.ini` uses
   `model = layers` with outdoor noise and discrete echoes; `realistic` is
   kept only for old configurations.)
2. **Simple (flat-noise-floor) model works excellently.** GCC-PHAT on it
   reproduces true TDOA to single-digit microseconds — this is what
   validated the detection pipeline's own logic is correct, independent of
   the reverb issue above.
3. **Shared-reference windowing is mandatory for GCC-PHAT.** Windowing each
   channel around its OWN independently-detected onset silently erases the
   very TDOA being measured (each window ends up centered on that channel's
   own arrival, cancelling the delay). All working code windows every
   channel around ONE shared reference channel's onset.
4. **Bearing (TDOA-only) is far more robust than range/Mach (amplitude-
   based).** Validated on the blind detection pipeline: bearing came out
   accurate to ~0.5° regardless of hardware calibration assumptions, since
   TDOA doesn't depend on absolute amplitude at all. Range and Mach
   (Methods B/C) depend on the detector's own assumed sensitivity
   calibration matching reality — get that assumption wrong (which, under
   the no-parameter-passing rule, it generally will vs. whatever the
   generation side used) and amplitude-based estimates are systematically
   biased even with perfect timing measurements.
5. **Fixed `dt_min` search gaps silently fail depending on RT60.** A blind
   detector can't assume it knows the true reverb time (that's generation
   ground truth). The fix: an adaptive detector that waits for the post-SW
   energy to actually settle near the noise floor before starting the MB
   search, rather than a fixed time gap.
6. **Shockwave amplitude/duration don't change with range** — only with
   perpendicular miss distance, a real feature of the far-field Whitham
   model (the shockwave forms along the whole trajectory; what matters is
   closest approach, not total range to the shooter). Muzzle blast
   amplitude does scale correctly with range (validated: 7.5× range
   increase → 7.5× amplitude decrease, matching the 1/r law exactly).
7. **Geometric Mach-cone validity**: `b ≤ a·β` where `β=√(M²−1)` — higher
   Mach numbers tolerate LARGER miss distances geometrically (narrower cone
   angle means larger β), counter to naive intuition.
8. **Channel quality issues found in the real datasets**: in the earlier
   per-channel Zenodo files (Ruger and Glock; removed from `/sound` in Nov
   2026 because their inter-channel timing is not physical) ch1 showed a
   consistent polarity inversion, suggesting a fixed wiring characteristic
   at that array position rather than a one-off fault. `util_get_accoustic_param.py
   --check-channels` catches this automatically via correlation against
   the array's own consensus mean.

9. **Echoes and bursts overlap in time** (Nov 2026). Outdoor echoes arrive
   25–200 ms after a shot (Zenodo: median 95 ms); rounds at the project's
   maximum of 600 rpm are 100 ms apart. `detect_shots` separates them by SW
   pairing (supersonic), by direction (array), and, for blast-only shots, by
   the 600 rpm limit (of two blasts closer than 90 ms the weaker is dropped).
   Generator bursts (5 rounds, ±4 dB, 5 seeds, PILAR-V 2.0 m above ground,
   gun 1.5 m, 200 m; sparse and dense echo sets, 300/450/600 rpm):
   supersonic 24–25/25 rounds and no false shots in every case;
   subsonic with sparse echoes 25/25 rounds and no false shots with an array
   (5–10 false without); subsonic with dense echoes 20–25/25 rounds and
   0–9 false shots. Zenodo (mono, blast only): 218 of 229 shots found,
   160 extra events (before the 90 ms rule: 930), 106 of them flagged
   `echo_risk`. Cipatat: shot-04 5/5 rounds, shot-17 5/5 shots.
   Near shots: a blast inside the SW's echo tail (50–100 m) is accepted when
   the array shows the SW and the blast at most `mb_tail_max_angle_deg`
   (75°) apart; generator shots at 45–100 m located within 0.2° and 1–3 %
   in range, 45 m flagged `friendly`.
10. **shot-17 range is +20–25 %** (88 m reference) although all 5 shots are
   found. The fitted miss distance is 15–35 m instead of 4 m: the SW
   direction is biased (the top microphone arrives 70–160 µs early, elevation
   2–4° instead of 0°), and at 88 m the range fit is very sensitive to the SW
   angle (each degree ≈ 3 m). The measured SW→MB delay is also 7 ms longer
   than the reference geometry gives. The calibre is now picked right
   (7.62×39, since the new T_N estimator and the amplitude term).

11. **The shockwave range only works while the bullet is still supersonic
   at the array** (computed library, 15 °C, Mach ≥ 1.1): 5.56 to ~650 m,
   7.62×39 to ~410 m, 7.62×51 and 7.62×54R to ~760 m, .338 and 12.7 mm to
   ~1300 m. Farther shots within the project's 2000 m give a muzzle blast
   only: bearing, but no range from the SW → MB delay. `test_det_range_drag.py`
   marks those cases. Cipatat with the computed library and temperature
   scaling: range RMS 4.5 % on the 8 development shots (shot-03 against 431 m;
   7.8 % with the old library), calibre right in 5 of 8 (4 before).
12. **`ammo_data.csv` sources** (checked 13 Nov 2026, in-app browser):
   5.56 SS109/M855 948 m/s, G7 0.151 (US Army BRL, Wikipedia); 7.62×51 M80
   838 m/s, G7 0.200 (BRL, Wikipedia); 7.62×39 57-N-231 718 m/s, G7 0.152
   and 7.62×54R 7N1 823 m/s, G7 0.206 (Wikipedia, Russian service-load
   data); .338 Lapua 250 gr Scenar GB488 905 m/s, G7 0.322, 39.4 mm
   (manufacturer data); 12.7 M33 ball 887 m/s, G1 0.70 (manufacturer G1 of
   M33-type loads 0.670–0.707). Only 9 mm (not supersonic, unused) is still
   unsourced. Check against the Cipatat SW → MB delays at the reference
   geometry: 5.56 within ±1 %, 7.62×39 −3.5 to −7 %, 12.7 −5.0 % (was
   −7.6 %); RMS 5.1 % over 8 shots (03 left out); shots 03, 04, 05 have
   doubtful reference geometry. PINDAD MU-series data, when available,
   should replace these (`util_ballistic_library.py --ini`).

13. **Calibre identification** (Nov 2026, Cipatat, 9 cuts). Project rule:
   three levels, coarse to fine. **Category**: small (< 7.62 mm), medium
   (7.62 to < 12.7 mm), big (12.7 to 40 mm), RPG, other. **Class**: the
   nominal bore number covering all its variants (5.56, 7.62, 8.6, 12.7, 20,
   …). **Variant**: one `[bullet_reference_library]` entry with its own
   properties (7.62×39, 7.62×51, 7.62×54R, .338 Lapua = 8.6×70, …). Class and
   category come from the name (`caliber_class_of`, `caliber_category_of` in
   `gs_det_signal_prepare.py`; override `<cal>_class`, `<cal>_category`).
   Output rule: **category mandatory** (`caliber_category`, null with
   `caliber_category_reason` when it cannot be determined), **class
   optional** (`caliber_class`, `[output] report_caliber_class`), **variant
   not reported** (only `range_ballistics`, the variant whose ballistics gave
   the range, and the fit details). Stage 4: `best_guess_caliber_category`
   and optional `best_guess_caliber_class`.
   Results (with the sourced bullet library): category right 9/9 (Stage 3)
   and 9/9 (Stage 4); class 9/9 and 8/9 (Stage 4: 05 .338 → 7.62, still
   medium); variant: 7.62×51 and 7.62×54R swap on 03 and 04 (almost the
   same bullet). Shot-06 (12.7) became right with the M33 coefficient
   0.70 (before: 0.62, taken as .338). Margins are small for 01 (small 1.98 vs medium 2.07) and 03.
   Generator checks (5.56, 7.62×51, 7.62×39, .338 at 150–600 m): all right;
   measured T_N there is 11–16 µs above the generator's, against 35 µs on
   PILAR-V.
14. **Generator shots near the limits** (Nov 2026). Weak shocks near Mach 1.1
   rise slower: generated 7.62×39 at 350 m 0.095 ms, 12.7 mm at 800 m
   0.085 ms (10–90 %), so `sw_max_rise_ms` was raised from 0.08 to 0.10 ms
   (Cipatat unchanged, Zenodo blasts taken as SW still 1 of 229). With it:
   7.62×39 at 350 m, 5.56 at 600 m and 12.7 mm at 800 m are found and
   located within 1 % and 0.3°; 7.62×39 at 400 m has no shockwave left
   (bullet below Mach 1.1) and is correctly not reported. Still failing:
   12.7 mm at 1200 m (the weak blast is confused with an echo, bearing
   wrong, no range); 7.62×51 at 700 m is found but 7 % short and taken as
   class 5.56 (bullet near Mach 1.1). Files must start ≥ 0.1 s before the
   first event (`[timing] pre_roll`): a shockwave in the first 10 ms is not
   found.

15. **Real-time timing and the shockwave-only category** (Nov 2026). The MB
   deadline bound holds for any miss distance: a shock heard at the array was
   emitted at s_e with the bullet at Mach ≥ 1.1, the shooter is at most
   s_e + R away, so SW → MB ≤ s_e/c − t_bullet(s_e), largest where the bullet
   slows to Mach 1.1. Library at 5 °C: 12.7 1.76 s, .338 1.67 s, 5.56 and
   7.62×51 ~0.87 s, 7.62×39 0.42 s; deadline 1.84 s, so a final report is
   always issued within ~1.9 s of the SW. The category from the SW alone
   (Whitham T_N and peak, Mach unknown between 1.1 and the muzzle Mach) is
   only a hint: right on 7 of 10 Cipatat cuts (first shot) and 8 of 18
   shots; small and medium overlap (03, 04 and 17 come out small, often
   flagged ambiguous). It is therefore NOT used to shorten the deadline.

16. **Clipped shockwaves (UMA-8)** (Nov 2026). The UMA-8 capsules overload at
   about 36 Pa peak (122 dB SPL, 10 % THD), so bullets passing within ~60–180 m
   (5.56 … 12.7 mm) clip the shockwave. Stage 1 flags a clipped SW (`sw_clipped`:
   flat top ≥ `clip_flat_samples` samples, or peak ≥ 90 % of `[calibration]
   sw_clip_pa[.UMA8]`) and estimates the true peak from the flat-top length
   (N-wave: P = C / (1 − 2τ/T_N), `dP_sw_unclipped_est`). Stage 3 then uses the
   measured peak as a lower bound plus that estimate with a looser uncertainty
   (`dp_clip_sigma_log` = 0.5); Stage 4 warns. Simulation, UMA8 geometry,
   48 kHz, hard clip at 36 Pa, 30 cases (5.56, 7.62×51, 7.62×39, .338, 12.7;
   200/400 m; miss 10/30/60 m): with the new handling the category is right in
   24 of 25 located cases and range RMS 4.3 %; taking the clipped peak at face
   value: 10 of 25 and 5.7 % (it calls most 7.62 and 12.7 shots one category
   too small). Bearing is unaffected by clipping (RMS 2.4°). Five cases fail in
   every mode because of the UMA-8 itself: the weak blast at 200–400 m gives a
   wrong bearing (73.6 mm aperture). Real capsule overload (sigma-delta
   instability, decimation ringing) may be worse than a clean clip: check with
   real shots. Also found: unclipped generated shockwaves resampled to 48 kHz
   rise in 0.104 ms (> `sw_max_rise_ms` 0.10) and are missed; the generator's
   fronts are slower than real ones (Cipatat 0.049 ms), to be looked at.

## Validated Reference Numbers

For the default `gen_config.ini` (range_miss 200 m, y_miss 10 m, Mach 2.5,
no drag, 5.56 NATO, Tetrahedral-0,3 array 2.0 m above ground, gun 1.5 m above
ground, 20 °C; Stage-1 ground truth, Nov 2026):
- Shockwave: `T_N ≈ 0.181 ms`, `dP ≈ 168 Pa` at b = 10 m (Whitham; identical at
  1500 m — see Known Limitations #6)
- Muzzle blast @ 200 m: `dP ≈ 52.8 Pa`, `t_pos ≈ 0.53 ms`; @ 1500 m: `dP ≈ 7.0 Pa`
  (Kinney-Graham, extrapolated beyond its range: far too strong against Cipatat,
  see `sound/Cipatat_clean/metadata.md`)
- `Dt` @ 200 m: `≈ 324 ms`; @ 1500 m: `≈ 2.60 s`
- True SW TDOA (pairs 0-1, 0-2, 0-3): `[-498.9, +19.5, +368.9] µs` (gun 0.5 m below the array; `[-466.7, +52.3, +401.8] µs` with both at the same height)
- True MB TDOA (pairs 0-1, 0-2, 0-3): `[-27.0, -425.6, +446.9] µs`
- Ground reflection on (array 2.0 m above ground): reflected shock 1.50–1.69 ms
  after the direct one, `dP ≈ 53 Pa` (with the array at 1.3 m: 0.84–1.01 ms)

## On the Horizon

- Hardware: DE10-Nano FPGA target. Onboard LTC2308 ADC is **multiplexed**
  (sequential channel sampling), which the simulation now models
  explicitly (`adc_multi_channel.type=multiplexed`) — confirmed via
  simulation that this introduces real, growing-with-channel-index TDOA
  error (~2-7µs per channel at 2µs/channel conversion time). An external
  simultaneous-sampling ADC (AD7606B / ADS8688) remains the recommended
  fix, now with simulated evidence quantifying why.
- Verilog/RTL prototyping of TDOA estimation and localization (in
  `/verilog`) — drafted in earlier sessions, not yet re-validated against
  the current Python pipeline's numbers.
- A far-field muzzle-blast model calibrated on real recordings (the
  Kinney-Graham blast is 3–90× too strong at 88–764 m) and blast directivity
  from the NIJ data.
- Making `mic_bandpass()` support genuinely different frequency responses
  per channel, not just per-channel gain/sensitivity (currently one shared
  filter for the whole array).

## Quick Usage

```bash
# Generate a simulated signal
python gs_gen_clean_signal.py gen_config.ini -o shot
python gs_gen_add_noise.py shot_clean.wav gen_config.ini
python gs_gen_add_noise.py shot_clean.wav gen_config.ini -o shot20 -layers colored,wind,tones,echoes -snr 20
python gs_gen_apply_adc.py shot_noisy.wav gen_config.ini

# Run detection on it (or on a real recording)
python gs_det_signal_prepare.py shot.wav det_config.ini -o result
python gs_det_tdoa.py result_prepared.json det_config.ini -o result
python gs_det_shooter_locator.py result_tdoa.json det_config.ini -o result
python gs_det_classify_bullet.py result_tdoa.json det_config.ini -o result

# Estimate real-world calibration parameters from an actual recording
python util_get_accoustic_param.py real_recording.wav --config-out estimated.ini
python util_get_accoustic_param.py --check-channels rec_ch0.wav rec_ch1.wav ... rec_ch6.wav

# Visualize / compare
python util_plot_signal.py shot.wav
python util_signal_similarity.py real.wav sim_clean.wav --synth-json sim.json --meas-scale 1520 --plot out/
python util_gen_cipatat_clean.py                # clean Cipatat twins -> ../sound/Cipatat_clean
```
