# Gunshot-Locator — Project Memory

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
- The time gap between them (`Dt`) encodes **range**, via the PILAR-V
  formula.

**Main publication reference:**
https://pub.dega-akustik.de/DAGA_1999-2008/data/articles/001903.pdf
(additional references in `/reference`)

**Gunshot recording datasets:**
- https://zenodo.org/record/7004819 — Publication: https://doi.org/10.1016/j.dib.2023.109091
- https://cadreforensics.com/audio/
- Local copies in `/sound`. Current datasets (each 7 channels, ch0-ch6 +
  mean): `glock_a283`, `ruger_04_1s`, `ruger_ar556_223cal_97V0`,
  `SmithWesson38cal_47V0`, `remington12g_0aV1`.

## Repository Layout

```
/python          -- all simulation, generation, and detection code (this file lives here)
/sound            -- real gunshot recordings (7ch per dataset + mean)
/reference        -- supporting publications
/verilog          -- FPGA/RTL prototyping (TDOA estimation, in progress)
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
| `gs_det_signal_prepare.py` | Stage 1. Blind SW/MB onset detection, two methods (`[detection] method`): **`classify`** (default, Nov 2026) labels every event SW or MB by frequency content and time/duration via `gs_det_events.py` — SW and MB are each optional and the MB need not follow a SW (blast-only and SW-only events supported) — then refines each channel's onset; **`threshold`** is the original per-channel detector (amplitude threshold for SW, adaptive energy-settling search for MB). Normalizes signal. Outputs `*_prepared.wav` + `*_prepared.json` (per-channel onset times/amplitudes, MB reliability flag, `event_class`, array-level events). |
| `gs_det_events.py` | Library: array-level SW/MB event classifier used by Stage 1 `classify` (and by `ana_sw_mb_range.py`). Thresholds tuned on real recordings (Cipatat). Stages 2–4 handle `SW only` (bullet direction, no bearing) and `MB only` (bearing, no range/Mach). |
| `gs_det_tdoa.py` | Stage 2. TDOA across all 6 mic pairs (SW and MB separately). `[tdoa] method` (Nov 2026): **`xcorr`** (default) — plain cross-correlation of the band-passed event (SW 1–8 kHz, MB 150–1000 Hz) in short windows, sinc-upsampled; **`gcc_phat`** — the original method, which returns ~0 µs on recordings whose noise is common to all channels (Cipatat). Both use **shared-reference windowing** (all channels windowed around ONE reference channel's onset — critical, see Known Limitations). Outputs `*_tdoa.json` (6 TDOA pairs × 2 events, plus `Dt`). |
| `gs_det_shooter_locator.py` | Stage 3. Solves DOA (bearing) from MB TDOA; solves the SW direction `k_hat_sw` from SW TDOA. Range (`[range] model`, Nov 2026): **`drag`** (default) fits shooter distance and trajectory angle per caliber of the detector's own `[bullet_reference_library]` (every calibre with `v0`; the patch's separate `[ballistics]` table was merged into it) (v(s) = v0·e^(−s/λ), local Mach-angle shock emission) to the SW angle, the SW→MB delay and the N-wave duration, and reports the best caliber, range, miss distance, bullet speed and 3D trajectory; **`constant`** is the old constant-speed PILAR-V solve. Speed of sound from `[physics] temperature_c` / `--temperature-c`. Doesn't depend on Stage 4. Validated against the generator's shock model by `test_det_range_drag.py`. |
| `gs_det_classify_bullet.py` | Stage 4. Three independent Mach-estimation methods (A: angle-only, calibration-independent; B: duration+amplitude, needs the detector's own assumed hardware calibration; C: geometric correction using B's result). Tests all calibers in the reference library (`det_config.ini` `[bullet_reference_library]`: the same six calibres, kept as an independent copy). Method B fits Whitham's N-wave equations (the generator's physics; checked on Cipatat, Nov 2026). Amplitude = Stage-1 normalized dP × `normalization_peak` × `assumed_peak_pa_at_full_scale` (1520 Pa, the PILAR-V estimate; per-array override `assumed_peak_pa_at_full_scale.<mic_config>`). |
| `det_config.ini` | Array geometry (its own installed hardware, legitimate to know), detection thresholds, TDOA settings (`xcorr` / `gcc_phat`), an **independent copy** of the bullet reference library, and the detector's own assumed sensitivity calibration (`assumed_peak_pa_at_full_scale`). |

Run order: `gs_det_signal_prepare.py` → `gs_det_tdoa.py` → `gs_det_shooter_locator.py` → `gs_det_classify_bullet.py` (last two can run independently once Stage 2 is done).

### Utilities (`util_*.py`)

| File | Role |
|---|---|
| `util_plot_signal.py` | Plots every channel of a multichannel `.wav`: full-duration + auto-zoomed view around the loudest event, overlay mode, handles any dtype (int PCM normalized, float kept in raw units — important since generation-block float WAVs hold real Pa values). |
| `util_get_accoustic_param.py` | Estimates noise/reverb/calibration parameters (RT60 via Schroeder backward integration with multi-window fit-quality reporting, noise floor, PSD slope, approximate mic passband) directly from a **real recording** — this is what should have produced `gen_config.ini`'s `[calibration]` constants in the first place. Also has `--check-channels` to flag bad channels (inverted polarity, low array-consensus correlation) across a multi-mic dataset before using any of them as a reference. Self-contained (no external project dependency). |
| `util_signal_similarity.py` | Measured vs synthesized gunshots, per Paper-Generator Section 4.3: Pearson r, SNR_synth, peak and duration errors (T_N / T_mb measured from the signals), event-PSD LSD, blast f_pk and N-wave f_null, pooled coherence (only with ≥ 8 segments), multichannel TDOA MAE/RMSE/Bias (vs generator labels with `--synth-json`), shock-to-blast delay error. Pools several recordings. See `util_signal_similarity.md`. (`util_signal_similarity_paper.py` forwards to it.) |
| `ana_ground_reflection.py` | Measures the ground-reflected shockwave in the Cipatat recordings (residual after removing the clean twin's direct wave, template match): delay, reflected/direct ratio, and the array height from the ring/top delay ratio, independent of the trajectory. |
| `util_gen_cipatat_clean.py` | Clean twins of the Cipatat development recordings (01–06, 06a, 15, 16, 17) in `sound/Cipatat_clean/`: the reference scene rebuilt with Stage 1 (PILAR-V, exponential drag, per-shot temperature), bullet side chosen from the SW TDOAs, SW and MB levels and durations fitted to each event, events placed at the recording's per-channel arrivals. No noise, reverberation or ground reflection. `--align array` keeps the model TDOAs, `--no-fit` the model levels. See `sound/Cipatat_clean/metadata.md`. |

## Key Physics & Algorithms Implemented

- **Whitham weak-shock N-wave model**: `dP = dP_ref·√(b_ref/b)`,
  `T_N = (2/Mc)·√(bL/β)`, `β = √(M²−1)`. Calibration constants
  (`dP0_sw=7.5 Pa @ b0_sw=50m`) are fit for **7.62×51 NATO at Mach 3**
  specifically — applying them at other Mach numbers extrapolates outside
  the fitted point (the amplitude law currently has no explicit Mach
  dependence at all, a known simplification).
- **Friedlander muzzle blast** with Hopkinson-Cranz cube-root scaling.
- **PILAR-V range formula**: solves `Dt = t_MB − t_SW` for range given
  Mach and miss distance, via root-finding (`scipy.optimize.brentq`).
- **Cross-correlation** of the band-passed event for TDOA (default since Nov 2026; GCC-PHAT optional), frequency-domain zero-padded for sub-sample
  resolution (validated to <10µs on clean signals).
- **ISO 9613-1 atmospheric absorption** — validated against published
  reference values (109.8 dB/km @ 4kHz/10%RH vs. reference 109; 23.1 vs.
  23 @ 70%RH, both <1% error). Applied per-pulse using a single
  characteristic frequency (`1/T_N` or `1/t_pos`), not full spectral
  reshaping of the waveform.
- **Analog gain chain**: `PA_TO_NORM = mic_sensitivity_V/Pa · 10^(gain_dB/20) / V_ref`.
  Supports per-channel mic/preamp variation and shared-ADC effects
  including fractional-sample timing skew for multiplexed ADCs (FFT-based
  phase-shift delay, validated to <0.001 sample error).

## Known Limitations & Important Findings

These are validated, reproducible findings from extensive testing this
project's own tools — not guesses:

1. **The `realistic` noise model (RT60=1.25s default) badly breaks GCC-PHAT
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
   every subsequent live detection.
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
8. **Channel quality issues found in the real datasets**: ch1 shows a
   consistent polarity inversion across multiple datasets (Ruger and
   Glock), suggesting a fixed wiring characteristic at that array position
   rather than a one-off fault. `util_get_accoustic_param.py
   --check-channels` catches this automatically via correlation against
   the array's own consensus mean.

## Validated Reference Numbers

For the default config (200m range, Mach 2.5, Y_MISS=10m, tetrahedral
array L=0.30m, 5.56 NATO):
- Shockwave: `T_N ≈ 0.739 ms`, `dP ≈ 16.77 Pa` (identical at 1500m too —
  see Known Limitations #6)
- Muzzle blast @ 200m: `dP ≈ 9.99 Pa`, `t_pos ≈ 2.04 ms`
- Muzzle blast @ 1500m: `dP ≈ 1.33 Pa`, `t_pos ≈ 3.98 ms`
- `Dt` @ 200m: `≈ 324 ms`; @ 1500m: `≈ 2.60 s`
- True SW TDOA (pairs 0-1, 0-2, 0-3): `[-467.0, +52.4, +402.1] µs`
- True MB TDOA (pairs 0-1, 0-2, 0-3): `[-25.2, -424.1, +449.0] µs`

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
- Extending atmospheric absorption from single-characteristic-frequency
  peak attenuation to full waveform reshaping (real long-range blast
  propagation preferentially strips high-frequency content from the sharp
  leading edge, rounding the N-wave shape — not currently modeled).
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
