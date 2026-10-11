# Gunshot Locator — Timing Budget

- **Target platform:** Terasic DE10-Standard (Cyclone V SX 5CSXFC6D6F31C6: FPGA fabric + dual-core Cortex-A9 HPS)
- **Reference product:** Metravib PILAR-V
- **Derived from:** `python/det_config.ini`, `python/gen_config.ini`, `python/project.md`, `verilog/`, `reference/Pilar Accoustic Gunshot 2004.pdf`
- **Version:** 2026-10-09 (v2: measured values from Milestones A and B; v1 2026-10-08 was the pre-implementation estimate)

Items marked **(verify)** are assumptions, not sourced figures. Values marked **measured** come from the bit-exact simulations in `fft2048/`, `xspec_phat/`, `gcc_engine/` (Milestone A) and `stream/` (Milestone B).

---

## 1. Summary

| Item | Value |
|---|---|
| Sample rate | 100 kHz per channel, 4 channels, simultaneous sampling |
| TDOA error budget | ≤ 10 µs per mic pair (95 %) → ≈ 0.9° azimuth. **Measured:** fixed-point part ≤ 0.2 µs; bearing ≤ 0.75° over 50–1500 m (§4.3) |
| Δt (SW → MB) error budget | ≤ max(0.5 ms, 0.5 % of Δt) (v1: 0.5 ms; §4.4). **Measured** worst: 0.77 ms at 1500 m = 0.03 % |
| FFT size | 2048 for both paths: SW at 100 kHz, MB decimated to 20 kHz |
| FPGA compute per event | **Measured** 65.2k cycles = 0.65 ms @ 100 MHz (v1 estimate 88k; budget 2 ms) |
| Processing after MB window closes | FPGA 0.65 ms (measured) + processor ≈ 67 ms (estimate, Milestone C) |
| Burst capability | **Measured** 900 rpm: 8 of 8 shots, 2.0 % engine load, no late jobs |
| Range coverage | **Measured** in simulation: 50–1500 m (§10) |
| Range accuracy | Dominated by bullet drag, not timing: a drag-aware range solver is required (§4.4) |

**Partition rule.** Anything that touches samples or timestamps runs in the FPGA (hard real time). The HPS only handles per-event math (soft real time). Linux scheduling jitter can therefore delay a report but never corrupt a measurement. All event times (t_SW, t_MB) come from the FPGA sample counter; the HPS never timestamps.

---

## 2. Performance targets

| Parameter | Target | Source |
|---|---|---|
| Azimuth | ± 2° | PILAR 2004 paper, Table 1 |
| Elevation | ± 5° | PILAR 2004 paper, Table 1 |
| Range | ± 10 % ≤ 200 m, ± 20 % ≤ 600 m, ± 30 % ≤ 1200 m | PILAR 2004 paper, Table 1 |
| Processing latency after MB arrival | ≤ 0.5 s | **(verify)** against current PILAR-V datasheet |
| Automatic fire | 900 rpm (one shot / 66.7 ms) | **(verify)** |

---

## 3. Reference parameters (from the repo)

| Parameter | Value | Source |
|---|---|---|
| Array | Tetrahedron, edge L = 0.30 m | `det_config.ini [geometry]` |
| Speed of sound c | 343 m/s | `det_config.ini [physics]` |
| Max pair TDOA | L / c = 875 µs | derived |
| SW window | 3 ms pre + 15 ms post | `det_config.ini [tdoa]` |
| MB window | 6 ms pre + 80 ms post | `det_config.ini [tdoa]` |
| Max lag | 1 ms | `det_config.ini [tdoa]` |
| Sub-sample interpolation | 16× | `det_config.ini [tdoa]` |
| MB search range | dt_min 10 ms … dt_max 5 s | `det_config.ini [detection]` |
| Default-case Δt | 324 ms @ 200 m, 2.60 s @ 1500 m | `project.md` |
| MB positive phase t_pos | 2.0 ms @ 200 m, 4.0 ms @ 1500 m | `project.md` |

---

## 4. Timing-accuracy budget

### 4.1 DOA sensitivity of the array

Plane-wave model: τ_ij = (r_j − r_i)·(−u) / c. The slowness vector is solved by least squares over all 6 pairs and then normalised, so c cancels out of the direction. Monte Carlo: 20,000 trials per direction, azimuth swept 0–110°, worst-case 95th-percentile error (script in Appendix A).

| σ per pair | Azimuth 95 % (el 0°) | Azimuth 95 % (el 20°) | Elevation 95 % |
|---|---|---|---|
| 1 µs | 0.09° | 0.10° | 0.09° |
| 2 µs | 0.18° | 0.20° | 0.18° |
| 5 µs | 0.46° | 0.49° | 0.46° |
| 10 µs | 0.91° | 0.98° | 0.92° |
| 20 µs | 1.84° | 1.95° | 1.83° |

**Rule of thumb:** ≈ 0.09–0.10° per µs of pair TDOA error. Elevation (± 5°) is never the binding constraint; azimuth is.

The Monte Carlo treats pair errors as independent and random. Systematic errors (mic position, AFE phase) don't average down the same way, so treat the allocations below as combined bias + random.

### 4.2 Azimuth error allocation

| Source | Allocation (95 %) | Notes |
|---|---|---|
| TDOA timing (this document) | ≈ 0.9° | ≤ 10 µs per pair |
| Vehicle heading reference (GNSS/IMU) | ≈ 1° | adds directly to azimuth on a moving platform |
| Wind / refraction | ≈ 0.8° | 5 m/s crosswind → atan(5/343) = 0.84°; reducible with an anemometer |
| **RSS** | **≈ 1.6°** | inside ± 2° |

### 4.3 TDOA error allocation (per pair, 95 %)

| Contributor | Allocation | Requirement |
|---|---|---|
| GCC-PHAT estimator (noise, reverb) | 6 µs | fs = 100 kHz, 16× sub-sample refinement; Python already reaches single-digit µs on clean signals |
| AFE phase mismatch between channels | 3 µs | see 4.3.1 |
| Mic position tolerance | 3 µs | ± 1 mm per capsule (1 mm = 2.92 µs) |
| Fixed-point FFT / PHAT | 1 µs | bit-true comparison against `gs_det_tdoa.py` |
| ADC inter-channel skew | 0.5 µs | simultaneous-sampling ADC (see 4.3.2) |
| Sample-clock jitter | ≈ 0 | CONVST generated from FPGA PLL |
| **RSS** | **7.4 µs** | ≈ 0.7° azimuth; 2.6 µs margin to the 10 µs target |

**Measured (v2).**
- **Fixed-point FFT / PHAT:** ≤ 0.2 µs against a double-precision reference, inside its 1 µs allocation. One step of the 1/256-sample measurement grid at 20 kHz.
- **Total per-pair error on generator scenarios:** 2–7 µs on SW pairs, 4–12 µs on MB pairs. The estimator dominates.
- **MB pairs exceed 10 µs in some cases:**
  - 12.4 µs at 50 m
  - 12.2 µs at 1500 m / 0.2 Pa
  - one 101 µs outlier at 1500 m / 0.5 Pa (13 dB blast SNR)

  The least-squares fit over 6 pairs keeps bearing within 0.75°, except 2.0° elevation in the 0.5 Pa case, which is still inside ±5°. The processor should reject outlier pairs (closure check τ_ab + τ_bc − τ_ac ≈ 0).
- **Mic positions off by up to 1 mm** (generator only, detector unaware): bearing error +0.5° elevation, consistent with the 3 µs allocation.

#### 4.3.1 AFE phase mismatch

**High-pass (coupling) corner.** A first-order HP at f_c leads the signal by τ = atan(f_c / f) / (2πf). The MB carries most of its energy below 500 Hz, so this mainly hits the MB TDOA.

| f_c | τ @ 300 Hz | ± 20 % f_c spread | τ @ 500 Hz | ± 20 % f_c spread |
|---|---|---|---|---|
| 16 Hz (current `devices/readme.md`) | 28.3 µs | ± 5.6 µs | 10.2 µs | ± 2.0 µs |
| 2 Hz (recommended) | 3.5 µs | ± 0.7 µs | 1.3 µs | ± 0.25 µs |

A ± 20 % spread is realistic for X7R or electrolytic coupling caps. Either drop the corner to ≤ 2 Hz (film caps where possible) or calibrate each channel's response.

**Anti-alias filter.** Low-frequency group delay is Σ 1 / (Q_i · 2πf_0) over the filter sections.

| Filter | Group delay | 1 % mismatch | 5 % mismatch |
|---|---|---|---|
| 2nd-order Butterworth, 7 kHz (current `devices/readme.md`) | 32.2 µs | 0.32 µs | 1.6 µs |
| 4th-order Butterworth, 25 kHz (recommended) | 16.6 µs | 0.17 µs | 0.83 µs |

**Gain stage.** MCP6002 has 1 MHz GBW. At gain 11 the closed-loop bandwidth is ≈ 90 kHz, a delay of ≈ 1.8 µs. A ± 30 % part-to-part GBW spread gives ≈ ± 0.5 µs. An op-amp with ≥ 10 MHz GBW makes this negligible.

**Capsule.** The electret capsule and its vent have their own low-frequency rolloff, which is unknown. A one-time per-channel calibration against a common source is recommended (the Akman 2017 thesis in `/reference` describes a TDoA-based calibration of system parameters).

#### 4.3.2 ADC skew

- **Simultaneous ADC:** aperture matching is in the ns range, so skew is negligible.
- **Multiplexed ADC (onboard):** channel i lags channel 0 by i × t_conv. The `project.md` simulation shows ~2–7 µs. This is deterministic, so it can be subtracted as a constant bias, but it must be calibrated to < 0.5 µs. The 12-bit dynamic range is a separate problem (Section 9).

### 4.4 Range and Δt budget

PILAR range formula (eq. 4 of the 2004 paper), which needs no Mach number:

    d = c·Δt / (1 − m·s)

Here m is the MB wave vector, s is the SW wave vector, and α is the angle between them. The sensitivity is:

    δd/d ≈ δc/c ⊕ δΔt/Δt ⊕ δα·cot(α/2)

Angle term with δα = 1.3° (m and s both at the 10 µs limit, combined):

| α | Angle term | Δt @ 50 m | Δt @ 200 m |
|---|---|---|---|
| 90° | 2.3 % | 146 ms | 583 ms |
| 60° | 3.9 % | 73 ms | 292 ms |
| 45° | 5.5 % | 43 ms | 171 ms |
| 30° | 8.5 % | 20 ms | 78 ms |
| 25° | 10.2 % | 14 ms | 55 ms |

Below α ≈ 25° the formula is ill-conditioned regardless of hardware: the shooter is firing nearly along the line to the array.

**Δt allocation (v2): ≤ max(0.5 ms, 0.5 % of Δt).** The 0.5 ms floor keeps < 3 % at 50 m for α ≥ 30°. At long range a fixed 0.5 ms is needlessly strict: Δt enters range as δd/d = δΔt/Δt, so 0.5 % of Δt stays about 60× inside the ±30 % spec. v1 used a flat 0.5 ms, which the 1000–1500 m scenarios "failed" with range errors of only 0.03 %.

Constant delays are removed in the FPGA: `fe_cond` reports `x_ts` (4-sample latency removed) and `mb_decim` reports `y_ts` (31-sample group delay removed), so SW and MB timestamps compare directly.

**Measured MB onset bias vs range** (weaker, smoother blasts are detected later; SW onsets are within −50 … 0 µs):

| Range | 50 m | 100 m | 200 m | ~430 m | 600 m | 800 m | 1000 m | 1200 m | 1500 m |
|---|---|---|---|---|---|---|---|---|---|
| MB onset error | −0.03 ms | −0.08 ms | +0.09 ms | +0.31 ms | +0.39 ms | +0.49 ms | +0.53 ms | +0.65 ms | +0.75 ms |
| As % of Δt | 0.05 % | 0.05 % | 0.03 % | 0.04 % | 0.04 % | 0.04 % | 0.03 % | 0.03 % | 0.03 % |

**Bullet drag dominates the range budget (v2 finding).** With the generator's drag model (`python/drag_model.md`), the constant-speed formula above **overestimates** range even with perfect inputs:

| Calibre | 200 m | 400 m | 600 m | 800 m | 1000 m | 1500 m |
|---|---|---|---|---|---|---|
| 5.56 | +8 % | +21 % | +48 % | +118 % | no SW | no SW |
| 7.62 | +8 % | +19 % | +41 % | +91 % | +270 % | no SW |
| 12.7 | +3 % | +7 % | +11 % | +16 % | +24 % | +57 % |

A decelerating bullet delivers its shockwave later, which shrinks Δt. Range therefore needs a drag-aware solver with the **detector's own** ballistic table (v₀ and decay length per calibre class, from independent data; never the generator's library). That is Milestone C software. Timing errors (above) are now a minor term.

**Speed of sound.** c changes by ≈ 0.17 % per °C. A ± 2 °C air-temperature sensor costs ≈ 0.35 % of range.

---

## 5. Sample rate and window sizing

Reasons for fs = 100 kHz:

- It matches `gen_config.ini [adc] fs`, so simulated results transfer directly to hardware.
- SW energy extends to ~10 kHz. At 100 kHz the anti-alias filter can sit at ~25 kHz with a gentle 4th-order slope, which keeps its group delay and mismatch small (4.3.1).
- Akman 2018 reports 51.2 kHz as sufficient, so 100 kHz gives 2× margin.
- The raw data rate is only 4 × 100 kS/s × 2 B = 800 kB/s.

| Path | Rate | Window | Samples | + max lag | NFFT |
|---|---|---|---|---|---|
| SW | 100 kHz | 3 + 15 ms | 1800 | 100 | 2048 |
| MB (÷ 5) | 20 kHz | 6 + 80 ms | 1720 | 20 | 2048 |

Decimating the MB path by 5 is safe because MB energy is below ~1 kHz (PILAR: < 500 Hz at hundreds of metres). One FFT engine size (2048) serves both events. Avoiding circular wrap within ± max_lag requires NFFT ≥ window + max lag, which both paths satisfy.

**Masks (v2).** The windows stay at 18 ms (SW) and 86 ms (MB). Overlapping events are cut out instead of shortening windows:
- **A shockwave inside an MB window** (burst fire): −1 … +3 ms around its onset.
- **A later blast inside an SW window** (short range): from 1 ms before it to the window end.

The masked samples are zeroed on all 4 channels by `gcc_loader`.

---

## 6. Block budgets

**Clock assumption:** 100 MHz fabric clock from a PLL (the board oscillator is 50 MHz). At 50 MHz every compute time below doubles and still fits.

### 6.1 FPGA — stream rate (hard real time) — measured

One 4-channel sample frame arrives every 10 µs = 1000 clock cycles.

| Block | Cycles used | Latency | Notes |
|---|---|---|---|
| `adc_if` (Milestone D) | ≤ 500 per frame (budget) | conversion + readout ≤ 5 µs | Replaced by `replay_src` in simulation |
| `timebase` | — | 0 | 64-bit frame index |
| `fe_cond` | 49 per frame | 4 samples, constant | 9-tap delay trim, DC blocker; `x_ts` corrected |
| `sw_det` | 6 per frame (+4 per noise block) | candidate at the crossing sample; confirmed 1 ms later; record 2.56 ms later | N-wave qualification |
| `mb_decim` | 263 per 5 frames | 31 samples (0.31 ms), constant | `y_ts` corrected; one sequential MAC |
| `mb_det` | 9 per 5 frames (+32 scan, +12 thresholds) | trigger on 16-sample energy; held ≤ 1 ms during SW qualification | onset back-tracking; delay lines in RAM |
| `sw_ring` / `mb_ring` | 1 write per frame / per 5 frames | — | 41 ms / 205 ms of history |
| `evt_sched` | every cycle | release at window end; dispatch within a few cycles | deadline margin 16 samples |

### 6.2 FPGA — event rate (per SW or MB window) — measured

| Step | Cycles | @ 100 MHz |
|---|---|---|
| `gcc_loader` (ring → Z0, Z1, masks) | 2,052 | 21 µs |
| 2 × FFT (`fft2048`, DIT, bit-reversed input) | 2 × 11,321 | 226 µs |
| `xspec_phat` | 6,167 | 62 µs |
| 3 × IFFT | 3 × 11,321 | 340 µs |
| `peak_search` | 406 (SW) / 86 (MB) | 4 / 1 µs |
| **Job total** | **65,239 (SW) / 64,919 (MB)** | **0.65 ms** |

v1 estimated 88k cycles. The FFT avoids the bit-reverse pass (producers write bit-reversed), the loader writes both packed buffers per cycle, and peak search reads all three correlation buffers in parallel. Two back-to-back jobs: 1.30 ms.

### 6.3 FPGA resources — measured (yosys, indicative)

| Resource | Milestone B, area-optimised | (before optimisation) | Available |
|---|---|---|---|
| LUTs | ≈ 19,800 (≈ 10–13k ALMs, 25–32 %) | ≈ 34,300 | 41,910 ALMs |
| FFs | ≈ 8,500 | ≈ 13,300 | 166,036 |
| M10K | 133 (24 %) | 131 | 553 |
| 18×18 multipliers | 32 (≈ 16 DSP blocks) | 66 | 112 DSP blocks |

v1 estimated about 115 M10K, 48 DSP and about 9k ALMs. The optimisation moved delay lines and histories into RAM, shared multipliers, and serialised the decimator, with byte-identical behaviour (`fpga_blocks.md` §5).

### 6.4 FPGA → HPS transfer

| Data | Size per event | Path | Budget |
|---|---|---|---|
| Correlation lags around peak — SW | 6 × 201 × 32 bit = 4.8 kB | LW H2F bridge, CPU reads | 1 ms |
| Correlation lags around peak — MB | 6 × 41 × 32 bit ≈ 1 kB | LW H2F bridge | 0.2 ms |
| Event metadata (t_SW / t_MB, peak values, flags) | < 64 B | LW H2F bridge | — |
| Raw stream (optional, for logging) | 800 kB/s continuous | F2SDRAM DMA into HPS DDR3 | background |

The FPGA raises one F2H interrupt per completed event.

### 6.5 HPS — per event (soft real time)

| Block | Budget | Notes |
|---|---|---|
| IRQ → user space | 10 ms worst case | UIO on a standard kernel, typically < 0.2 ms; PREEMPT_RT not required |
| Read lags from FPGA | 1 ms | see 6.4 |
| 16× windowed-sinc peak refinement | 1 ms | 6 pairs × 33 fine lags × 32 taps ≈ 6.3k MAC |
| SW: trajectory direction s | 1 ms | 6 × 3 least squares |
| MB: bearing m, range, Mach / classification | 3 ms | `gs_det_shooter_locator.py` + `gs_det_classify_bullet.py` logic |
| SW ↔ MB association | 1 ms | queue of pending SW events |
| Report / network / display | 50 ms | |
| **Total** | **≈ 67 ms** | |

The windowed-sinc refinement replaces the 16× zero-padded IFFT in `gs_det_tdoa.py`. It is valid because the correlation is band-limited. With ± max_lag samples available, a 32-tap kernel centred on the coarse peak stays well away from the edges.

---

## 7. End-to-end timeline

Times from SW arrival at the array (t = 0). FPGA entries are **measured** (Milestone B, 200 m scenario); processor entries are still estimates (Milestone C).

| Event | 200 m (Δt ≈ 324 ms) | 1500 m (Δt ≈ 2.56 s) |
|---|---|---|
| SW candidate (first crossing) | 0 (onset −50 … 0 µs) | 0 |
| SW confirmed (N-wave qualification) | 1 ms | 1 ms |
| SW window closed | 15 ms | 15 ms |
| SW TDOAs (FPGA) | 15.7 ms | 15.7 ms |
| Early alert: shot detected + SW direction (processor estimate) | ≈ 85 ms | ≈ 85 ms |
| MB arrives | 324 ms | 2.56 s |
| MB detected (onset bias) | +0.09 ms | +0.75 ms |
| MB window closed (+ 80 ms) | 404 ms | 2.64 s |
| MB TDOAs (FPGA) | 404.7 ms | 2.64 s |
| Final alert: bearing + range (processor estimate) | ≈ 475 ms | ≈ 2.71 s |

- The FPGA adds 0.65 ms per window. Everything else is acoustic propagation, the window lengths, and the processor (estimated ≈ 70 ms).
- A subsonic bullet produces only the MB events: bearing at ≈ MB + 81 ms, no range from Δt.

---

## 8. Burst fire (900 rpm) — measured

8 shots, 6667 samples apart, 200 m (Milestone B `burst` scenario):

- **Detection:** 8 of 8 SW and 8 of 8 MB events, every block bit-exact, no late jobs or overflow.
- **Compute:** 2 jobs × 0.65 ms per 66.7 ms ≈ 2.0 % of the engine.
- **Window overlap:** at 200 m, shot k+5's shockwave lands 9.5 ms after shot k's blast, inside its MB window. Without a mask that gave a 10° bearing error. **Implemented:** SW masks in MB windows (`evt_sched` + `gcc_loader`), bringing it to ≤ 0.33°. 5 masks were applied in the burst. Trimming the MB window (v1's suggestion) can't fix it at this range.
- **MB detector under burst:** `mb_det` arms once settled and stays armed (latched), and it runs continuously. Without the latch, 1 of 8 blasts was missed.
- **Association** (up to `dt_max_s` / 66.7 ms ≈ 75 pending SW events) is still a processor task (Milestone C).

---

## 9. Repo items that conflict with this budget — status

- [ ] `devices/readme.md`: 16 kHz sampling and 7 kHz LPF → 100 kHz and a 4th-order ~25 kHz anti-alias filter. *(Milestone D hardware)*
- [ ] `devices/readme.md`: 1 µF / 10 kΩ coupling (16 Hz corner) → ≤ 2 Hz, or per-channel calibration (4.3.1). *(D)*
- [ ] `devices/readme.md`: MCP6002 (1 MHz GBW) → ≥ 10 MHz GBW op-amp (4.3.1). *(D)*
- [ ] Onboard ADC: multiplexed 12-bit → external 16-bit simultaneous-sampling ADC. Milestone B also showed a ~40× level range between a 50 m and a 1500 m blast. *(D)*
- [ ] `python/project.md`: ADS8688 is listed as simultaneous-sampling, but it is multiplexed. Simultaneous options: AD7606B, ADS8588S.
- [x] `gcc_phat_top.v` parameters, `peak_detector_tdoa.v` `tdoa_calc`, `fft1024`: superseded by `fft2048/`, `xspec_phat/` and `gcc_engine/` (old RTL moved to `old/`).
- [ ] `python/gs_det_tdoa.py`: add the windowed-sinc refinement and the PHAT band masks (SW 20–512, MB 10–205 bins), so Python and FPGA use the same estimator.
- [x] `python/det_config.ini` `mb_window_post_s = 0.080`: kept; burst overlap is handled by masks instead (§8).
- [x] Constant delays before Δt: compensated in the FPGA (`x_ts`, `y_ts`). *(Python side still open)*
- [ ] `python/gs_det_signal_prepare.py`: add a minimum noise floor (`noise_min`). With ADC noise below 1 LSB its thresholds collapse and onsets come out up to 322 ms early.
- [ ] `python/det_config.ini [bullet_reference_library]`: give the detector its own bullet values (they currently equal the generator's), per the separation rule in `SIMULATION.md`.
- [ ] Processor (Milestone C): drag-aware range solver with a detector-owned ballistic table (§4.4); blast-only (subsonic) path; outlier-pair rejection (§4.3).

---

## 10. Range coverage (v2, measured in simulation)

| Range | Status | Limiting factor |
|---|---|---|
| < 50 m | Out of scope by decision | Below ~20 m the blast falls inside the 10 ms MB blanking. 25 m passes with defaults; 15 m passes with 3 ms blanking and 2 ms settling, untested against ground reflections. |
| 50 m | ✓ | Blast front as steep as a shockwave → N-wave symmetry classification |
| 50–1500 m | ✓ | MB onset bias grows to 0.75 ms (0.03 % of Δt) |
| 1500 m | ✓ at 0.05 / 0.2 / 0.5 Pa noise | 0.5 Pa (13 dB blast SNR) is the edge: one outlier MB pair, 2.0° elevation error |
| Subsonic at the array | ✓ blast-only event | 5.56 beyond ≈ 900 m, 7.62 beyond ≈ 1000 m (generator drag model): bearing only, range needs another source |
| Mismatch (detector unaware) | ✓ | Mics off by up to 1 mm; different SW / MB shapes; 40 °C air |

---

## Appendix A — DOA sensitivity script

```python
import numpy as np

L, c = 0.30, 343.0
# Tetrahedron: base triangle in z = 0, apex up, centred on origin
R, h = L / np.sqrt(3), L * np.sqrt(2 / 3)
P = np.array([[R * np.cos(a), R * np.sin(a), 0] for a in np.deg2rad([90, 210, 330])]
             + [[0, 0, h]])
P -= P.mean(0)
pairs = [(i, j) for i in range(4) for j in range(i + 1, 4)]
D = np.array([P[j] - P[i] for i, j in pairs])

def doa_error(az, el, sigma_pair, n=20000, rng=np.random.default_rng(0)):
    u = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    tau = -(D @ u) / c
    T = tau + rng.normal(0, sigma_pair, (n, 6))
    S = np.linalg.lstsq(D, T.T, rcond=None)[0].T        # slowness vectors
    U = -S / np.linalg.norm(S, axis=1, keepdims=True)    # unit directions
    az_hat = np.arctan2(U[:, 1], U[:, 0])
    el_hat = np.arcsin(np.clip(U[:, 2], -1, 1))
    d_az = np.rad2deg(np.angle(np.exp(1j * (az_hat - az))))
    d_el = np.rad2deg(el_hat - el)
    return np.percentile(np.abs(d_az), 95), np.percentile(np.abs(d_el), 95)

for s_us in [1, 2, 5, 10, 20]:
    for el in [0, 20]:
        res = [doa_error(np.deg2rad(a), np.deg2rad(el), s_us * 1e-6) for a in range(0, 120, 10)]
        print(f"sigma={s_us:>2} us  el={el:>2}  az95={max(r[0] for r in res):.2f}  "
              f"el95={max(r[1] for r in res):.2f} deg")
```

---

## Appendix B — Key formulas

| Quantity | Formula |
|---|---|
| Max pair TDOA | L / c |
| HP phase lead | τ = atan(f_c / f) / (2πf) |
| LP group delay at DC | Σ 1 / (Q_i · 2πf_0) |
| Mic position → TDOA | δτ = δx / c |
| Wind → azimuth | δθ ≈ atan(v_cross / c) |
| Range (PILAR eq. 4) | d = c·Δt / (1 − m·s) |
| Range sensitivity | δd/d ≈ δc/c ⊕ δΔt/Δt ⊕ δα·cot(α/2) |
| FFT cycles (current radix-2 core) | 2N + (N/2)·log₂N |
| Sample frame budget | f_clk / fs = 1000 cycles @ 100 MHz, 100 kHz |
