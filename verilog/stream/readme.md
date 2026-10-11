# stream — continuous detection chain (Milestone B)

Turns a continuous 4-channel 100 kHz recording into timestamped SW / MB events and GCC-PHAT results, in real time, with the Milestone A engine.

```
replay_src -> timebase -> fe_cond -+-> sw_ring --------------------+
                                   +-> sw_det ------+               |
                                   +-> mb_decim -+-> mb_ring -------+
                                                 +-> mb_det -+      |
                                 sw_det / mb_det -> evt_sched -> gcc_engine
```

- **Status: Milestone B passed over the full range, 50–1500 m** (Verilator 5.020), and **area-optimised** (−42% LUTs, −52% multipliers) with byte-identical behaviour; see Resources. Fifteen scenarios from your generator:
  - single shots at 50, 200, ~430, 600 and 1500 m (three noise levels at 1500 m)
  - three shots with **bullet drag** (5.56 and 7.62 at 1000 m, 12.7 at 1500 m), one of them subsonic, so blast-only
  - four **mismatch** scenarios, where the generator's world differs from the detector's assumptions and the detector is not re-tuned
  - an 8-shot 900 rpm burst
- **Results:** every block is bit-exact against its Python model, every shot gives exactly one SW and one MB event (MB only when subsonic), and no run has late jobs, overflows or underflows. Not yet through Quartus.
- **Depends on:** `../fft2048`, `../xspec_phat`, `../gcc_engine`, your generator in `python/` (found automatically), and **Verilator**.
- **Version:** 2026-10-11 (area-optimised)

---

## Files

| File | Block | Spec (`fpga_blocks.md`) |
|---|---|---|
| `replay_src.v` | Frame source: replay FIFO (or ADC pass-through) + frame cadence | 2.2 |
| `timebase.v` | 64-bit sample index | 2.3 |
| `fe_cond.v` | Offset, gain, 9-tap fractional-delay trim, DC blocker | 2.4 |
| `sw_det.v` | Shockwave detector: candidate, N-wave qualification, features | 2.5 |
| `mb_decim.v` | ÷5 decimator, 63-tap linear phase, group delay removed from timestamps | 2.6 |
| `mb_det.v` | Muzzle-blast detector: sliding energy, latched settling, held events, onset refinement | 2.7 |
| `evt_sched.v` | Event ids, window release, SW/MB masks, dispatch, deadlines | 3.1 |
| `gsl_core.v` | Milestone B top: everything above + `cap_ring` ×2 + `gcc_engine` | — |
| `gen_stream_coef.py` | Generates `fe_fdelay.hex` and `decim63.hex` | — |
| `tb_gsl_core.v` | Testbench: plays a recording, logs streams / events / results | — |
| `milestone_b.py` | `prepare` (scenarios), `check` (bit-true models + accuracy) | — |
| `Makefile` | `make` = prepare → build → simulate all → check | — |

## Quick start

```bash
sudo apt install verilator          # once (macOS: brew install verilator)
cd GSL/verilog/stream
make                                 # ~25 min: prepare, 15 simulations (~8 min), check (~15 min)
```

- **`make FRAME_DIV=1000`:** real-time cadence (1000 clocks per sample at 100 MHz). The default of 250 gives identical results 4× faster, and is a harsher scheduling test, since the engine gets 4× less time per sample.
- **`make SIM=icarus`:** uses Icarus. Results are identical (checked bit-for-bit), but it is about 125× slower.
- **`python3 milestone_b.py check --only burst`:** checks one scenario. `--set key=value` overrides a configuration value.

---

## What each block does

| Block | Per-sample work | Multipliers | Notes |
|---|---|---|---|
| `replay_src` | Pop 1 frame per tick | 0 | FIFO 1024 frames; `rp_underflow` sticky = broken timeline |
| `fe_cond` | 49 cycles (4 ch × 12 + 1) | 1 | dtrim ±1 sample in 1/128 steps, exact at 0; DC corner 1.94 Hz; latency 4 samples |
| `sw_det` | 6 cycles | 2 (shared) | Block noise floor, 10× RMS candidate, 1 ms N-wave qualification |
| `mb_decim` | 263 cycles every 5 samples | 1 | 0.02 dB ripple to 4 kHz, < −71 dB from 16 kHz; `y_ts` = filter-centre time |
| `mb_det` | 9 cycles every 5 samples (+32 scan) | 2 | 16-sample energy, 20× trigger, 9× / 5 ms settle (latched), onset back-tracking; delay lines in RAM |
| `evt_sched` | Every cycle | 0 | Release when the window is complete; up to 4 masks per job |

## Configuration (ports on `gsl_core`; defaults in `milestone_b.py`)

| Group | Setting | Default | Meaning |
|---|---|---|---|
| fe | `gain`, `offset`, `dtrim`, `dcblk_en` | 16384, 0, 0, 1 | Per channel: Q2.14 gain, s16 offset, 1/128-sample trim |
| sw | `thr_mult2`, `noise_min` | 100, 16 | Candidate at (10× RMS)²; minimum noise mean square (LSB²) |
| sw | `holdoff`, `dt_min`, `feat_win` | 2000, 1000, 256 | Samples (20 ms, 10 ms MB blanking, 2.56 ms features) |
| sw | `sym`, `tn_max`, `pk_mult` | 32, 60, 9 | N-wave tests: swing ratio ≥ 0.5, peaks ≤ 0.6 ms apart, peak ≥ 3× threshold |
| mb | `thr_mult`, `settle_mult` | 20, 9 | Energy ratios (Python: 20× noise², settle at 3× RMS) |
| mb | `dwell`, `holdoff`, `settle_timeout`, `noise_min` | 100, 400, 4000, 16 | 20 kHz samples (5 ms, 20 ms, 200 ms); LSB² |
| windows | `sw_pre/len/maxlag`, `mb_pre/len/maxlag` | 300/1800/100, 120/1720/20 | As `det_config.ini` |
| PHAT | `sw_band`, `mb_band`, floors | 20–512, 10–205, 0 | Recommended band masks |
| masks | `mb_mask_pre/post` | 20, 60 | SW masked in MB windows: −1 … +3 ms (20 kHz samples) |
| masks | `sw_mask_pre/post` | 100, 2000 | Later MB masked in SW windows: −1 ms … window end (100 kHz samples) |

---

## Milestone B results

| Scenario | Events SW / MB | Bit-exact | Δt error / budget | Worst TDOA error (SW / MB) | MB bearing error (az / el) |
|---|---|---|---|---|---|
| r50: 50 m (minimum) | 1 / 1 | ✓ | 16 / 500 µs | 4.7 / 12.4 µs | 0.00° / 0.00° |
| s1: 200 m | 1 / 1 | ✓ | 142 / 1,618 µs | 4.1 / 7.4 µs | −0.19° / −0.33° |
| s2: shooter (350, 250, 20) m, 0.05 Pa | 1 / 1 | ✓ | 357 / 3,599 µs | 4.6 / 5.7 µs | −0.30° / −0.22° |
| s3: 600 m, 0.02 Pa | 1 / 1 | ✓ | 400 / 4,919 µs | 5.3 / 5.0 µs | −0.19° / +0.25° |
| L1500_n05: 1500 m, 0.05 Pa | 1 / 1 | ✓ | 767 / 12,781 µs | 5.7 / 4.5 µs | −0.11° / +0.25° |
| L1500_n20: 1500 m, 0.2 Pa | 1 / 1 | ✓ | 717 / 12,781 µs | 6.4 / 12.2 µs | −0.28° / +0.75° |
| L1500_n50: 1500 m, 0.5 Pa | 1 / 1 | ✓ | 747 / 12,781 µs | 7.0 / 101.5 µs | +0.41° / +2.00° |
| burst: 8 shots, 900 rpm (200 m) | 8 / 8 | ✓ | 182 / 1,618 µs | 4.7 / 7.4 µs | ≤ 0.19° / ≤ 0.33° |
| d556_1000: 5.56 with drag, 1000 m (**subsonic**) | 0 / 1 | ✓ | n/a (no SW) | n/a / 10.2 µs | +0.14° / −0.58° |
| d762_1000: 7.62 with drag, 1000 m (Mach 1.12) | 1 / 1 | ✓ | 521 / 5,237 µs | 2.4 / 7.0 µs | +0.02° / +0.25° |
| d127_1500: 12.7 with drag, 1500 m (Mach 1.44) | 1 / 1 | ✓ | 774 / 10,048 µs | 6.2 / 11.5 µs | +0.60° / +0.50° |
| m_tol: mics off by up to 1 mm (generator only) | 1 / 1 | ✓ | 145 / 1,618 µs | 5.0 / 7.8 µs | −0.22° / +0.50° |
| m_sw: shockwave 1.8× longer, 0.5× amplitude | 1 / 1 | ✓ | 142 / 1,618 µs | 3.3 / 7.3 µs | −0.16° / −0.17° |
| m_mb: blast 0.4× duration (sharper) | 1 / 1 | ✓ | 142 / 1,618 µs | 4.7 / 5.6 µs | −0.08° / −0.17° |
| m_hot: air 40 °C, detector assumes 343 m/s | 1 / 1 | ✓ | 159 / 1,566 µs | 4.7 / 8.9 µs | −0.18° / −0.17° |

Also verified during development but not in the standard set: 25 m, 100 m, 800 m, 1000 m and 1200 m all pass. 15 m passes only with shorter timing; see "Range limits".

- **Separation:** bearing uses the detector's own geometry and speed of sound (`python/det_config.ini`); generator ground truth is used only to score. See the test rule in `../SIMULATION.md`.
- **Bit-exact** covers the x and y streams, every SW / MB record, every job (id, type, ring start, number of masks) and every engine result including lag windows.
- **Δt budget** = max(500 µs, 0.5% of the true SW→MB gap). 500 µs comes from `timing_budget.md` §4.4, where it is set by the 50 m case. At long range the same relative error matters far less: Δt error enters range as d_range/range = d_Δt/Δt, so 0.5% stays about 60× inside the ±30% spec. Consider updating §4.4 the same way.
- **1500 m at 0.5 Pa** (13 dB blast SNR) is the edge of the envelope: one of the six MB pairs is off by 101 µs. Bearing is still inside spec (±2° az, ±5° el), but the processor should reject outlier pairs with a closure check (τ_ab + τ_bc − τ_ac ≈ 0) in Milestone C.

---

## Range limits (design target: 50–1500 m)

| Range | What limits it | Status |
|---|---|---|
| < 50 m | SW→MB gap shrinks: 23 ms at 25 m, 8.5 ms at 15 m. Below ~20 m the blast falls inside the 10 ms MB blanking (`dt_min`, from your `det_config.ini`). | 25 m passes with defaults. 15 m passes with `sw_dt_min = 300` and `mb_dwell = 40` (3 ms / 2 ms). Your generator is free-field, so whether short blanking survives real ground reflections is untested. Out of scope by decision. |
| 50 m | Blast front as steep as a shockwave (steepness 0.21 vs 0.34) | Solved by N-wave symmetry classification |
| 50–1500 m | MB onset detected progressively late: +0.09 ms at 200 m, +0.39 ms at 600 m, +0.75 ms at 1500 m | Fine relative to the gap (≤ 0.03% at 1500 m) |
| 1500 m | Blast SNR: 2.24 Pa peak. 33 dB at 0.05 Pa noise, 13 dB at 0.5 Pa | Works; 0.5 Pa is the edge (see above) |
| > ~900 m (5.56), > ~1000 m (7.62) | **Subsonic bullets make no shockwave** (generator drag model, `python/drag_model.md`) | The FPGA handles it: blast-only event, bearing within 0.6°. Range needs another source (Milestone C) |
| Transonic (Mach 1.1–1.2) | Weak, slightly longer N-wave | Still detected and classified (7.62 at 1000 m, Mach 1.12) |

---

## Design additions beyond the Python detector

Each came out of a simulation, is configurable, and has a measured reason.

1. **N-wave qualification (`sw_det`).** In continuous operation the SW detector also fires on the muzzle blast. A candidate is confirmed as a shockwave only if, within 1 ms:
   - **Symmetry:** the smaller of the positive and negative swings is ≥ 0.5 × the larger. N-waves measure 0.89–1.04 at 15–1500 m; blasts measure 0.15–0.17.
   - **Duration:** the two peaks are ≤ 0.6 ms apart. N-waves measure 0.18 ms.
   - **Margin:** the peak is ≥ 3× the trigger threshold.

   An earlier front-steepness test failed below about 50 m, where blast fronts become as steep as shockwaves (0.21–0.30 vs 0.34). The duration and margin checks were added after a weak 1500 m blast in 0.2 Pa noise passed the symmetry test alone.
2. **Held MB events (`mb_det`).** An MB trigger during shockwave qualification is held. If the shockwave is confirmed, the event was the shockwave's own energy and is dropped, along with its holdoff. If the candidate is rejected, the event was a real blast and is released with its original onset.
3. **Latched MB arming (`mb_det`).** Once settled, the detector stays armed until it triggers. Without the latch, 1 of 8 burst blasts was missed when its front crossed the settle level one sample before the trigger level.
4. **Masks (`evt_sched` + `gcc_loader`).**
   - **SW masked in MB windows:** at 900 rpm a later shot's shockwave lands in a blast window. That gave a 10° bearing error; with the mask it's 0.13°.
   - **Later MB masked in SW windows:** at short range the blast lands in the shockwave window. That gave a 459 µs SW TDOA error at 15 m; with the mask it's 4.1 µs. Only blasts *after* the SW onset are masked, because in a burst an earlier shot's blast precedes the shockwave by 9.3 ms, and masking it would cut the shockwave.
5. **Noise floor minimum (`noise_min`).** The generator's ADC noise is below 1 LSB, so a pure 10× RMS threshold collapses to about 0.

## Findings

- **Bullet drag changes the software design (most important).** With the new generator drag model (`python/drag_model.md`):
  - **Subsonic bullets:** 5.56 is subsonic beyond about 900 m and 7.62 beyond about 1000 m, so there is no shockwave. The FPGA handles it (scenario `d556_1000`), but the processor needs a blast-only path: bearing, with range from another source.
  - **Range formula:** the constant-speed PILAR formula d = c·Δt / (1 − m·s) **overestimates** range with a slowing bullet, even with perfect inputs. For 5.56: +8% at 200 m, +48% at 600 m. For 12.7: +57% at 1500 m. Milestone C needs a drag-aware range solver.
  - **No RTL change needed:** the hardware results (events, TDOAs, bearing) are unaffected.
- **Your Python detector misfires on these recordings.** Its noise RMS ≈ 0 makes its thresholds ≈ 0, giving SW onsets 3–8 ms early and MB onsets up to 322 ms early. A minimum noise floor would align it with the hardware.
- **Absolute levels.** The generator's ADC stage scales each recording to its own peak (`peak_match`), so these tests can't reveal clipping. A real fixed-gain front end must hold a 50 m blast (≈ 4× the 200 m level) without clipping and still resolve a 2 Pa blast at 1500 m. That is roughly a 40× level range, comfortable for a 16-bit ADC but not for the onboard 12-bit one.

## Not covered by Milestone B

- **Real ADC and microphones** (Milestone D).
- **`csr`, `result_mbox` and the processor** (Milestone C).
- **The generator's `realistic` noise model,** which is unusable for any TDOA method (`gcc_engine/readme.md`, finding 1).

---

## Resources (indicative, yosys `synth_intel_alm -family cyclonev`)

| Block | LUTs before → after | FFs before → after | M10K | 18×18 multipliers before → after |
|---|---|---|---|---|
| `replay_src` | 292 | 106 | 7 | 0 |
| `fe_cond` | 3,992 → **2,174** | 1,282 → 638 | 3 → 4 | 6 → 4 |
| `mb_decim` | 2,410 → **1,056** | 449 → 325 | 3 | 16 → 2 |
| `sw_det` | 2,970 → 3,048 | 1,360 → 1,596 | 0 | 18 → 8 |
| `mb_det` | 10,828 → **2,546** | 3,677 → 1,463 | 0 → 1 | 14 → 6 |
| `evt_sched` | 8,238 → **5,142** | 4,001 → 2,029 | 0 | 0 |
| `gcc_engine` (loader, bank, FFT, xspec, peak) | 5,590 | 2,390 | 46 | 12 |
| `cap_ring` × 2 | — | — | 72 | 0 |
| **Total** | **34,320 → 19,848 (−42%)** | **13,265 → 8,547 (−36%)** | **131 → 133 / 553 (24%)** | **66 → 32 (≈ 16 / 112 DSP)** |

**Area optimisation (done).** About 10–13k ALMs, or 25–32% of the device, down from 40–55%. Behaviour is unchanged: the optimised RTL produces byte-identical logs to the original on all 15 Milestone B scenarios, plus a 15 m case with masks in both directions and `fe_cond` with non-zero trims, gains and offsets. It also matches in Icarus at 100 clocks per sample. The bit-true Python models need no changes.

| Block | Change |
|---|---|
| `mb_det` | 4 × 16 delay lines of squares moved from registers into a 64 × 36 RAM (circular buffers, cleared after reset); scan reads the RAM; one shared threshold multiplier |
| `evt_sched` | Job queue stores mask centres, not bounds (formed at dispatch); one set of overlap units for both job types; MB window time recomputed at dispatch; **queue depth 8 → 4** (occupancy never exceeds 1, even at 900 rpm) |
| `fe_cond` | Filter history in a 64 × 18 RAM read in step with the coefficient ROM; gain and taps share one multiplier |
| `mb_decim` | One multiply-accumulator, channels in sequence (263 cycles per output); history 64 → 128 slots so new samples can't overwrite taps still to be read at any frame cadence |
| `sw_det` | Two shared multipliers; the margin threshold is stored per channel. Saves DSP, not LUTs (its logic is mostly feature tracking) |

**Minimum frame cadence:** every block needs at least about 64 clocks per sample (`fe_cond` uses 49). Real time is 1000; simulation uses 250 (Verilator) or 100 (Icarus).

## Interface changes to Milestone A blocks

- **`gcc_engine`** gains `job_ts0`, `job_mask_en`, `job_mask_lo` and `job_mask_hi`, passed to the loader. `tb_gcc_engine.v` ties them off.
- **`gcc_loader`** zeros all 4 channels of window sample n when its time `ts0 + n·step` falls inside an enabled mask interval (step 5 for MB, 1 for SW).
- **Milestone A re-verified** after these changes: 10 jobs, bit-exact, PASSED.
