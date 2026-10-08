# gcc_engine — GCC-PHAT core (Milestone A)

One job = one 4-channel window (SW or MB) in, 6 coarse TDOAs plus lag windows out. This folder adds the last three Milestone A blocks and the core's top level:

| Block | Spec section (`fpga_blocks.md`) | Role |
|---|---|---|
| `cap_ring` | 2.8 | 4096 × 72-bit capture ring (used for both `sw_ring` and `mb_ring`) |
| `gcc_loader` | 3.2 | Ring → Z0 / Z1, packed, bit-reversed, zero-padded |
| `peak_search` | 3.5 | Coarse peak, peak value and sidelobe per pair; lag-window stream |
| `gcc_engine` | 3.1–3.5 | Sequencer + buffer-port mux; instantiates loader, bank, `fft2048`, `xspec_phat`, `peak_search` |

- **Status: Milestone A passed** (Icarus Verilog 12): 10 jobs, every stage bit-exact against the Python model, worst fixed-point TDOA error 0.195 µs, 0.65 ms per job. Not yet through Quartus.
- **Depends on:** `../fft2048`, `../xspec_phat`, and your generator in `../../python` (for the realistic jobs).
- **Version:** 2026-10-08

---

## Files

| File | Description |
|---|---|
| `gcc_engine.v` | GCC-PHAT core top: job sequencer, client-port mux, debug read port |
| `gcc_loader.v` | Window loader |
| `peak_search.v` | Peak search (2 passes over ± max_lag) |
| `cap_ring.v` | Capture ring RAM (2-cycle read latency) |
| `tb_gcc_engine.v` | Data-driven Milestone A testbench |
| `milestone_a.py` | Builds the jobs (`prepare`) and verifies the results (`check`) |
| `Makefile` | `make` = prepare → simulate → check |

## Quick start

```bash
cd verilog/gcc_engine
make            # ~3 min: prepare jobs, simulate 10 jobs, check
```

Needs numpy **and scipy** (your generator uses scipy, and so does the MB decimator stand-in). See `../SIMULATION.md`.

---

## Interface (`gcc_engine`)

| Port | Dir | Width | Description |
|---|---|---|---|
| `job_valid`, `job_ready` | in, out | 1 | Job handshake (accepted when `job_ready`) |
| `job_type` | in | 1 | 0 = SW (`sw_ring`), 1 = MB (`mb_ring`) |
| `job_id` | in | 16 | Returned in `res_id` |
| `job_start` | in | 12 | Ring address of window sample 0 (onset − pre) |
| `job_len` | in | 12 | Window length (1800 SW, 1720 MB) |
| `job_maxlag` | in | 7 | ± lag range (100 SW, 20 MB) |
| `cfg_{sw,mb}_band_lo/hi` | in | 11 | PHAT band mask per job type |
| `cfg_{sw,mb}_floor` | in | 7 | PHAT floor per job type (0 = off) |
| `sw_ring_raddr/rdata`, `mb_ring_raddr/rdata` | out / in | 12 / 72 | Capture-ring read ports |
| `res_valid` | out | 1 | Result strobe |
| `res_type`, `res_id` | out | 1, 16 | Job echo |
| `res_lag` | out | 6 × s8 | Coarse lag per pair, samples (lag > 0: second channel lags) |
| `res_pk` | out | 6 × s18 | Correlation value at the peak |
| `res_side` | out | 6 × u18 | Max \|r\| more than 2 samples from the peak |
| `res_exp_z`, `res_exp_r` | out | 2 × 5, 3 × 5 | BFP exponents (forward FFTs, IFFTs) |
| `res_err` | out | 1 | `fft_err` or buffer-port conflict during the job |
| `lw_we`, `lw_idx`, `lw_data` | out | 1, 8, 6 × s18 | Lag-window stream: one word per lag (−L … L), all 6 pairs |
| `dbg_re`, `dbg_sel`, `dbg_raddr`, `dbg_rdata` | in / out | — | Read any buffer while idle (testbench / processor debug) |

Pair order everywhere: (0,1) (0,2) (0,3) (1,2) (1,3) (2,3). Pair j sits at bits `[j*W +: W]`.

Peak-to-sidelobe ratio = `res_pk / res_side` is left to the processor (no hardware divider). A peak ≤ 2 marks a suppressed (dead) channel.

---

## Job timing (measured, 100 MHz)

| Step | Cycles |
|---|---|
| LOAD (`gcc_loader`) | 2,052 |
| FFT Z0, FFT Z1 | 2 × 11,321 |
| XSPEC | 6,167 |
| IFFT W0, W1, W2 | 3 × 11,321 |
| PEAK | 406 (SW) / 86 (MB) |
| Sequencer overhead | ~10 |
| **SW job** | **65,239 = 0.652 ms** |
| **MB job** | **64,919 = 0.649 ms** |
| Two jobs back to back | 129,838 = 1.298 ms |

Budget: 2 ms per job. At 900 rpm (2 jobs per 66.7 ms) the engine is busy about 2% of the time.

---

## Milestone A results

10 jobs, ring start addresses spread over the ring (including three that wrap around address 4095 → 0):

| # | Job | Bit-exact | Fixed-point error | Worst \|TDOA − truth\| | Bearing error (az / el) |
|---|---|---|---|---|---|
| 1 | Synthetic noise, integer delays | ✓ | 0.000 µs | 0.0 µs | — |
| 2 | Synthetic N-wave, az 37° el 10°, ring wrap | ✓ | 0.000 µs | 0.8 µs | −0.04° / +0.01° |
| 3 | Synthetic Friedlander, MB path, all bins | ✓ | 0.195 µs | 23.8 µs | +0.01° / +1.37° |
| 4 | Dead mic + band 20–500 + floor 27 | ✓ | 0.000 µs | 0.0 µs (live pairs); dead pairs flagged | — |
| 5 | s1 SW (200 m, miss 10 m) | ✓ | 0.039 µs | 4.7 µs | — |
| 6 | s1 MB | ✓ | 0.000 µs | 5.6 µs | −0.08° / −0.17° |
| 7 | s2 SW (shooter at 350, 250, 20 m; noise 0.05 Pa) | ✓ | 0.078 µs | 7.2 µs | — |
| 8 | s2 MB | ✓ | 0.195 µs | 3.9 µs | −0.25° / −0.14° |
| 9 | s3 SW (600 m, miss 25 m) | ✓ | 0.156 µs | 4.8 µs | — |
| 10 | s3 MB | ✓ | 0.195 µs | 5.0 µs | −0.19° / +0.25° |

- **Bit-exact** covers the loader addressing (rebuilt from the ring image), both forward FFTs and exponents, `xspec_phat` + the three IFFTs, and `peak_search` lags, peaks and sidelobes. The testbench also checks that the lag-window stream equals the correlation buffers.
- **Fixed-point error:** RTL correlations vs double-precision GCC-PHAT on the same windows and band, both refined with the planned windowed-sinc method on a 1/256-sample grid. 0.195 µs is one grid step at 20 kHz. Budget: 1 µs.
- **Realistic jobs** (5–10) come from your generator (`gs_gen_clean_signal` → `gs_gen_add_noise` → `gs_gen_apply_adc`, int16 ADC codes) with the recommended band masks (SW 20–512, MB 10–205). Bearing is the least-squares direction from the 6 MB TDOAs vs the true shooter direction.
- The `repo-truth` column in `milestone_a.py check` shows `gs_det_tdoa.gcc_phat` (all bins, 100 kHz) on the same spans for comparison; it is up to 22.6 µs off on MB pairs where the band-limited FPGA path is within 4 µs.

### What Milestone A does not yet cover

- **Onsets come from ground truth**, not from a detector: windows are placed around the generator's true arrival times. The detectors are Milestone B.
- **MB decimation is a Python stand-in** (63-tap linear-phase FIR, ÷5, group delay removed) until `mb_decim` exists (Milestone B).
- **Only the `simple` noise model** is used; see finding 1.
- Bearing errors exclude vehicle heading, wind and array-position errors (see `timing_budget.md` §4.2).

---

## Findings

1. **The generator's `realistic` noise model makes TDOA impossible.** `make_room_ir` gives a diffuse tail starting at sample 1 with about 40 dB more energy than the direct sound (its own docstring says so), and each mic gets an independent tail. Every channel then sees a different random waveform: float GCC-PHAT returns ≈ 0 lag for all pairs even with perfect windows, and `gs_det_signal_prepare.py` reports a 5 ms SW onset spread. Open-field gunfire is direct-path dominated. Consider adding a direct-to-reverberant ratio parameter and an initial delay gap before the tail.
2. **Band masks matter on the MB path.** Job 3 (all bins) has up to 24 µs TDOA error vs ≤ 5.6 µs for the band-limited realistic MB jobs. Keep SW 20–512 and MB 10–205 as defaults.
3. **PSR needs per-type thresholds.** The MB correlation peak is several samples wide at 20 kHz, so its PSR (pk / side, side excluding ±2 samples) is only 1.7–2.1 on the realistic MB jobs, while SW jobs give 2.7–5.4 and integer-delay noise 50–90. Either widen the exclusion zone for MB or threshold per type on the processor.

---

## Resources (indicative, yosys `synth_intel_alm -family cyclonev`)

| Block | M10K | 18×18 multipliers | LUTs | FFs |
|---|---|---|---|---|
| `gcc_engine` (whole core: loader, bank, FFT, xspec, peak) | 46 | 12 (= 6 DSP blocks) | ≈ 5,600 | ≈ 2,400 |
| `cap_ring` (each; 2 needed) | 36 | 0 | — | — |
| **Total GCC-PHAT core + 2 rings** | **118 / 553 (21%)** | **6 / 112 DSP** | **≈ 5,600** | **≈ 2,400** |

## Changes vs `fpga_blocks.md`

| Spec | Implementation | Reason |
|---|---|---|
| Loader ≤ 4,096 cycles | 2,052 | Bank v2 writes Z0 and Z1 in the same cycle |
| `pk_psr` Q8.8 | `res_side` (u18); PSR computed by the processor | No hardware divider |
| Peak search ≤ 1,206 cycles | 406 (SW) / 86 (MB) | All three W buffers read in parallel |
| Job ≤ 2 ms | 0.65 ms | — |
| `result_mbox` | Not yet: results on ports + lag-window stream | Milestone C |
