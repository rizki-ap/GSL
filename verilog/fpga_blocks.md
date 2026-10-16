# Gunshot Locator — FPGA Block Specification (v2, as built)

- **Target:** DE10-Standard (Cyclone V SX 5CSXFC6D6F31C6), fabric + HPS, or Nios V
- **Companion to:** `timing_budget.md` (budgets) and `SIMULATION.md` (how to run the checks)
- **Status:** Milestones A and B built and verified in simulation (Verilator / Icarus), and area-optimised (§5). Blocks marked *planned* are Milestone C or D. Nothing has been through Quartus yet.
- **Version:** 2026-10-11 (v2 as built and area-optimised; v1, 2026-10-08, was the pre-implementation proposal)

v1 → v2 changes are listed per block under **Changes vs v1**.

---

## 0. Conventions

| Item | Definition |
|---|---|
| Clock | Single domain `clk_sys` = 100 MHz. All blocks are synchronous to it. |
| Reset | `rst_n`, active-low, asynchronous assert. Omitted from the port tables. |
| Strobes | Every `*_valid` / `*_done` / `*_evt` is a 1-cycle pulse. |
| `s16` / `s18` | 16-bit / 18-bit two's complement. Stream samples are s18 with \|x\| < 2^15. |
| `x_ts` | 64-bit 100 kHz sample index, already corrected for `fe_cond` latency. Every SW-side time is in this base. |
| `y_idx` / `y_ts` | 20 kHz MB-stream index; `y_ts` is the same sample's time in the 100 kHz base (decimator group delay removed). |
| Frame | One 4-channel sample every `FRAME_DIV` clocks: 1000 at 100 MHz = 100 kHz. Simulation runs at 250 by default; results are identical. Minimum about 64. |
| Ring word | 72 bit = `{ch0, ch1, ch2, ch3}`, s18 each, ch0 in the top bits. |
| Pair order | (0,1) (0,2) (0,3) (1,2) (1,3) (2,3). Lag > 0: second channel lags. Matches `gs_det_tdoa.py`. |
| Configuration | Ports on `gsl_core` for now (43 values); `csr` will drive them in Milestone C. |

---

## 1. Block overview

```
replay_src -> timebase -> fe_cond -+-> sw_ring (cap_ring) ----------------+
 (adc_if, D)                       +-> sw_det -----------+                 |
                                   +-> mb_decim -+-> mb_ring (cap_ring) ---+
                                                 +-> mb_det -+             |
                                 sw_det / mb_det -> evt_sched -> gcc_engine (loader, fft2048,
                                                                 xspec_phat, peak_search)
                                                              -> result_mbox (C) -> processor
```

| # | Block | Folder | Status | Per-sample / per-job cost (measured) |
|---|---|---|---|---|
| 1 | `adc_if` | — | planned (D) | — |
| 2 | `replay_src` | `stream/` | built, verified (B) | 1 pop per frame |
| 3 | `timebase` | `stream/` | built (B); PPS planned (D) | — |
| 4 | `fe_cond` | `stream/` | built, verified (B) | 49 of 1000 cycles per frame; latency 4 samples |
| 5 | `sw_det` | `stream/` | built, verified (B) | 6 cycles per frame (+4 per noise block) |
| 6 | `mb_decim` | `stream/` | built, verified (B) | 263 cycles per 5 frames; delay 31 samples |
| 7 | `mb_det` | `stream/` | built, verified (B) | 9 cycles per 5 frames (+32 scan, +12 thresholds) |
| 8 | `cap_ring` ×2 | `gcc_engine/` | built, verified (A, B) | 1 write per frame / per 5 frames |
| 9 | `evt_sched` | `stream/` | built, verified (B) | every cycle |
| 10 | `gcc_loader` | `gcc_engine/` | built, verified (A, B) | 2,052 cycles per job |
| 11 | `fft2048` | `fft2048/` | built, verified (A) | 11,321 cycles per transform |
| 12 | `xspec_phat` | `xspec_phat/` | built, verified (A) | 6,167 cycles per job |
| 13 | `peak_search` | `gcc_engine/` | built, verified (A, B) | 406 (SW) / 86 (MB) cycles per job |
| — | `gcc_engine` | `gcc_engine/` | built, verified (A, B) | 65,239 (SW) / 64,919 (MB) cycles = 0.65 ms per job |
| — | `gsl_core` | `stream/` | built, verified (B) | Milestone B top |
| 14 | `result_mbox` | `gsl_mm/` | built, verified (C1) as `gsl_mbox` | ≤ 1,430 cycles per message |
| 15 | `csr` | `gsl_mm/` | built, verified (C1) in `gsl_mm` | — |
| 16 | `raw_dma` | — | planned (C, optional) | — |

---

## 2. Stream-rate blocks

### 2.1 `adc_if` — AD7606B controller *(planned, Milestone D)*

Unchanged from v1: one CONVST per frame from a free-running divider, readout ≤ 5 µs, all 4 channels from one conversion, `adc_err` on BUSY timeout. In Milestone B its role is played by `replay_src`, which already has the ADC pass-through port (`adc_valid`, `adc_data`). In Milestone D the frame cadence moves from `replay_src` to `adc_if`.

### 2.2 `replay_src` *(built)*

| Port | Dir | Width | Description |
|---|---|---|---|
| `run` | in | 1 | Enables the frame cadence |
| `src_sel` | in | 1 | 1 = replay FIFO, 0 = ADC pass-through |
| `rp_wr_en`, `rp_wr_data` | in | 1, 64 | Frame writes from the host, `{ch0..ch3}` s16 |
| `rp_level`, `rp_full` | out | 11, 1 | FIFO fill level |
| `adc_valid`, `adc_data` | in | 1, 64 | Live ADC frames (Milestone D) |
| `frame_tick` | out | 1 | Cadence strobe, every `FRAME_DIV` clocks while `run` |
| `src_valid`, `src_data` | out | 1, 64 | Frame to `fe_cond` |
| `rp_underflow` | out | 1 | Sticky: a tick found the FIFO empty (timeline broken) |

| Parameter | Value |
|---|---|
| FIFO | 1024 frames (10 ms at 100 kHz) |
| Cadence | `FRAME_DIV` clocks (parameter; 1000 = 100 kHz at 100 MHz) |

**Changes vs v1:** FIFO 1024 frames (v1 ≥ 200); the frame cadence is generated here until `adc_if` exists.

### 2.3 `timebase` *(built)*

64-bit frame index `ts64`; `src_ts` (valid with `src_valid`) carries the index of the frame being delivered. `ts_clear` resets it. PPS latching (`pps_ts`, `pps_frac`) is planned for Milestone D.

### 2.4 `fe_cond` *(built)*

Per channel: `a = x + offset` → `g = sat18(round(a·gain, 14))` → 9-tap fractional-delay FIR (Q1.16) → DC blocker → output.

| Port | Dir | Width | Description |
|---|---|---|---|
| `src_valid`, `src_data`, `src_ts` | in | 1, 64, 64 | Frame in |
| `cfg_gain` | in | 4 × 16 | Q2.14, 16384 = 1.0 |
| `cfg_offset` | in | 4 × s16 | |
| `cfg_dtrim` | in | 4 × s9 | Fractional delay, 1/128 sample, −128 … +128 |
| `cfg_dcblk_en` | in | 1 | DC blocker enable |
| `x_valid`, `x_data`, `x_ts` | out | 1, 72, 64 | `x_ts = src_ts − 4` |

| Parameter | Value |
|---|---|
| Latency | 4 samples, identical on all channels; `x_ts` already corrected |
| Delay trim | dtrim = 0 is an exact unit impulse (transparent); error ≤ 0.025 sample (0.25 µs) below 20 kHz |
| DC blocker | Shift-only, `acc += ((y − y_prev) << 13) − (acc >>> 13)`; corner 1.94 Hz |
| Rounding | Convergent (round-half-to-even) |
| Timing | Sequential over channels: 4 × 12 + 1 = 49 cycles per frame |
| ROM | `fe_fdelay.hex`: 129 phases × 9 taps; negative trims use mirrored taps |

**Changes vs v1:** 9 taps (v1 8); the coefficient ROM has 2-cycle latency (a 1-sample bug from assuming 1 cycle was caught by the bit-exact check and fixed).

### 2.5 `sw_det` *(built; substantially extended)*

| Port | Dir | Width | Description |
|---|---|---|---|
| `x_valid`, `x_data`, `x_ts` | in | 1, 72, 64 | Conditioned stream |
| `cfg_thr_mult2` | in | 16 | Candidate threshold (mult²): 100 = 10× noise RMS |
| `cfg_noise_min` | in | 32 | Minimum noise mean square, LSB² (default 16) |
| `cfg_holdoff`, `cfg_dt_min` | in | 16, 16 | Samples: 2000 (20 ms), 1000 (10 ms MB blanking) |
| `cfg_feat_win` | in | 9 | Feature window, samples (256), > 100 |
| `cfg_sym` | in | 8 | N-wave symmetry, 1/64 units (32 = 0.5) |
| `cfg_tn_max` | in | 8 | Max N-wave duration, samples (60 = 0.6 ms) |
| `cfg_pk_mult` | in | 8 | Peak margin over threshold, power ratio (9 = 3×) |
| `rec_valid`, `rec_ts`, `rec_ref` | out | 1, 64, 2 | SW record (end of feature window) |
| `rec_feat` | out | 4 × 54 | Per channel `{pk_pos, pos_off, pk_neg, neg_off}` |
| `conf_valid`, `conf_ts` | out | 1, 64 | SW confirmed (1 ms after onset) |
| `sw_qual` | out | 1 | Level: a candidate is being qualified |
| `sw_blank` | out | 1 | Level: `cfg_dt_min` samples after a confirmed SW |
| `noise_valid` | out | 1 | First clean noise block seen |

Behaviour:
- **Noise floor:** mean square per channel over 512-sample blocks, accepted only if the block had no candidate, qualification or holdoff.
- **Candidate:** a rising crossing of x² > `thr_mult2 × max(noise, noise_min)` on any channel. `rec_ts` is the crossing sample; `rec_ref` is the lowest channel that crossed.
- **Qualification (1 ms, ref channel):** with P = max(max x, 0) and Q = max(−min x, 0), it is a shockwave only if all three hold:
  - **Symmetry:** 64·min(P, Q) ≥ `cfg_sym`·max(P, Q)
  - **Duration:** |offset(min) − offset(max)| ≤ `cfg_tn_max`
  - **Margin:** max(P, Q)² ≥ `cfg_pk_mult`·threshold
- **Rejected candidate** (typically the muzzle blast itself): no record, no holdoff, no blanking.

**Changes vs v1:**
- **Qualification added.** Without it the detector fires on muzzle blasts. A front-steepness test was tried first and failed below 50 m.
- **New ports:** `conf_valid` / `conf_ts` and `sw_qual`.
- **`noise_min`:** needed because the generator's ADC noise is below 1 LSB.
- **Blanking** now starts at confirmation.

### 2.6 `mb_decim` *(built)*

| Port | Dir | Width | Description |
|---|---|---|---|
| `x_valid`, `x_data`, `x_ts` | in | 1, 72, 64 | 100 kHz stream |
| `y_valid`, `y_data` | out | 1, 72 | 20 kHz stream |
| `y_ts` | out | 64 | Filter-centre time: newest `x_ts` − 31 |
| `y_idx` | out | 32 | Output index (mb_ring address) |

| Parameter | Value |
|---|---|
| Filter | 63-tap linear phase (Kaiser β 6, fc 7 kHz), Q1.17, DC gain exactly 1 (`decim63.hex`) |
| Response | 0.02 dB ripple to 4 kHz; < −71 dB from 16 kHz (spec 0.1 dB / −60 dB) |
| Output | Every 5th input, once 63 inputs have been seen |
| Timing | One multiply-accumulator, channels in sequence: 263 cycles per output; 128-slot history |

### 2.7 `mb_det` *(built; extended)*

| Port | Dir | Width | Description |
|---|---|---|---|
| `y_valid`, `y_data`, `y_ts`, `y_idx` | in | 1, 72, 64, 32 | Decimated stream |
| `sw_blank`, `sw_qual`, `sw_conf` | in | 1 | From `sw_det` |
| `cfg_thr_mult`, `cfg_settle_mult` | in | 8, 8 | Energy ratios: 20, 9 (= 3× RMS) |
| `cfg_dwell`, `cfg_holdoff`, `cfg_settle_timeout` | in | 16 each | 20 kHz samples: 100 (5 ms), 400 (20 ms), 4000 (200 ms) |
| `cfg_noise_min` | in | 32 | LSB² |
| `rec_valid`, `rec_ts`, `rec_idx`, `rec_ref` | out | 1, 64, 32, 2 | MB record |
| `rec_reliable` | out | 1 | Armed by settling (0 = armed by timeout) |

Behaviour:
- **Energy:** 16-sample sliding energy per channel. Trigger when E > `thr_mult`·16·noise on any channel (rising edge) while armed.
- **Arming is latched:** once all channels have been quiet (E < `settle_mult`·16·noise) for `dwell` samples, the detector stays armed until it triggers or a new SW blanking starts.
- **Held events:** a trigger while `sw_qual` is high is held. If `sw_conf` follows, it is dropped together with its holdoff (it was the shockwave's own energy). Otherwise it is released.
- **Onset:** the oldest sample in the 16-sample window whose square exceeds `thr_mult`·noise; `rec_ts = y_ts − 5·back`.
- **Continuous operation:** SW↔MB pairing and the `dt_max` timeout are done by the processor.

**Changes vs v1:**
- **Latched arming:** without it, 1 of 8 burst blasts was missed.
- **Held events:** needed because SW qualification takes 1 ms.
- **`settle_timeout` / `rec_reliable`** added.

### 2.8 `cap_ring` (`sw_ring`, `mb_ring`) *(built)*

4096 × 72-bit dual-port RAM, 2-cycle read latency. `sw_ring` is written at `x_ts[11:0]` (41 ms of history), `mb_ring` at `y_idx[11:0]` (205 ms). Writes never stall; `evt_sched` enforces the read deadline.

---

## 3. Event-rate blocks

### 3.1 `evt_sched` *(built; extended)*

| Port | Dir | Width | Description |
|---|---|---|---|
| `x_valid`, `x_ts`, `y_valid`, `y_idx` | in | — | Last written ring positions |
| `sw_rec_*`, `sw_conf_*`, `mb_rec_*` | in | — | Detector records |
| `cfg_{sw,mb}_{pre,len,maxlag}` | in | 12, 12, 7 | Windows: SW 300 / 1800 / 100, MB 120 / 1720 / 20 |
| `cfg_mb_mask_pre/post` | in | 8, 8 | SW masked in MB windows, 20 kHz samples: 20 / 60 (−1 … +3 ms) |
| `cfg_sw_mask_pre/post` | in | 16, 16 | Later MB masked in SW windows, 100 kHz samples: 100 / 2000 |
| `job_valid`, `job_ready` | out / in | 1 | Engine handshake |
| `job_type`, `job_id`, `job_start`, `job_len`, `job_maxlag` | out | — | Job |
| `job_ts0`, `job_mask_en`, `job_mask_lo`, `job_mask_hi` | out | 32, 4, 4 × 32, 4 × 32 | Window time base and mask intervals |
| `sw_id`, `mb_id` | out | 16 | Event id for the record valid this cycle |
| `late_valid`, `late_id`, `late_cnt`, `ovf_cnt` | out | — | Health |
| `busy` | out | 1 | Anything pending |

Behaviour:
- **Event ids:** one counter shared by SW and MB, assigned on record arrival (SW first if simultaneous).
- **Pending FIFOs:** SW and MB, 4 entries each. Release happens once the last window sample is written to the ring.
- **Job queue:** 4 entries (v1: 8; occupancy never exceeds 1).
- **Dispatch deadline:** a job is dropped and counted as late if its window start could be overwritten during the load: `last − start ≥ 4096 − 16`.
- **SW masks in MB windows:** at MB release, every remembered SW onset (last 4) whose −1 … +3 ms interval overlaps the window.
- **MB masks in SW windows:** at SW release, every remembered MB onset (last 4) that lies after the SW onset and overlaps the window.

**Changes vs v1:**
- **Masks** replace v1's "trim the MB window to 40 ms" idea, which can't work at 200 m: the next shockwave lands 9.5 ms after the blast.
- **Event ids** are assigned here.
- **Mask outputs** added.

### 3.2 `gcc_loader` *(built)*

Ring → Z0 = ch0 + j·ch1, Z1 = ch2 + j·ch3 at bit-reversed addresses, zero-padded to 2048. All 4 channels come from the same ring address. Mask inputs `ld_ts0`, `ld_step` (1 SW, 5 MB) and 4 intervals: window sample n is zeroed on all channels if `ts0 + n·step` lies inside an enabled interval. **2,052 cycles.**

### 3.3 `fft2048` *(built)*

Radix-2 DIT, in place, s18 data and twiddles, block-floating-point. Input bit-reversed, output natural. Input contract \|x\| < 0.5 FS. IFFT via re/im swap. **11,321 cycles per transform.** Full interface: `fft2048/readme.md`.

### 3.4 `xspec_phat` *(built)*

Unpack, 6 cross-spectra G_ab = X_b·conj(X_a), PHAT normalisation to magnitude 32767 (rsqrt ROM, phase exact), band mask, `cfg_floor`, packed IFFT inputs. **6,167 cycles.** Full interface: `xspec_phat/readme.md`.

### 3.5 `peak_search` *(built)*

Two passes over ±max_lag: per-pair argmax (first maximum, as numpy) and sidelobe (max \|r\| more than 2 samples from the peak); lag-window stream for the processor's sub-sample refinement. PSR = pk / side is computed by the processor. **406 (SW) / 86 (MB) cycles.**

### 3.6 `result_mbox` *(built in Milestone C1 as `gsl_mbox`; see `gsl_mm/c1_interface.md`)*

As in v1: descriptor ring to the processor with IRQ. Contents now known from `gsl_core`'s outputs:
- **SW / MB records:** id, timestamps, reference channel, features, reliability.
- **Job results:** lags, peaks, sidelobes, exponents, error flag.
- **Lag windows:** 6 × (2·max_lag + 1).
- **Health:** late / overflow / underflow counters.

---

## 4. Infrastructure *(`csr` built in Milestone C1 inside `gsl_mm`; register map in `gsl_mm/c1_interface.md`)*

### 4.1 `csr`

Drives the 43 configuration ports of `gsl_core` and reads its status. Register map to be fixed in Milestone C. Configuration groups:
- `fe` (gain, offset, dtrim, DC enable)
- `sw` (thresholds, holdoff, blanking, features, N-wave tests)
- `mb` (thresholds, settle, holdoff, timeout)
- windows
- PHAT bands and floors
- masks

Status:
- `rp_underflow`, `late_cnt`, `ovf_cnt`, `res_err`
- `noise_valid` (SW, MB)
- event counters

### 4.2 `raw_dma` (optional)

Unchanged from v1.

---

## 5. Resources (yosys `synth_intel_alm -family cyclonev`, indicative)

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

---

## 6. Verification (as done)

| Level | Where | Method | Result |
|---|---|---|---|
| `fft2048` | `fft2048/` | DC, tones, complex, IFFT round trip vs double DFT; bit-true model | bit-exact; SNR 82–101 dB |
| `xspec_phat` | `xspec_phat/` | 4 signal cases through FFT → xspec → IFFT | bit-exact; fixed-point TDOA error ≤ 0.2 µs |
| **Milestone A** | `gcc_engine/` | 10 jobs incl. generator scenarios, ring wrap-around | bit-exact; bearing ≤ 0.25° |
| **Milestone B** | `stream/` | 15 scenarios played continuously: 50–1500 m, 3 noise levels, bullet drag (incl. subsonic), 4 mismatch scenarios, 900 rpm burst | every block bit-exact; Δt within budget; bearing ≤ 0.75° (2.0° el at 1500 m / 0.5 Pa) |
| Cadence | `stream/` | 250 vs 1000 clocks per sample on 6 scenarios | byte-identical logs |
| Simulators | `stream/` | Icarus vs Verilator (50 m scenario) | byte-identical logs |

**Separation rule:** detector-side computations never use generator parameters. See `SIMULATION.md`.

---

## 7. Folder map

| Folder | Contents |
|---|---|
| `fft2048/` | `fft2048`, `fft2048_bfly`, twiddle ROM, `fft_buf_bank`, `fft_sdp_ram` |
| `xspec_phat/` | `xspec_phat`, `xspec_phat_unit`, rsqrt ROM |
| `gcc_engine/` | `gcc_engine`, `gcc_loader`, `peak_search`, `cap_ring`; Milestone A testbench |
| `stream/` | `replay_src`, `timebase`, `fe_cond`, `sw_det`, `mb_decim`, `mb_det`, `evt_sched`, `gsl_core`; Milestone B testbench |
| `old/` | Earlier `fft1024` / `gcc_phat_top` RTL, superseded |
