# Gunshot Locator — FPGA Block Specification

- **Target:** DE10-Standard (Cyclone V SX 5CSXFC6D6F31C6), fabric + HPS
- **Companion to:** `timing_budget.md` (all budgets here are derived there)
- **Status:** proposal — port names, register map and descriptor layout are a starting point, not frozen
- **Version:** 2026-10-08

---

## 0. Conventions

| Item | Definition |
|---|---|
| Clock | Single domain `clk_sys` = 100 MHz from PLL (50 MHz board oscillator). The LW H2F bridge's FPGA side is also clocked by `clk_sys`. |
| Reset | `rst_n`, active-low, synchronous to `clk_sys`. Omitted from the port tables. |
| Strobes | Every `*_valid` / `*_evt` / `*_done` is a 1-cycle pulse. |
| `s16` / `s18` | 16-bit / 18-bit two's complement. `s18` is Q1.17. |
| `ts64` | 64-bit sample index at 100 kHz. Every timestamp in the system is a `ts64`, never a clock-cycle count. |
| Frame | One 4-channel sample set, every 10 µs = 1000 `clk_sys` cycles. |
| MB frame | One 4-channel decimated sample set, every 50 µs (20 kHz). |
| Pair order | (0,1) (0,2) (0,3) (1,2) (1,3) (2,3) — must match `gs_det_tdoa.py`. |
| Lag sign | lag > 0 ⇒ second channel of the pair lags the first — must match `gs_det_tdoa.py`. |

---

## 1. Block overview

```
Stream (every 10 µs, hard real time):
 adc_if ─┐
 replay_src ─┴─► fe_cond ─┬─► sw_det ─────┐
                          ├─► sw_ring     │
                          └─► mb_decim ─┬─► mb_det ─┐
                                        └─► mb_ring │
Event (per window, ≤ 2 ms):                         │
 evt_sched ◄────────────────────────────────────────┘
   └─► gcc_loader ─► fft2048 ⇄ xspec_phat
                     └─► peak_search ─► result_mbox ─► HPS (LW bridge + IRQ)
Shared: timebase, csr, raw_dma
```

| # | Block | Class | Throughput | Latency / deadline | Key spec |
|---|---|---|---|---|---|
| 1 | `adc_if` | stream | 4 ch / 10 µs | ≤ 5 µs CONVST → data | one CONVST for all channels |
| 2 | `replay_src` | stream | 4 ch / 10 µs | same cadence as ADC | bit-exact replay of Python WAVs |
| 3 | `timebase` | stream | 1 / frame | 0 | 64-bit, never wraps |
| 4 | `fe_cond` | stream | 4 ch / 10 µs | constant ≤ 0.5 ms | delay trim 1/128 sample |
| 5 | `sw_det` | stream | 4 ch / 10 µs | trigger ≤ 0.5 ms | onset ±1 sample |
| 6 | `mb_decim` | stream | 4 ch / 10 µs in | constant ≤ 0.31 ms | ÷ 5, linear phase |
| 7 | `mb_det` | stream | 4 ch / 50 µs | trigger ≤ 1 ms | onset ≤ 0.4 ms error |
| 8 | `sw_ring` / `mb_ring` | stream | 1 write / frame | read deadline 23 / 119 ms | 4096 deep |
| 9 | `evt_sched` | event | ≥ 1 job / 33 ms | dispatch ≤ 10 cycles | queue ≥ 4 |
| 10 | `gcc_loader` | event | — | ≤ 4,096 cycles | shared-reference window |
| 11 | `fft2048` | event | — | ≤ 15,360 cycles / transform | BFP, ≤ 1 µs TDOA error |
| 12 | `xspec_phat` | event | — | ≤ 6,200 cycles | rsqrt error ≤ 2⁻¹² |
| 13 | `peak_search` | event | — | ≤ 1,206 cycles | ± max_lag only |
| 14 | `result_mbox` | event | — | IRQ ≤ 1 µs after write | 8 slots, never overwrite unread |
| 15 | `csr` | infra | — | applied at frame boundary | atomic 64-bit reads |
| 16 | `raw_dma` | infra | 800 kB/s | FIFO ≥ 1.28 ms | optional |

Per-job engine total (blocks 10–13): ≈ 88k cycles = 0.88 ms, budget 2 ms.

---

## 2. Stream-rate blocks

These run every frame and may never stall. A missed frame is a hard error.

### 2.1 `adc_if` — AD7606B controller

| Port | Dir | Width | Description |
|---|---|---|---|
| `adc_convst` | out (pin) | 1 | Conversion start, rising edge once per frame |
| `adc_cs_n` | out (pin) | 1 | Chip select |
| `adc_sclk` | out (pin) | 1 | Serial clock, 25 MHz (`clk_sys` / 4) |
| `adc_reset` | out (pin) | 1 | Power-up reset pulse |
| `adc_busy` | in (pin) | 1 | Conversion in progress (2-FF synchronised) |
| `adc_dout` | in (pin) | 1–2 | Serial data line(s) |
| `frame_tick` | out | 1 | Strobe at each CONVST; master cadence for the whole design |
| `adc_ch0..3` | out | s16 | One sample per channel |
| `adc_valid` | out | 1 | Strobe when all 4 channels are read |
| `adc_err` | out | 1 | BUSY timeout or framing error (sticky in `csr`) |

| Parameter | Value |
|---|---|
| Frame period | exactly 1000 cycles = 10.000 µs, from a free-running counter, never stalled |
| Inter-channel skew | single CONVST for all 4 channels; skew = ADC aperture matching (ns) |
| BUSY timeout | `adc_err` if BUSY is not low within datasheet t_CONV(max) |
| Readout | e.g. 64 SCLK at 25 MHz ≈ 2.6 µs on one DOUT line |
| CONVST → `adc_valid` | ≤ 5 µs |

Check the AD7606B serial-mode section for DOUT line count and which channels appear on which line before fixing the pin map.

### 2.2 `replay_src` — test-vector source (strongly recommended)

| Port | Dir | Width | Description |
|---|---|---|---|
| `rp_wr_data` | in | 64 | One frame (4 × s16) from the HPS (Avalon-MM write or F2SDRAM read master) |
| `rp_wr_en` | in | 1 | FIFO write |
| `rp_fifo_lvl` | out | 12 | FIFO fill level, to `csr` |
| `frame_tick` | in | 1 | Cadence from `adc_if` |
| `src_sel` | in | 1 | 0 = ADC, 1 = replay (from `csr`) |
| `adc_ch0..3`, `adc_valid` | in | s16, 1 | Live ADC stream |
| `src_ch0..3` | out | s16 | Muxed stream into `fe_cond` |
| `src_valid` | out | 1 | Muxed valid |
| `rp_underflow` | out | 1 | Sticky, to `csr` |

| Parameter | Value |
|---|---|
| Cadence | Identical to the ADC (driven by `frame_tick`), so downstream blocks cannot tell the difference |
| FIFO depth | ≥ 200 frames (2 ms) |
| Purpose | Play `gs_gen_apply_adc.py` output through the FPGA and compare against `gs_det_*.py` bit-for-bit |

### 2.3 `timebase`

| Port | Dir | Width | Description |
|---|---|---|---|
| `src_valid` | in | 1 | Frame strobe |
| `pps_in` | in (pin) | 1 | Optional GNSS PPS (2-FF synchronised) |
| `ts64` | out | 64 | Sample index, increments on `src_valid` |
| `pps_ts` | out | 64 | `ts64` latched at each PPS rising edge |
| `pps_frac` | out | 10 | Cycles since last frame at the PPS edge (10 ns resolution) |
| `pps_valid` | out | 1 | Strobe |

| Parameter | Value |
|---|---|
| Increment | Same cycle as `src_valid` |
| Wrap | Never in practice (2⁶⁴ samples ≈ 5.8 × 10⁶ years) |
| Clear | Only on reset or `CTRL.ts_clear` |

### 2.4 `fe_cond` — DC block, gain, offset, fractional-delay trim

| Port | Dir | Width | Description |
|---|---|---|---|
| `src_ch0..3`, `src_valid` | in | s16, 1 | Raw frame |
| `ts64` | in | 64 | Current sample index |
| `cfg_gain[0..3]` | in | 16 | Q2.14 per-channel gain |
| `cfg_offset[0..3]` | in | s16 | Per-channel offset |
| `cfg_dtrim[0..3]` | in | s9 | Fractional delay in 1/128-sample units |
| `cfg_dcblk_en` | in | 1 | DC blocker enable |
| `x_ch0..3` | out | s18 | Conditioned samples |
| `x_valid` | out | 1 | Strobe |
| `x_ts` | out | 64 | `ts64` of the input sample this output corresponds to (latency-corrected) |
| `lat_fe` | out | 16 | Constant latency in samples, to `csr` |

| Parameter | Value |
|---|---|
| Latency | Constant, identical on all 4 channels, ≤ 50 samples (0.5 ms) |
| DC blocker | 1st-order IIR, corner ≤ 2 Hz, identical coefficients on all channels |
| Delay trim | ± 1 sample used, step 1/128 sample = 0.078 µs (8-tap × 128-phase polyphase, or cubic Farrow) |
| Order | offset → gain → delay trim → DC block; saturate to s18 |
| Compute | ≤ 250 MAC / channel / frame (one time-shared DSP) |

Digital filters are identical across channels, so they cost no TDOA budget; only `dtrim` changes relative timing. `dtrim` is the calibration knob for AFE mismatch and, if the onboard multiplexed ADC is kept, its channel skew.

### 2.5 `sw_det` — shockwave onset detector

| Port | Dir | Width | Description |
|---|---|---|---|
| `x_ch0..3`, `x_valid`, `x_ts` | in | s18, 1, 64 | Conditioned stream |
| `cfg_thr_mult` | in | 8 | Threshold multiplier (default 10 = `sw_threshold_mult`) |
| `cfg_noise_win` | in | 12 | Noise-floor window, samples (default 500 = 5 ms) |
| `cfg_sta_win` | in | 6 | Short-term window, samples (default 10 = 0.1 ms) |
| `cfg_holdoff` | in | 16 | Re-trigger holdoff, samples (default 2000 = 20 ms) |
| `cfg_tn_max` | in | 8 | Feature window, samples (default 150 = 1.5 ms) |
| `cfg_dt_min` | in | 16 | MB blanking after SW, samples (default 1000 = 10 ms) |
| `sw_evt` | out | 1 | Onset strobe |
| `sw_ts` | out | 64 | `x_ts` of first threshold crossing |
| `sw_ref_ch` | out | 2 | First channel to cross |
| `sw_pk_pos[0..3]` | out | s18 | Positive peak per channel |
| `sw_pk_neg[0..3]` | out | s18 | Negative peak per channel |
| `sw_tn[0..3]` | out | 8 | N-wave duration (positive-to-negative peak spacing), samples |
| `sw_feat_valid` | out | 1 | Strobe when features are complete |
| `sw_blank` | out | 1 | Level, high for `cfg_dt_min` after each `sw_evt` (to `mb_det`) |

| Parameter | Value |
|---|---|
| Statistic | Short-term amplitude over `sta_win` vs noise floor over `noise_win`, any channel |
| Onset timestamp | First crossing on `sw_ref_ch`, 1-sample (10 µs) resolution |
| Trigger latency | ≤ 50 samples (0.5 ms) after onset |
| Feature latency | `sw_feat_valid` ≤ 2 ms after onset |
| Holdoff | ≤ 20 ms (must stay below the 66.7 ms burst interval) |
| Noise floor | Frozen from onset to end of holdoff |
| Event rate | ≥ 1 per 66.7 ms sustained |

### 2.6 `mb_decim` — ÷ 5 decimator for the muzzle-blast path

| Port | Dir | Width | Description |
|---|---|---|---|
| `x_ch0..3`, `x_valid`, `x_ts` | in | s18, 1, 64 | 100 kHz stream |
| `y_ch0..3` | out | s18 | 20 kHz stream |
| `y_valid` | out | 1 | Strobe, once per 5 frames |
| `y_ts` | out | 64 | Timestamp in 100 kHz units, corrected for group delay |
| `y_idx` | out | 32 | 20 kHz sample index (ring addressing) |
| `lat_dec` | out | 8 | Group delay in 100 kHz samples, to `csr` |

| Parameter | Value |
|---|---|
| Filter | Symmetric (linear-phase) FIR, ≤ 63 taps, polyphase |
| Passband | 0–4 kHz, ripple ≤ 0.1 dB |
| Stopband | ≥ 16 kHz, ≥ 60 dB |
| Group delay | (N_taps − 1)/2 ≤ 31 samples = 0.31 ms, constant, removed in `y_ts` |
| Compute | ≈ 13 MAC / input sample / channel |

Correcting `y_ts` here means MB and SW timestamps compare directly when forming Δt (budget ≤ 0.5 ms).

### 2.7 `mb_det` — muzzle-blast onset detector

| Port | Dir | Width | Description |
|---|---|---|---|
| `y_ch0..3`, `y_valid`, `y_ts`, `y_idx` | in | s18, 1, 64, 32 | Decimated stream |
| `sw_blank` | in | 1 | Blanking from `sw_det` |
| `cfg_thr_mult` | in | 8 | Default 20 (`mb_energy_threshold_mult`) |
| `cfg_noise_win` | in | 12 | Noise-floor window, 20 kHz samples |
| `cfg_sta_win` | in | 6 | Short-term window (default 20 = 1 ms) |
| `cfg_settle_mult` | in | 8 | Settling threshold relative to noise floor |
| `cfg_holdoff` | in | 16 | Re-trigger holdoff |
| `mb_evt` | out | 1 | Onset strobe |
| `mb_ts` | out | 64 | Onset, 100 kHz units |
| `mb_idx` | out | 32 | Onset, 20 kHz index (for `mb_ring`) |
| `mb_ref_ch` | out | 2 | First channel to cross |
| `mb_reliable` | out | 1 | Low if MB arrived before post-SW energy settled |

| Parameter | Value |
|---|---|
| Statistic | Short-term energy vs adaptive noise floor |
| Arming | After each SW: wait for `sw_blank` low **and** energy to settle within `settle_mult` × floor (`project.md`, limitation 5) |
| Onset accuracy | ≤ 0.4 ms (≤ 8 samples at 20 kHz); back-track from the crossing to the start of the energy rise |
| Trigger latency | ≤ 1 ms (20 samples) after onset |
| Mode | Continuous; SW↔MB pairing and the 5 s timeout live on the HPS |

### 2.8 `sw_ring` / `mb_ring` — capture rings

| Port | Dir | Width | Description |
|---|---|---|---|
| `wr_data` | in | 72 | 4 × s18 |
| `wr_en` | in | 1 | `x_valid` (SW) / `y_valid` (MB) |
| `wr_addr` | in | 12 | `x_ts[11:0]` (SW) / `y_idx[11:0]` (MB) |
| `rd_addr` | in | 12 | From `gcc_loader` |
| `rd_data` | out | 72 | 4 × s18, 2-cycle read latency |

| Parameter | `sw_ring` | `mb_ring` |
|---|---|---|
| Rate | 100 kHz | 20 kHz |
| Depth | 4096 (41 ms) | 4096 (205 ms) |
| Window | 300 pre + 1500 post | 120 pre + 1600 post |
| Read deadline after window closes | 2296 samples = 23 ms | 2376 samples = 119 ms |
| M10K | ≈ 32 | ≈ 32 |

Writes never stall. Overlapping events are just different start addresses into the same ring.

---

## 3. Event-rate blocks

One shared engine processes one job (one SW or one MB window) at a time. Budget per job: 2 ms.

### 3.1 `evt_sched` — job scheduler

| Port | Dir | Width | Description |
|---|---|---|---|
| `sw_evt`, `sw_ts`, `sw_ref_ch` | in | 1, 64, 2 | SW onset |
| `sw_pk_*`, `sw_tn`, `sw_feat_valid` | in | — | SW features (latched into job metadata) |
| `mb_evt`, `mb_ts`, `mb_idx`, `mb_ref_ch`, `mb_reliable` | in | — | MB onset |
| `x_ts`, `y_idx` | in | 64, 32 | Current write positions |
| `eng_busy` | in | 1 | Engine occupied |
| `job_valid` | out | 1 | Dispatch strobe |
| `job_type` | out | 1 | 0 = SW, 1 = MB |
| `job_id` | out | 16 | Incrementing ID |
| `job_start` | out | 12 | Ring address of window start (onset − pre) |
| `job_len` | out | 11 | 1800 (SW) / 1720 (MB) |
| `job_maxlag` | out | 7 | 100 (SW) / 20 (MB) |
| `job_meta` | out | — | ts, ref_ch, features, reliable flag → `result_mbox` |

| Parameter | Value |
|---|---|
| SW release | when `x_ts` ≥ `sw_ts` + 1500 |
| MB release | when `y_idx` ≥ `mb_idx` + 1600 |
| Queue | ≥ 4 jobs, FIFO |
| Dispatch | ≤ 10 cycles after release when engine idle |
| Deadline | Start ≤ 23 ms (SW) / 119 ms (MB) after release, else drop and set `job_late` |
| Overflow | Sticky counter in `csr` |

### 3.2 `gcc_loader` — window loader

| Port | Dir | Width | Description |
|---|---|---|---|
| `job_*` | in | — | From `evt_sched` |
| `sw_rd_addr` / `mb_rd_addr` | out | 12 | Ring read address |
| `sw_rd_data` / `mb_rd_data` | in | 72 | Ring read data |
| `z0_wr_*`, `z1_wr_*` | out | 11 addr, 36 data | Write ports of FFT buffers Z0, Z1 |
| `ld_done` | out | 1 | Strobe |

| Parameter | Value |
|---|---|
| Windowing | **Same start address for all 4 channels** (shared reference, `project.md` limitation 3) |
| Packing | Z0 = ch0 + j·ch1, Z1 = ch2 + j·ch3 |
| Length | `job_len` samples, zero-filled to 2048 |
| Taper | Optional (CSR), default matches `gs_det_tdoa.py` |
| Cycles | ≤ 4,096 (41 µs) |
| Optimisation | Writing in bit-reversed address order removes the FFT's 2N bit-reverse pass (≈ 27 % fewer cycles) |

### 3.3 `fft2048` — FFT / IFFT engine

| Port | Dir | Width | Description |
|---|---|---|---|
| `fft_start` | in | 1 | Strobe |
| `fft_inverse` | in | 1 | 0 = FFT, 1 = IFFT |
| `fft_buf_sel` | in | 3 | Buffer: Z0, Z1, W0, W1, W2 |
| buffer ports | in/out | 11 addr, 36 data | 2 reads + 2 writes per cycle (true dual-port) |
| `fft_done` | out | 1 | Strobe |
| `fft_bfp_exp` | out | 5 | Total right shifts applied |

| Parameter | Value |
|---|---|
| Size | 2048, radix-2 DIT, in place |
| Data / twiddles | s18 complex / 18-bit |
| Scaling | Block-floating-point (shift a stage only when overflow is possible), total in `fft_bfp_exp` |
| Cycles | ≤ 15,360 per transform (2N + N/2·log₂N); ≤ 11,264 with bit-reversed input |
| Per job | 2 forward + 3 inverse ≤ 76,800 cycles |
| IFFT scaling | No 1/N: only argmax and relative values are used downstream |
| Accuracy | ≤ 1 µs TDOA error (95 %) vs double precision on the replay test set |

Buffers: 5 × 2048 × 36 bit (Z0, Z1, W0, W1, W2) ≈ 40 M10K.

### 3.4 `xspec_phat` — unpack, cross-spectrum, PHAT, IFFT packing

| Port | Dir | Width | Description |
|---|---|---|---|
| `xp_start` | in | 1 | Strobe |
| `z0_rd_*`, `z1_rd_*` | in/out | 11 addr, 36 data | Two reads per cycle: bins k and N−k |
| `bfp_exp0/1` | in | 5 | Exponents of Z0, Z1 |
| `cfg_band_lo/hi` | in | 11 + 11 | Per job type (SW / MB) |
| `cfg_eps` | in | 18 | PHAT floor |
| `w0..w2_wr_*` | out | 11 addr, 36 data | Two writes per cycle: bins k and N−k |
| `xp_done` | out | 1 | Strobe |

| Step | Definition |
|---|---|
| Unpack | X_a[k] = (Z[k] + Z*[N−k]) / 2, X_b[k] = (Z[k] − Z*[N−k]) / 2j |
| Cross-spectrum | G = X_a · X_b*, exponents aligned |
| PHAT | G / max(\|G\|, eps); 1/√\|G\|² from ROM seed + 1 Newton step, relative error ≤ 2⁻¹² |
| Band mask | Bins outside [lo, hi] set to 0; default all bins (comparable with Python) |
| Packing | W_p[k] = G_2p[k] + j·G_2p+1[k]; W_p[N−k] = G_2p[k]* + j·G_2p+1[k]* |

| Parameter | Value |
|---|---|
| Cycles | ≤ 6,200 (1 bin / pair / cycle over k = 0…1024) |
| Pipeline latency | < 32 cycles |
| DSP | ≈ 6 (complex multiply + Newton step) |

### 3.5 `peak_search` — coarse TDOA peak

| Port | Dir | Width | Description |
|---|---|---|---|
| `ps_start` | in | 1 | Strobe |
| `w0..w2_rd_*` | in/out | 11 addr, 36 data | IFFT outputs: re = r_2p, im = r_2p+1 |
| `cfg_maxlag` | in | 7 | 100 (SW) / 20 (MB) |
| `pk_lag[0..5]` | out | s8 | Coarse lag, samples |
| `pk_val[0..5]` | out | s18 | Peak value |
| `pk_psr[0..5]` | out | 16 | Peak-to-sidelobe ratio, Q8.8 |
| `lagwin_wr_*` | out | 12 addr, 32 data | All lags in ± max_lag, to `result_mbox` |
| `ps_done` | out | 1 | Strobe |

| Parameter | Value |
|---|---|
| Lag range | ± max_lag only: addresses N−L … N−1, then 0 … L |
| Coarse peak | Argmax per pair, 1-sample resolution |
| PSR | Peak / largest \|r\| outside ± 2 samples of the peak (within ± max_lag) |
| Lag window | 2·max_lag + 1 values per pair (201 SW / 41 MB) for HPS sinc refinement |
| Cycles | ≤ 1,206 |

### 3.6 `result_mbox` — event descriptors to the HPS

| Port | Dir | Width | Description |
|---|---|---|---|
| `job_meta`, `bfp_exp`, `pk_*`, `lagwin_wr_*` | in | — | Per-job results |
| `avs_address` | in | 14 | Avalon-MM slave on LW H2F bridge (word address, 64 KB) |
| `avs_read`, `avs_readdata` | in, out | 1, 32 | Read port |
| `avs_waitrequest` | out | 1 | |
| `irq` | out | 1 | To `f2h_irq0` |
| `mbox_wr_ptr` | out | 3 | To `csr` |
| `mbox_rd_ptr` | in | 3 | From `csr` (written by HPS) |

| Parameter | Value |
|---|---|
| Slots | 8 × 8 KB ring (≈ 5 kB used per SW descriptor) |
| IRQ | Level, high while `wr_ptr` ≠ `rd_ptr`; rises ≤ 1 µs (100 cycles) after a slot completes |
| Overflow | If all slots are full, drop the new job and count it; never overwrite an unread slot |
| HPS read time | ≤ 1 ms per SW descriptor |

Proposed descriptor layout (32-bit words, offsets in bytes):

| Offset | Field |
|---|---|
| 0x000 | magic[15:0], type[16], job_id[31:17] |
| 0x004 / 0x008 | ts lo / hi (onset, 100 kHz units) |
| 0x00C | ref_ch[1:0], reliable[2], late[3], saturated[4], max_lag[14:8], bfp_exp[20:16] |
| 0x010 – 0x02C | SW features: pk_pos, pk_neg, tn per channel (packed) |
| 0x040 – 0x054 | pk_lag[6] (s8) + pk_psr[6] (Q8.8), packed |
| 0x058 – 0x06C | pk_val[6] |
| 0x080 – … | lag window: 6 pairs × (2·max_lag + 1) words, pair-major |

---

## 4. Infrastructure

### 4.1 `csr` — control / status registers

Avalon-MM slave on the LW H2F bridge, 32-bit.

| Parameter | Value |
|---|---|
| Config writes | Shadowed; applied together at the next frame boundary |
| 64-bit reads | Reading the low word latches the high word |
| Status bits | Sticky, write-1-to-clear |

Proposed register map:

| Offset | Name | Access | Content |
|---|---|---|---|
| 0x000 | `VERSION` | RO | Build ID |
| 0x004 | `CTRL` | RW | enable, `src_sel`, `ts_clear`, soft reset |
| 0x008 | `STATUS` | RO / W1C | `adc_err`, `rp_underflow`, `job_late`, `queue_ovf`, `mbox_ovf`, `dma_ovr` |
| 0x010 / 0x014 | `TS_LO` / `TS_HI` | RO | Current `ts64` |
| 0x018 – 0x020 | `PPS_TS_LO/HI`, `PPS_FRAC` | RO | Last PPS timestamp |
| 0x040 – 0x05C | `GAIN[4]`, `OFFSET[4]` | RW | `fe_cond` |
| 0x060 – 0x06C | `DTRIM[4]` | RW | `fe_cond` |
| 0x070 – 0x07C | `SW_THR`, `SW_WIN`, `SW_HOLDOFF`, `DT_MIN` | RW | `sw_det` |
| 0x080 – 0x08C | `MB_THR`, `MB_WIN`, `MB_SETTLE`, `MB_HOLDOFF` | RW | `mb_det` |
| 0x090 – 0x09C | `WIN_SW`, `WIN_MB`, `MAXLAG`, `TAPER` | RW | Windows and lags |
| 0x0A0 – 0x0AC | `BAND_SW`, `BAND_MB`, `PHAT_EPS` | RW | `xspec_phat` |
| 0x0B0 | `LATENCY` | RO | `lat_fe`, `lat_dec` |
| 0x0C0 / 0x0C4 | `MBOX_WR_PTR` / `MBOX_RD_PTR` | RO / RW | Mailbox ring |
| 0x0D0 – 0x0DC | `DMA_BASE`, `DMA_SIZE`, `DMA_WR_PTR`, `DMA_CTRL` | RW / RO | `raw_dma` |
| 0x0E0 – 0x0F4 | Counters: frames, sw_evts, mb_evts, jobs, drops, adc_errs | RO | Diagnostics |
| 0x100 – 0x1FF | `RP_FIFO` | WO | `replay_src` frame writes |

### 4.2 `raw_dma` — continuous raw logging (optional)

| Port | Dir | Width | Description |
|---|---|---|---|
| `src_ch0..3`, `src_valid` | in | s16, 1 | Raw (pre-`fe_cond`) frames |
| `cfg_base`, `cfg_size`, `cfg_en` | in | 32, 32, 1 | DDR3 circular buffer |
| `avm_*` | out | 32 addr, 64 data, burst 8 | Write master on an F2SDRAM port |
| `dma_wr_ptr` | out | 32 | To `csr` |
| `dma_ovr` | out | 1 | Sticky |

| Parameter | Value |
|---|---|
| Rate | 800 kB/s sustained (one 64-bit frame per 10 µs) |
| Burst | 8 frames = 64 B |
| FIFO | ≥ 128 frames (1.28 ms) |
| Buffer | e.g. 8 MB = 10 s, reserved on the HPS via device-tree `reserved-memory` |

---

## 5. Resource estimate (5CSXFC6, estimate only)

| Resource | Use | Available |
|---|---|---|
| M10K — rings | ≈ 64 | 553 |
| M10K — FFT buffers (5) | ≈ 40 | |
| M10K — mailbox (64 KB) | ≈ 52 | |
| M10K — FIFOs, ROMs | ≈ 10 | |
| **M10K total** | **≈ 165 (30 %)** | |
| DSP — FFT | ≈ 36 | 112 |
| DSP — `fe_cond`, decimator, detectors, PHAT | ≈ 12 | |
| **DSP total** | **≈ 48 (43 %)** | |

---

## 6. Verification plan

| Level | Method | Pass criterion |
|---|---|---|
| Block | Self-checking testbench per block | Cycle counts within the budgets above |
| `fft2048` | Random + tone vectors vs NumPy | Output SNR consistent with ≤ 1 µs TDOA error |
| Engine (10–13) | Windows from `gs_det_signal_prepare.py` output vs `gs_det_tdoa.py` | Coarse lags identical; refined TDOA within 1 µs |
| System (sim) | Full RTL with `replay_src` fed from `gs_gen_apply_adc.py` WAVs | Onsets and TDOAs match the Python pipeline |
| System (board) | Same WAVs through `replay_src` on hardware | Same as above, plus no `STATUS` errors over a 1 h soak |
| Burst | Synthetic 900 rpm sequences | No `job_late`, `queue_ovf` or `mbox_ovf` |
| Live ADC | Calibration source at known bearing | Residual skew after `dtrim` ≤ 0.5 µs |

---

## 7. Mapping to existing RTL

| Existing | Becomes | Changes |
|---|---|---|
| `fft1024-claude/fft1024.v` | `fft2048` | N = 2048, block-floating-point, multi-buffer select |
| `fromClaude/gcc_phat_top.v` → `rfft1024_packed` | `gcc_loader` + unpack in `xspec_phat` | 4 channels in 2 packed FFTs |
| `fromClaude/gcc_phat_core.v` | `xspec_phat` | 6 pairs, band mask, IFFT packing |
| `fromClaude/peak_detector_tdoa.v` | `peak_search` | ± max_lag scan, PSR, lag window out; drop `tdoa_calc` |
| `gcc_phat_top.v` parameters | — | FS 16000 → 100000, NFFT 1024 → 2048, D_MM 150 → 300 |
