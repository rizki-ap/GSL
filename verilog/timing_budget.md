# Gunshot Locator — Timing Budget

- **Target platform:** Terasic DE10-Standard (Cyclone V SX 5CSXFC6D6F31C6: FPGA fabric + dual-core Cortex-A9 HPS)
- **Reference product:** Metravib PILAR-V
- **Derived from:** `python/det_config.ini`, `python/gen_config.ini`, `python/project.md`, `verilog/`, `reference/Pilar Accoustic Gunshot 2004.pdf`
- **Version:** 2026-10-08

Items marked **(verify)** are assumptions, not sourced figures.

---

## 1. Summary

| Item | Value |
|---|---|
| Sample rate | 100 kHz per channel, 4 channels, simultaneous sampling |
| TDOA error budget | ≤ 10 µs per mic pair (95 %) → ≈ 0.9° azimuth |
| Δt (SW → MB) error budget | ≤ 0.5 ms |
| FFT size | 2048 for both paths: SW at 100 kHz, MB decimated to 20 kHz |
| FPGA compute per event | ≈ 88k cycles = 0.88 ms @ 100 MHz (budget 2 ms) |
| Processing after MB window closes | ≈ 70 ms worst case |
| Burst capability | 900 rpm (66.7 ms/shot) at ≈ 3 % FPGA engine load |

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

**Δt allocation: ≤ 0.5 ms.** That is < 3 % at 50 m for α ≥ 30°, and < 1 % at 200 m. It is easy for an energy onset detector. However, every constant filter delay in the SW and MB paths (FIR group delay, decimator delay) must be added back to t_SW and t_MB before forming Δt. An uncompensated 0.3 ms decimator delay alone would cost 1.5 % at α = 30°, 50 m.

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

---

## 6. Block budgets

**Clock assumption:** 100 MHz fabric clock from a PLL (the board oscillator is 50 MHz). At 50 MHz every compute time below doubles and still fits.

### 6.1 FPGA — stream rate (hard real time)

One 4-channel sample frame arrives every 10 µs = 1000 clock cycles.

| Block | Throughput | Latency budget | Notes |
|---|---|---|---|
| ADC interface (e.g. AD7606B) | 4 ch / 10 µs | conversion + readout ≤ 5 µs | CONVST from PLL; serial 2-line readout at 25 MHz ≈ 1.3 µs |
| Sample counter / timestamp | 1 / frame | 0 | 64-bit; latched on each detected onset |
| HP + per-channel gain/delay trim | 4 ch / frame | ≤ 0.5 ms, constant | one shared DSP gives 250 MAC / ch / frame |
| SW onset detector | 4 ch / frame | trigger ≤ 0.5 ms after onset | amplitude/energy vs 5 ms noise floor (`sw_threshold_mult = 10`) |
| MB decimator ÷ 5 | 4 ch / 10 µs in, 4 ch / 50 µs out | ≈ 0.3 ms group delay, constant | ~61-tap polyphase FIR; compensate in t_MB |
| MB onset detector | 4 ch / 50 µs | trigger ≤ 1 ms after onset | armed after SW; adaptive energy-settling FSM |
| SW ring buffer | continuous | — | 4096 × 18 bit per channel = 41 ms |
| MB ring buffer | continuous | — | 4096 × 18 bit per channel = 205 ms |

The ring buffers replace ping-pong buffers. Each holds ≥ 2× its window, so a window can be read into the FFT right after it closes while acquisition continues. Overlapping events are just multiple start pointers into the same ring.

### 6.2 FPGA — event rate (per SW or MB window)

Cycle model of the current radix-2 core (`verilog/fft1024-claude`): 2N (bit-reverse) + (N/2)·log₂N (butterflies). For N = 2048 that is 4096 + 11,264 = 15,360 cycles.

| Step | Operations | Cycles | @ 100 MHz | @ 50 MHz |
|---|---|---|---|---|
| Load 4 ch into 2 packed complex inputs | 4 × 2048 samples, 2 / cycle | 4.1k | 41 µs | 82 µs |
| Forward FFT | 2 × FFT-2048 | 30.7k | 307 µs | 614 µs |
| Unpack + cross-spectrum + PHAT | 6 pairs × 1025 bins | 6.2k | 62 µs | 123 µs |
| Inverse FFT | 3 × IFFT-2048 (2 pairs each) | 46.1k | 461 µs | 922 µs |
| Peak search (± max_lag only) | 6 × 201 lags | 1.2k | 12 µs | 24 µs |
| **Total** | | **88.3k** | **0.88 ms** | **1.77 ms** |
| **Budget** | | | **2 ms** | **4 ms** |

**Packing tricks:**

- **Forward:** two real channels per complex FFT (re = ch_a, im = ch_b), as in `rfft1024_packed`. Unpack with X_a[k] = (Z[k] + Z*[N−k]) / 2 and X_b[k] = (Z[k] − Z*[N−k]) / 2j.
- **Inverse:** each GCC-PHAT correlation is real, so IFFT(G₁ + j·G₂) = r₁ + j·r₂. That puts 6 pairs into 3 IFFTs.

**Fixed point.** The current core is 18-bit Q1.17 with an unconditional ÷ 2 per stage. That costs dynamic range on short, low-energy windows (a ~1 ms SW pulse inside an 18 ms window). Use block-floating-point (per-stage overflow detect), and confirm the 1 µs allocation by bit-true comparison with `gs_det_tdoa.py`.

### 6.3 FPGA resources (estimate)

| Resource | Use | Available (5CSXFC6) |
|---|---|---|
| M10K — ring buffers | 4 ch × 2 paths × 8 blocks = 64 | 553 |
| M10K — FFT work memory + spectra | ≈ 50 | |
| M10K — total | ≈ 115 (≈ 21 %) | |
| DSP — FFT | ≈ 36 (per `fft1024` readme) | 112 |
| DSP — filters, detectors | ≈ 8 | |

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

Times are measured from SW arrival at the array (t = 0), using worst-case budgets.

| Event | 200 m (Δt ≈ 324 ms) | 1500 m (Δt ≈ 2.60 s) |
|---|---|---|
| SW arrives | 0 | 0 |
| SW trigger | 0.5 ms | 0.5 ms |
| SW window closed | 15 ms | 15 ms |
| SW TDOAs (FPGA) | 17 ms | 17 ms |
| Early alert: shot detected + SW direction | ≈ 85 ms | ≈ 85 ms |
| MB arrives | 324 ms | 2.60 s |
| MB window closed (+ 80 ms) | 404 ms | 2.68 s |
| MB TDOAs (FPGA) | 406 ms | 2.68 s |
| Final alert: bearing + range | ≈ 475 ms | ≈ 2.75 s |

- Processing adds ≈ 70 ms after the MB window closes. Everything else is acoustic propagation.
- With a single array, the early alert's trajectory has a two-fold ambiguity until the MB arrives (PILAR 2004).
- If no MB is detected (suppressed weapon, or MB below the noise floor), the SW-only report waits for `dt_max_s` = 5 s. Consider issuing a provisional report earlier.

---

## 8. Burst fire (900 rpm)

- **Compute:** 2 events × 0.88 ms per 66.7 ms shot ≈ 2.6 % of one FFT engine at 100 MHz. Throughput is not a constraint.
- **Window overlap:** the 86 ms MB window is longer than the shot interval. It will contain the next shot's MB and, at 200 m, the SW of shot n+5. Trim `mb_window_post_s` to ≈ 40 ms (the Friedlander positive phase is only 2–4 ms), and/or blank intervals around detected SWs before GCC-PHAT.
- **Association:** up to `dt_max_s` / 66.7 ms ≈ 75 SW events can be waiting for an MB at once. The HPS queue must hold that many.
- **MB detector:** the adaptive energy-settling search (`project.md`, Known Limitation 5) assumes an isolated event. Under burst fire the energy never settles, so it needs a burst-aware mode.

---

## 9. Repo items that conflict with this budget

- [ ] `devices/readme.md`: 16 kHz sampling and 7 kHz LPF → 100 kHz and a 4th-order ~25 kHz anti-alias filter.
- [ ] `devices/readme.md`: 1 µF / 10 kΩ coupling (16 Hz corner) → ≤ 2 Hz, or per-channel calibration (4.3.1).
- [ ] `devices/readme.md`: MCP6002 (1 MHz GBW) → ≥ 10 MHz GBW op-amp to remove the gain-stage phase spread (4.3.1).
- [ ] Onboard ADC: single multiplexed 12-bit converter (`devices/readme.md` assumes ADC128S022 — confirm the part on your board revision). The skew is calibratable; the ~74 dB dynamic range is the larger problem → external 16-bit simultaneous-sampling ADC.
- [ ] `python/project.md`: ADS8688 is listed as simultaneous-sampling, but it is a multiplexed 8-channel SAR. Simultaneous options: AD7606B, ADS8588S.
- [ ] `verilog/fromClaude/gcc_phat_top.v`: FS = 16000, NFFT = 1024, D_MM = 150 → 100000 / 2048 / 300.
- [ ] `verilog/fromClaude/peak_detector_tdoa.v`: `tdoa_calc` computes a single-pair 2D angle → output raw lags only and do the 3D solve on the HPS.
- [ ] `verilog/fft1024-claude`: generalise to N = 2048, add block-floating-point scaling, and pack 2 pairs per IFFT.
- [ ] `python/gs_det_tdoa.py`: add the windowed-sinc refinement alongside the 16× zero-padded IFFT, so hardware and simulation use the same estimator.
- [ ] `python/det_config.ini`: `mb_window_post_s = 0.080` → ≈ 0.040 for burst operation (Section 8).
- [ ] `python/gs_det_*`: compensate constant filter/decimator delays before forming Δt (4.4).

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
