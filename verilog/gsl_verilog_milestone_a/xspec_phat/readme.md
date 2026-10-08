# xspec_phat — channel unpack, 6 cross-spectra, PHAT, IFFT packing

Implements block 12 of `verilog/fpga_blocks.md` (section 3.4). Sits between the two forward FFTs and the three IFFTs of a GCC-PHAT job.

- **Status:** simulated in Icarus Verilog 12 together with `fft2048`; every stage is bit-exact against the Python model, and the fixed-point TDOA error is ≤ 0.2 µs. Not yet through Quartus.
- **Depends on:** `../fft2048` (engine, twiddles and the v2 buffer bank).
- **Version:** 2026-10-08

---

## Files

| File | Description |
|---|---|
| `xspec_phat.v` | Top: Z reads, unpack, pair scheduling, band mask, IFFT packing, W writes |
| `xspec_phat_unit.v` | Pipelined cross-spectrum + PHAT normalisation (1 pair / cycle) |
| `gen_rsqrt_rom.py` | Generates `rsqrt1024.hex` (reciprocal-sqrt ROM) |
| `rsqrt1024.hex` | ROM init file (pre-generated) |
| `xspec_phat_model.py` | Bit-true model of `xspec_phat` + the full fixed-point GCC-PHAT chain, float references, checker |
| `tb_xspec_phat.v` | System testbench: fft2048 ×2 → xspec_phat → IFFT ×3 → coarse peaks |
| `Makefile` | Icarus flow (uses `../fft2048`) |

## Quick start

```bash
cd verilog/xspec_phat
make          # ROM -> simulate -> bit-true + TDOA accuracy check (needs numpy)
```

---

## Interface

| Port | Dir | Width | Description |
|---|---|---|---|
| `xp_start` | in | 1 | Strobe |
| `cfg_band_lo`, `cfg_band_hi` | in | 11 | Keep bins k in [lo, hi] (and mirrors); 0 / 1024 = all |
| `cfg_floor` | in | 7 | Zero G when floor(log₂\|G\|) < floor; 0 = off |
| `z0_exp`, `z1_exp` | in | 5 | `fft_bfp_exp` of the Z0 and Z1 transforms (used only by the floor) |
| `xp_busy`, `xp_done` | out | 1 | Busy level, done strobe |
| `z_re`, `z_raddr` | out | 1, 11 | Read Z0 and Z1 at the same address (natural order) |
| `z0_rdata`, `z1_rdata` | in | 36 | From the bank client ports of buffers 0 and 1 |
| `w_we` | out | 3 | One-hot write enable for W0 / W1 / W2 |
| `w_waddr`, `w_wdata` | out | 11, 36 | Shared write address (bit-reversed) and data |

## Data contract

| Item | Rule |
|---|---|
| Input | Z0 = FFT(ch0 + j·ch1), Z1 = FFT(ch2 + j·ch3), natural order, straight from `fft2048` |
| Output | W0 = G01 + j·G02, W1 = G03 + j·G12, W2 = G13 + j·G23, at **bit-reversed** addresses (`fft2048` input contract) |
| Output level | \|G\| ≈ 32767 (0.25 FS), so W components ≤ 65534 < 0.5 FS |
| Lag sign | G_ab = X_b · X_a*, so lag > 0 means channel b lags channel a — same as `gs_det_tdoa.py` (`G = X2 · conj(X1)`) |
| IFFT result | IFFT(W_p) = r_A + j·r_B: real part is the first pair's correlation, imaginary part the second's |

---

## Design

### Unpack

    X_A[k] = (Z[k] + conj(Z[N−k])) / 2
    X_B[k] = (Z[k] − conj(Z[N−k])) / 2j

The /2 uses round-half-up so X stays s18 and fits the 18×18 multipliers.

### Cross-spectrum and PHAT (`xspec_phat_unit`, latency 7)

| Stage | Operation |
|---|---|
| S1 | G = X_b · conj(X_a), 37 bits |
| S2 | One's-complement magnitude OR → leading-one position; zero / floor test |
| S3 | Normalise G so its largest component is in [2^16, 2^17) |
| S4 | q = \|G_n\|², in [2^32, 2^35] |
| S5 | q = m · 2^(2e), e ∈ {16, 17}; ROM lookup y = 1/√m (10-bit index) |
| S6 | G_n · y |
| S7 | Convergent rounding by e + 2, clamp to ±32767 |

Both components are scaled by the same real factor, so **phase comes only from G_n and is exact**. The ROM's ~0.1% error affects only magnitude.

### Exponents and floor

The BFP exponents of Z0 and Z1 don't change the result (PHAT keeps only phase). They are used only by the optional floor:

    zero G  if  msb(|G|) + exp(Xa's spectrum) + exp(Xb's spectrum) < cfg_floor

### Packing (Hermitian completion)

    W[k]   = (GA.re − GB.im) + j(GA.im + GB.re)
    W[N−k] = (GA.re + GB.im) + j(GB.re − GA.im)

### Schedule and timing

One 6-cycle slot per bin k = 0 … 1024:

| Slot s, cycle | Operation |
|---|---|
| 0 / 1 | Read Z0, Z1 at k = s, then at N − s |
| 2 | Unpack → X_next |
| (slot s+1) 0 … 5 | Issue pairs (0,1) (0,2) (0,3) (1,2) (1,3) (2,3) of bin s |
| unit output, odd pair | Write W[k]; one cycle later write W[N−k] |

At most one W write per cycle.

    1026 slots × 6 + 11 drain = 6,167 cycles  (61.7 µs @ 100 MHz)

---

## Test results (Icarus Verilog 12)

| Test | Path | Coarse lags | Stages bit-exact | Fixed-point TDOA error |
|---|---|---|---|---|
| 1. Noise, integer delays (0, 7, −12, 25) | SW, 100 kHz | all exact | Z, W, R ✓ | 0.000 µs |
| 2. Plane-wave N-wave, az 37° el 10°, noise A/100 | SW, 100 kHz | within 1 sample | Z, W, R ✓ | 0.000 µs |
| 3. Plane-wave Friedlander, az 200° el 5° | MB, 20 kHz | within 1 sample | Z, W, R ✓ | 0.195 µs |
| 4. Dead mic (ch3 idle ±1 LSB), band 20–500, floor 27 | SW, 100 kHz | live pairs exact; dead pairs max\|r\| ≤ 1 | Z, W, R ✓ | 0.000 µs |

- **Cycles:** `xspec_phat` 6,167 every test (budget 6,200); `fft2048` 11,321 per transform.
- **Fixed-point error:** RTL correlations vs float GCC-PHAT on the same windows, both refined with the planned windowed-sinc method on a 1/256-sample grid. 0.195 µs is one grid step at 20 kHz, i.e. the measurement limit. Budget: 1 µs.
- **Band mask:** test 4 has zero non-zero bins outside [20, 500].
- **No client-port conflicts** in any test.

---

## Findings to act on

### 1. Band-limit PHAT (big accuracy gain)

With all bins enabled (what `gs_det_tdoa.py` does today), PHAT gives noise-only bins the same weight as signal bins. FPGA path, same test signals, |TDOA − truth|:

| Signal | All bins | Band-limited | Band |
|---|---|---|---|
| MB Friedlander (test 3) | 16.1 µs worst, 11.6 mean | **4.7 µs worst, 2.2 mean** | 10–205 (≈ 0.1–2 kHz @ 20 kHz) |
| SW N-wave (test 2) | 1.3 µs worst, 0.6 mean | **0.7 µs worst, 0.3 mean** | 20–512 (≈ 1–25 kHz @ 100 kHz) |

The MB all-bins case exceeds the 6 µs estimator allocation in `timing_budget.md`; band-limited, it fits. Recommended defaults: SW `lo = 20, hi = 512`; MB `lo = 10, hi = 205`. Add the same mask to `gs_det_tdoa.gcc_phat` so Python and FPGA stay comparable.

### 2. Dead or unplugged microphone

PHAT normalises any non-zero bin to full magnitude, so an idle channel's ±1 LSB noise produces correlations about 40% as strong as real ones (spurious TDOAs). With `cfg_floor = 27` they drop to ≤ 1 LSB. The right floor depends on the real ADC noise level, so calibrate it on hardware; 0 (off) reproduces Python behaviour.

---

## Resources (indicative, yosys `synth_intel_alm -family cyclonev`)

| Block | M10K | 18×18 multipliers | LUTs | FFs |
|---|---|---|---|---|
| `xspec_phat` | 2 (rsqrt ROM) | 8 (= 4 DSP blocks) | ≈ 1,570 | ≈ 1,080 |
| `fft_buf_bank` v2 (5 buffers) | 40 | 0 | — | — |

## Changes vs `fpga_blocks.md` section 3.4

| Spec | Implementation | Reason |
|---|---|---|
| G = X_a · X_b* | G = X_b · X_a* | Matches `gs_det_tdoa.py` lag sign |
| rsqrt: ROM seed + 1 Newton step, error ≤ 2⁻¹² | ROM only, ≤ 0.1% | Error is magnitude-only (phase exact); saves 3 multipliers |
| `cfg_eps` | `cfg_floor` (log₂ threshold using BFP exponents) | Scale-independent; needs `z0_exp`, `z1_exp` |
| Unit-magnitude PHAT, scaled ≤ 0.25 FS | Magnitude 32767 | Keeps W < 0.5 FS for the IFFT |
| ≤ 6,200 cycles, ≈ 6 DSP | 6,167 cycles, 4 DSP blocks | — |

## Integration notes (for `gcc_engine`)

- **Port mux:** while `xp_busy`, route `z_raddr` / `z_re` to client ports 0 and 1, and `w_we` / `w_waddr` / `w_wdata` to client ports 2–4. `tb_xspec_phat.v` shows the exact mux.
- **Exponents:** latch `fft_bfp_exp` after each forward FFT and feed them to `z0_exp` / `z1_exp`.
- **Downstream:** after the 3 IFFTs, `peak_search` reads W0–W2 in natural order; the real part is pair A, the imaginary part pair B.
