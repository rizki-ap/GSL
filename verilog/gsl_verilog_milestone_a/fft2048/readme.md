# fft2048 — 2048-point FFT/IFFT engine with block-floating-point

Implements block 11 of `verilog/fpga_blocks.md` (section 3.3), plus the buffer memory it runs in.

- **Status:** simulated in Icarus Verilog 12; RTL is bit-exact against the Python model on all test vectors. Not yet through Quartus (timing closure at 100 MHz unverified).
- **Version:** 2026-10-08

---

## Files

| File | Description |
|---|---|
| `fft2048.v` | Engine: address generation, pipeline, BFP control, FSM |
| `fft2048_bfly.v` | Radix-2 butterfly: complex multiply, convergent rounding, BFP shift, saturation |
| `fft2048_twiddle_rom.v` | 1024-entry twiddle ROM (W = cos − j·sin, Q1.17) |
| `fft_buf_bank.v` | NBUF (default 5) complex buffers, each split into 2 banks; engine port + per-buffer client ports (v2) |
| `fft_sdp_ram.v` | 1024 × 36 simple dual-port RAM (M10K) |
| `gen_twiddle2048.py` | Generates `tw2048_re.hex`, `tw2048_im.hex` |
| `tw2048_re.hex`, `tw2048_im.hex` | Twiddle ROM init files (pre-generated) |
| `fft2048_model.py` | Bit-true Python model + RTL checker |
| `tb_fft2048.v` | Self-checking testbench (DC, 1 tone, 2 tones, arbitrary complex, IFFT round trip) |
| `Makefile` | Icarus flow |
| `do_fft2048.do` | ModelSim / Questa script |

---

## Quick start

```bash
cd verilog/fft2048
make          # twiddles -> simulate -> bit-true check
make wave     # same, with tb_fft2048.vcd for GTKWave
make model    # Python model vs numpy (needs numpy)
```

ModelSim / Questa: `python3 gen_twiddle2048.py` then `vsim -c -do do_fft2048.do`.

---

## Interface

### `fft2048`

| Port | Dir | Width | Description |
|---|---|---|---|
| `clk`, `rst_n` | in | 1 | 100 MHz, active-low async reset |
| `fft_start` | in | 1 | Strobe; ignored while busy |
| `fft_inverse` | in | 1 | 0 = FFT, 1 = IFFT (unnormalised, no 1/N) |
| `fft_buf_sel` | in | 3 | Buffer to transform in place (Z0, Z1, W0, W1, W2) |
| `fft_busy` | out | 1 | High while the engine owns the buffer |
| `fft_done` | out | 1 | 1-cycle strobe |
| `fft_bfp_exp` | out | 5 | Total right shifts; true X[k] = out[k] · 2^exp |
| `fft_err` | out | 1 | Input contract violated or saturation occurred (valid at `fft_done`) |
| `mem_*` | out / in | — | Bank-level memory master → `fft_buf_bank` engine port |

### `fft_buf_bank` (v2)

| Port | Dir | Width | Description |
|---|---|---|---|
| `eng_*` | in / out | — | Engine port (connect to `fft2048` `mem_*`) |
| `cl_re[i]` | in | 1 | Client read strobe for buffer i (used for conflict detection) |
| `cl_raddr[i]` | in | 11 | Client read word address |
| `cl_rdata[i]` | out | 36 | Read data, 1 cycle after `cl_raddr` |
| `cl_we[i]` | in | 1 | Client write enable |
| `cl_waddr[i]`, `cl_wdata[i]` | in | 11, 36 | Client write address / data `{re, im}` |
| `cl_conflict` | out | 1 | A client touched the buffer the engine currently owns |

Client ports are flattened (buffer i at `[i*W +: W]`) and give 1 read + 1 write per buffer per cycle, so different buffers can be accessed in parallel (e.g. `xspec_phat` reads Z0 and Z1 while writing W0–W2). Which block drives each client port is muxed outside the bank, in `gcc_engine`.

---

## Data contract

| Item | Rule |
|---|---|
| Word format | `{re[17:0], im[17:0]}`, each s18 Q1.17 |
| Input order | **Bit-reversed**: write x[n] at address `bitrev11(n)` |
| Output order | Natural: X[k] at address k |
| Input range | \|re\|, \|im\| < 2^16 (< 0.5 FS). Violations set `fft_err` |
| Output scaling | X[k] = out[k] · 2^`fft_bfp_exp` |
| IFFT | Unnormalised: IFFT(FFT(x)) = N · x |

Producers get bit-reversed order for free by wiring their write address reversed (`gcc_loader` for the forward FFT, `xspec_phat` for the IFFT).

---

## Design

### Memory banking

Each buffer is two 1024 × 36 simple dual-port RAMs:

    bank(addr) = XOR of all 11 address bits
    row(addr)  = addr[10:1]

The two operands of every radix-2 butterfly differ in exactly one address bit, so they always land in different banks. Each bank sees one read (butterfly n) and one write (butterfly n − 4) per cycle, which a single M10K handles.

### Butterfly and rounding

    P = (A + W·B) / 2^sh,   Q = (A − W·B) / 2^sh

A·2^17 ± W·B is kept at full 38-bit precision, then rounded once (convergent, round-half-to-even) by 17 + sh bits and saturated to s18. For k = 0 (W = 1, not representable in Q1.17) the multiplier is bypassed, so W = 1 and W = −j are both exact.

### Block-floating-point

After each stage, the engine ORs the one's-complement magnitudes of every value written. With p = MSB position of that OR, the next stage shifts by sh = max(0, p − 14). Per-component radix-2 growth is at most 1 + √2, so this can never overflow. Stage 0 has no BFP decision, which is why the input contract requires < 0.5 FS (stage 0 has W = 1, growth ≤ 2).

58 adversarial full-scale inputs (constant max, alternating, Nyquist, random-sign max, 45° tones), FFT and IFFT, produced zero saturations in the bit-true model.

### Inverse

IDFT(x) = swap(DFT(swap(x))), where swap exchanges re and im. The swap is applied on stage-0 reads and last-stage writes, so the twiddle ROM is forward-only and costs nothing extra.

### Pipeline and timing

| Cycle | Operation |
|---|---|
| t | Address generation → RAM / ROM read addresses |
| t+1 | RAM / ROM data, bank un-swap, inverse swap → input registers |
| t+2 | Butterfly M (multiply) |
| t+3 | Butterfly A (add, round, saturate) |
| t+4 | Bank swap → RAM write, BFP max update |

One butterfly per cycle, 1024 per stage, plus a 5-cycle drain per stage (the BFP decision needs the whole stage):

    11 stages × 1029 + 2 = 11,321 cycles = 113.2 µs @ 100 MHz

---

## Test results (Icarus Verilog 12)

| Test | Buffer | Cycles | bfp_exp | SNR vs double DFT | Max bin error | Checks |
|---|---|---|---|---|---|---|
| 1. DC, 0.25 FS | 0 | 11,321 | 11 | 294.8 dB (exact) | 0.00 LSB | peak at bin 0, X[0] = N·A |
| 2. 0.40 cos, bin 37 | 1 | 11,321 | 10 | 101.4 dB | 0.20 LSB | peaks at 37 / 2011, amplitude |
| 3. 0.25 cos bin 100 + 0.15 sin bin 333 | 2 | 11,321 | 10 | 98.5 dB | 0.21 LSB | 4 peaks, amplitude, sin → −j phase |
| 4. Complex chirp 20→600 + noise | 3 | 11,321 | 7 | 82.1 dB | 3.57 LSB | SNR |
| 5. IFFT of test 4 spectrum | 4 | 11,321 | 4 | 82.8 dB | 3.97 LSB | round trip 78.0 dB vs original x |

All five outputs are bit-exact against `fft2048_model.py`. A mutation test (wrong twiddle sign in the butterfly) fails tests 2–5, so the checks are not vacuous; DC still passes because k = 0 bypasses the multiplier.

**SNR threshold (75 dB):** it sits ≥ 15 dB below the best input SNR the array can deliver (16-bit ADC ≈ 90 dB; acoustic SNR is far lower), so FFT rounding never limits the TDOA estimate. Broadband full-scale input measures ~79–82 dB; tones ~100 dB.

---

## Resources (indicative, yosys `synth_intel_alm -family cyclonev`)

| Block | M10K | 18×18 multipliers | LUTs | FFs |
|---|---|---|---|---|
| `fft2048` (engine + twiddle ROM) | 4 | 4 (= 2 DSP blocks) | ≈ 1,350 | ≈ 420 |
| `fft_buf_bank` (5 buffers) | 40 | 0 | ≈ 1,100 (muxes) | — |

Quartus numbers will differ; check after the first fit.

---

## Changes vs `fpga_blocks.md` section 3.3

| Spec | Implementation | Reason |
|---|---|---|
| ≤ 15,360 cycles (2N bit-reverse + N/2·log₂N) | 11,321 cycles | No bit-reverse pass; producers write bit-reversed |
| Buffers "true dual-port, 2R + 2W per cycle" | 2 banks per buffer, each simple dual-port | One M10K has only 2 ports; parity banking gives conflict-free 2R + 2W |
| Input range unspecified | \|x\| < 0.5 FS | Stage 0 has no BFP decision |
| IFFT via conjugate twiddles | IFFT via re/im swap | Forward-only ROM; −sin = −1.0 is representable but +1.0 is not |

## Integration notes

- **`xspec_phat` scaling:** each packed IFFT input W_p = G_a + j·G_b must have components < 0.5 FS, so scale the unit-magnitude PHAT output to ≤ 0.25 FS.
- **Client ports (v2):** the original single external port was replaced by per-buffer client ports so `xspec_phat` can access Z0, Z1 and W0–W2 in the same cycle. `tb_fft2048.v` was updated accordingly; all results are unchanged and bit-exact.
- **Timing closure:** if 100 MHz fails, the likely critical path is the butterfly A stage (38-bit add + rounding + saturate). Split it into two registers and add one stage to the control pipeline (v1..v5); the drain grows to 6 cycles (+11 cycles per transform).
