# Milestone C1 — processor interface specification (Nios V)

- **Processor:** Nios V/m (bare-metal C), Quartus Prime Lite 22.1 or newer, DE10-Standard
- **Status:** implemented and verified in simulation (`gsl_mm.v`, `gsl_mbox.v`, `gsl_replay_dma.v`; see `readme.md`)
- **Version:** 2026-10-12

The FPGA side stays processor-neutral (Avalon-MM), so this interface would also work from the HPS lightweight bridge if that is ever wanted.

---

## 1. System (Platform Designer)

```
CLOCK_50 ─ PLL ─┬─ clk_sys 100 MHz ─────────────────────────────────────────┐
                └─ sdram_clk 100 MHz (phase per Terasic DE10-Standard demos)│
                                                                             │
 Nios V/m ──┐                                                                │
 JTAG master┼── Avalon-MM interconnect ──┬── on-chip RAM 128 KB (code, data) │
 (System    │                            ├── JTAG UART (results, logs)       │
  Console)  │                            ├── interval timer (profiling)      │
            │                            ├── sysid                           │
            │                            ├── SDRAM controller ── 64 MB SDRAM │
            │                            └── gsl_mm slave (64 KB window)     │
            │                                                                │
 gsl_mm DMA master ──────────────────────── SDRAM controller (recordings)    │
 gsl_mm IRQ ─────────────────────────────── Nios V                           │
```

| Component | Purpose |
|---|---|
| Nios V/m | Runs the processing software (C2) |
| On-chip RAM 128 KB | Code and data (Nios V/m has no cache, so running from SDRAM would be slow) |
| SDRAM controller, 64 MB | Recordings for replay (all 15 scenarios ≈ 11 MB) |
| JTAG-to-Avalon master | System Console loads recordings, reads and writes registers for debug |
| JTAG UART | Result and log output to the host |
| Interval timer | Measures per-event processing time |
| `gsl_mm` | `gsl_core` + registers + result mailbox + replay DMA (this document) |

The concrete system for the DE10-Standard (instance names, connections, base addresses: `gsl_mm` at 0x0800_0000, SDRAM at 0x0C00_0000) is in `../board/readme.md`, step 2.

---

## 2. `gsl_mm` — what is new in RTL

| Block | Function |
|---|---|
| `gsl_csr` | Control, status, counters, the 43 configuration values, replay and mailbox registers |
| `gsl_mbox` | Serialises SW / MB records, job dispatches, late jobs and engine results (with lag windows) into a ring of 32-bit words the processor reads; IRQ while not empty |
| `gsl_replay_dma` | Avalon-MM read master: streams frames from SDRAM into `replay_src`, with loop mode for the soak test; stops the frame cadence cleanly at the end of a non-looping recording |
| `gsl_core` | Unchanged from Milestone B |

---

## 3. Register map (Avalon-MM slave, 32-bit words, byte offsets)

### 3.1 Control, status, counters (0x000 – 0x0FF)

| Offset | Name | Access | Bits |
|---|---|---|---|
| 0x000 | `ID` | RO | `0x47534C31` ("GSL1") |
| 0x004 | `VERSION` | RO | Build date, BCD `0xYYYYMMDD` |
| 0x008 | `CTRL` | RW | [0] run, [1] src_sel (1 = replay), [2] soft_reset (self-clearing: clears `gsl_core` state, mailbox, counters), [3] irq_en, [4] replay_loop |
| 0x00C | `STATUS` | RO / W1C | [0] rp_underflow (sticky), [1] sw_noise_valid, [2] mb_noise_valid, [3] mbox_drop (sticky), [4] res_err seen (sticky), [5] engine busy, [6] replay_done, [7] dma_error (sticky) |
| 0x010 | `TS_LO` | RO | Current sample index, low word (reading it latches `TS_HI`) |
| 0x014 | `TS_HI` | RO | High word |
| 0x018 | `LATE_CNT` | RO | Jobs dropped as late (`evt_sched`) |
| 0x01C | `OVF_CNT` | RO | Events dropped, pending FIFO full (`evt_sched`) |
| 0x020 | `MBOX_DROP_CNT` | RO | Messages dropped, mailbox full |
| 0x024 | `EVENT_CNT` | RO | [15:0] SW records, [31:16] MB records |
| 0x028 | `RESULT_CNT` | RO | Engine results |
| 0x030 | `MBOX_WR` | RO | Write pointer (word index, 0 … 8191) |
| 0x034 | `MBOX_RD` | RW | Read pointer; the processor writes it after consuming a message |
| 0x038 | `MBOX_LEVEL` | RO | Words waiting |
| 0x040 | `RP_BASE` | RW | SDRAM byte address of frame 0 (8-byte aligned) |
| 0x044 | `RP_NFRAMES` | RW | Frames in the recording |
| 0x048 | `RP_POS` | RO | Frames delivered in the current pass |
| 0x04C | `RP_LOOPS` | RO | Completed passes (loop mode) |
| 0x050 | `RP_LEVEL` | RO | Replay FIFO level |

### 3.2 Configuration (0x100 – 0x1AB)

One register per value, **in the same order as `mb_cfg.txt`** (Milestone B), at `0x100 + 4·k` for k = 0 … 42. Reset values are the Milestone B defaults (`stream/milestone_b.py`, `CFG`). Signed values are written sign-extended.

| k | Offset | Value | | k | Offset | Value |
|---|---|---|---|---|---|---|
| 0–3 | 0x100–0x10C | `fe_gain[0..3]` | | 24 | 0x160 | `mb_settle_timeout` |
| 4–7 | 0x110–0x11C | `fe_offset[0..3]` | | 25–27 | 0x164–0x16C | `sw_pre`, `sw_len`, `sw_maxlag` |
| 8–11 | 0x120–0x12C | `fe_dtrim[0..3]` | | 28–30 | 0x170–0x178 | `mb_pre`, `mb_len`, `mb_maxlag` |
| 12 | 0x130 | `fe_dcblk_en` | | 31–33 | 0x17C–0x184 | `sw_band_lo`, `sw_band_hi`, `sw_floor` |
| 13–18 | 0x134–0x148 | `sw_thr_mult2`, `sw_noise_min`, `sw_holdoff`, `sw_dt_min`, `sw_feat_win`, `sw_sym` | | 34–36 | 0x188–0x190 | `mb_band_lo`, `mb_band_hi`, `mb_floor` |
| 19–23 | 0x14C–0x15C | `mb_thr_mult`, `mb_settle_mult`, `mb_dwell`, `mb_holdoff`, `mb_noise_min` | | 37–40 | 0x194–0x1A0 | `mb_mask_pre`, `mb_mask_post`, `sw_mask_pre`, `sw_mask_post` |
| | | | | 41–42 | 0x1A4–0x1A8 | `sw_tn_max`, `sw_pk_mult` |

Change configuration only while `CTRL.run = 0` (followed by a soft reset), so the detectors never see a half-applied set.

### 3.3 Mailbox RAM (0x8000 – 0xFFFF)

8192 words, read-only: word i at `0x8000 + 4·i`. Messages wrap around the end of the ring.

---

## 4. Mailbox messages

Every message starts with a header word; a message is either written completely or dropped (`MBOX_DROP_CNT`), never partially.

**Header:** `[31:28]` type, `[27:16]` length in words (header included), `[15:0]` event or job id.

| Type | Message | Length | Payload after the header |
|---|---|---|---|
| 1 | SW record | 16 | ts_lo, ts_hi, ref_ch; then per channel c = 0…3: pk_pos (s18, sign-extended), pk_neg (s18), `{neg_off[8:0] << 16 \| pos_off[8:0]}` |
| 2 | MB record | 5 | ts_lo, ts_hi (100 kHz units, group delay removed), idx (20 kHz index), `{reliable << 8 \| ref_ch}` |
| 3 | Job dispatched | 2 | `{type << 31 \| mask_en << 16 \| start}` |
| 4 | Job late (dropped) | 1 | — (id in the header) |
| 5 | Engine result | 17 + 6·(2L+1) | see below |

**Engine result (type 5):**

| Word | Content |
|---|---|
| 1 | `{type << 31 \| err << 30 \| maxlag L << 16}` |
| 2 | Exponents: `{er2 << 20 \| er1 << 15 \| er0 << 10 \| ez1 << 5 \| ez0}` |
| 3–4 | Coarse lags, s8 each: pairs 0–3 in word 3 (pair 0 in [7:0]), pairs 4–5 in word 4 |
| 5–10 | Peak values, pairs 0–5 (s18, sign-extended) |
| 11–16 | Sidelobes, pairs 0–5 (u18) |
| 17 … | Lag window, lag-major: for m = −L … +L, the 6 pair values (s18, sign-extended) |

Lengths: SW result (L = 100): 1223 words; MB result (L = 20): 263 words. The 8192-word ring holds about 6 SW results.

**Ordering:** messages appear in the order the events occur, with one exception. When several sources are waiting at the same moment, they are written in a fixed priority: late job, result, job dispatch, SW record, MB record. So events within the same few microseconds can swap. In all 15 Milestone B scenarios the order matched Milestone B exactly. A result is emitted when the engine finishes; its lag window is staged inside `gsl_mbox` during the peak search.

**Pair order:** (0,1) (0,2) (0,3) (1,2) (1,3) (2,3); lag > 0 means the second channel lags.

---

## 5. Replay DMA

- **Frame format in SDRAM:** 64-bit little-endian words, the same `{ch0, ch1, ch2, ch3}` s16 frames as `mb_*.hex` (ch0 in bits 63:48), stored consecutively from `RP_BASE`.
- **Refill:** the DMA keeps the replay FIFO between 512 and 1000 frames, reading as 32-bit Avalon transfers (two per frame).
- **Non-looping:** after `RP_NFRAMES` frames it stops fetching; when the FIFO drains it stops the frame cadence and sets `replay_done`, so the end of a recording is not an underflow.
- **Looping** (soak test): it wraps to `RP_BASE` and increments `RP_LOOPS`. Timestamps keep counting, so the same shots repeat every pass, offset by the recording length.
- **What a pass reproduces.** Every pass gives the same events with near-identical values, but not necessarily bit-identical ones: the noise blocks, the ÷5 decimation phase and the filter states don't line up across the wrap. Measured on the 200 m recording: identical SW onsets and SW results, MB onsets shifted by up to 4 samples (within one decimated sample), equal coarse lags.
- **Bandwidth:** 800 kB/s, a tiny fraction of the SDRAM's ~200 MB/s.

---

## 6. Verification (C1, simulation)

| Check | Method | Result |
|---|---|---|
| All 15 Milestone B scenarios | Replayed from an SDRAM model (random wait states) through the DMA; a bus-model processor reads every message over Avalon-MM and rebuilds the Milestone B log lines | **Identical** to Milestone B: every record, job, result and lag window, in the same order; x and y streams byte-identical |
| Overflow | Processor 8000 cycles slower per word, 900 rpm burst | 5 messages dropped and counted; the 11 results that arrived are complete and identical |
| Loop mode | 3 passes of the 200 m recording | Same events every pass (see §5) |
| Registers | ID read, configuration written and read back over the bus | Used by every run above |

## 7. Needed from you later (C3, board stage)

1. **The exact Quartus version** (for example 22.1std or 23.1std), and whether the IP Catalog shows **"Nios V/m Microcontroller"** (and **"Nios V/g"**, which has an FPU).
2. **The Terasic DE10-Standard System CD** pin assignments (`DE10_Standard_golden_top` `.qsf`), so the top level uses the board's real pin names (CLOCK_50, DRAM_*, LEDR, KEY, …).
