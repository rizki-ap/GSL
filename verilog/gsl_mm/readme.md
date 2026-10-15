# gsl_mm — processor interface (Milestone C1)

`gsl_core` wrapped for a processor: Avalon-MM registers, a result mailbox and a replay DMA. Built for Nios V; the interface is processor-neutral.

```
cd verilog/stream && make        # Milestone B first: C1 compares against its outputs
cd ../gsl_mm && make             # ~12 min: build, 15 scenarios + 2 stress tests
                                 # ends: Milestone C1: PASSED
```

- **Status: Milestone C1 passed in simulation** (Verilator 5.020). Not yet on the board.
- **Depends on:** `../stream` (with its scenario files and outputs), `../gcc_engine`, `../fft2048`, `../xspec_phat`, **Verilator**.
- **Interface specification:** `c1_interface.md` (register map, message format, replay DMA, Platform Designer system).
- **Version:** 2026-10-12

## Files

| File | What it is |
|---|---|
| `gsl_mm.v` | Top: Avalon-MM slave (64 KB window, read latency 2), registers, configuration, soft reset, IRQ; instantiates `gsl_core`, `gsl_mbox`, `gsl_replay_dma` |
| `gsl_mbox.v` | Result mailbox: serialises records, jobs, late jobs and results (with lag windows) into an 8192-word ring; whole messages or nothing |
| `gsl_replay_dma.v` | Avalon-MM read master: streams a recording from SDRAM into `replay_src`; loop mode for the soak test |
| `tb_gsl_mm.v` | Testbench: SDRAM model with random wait states; a processor model that uses only the bus |
| `c1_check.py` | Runs every Milestone B scenario through `gsl_mm` and compares with Milestone B; overflow and loop tests |
| `c1_interface.md` | The specification |

`../stream/gsl_core.v` gained three status outputs for this (`sw_noise_valid`, `mb_noise_valid`, `job_mask_en_o`); Milestone B is unchanged.

## Results

| Check | Result |
|---|---|
| 15 Milestone B scenarios | Every record, job, result and lag window identical to Milestone B, **in the same order**; x and y streams byte-identical, so the DMA delivered every frame |
| Overflow (processor 8000 cycles slower per word, 900 rpm burst) | 5 messages dropped and counted; the 11 results that arrived are complete and identical |
| Loop mode (3 passes, 200 m) | Same events every pass; SW onsets and SW results identical; MB onsets within one decimated sample (decimation phase shifts at the wrap) |
| Errors | No late jobs, overflows, underflows or result errors in any run |

## Resources (yosys, indicative)

| Block | LUTs | FFs | M10K |
|---|---|---|---|
| `gsl_mbox` | 1,820 | 1,009 | 35 |
| `gsl_replay_dma` | 956 | 267 | 0 |
| Registers in `gsl_mm` | ≈ 1–2k (estimate) | — | 0 |

With the area-optimised `gsl_core` that is about 24k LUTs and 168 M10K before the Nios V system (Nios V/m plus 128 KB code RAM adds about 105 M10K: ~50% of the device).

## Notes

- **Messages and ordering:** see `c1_interface.md` §4. When several sources wait at the same moment a fixed priority applies, so events within microseconds of each other could swap; it never happened in the 15 scenarios.
- **Loop passes are near-identical, not bit-identical** (§5 of the spec), so the soak test compares events per pass with small tolerances.
- **Configuration** should change only while stopped, followed by a soft reset (`CTRL` bit 2).
