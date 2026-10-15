# board — DE10-Standard bring-up (Milestone C3)

Goal of this first board step: **replay the Milestone B scenarios on the real FPGA and get results byte-identical to the simulation.** Your PC loads a recording into SDRAM over JTAG, a small Nios V program replays it and logs every detector message into SDRAM, your PC reads the log back and `board_check.py` compares it with the simulation.

- **Tools:** Quartus Prime Lite 24.1 (with Nios V), the DE10-Standard on USB-Blaster, Python 3.
- **Status:** written and checked as far as possible without Quartus. The top level, C program and PC scripts were compile-checked or tested against the C1 simulation; the Platform Designer steps and Quartus scripts are untested until you run them. Expect a round or two of fixes: send me the exact messages.
- **Version:** 2026-10-12

## Files

| File | Purpose |
|---|---|
| `create_project.tcl` | Creates the Quartus project `gsl_board` (device, files, 106 pins from Terasic's `.qsf`) and copies the coefficient files |
| `gsl_top.v` | Top level with the board's pin names |
| `gsl_top.sdc` | Timing constraints (100 MHz system, SDRAM I/O) |
| `software/main.c` | Nios V program: ID check, SDRAM test, replay and logging |
| `gsl_sc.tcl` | System Console commands: load a recording, save the log |
| `board_check.py` | Converts recordings for SDRAM; compares board logs with the simulation |
| `../gsl_mm/gsl_mm_hw.tcl` | Makes the detector a Platform Designer component |

You create `gsl_sys.qsys` in step 2; **commit it**, it is part of the design.

---

## 1. Create the project

```bash
cd verilog/board
quartus_sh -t create_project.tcl
```

Then open `gsl_board.qpf` in Quartus.

## 2. Build the system in Platform Designer

**Tools > Platform Designer**, then **File > New System**, and save it as `gsl_sys.qsys` in `verilog/board`.

**Make the detector visible:** **Tools > Options > IP Search Locations**, add the folder `verilog/gsl_mm` (full path), then **File > Refresh System**. "GSL detector (gsl_mm)" now appears in the IP Catalog under **GSL**.

### 2a. Components

Use exactly these instance names (right-click > Rename); the C program depends on them.

| Instance | IP Catalog name | Settings |
|---|---|---|
| `clk` | (already there: Clock Source) | 50 MHz |
| `pll` | PLL Intel FPGA IP | Reference clock 50 MHz; **2 clocks**: outclk0 100 MHz phase 0; outclk1 100 MHz, phase shift **7000 ps** (= −3 ns); locked output off |
| `cpu` | Nios V/g Processor | **Floating-point unit on**. Reset agent `ram.s1`, offset 0. **Peripheral region A: base 0x0800_0000, size 128 MB** (all I/O and SDRAM uncached) |
| `ram` | On-Chip Memory (RAM or ROM) Intel FPGA IP | RAM, 32 bits, **131072 bytes**, "Initialize memory content" **off** |
| `jtag_uart` | JTAG UART Intel FPGA IP | defaults |
| `sysid` | System ID Peripheral Intel FPGA IP | ID 0x20261012 |
| `led` | PIO (Parallel I/O) Intel FPGA IP | 8 bits, output |
| `sdram` | SDRAM Controller Intel FPGA IP | data width 16; chip selects 1, banks 4; row 13, column 10; CAS latency 3; 2 initialization refreshes; refresh every 7.8125 µs; power-up delay 100 µs; t_rfc 70 ns, t_rp 15 ns, t_rcd 15 ns, t_ac 5.4 ns, t_wr 14 ns |
| `jtag_master` | JTAG to Avalon Master Bridge Intel FPGA IP | defaults |
| `gsl` | GSL detector (gsl_mm) | FRAME_DIV **1000** (default) |

### 2b. Connections

| From | To |
|---|---|
| `clk.clk` | `pll.refclk` |
| `clk.clk_reset` | `pll.reset` |
| `pll.outclk0` | the clock input of **every** other component |
| `cpu.instruction_manager` | `ram.s1`, `cpu.dm_agent` |
| `cpu.data_manager` | `ram.s1`, `cpu.dm_agent`, `cpu.timer_sw_agent`, `jtag_uart.avalon_jtag_slave`, `sysid.control_slave`, `led.s1`, `sdram.s1`, `gsl.s0` |
| `jtag_master.master` | `ram.s1`, `sdram.s1`, `gsl.s0` |
| `gsl.m0` | `sdram.s1` |
| Interrupts (`cpu.platform_irq_rx`) | `jtag_uart.irq` = 0, `gsl.irq` = 1 |
| Resets | **System > Create Global Reset Network** |

### 2c. Exports (Export column, double-click)

| Interface | Export name |
|---|---|
| `clk.clk_in` | `clk` |
| `clk.clk_in_reset` | `reset` |
| `pll.outclk1` | `sdram_clk` |
| `sdram.wire` | `sdram` |
| `led.external_connection` | `led` |
| `gsl.adc` | `adc` |

### 2d. Base addresses (then lock them)

| Interface | Base | Size |
|---|---|---|
| `ram.s1` | 0x0000_0000 | 128 KB |
| `gsl.s0` | 0x0800_0000 | 64 KB |
| `jtag_uart.avalon_jtag_slave` | 0x0801_0000 | |
| `sysid.control_slave` | 0x0801_0040 | |
| `led.s1` | 0x0801_0050 | |
| `cpu.timer_sw_agent` | 0x0801_0080 | |
| `cpu.dm_agent` | 0x0802_0000 | |
| `sdram.s1` | 0x0C00_0000 | 64 MB |

The Messages panel must show **no errors** (warnings are fine). **Generate > Generate HDL**: Verilog, no simulation model. Close Platform Designer.

## 3. Compile

**Processing > Start Compilation** (about 20–40 minutes).

**Send me:**
1. **Fitter > Resource Section > Resource Usage Summary:** ALMs, M10K blocks, DSP blocks.
2. **Timing Analyzer > Slow 1100mV 85C Model > Fmax Summary**, plus the setup slack per clock (is anything red?).
3. Every **critical warning**, and whether this one appears: `gsl_top.sdc: SDRAM clock pin not found`.

## 4. Program the FPGA

**Tools > Programmer**, file `output_files/gsl_board.sof`, **Start**. `LEDR[9]` starts blinking.

## 5. Build and start the Nios V program

In the **Nios V Command Shell** (Windows Start menu; on Linux `niosv-shell` in the Quartus installation):

```bash
cd verilog/board
niosv-bsp -c -t=hal --sopcinfo=gsl_sys.sopcinfo software/bsp/settings.bsp
niosv-app --bsp-dir=software/bsp --app-dir=software/app --srcs=software/main.c
cmake -S software/app -B software/app/build -G "Unix Makefiles" -DCMAKE_BUILD_TYPE=Release
make -C software/app/build
niosv-download -g software/app/build/*.elf
juart-terminal
```

Expected in the terminal:

```
GSL board bring-up (Milestone C3)
gsl_mm ID 0x47534c31 (ok), VERSION 0x20261012
SDRAM test (1 MB at offset 0x03f00000): ok (0 errors)
Waiting for a recording (System Console: gsl_load <file.bin>) ...
```

`LEDR[0]` blinks while it waits. Leave the terminal open.

## 6. Replay the scenarios

Recordings for SDRAM (once):

```bash
python3 board_check.py prepare        # needs ../stream prepared: make prepare in ../stream
```

**Tools > System Debugging Tools > System Console**; in its Tcl console:

```tcl
cd <full path to>/verilog/board
source gsl_sc.tcl
gsl_info                 ;# gsl_mm ID 0x47534C31 (ok) ...
gsl_run s1               ;# load, replay in real time, save log/s1.bin
```

The Nios V terminal shows the replay and a summary. Then:

```bash
python3 board_check.py check s1       # ends: s1: PASS - identical to simulation
```

All 15 scenarios: in System Console `foreach n {s1 s2 s3 r50 L1500_n05 L1500_n20 L1500_n50 d556_1000 d762_1000 d127_1500 m_tol m_sw m_mb m_hot burst} { gsl_run $n }`, then `for n in s1 s2 s3 r50 L1500_n05 L1500_n20 L1500_n50 d556_1000 d762_1000 d127_1500 m_tol m_sw m_mb m_hot burst; do python3 board_check.py check $n; done`. The checks need the Milestone B outputs (`make` in `../stream`).

## 7. What to send me

- Step 3: resources, Fmax and slack, critical warnings.
- Step 5: the terminal output.
- Step 6: the `board_check.py` output for `s1` (and for all 15 when that works).
- Any error, **copied exactly**.

## Troubleshooting

| Symptom | Likely cause |
|---|---|
| `GSL_BASE` (or `SDRAM_BASE`, `LED_BASE`) undeclared when building the program | An instance name differs from the table in 2a |
| ID reads WRONG | `gsl.s0` not at 0x0800_0000, or `gsl` not connected to `cpu.data_manager` |
| SDRAM test FAILED | SDRAM clock phase: try 6000 ps or 8000 ps for `pll.outclk1`, regenerate, recompile |
| `Waiting for a recording` never ends after `gsl_load` | Data cache: check the peripheral region (2a); `gsl_info` in System Console still works |
| System Console: "no JTAG master found" | Board not programmed, or `jtag_master` missing |
| `board_check.py`: frames n/m not equal, or underflow | Replay DMA could not keep up: send me the log summary |
