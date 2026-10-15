# ============================================================
# gsl_mm_hw.tcl — Platform Designer component: GSL detector (gsl_mm)
#
# Add verilog/gsl_mm to Platform Designer's IP search path; the component
# then appears as "GSL detector (gsl_mm)" in the IP Catalog under "GSL".
#
# The coefficient files (.hex) are read with $readmemh, which Quartus
# resolves relative to the PROJECT directory: create_project.tcl copies
# them into verilog/board/.
# ============================================================
package require -exact qsys 16.1

set_module_property NAME                         gsl_mm
set_module_property DISPLAY_NAME                 "GSL detector (gsl_mm)"
set_module_property DESCRIPTION                  "Gunshot locator detection core with registers, result mailbox and replay DMA"
set_module_property VERSION                      1.0
set_module_property GROUP                        "GSL"
set_module_property AUTHOR                       "GSL"
set_module_property INSTANTIATE_IN_SYSTEM_MODULE true
set_module_property EDITABLE                     false
set_module_property ELABORATION_CALLBACK         ""

# ---------------- files ----------------
add_fileset QUARTUS_SYNTH QUARTUS_SYNTH "" ""
set_fileset_property QUARTUS_SYNTH TOP_LEVEL gsl_mm
foreach {name path} {
    gsl_mm.v              gsl_mm.v
    gsl_mbox.v            gsl_mbox.v
    gsl_replay_dma.v      gsl_replay_dma.v
    gsl_core.v            ../stream/gsl_core.v
    replay_src.v          ../stream/replay_src.v
    timebase.v            ../stream/timebase.v
    fe_cond.v             ../stream/fe_cond.v
    mb_decim.v            ../stream/mb_decim.v
    sw_det.v              ../stream/sw_det.v
    mb_det.v              ../stream/mb_det.v
    evt_sched.v           ../stream/evt_sched.v
    cap_ring.v            ../gcc_engine/cap_ring.v
    gcc_loader.v          ../gcc_engine/gcc_loader.v
    peak_search.v         ../gcc_engine/peak_search.v
    gcc_engine.v          ../gcc_engine/gcc_engine.v
    fft_sdp_ram.v         ../fft2048/fft_sdp_ram.v
    fft_buf_bank.v        ../fft2048/fft_buf_bank.v
    fft2048_twiddle_rom.v ../fft2048/fft2048_twiddle_rom.v
    fft2048_bfly.v        ../fft2048/fft2048_bfly.v
    fft2048.v             ../fft2048/fft2048.v
    xspec_phat_unit.v     ../xspec_phat/xspec_phat_unit.v
    xspec_phat.v          ../xspec_phat/xspec_phat.v
} {
    if {$name == "gsl_mm.v"} {
        add_fileset_file $name VERILOG PATH $path TOP_LEVEL_FILE
    } else {
        add_fileset_file $name VERILOG PATH $path
    }
}

# ---------------- parameters ----------------
add_parameter FRAME_DIV INTEGER 1000 "Clock cycles per 4-channel sample (1000 = 100 kHz at 100 MHz)"
set_parameter_property FRAME_DIV HDL_PARAMETER true
set_parameter_property FRAME_DIV ALLOWED_RANGES 64:65535
foreach {p v} {
    FD_INIT    fe_fdelay.hex
    DEC_INIT   decim63.hex
    TW_INIT_RE tw2048_re.hex
    TW_INIT_IM tw2048_im.hex
    RSQRT_INIT rsqrt1024.hex
} {
    add_parameter $p STRING $v "Coefficient file (relative to the Quartus project directory)"
    set_parameter_property $p HDL_PARAMETER true
}

# ---------------- clock and reset ----------------
add_interface clock clock end
set_interface_property clock clockRate 0
add_interface_port clock clk clk Input 1

add_interface reset reset end
set_interface_property reset associatedClock clock
set_interface_property reset synchronousEdges DEASSERT
add_interface_port reset rst_n reset_n Input 1

# ---------------- s0: registers + mailbox (64 KB) ----------------
add_interface s0 avalon end
set_interface_property s0 addressUnits                    WORDS
set_interface_property s0 associatedClock                 clock
set_interface_property s0 associatedReset                 reset
set_interface_property s0 readLatency                     0
set_interface_property s0 maximumPendingReadTransactions  4
set_interface_property s0 readWaitTime                    0
set_interface_property s0 writeWaitTime                   0
set_interface_property s0 bitsPerSymbol                   8
set_interface_property s0 explicitAddressSpan             65536
add_interface_port s0 avs_address       address       Input  14
add_interface_port s0 avs_read          read          Input  1
add_interface_port s0 avs_write         write         Input  1
add_interface_port s0 avs_writedata     writedata     Input  32
add_interface_port s0 avs_readdata      readdata      Output 32
add_interface_port s0 avs_readdatavalid readdatavalid Output 1

# ---------------- m0: replay DMA (reads SDRAM) ----------------
add_interface m0 avalon start
set_interface_property m0 addressUnits    SYMBOLS
set_interface_property m0 associatedClock clock
set_interface_property m0 associatedReset reset
set_interface_property m0 bitsPerSymbol   8
add_interface_port m0 avm_address       address       Output 32
add_interface_port m0 avm_read          read          Output 1
add_interface_port m0 avm_waitrequest   waitrequest   Input  1
add_interface_port m0 avm_readdata      readdata      Input  32
add_interface_port m0 avm_readdatavalid readdatavalid Input  1

# ---------------- interrupt: mailbox not empty ----------------
add_interface irq interrupt end
set_interface_property irq associatedAddressablePoint s0
set_interface_property irq associatedClock            clock
set_interface_property irq associatedReset            reset
add_interface_port irq irq irq Output 1

# ---------------- adc: live ADC frames (Milestone D; tie off for now) ----------------
add_interface adc conduit end
add_interface_port adc adc_valid valid Input 1
add_interface_port adc adc_data  data  Input 64
