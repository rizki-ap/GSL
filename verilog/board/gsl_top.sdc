# ============================================================
# gsl_top.sdc — timing constraints, DE10-Standard (Milestone C3)
# ============================================================

# 50 MHz board clock; the PLL derives 100 MHz system and SDRAM clocks
create_clock -name CLOCK_50 -period 20.000 [get_ports CLOCK_50]
derive_pll_clocks
derive_clock_uncertainty

# JTAG (Nios V debug, JTAG UART, JTAG-to-Avalon master)
create_clock -name altera_reserved_tck -period 41.667 [get_ports -nowarn altera_reserved_tck]
set_clock_groups -asynchronous -group [get_clocks -nowarn altera_reserved_tck]

# Asynchronous board I/O
set_false_path -from [get_ports {KEY[*] SW[*]}]
set_false_path -to   [get_ports {LEDR[*] HEX0[*] HEX1[*] HEX2[*] HEX3[*] HEX4[*] HEX5[*]}]

# ------------------------------------------------------------
# SDRAM (IS42S16320, 100 MHz, CAS latency 3)
# DRAM_CLK comes from PLL output 1 (shifted -3 ns). Its internal pin name
# exists only after Platform Designer generation, so it is looked up by
# pattern; if the pattern does not match, the SDRAM I/O stays unconstrained
# and a warning says so (send the name from the PLL report and it gets fixed).
# ------------------------------------------------------------
set dram_src [get_pins -nowarn -compatibility_mode {*general[1]*|divclk}]
if {[get_collection_size $dram_src] == 1} {
    create_generated_clock -name dram_clk -source $dram_src [get_ports DRAM_CLK]
    set dram_out [get_ports {DRAM_ADDR[*] DRAM_BA[*] DRAM_CAS_N DRAM_CKE DRAM_CS_N DRAM_DQ[*] DRAM_LDQM DRAM_UDQM DRAM_RAS_N DRAM_WE_N}]
    # SDRAM data valid tAC = 5.4 ns after its clock, held tOH = 2.7 ns
    set_input_delay  -clock dram_clk -max 5.4 [get_ports {DRAM_DQ[*]}]
    set_input_delay  -clock dram_clk -min 2.7 [get_ports {DRAM_DQ[*]}]
    # SDRAM setup 1.5 ns, hold 0.8 ns
    set_output_delay -clock dram_clk -max 1.5  $dram_out
    set_output_delay -clock dram_clk -min -0.8 $dram_out
} else {
    post_message -type warning "gsl_top.sdc: SDRAM clock pin not found ([get_collection_size $dram_src] matches); SDRAM I/O is unconstrained"
}
