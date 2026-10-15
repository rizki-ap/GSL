// ============================================================
// gsl_top — DE10-Standard top level (Milestone C3)
//
//  Pin names follow Terasic's DE10_Standard_default.qsf. The Platform
//  Designer system gsl_sys (see readme.md) holds the PLL, Nios V/g,
//  on-chip RAM, JTAG UART, LED PIO, SDRAM controller, JTAG-to-Avalon
//  master and the detector (gsl_mm).
//
//  KEY[0]  reset (press to reset the system)
//  LEDR[9] heartbeat from CLOCK_50: the FPGA is configured
//  LEDR[7:0] driven by the Nios V program
// ============================================================
module gsl_top (
    input  wire        CLOCK_50,
    input  wire [3:0]  KEY,
    input  wire [9:0]  SW,
    output wire [9:0]  LEDR,
    output wire [6:0]  HEX0, HEX1, HEX2, HEX3, HEX4, HEX5,

    output wire [12:0] DRAM_ADDR,
    output wire [1:0]  DRAM_BA,
    output wire        DRAM_CAS_N,
    output wire        DRAM_CKE,
    output wire        DRAM_CLK,
    output wire        DRAM_CS_N,
    inout  wire [15:0] DRAM_DQ,
    output wire        DRAM_LDQM,
    output wire        DRAM_RAS_N,
    output wire        DRAM_UDQM,
    output wire        DRAM_WE_N
);
    // ---- heartbeat: ~1.5 Hz from the 50 MHz board clock ----
    reg [24:0] hb = 25'd0;
    always @(posedge CLOCK_50) hb <= hb + 1'b1;
    assign LEDR[9] = hb[24];
    assign LEDR[8] = 1'b0;

    // seven-segment displays off (active low)
    assign HEX0 = 7'h7F; assign HEX1 = 7'h7F; assign HEX2 = 7'h7F;
    assign HEX3 = 7'h7F; assign HEX4 = 7'h7F; assign HEX5 = 7'h7F;

    gsl_sys u_sys (
        .clk_clk        (CLOCK_50),
        .reset_reset_n  (KEY[0]),
        .sdram_clk_clk  (DRAM_CLK),
        .sdram_addr     (DRAM_ADDR),
        .sdram_ba       (DRAM_BA),
        .sdram_cas_n    (DRAM_CAS_N),
        .sdram_cke      (DRAM_CKE),
        .sdram_cs_n     (DRAM_CS_N),
        .sdram_dq       (DRAM_DQ),
        .sdram_dqm      ({DRAM_UDQM, DRAM_LDQM}),
        .sdram_ras_n    (DRAM_RAS_N),
        .sdram_we_n     (DRAM_WE_N),
        .led_export     (LEDR[7:0]),
        .adc_valid      (1'b0),
        .adc_data       (64'd0)
    );
endmodule
