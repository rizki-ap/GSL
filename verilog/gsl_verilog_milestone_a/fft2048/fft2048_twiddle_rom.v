// ============================================================
//  fft2048_twiddle_rom.v
//  N/2 = 1024-entry twiddle ROM, W_N^k = cos - j*sin, Q1.17.
//  Registered output, 1-cycle latency. Init files from
//  gen_twiddle2048.py (Quartus infers M10K ROM from $readmemh).
// ============================================================
`timescale 1ns/1ps

module fft2048_twiddle_rom #(
    parameter TW      = 18,
    parameter AW      = 10,
    parameter INIT_RE = "tw2048_re.hex",
    parameter INIT_IM = "tw2048_im.hex"
)(
    input  wire                 clk,
    input  wire [AW-1:0]        addr,
    output reg  signed [TW-1:0] w_re,
    output reg  signed [TW-1:0] w_im
);
    (* romstyle = "M10K" *) reg [TW-1:0] rom_re [0:(1<<AW)-1];
    (* romstyle = "M10K" *) reg [TW-1:0] rom_im [0:(1<<AW)-1];

    initial begin
        $readmemh(INIT_RE, rom_re);
        $readmemh(INIT_IM, rom_im);
    end

    always @(posedge clk) begin
        w_re <= rom_re[addr];
        w_im <= rom_im[addr];
    end
endmodule
