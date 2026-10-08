// ============================================================
//  fft_sdp_ram.v
//  Simple dual-port RAM: 1 write port + 1 read port.
//  Registered read, 1-cycle latency. Infers Cyclone V M10K.
//
//  Read-during-write to the SAME address is "don't care":
//  the FFT schedule never does it (see fft2048.v).
// ============================================================
`timescale 1ns/1ps

module fft_sdp_ram #(
    parameter DW = 36,
    parameter AW = 10
)(
    input  wire          clk,
    input  wire          we,
    input  wire [AW-1:0] waddr,
    input  wire [DW-1:0] wdata,
    input  wire [AW-1:0] raddr,
    output reg  [DW-1:0] rdata
);
    (* ramstyle = "M10K, no_rw_check" *)
    reg [DW-1:0] mem [0:(1<<AW)-1];

    always @(posedge clk) begin
        if (we) mem[waddr] <= wdata;
        rdata <= mem[raddr];
    end
endmodule
