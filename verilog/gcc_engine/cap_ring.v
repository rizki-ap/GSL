`timescale 1ns/1ps

// ============================================================
//  cap_ring.v
//  Capture ring: 4096 x 72 bit = one 4-channel s18 frame per word,
//  {ch0[71:54], ch1[53:36], ch2[35:18], ch3[17:0]}.
//  Used twice: sw_ring (100 kHz, 41 ms) and mb_ring (20 kHz, 205 ms).
//
//  Write port: stream path (fe_cond / mb_decim), address = sample
//  index mod 4096, never stalls. Read port: gcc_loader.
//  Read latency 2 cycles (M10K with output register).
// ============================================================
module cap_ring #(
    parameter DW = 72,
    parameter AW = 12
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
    reg [DW-1:0] q;

    always @(posedge clk) begin
        if (we) mem[waddr] <= wdata;
        q     <= mem[raddr];
        rdata <= q;
    end
endmodule
