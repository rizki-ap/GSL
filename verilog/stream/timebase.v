`timescale 1ns/1ps

// ============================================================
//  timebase.v
//  64-bit sample index. Frame k (the k-th src_valid after reset or
//  ts_clear) carries index k: src_ts is valid together with src_valid,
//  and ts64 increments right after it.
//
//  Spec: verilog/fpga_blocks.md section 2.3  (PPS latch: Milestone D)
// ============================================================
module timebase (
    input  wire        clk,
    input  wire        rst_n,
    input  wire        src_valid,
    input  wire        ts_clear,
    output wire [63:0] src_ts,
    output reg  [63:0] ts64
);
    assign src_ts = ts64;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n)          ts64 <= 64'd0;
        else if (ts_clear)   ts64 <= 64'd0;
        else if (src_valid)  ts64 <= ts64 + 64'd1;
    end
endmodule
