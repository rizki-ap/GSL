`timescale 1ns/1ps

// ============================================================
//  replay_src.v
//  Sample source for the stream path: replays 4-channel frames written
//  by the host (testbench now, HPS / Nios later), or passes the live
//  ADC through. Also generates the frame cadence.
//
//  Spec: verilog/fpga_blocks.md section 2.2
//
//  Frame format (64 bit): {ch0[63:48], ch1[47:32], ch2[31:16], ch3[15:0]}, s16
//
//  Cadence: frame_tick every FRAME_DIV clocks while run = 1
//           (1000 @ 100 MHz = 100 kHz). In Milestone D, adc_if takes over
//           the cadence; the replay path then follows adc's frame_tick.
//  On each tick (src_sel = 1) one frame is popped; if the FIFO is empty
//  no sample is produced and rp_underflow (sticky) is set -- the
//  timeline would be broken, so treat it as a test failure.
// ============================================================
module replay_src #(
    parameter FRAME_DIV = 1000,
    parameter FIFO_AW   = 10                 // 1024 frames = 10 ms
)(
    input  wire              clk,
    input  wire              rst_n,

    input  wire              run,
    input  wire              src_sel,        // 1 = replay, 0 = ADC

    input  wire              rp_wr_en,
    input  wire [63:0]       rp_wr_data,
    output wire [FIFO_AW:0]  rp_level,
    output wire              rp_full,

    input  wire              adc_valid,
    input  wire [63:0]       adc_data,

    output reg               frame_tick,
    output reg               src_valid,
    output reg  [63:0]       src_data,
    output reg               rp_underflow
);
    // ---------------- frame divider ----------------
    reg [15:0] div;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            div        <= 16'd0;
            frame_tick <= 1'b0;
        end else begin
            frame_tick <= 1'b0;
            if (run) begin
                if (div == FRAME_DIV - 1) begin
                    div        <= 16'd0;
                    frame_tick <= 1'b1;
                end else
                    div <= div + 1'b1;
            end else
                div <= 16'd0;
        end
    end

    // ---------------- FIFO ----------------
    reg [63:0]      mem [0:(1<<FIFO_AW)-1];
    reg [FIFO_AW:0] wptr, rptr;
    assign rp_level = wptr - rptr;
    assign rp_full  = (rp_level == (1 << FIFO_AW));
    wire   empty    = (wptr == rptr);
    wire   pop      = frame_tick && src_sel && !empty;

    reg [63:0] q;
    always @(posedge clk) begin
        if (rp_wr_en && !rp_full) mem[wptr[FIFO_AW-1:0]] <= rp_wr_data;
        q <= mem[rptr[FIFO_AW-1:0]];
    end

    reg pop_d;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            wptr <= 0; rptr <= 0; pop_d <= 1'b0;
            src_valid <= 1'b0; src_data <= 64'd0; rp_underflow <= 1'b0;
        end else begin
            if (rp_wr_en && !rp_full) wptr <= wptr + 1'b1;
            if (pop) rptr <= rptr + 1'b1;
            pop_d <= pop;
            if (frame_tick && src_sel && empty) rp_underflow <= 1'b1;

            if (src_sel) begin
                src_valid <= pop_d;                 // data of the popped address
                if (pop_d) src_data <= q;
            end else begin
                src_valid <= adc_valid;
                if (adc_valid) src_data <= adc_data;
            end
        end
    end
endmodule
