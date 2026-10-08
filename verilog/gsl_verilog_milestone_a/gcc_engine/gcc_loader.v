`timescale 1ns/1ps

// ============================================================
//  gcc_loader.v
//  Copy one 4-channel window from a capture ring into the two FFT
//  input buffers, packed and bit-reversed, zero-padded to 2048.
//
//  Spec: verilog/fpga_blocks.md section 3.2
//
//    Z0[bitrev(n)] = ch0[n] + j*ch1[n]
//    Z1[bitrev(n)] = ch2[n] + j*ch3[n]      n = 0 .. 2047
//    ch_i[n] = ring[(ld_base + n) mod 4096]  for n < ld_len, else 0
//
//  All 4 channels come from the SAME ring address (shared-reference
//  windowing, project.md limitation 3): one ring word holds all 4.
//  Z0 and Z1 are different buffers, so both are written every cycle.
//
//  Input range: ring samples come from fe_cond with |x| < 2^15,
//  inside the fft2048 input contract (|x| < 2^16).
//
//  Cycles: 2048 + RL + 2 = 2,052   (20.5 us @ 100 MHz)
// ============================================================
module gcc_loader #(
    parameter N     = 2048,
    parameter LOG2N = 11,
    parameter RAW   = 12,            // ring address width
    parameter RL    = 2              // ring read latency
)(
    input  wire              clk,
    input  wire              rst_n,

    input  wire              ld_start,
    input  wire [RAW-1:0]    ld_base,      // ring address of window sample 0
    input  wire [LOG2N:0]    ld_len,       // window length, 1 .. 2048
    output reg               ld_busy,
    output reg               ld_done,

    output wire [RAW-1:0]    ring_raddr,
    input  wire [71:0]       ring_rdata,

    output wire              z_we,
    output wire [LOG2N-1:0]  z_waddr,      // bit-reversed
    output wire [35:0]       z0_wdata,
    output wire [35:0]       z1_wdata
);
    function [LOG2N-1:0] bitrev;
        input [LOG2N-1:0] v;
        integer i;
        for (i = 0; i < LOG2N; i = i + 1)
            bitrev[i] = v[LOG2N-1-i];
    endfunction

    reg [LOG2N:0]   n;                      // 0 .. N
    reg             issuing;
    wire            iss_v   = issuing;
    wire            in_win  = (n < ld_len_r);
    reg [LOG2N:0]   ld_len_r;
    reg [RAW-1:0]   base_r;

    assign ring_raddr = base_r + n[RAW-1:0];

    // delay line for (valid, in-window, n) to match ring latency
    reg [RL-1:0]        v_sr, w_sr;
    reg [LOG2N-1:0]     n_sr [0:RL-1];
    integer k;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            v_sr <= {RL{1'b0}};
        end else begin
            v_sr <= {v_sr[RL-2:0], iss_v};
        end
    end

    always @(posedge clk) begin
        w_sr    <= {w_sr[RL-2:0], in_win};
        n_sr[0] <= n[LOG2N-1:0];
        for (k = 1; k < RL; k = k + 1)
            n_sr[k] <= n_sr[k-1];
    end

    assign z_we     = v_sr[RL-1];
    assign z_waddr  = bitrev(n_sr[RL-1]);
    assign z0_wdata = w_sr[RL-1] ? ring_rdata[71:36] : 36'd0;
    assign z1_wdata = w_sr[RL-1] ? ring_rdata[35:0]  : 36'd0;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            n        <= {(LOG2N+1){1'b0}};
            issuing  <= 1'b0;
            ld_busy  <= 1'b0;
            ld_done  <= 1'b0;
            ld_len_r <= {(LOG2N+1){1'b0}};
            base_r   <= {RAW{1'b0}};
        end else begin
            ld_done <= 1'b0;
            if (ld_start && !ld_busy) begin
                base_r   <= ld_base;
                ld_len_r <= ld_len;
                n        <= {(LOG2N+1){1'b0}};
                issuing  <= 1'b1;
                ld_busy  <= 1'b1;
            end else if (issuing) begin
                if (n == N - 1) issuing <= 1'b0;
                n <= n + 1'b1;
            end else if (ld_busy && v_sr == {RL{1'b0}}) begin
                ld_busy <= 1'b0;
                ld_done <= 1'b1;
            end
        end
    end
endmodule
