// ============================================================
//  fft2048_bfly.v
//  Pipelined radix-2 DIT butterfly with block-floating-point
//  output scaling and convergent rounding.
//
//      P = (A + W*B) / 2^sh
//      Q = (A - W*B) / 2^sh          sh = 0, 1 or 2
//
//  Full precision is kept until the single rounding step:
//      A*2^17 +/- W*B   (38-bit)  -> round(>> 17+sh) -> sat s18
//
//  bypass = 1 means W = 1 (k = 0), which Q1.17 cannot represent:
//  the multiplier is skipped and B is used directly (exact).
//
//  Latency: 2 cycles (M = multiply, A = add/round/saturate).
//  Cyclone V mapping: each of W*B re / im is an 18x18
//  "sum of two" -> 2 variable-precision DSP blocks per butterfly.
// ============================================================
`timescale 1ns/1ps

module fft2048_bfly #(
    parameter DW = 18,              // data width (Q1.17)
    parameter TW = 18               // twiddle width (Q1.17)
)(
    input  wire                 clk,
    input  wire signed [DW-1:0] a_re,
    input  wire signed [DW-1:0] a_im,
    input  wire signed [DW-1:0] b_re,
    input  wire signed [DW-1:0] b_im,
    input  wire signed [TW-1:0] w_re,
    input  wire signed [TW-1:0] w_im,
    input  wire                 bypass,
    input  wire [1:0]           sh,
    output reg  signed [DW-1:0] p_re,
    output reg  signed [DW-1:0] p_im,
    output reg  signed [DW-1:0] q_re,
    output reg  signed [DW-1:0] q_im,
    output reg                  sat
);
    localparam integer FRAC = TW - 1;          // 17
    localparam integer PW   = DW + TW + 1;     // 37: W*B width
    localparam integer SW   = PW + 1;          // 38: A*2^17 +/- W*B
    localparam integer MAXV =  (1 << (DW-1)) - 1;
    localparam integer MINV = -(1 << (DW-1));

    // ---------------- M stage: complex multiply ----------------
    reg signed [PW-1:0] m_re, m_im;
    reg signed [DW-1:0] a_re_d, a_im_d;
    reg [1:0]           sh_d;

    always @(posedge clk) begin
        if (bypass) begin
            m_re <= $signed({b_re, {FRAC{1'b0}}});
            m_im <= $signed({b_im, {FRAC{1'b0}}});
        end else begin
            m_re <= b_re * w_re - b_im * w_im;
            m_im <= b_re * w_im + b_im * w_re;
        end
        a_re_d <= a_re;
        a_im_d <= a_im;
        sh_d   <= sh;
    end

    // ---------------- A stage: add, round, saturate ------------
    wire signed [SW-1:0] a_re_x = $signed({a_re_d, {FRAC{1'b0}}});
    wire signed [SW-1:0] a_im_x = $signed({a_im_d, {FRAC{1'b0}}});

    wire signed [SW-1:0] s_re = a_re_x + m_re;
    wire signed [SW-1:0] s_im = a_im_x + m_im;
    wire signed [SW-1:0] d_re = a_re_x - m_re;
    wire signed [SW-1:0] d_im = a_im_x - m_im;

    // Convergent (round-half-to-even) right shift by s, then
    // saturate to DW bits. Returns {sat, value}.
    function [DW:0] rnd_sat;
        input signed [SW-1:0] x;
        input integer         s;
        reg   signed [SW:0]   t;
        begin
            t = x;
            t = t + ((1 <<< (s-1)) - 1) + ((x >>> s) & 1);
            t = t >>> s;
            if (t > MAXV)      rnd_sat = {1'b1, MAXV[DW-1:0]};
            else if (t < MINV) rnd_sat = {1'b1, MINV[DW-1:0]};
            else               rnd_sat = {1'b0, t[DW-1:0]};
        end
    endfunction

    wire [DW:0] r_pr = rnd_sat(s_re, FRAC + sh_d);
    wire [DW:0] r_pi = rnd_sat(s_im, FRAC + sh_d);
    wire [DW:0] r_qr = rnd_sat(d_re, FRAC + sh_d);
    wire [DW:0] r_qi = rnd_sat(d_im, FRAC + sh_d);

    always @(posedge clk) begin
        p_re <= r_pr[DW-1:0];
        p_im <= r_pi[DW-1:0];
        q_re <= r_qr[DW-1:0];
        q_im <= r_qi[DW-1:0];
        sat  <= r_pr[DW] | r_pi[DW] | r_qr[DW] | r_qi[DW];
    end
endmodule
