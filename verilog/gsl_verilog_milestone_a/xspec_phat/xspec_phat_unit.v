`timescale 1ns/1ps

// ============================================================
//  xspec_phat_unit.v
//  One cross-spectrum + PHAT normalisation per cycle (pipelined).
//
//      G   = Xb * conj(Xa)                 (matches gs_det_tdoa.py:
//                                           G = X2 * conj(X1))
//      out = 32767 * G / |G|               (0 if G == 0 or below floor)
//
//  Floor: G is zeroed when  msb(|G|) + in_esum < cfg_floor, where
//  in_esum = BFP exponent of Xa's spectrum + that of Xb's spectrum.
//  msb + esum is floor(log2 |G|) in Q1.17 units of the unscaled
//  spectra. cfg_floor = 0 disables it (matches gs_det_tdoa.py).
//  Use it to stop PHAT from normalising rounding residue (e.g. a
//  dead microphone) up to full magnitude.
//
//  Pipeline (latency 7):
//    S1  complex multiply                  37-bit G
//    S2  one's-complement magnitude OR, leading-one position
//    S3  normalise G so its largest component sits in [2^16, 2^17)
//    S4  q = |Gn|^2                        in [2^32, 2^35]
//    S5  q = m * 2^(2e), e in {16,17}; ROM y = 1/sqrt(m)  (rsqrt1024.hex)
//    S6  Gn * y
//    S7  round (convergent) by e + 2, clamp to +/-32767
//
//  Phase comes only from Gn (both components scaled by the same
//  real factor), so the ~0.1 % ROM error affects magnitude only.
//  Output magnitude is ~32767 (0.25 FS), so a packed IFFT input
//  W = G_a + j*G_b stays below 0.5 FS (fft2048 input contract).
//
//  DSP use: S1 2 x "sum of two 18x18", S4 1 x "sum of two",
//           S6 2 x 18x18  -> 4 Cyclone V DSP blocks.
// ============================================================
module xspec_phat_unit #(
    parameter DW       = 18,
    parameter TAGW     = 14,
    parameter ROM_INIT = "rsqrt1024.hex"
)(
    input  wire                 clk,
    input  wire                 rst_n,

    input  wire [6:0]           cfg_floor,

    input  wire                 in_v,
    input  wire [TAGW-1:0]      in_tag,
    input  wire [5:0]           in_esum,
    input  wire signed [DW-1:0] xa_re,
    input  wire signed [DW-1:0] xa_im,
    input  wire signed [DW-1:0] xb_re,
    input  wire signed [DW-1:0] xb_im,

    output wire                 out_v,
    output wire [TAGW-1:0]      out_tag,
    output reg  signed [DW-1:0] g_re,
    output reg  signed [DW-1:0] g_im,
    output wire                 busy
);
    localparam integer GW   = 2*DW + 1;        // 37: G width
    localparam integer MW   = 2*DW;            // 36: magnitude width
    localparam integer LAT  = 7;
    localparam integer OMAX = 32767;

    // ---------------- valid / tag pipeline ----------------
    reg [LAT-1:0]  v_sr;
    reg [TAGW-1:0] tag_sr [0:LAT-1];
    integer j;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) v_sr <= {LAT{1'b0}};
        else        v_sr <= {v_sr[LAT-2:0], in_v};
    end

    always @(posedge clk) begin
        tag_sr[0] <= in_tag;
        for (j = 1; j < LAT; j = j + 1)
            tag_sr[j] <= tag_sr[j-1];
    end

    assign out_v   = v_sr[LAT-1];
    assign out_tag = tag_sr[LAT-1];
    assign busy    = |v_sr;

    // ---------------- S1: G = Xb * conj(Xa) ----------------
    reg signed [GW-1:0] s1_re, s1_im;
    reg [5:0]           s1_esum;
    always @(posedge clk) begin
        s1_re   <= xb_re * xa_re + xb_im * xa_im;
        s1_im   <= xb_im * xa_re - xb_re * xa_im;
        s1_esum <= in_esum;
    end

    // ---------------- S2: magnitude, leading one ----------------
    function [MW-1:0] mag1;                     // one's-complement magnitude
        input [GW-1:0] v;
        mag1 = v[GW-1] ? ~v[MW-1:0] : v[MW-1:0];
    endfunction

    function [5:0] msb_pos;
        input [MW-1:0] v;
        integer i;
        begin
            msb_pos = 6'd0;
            for (i = 0; i < MW; i = i + 1)
                if (v[i]) msb_pos = i;
        end
    endfunction

    wire [MW-1:0] s1_mor = mag1(s1_re) | mag1(s1_im);
    wire [5:0]    s1_msb = msb_pos(s1_mor);
    wire [6:0]    s1_lvl = s1_msb + s1_esum;          // floor(log2 |G|)

    reg signed [GW-1:0] s2_re, s2_im;
    reg [5:0]           s2_msb;
    reg                 s2_zero;
    always @(posedge clk) begin
        s2_re   <= s1_re;
        s2_im   <= s1_im;
        s2_msb  <= s1_msb;
        s2_zero <= ((s1_re == 0) && (s1_im == 0)) || (s1_lvl < cfg_floor);
    end

    // ---------------- S3: normalise to [2^16, 2^17) ----------------
    wire signed [GW-1:0] s2_nre = (s2_msb >= 6'd16) ? (s2_re >>> (s2_msb - 6'd16))
                                                    : (s2_re <<< (6'd16 - s2_msb));
    wire signed [GW-1:0] s2_nim = (s2_msb >= 6'd16) ? (s2_im >>> (s2_msb - 6'd16))
                                                    : (s2_im <<< (6'd16 - s2_msb));

    reg signed [DW-1:0] s3_re, s3_im;
    reg                 s3_zero;
    always @(posedge clk) begin
        s3_re   <= s2_nre[DW-1:0];
        s3_im   <= s2_nim[DW-1:0];
        s3_zero <= s2_zero;
    end

    // ---------------- S4: q = |Gn|^2 ----------------
    reg [MW-1:0]        s4_q;
    reg signed [DW-1:0] s4_re, s4_im;
    reg                 s4_zero;
    always @(posedge clk) begin
        s4_q    <= s3_re * s3_re + s3_im * s3_im;
        s4_re   <= s3_re;
        s4_im   <= s3_im;
        s4_zero <= s3_zero;
    end

    // ---------------- S5: rsqrt ROM ----------------
    // q in [2^32, 2^35]; clamp 2^35 (both components = -2^17).
    wire [MW-1:0] qc  = s4_q[MW-1] ? {1'b0, {(MW-1){1'b1}}} : s4_q;
    wire          e17 = qc[34];                       // e = 16 + e17
    wire [9:0]    idx = e17 ? qc[35:26] : qc[33:24];  // floor(m * 256)

    (* romstyle = "M10K" *) reg [DW-1:0] rom [0:1023];
    initial $readmemh(ROM_INIT, rom);

    reg [DW-1:0]        s5_y;
    reg signed [DW-1:0] s5_re, s5_im;
    reg                 s5_e17, s5_zero;
    always @(posedge clk) begin
        s5_y    <= rom[idx];
        s5_re   <= s4_re;
        s5_im   <= s4_im;
        s5_e17  <= e17;
        s5_zero <= s4_zero;
    end

    // ---------------- S6: Gn * y ----------------
    wire signed [DW-1:0] y_s = s5_y;                  // < 2^17, positive

    reg signed [MW-1:0] s6_pr, s6_pi;
    reg                 s6_e17, s6_zero;
    always @(posedge clk) begin
        s6_pr   <= s5_re * y_s;
        s6_pi   <= s5_im * y_s;
        s6_e17  <= s5_e17;
        s6_zero <= s5_zero;
    end

    // ---------------- S7: round, clamp ----------------
    // Convergent rounding of x by s bits, clamp to +/-OMAX.
    function [DW-1:0] rnd_clamp;
        input signed [MW-1:0] x;
        input integer         s;
        reg   signed [MW:0]   t;
        begin
            t = x;
            t = t + ((1 <<< (s-1)) - 1) + ((x >>> s) & 1);
            t = t >>> s;
            if (t >  OMAX)     rnd_clamp =  OMAX;
            else if (t < -OMAX) rnd_clamp = -OMAX;
            else               rnd_clamp = t[DW-1:0];
        end
    endfunction

    always @(posedge clk) begin
        g_re <= s6_zero ? {DW{1'b0}} : rnd_clamp(s6_pr, 18 + s6_e17);
        g_im <= s6_zero ? {DW{1'b0}} : rnd_clamp(s6_pi, 18 + s6_e17);
    end
endmodule
