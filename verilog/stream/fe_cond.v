`timescale 1ns/1ps

// ============================================================
//  fe_cond.v
//  Front-end conditioning, per channel:
//      a = x + offset                               (s17)
//      g = sat18( round(a * gain, 14) )             gain Q2.14 (16384 = 1.0)
//      y = sat18( round(sum_k c[k] * g[n-k], 16) )  9-tap fractional delay, Q1.16
//      z = DC blocker( y ),  corner fs / (2*pi*2^13) = 1.94 Hz
//  Output x = z (or y when cfg_dcblk_en = 0), s18; |x| < 2^15 for s16 input.
//
//  Spec: verilog/fpga_blocks.md section 2.4
//
//  Fractional delay: cfg_dtrim in 1/128-sample units, -128 .. +128.
//    dtrim = 0 -> c = unit impulse at tap 4: exact 4-sample delay
//    (transparent). Positive dtrim delays that channel further.
//    Delay error <= 0.025 sample (0.25 us) below 20 kHz.
//  DC blocker (shift-only, no multiplier):
//      acc += ((y - y_prev) << 13) - (acc >>> 13);  z = round(acc, 13)
//  Latency: constant 4 samples on every channel. x_ts = src_ts - 4 is the
//  input-sample index the output corresponds to.
//
//  Timing: sequential over the 4 channels, 12 cycles each. The coefficient
//  ROM has a registered address and a registered output (2-cycle latency):
//    t = 0      gain, shift history, address tap 0
//    t = 1      address tap 1
//    t = 2..10  acc += c[t-2] * g[n-(t-2)], address tap t (t <= 8)
//    t = 11     round, DC block, store
//  -> 48 cycles + 1 per frame (frame = 1000 cycles).
//  Rounding: convergent (round-half-to-even) everywhere.
//
//  Implementation (area-optimised, v2): the 4 x 9 sample history lives in a
//  64 x 18-bit RAM (one 16-slot circular buffer per channel, shared write
//  pointer), read one cycle after its address so it lines up with the
//  coefficient ROM; gain and filter taps share one multiplier (gain at t = 0,
//  taps at t = 2..10). Behaviour is identical to v1 (byte-identical logs).
// ============================================================
module fe_cond #(
    parameter FD_INIT = "fe_fdelay.hex"
)(
    input  wire          clk,
    input  wire          rst_n,

    input  wire          src_valid,
    input  wire [63:0]   src_data,       // {ch0..ch3} s16
    input  wire [63:0]   src_ts,

    input  wire [63:0]   cfg_gain,       // 4 x u16 Q2.14, ch c at [16c +: 16]
    input  wire [63:0]   cfg_offset,     // 4 x s16
    input  wire [35:0]   cfg_dtrim,      // 4 x s9, 1/128 sample, -128 .. 128
    input  wire          cfg_dcblk_en,

    output reg           x_valid,
    output reg  [71:0]   x_data,         // {ch0..ch3} s18
    output reg  [63:0]   x_ts
);
    localparam integer LAT = 4;

    function signed [17:0] sat18;
        input signed [47:0] v;
        if (v > 48'sd131071)       sat18 = 18'sd131071;
        else if (v < -48'sd131072) sat18 = -18'sd131072;
        else                       sat18 = v[17:0];
    endfunction

    function signed [47:0] rconv;            // convergent round of v / 2^s
        input signed [47:0] v;
        input integer       s;
        reg   signed [48:0] tt;
        begin
            tt = v;
            tt = tt + ((49'sd1 <<< (s-1)) - 49'sd1) + ((v >>> s) & 49'sd1);
            rconv = tt >>> s;
        end
    endfunction

    // ---------------- coefficient ROM (p*9 + tap) ----------------
    (* romstyle = "M10K" *) reg [17:0] rom [0:1160];
    initial $readmemh(FD_INIT, rom);
    reg  [10:0]        rom_addr;
    reg  signed [17:0] coef;
    always @(posedge clk) coef <= rom[rom_addr];

    // ---------------- history RAM: address {channel, slot} ----------------
    reg [17:0] hram [0:63];
    reg [3:0]  wp;                           // slot of the newest sample (all channels)
    reg [5:0]  hra;                          // read address (combinational below)
    reg signed [17:0] hq;                    // read data, 1 cycle later
    wire       hwe;
    wire [5:0] hwa;
    wire [17:0] hwd;
    always @(posedge clk) begin
        if (hwe) hram[hwa] <= hwd;
        hq <= hram[hra];
    end

    // ---------------- state ----------------
    reg signed [17:0] yprev [0:3];
    reg signed [39:0] dacc  [0:3];

    reg        busy, fin;
    reg [1:0]  c;
    reg [3:0]  t;
    reg [63:0] in_data, in_ts;
    reg signed [43:0] acc;
    reg [71:0] outv;

    wire signed [15:0] xin   = in_data[63 - 16*c -: 16];
    wire        [15:0] gain  = cfg_gain[16*c +: 16];
    wire signed [15:0] offs  = cfg_offset[16*c +: 16];
    wire signed [8:0]  dtr   = cfg_dtrim[9*c +: 9];
    wire               neg   = dtr[8];
    wire        [8:0]  pabs  = neg ? -dtr : dtr;
    wire        [7:0]  phase = (pabs > 9'd128) ? 8'd128 : pabs[7:0];

    function [10:0] taddr;                   // mirrored taps for dtrim < 0
        input [7:0] p;
        input       ng;
        input [3:0] k;
        taddr = p * 9 + (ng ? (4'd8 - k) : k);
    endfunction

    // one shared multiplier: gain at t = 0, filter taps at t = 2..10
    wire signed [16:0] a_s   = xin + offs;
    wire signed [17:0] m_a   = (t == 4'd0) ? {a_s[16], a_s} : coef;
    wire signed [17:0] m_b   = (t == 4'd0) ? $signed({2'b00, gain}) : hq;
    wire signed [35:0] m_p   = m_a * m_b;
    wire signed [17:0] g_new = sat18(rconv(m_p, 14));

    // history: write the new sample at slot wp+1 at t = 0; read tap k at t = k+1
    assign hwe = busy && (t == 4'd0);
    assign hwa = {c, wp + 4'd1};
    assign hwd = g_new;
    always @* hra = {c, wp + 4'd1 - (t - 4'd1)};

    // final step (t = 11), combinational
    wire signed [17:0] y_fin = sat18(rconv(acc, 16));
    wire signed [39:0] dy    = y_fin - yprev[c];
    wire signed [39:0] d_fin = dacc[c] + (dy <<< 13) - (dacc[c] >>> 13);
    wire signed [17:0] z_fin = sat18(rconv(d_fin, 13));

    integer i;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            busy <= 1'b0; fin <= 1'b0; x_valid <= 1'b0;
            c <= 2'd0; t <= 4'd0; wp <= 4'd0;
            for (i = 0; i < 4; i = i + 1) begin
                yprev[i] <= 18'sd0;
                dacc[i]  <= 40'sd0;
            end
        end else begin
            x_valid <= 1'b0;
            fin     <= 1'b0;
            if (fin) begin
                x_valid <= 1'b1;
                x_data  <= outv;
                x_ts    <= in_ts - LAT;
            end
            if (!busy) begin
                if (src_valid) begin
                    in_data <= src_data;
                    in_ts   <= src_ts;
                    busy    <= 1'b1;
                    c       <= 2'd0;
                    t       <= 4'd0;
                end
            end else if (t == 4'd0) begin            // gain; history write; ROM tap 0
                rom_addr   <= taddr(phase, neg, 4'd0);
                acc        <= 44'sd0;
                t          <= 4'd1;
            end else if (t == 4'd1) begin
                rom_addr <= taddr(phase, neg, 4'd1);
                t        <= 4'd2;
            end else if (t == 4'd11) begin
                dacc[c]  <= d_fin;
                yprev[c] <= y_fin;
                outv[71 - 18*c -: 18] <= cfg_dcblk_en ? z_fin : y_fin;
                if (c == 2'd3) begin
                    busy <= 1'b0;
                    fin  <= 1'b1;
                    wp   <= wp + 4'd1;
                end
                c <= c + 1'b1;
                t <= 4'd0;
            end else begin                           // t = 2 .. 10
                acc <= acc + m_p;
                if (t <= 4'd8) rom_addr <= taddr(phase, neg, t);
                t <= t + 1'b1;
            end
        end
    end

    // history starts at zero (simulation); in hardware the first 8
    // outputs after reset are a start-up transient
    integer kk;
    initial for (kk = 0; kk < 64; kk = kk + 1) hram[kk] = 18'sd0;
endmodule
