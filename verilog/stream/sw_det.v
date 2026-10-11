`timescale 1ns/1ps

// ============================================================
//  sw_det.v
//  Shockwave onset detector (streaming version of the SW stage in
//  python/gs_det_signal_prepare.py).
//
//  Spec: verilog/fpga_blocks.md section 2.5
//
//  Per channel c, every sample:
//    sq      = x^2
//    above   = noise_valid && sq > thr[c]       thr = mult2 * max(noise_ms, noise_min)
//    rise[c] = above && !above_prev[c]          (rising crossing, as in Python)
//  Candidate: noise_valid && no holdoff && any rise.
//    rec_ts  = x_ts of the candidate sample (first crossing on any channel)
//    rec_ref = lowest channel index that crossed
//  Qualification (SW vs MB): over QW = 100 samples (1 ms) from the
//  candidate, on the ref channel, with P = max(max x, 0), Q = max(-min x, 0):
//      confirmed as a shockwave only if  64 * min(P, Q) >= cfg_sym * max(P, Q)
//  i.e. an N-wave: the rear swing is as large as the front (ratio ~0.9-1.0
//  at 15-200 m), while a muzzle blast's negative phase is weak (~0.15-0.17)
//  -- independent of range. Default cfg_sym = 32 (0.5).
//  Two further checks reject noise-level candidates (e.g. a weak blast at
//  1500 m in 0.2 Pa noise):
//    duration : |offset(min x) - offset(max x)| <= cfg_tn_max  (60 = 0.6 ms;
//               rifle N-waves 0.15-0.5 ms, measured 18 samples)
//    margin   : max(P, Q)^2 >= cfg_pk_mult * thr  (9 = 3x the trigger
//               amplitude; real shockwaves are ~60x)
//  (An earlier front-steepness test failed below ~50 m, where the blast
//  front becomes as steep as a shockwave.)
//  A rejected candidate (the muzzle blast itself) produces no record, no
//  holdoff and no MB blanking.
//  Python compares |x| > 10 * noise_rms; here x^2 > 100 * noise_ms
//  (cfg_thr_mult2 = mult^2 = 100). noise_min keeps the threshold sane when
//  the ADC noise is below 1 LSB.
//
//  Noise floor: mean square per channel over 512-sample blocks (5.12 ms),
//  accepted only if the whole block had no trigger / holdoff ("clean").
//  The first clean block arms the detector (noise_valid).
//
//  Features (for classification, Milestone C) over cfg_feat_win samples
//  from the trigger sample: per channel positive / negative peak value and
//  offset; T_N ~ neg_off - pos_off, onset refinement ~ pos_off.
//  rec_valid fires at the end of the feature window.
//
//  sw_blank: high for cfg_dt_min samples after each CONFIRMED shockwave (MB guard).
//
//  Timing: per sample, 1 latch + 4 channel cycles + 1 decision cycle;
//  after a clean block, 8 cycles to recompute thresholds.
//
//  Implementation (area-optimised, v2): two shared multipliers. M1 squares x
//  in the channel cycles and the qualification peak in the decision cycle; M2
//  forms the threshold and the margin threshold (pk_mult x thr, stored per
//  channel) after a noise block, and sym x peak in the decision cycle.
//  Behaviour is identical to v1 (byte-identical simulation logs).
// ============================================================
module sw_det (
    input  wire          clk,
    input  wire          rst_n,

    input  wire          x_valid,
    input  wire [71:0]   x_data,
    input  wire [63:0]   x_ts,

    input  wire [15:0]   cfg_thr_mult2,   // mult^2 (100 = 10x RMS)
    input  wire [31:0]   cfg_noise_min,   // minimum noise mean square, LSB^2
    input  wire [15:0]   cfg_holdoff,     // samples
    input  wire [15:0]   cfg_dt_min,      // samples (MB blanking)
    input  wire [8:0]    cfg_feat_win,    // samples, > QW (100)
    input  wire [7:0]    cfg_sym,         // N-wave symmetry threshold, 1/64 units (32 = 0.5)
    input  wire [7:0]    cfg_tn_max,      // max N-wave duration, samples (60)
    input  wire [7:0]    cfg_pk_mult,     // peak margin over threshold, power ratio (9)

    output reg           rec_valid,
    output reg  [63:0]   rec_ts,
    output reg  [1:0]    rec_ref,
    output reg  [215:0]  rec_feat,        // ch c at [54c +: 54]: {pk_pos, pos_off, pk_neg, neg_off}
    output reg           conf_valid,      // SW confirmed (1 ms = 100 samples after onset)
    output reg  [63:0]   conf_ts,         // its onset (same as rec_ts later)
    output wire          sw_blank,
    output wire          sw_qual,         // qualifying a candidate (level)
    output reg           noise_valid
);
    localparam [2:0] S_IDLE = 3'd0, S_CH = 3'd1, S_DEC = 3'd2, S_THR = 3'd3;

    reg [2:0]  st;
    reg [1:0]  c;
    reg [71:0] xl;
    reg [63:0] tsl;

    reg [44:0] bsum [0:3];
    reg [35:0] nms  [0:3];
    reg [51:0] thr  [0:3];
    reg [59:0] thrpk[0:3];                 // cfg_pk_mult * thr
    reg        tp;                         // threshold sub-step
    reg [3:0]  ap, rise;
    reg [8:0]  bcnt;
    reg        clean;
    reg [15:0] hcnt, bkcnt;
    reg        fa;
    reg [8:0]  fi;
    reg signed [17:0] ppos [0:3], pneg [0:3];
    reg        qa;                         // qualifying a candidate
    reg [6:0]  qc;
    reg signed [17:0] qpos, qneg;          // max x, min x on the ref channel
    reg [6:0]  qpo, qno;                   // their offsets in the window
    localparam [6:0] QW1 = 7'd99;          // QW - 1
    reg [8:0]  opos [0:3], oneg [0:3];

    assign sw_blank = (bkcnt != 16'd0);
    assign sw_qual  = qa;

    wire signed [17:0] xc   = xl[71 - 18*c -: 18];
    wire        [17:0] pe    = qpos[17] ? 18'd0 : qpos;              // P = max(max x, 0)
    wire        [18:0] qn    = -{qneg[17], qneg};
    wire        [17:0] ne_   = qneg[17] ? qn[17:0] : 18'd0;          // Q = max(-min x, 0)
    wire        [17:0] mn    = (pe < ne_) ? pe : ne_;
    wire        [17:0] mx    = (pe < ne_) ? ne_ : pe;

    // ---- M1: x^2 (channel cycles) / peak^2 (decision cycle) ----
    wire signed [18:0] m1a  = (st == S_CH) ? {xc[17], xc} : {1'b0, mx};
    wire        [37:0] m1p  = m1a * m1a;
    wire        [35:0] sq   = m1p[35:0];

    wire               abv  = noise_valid && (sq > thr[c]);
    wire        [35:0] ne   = (nms[c] > {4'd0, cfg_noise_min}) ? nms[c] : {4'd0, cfg_noise_min};


    wire               hact = (hcnt != 16'd0);
    wire               trig = noise_valid && !hact && !qa && (rise != 4'd0);

    wire        [25:0] sy_l  = {2'b00, mn, 6'd0};
    wire        [25:0] sy_r  = m2p[25:0];                             // cfg_sym * mx
    wire        [6:0]  tn    = (qno > qpo) ? (qno - qpo) : (qpo - qno);
    wire        [35:0] mx2   = m1p[35:0];                             // mx * mx
    wire               sharp = (mx != 18'd0) && (sy_l >= sy_r)          // N-wave confirmed
                               && ({1'b0, tn} <= cfg_tn_max)
                               && ({24'd0, mx2} >= thrpk[rec_ref]);

    // ---- M2: thresholds (after a noise block) / sym * peak (decision cycle) ----
    reg         [15:0] m2a;
    reg         [51:0] m2b;
    always @* begin
        if (st == S_THR) begin
            m2a = tp ? {8'd0, cfg_pk_mult} : cfg_thr_mult2;
            m2b = tp ? thr[c] : {16'd0, ne};
        end else begin
            m2a = {8'd0, cfg_sym};
            m2b = {34'd0, mx};
        end
    end
    wire        [67:0] m2p   = m2a * m2b;
    wire        [1:0]  rsel = rise[0] ? 2'd0 : rise[1] ? 2'd1 : rise[2] ? 2'd2 : 2'd3;

    integer i;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            st <= S_IDLE; c <= 2'd0; rec_valid <= 1'b0; noise_valid <= 1'b0; conf_valid <= 1'b0;
            qa <= 1'b0; qc <= 7'd0; qpos <= 18'sd0; qneg <= 18'sd0; tp <= 1'b0;
            ap <= 4'd0; rise <= 4'd0; bcnt <= 9'd0; clean <= 1'b1;
            hcnt <= 16'd0; bkcnt <= 16'd0; fa <= 1'b0; fi <= 9'd0;
            for (i = 0; i < 4; i = i + 1) begin
                bsum[i] <= 45'd0; nms[i] <= 36'd0; thr[i] <= {52{1'b1}}; thrpk[i] <= {60{1'b1}};

            end
        end else begin
            rec_valid  <= 1'b0;
            conf_valid <= 1'b0;
            case (st)
                S_IDLE: if (x_valid) begin
                    xl <= x_data; tsl <= x_ts; c <= 2'd0; st <= S_CH;
                end

                S_CH: begin
                    rise[c]  <= abv && !ap[c];
                    ap[c]    <= abv;
                    bsum[c]  <= bsum[c] + sq;
                    if (qa && c == rec_ref) begin
                        if (xc > qpos) begin qpos <= xc; qpo <= qc; end
                        if (xc < qneg) begin qneg <= xc; qno <= qc; end
                    end
                    if (fa) begin
                        if (xc > ppos[c]) begin ppos[c] <= xc; opos[c] <= fi; end
                        if (xc < pneg[c]) begin pneg[c] <= xc; oneg[c] <= fi; end
                    end
                    c <= c + 1'b1;
                    if (c == 2'd3) st <= S_DEC;
                end

                S_DEC: begin
                    if (trig) begin
                        rec_ts  <= tsl;
                        rec_ref <= rsel;
                        qa      <= 1'b1;
                        qc      <= 7'd1;
                        qpos    <= xl[71 - 18*rsel -: 18];
                        qneg    <= xl[71 - 18*rsel -: 18];
                        qpo     <= 7'd0;
                        qno     <= 7'd0;
                        fa      <= 1'b1;
                        fi      <= 9'd1;
                        for (i = 0; i < 4; i = i + 1) begin
                            ppos[i] <= xl[71 - 18*i -: 18]; opos[i] <= 9'd0;
                            pneg[i] <= xl[71 - 18*i -: 18]; oneg[i] <= 9'd0;
                        end
                    end else begin
                        if (hcnt  != 16'd0) hcnt  <= hcnt  - 1'b1;
                        if (bkcnt != 16'd0) bkcnt <= bkcnt - 1'b1;
                        if (qa) begin
                            if (qc == QW1) begin
                                qa <= 1'b0;
                                if (sharp) begin      // confirmed shockwave
                                    hcnt  <= cfg_holdoff;
                                    bkcnt <= cfg_dt_min;
                                    conf_valid <= 1'b1;
                                    conf_ts    <= rec_ts;
                                end else
                                    fa <= 1'b0;       // rejected: drop features, no record
                            end else
                                qc <= qc + 1'b1;
                        end
                        if (fa && !(qa && qc == QW1 && !sharp)) begin
                            if (fi == cfg_feat_win - 1'b1) begin
                                fa        <= 1'b0;
                                rec_valid <= 1'b1;
                                for (i = 0; i < 4; i = i + 1)
                                    rec_feat[54*i +: 54] <= {ppos[i], opos[i], pneg[i], oneg[i]};
                            end else
                                fi <= fi + 1'b1;
                        end
                    end
                    // noise block
                    if (bcnt == 9'd511) begin
                        bcnt  <= 9'd0;
                        clean <= 1'b1;
                        for (i = 0; i < 4; i = i + 1) bsum[i] <= 45'd0;
                        if (clean && !trig && !hact && !qa) begin
                            for (i = 0; i < 4; i = i + 1) nms[i] <= bsum[i][44:9];
                            noise_valid <= 1'b1;
                            c  <= 2'd0;
                            st <= S_THR;
                        end else
                            st <= S_IDLE;
                    end else begin
                        bcnt <= bcnt + 1'b1;
                        if (trig || hact || qa) clean <= 1'b0;
                        st <= S_IDLE;
                    end
                end

                S_THR: begin                      // per channel: thr, then pk_mult * thr
                    if (!tp) begin
                        thr[c] <= m2p[51:0];
                        tp     <= 1'b1;
                    end else begin
                        thrpk[c] <= m2p[59:0];
                        tp <= 1'b0;
                        c  <= c + 1'b1;
                        if (c == 2'd3) st <= S_IDLE;
                    end
                end

                default: st <= S_IDLE;
            endcase
        end
    end
endmodule
