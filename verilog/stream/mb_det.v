`timescale 1ns/1ps

// ============================================================
//  mb_det.v
//  Muzzle-blast onset detector on the 20 kHz stream (streaming version
//  of the MB stage in python/gs_det_signal_prepare.py).
//
//  Spec: verilog/fpga_blocks.md section 2.7
//
//  Per channel c, every sample:
//    E[c]   = sum of the last W = 16 squares (0.8 ms sliding energy)
//    above  = E > thr_mult * W * ne             (Python: 1 ms energy > 20 * noise^2)
//    quiet  = E < settle_mult * W * ne          (Python: settle at 3x RMS -> 9x power)
//    ne     = max(noise_ms, noise_min)
//  Armed  = noise_valid && !sw_blank && !holdoff && (settled_latch || timeout)
//    settled : all channels quiet for >= cfg_dwell consecutive samples
//              (Python: 5 ms window below 3x RMS). It is LATCHED: once
//              settled, the detector stays armed until it triggers or a new
//              SW blanking starts -- like Python, which finds the settled
//              point once and then searches. (Without the latch, an MB whose
//              front crosses the 9x settle level one sample before the 20x
//              trigger level disarms the detector and is missed.)
//    timeout : >= cfg_settle_timeout samples since sw_blank ended without
//              settling -> arm anyway, event marked unreliable
//  Trigger: armed && any rising crossing of 'above'.
//  Onset refinement: the oldest sample in the 16-sample window of the
//  trigger channel whose square exceeds thr_mult * ne (one always does).
//    rec_idx = y_idx - back,  rec_ts = y_ts - 5 * back  (100 kHz units)
//  rec_reliable = armed by settling (0 = armed by timeout).
//
//  Noise floor: mean square per channel over 128-sample blocks (6.4 ms),
//  accepted only if the whole block had no trigger / holdoff / sw_blank.
//
//  Shockwave qualification (sw_det needs 1 ms to classify a candidate):
//  an MB trigger while sw_qual is high is HELD. If sw_det then confirms a
//  shockwave (sw_conf), the held event was the shockwave's own energy: it is
//  dropped and its holdoff cancelled. If the candidate is rejected (it was
//  the blast itself), the held event is released with its original onset.
//
//  Runs continuously (no one-shot pairing); SW<->MB pairing and the
//  dt_max timeout are done by the processor.
//  Implementation (area-optimised, v2): the 16-sample delay lines of
//  squares live in one 64 x 36-bit RAM (4 circular buffers, shared write
//  pointer), cleared after reset. Per sample: 2 cycles per channel (read the
//  oldest square, then write the newest) + 1 decision cycle; up to 32 scan
//  cycles after a trigger; 12 threshold cycles after a clean block, with one
//  shared 8 x 36 multiplier (x16 is a shift). Sample period = 5000 cycles.
//  Behaviour is identical to v1 (byte-identical simulation logs).
// ============================================================
module mb_det (
    input  wire          clk,
    input  wire          rst_n,

    input  wire          y_valid,
    input  wire [71:0]   y_data,
    input  wire [63:0]   y_ts,
    input  wire [31:0]   y_idx,
    input  wire          sw_blank,
    input  wire          sw_qual,             // sw_det is qualifying a candidate
    input  wire          sw_conf,             // sw_det confirmed a shockwave (pulse)

    input  wire [7:0]    cfg_thr_mult,        // energy ratio (20)
    input  wire [7:0]    cfg_settle_mult,     // power ratio (9 = 3x RMS)
    input  wire [15:0]   cfg_dwell,           // samples (100 = 5 ms)
    input  wire [15:0]   cfg_holdoff,         // samples (400 = 20 ms)
    input  wire [31:0]   cfg_noise_min,       // LSB^2
    input  wire [15:0]   cfg_settle_timeout,  // samples (4000 = 200 ms)

    output reg           rec_valid,
    output reg  [63:0]   rec_ts,
    output reg  [31:0]   rec_idx,
    output reg  [1:0]    rec_ref,
    output reg           rec_reliable,
    output reg           noise_valid
);
    localparam integer W = 16;
    localparam [2:0] S_IDLE = 3'd0, S_CH = 3'd1, S_DEC = 3'd2, S_SCAN = 3'd3,
                     S_THR = 3'd4, S_CLR = 3'd5;

    reg [2:0]  st;
    reg [1:0]  c;
    reg        ph;                         // sub-phase within a channel / scan step
    reg [3:0]  tc;                         // threshold-step counter
    reg [71:0] yl;
    reg [63:0] tsl;
    reg [31:0] idl;
    reg        blk;
    reg        qual_l, pend;

    // ---- delay-line RAM: address {channel, slot}; slot wp holds the newest ----
    reg [35:0] dlram [0:63];
    reg [5:0]  ra;                         // read address (combinational below)
    reg [35:0] rq;                         // read data, 1 cycle after ra
    reg        we;
    reg [5:0]  wa;
    reg [35:0] wd;
    always @(posedge clk) begin
        if (we) dlram[wa] <= wd;
        rq <= dlram[ra];
    end
    reg [3:0]  wp;                         // slot of the newest sample
    reg [5:0]  clr;                        // reset clear counter

    reg [39:0] E    [0:3];
    reg [41:0] bsum [0:3];
    reg [35:0] nms  [0:3];
    reg [47:0] thrE [0:3], thrS [0:3];
    reg [43:0] thrI [0:3];
    reg [3:0]  ap, rise, quiet;
    reg [6:0]  bcnt;
    reg        clean;
    reg [15:0] hcnt, qcnt, age;
    reg        arm_l;

    reg [1:0]  sref;
    reg [4:0]  sa;                         // scan age (15 .. 0)
    reg [3:0]  back;
    reg        found;
    reg        srel;

    wire signed [17:0] yc   = yl[71 - 18*c -: 18];
    wire        [35:0] sq   = yc * yc;
    reg         [35:0] sq_r;
    wire        [39:0] Enew = E[c] + sq_r - rq;            // rq = oldest square of channel c

    wire        [35:0] ne   = (nms[c] > {4'd0, cfg_noise_min}) ? nms[c] : {4'd0, cfg_noise_min};
    wire        [7:0]  tmul = (tc[1:0] == 2'd1) ? cfg_settle_mult : cfg_thr_mult;
    wire        [43:0] tprod = tmul * ne;                   // one shared multiplier

    wire               hact    = (hcnt != 16'd0);
    wire               settled = (qcnt >= cfg_dwell);
    wire               tmo     = (age >= cfg_settle_timeout);
    wire               armed   = noise_valid && !blk && !hact && (arm_l || settled || tmo);
    wire               trig    = armed && (rise != 4'd0);
    wire        [1:0]  rsel    = rise[0] ? 2'd0 : rise[1] ? 2'd1 : rise[2] ? 2'd2 : 2'd3;

    // RAM read address
    always @* begin
        case (st)
            S_SCAN:  ra = {sref, wp - sa[3:0]};
            default: ra = {c, wp + 4'd1};               // S_CH phase 0: oldest of channel c
        endcase
    end

    integer i;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            st <= S_CLR; clr <= 6'd0; c <= 2'd0; ph <= 1'b0; wp <= 4'd0; we <= 1'b0;
            rec_valid <= 1'b0; noise_valid <= 1'b0;
            ap <= 4'd0; rise <= 4'd0; quiet <= 4'd0; bcnt <= 7'd0; clean <= 1'b1;
            hcnt <= 16'd0; qcnt <= 16'd0; age <= 16'd0; arm_l <= 1'b0; pend <= 1'b0;
            for (i = 0; i < 4; i = i + 1) begin
                E[i] <= 40'd0; bsum[i] <= 42'd0; nms[i] <= 36'd0;
                thrE[i] <= {48{1'b1}}; thrS[i] <= 48'd0; thrI[i] <= {44{1'b1}};
            end
        end else begin
            rec_valid <= 1'b0;
            we        <= 1'b0;
            case (st)
                S_CLR: begin                     // zero the delay lines after reset
                    we <= 1'b1; wa <= clr; wd <= 36'd0;
                    clr <= clr + 1'b1;
                    if (clr == 6'd63) st <= S_IDLE;
                end

                S_IDLE: if (y_valid) begin
                    yl <= y_data; tsl <= y_ts; idl <= y_idx; blk <= sw_blank; qual_l <= sw_qual;
                    c <= 2'd0; ph <= 1'b0; st <= S_CH;
                end

                S_CH: begin
                    if (!ph) begin               // read of the oldest square issued (ra)
                        sq_r <= sq;
                        ph   <= 1'b1;
                    end else begin               // rq valid: update, write the newest
                        we <= 1'b1; wa <= {c, wp + 4'd1}; wd <= sq_r;
                        E[c]     <= Enew;
                        rise[c]  <= noise_valid && (Enew > thrE[c]) && !ap[c];
                        ap[c]    <= noise_valid && (Enew > thrE[c]);
                        quiet[c] <= noise_valid && (Enew < thrS[c]);
                        bsum[c]  <= bsum[c] + sq_r;
                        ph <= 1'b0;
                        c  <= c + 1'b1;
                        if (c == 2'd3) begin
                            wp <= wp + 1'b1;
                            st <= S_DEC;
                        end
                    end
                end

                S_DEC: begin
                    qcnt <= (quiet == 4'hF) ? ((qcnt == 16'hFFFF) ? qcnt : qcnt + 1'b1) : 16'd0;
                    age  <= blk ? 16'd0 : ((age == 16'hFFFF) ? age : age + 1'b1);
                    arm_l <= (trig || blk || !noise_valid) ? 1'b0 : (arm_l || settled);
                    if (trig) begin
                        hcnt  <= cfg_holdoff;
                        sref  <= rsel;
                        srel  <= arm_l || settled;
                        sa    <= 5'd15;
                        ph    <= 1'b0;
                        found <= 1'b0;
                        back  <= 4'd0;
                    end else if (hcnt != 16'd0)
                        hcnt <= hcnt - 1'b1;
                    // noise block
                    if (bcnt == 7'd127) begin
                        bcnt  <= 7'd0;
                        clean <= 1'b1;
                        for (i = 0; i < 4; i = i + 1) bsum[i] <= 42'd0;
                        if (clean && !trig && !hact && !blk) begin
                            for (i = 0; i < 4; i = i + 1) nms[i] <= bsum[i][41:7];
                            noise_valid <= 1'b1;
                        end
                    end else begin
                        bcnt <= bcnt + 1'b1;
                        if (trig || hact || blk) clean <= 1'b0;
                    end
                    if (trig)
                        st <= S_SCAN;
                    else if (bcnt == 7'd127 && clean && !hact && !blk) begin
                        c <= 2'd0; tc <= 4'd0; st <= S_THR;
                    end else
                        st <= S_IDLE;
                end

                S_SCAN: begin                    // ph 0: read age sa issued; ph 1: compare
                    if (!ph)
                        ph <= 1'b1;
                    else begin
                        ph <= 1'b0;
                        if (!found && rq > thrI[sref]) begin
                            found <= 1'b1;
                            back  <= sa[3:0];
                        end
                        if (sa == 5'd0) begin
                            if (qual_l) pend <= 1'b1;
                            else        rec_valid <= 1'b1;
                            rec_ref      <= sref;
                            rec_reliable <= srel;
                            begin : emit
                                reg [3:0] b;
                                b = (!found && rq > thrI[sref]) ? 4'd0 : back;
                                rec_idx <= idl - b;
                                rec_ts  <= tsl - 64'd5 * b;
                            end
                            st <= S_IDLE;
                        end else
                            sa <= sa - 1'b1;
                    end
                end

                S_THR: begin                     // tc 0..11: (thrE, thrS, thrI) for c = 0..3
                    case (tc[1:0])
                        2'd0: thrE[c] <= {tprod, 4'd0};
                        2'd1: thrS[c] <= {tprod, 4'd0};
                        default: begin
                            thrI[c] <= tprod;
                            c <= c + 1'b1;
                        end
                    endcase
                    tc <= (tc[1:0] == 2'd2) ? {tc[3:2] + 1'b1, 2'b00} : tc + 1'b1;
                    if (tc[1:0] == 2'd2 && c == 2'd3) st <= S_IDLE;
                end

                default: st <= S_IDLE;
            endcase

            // resolve a held event (after the case: overrides hcnt)
            if (pend) begin
                if (sw_conf) begin
                    pend <= 1'b0;
                    hcnt <= 16'd0;
                end else if (!sw_qual) begin
                    pend      <= 1'b0;
                    rec_valid <= 1'b1;
                end
            end
        end
    end

    integer k;
    initial for (k = 0; k < 64; k = k + 1) dlram[k] = 36'd0;
endmodule
