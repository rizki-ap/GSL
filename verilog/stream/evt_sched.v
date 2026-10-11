`timescale 1ns/1ps

// ============================================================
//  evt_sched.v
//  Event scheduler: assigns event ids, waits until each event's window
//  is completely in its capture ring, then hands it to gcc_engine.
//
//  Spec: verilog/fpga_blocks.md section 3.1
//
//  SW event (from sw_det): window start = rec_ts - sw_pre (x_ts domain,
//    sw_ring address = index[11:0]); released when the last written
//    sw_ring index >= start + sw_len - 1.
//  MB event (from mb_det): window start = rec_idx - mb_pre (20 kHz y_idx
//    domain, mb_ring address); released when last y_idx >= start + mb_len - 1.
//  Released events queue for the engine (FIFO, depth 4). At dispatch the
//  ring deadline is checked: the window start must not be overwritten
//  during the load, i.e. (last - start) < 4096 - 16. Otherwise the job is
//  dropped and counted as late (late_valid / late_cnt).
//
//  Shockwave masks (burst fire): the onsets of the last 4 confirmed SWs
//  (sw_det conf_valid, 1 ms after onset) are kept. When an MB window is
//  released, every remembered SW whose interval
//      [onset - 5*mask_pre, onset + 5*mask_post]   (mask_* in 20 kHz samples)
//  overlaps the window becomes a mask interval for gcc_loader, which zeros
//  those samples on all 4 channels. Without it, at 900 rpm and 200 m the
//  shockwave of shot k+5 lands 9.5 ms into shot k's MB window (10 deg
//  bearing error in simulation; 0.1 deg with the mask).
//
//  Muzzle-blast masks (short range): the onsets of the last 4 MB records
//  are kept too. When an SW window is released, every remembered MB whose
//  onset lies AFTER the SW onset and whose interval
//      [onset - sw_mask_pre, onset + sw_mask_post]   (100 kHz samples)
//  overlaps the window is masked. At 15 m the blast arrives 8.5 ms after
//  the shockwave, inside its window, and is 4x louder (459 us SW TDOA
//  error unmasked). Blasts before the SW onset (earlier shots in a burst)
//  are not masked: that would cut the shockwave itself.
//
//  Event ids: one counter shared by SW and MB, assigned when the record
//  arrives (SW first if both arrive in the same cycle). The engine
//  result carries the same id (res_id).
//  Overflow (a pending FIFO full) drops the event and counts ovf_cnt.
// ============================================================
module evt_sched (
    input  wire          clk,
    input  wire          rst_n,

    input  wire          x_valid,
    input  wire [63:0]   x_ts,
    input  wire          y_valid,
    input  wire [31:0]   y_idx,

    input  wire          sw_rec_valid,
    input  wire [63:0]   sw_rec_ts,
    input  wire          sw_conf_valid,
    input  wire [63:0]   sw_conf_ts,
    input  wire          mb_rec_valid,
    input  wire [63:0]   mb_rec_ts,
    input  wire [31:0]   mb_rec_idx,

    input  wire [11:0]   cfg_sw_pre,
    input  wire [11:0]   cfg_sw_len,
    input  wire [6:0]    cfg_sw_maxlag,
    input  wire [11:0]   cfg_mb_pre,
    input  wire [11:0]   cfg_mb_len,
    input  wire [6:0]    cfg_mb_maxlag,
    input  wire [7:0]    cfg_mb_mask_pre,     // 20 kHz samples before a SW onset (20 = 1 ms)
    input  wire [7:0]    cfg_mb_mask_post,    // 20 kHz samples after it (60 = 3 ms)
    input  wire [15:0]   cfg_sw_mask_pre,     // 100 kHz samples before an MB onset (100 = 1 ms)
    input  wire [15:0]   cfg_sw_mask_post,    // 100 kHz samples after it (2000 = to window end)

    output reg           job_valid,
    input  wire          job_ready,
    output reg           job_type,
    output reg  [15:0]   job_id,
    output reg  [11:0]   job_start,
    output reg  [11:0]   job_len,
    output reg  [6:0]    job_maxlag,
    output reg  [31:0]   job_ts0,
    output reg  [3:0]    job_mask_en,
    output reg  [127:0]  job_mask_lo,
    output reg  [127:0]  job_mask_hi,

    output wire [15:0]   sw_id,           // id of the SW record valid this cycle
    output wire [15:0]   mb_id,           // id of the MB record valid this cycle

    output reg           late_valid,
    output reg  [15:0]   late_id,
    output reg  [15:0]   late_cnt,
    output reg  [15:0]   ovf_cnt,
    output wire          busy
);
    localparam integer MARGIN = 16;

    // Implementation (area-optimised, v2): the job queue stores mask CENTRES
    // (onsets) and the bounds are formed at dispatch; one set of 4 overlap units
    // serves both job types (one release per cycle); the MB window time base is
    // recomputed at dispatch from its ring start (MB time - 5 * MB index is a
    // constant, taken from the MB records).
    // Behaviour is identical to v1 (byte-identical simulation logs).

    reg [15:0] next_id;
    assign sw_id = next_id;
    assign mb_id = next_id + (sw_rec_valid ? 16'd1 : 16'd0);

    reg [31:0] x_last, y_last, yoff;       // yoff = MB time - 5 * MB index (constant)
    reg        x_seen, y_seen;

    // ---- pending FIFOs (depth 4); for SW the window time base equals the start ----
    reg [15:0] sp_id [0:3]; reg [31:0] sp_st [0:3];                         reg [2:0] sp_w, sp_r;
    reg [15:0] mp_id [0:3]; reg [31:0] mp_st [0:3]; reg [31:0] mp_t0 [0:3]; reg [2:0] mp_w, mp_r;

    // ---- histories: last 4 confirmed SW onsets, last 4 MB onsets ----
    reg [31:0] hs [0:3]; reg [3:0] hv;  reg [1:0] hw;
    reg [31:0] hm [0:3]; reg [3:0] hmv; reg [1:0] hmw;
    wire sp_empty = (sp_w == sp_r), sp_full = (sp_w - sp_r == 3'd4);
    wire mp_empty = (mp_w == mp_r), mp_full = (mp_w - mp_r == 3'd4);

    // ---- job queue (depth 4): type, id, ring start, mask enables + centres ----
    // Depth 4 (v1: 8): a job takes 0.65 ms and windows are released tens of ms
    // apart even at 900 rpm, so occupancy stays <= 1; overflow is held off by
    // jq_full (the release simply waits) and pending overflow is counted.
    reg        jq_t  [0:3]; reg [15:0] jq_id [0:3]; reg [31:0] jq_st [0:3];
    reg [3:0]  jq_me [0:3]; reg [127:0] jq_mc [0:3];
    reg [3:0]  jq_w, jq_r;
    wire jq_empty = (jq_w == jq_r), jq_full = (jq_w - jq_r == 4'd4);

    wire [31:0] sp_age = x_last - sp_st[sp_r[1:0]];
    wire [31:0] mp_age = y_last - mp_st[mp_r[1:0]];
    wire sw_rel = !sp_empty && x_seen && !sp_age[31] && (sp_age >= {20'd0, cfg_sw_len} - 32'd1);
    wire mb_rel = !mp_empty && y_seen && !mp_age[31] && (mp_age >= {20'd0, cfg_mb_len} - 32'd1);
    wire do_sw  = !jq_full && sw_rel;                  // SW has priority
    wire do_mb  = !jq_full && !sw_rel && mb_rel;

    wire        hd_t   = jq_t [jq_r[1:0]];
    wire [31:0] hd_st  = jq_st[jq_r[1:0]];
    wire [31:0] hd_age = (hd_t ? y_last : x_last) - hd_st;
    reg         cool;

    // ---- per-type mask parameters (100 kHz ticks) ----
    wire [31:0] mpre5  = 32'd5 * {24'd0, cfg_mb_mask_pre};
    wire [31:0] mpost5 = 32'd5 * {24'd0, cfg_mb_mask_post};

    // ---- shared overlap units: en = valid && (h + P1 - L >= 0) && (H - (h - P2) >= 0) ----
    //   MB release (SW history): P1 = 5*mb_post, L = window start, H = window end, P2 = 5*mb_pre
    //   SW release (MB history): P1 = 0, L = SW onset + 1, H = window end, P2 = sw_mask_pre
    wire [31:0] w_lo_mb = mp_t0[mp_r[1:0]];
    wire [31:0] w_hi_mb = w_lo_mb + 32'd5 * ({20'd0, cfg_mb_len} - 32'd1);
    wire [31:0] w_lo_sw = sp_st[sp_r[1:0]];
    wire [31:0] sw_on   = w_lo_sw + {20'd0, cfg_sw_pre};
    wire [31:0] w_hi_sw = w_lo_sw + {20'd0, cfg_sw_len} - 32'd1;
    wire [31:0] P1 = do_sw ? 32'd0                     : mpost5;
    wire [31:0] P2 = do_sw ? {16'd0, cfg_sw_mask_pre} : mpre5;
    wire [31:0] Lb = do_sw ? sw_on + 32'd1            : w_lo_mb;
    wire [31:0] Hb = do_sw ? w_hi_sw                   : w_hi_mb;
    wire [3:0]  m_en;
    wire [127:0] m_c;
    genvar gk;
    generate
        for (gk = 0; gk < 4; gk = gk + 1) begin : g_mask
            wire [31:0] h  = do_sw ? hm[gk] : hs[gk];
            wire        hvk = do_sw ? hmv[gk] : hv[gk];
            wire [31:0] d1 = h + P1 - Lb;
            wire [31:0] d2 = Hb - (h - P2);
            assign m_en[gk]         = hvk && !d1[31] && !d2[31];
            assign m_c[32*gk +: 32] = h;
        end
    endgenerate

    // ---- bounds at dispatch ----
    wire [31:0] d_pre  = hd_t ? mpre5  : {16'd0, cfg_sw_mask_pre};
    wire [31:0] d_post = hd_t ? mpost5 : {16'd0, cfg_sw_mask_post};
    wire [127:0] hd_mc = jq_mc[jq_r[1:0]];
    wire [127:0] hd_lo, hd_hi;
    generate
        for (gk = 0; gk < 4; gk = gk + 1) begin : g_bnd
            assign hd_lo[32*gk +: 32] = hd_mc[32*gk +: 32] - d_pre;
            assign hd_hi[32*gk +: 32] = hd_mc[32*gk +: 32] + d_post;
        end
    endgenerate

    assign busy = !sp_empty || !mp_empty || !jq_empty || job_valid || cool;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            next_id <= 16'd0; x_last <= 32'd0; y_last <= 32'd0; yoff <= 32'd0;
            x_seen <= 1'b0; y_seen <= 1'b0;
            sp_w <= 3'd0; sp_r <= 3'd0; mp_w <= 3'd0; mp_r <= 3'd0; jq_w <= 4'd0; jq_r <= 4'd0;
            hv <= 4'd0; hw <= 2'd0; hmv <= 4'd0; hmw <= 2'd0;
            job_valid <= 1'b0; cool <= 1'b0; late_valid <= 1'b0;
            late_cnt <= 16'd0; ovf_cnt <= 16'd0;
            job_type <= 1'b0; job_id <= 16'd0; job_start <= 12'd0; job_len <= 12'd0; job_maxlag <= 7'd0;
            job_ts0 <= 32'd0; job_mask_en <= 4'd0; job_mask_lo <= 128'd0; job_mask_hi <= 128'd0;
            late_id <= 16'd0;
        end else begin
            job_valid  <= 1'b0;
            late_valid <= 1'b0;
            cool       <= job_valid;

            if (x_valid) begin x_last <= x_ts[31:0]; x_seen <= 1'b1; end
            if (y_valid) begin y_last <= y_idx; y_seen <= 1'b1; end

            // ---- histories ----
            if (sw_conf_valid) begin hs[hw] <= sw_conf_ts[31:0]; hv[hw] <= 1'b1; hw <= hw + 1'b1; end
            if (mb_rec_valid)  begin
                hm[hmw] <= mb_rec_ts[31:0]; hmv[hmw] <= 1'b1; hmw <= hmw + 1'b1;
                yoff    <= mb_rec_ts[31:0] - 32'd5 * mb_rec_idx;
            end

            // ---- records in ----
            next_id <= next_id + (sw_rec_valid ? 16'd1 : 16'd0) + (mb_rec_valid ? 16'd1 : 16'd0);
            if (sw_rec_valid) begin
                if (sp_full) ovf_cnt <= ovf_cnt + 1'b1;
                else begin
                    sp_id[sp_w[1:0]] <= sw_id;
                    sp_st[sp_w[1:0]] <= sw_rec_ts[31:0] - {20'd0, cfg_sw_pre};
                    sp_w <= sp_w + 1'b1;
                end
            end
            if (mb_rec_valid) begin
                if (mp_full) ovf_cnt <= ovf_cnt + 1'b1 + (sw_rec_valid && sp_full);
                else begin
                    mp_id[mp_w[1:0]] <= mb_id;
                    mp_st[mp_w[1:0]] <= mb_rec_idx - {20'd0, cfg_mb_pre};
                    mp_t0[mp_w[1:0]] <= mb_rec_ts[31:0] - 32'd5 * {20'd0, cfg_mb_pre};
                    mp_w <= mp_w + 1'b1;
                end
            end

            // ---- release (one per cycle, SW first) ----
            if (do_sw || do_mb) begin
                jq_t [jq_w[1:0]] <= do_mb;
                jq_id[jq_w[1:0]] <= do_sw ? sp_id[sp_r[1:0]] : mp_id[mp_r[1:0]];
                jq_st[jq_w[1:0]] <= do_sw ? sp_st[sp_r[1:0]] : mp_st[mp_r[1:0]];
                jq_me[jq_w[1:0]] <= m_en;
                jq_mc[jq_w[1:0]] <= m_c;
                jq_w <= jq_w + 1'b1;
                if (do_sw) sp_r <= sp_r + 1'b1;
                else       mp_r <= mp_r + 1'b1;
            end

            // ---- dispatch ----
            if (!jq_empty && job_ready && !job_valid && !cool) begin
                jq_r <= jq_r + 1'b1;
                if (hd_age < 32'd4096 - MARGIN) begin
                    job_valid   <= 1'b1;
                    job_type    <= hd_t;
                    job_id      <= jq_id[jq_r[1:0]];
                    job_start   <= hd_st[11:0];
                    job_len     <= hd_t ? cfg_mb_len    : cfg_sw_len;
                    job_maxlag  <= hd_t ? cfg_mb_maxlag : cfg_sw_maxlag;
                    job_ts0     <= hd_t ? (32'd5 * hd_st + yoff) : hd_st;
                    job_mask_en <= jq_me[jq_r[1:0]];
                    job_mask_lo <= hd_lo;
                    job_mask_hi <= hd_hi;
                end else begin
                    late_valid <= 1'b1;
                    late_id    <= jq_id[jq_r[1:0]];
                    late_cnt   <= late_cnt + 1'b1;
                end
            end
        end
    end
endmodule
