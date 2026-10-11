`timescale 1ns/1ps

// ============================================================
//  gsl_core.v  -- Milestone B top
//  Continuous stream path + GCC-PHAT engine:
//
//   replay_src -> timebase -> fe_cond -+-> sw_ring (cap_ring) ------+
//                                      +-> sw_det ---------+        |
//                                      +-> mb_decim -+-> mb_ring ---+
//                                                    +-> mb_det -+  |
//                                     sw_det/mb_det -> evt_sched -> gcc_engine
//
//  Outputs: detector records (with event ids), dispatched / late jobs,
//  engine results + lag windows, and the x / y streams for verification.
//  Configuration is on ports here; csr (Milestone C) will drive them.
// ============================================================
module gsl_core #(
    parameter FRAME_DIV  = 1000,
    parameter FD_INIT    = "fe_fdelay.hex",
    parameter DEC_INIT   = "decim63.hex",
    parameter TW_INIT_RE = "tw2048_re.hex",
    parameter TW_INIT_IM = "tw2048_im.hex",
    parameter RSQRT_INIT = "rsqrt1024.hex"
)(
    input  wire          clk,
    input  wire          rst_n,

    // ---- source ----
    input  wire          run,
    input  wire          src_sel,
    input  wire          rp_wr_en,
    input  wire [63:0]   rp_wr_data,
    output wire [10:0]   rp_level,
    output wire          rp_underflow,
    input  wire          adc_valid,
    input  wire [63:0]   adc_data,

    // ---- configuration ----
    input  wire [63:0]   cfg_fe_gain,
    input  wire [63:0]   cfg_fe_offset,
    input  wire [35:0]   cfg_fe_dtrim,
    input  wire          cfg_fe_dcblk_en,
    input  wire [15:0]   cfg_sw_thr_mult2,
    input  wire [31:0]   cfg_sw_noise_min,
    input  wire [15:0]   cfg_sw_holdoff,
    input  wire [15:0]   cfg_sw_dt_min,
    input  wire [8:0]    cfg_sw_feat_win,
    input  wire [7:0]    cfg_sw_sym,
    input  wire [7:0]    cfg_sw_tn_max,
    input  wire [7:0]    cfg_sw_pk_mult,
    input  wire [7:0]    cfg_mb_thr_mult,
    input  wire [7:0]    cfg_mb_settle_mult,
    input  wire [15:0]   cfg_mb_dwell,
    input  wire [15:0]   cfg_mb_holdoff,
    input  wire [31:0]   cfg_mb_noise_min,
    input  wire [15:0]   cfg_mb_settle_timeout,
    input  wire [11:0]   cfg_sw_pre,
    input  wire [11:0]   cfg_sw_len,
    input  wire [6:0]    cfg_sw_maxlag,
    input  wire [11:0]   cfg_mb_pre,
    input  wire [11:0]   cfg_mb_len,
    input  wire [6:0]    cfg_mb_maxlag,
    input  wire [7:0]    cfg_mb_mask_pre,
    input  wire [7:0]    cfg_mb_mask_post,
    input  wire [15:0]   cfg_sw_mask_pre,
    input  wire [15:0]   cfg_sw_mask_post,
    input  wire [10:0]   cfg_sw_band_lo,
    input  wire [10:0]   cfg_sw_band_hi,
    input  wire [6:0]    cfg_sw_floor,
    input  wire [10:0]   cfg_mb_band_lo,
    input  wire [10:0]   cfg_mb_band_hi,
    input  wire [6:0]    cfg_mb_floor,

    // ---- streams (verification / debug) ----
    output wire          x_valid,
    output wire [71:0]   x_data,
    output wire [63:0]   x_ts,
    output wire          y_valid,
    output wire [71:0]   y_data,
    output wire [63:0]   y_ts,
    output wire [31:0]   y_idx,

    // ---- detector records ----
    output wire          sw_rec_valid,
    output wire [15:0]   sw_rec_id,
    output wire [63:0]   sw_rec_ts,
    output wire [1:0]    sw_rec_ref,
    output wire [215:0]  sw_rec_feat,
    output wire          mb_rec_valid,
    output wire [15:0]   mb_rec_id,
    output wire [63:0]   mb_rec_ts,
    output wire [31:0]   mb_rec_idx,
    output wire [1:0]    mb_rec_ref,
    output wire          mb_rec_reliable,

    // ---- scheduling ----
    output wire          job_valid,
    output wire          job_type,
    output wire [15:0]   job_id,
    output wire [11:0]   job_start,
    output wire          late_valid,
    output wire [15:0]   late_id,
    output wire [15:0]   late_cnt,
    output wire [15:0]   ovf_cnt,

    // ---- engine results ----
    output wire          res_valid,
    output wire          res_type,
    output wire [15:0]   res_id,
    output wire [47:0]   res_lag,
    output wire [107:0]  res_pk,
    output wire [107:0]  res_side,
    output wire [9:0]    res_exp_z,
    output wire [14:0]   res_exp_r,
    output wire          res_err,
    output wire          lw_we,
    output wire [7:0]    lw_idx,
    output wire [107:0]  lw_data,

    output wire          busy               // scheduler or engine still working
);
    // ---------------- source + timebase ----------------
    wire        frame_tick, src_valid;
    wire [63:0] src_data, src_ts, ts64;
    wire        rp_full;

    replay_src #(.FRAME_DIV(FRAME_DIV), .FIFO_AW(10)) u_src (
        .clk(clk), .rst_n(rst_n), .run(run), .src_sel(src_sel),
        .rp_wr_en(rp_wr_en), .rp_wr_data(rp_wr_data), .rp_level(rp_level), .rp_full(rp_full),
        .adc_valid(adc_valid), .adc_data(adc_data),
        .frame_tick(frame_tick), .src_valid(src_valid), .src_data(src_data),
        .rp_underflow(rp_underflow)
    );

    timebase u_tb (.clk(clk), .rst_n(rst_n), .src_valid(src_valid), .ts_clear(1'b0),
                   .src_ts(src_ts), .ts64(ts64));

    // ---------------- conditioning ----------------
    fe_cond #(.FD_INIT(FD_INIT)) u_fe (
        .clk(clk), .rst_n(rst_n),
        .src_valid(src_valid), .src_data(src_data), .src_ts(src_ts),
        .cfg_gain(cfg_fe_gain), .cfg_offset(cfg_fe_offset), .cfg_dtrim(cfg_fe_dtrim),
        .cfg_dcblk_en(cfg_fe_dcblk_en),
        .x_valid(x_valid), .x_data(x_data), .x_ts(x_ts)
    );

    // ---------------- SW path ----------------
    wire [11:0] sw_raddr, mb_raddr;
    wire [71:0] sw_rdata, mb_rdata;

    cap_ring u_sw_ring (.clk(clk), .we(x_valid), .waddr(x_ts[11:0]), .wdata(x_data),
                        .raddr(sw_raddr), .rdata(sw_rdata));

    wire sw_blank, sw_nv, mb_nv;
    wire sw_conf_valid, sw_qual;
    wire [63:0] sw_conf_ts;
    sw_det u_swd (
        .clk(clk), .rst_n(rst_n),
        .x_valid(x_valid), .x_data(x_data), .x_ts(x_ts),
        .cfg_thr_mult2(cfg_sw_thr_mult2), .cfg_noise_min(cfg_sw_noise_min),
        .cfg_holdoff(cfg_sw_holdoff), .cfg_dt_min(cfg_sw_dt_min), .cfg_feat_win(cfg_sw_feat_win),
        .cfg_sym(cfg_sw_sym), .cfg_tn_max(cfg_sw_tn_max), .cfg_pk_mult(cfg_sw_pk_mult),
        .rec_valid(sw_rec_valid), .rec_ts(sw_rec_ts), .rec_ref(sw_rec_ref), .rec_feat(sw_rec_feat),
        .conf_valid(sw_conf_valid), .conf_ts(sw_conf_ts), .sw_qual(sw_qual),
        .sw_blank(sw_blank), .noise_valid(sw_nv)
    );

    // ---------------- MB path ----------------
    mb_decim #(.TAP_INIT(DEC_INIT)) u_dec (
        .clk(clk), .rst_n(rst_n),
        .x_valid(x_valid), .x_data(x_data), .x_ts(x_ts),
        .y_valid(y_valid), .y_data(y_data), .y_ts(y_ts), .y_idx(y_idx)
    );

    cap_ring u_mb_ring (.clk(clk), .we(y_valid), .waddr(y_idx[11:0]), .wdata(y_data),
                        .raddr(mb_raddr), .rdata(mb_rdata));

    mb_det u_mbd (
        .clk(clk), .rst_n(rst_n),
        .y_valid(y_valid), .y_data(y_data), .y_ts(y_ts), .y_idx(y_idx), .sw_blank(sw_blank),
        .sw_qual(sw_qual), .sw_conf(sw_conf_valid),
        .cfg_thr_mult(cfg_mb_thr_mult), .cfg_settle_mult(cfg_mb_settle_mult),
        .cfg_dwell(cfg_mb_dwell), .cfg_holdoff(cfg_mb_holdoff),
        .cfg_noise_min(cfg_mb_noise_min), .cfg_settle_timeout(cfg_mb_settle_timeout),
        .rec_valid(mb_rec_valid), .rec_ts(mb_rec_ts), .rec_idx(mb_rec_idx),
        .rec_ref(mb_rec_ref), .rec_reliable(mb_rec_reliable), .noise_valid(mb_nv)
    );

    // ---------------- scheduler ----------------
    wire        job_ready, eng_busy, sch_busy;
    wire [11:0] job_len;
    wire [6:0]  job_maxlag;
    wire [31:0] job_ts0;
    wire [3:0]  job_mask_en;
    wire [127:0] job_mask_lo, job_mask_hi;

    evt_sched u_sch (
        .clk(clk), .rst_n(rst_n),
        .x_valid(x_valid), .x_ts(x_ts), .y_valid(y_valid), .y_idx(y_idx),
        .sw_rec_valid(sw_rec_valid), .sw_rec_ts(sw_rec_ts),
        .sw_conf_valid(sw_conf_valid), .sw_conf_ts(sw_conf_ts),
        .mb_rec_valid(mb_rec_valid), .mb_rec_ts(mb_rec_ts), .mb_rec_idx(mb_rec_idx),
        .cfg_sw_pre(cfg_sw_pre), .cfg_sw_len(cfg_sw_len), .cfg_sw_maxlag(cfg_sw_maxlag),
        .cfg_mb_pre(cfg_mb_pre), .cfg_mb_len(cfg_mb_len), .cfg_mb_maxlag(cfg_mb_maxlag),
        .cfg_mb_mask_pre(cfg_mb_mask_pre), .cfg_mb_mask_post(cfg_mb_mask_post),
        .cfg_sw_mask_pre(cfg_sw_mask_pre), .cfg_sw_mask_post(cfg_sw_mask_post),
        .job_valid(job_valid), .job_ready(job_ready), .job_type(job_type), .job_id(job_id),
        .job_start(job_start), .job_len(job_len), .job_maxlag(job_maxlag),
        .job_ts0(job_ts0), .job_mask_en(job_mask_en), .job_mask_lo(job_mask_lo), .job_mask_hi(job_mask_hi),
        .sw_id(sw_rec_id), .mb_id(mb_rec_id),
        .late_valid(late_valid), .late_id(late_id), .late_cnt(late_cnt), .ovf_cnt(ovf_cnt),
        .busy(sch_busy)
    );

    // ---------------- engine ----------------
    wire [35:0] dbg_rdata;
    gcc_engine #(.TW_INIT_RE(TW_INIT_RE), .TW_INIT_IM(TW_INIT_IM), .RSQRT_INIT(RSQRT_INIT)) u_eng (
        .clk(clk), .rst_n(rst_n),
        .job_valid(job_valid), .job_ready(job_ready), .job_type(job_type), .job_id(job_id),
        .job_start(job_start), .job_len(job_len), .job_maxlag(job_maxlag),
        .job_ts0(job_ts0), .job_mask_en(job_mask_en), .job_mask_lo(job_mask_lo), .job_mask_hi(job_mask_hi),
        .cfg_sw_band_lo(cfg_sw_band_lo), .cfg_sw_band_hi(cfg_sw_band_hi), .cfg_sw_floor(cfg_sw_floor),
        .cfg_mb_band_lo(cfg_mb_band_lo), .cfg_mb_band_hi(cfg_mb_band_hi), .cfg_mb_floor(cfg_mb_floor),
        .sw_ring_raddr(sw_raddr), .sw_ring_rdata(sw_rdata),
        .mb_ring_raddr(mb_raddr), .mb_ring_rdata(mb_rdata),
        .res_valid(res_valid), .res_type(res_type), .res_id(res_id),
        .res_lag(res_lag), .res_pk(res_pk), .res_side(res_side),
        .res_exp_z(res_exp_z), .res_exp_r(res_exp_r), .res_err(res_err),
        .lw_we(lw_we), .lw_idx(lw_idx), .lw_data(lw_data),
        .dbg_re(1'b0), .dbg_sel(3'd0), .dbg_raddr(11'd0), .dbg_rdata(dbg_rdata),
        .busy(eng_busy)
    );

    assign busy = sch_busy || eng_busy;
endmodule
