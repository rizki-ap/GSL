// ============================================================
//  tb_gsl_core.v  -- Milestone B testbench
//
//  Plays a 4-channel recording through gsl_core in real time
//  (FRAME_DIV clocks per sample) and logs everything for
//  milestone_b.py check.
//
//  Plusargs (written by milestone_b.py prepare):
//    +stim=<file>     hex, one 64-bit frame per line {ch0..ch3} s16
//    +nframes=<n>     number of frames in the file
//    +cfg=<file>      43 integers, order documented in milestone_b.py
//    +out=<prefix>    output files <prefix>_x.txt, _y.txt, _ev.txt
//  Parameter FRAME_DIV: 1000 = real time @ 100 MHz (Verilator default);
//  Icarus fallback uses a smaller value for speed (see Makefile).
// ============================================================
`timescale 1ns/1ps

module tb_gsl_core;
    parameter FRAME_DIV = 1000;
    localparam integer MAXF = 524288;              // 5.2 s of audio

    reg clk = 1'b0;
    always #5 clk = ~clk;
    reg rst_n = 1'b0;

    // ---------------- stimulus ----------------
    reg  [63:0] stim [0:MAXF-1];
    integer     nframes, wp;
    reg         run = 1'b0, rp_wr_en = 1'b0;
    reg  [63:0] rp_wr_data = 64'd0;
    wire [10:0] rp_level;
    wire        rp_underflow;

    // ---------------- configuration ----------------
    reg [63:0] fe_gain, fe_offset;
    reg [35:0] fe_dtrim;
    reg        fe_dcblk_en;
    reg [15:0] sw_thr_mult2, sw_holdoff, sw_dt_min;
    reg [31:0] sw_noise_min, mb_noise_min;
    reg [8:0]  sw_feat_win;
    reg [7:0]  sw_sym, sw_tn_max, sw_pk_mult;
    reg [7:0]  mb_thr_mult, mb_settle_mult;
    reg [15:0] mb_dwell, mb_holdoff, mb_settle_timeout;
    reg [11:0] sw_pre, sw_len, mb_pre, mb_len;
    reg [6:0]  sw_maxlag, mb_maxlag, sw_floor, mb_floor;
    reg [7:0]  mb_mask_pre, mb_mask_post;
    reg [15:0] sw_mask_pre, sw_mask_post;
    reg [10:0] sw_band_lo, sw_band_hi, mb_band_lo, mb_band_hi;

    // ---------------- DUT ----------------
    wire         x_valid, y_valid;
    wire [71:0]  x_data, y_data;
    wire [63:0]  x_ts, y_ts;
    wire [31:0]  y_idx;
    wire         sw_rec_valid, mb_rec_valid, mb_rec_reliable;
    wire [15:0]  sw_rec_id, mb_rec_id;
    wire [63:0]  sw_rec_ts, mb_rec_ts;
    wire [31:0]  mb_rec_idx;
    wire [1:0]   sw_rec_ref, mb_rec_ref;
    wire [215:0] sw_rec_feat;
    wire         job_valid, job_type, late_valid;
    wire [15:0]  job_id, late_id, late_cnt, ovf_cnt;
    wire [11:0]  job_start;
    wire         res_valid, res_type, res_err, lw_we, busy;
    wire [15:0]  res_id;
    wire [47:0]  res_lag;
    wire [107:0] res_pk, res_side, lw_data;
    wire [9:0]   res_exp_z;
    wire [14:0]  res_exp_r;
    wire [7:0]   lw_idx;

    gsl_core #(
        .FRAME_DIV(FRAME_DIV),
        .FD_INIT("fe_fdelay.hex"), .DEC_INIT("decim63.hex"),
        .TW_INIT_RE("../fft2048/tw2048_re.hex"), .TW_INIT_IM("../fft2048/tw2048_im.hex"),
        .RSQRT_INIT("../xspec_phat/rsqrt1024.hex")
    ) dut (
        .clk(clk), .rst_n(rst_n),
        .run(run), .src_sel(1'b1), .rp_wr_en(rp_wr_en), .rp_wr_data(rp_wr_data),
        .rp_level(rp_level), .rp_underflow(rp_underflow), .adc_valid(1'b0), .adc_data(64'd0),
        .cfg_fe_gain(fe_gain), .cfg_fe_offset(fe_offset), .cfg_fe_dtrim(fe_dtrim),
        .cfg_fe_dcblk_en(fe_dcblk_en),
        .cfg_sw_thr_mult2(sw_thr_mult2), .cfg_sw_noise_min(sw_noise_min),
        .cfg_sw_holdoff(sw_holdoff), .cfg_sw_dt_min(sw_dt_min), .cfg_sw_feat_win(sw_feat_win), .cfg_sw_sym(sw_sym), .cfg_sw_tn_max(sw_tn_max), .cfg_sw_pk_mult(sw_pk_mult),
        .cfg_mb_thr_mult(mb_thr_mult), .cfg_mb_settle_mult(mb_settle_mult),
        .cfg_mb_dwell(mb_dwell), .cfg_mb_holdoff(mb_holdoff),
        .cfg_mb_noise_min(mb_noise_min), .cfg_mb_settle_timeout(mb_settle_timeout),
        .cfg_sw_pre(sw_pre), .cfg_sw_len(sw_len), .cfg_sw_maxlag(sw_maxlag),
        .cfg_mb_pre(mb_pre), .cfg_mb_len(mb_len), .cfg_mb_maxlag(mb_maxlag),
        .cfg_mb_mask_pre(mb_mask_pre), .cfg_mb_mask_post(mb_mask_post),
        .cfg_sw_mask_pre(sw_mask_pre), .cfg_sw_mask_post(sw_mask_post),
        .cfg_sw_band_lo(sw_band_lo), .cfg_sw_band_hi(sw_band_hi), .cfg_sw_floor(sw_floor),
        .cfg_mb_band_lo(mb_band_lo), .cfg_mb_band_hi(mb_band_hi), .cfg_mb_floor(mb_floor),
        .x_valid(x_valid), .x_data(x_data), .x_ts(x_ts),
        .y_valid(y_valid), .y_data(y_data), .y_ts(y_ts), .y_idx(y_idx),
        .sw_rec_valid(sw_rec_valid), .sw_rec_id(sw_rec_id), .sw_rec_ts(sw_rec_ts),
        .sw_rec_ref(sw_rec_ref), .sw_rec_feat(sw_rec_feat),
        .mb_rec_valid(mb_rec_valid), .mb_rec_id(mb_rec_id), .mb_rec_ts(mb_rec_ts),
        .mb_rec_idx(mb_rec_idx), .mb_rec_ref(mb_rec_ref), .mb_rec_reliable(mb_rec_reliable),
        .job_valid(job_valid), .job_type(job_type), .job_id(job_id), .job_start(job_start),
        .late_valid(late_valid), .late_id(late_id), .late_cnt(late_cnt), .ovf_cnt(ovf_cnt),
        .res_valid(res_valid), .res_type(res_type), .res_id(res_id),
        .res_lag(res_lag), .res_pk(res_pk), .res_side(res_side),
        .res_exp_z(res_exp_z), .res_exp_r(res_exp_r), .res_err(res_err),
        .lw_we(lw_we), .lw_idx(lw_idx), .lw_data(lw_data),
        .busy(busy)
    );

    // ---------------- helpers ----------------
    function integer s18;  input [17:0] v; s18 = $signed(v); endfunction
    function integer s8;   input [7:0]  v; s8  = $signed(v); endfunction

    // ---------------- logging ----------------
    integer fx, fy, fe;
    integer n_sw, n_mb, n_res, n_err;
    reg [15:0] cur_job;
    integer j;

    always @(posedge clk) if (rst_n) begin
        if (x_valid)
            $fdisplay(fx, "%0d %0d %0d %0d %0d", $signed(x_ts),
                      s18(x_data[71:54]), s18(x_data[53:36]), s18(x_data[35:18]), s18(x_data[17:0]));
        if (y_valid)
            $fdisplay(fy, "%0d %0d %0d %0d %0d %0d", $signed(y_ts), y_idx,
                      s18(y_data[71:54]), s18(y_data[53:36]), s18(y_data[35:18]), s18(y_data[17:0]));
        if (sw_rec_valid) begin
            n_sw = n_sw + 1;
            $fwrite(fe, "S %0d %0d %0d", sw_rec_id, $signed(sw_rec_ts), sw_rec_ref);
            for (j = 0; j < 4; j = j + 1)
                $fwrite(fe, " %0d %0d %0d %0d",
                        s18(sw_rec_feat[54*j+36 +: 18]), sw_rec_feat[54*j+27 +: 9],
                        s18(sw_rec_feat[54*j+9  +: 18]), sw_rec_feat[54*j    +: 9]);
            $fwrite(fe, "\n");
        end
        if (mb_rec_valid) begin
            n_mb = n_mb + 1;
            $fdisplay(fe, "M %0d %0d %0d %0d %0d", mb_rec_id, $signed(mb_rec_ts), mb_rec_idx,
                      mb_rec_ref, mb_rec_reliable);
        end
        if (job_valid) begin
            cur_job = job_id;
            $fdisplay(fe, "J %0d %0d %0d %0d", job_id, job_type, job_start, dut.u_sch.job_mask_en);
        end
        if (late_valid)
            $fdisplay(fe, "X %0d", late_id);
        if (lw_we) begin
            $fwrite(fe, "L %0d %0d", cur_job, lw_idx);
            for (j = 0; j < 6; j = j + 1) $fwrite(fe, " %0d", s18(lw_data[18*j +: 18]));
            $fwrite(fe, "\n");
        end
        if (res_valid) begin
            n_res = n_res + 1;
            if (res_err) n_err = n_err + 1;
            $fwrite(fe, "R %0d %0d %0d %0d %0d %0d %0d %0d", res_id, res_type, res_err,
                    res_exp_z[4:0], res_exp_z[9:5], res_exp_r[4:0], res_exp_r[9:5], res_exp_r[14:10]);
            for (j = 0; j < 6; j = j + 1) $fwrite(fe, " %0d", s8(res_lag[8*j +: 8]));
            for (j = 0; j < 6; j = j + 1) $fwrite(fe, " %0d", s18(res_pk[18*j +: 18]));
            for (j = 0; j < 6; j = j + 1) $fwrite(fe, " %0d", res_side[18*j +: 18]);
            $fwrite(fe, "\n");
        end
    end

    // ---------------- main ----------------
    reg [8*256:1] stim_file, cfg_file, out_pfx, fname;
    integer cf, r, k, v, idle;
    integer cfgv [0:42];
    reg [63:0] cycles;

    always @(posedge clk) cycles <= cycles + 64'd1;

    initial begin
        cycles = 0;
        n_sw = 0; n_mb = 0; n_res = 0; n_err = 0; cur_job = 0;
        if (!$value$plusargs("stim=%s", stim_file)) stim_file = "mb_stim.hex";
        if (!$value$plusargs("nframes=%d", nframes)) nframes = 0;
        if (!$value$plusargs("cfg=%s", cfg_file))   cfg_file  = "mb_cfg.txt";
        if (!$value$plusargs("out=%s", out_pfx))    out_pfx   = "mb_out";

        $readmemh(stim_file, stim, 0, nframes - 1);
        cf = $fopen(cfg_file, "r");
        if (cf == 0) begin $display("FATAL: cannot open cfg %0s", cfg_file); $finish; end
        for (k = 0; k < 43; k = k + 1) r = $fscanf(cf, "%d", cfgv[k]);
        $fclose(cf);
        for (k = 0; k < 4; k = k + 1) begin
            v = cfgv[k];     fe_gain  [16*k +: 16] = v[15:0];
            v = cfgv[4+k];   fe_offset[16*k +: 16] = v[15:0];
            v = cfgv[8+k];   fe_dtrim [9*k  +: 9]  = v[8:0];
        end
        fe_dcblk_en = cfgv[12];
        sw_thr_mult2 = cfgv[13]; sw_noise_min = cfgv[14]; sw_holdoff = cfgv[15];
        sw_dt_min = cfgv[16]; sw_feat_win = cfgv[17]; sw_sym = cfgv[18];
        mb_thr_mult = cfgv[19]; mb_settle_mult = cfgv[20]; mb_dwell = cfgv[21];
        mb_holdoff = cfgv[22]; mb_noise_min = cfgv[23]; mb_settle_timeout = cfgv[24];
        sw_pre = cfgv[25]; sw_len = cfgv[26]; sw_maxlag = cfgv[27];
        mb_pre = cfgv[28]; mb_len = cfgv[29]; mb_maxlag = cfgv[30];
        sw_band_lo = cfgv[31]; sw_band_hi = cfgv[32]; sw_floor = cfgv[33];
        mb_band_lo = cfgv[34]; mb_band_hi = cfgv[35]; mb_floor = cfgv[36];
        mb_mask_pre = cfgv[37]; mb_mask_post = cfgv[38];
        sw_mask_pre = cfgv[39]; sw_mask_post = cfgv[40];
        sw_tn_max = cfgv[41]; sw_pk_mult = cfgv[42];

        $sformat(fname, "%0s_x.txt",  out_pfx); fx = $fopen(fname, "w");
        $sformat(fname, "%0s_y.txt",  out_pfx); fy = $fopen(fname, "w");
        $sformat(fname, "%0s_ev.txt", out_pfx); fe = $fopen(fname, "w");

        $display("============================================================");
        $display(" Milestone B: gsl_core, %0d frames (%0.3f s), FRAME_DIV %0d",
                 nframes, nframes / 100000.0, FRAME_DIV);
        $display("============================================================");

        repeat (5) @(negedge clk);
        rst_n = 1'b1;
        repeat (5) @(negedge clk);

        // feed the replay FIFO; start the cadence once it is primed
        wp = 0;
        while (wp < nframes || rp_level != 0) begin
            @(negedge clk);
            rp_wr_en = 1'b0;
            if (wp < nframes && rp_level < 11'd1000) begin
                rp_wr_en   = 1'b1;
                rp_wr_data = stim[wp];
                wp = wp + 1;
            end
            if (!run && (rp_level >= 11'd256 || wp == nframes)) run = 1'b1;
            if ((wp % 20000) == 0 && rp_wr_en) $display("  %0d frames replayed", wp);
        end
        // last frame popped: stop the cadence before the next tick
        @(negedge clk);
        rp_wr_en = 1'b0;
        run = 1'b0;

        // drain: let the scheduler / engine finish
        idle = 0;
        while (idle < 20000) begin
            @(negedge clk);
            idle = busy ? 0 : idle + 1;
        end

        $fdisplay(fe, "E %0d %0d %0d %0d %0d", rp_underflow, late_cnt, ovf_cnt, nframes, cycles);
        $fclose(fx); $fclose(fy); $fclose(fe);
        $display("  SW records %0d, MB records %0d, results %0d (res_err %0d), late %0d, overflow %0d, underflow %0d",
                 n_sw, n_mb, n_res, n_err, late_cnt, ovf_cnt, rp_underflow);
        $display("  simulated %0d cycles", cycles);
        $finish;
    end
endmodule
