// ============================================================
// gsl_mm — gsl_core with its processor interface (Milestone C1)
//
//  Avalon-MM slave (32-bit, word address, 64 KB window, fixed read
//  latency 2, no waitrequest), Avalon-MM read master for replay (32-bit,
//  pipelined), level interrupt. Register map and message format:
//  c1_interface.md.
//
//    0x0000 - 0x00FF  control, status, counters
//    0x0100 - 0x01AB  configuration (43 values, mb_cfg.txt order)
//    0x8000 - 0xFFFF  mailbox ring (8192 words, read-only)
// ============================================================
module gsl_mm #(
    parameter FRAME_DIV  = 1000,
    parameter FD_INIT    = "fe_fdelay.hex",
    parameter DEC_INIT   = "decim63.hex",
    parameter TW_INIT_RE = "tw2048_re.hex",
    parameter TW_INIT_IM = "tw2048_im.hex",
    parameter RSQRT_INIT = "rsqrt1024.hex",
    parameter [31:0] VERSION = 32'h2026_1012
)(
    input  wire          clk,
    input  wire          rst_n,

    // Avalon-MM slave (processor)
    input  wire [13:0]   avs_address,          // word address
    input  wire          avs_read,
    input  wire          avs_write,
    input  wire [31:0]   avs_writedata,
    output reg  [31:0]   avs_readdata,
    output reg           avs_readdatavalid,

    // Avalon-MM read master (replay from SDRAM)
    output wire [31:0]   avm_address,
    output wire          avm_read,
    input  wire          avm_waitrequest,
    input  wire [31:0]   avm_readdata,
    input  wire          avm_readdatavalid,

    output wire          irq,

    // live ADC frames (Milestone D)
    input  wire          adc_valid,
    input  wire [63:0]   adc_data
);
    // ---------------- control registers ----------------
    reg        c_run, c_src, c_irq_en, c_loop;
    reg [3:0]  srst_cnt;                       // soft reset in progress
    wire       clr = (srst_cnt != 4'd0);
    reg        core_rst_n;
    reg [31:0] rp_base, rp_nframes;
    reg [12:0] mb_rd;                          // processor's mailbox read pointer

    // ---------------- configuration (43 values) ----------------
    reg [31:0] cfg [0:42];

    function integer cw;                       // field width
        input integer k;
        begin
            if      (k < 4)  cw = 16;   else if (k < 8)  cw = 16;  else if (k < 12) cw = 9;
            else if (k == 12) cw = 1;   else if (k == 13) cw = 16; else if (k == 14) cw = 32;
            else if (k == 15) cw = 16;  else if (k == 16) cw = 16; else if (k == 17) cw = 9;
            else if (k == 18) cw = 8;   else if (k == 19) cw = 8;  else if (k == 20) cw = 8;
            else if (k == 21) cw = 16;  else if (k == 22) cw = 16; else if (k == 23) cw = 32;
            else if (k == 24) cw = 16;  else if (k == 25) cw = 12; else if (k == 26) cw = 12;
            else if (k == 27) cw = 7;   else if (k == 28) cw = 12; else if (k == 29) cw = 12;
            else if (k == 30) cw = 7;   else if (k == 31) cw = 11; else if (k == 32) cw = 11;
            else if (k == 33) cw = 7;   else if (k == 34) cw = 11; else if (k == 35) cw = 11;
            else if (k == 36) cw = 7;   else if (k == 37) cw = 8;  else if (k == 38) cw = 8;
            else if (k == 39) cw = 16;  else if (k == 40) cw = 16; else cw = 8;
        end
    endfunction

    function csigned;                          // fe_offset, fe_dtrim
        input integer k;
        csigned = (k >= 4 && k < 12);
    endfunction

    function [31:0] cdef;                      // Milestone B defaults
        input integer k;
        begin
            case (k)
                0, 1, 2, 3: cdef = 16384;
                12: cdef = 1;    13: cdef = 100;  14: cdef = 16;   15: cdef = 2000;
                16: cdef = 1000; 17: cdef = 256;  18: cdef = 32;   19: cdef = 20;
                20: cdef = 9;    21: cdef = 100;  22: cdef = 400;  23: cdef = 16;
                24: cdef = 4000; 25: cdef = 300;  26: cdef = 1800; 27: cdef = 100;
                28: cdef = 120;  29: cdef = 1720; 30: cdef = 20;   31: cdef = 20;
                32: cdef = 512;  33: cdef = 0;    34: cdef = 10;   35: cdef = 205;
                36: cdef = 0;    37: cdef = 20;   38: cdef = 60;   39: cdef = 100;
                40: cdef = 2000; 41: cdef = 60;   42: cdef = 9;
                default: cdef = 0;
            endcase
        end
    endfunction

    function [31:0] cmask;
        input integer k;
        integer w;
        begin
            w = cw(k);
            cmask = (w >= 32) ? 32'hFFFF_FFFF : ((32'd1 << w) - 32'd1);
        end
    endfunction

    function [31:0] cread;                     // sign-extend signed fields
        input integer k;
        input [31:0] v;
        integer w;
        begin
            w = cw(k);
            if (csigned(k) && v[w-1]) cread = v | ~cmask(k);
            else                      cread = v;
        end
    endfunction

    // ---------------- core ----------------
    wire        rp_wr_en;
    wire [63:0] rp_wr_data;
    wire [10:0] rp_level;
    wire        rp_underflow;
    wire        x_valid, y_valid;
    wire [71:0] x_data, y_data;
    wire [63:0] x_ts, y_ts;
    wire [31:0] y_idx;
    wire        sw_rec_valid, mb_rec_valid, job_valid, late_valid, res_valid, lw_we;
    wire [15:0] sw_rec_id, mb_rec_id, job_id, late_id, res_id, late_cnt, ovf_cnt;
    wire [63:0] sw_rec_ts, mb_rec_ts;
    wire [1:0]  sw_rec_ref, mb_rec_ref;
    wire [215:0] sw_rec_feat;
    wire [31:0] mb_rec_idx;
    wire        mb_rec_reliable, job_type, res_type, res_err, core_busy;
    wire [11:0] job_start;
    wire [47:0] res_lag;
    wire [107:0] res_pk, res_side, lw_data;
    wire [9:0]  res_exp_z;
    wire [14:0] res_exp_r;
    wire [7:0]  lw_idx;
    wire        sw_nv, mb_nv;
    wire [3:0]  job_men;

    wire dma_done, dma_primed;
    wire replay_done = c_src && dma_done && (rp_level == 11'd0);
    wire core_run    = c_run && (c_src ? (dma_primed && !replay_done) : 1'b1);

    gsl_core #(
        .FRAME_DIV(FRAME_DIV), .FD_INIT(FD_INIT), .DEC_INIT(DEC_INIT),
        .TW_INIT_RE(TW_INIT_RE), .TW_INIT_IM(TW_INIT_IM), .RSQRT_INIT(RSQRT_INIT)
    ) u_core (
        .clk(clk), .rst_n(core_rst_n), .run(core_run), .src_sel(c_src),
        .rp_wr_en(rp_wr_en), .rp_wr_data(rp_wr_data), .rp_level(rp_level), .rp_underflow(rp_underflow),
        .adc_valid(adc_valid), .adc_data(adc_data),
        .cfg_fe_gain({cfg[3][15:0], cfg[2][15:0], cfg[1][15:0], cfg[0][15:0]}),
        .cfg_fe_offset({cfg[7][15:0], cfg[6][15:0], cfg[5][15:0], cfg[4][15:0]}),
        .cfg_fe_dtrim({cfg[11][8:0], cfg[10][8:0], cfg[9][8:0], cfg[8][8:0]}),
        .cfg_fe_dcblk_en(cfg[12][0]),
        .cfg_sw_thr_mult2(cfg[13][15:0]), .cfg_sw_noise_min(cfg[14]), .cfg_sw_holdoff(cfg[15][15:0]),
        .cfg_sw_dt_min(cfg[16][15:0]), .cfg_sw_feat_win(cfg[17][8:0]), .cfg_sw_sym(cfg[18][7:0]),
        .cfg_sw_tn_max(cfg[41][7:0]), .cfg_sw_pk_mult(cfg[42][7:0]),
        .cfg_mb_thr_mult(cfg[19][7:0]), .cfg_mb_settle_mult(cfg[20][7:0]), .cfg_mb_dwell(cfg[21][15:0]),
        .cfg_mb_holdoff(cfg[22][15:0]), .cfg_mb_noise_min(cfg[23]), .cfg_mb_settle_timeout(cfg[24][15:0]),
        .cfg_sw_pre(cfg[25][11:0]), .cfg_sw_len(cfg[26][11:0]), .cfg_sw_maxlag(cfg[27][6:0]),
        .cfg_mb_pre(cfg[28][11:0]), .cfg_mb_len(cfg[29][11:0]), .cfg_mb_maxlag(cfg[30][6:0]),
        .cfg_mb_mask_pre(cfg[37][7:0]), .cfg_mb_mask_post(cfg[38][7:0]),
        .cfg_sw_mask_pre(cfg[39][15:0]), .cfg_sw_mask_post(cfg[40][15:0]),
        .cfg_sw_band_lo(cfg[31][10:0]), .cfg_sw_band_hi(cfg[32][10:0]), .cfg_sw_floor(cfg[33][6:0]),
        .cfg_mb_band_lo(cfg[34][10:0]), .cfg_mb_band_hi(cfg[35][10:0]), .cfg_mb_floor(cfg[36][6:0]),
        .x_valid(x_valid), .x_data(x_data), .x_ts(x_ts),
        .y_valid(y_valid), .y_data(y_data), .y_ts(y_ts), .y_idx(y_idx),
        .sw_rec_valid(sw_rec_valid), .sw_rec_id(sw_rec_id), .sw_rec_ts(sw_rec_ts),
        .sw_rec_ref(sw_rec_ref), .sw_rec_feat(sw_rec_feat),
        .mb_rec_valid(mb_rec_valid), .mb_rec_id(mb_rec_id), .mb_rec_ts(mb_rec_ts),
        .mb_rec_idx(mb_rec_idx), .mb_rec_ref(mb_rec_ref), .mb_rec_reliable(mb_rec_reliable),
        .job_valid(job_valid), .job_type(job_type), .job_id(job_id), .job_start(job_start),
        .late_valid(late_valid), .late_id(late_id), .late_cnt(late_cnt), .ovf_cnt(ovf_cnt),
        .res_valid(res_valid), .res_type(res_type), .res_id(res_id), .res_lag(res_lag),
        .res_pk(res_pk), .res_side(res_side), .res_exp_z(res_exp_z), .res_exp_r(res_exp_r),
        .res_err(res_err),
        .lw_we(lw_we), .lw_idx(lw_idx), .lw_data(lw_data),
        .busy(core_busy),
        .sw_noise_valid(sw_nv), .mb_noise_valid(mb_nv), .job_mask_en_o(job_men)
    );

    // ---------------- replay DMA ----------------
    wire [31:0] rp_pos, rp_loops;
    gsl_replay_dma u_dma (
        .clk(clk), .rst_n(rst_n), .clr(clr),
        .enable(c_run && c_src), .loop(c_loop), .base(rp_base), .nframes(rp_nframes),
        .avm_address(avm_address), .avm_read(avm_read), .avm_waitrequest(avm_waitrequest),
        .avm_readdata(avm_readdata), .avm_readdatavalid(avm_readdatavalid),
        .fifo_level(rp_level), .fifo_we(rp_wr_en), .fifo_wdata(rp_wr_data),
        .pos(rp_pos), .loops(rp_loops), .done(dma_done), .primed(dma_primed)
    );

    // ---------------- mailbox ----------------
    wire [31:0] mb_rdata;
    wire [12:0] mb_wr, mb_level;
    wire [31:0] mb_drops;
    wire        mb_drop_sticky;
    gsl_mbox #(.AW(13)) u_mbox (
        .clk(clk), .rst_n(rst_n), .clr(clr),
        .sw_v(sw_rec_valid), .sw_id(sw_rec_id), .sw_ts(sw_rec_ts), .sw_ref(sw_rec_ref), .sw_feat(sw_rec_feat),
        .mb_v(mb_rec_valid), .mb_id(mb_rec_id), .mb_ts(mb_rec_ts), .mb_idx(mb_rec_idx),
        .mb_ref(mb_rec_ref), .mb_rel(mb_rec_reliable),
        .job_v(job_valid), .job_type(job_type), .job_id(job_id), .job_start(job_start), .job_men(job_men),
        .late_v(late_valid), .late_id(late_id),
        .res_v(res_valid), .res_type(res_type), .res_id(res_id), .res_lag(res_lag),
        .res_pk(res_pk), .res_side(res_side), .res_exp_z(res_exp_z), .res_exp_r(res_exp_r),
        .res_err(res_err), .cfg_sw_maxlag(cfg[27][6:0]), .cfg_mb_maxlag(cfg[30][6:0]),
        .lw_we(lw_we), .lw_idx(lw_idx), .lw_data(lw_data),
        .rd_addr(avs_address[12:0]), .rd_data(mb_rdata), .rd_ptr(mb_rd),
        .wr_ptr(mb_wr), .level(mb_level), .drop_cnt(mb_drops), .drop_sticky(mb_drop_sticky)
    );

    assign irq = c_irq_en && (mb_level != 13'd0);

    // ---------------- status and counters ----------------
    reg [63:0] ts_cur;
    reg [31:0] ts_hi_l;
    reg [15:0] n_sw, n_mb;
    reg [31:0] n_res;
    reg        s_uf, s_drop, s_err;
    reg        uf_d;
    reg [31:0] drops_d;

    // ---------------- Avalon slave ----------------
    reg        rd1;
    reg [13:0] a1;
    wire [7:0] w1 = a1[7:0];

    integer k;
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            c_run <= 1'b0; c_src <= 1'b1; c_irq_en <= 1'b0; c_loop <= 1'b0;
            srst_cnt <= 4'd0; core_rst_n <= 1'b0;
            rp_base <= 32'd0; rp_nframes <= 32'd0; mb_rd <= 13'd0;
            for (k = 0; k < 43; k = k + 1) cfg[k] <= cdef(k);
            ts_cur <= 64'd0; ts_hi_l <= 32'd0; n_sw <= 16'd0; n_mb <= 16'd0; n_res <= 32'd0;
            s_uf <= 1'b0; s_drop <= 1'b0; s_err <= 1'b0; uf_d <= 1'b0; drops_d <= 32'd0;
            rd1 <= 1'b0; avs_readdatavalid <= 1'b0; avs_readdata <= 32'd0;
        end else begin
            core_rst_n <= !clr;
            if (srst_cnt != 4'd0) srst_cnt <= srst_cnt - 1'b1;

            // ---- status ----
            if (x_valid) ts_cur <= x_ts;
            uf_d    <= rp_underflow;
            drops_d <= mb_drops;
            if (clr) begin
                n_sw <= 16'd0; n_mb <= 16'd0; n_res <= 32'd0; mb_rd <= 13'd0;
                s_uf <= 1'b0; s_drop <= 1'b0; s_err <= 1'b0;
            end else begin
                if (sw_rec_valid) n_sw <= n_sw + 1'b1;
                if (mb_rec_valid) n_mb <= n_mb + 1'b1;
                if (res_valid)    n_res <= n_res + 1'b1;
                if (rp_underflow && !uf_d)  s_uf   <= 1'b1;
                if (mb_drops != drops_d)    s_drop <= 1'b1;
                if (res_valid && res_err)   s_err  <= 1'b1;
            end

            // ---- writes ----
            if (avs_write && !avs_address[13]) begin
                if (avs_address[7:6] == 2'b01) begin           // 0x100 + 4k
                    for (k = 0; k < 43; k = k + 1)
                        if (avs_address[5:0] == k) cfg[k] <= avs_writedata & cmask(k);
                end else case (avs_address[7:0])
                    8'd2: begin
                        c_run <= avs_writedata[0]; c_src <= avs_writedata[1];
                        c_irq_en <= avs_writedata[3]; c_loop <= avs_writedata[4];
                        if (avs_writedata[2]) srst_cnt <= 4'd8;
                    end
                    8'd3: begin                                  // W1C
                        if (avs_writedata[0]) s_uf   <= 1'b0;
                        if (avs_writedata[3]) s_drop <= 1'b0;
                        if (avs_writedata[4]) s_err  <= 1'b0;
                    end
                    8'd13: mb_rd      <= avs_writedata[12:0];
                    8'd16: rp_base    <= {avs_writedata[31:3], 3'b000};
                    8'd17: rp_nframes <= avs_writedata;
                    default: ;
                endcase
            end

            // ---- reads: latency 2 ----
            rd1 <= avs_read;
            a1  <= avs_address;
            avs_readdatavalid <= rd1;
            if (rd1) begin
                if (a1[13])
                    avs_readdata <= mb_rdata;
                else if (w1[7:6] == 2'b01) begin
                    avs_readdata <= 32'd0;
                    for (k = 0; k < 43; k = k + 1)
                        if (w1[5:0] == k) avs_readdata <= cread(k, cfg[k]);
                end else case (w1)
                    8'd0:  avs_readdata <= 32'h4753_4C31;           // "GSL1"
                    8'd1:  avs_readdata <= VERSION;
                    8'd2:  avs_readdata <= {27'd0, c_loop, c_irq_en, clr, c_src, c_run};
                    8'd3:  avs_readdata <= {24'd0, 1'b0, replay_done, core_busy, s_err, s_drop, mb_nv, sw_nv, s_uf};
                    8'd4:  begin avs_readdata <= ts_cur[31:0]; ts_hi_l <= ts_cur[63:32]; end
                    8'd5:  avs_readdata <= ts_hi_l;
                    8'd6:  avs_readdata <= {16'd0, late_cnt};
                    8'd7:  avs_readdata <= {16'd0, ovf_cnt};
                    8'd8:  avs_readdata <= mb_drops;
                    8'd9:  avs_readdata <= {n_mb, n_sw};
                    8'd10: avs_readdata <= n_res;
                    8'd12: avs_readdata <= {19'd0, mb_wr};
                    8'd13: avs_readdata <= {19'd0, mb_rd};
                    8'd14: avs_readdata <= {19'd0, mb_level};
                    8'd16: avs_readdata <= rp_base;
                    8'd17: avs_readdata <= rp_nframes;
                    8'd18: avs_readdata <= rp_pos;
                    8'd19: avs_readdata <= rp_loops;
                    8'd20: avs_readdata <= {21'd0, rp_level};
                    default: avs_readdata <= 32'd0;
                endcase
            end
        end
    end
endmodule
