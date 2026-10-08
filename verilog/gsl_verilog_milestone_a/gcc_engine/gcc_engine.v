`timescale 1ns/1ps

// ============================================================
//  gcc_engine.v
//  GCC-PHAT core: one job = one 4-channel window (SW or MB) in,
//  6 coarse TDOAs + lag windows out.
//
//  Contains: gcc_loader, fft_buf_bank (5 buffers), fft2048,
//            xspec_phat, peak_search.
//  Rings (cap_ring) live outside: they are written by the stream
//  path (fe_cond / mb_decim) and read here by gcc_loader.
//
//  Job sequence:
//    LOAD   ring -> Z0, Z1 (packed, bit-reversed)          ~2,052
//    FFT0   Z0 forward, latch exp                          11,321
//    FFT1   Z1 forward, latch exp                          11,321
//    XSPEC  Z0, Z1 -> W0, W1, W2                            6,167
//    IFFT0  W0 inverse                                     11,321
//    IFFT1  W1 inverse                                     11,321
//    IFFT2  W2 inverse                                     11,321
//    PEAK   W0..W2 -> lags, peaks, sidelobes, lag window      406
//    DONE   res_valid                                    ~65,250 total
//
//  Buffer client-port ownership (fft_buf_bank v2):
//    buffers 0,1  write: gcc_loader (LOAD)   read: xspec_phat (XSPEC)
//    buffers 2..4 write: xspec_phat (XSPEC)  read: peak_search (PEAK)
//    any buffer   read : debug port (IDLE only)
//  The FFT engine uses the bank's separate engine port.
//
//  Per-type configuration (SW / MB) for band mask and floor; window
//  start / length / max lag come with each job (from evt_sched).
// ============================================================
module gcc_engine #(
    parameter TW_INIT_RE = "tw2048_re.hex",
    parameter TW_INIT_IM = "tw2048_im.hex",
    parameter RSQRT_INIT = "rsqrt1024.hex"
)(
    input  wire          clk,
    input  wire          rst_n,

    // ---- job (from evt_sched) ----
    input  wire          job_valid,
    output wire          job_ready,
    input  wire          job_type,        // 0 = SW (sw_ring), 1 = MB (mb_ring)
    input  wire [15:0]   job_id,
    input  wire [11:0]   job_start,       // ring address of window sample 0
    input  wire [11:0]   job_len,         // 1 .. 2048
    input  wire [6:0]    job_maxlag,

    // ---- configuration (from csr) ----
    input  wire [10:0]   cfg_sw_band_lo,
    input  wire [10:0]   cfg_sw_band_hi,
    input  wire [6:0]    cfg_sw_floor,
    input  wire [10:0]   cfg_mb_band_lo,
    input  wire [10:0]   cfg_mb_band_hi,
    input  wire [6:0]    cfg_mb_floor,

    // ---- capture rings ----
    output wire [11:0]   sw_ring_raddr,
    input  wire [71:0]   sw_ring_rdata,
    output wire [11:0]   mb_ring_raddr,
    input  wire [71:0]   mb_ring_rdata,

    // ---- results (to result_mbox) ----
    output reg           res_valid,
    output reg           res_type,
    output reg  [15:0]   res_id,
    output wire [47:0]   res_lag,         // s8 per pair
    output wire [107:0]  res_pk,          // s18 per pair
    output wire [107:0]  res_side,        // u18 per pair
    output reg  [9:0]    res_exp_z,       // {exp Z1, exp Z0}
    output reg  [14:0]   res_exp_r,       // {exp W2, exp W1, exp W0}
    output reg           res_err,         // fft_err or port conflict during the job
    output wire          lw_we,           // lag window stream
    output wire [7:0]    lw_idx,
    output wire [107:0]  lw_data,

    // ---- debug read port (only while idle) ----
    input  wire          dbg_re,
    input  wire [2:0]    dbg_sel,
    input  wire [10:0]   dbg_raddr,
    output wire [35:0]   dbg_rdata,       // 1 cycle after dbg_raddr

    output wire          busy
);
    localparam [3:0] S_IDLE  = 4'd0, S_LOAD  = 4'd1, S_FFT0  = 4'd2, S_FFT1  = 4'd3,
                     S_XSPEC = 4'd4, S_IFFT0 = 4'd5, S_IFFT1 = 4'd6, S_IFFT2 = 4'd7,
                     S_PEAK  = 4'd8, S_DONE  = 4'd9;

    reg [3:0]  state;
    reg        launched;
    reg        j_type;
    reg [15:0] j_id;
    reg [11:0] j_start, j_len;
    reg [6:0]  j_maxlag;
    reg [4:0]  ez0, ez1, er0, er1, er2;
    reg        err_acc;

    assign job_ready = (state == S_IDLE);
    assign busy      = (state != S_IDLE);

    // ------------------------------------------------------------
    // Sub-blocks
    // ------------------------------------------------------------
    // loader
    wire        ld_busy, ld_done, ld_z_we;
    wire [10:0] ld_z_waddr;
    wire [35:0] ld_z0_wdata, ld_z1_wdata;
    wire [11:0] ld_raddr;
    wire        ld_start = (state == S_LOAD) && !launched;

    assign sw_ring_raddr = ld_raddr;
    assign mb_ring_raddr = ld_raddr;

    gcc_loader u_loader (
        .clk(clk), .rst_n(rst_n),
        .ld_start(ld_start), .ld_base(j_start), .ld_len(j_len),
        .ld_busy(ld_busy), .ld_done(ld_done),
        .ring_raddr(ld_raddr), .ring_rdata(j_type ? mb_ring_rdata : sw_ring_rdata),
        .z_we(ld_z_we), .z_waddr(ld_z_waddr), .z0_wdata(ld_z0_wdata), .z1_wdata(ld_z1_wdata)
    );

    // fft engine
    reg         fft_start;
    reg         fft_inverse;
    reg  [2:0]  fft_buf_sel;
    wire        fft_busy, fft_done, fft_err;
    wire [4:0]  fft_bfp_exp;
    wire        mem_active;
    wire [2:0]  mem_sel;
    wire [9:0]  mem_raddr0, mem_raddr1, mem_waddr0, mem_waddr1;
    wire [35:0] mem_rdata0, mem_rdata1, mem_wdata0, mem_wdata1;
    wire        mem_we0, mem_we1;

    fft2048 #(.TW_INIT_RE(TW_INIT_RE), .TW_INIT_IM(TW_INIT_IM)) u_fft (
        .clk(clk), .rst_n(rst_n),
        .fft_start(fft_start), .fft_inverse(fft_inverse), .fft_buf_sel(fft_buf_sel),
        .fft_busy(fft_busy), .fft_done(fft_done), .fft_bfp_exp(fft_bfp_exp), .fft_err(fft_err),
        .mem_active(mem_active), .mem_sel(mem_sel),
        .mem_raddr0(mem_raddr0), .mem_raddr1(mem_raddr1),
        .mem_rdata0(mem_rdata0), .mem_rdata1(mem_rdata1),
        .mem_we0(mem_we0), .mem_we1(mem_we1),
        .mem_waddr0(mem_waddr0), .mem_waddr1(mem_waddr1),
        .mem_wdata0(mem_wdata0), .mem_wdata1(mem_wdata1)
    );

    // xspec_phat
    wire        xp_busy, xp_done, xp_z_re;
    wire [10:0] xp_z_raddr, xp_w_waddr;
    wire [2:0]  xp_w_we;
    wire [35:0] xp_w_wdata;
    wire        xp_start = (state == S_XSPEC) && !launched;
    wire [179:0] cl_rdata;

    xspec_phat #(.ROM_INIT(RSQRT_INIT)) u_xspec (
        .clk(clk), .rst_n(rst_n),
        .xp_start(xp_start),
        .cfg_band_lo(j_type ? cfg_mb_band_lo : cfg_sw_band_lo),
        .cfg_band_hi(j_type ? cfg_mb_band_hi : cfg_sw_band_hi),
        .cfg_floor  (j_type ? cfg_mb_floor   : cfg_sw_floor),
        .z0_exp(ez0), .z1_exp(ez1),
        .xp_busy(xp_busy), .xp_done(xp_done),
        .z_re(xp_z_re), .z_raddr(xp_z_raddr),
        .z0_rdata(cl_rdata[0 +: 36]), .z1_rdata(cl_rdata[36 +: 36]),
        .w_we(xp_w_we), .w_waddr(xp_w_waddr), .w_wdata(xp_w_wdata)
    );

    // peak_search
    wire        ps_busy, ps_done, ps_w_re;
    wire [10:0] ps_w_raddr;
    wire        ps_start = (state == S_PEAK) && !launched;

    peak_search u_peak (
        .clk(clk), .rst_n(rst_n),
        .ps_start(ps_start), .ps_maxlag(j_maxlag),
        .ps_busy(ps_busy), .ps_done(ps_done),
        .w_re(ps_w_re), .w_raddr(ps_w_raddr),
        .w0_rdata(cl_rdata[72 +: 36]), .w1_rdata(cl_rdata[108 +: 36]), .w2_rdata(cl_rdata[144 +: 36]),
        .pk_lag(res_lag), .pk_val(res_pk), .pk_side(res_side),
        .lw_we(lw_we), .lw_idx(lw_idx), .lw_data(lw_data)
    );

    // ------------------------------------------------------------
    // Buffer bank + client-port mux
    // ------------------------------------------------------------
    wire in_load  = (state == S_LOAD);
    wire in_xspec = (state == S_XSPEC);
    wire in_peak  = (state == S_PEAK);
    wire in_idle  = (state == S_IDLE);

    wire [4:0]   cl_re, cl_we;
    wire [54:0]  cl_raddr, cl_waddr;
    wire [179:0] cl_wdata;
    wire         cl_conflict;

    genvar g;
    generate
        for (g = 0; g < 5; g = g + 1) begin : g_mux
            wire dbg_hit = in_idle && dbg_re && (dbg_sel == g);
            if (g < 2) begin : g_z
                assign cl_re[g]             = in_xspec ? xp_z_re    : dbg_hit;
                assign cl_raddr[g*11 +: 11] = in_xspec ? xp_z_raddr : dbg_raddr;
                assign cl_we[g]             = in_load && ld_z_we;
                assign cl_waddr[g*11 +: 11] = ld_z_waddr;
                assign cl_wdata[g*36 +: 36] = (g == 0) ? ld_z0_wdata : ld_z1_wdata;
            end else begin : g_w
                assign cl_re[g]             = in_peak ? ps_w_re    : dbg_hit;
                assign cl_raddr[g*11 +: 11] = in_peak ? ps_w_raddr : dbg_raddr;
                assign cl_we[g]             = in_xspec && xp_w_we[g-2];
                assign cl_waddr[g*11 +: 11] = xp_w_waddr;
                assign cl_wdata[g*36 +: 36] = xp_w_wdata;
            end
        end
    endgenerate

    fft_buf_bank #(.NBUF(5), .DW(36), .AW(11)) u_bank (
        .clk(clk),
        .eng_active(mem_active), .eng_sel(mem_sel),
        .eng_raddr0(mem_raddr0), .eng_raddr1(mem_raddr1),
        .eng_rdata0(mem_rdata0), .eng_rdata1(mem_rdata1),
        .eng_we0(mem_we0), .eng_we1(mem_we1),
        .eng_waddr0(mem_waddr0), .eng_waddr1(mem_waddr1),
        .eng_wdata0(mem_wdata0), .eng_wdata1(mem_wdata1),
        .cl_re(cl_re), .cl_raddr(cl_raddr), .cl_rdata(cl_rdata),
        .cl_we(cl_we), .cl_waddr(cl_waddr), .cl_wdata(cl_wdata),
        .cl_conflict(cl_conflict)
    );

    reg [2:0] dbg_sel_d;
    always @(posedge clk) dbg_sel_d <= dbg_sel;
    assign dbg_rdata = cl_rdata[dbg_sel_d*36 +: 36];

    // ------------------------------------------------------------
    // Sequencer
    // ------------------------------------------------------------
    // FFT step helper: buffer and direction for each FFT state
    always @* begin
        case (state)
            S_FFT0:  begin fft_buf_sel = 3'd0; fft_inverse = 1'b0; end
            S_FFT1:  begin fft_buf_sel = 3'd1; fft_inverse = 1'b0; end
            S_IFFT0: begin fft_buf_sel = 3'd2; fft_inverse = 1'b1; end
            S_IFFT1: begin fft_buf_sel = 3'd3; fft_inverse = 1'b1; end
            default: begin fft_buf_sel = 3'd4; fft_inverse = 1'b1; end
        endcase
        fft_start = !launched && (state == S_FFT0 || state == S_FFT1 ||
                                  state == S_IFFT0 || state == S_IFFT1 || state == S_IFFT2);
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state     <= S_IDLE;
            launched  <= 1'b0;
            res_valid <= 1'b0;
            err_acc   <= 1'b0;
            j_type    <= 1'b0;
            j_id      <= 16'd0;
            j_start   <= 12'd0;
            j_len     <= 12'd0;
            j_maxlag  <= 7'd0;
            ez0 <= 5'd0; ez1 <= 5'd0; er0 <= 5'd0; er1 <= 5'd0; er2 <= 5'd0;
            res_type  <= 1'b0;
            res_id    <= 16'd0;
            res_exp_z <= 10'd0;
            res_exp_r <= 15'd0;
            res_err   <= 1'b0;
        end else begin
            res_valid <= 1'b0;
            if (cl_conflict) err_acc <= 1'b1;

            case (state)
                S_IDLE: if (job_valid) begin
                    j_type   <= job_type;
                    j_id     <= job_id;
                    j_start  <= job_start;
                    j_len    <= job_len;
                    j_maxlag <= job_maxlag;
                    err_acc  <= 1'b0;
                    launched <= 1'b0;
                    state    <= S_LOAD;
                end

                S_LOAD:  if (!launched) launched <= 1'b1;
                         else if (ld_done) begin launched <= 1'b0; state <= S_FFT0; end

                S_FFT0, S_FFT1, S_IFFT0, S_IFFT1, S_IFFT2:
                    if (!launched) launched <= 1'b1;
                    else if (fft_done) begin
                        launched <= 1'b0;
                        if (fft_err) err_acc <= 1'b1;
                        case (state)
                            S_FFT0:  begin ez0 <= fft_bfp_exp; state <= S_FFT1;  end
                            S_FFT1:  begin ez1 <= fft_bfp_exp; state <= S_XSPEC; end
                            S_IFFT0: begin er0 <= fft_bfp_exp; state <= S_IFFT1; end
                            S_IFFT1: begin er1 <= fft_bfp_exp; state <= S_IFFT2; end
                            default: begin er2 <= fft_bfp_exp; state <= S_PEAK;  end
                        endcase
                    end

                S_XSPEC: if (!launched) launched <= 1'b1;
                         else if (xp_done) begin launched <= 1'b0; state <= S_IFFT0; end

                S_PEAK:  if (!launched) launched <= 1'b1;
                         else if (ps_done) begin launched <= 1'b0; state <= S_DONE; end

                S_DONE: begin
                    res_valid <= 1'b1;
                    res_type  <= j_type;
                    res_id    <= j_id;
                    res_exp_z <= {ez1, ez0};
                    res_exp_r <= {er2, er1, er0};
                    res_err   <= err_acc;
                    state     <= S_IDLE;
                end

                default: state <= S_IDLE;
            endcase
        end
    end
endmodule
