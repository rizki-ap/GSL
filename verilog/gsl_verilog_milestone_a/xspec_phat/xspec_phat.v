`timescale 1ns/1ps

// ============================================================
//  xspec_phat.v
//  Unpack 4 real channels from 2 packed spectra, form the 6 GCC-PHAT
//  cross-spectra, and write them packed 2-per-buffer as IFFT inputs.
//
//  Spec: verilog/fpga_blocks.md section 3.4
//
//  Inputs  (natural order, fft2048 outputs)
//    Z0[k] = FFT(ch0 + j*ch1),  Z1[k] = FFT(ch2 + j*ch3)
//  Outputs (BIT-REVERSED addresses, the fft2048 input contract)
//    W0 = G01 + j*G02,  W1 = G03 + j*G12,  W2 = G13 + j*G23
//    with Hermitian completion, so IFFT(Wp) = r_A + j*r_B.
//
//  Unpack (the /2 is applied with round-half-up to stay in s18):
//    X_A[k] = (Z[k] + conj(Z[N-k])) / 2
//    X_B[k] = (Z[k] - conj(Z[N-k])) / 2j
//  Cross-spectrum (lag > 0  <=>  channel b lags channel a):
//    G_ab = PHAT( X_b * conj(X_a) ),  |G_ab| ~ 32767 (0.25 FS)
//  Packing for pair p = (A, B):
//    W[k]   = (GA.re - GB.im) + j(GA.im + GB.re)
//    W[N-k] = (GA.re + GB.im) + j(GB.re - GA.im)
//
//  BFP exponents of Z0 / Z1 do not affect the result (PHAT keeps only
//  the phase, so each spectrum's scale cancels). They are used only
//  for the optional floor (cfg_floor, see xspec_phat_unit.v).
//
//  Band mask: bins k outside [cfg_band_lo, cfg_band_hi] (and their
//  mirrors) are written as 0. Default 0..1024 = all bins, which
//  matches gs_det_tdoa.py.
//
//  Schedule: one 6-cycle slot per bin k = 0..1024.
//    slot s, cycle 0/1 : read Z0,Z1 at k = s, then at N-s
//    slot s, cycle 2   : unpack -> X_next
//    slot s+1, cycles 0..5 : issue pairs 0..5 of bin s to the PHAT unit
//    unit out, odd pair : write W[k], then W[N-k] one cycle later
//  At most one W write per cycle (one-hot w_we).
//
//  Cycles: 1026 slots x 6 + ~11 drain = 6,167  (61.7 us @ 100 MHz)
// ============================================================
module xspec_phat #(
    parameter N        = 2048,
    parameter LOG2N    = 11,
    parameter DW       = 18,
    parameter ROM_INIT = "rsqrt1024.hex"
)(
    input  wire                 clk,
    input  wire                 rst_n,

    // ---- control ----
    input  wire                 xp_start,
    input  wire [LOG2N-1:0]     cfg_band_lo,     // 0 .. N/2
    input  wire [LOG2N-1:0]     cfg_band_hi,     // 0 .. N/2
    input  wire [6:0]           cfg_floor,       // 0 = off
    input  wire [4:0]           z0_exp,          // fft_bfp_exp of Z0
    input  wire [4:0]           z1_exp,          // fft_bfp_exp of Z1
    output reg                  xp_busy,
    output reg                  xp_done,

    // ---- Z0 / Z1 read (same address to both buffers) ----
    output wire                 z_re,
    output wire [LOG2N-1:0]     z_raddr,
    input  wire [2*DW-1:0]      z0_rdata,
    input  wire [2*DW-1:0]      z1_rdata,

    // ---- W0 / W1 / W2 write (one-hot enable, shared addr / data) ----
    output wire [2:0]           w_we,
    output wire [LOG2N-1:0]     w_waddr,
    output wire [2*DW-1:0]      w_wdata
);
    localparam integer AW   = LOG2N;
    localparam integer HALF = N / 2;              // 1024
    localparam integer TAGW = 3 + AW;             // {pair, k}

    localparam [1:0] S_IDLE = 2'd0, S_RUN = 2'd1, S_FLUSH = 2'd2;

    function [AW-1:0] bitrev;
        input [AW-1:0] v;
        integer i;
        for (i = 0; i < AW; i = i + 1)
            bitrev[i] = v[AW-1-i];
    endfunction

    // ------------------------------------------------------------
    // Slot / phase counters
    // ------------------------------------------------------------
    reg [1:0]     state;
    reg [AW:0]    slot;              // 0 .. HALF+1
    reg [2:0]     ph;                // 0 .. 5
    wire          rd_slot  = (state == S_RUN) && (slot <= HALF);

    assign z_re    = rd_slot && (ph <= 3'd1);
    assign z_raddr = (ph == 3'd0) ? slot[AW-1:0] : (N - slot[AW-1:0]);   // mod N

    // ------------------------------------------------------------
    // Capture Z[k] (cycle 1), unpack with Z[N-k] (cycle 2)
    // ------------------------------------------------------------
    reg [2*DW-1:0] zk0, zk1;

    function signed [DW-1:0] half_rnd;           // (s + 1) >>> 1
        input signed [DW:0] s;
        reg   signed [DW+1:0] t;
        begin
            t = s;
            t = (t + 1) >>> 1;
            half_rnd = t[DW-1:0];
        end
    endfunction

    // Unpack one packed spectrum: returns {XA_re, XA_im, XB_re, XB_im}
    function [4*DW-1:0] unpack;
        input [2*DW-1:0] zk, zm;
        reg signed [DW-1:0] zr, zi, mr, mi;
        begin
            zr = zk[2*DW-1:DW]; zi = zk[DW-1:0];
            mr = zm[2*DW-1:DW]; mi = zm[DW-1:0];
            unpack = { half_rnd(zr + mr), half_rnd(zi - mi),      // X_A
                       half_rnd(zi + mi), half_rnd(mr - zr) };    // X_B
        end
    endfunction

    reg [4*DW-1:0] xn01, xn23;       // X_next: {X0, X1}, {X2, X3}
    reg [4*DW-1:0] xc01, xc23;       // X_cur  (stable for a whole slot)
    reg [AW-1:0]   k_cur;
    reg            cur_v;

    always @(posedge clk) begin
        if (rd_slot && ph == 3'd1) begin
            zk0 <= z0_rdata;
            zk1 <= z1_rdata;
        end
        if (rd_slot && ph == 3'd2) begin
            xn01 <= unpack(zk0, z0_rdata);
            xn23 <= unpack(zk1, z1_rdata);
        end
    end

    // X_cur as 4 complex values
    wire signed [DW-1:0] x0r = xc01[4*DW-1:3*DW], x0i = xc01[3*DW-1:2*DW];
    wire signed [DW-1:0] x1r = xc01[2*DW-1:DW],   x1i = xc01[DW-1:0];
    wire signed [DW-1:0] x2r = xc23[4*DW-1:3*DW], x2i = xc23[3*DW-1:2*DW];
    wire signed [DW-1:0] x3r = xc23[2*DW-1:DW],   x3i = xc23[DW-1:0];

    // ------------------------------------------------------------
    // Pair issue: ph = pair index
    //   0:(0,1) 1:(0,2) 2:(0,3) 3:(1,2) 4:(1,3) 5:(2,3)
    // ------------------------------------------------------------
    reg signed [DW-1:0] ia_re, ia_im, ib_re, ib_im;
    reg [5:0]           i_esum;
    always @* begin
        case (ph)
            3'd0:    begin ia_re = x0r; ia_im = x0i; ib_re = x1r; ib_im = x1i; i_esum = z0_exp + z0_exp; end
            3'd1:    begin ia_re = x0r; ia_im = x0i; ib_re = x2r; ib_im = x2i; i_esum = z0_exp + z1_exp; end
            3'd2:    begin ia_re = x0r; ia_im = x0i; ib_re = x3r; ib_im = x3i; i_esum = z0_exp + z1_exp; end
            3'd3:    begin ia_re = x1r; ia_im = x1i; ib_re = x2r; ib_im = x2i; i_esum = z0_exp + z1_exp; end
            3'd4:    begin ia_re = x1r; ia_im = x1i; ib_re = x3r; ib_im = x3i; i_esum = z0_exp + z1_exp; end
            default: begin ia_re = x2r; ia_im = x2i; ib_re = x3r; ib_im = x3i; i_esum = z1_exp + z1_exp; end
        endcase
    end

    wire                 u_out_v, u_busy;
    wire [TAGW-1:0]      u_out_tag;
    wire signed [DW-1:0] u_g_re, u_g_im;

    xspec_phat_unit #(.DW(DW), .TAGW(TAGW), .ROM_INIT(ROM_INIT)) u_phat (
        .clk(clk), .rst_n(rst_n), .cfg_floor(cfg_floor),
        .in_v((state == S_RUN) && cur_v), .in_tag({ph, k_cur}), .in_esum(i_esum),
        .xa_re(ia_re), .xa_im(ia_im), .xb_re(ib_re), .xb_im(ib_im),
        .out_v(u_out_v), .out_tag(u_out_tag),
        .g_re(u_g_re), .g_im(u_g_im), .busy(u_busy)
    );

    // ------------------------------------------------------------
    // Output: band mask, pack pairs, schedule the two writes
    // ------------------------------------------------------------
    wire [2:0]    o_pair = u_out_tag[TAGW-1:AW];
    wire [AW-1:0] o_k    = u_out_tag[AW-1:0];
    wire          o_mask = (o_k < cfg_band_lo) || (o_k > cfg_band_hi);

    wire signed [DW-1:0] gb_re = o_mask ? {DW{1'b0}} : u_g_re;
    wire signed [DW-1:0] gb_im = o_mask ? {DW{1'b0}} : u_g_im;

    reg signed [DW-1:0] ga_re, ga_im;              // held pair A

    reg            wa_v, wb_v, wb2_v;
    reg [1:0]      wa_buf, wb_buf, wb2_buf;
    reg [AW-1:0]   wa_addr, wb_addr, wb2_addr;
    reg [2*DW-1:0] wa_data, wb_data, wb2_data;
    reg            last_seen;

    wire signed [DW-1:0] wk_re  = ga_re - gb_im;   // W[k]
    wire signed [DW-1:0] wk_im  = ga_im + gb_re;
    wire signed [DW-1:0] wm_re  = ga_re + gb_im;   // W[N-k]
    wire signed [DW-1:0] wm_im  = gb_re - ga_im;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            wa_v <= 1'b0; wb_v <= 1'b0; wb2_v <= 1'b0;
        end else begin
            wa_v  <= u_out_v && o_pair[0];
            wb_v  <= u_out_v && o_pair[0];
            wb2_v <= wb_v;
        end
    end

    always @(posedge clk) begin
        if (u_out_v && !o_pair[0]) begin
            ga_re <= gb_re;
            ga_im <= gb_im;
        end
        wa_buf  <= o_pair[2:1];
        wa_addr <= bitrev(o_k);
        wa_data <= {wk_re, wk_im};
        wb_buf  <= o_pair[2:1];
        wb_addr <= bitrev(N - o_k);                 // mod N
        wb_data <= {wm_re, wm_im};
        wb2_buf  <= wb_buf;
        wb2_addr <= wb_addr;
        wb2_data <= wb_data;
    end

    // wa (W[k]) and wb2 (W[N-k]) never coincide: unit outputs alternate
    // even / odd pairs every cycle.
    wire       w_sel_b = !wa_v && wb2_v;
    wire [1:0] w_buf   = w_sel_b ? wb2_buf : wa_buf;
    assign w_we    = (wa_v || wb2_v) ? (3'b001 << w_buf) : 3'b000;
    assign w_waddr = w_sel_b ? wb2_addr : wa_addr;
    assign w_wdata = w_sel_b ? wb2_data : wa_data;

    wire pipe_empty = !u_busy && !u_out_v && !wa_v && !wb_v && !wb2_v;

    // ------------------------------------------------------------
    // FSM
    // ------------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state   <= S_IDLE;
            slot    <= {(AW+1){1'b0}};
            ph      <= 3'd0;
            cur_v   <= 1'b0;
            xp_busy <= 1'b0;
            xp_done <= 1'b0;
        end else begin
            xp_done <= 1'b0;
            case (state)
                S_IDLE: begin
                    if (xp_start) begin
                        slot    <= {(AW+1){1'b0}};
                        ph      <= 3'd0;
                        cur_v   <= 1'b0;
                        xp_busy <= 1'b1;
                        state   <= S_RUN;
                    end
                end

                S_RUN: begin
                    if (ph == 3'd5) begin
                        ph <= 3'd0;
                        // hand bin 'slot' to the issue stage for the next slot
                        cur_v <= (slot <= HALF);
                        if (slot == HALF + 1) begin
                            cur_v <= 1'b0;
                            state <= S_FLUSH;
                        end
                        slot <= slot + 1'b1;
                    end else begin
                        ph <= ph + 1'b1;
                    end
                end

                S_FLUSH: begin
                    if (pipe_empty) begin
                        xp_busy <= 1'b0;
                        xp_done <= 1'b1;
                        state   <= S_IDLE;
                    end
                end

                default: state <= S_IDLE;
            endcase
        end
    end

    // Datapath hand-over to the issue stage (no reset needed)
    always @(posedge clk)
        if (state == S_RUN && ph == 3'd5) begin
            k_cur <= slot[AW-1:0];
            xc01  <= xn01;
            xc23  <= xn23;
        end
endmodule
