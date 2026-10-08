// ============================================================
//  fft2048.v
//  2048-point in-place radix-2 DIT FFT / IFFT engine with
//  block-floating-point (BFP) scaling.
//
//  Spec: verilog/fpga_blocks.md, section 3.3
//
//  Data      : complex s18 (Q1.17), word = {re[17:0], im[17:0]}
//  Twiddles  : s18 (Q1.17), N/2-entry ROM, W = cos - j*sin
//  Memory    : fft_buf_bank (2 banks per buffer, 2R + 2W / cycle)
//
//  ORDERING  : input in BIT-REVERSED order, output in NATURAL order.
//              x[n] must be stored at address bitrev11(n). Producers
//              (gcc_loader, xspec_phat) do this for free by wiring
//              their write address reversed.
//
//  INPUT CONTRACT: |re|, |im| < 2^16 (i.e. < 0.5 full scale).
//              Stage 0 (W = 1) then cannot overflow without scaling.
//              Violations set fft_err (and are saturated).
//
//  BFP RULE  : after each stage, OR the one's-complement magnitudes
//              of every value written. With MSB position p of that
//              OR, the next stage shifts by sh = max(0, p - 14).
//              Radix-2 growth is <= 1 + sqrt(2) per component, so
//              this guarantees no overflow. fft_bfp_exp = sum(sh).
//              True result:  X[k] = out[k] * 2^fft_bfp_exp  (Q1.17)
//
//  INVERSE   : unnormalised IDFT (no 1/N) via
//              IDFT(x) = swap(DFT(swap(x))), swap = exchange re/im.
//              The swap is applied on stage-0 reads and last-stage
//              writes, so the twiddle ROM is forward-only.
//
//  TIMING    : 1 butterfly / cycle, 1024 per stage, + 5-cycle drain
//              per stage (the BFP decision needs the whole stage):
//              11 x 1029 + 2 = 11,321 cycles  (113.2 us @ 100 MHz)
//
//  Pipeline (issue at cycle t):
//    t   : address gen -> RAM / ROM read addresses
//    t+1 : RAM / ROM data, bank un-swap, inverse swap -> input regs
//    t+2 : butterfly M (multiply)
//    t+3 : butterfly A (add / round / saturate)
//    t+4 : bank swap -> RAM write, BFP max update
// ============================================================
`timescale 1ns/1ps

module fft2048 #(
    parameter N          = 2048,
    parameter LOG2N      = 11,
    parameter DW         = 18,
    parameter TW         = 18,
    parameter TW_INIT_RE = "tw2048_re.hex",
    parameter TW_INIT_IM = "tw2048_im.hex"
)(
    input  wire                 clk,
    input  wire                 rst_n,

    // ---- control ----
    input  wire                 fft_start,     // 1-cycle strobe (ignored while busy)
    input  wire                 fft_inverse,   // 0 = FFT, 1 = IFFT
    input  wire [2:0]           fft_buf_sel,   // buffer to transform in place
    output reg                  fft_busy,
    output reg                  fft_done,      // 1-cycle strobe
    output reg  [4:0]           fft_bfp_exp,   // total right shifts, valid at fft_done
    output reg                  fft_err,       // input contract or saturation, valid at fft_done

    // ---- memory master (to fft_buf_bank engine port) ----
    output wire                 mem_active,
    output wire [2:0]           mem_sel,
    output wire [LOG2N-2:0]     mem_raddr0,
    output wire [LOG2N-2:0]     mem_raddr1,
    input  wire [2*DW-1:0]      mem_rdata0,
    input  wire [2*DW-1:0]      mem_rdata1,
    output wire                 mem_we0,
    output wire                 mem_we1,
    output wire [LOG2N-2:0]     mem_waddr0,
    output wire [LOG2N-2:0]     mem_waddr1,
    output wire [2*DW-1:0]      mem_wdata0,
    output wire [2*DW-1:0]      mem_wdata1
);
    localparam integer AW = LOG2N;        // word address width (11)
    localparam integer RW = LOG2N - 1;    // bank row / counter / twiddle width (10)

    localparam [1:0] S_IDLE  = 2'd0,
                     S_RUN   = 2'd1,
                     S_DRAIN = 2'd2,
                     S_DONE  = 2'd3;

    // ------------------------------------------------------------
    // Control state
    // ------------------------------------------------------------
    reg [1:0]     state;
    reg [3:0]     stage;          // 0 .. LOG2N-1
    reg [RW-1:0]  ctr;            // butterfly index within stage
    reg           inv_r;
    reg [2:0]     sel_r;
    reg [1:0]     sh_cur;         // BFP shift used in the current stage
    reg [4:0]     exp_acc;
    reg [DW-2:0]  max_or;         // OR of magnitudes written this stage
    reg           err_in, err_sat;

    wire issue = (state == S_RUN);

    assign mem_active = fft_busy;
    assign mem_sel    = sel_r;

    // ------------------------------------------------------------
    // Address generation (cycle t)
    //   idx_a = ctr with a 0 inserted at bit 'stage'
    //   idx_b = idx_a + 2^stage
    //   k     = (ctr mod 2^stage) << (LOG2N-1-stage)
    // ------------------------------------------------------------
    wire [AW-1:0] one      = {{(AW-1){1'b0}}, 1'b1};
    wire [AW-1:0] c_ext    = {1'b0, ctr};
    wire [AW-1:0] low_mask = (one << stage) - one;
    wire [AW-1:0] idx_a    = ((c_ext & ~low_mask) << 1) | (c_ext & low_mask);
    wire [AW-1:0] idx_b    = idx_a | (one << stage);
    wire [RW-1:0] tw_k     = (ctr & low_mask[RW-1:0]) << (LOG2N - 1 - stage);

    wire          bank_a   = ^idx_a;           // bank of operand A (B is in ~bank_a)
    wire [RW-1:0] row_a    = idx_a[AW-1:1];
    wire [RW-1:0] row_b    = idx_b[AW-1:1];

    assign mem_raddr0 = bank_a ? row_b : row_a;
    assign mem_raddr1 = bank_a ? row_a : row_b;

    wire signed [TW-1:0] w_re, w_im;
    fft2048_twiddle_rom #(
        .TW(TW), .AW(RW), .INIT_RE(TW_INIT_RE), .INIT_IM(TW_INIT_IM)
    ) u_rom (
        .clk(clk), .addr(tw_k), .w_re(w_re), .w_im(w_im)
    );

    // ------------------------------------------------------------
    // Control pipeline
    // ------------------------------------------------------------
    reg           v1, v2, v3, v4;
    reg           ba1, ba2, ba3, ba4;
    reg [RW-1:0]  ra1, ra2, ra3, ra4;
    reg [RW-1:0]  rb1, rb2, rb3, rb4;
    reg           byp1;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            v1 <= 1'b0; v2 <= 1'b0; v3 <= 1'b0; v4 <= 1'b0;
        end else begin
            v1 <= issue; v2 <= v1; v3 <= v2; v4 <= v3;
        end
    end

    always @(posedge clk) begin
        ba1 <= bank_a; ra1 <= row_a; rb1 <= row_b; byp1 <= (tw_k == {RW{1'b0}});
        ba2 <= ba1;    ra2 <= ra1;   rb2 <= rb1;
        ba3 <= ba2;    ra3 <= ra2;   rb3 <= rb2;
        ba4 <= ba3;    ra4 <= ra3;   rb4 <= rb3;
    end

    // ------------------------------------------------------------
    // Cycle t+1: un-swap banks, inverse swap on stage 0, range check
    // ------------------------------------------------------------
    wire [2*DW-1:0] rd_a = ba1 ? mem_rdata1 : mem_rdata0;
    wire [2*DW-1:0] rd_b = ba1 ? mem_rdata0 : mem_rdata1;

    wire swap_in = inv_r && (stage == 4'd0);

    wire signed [DW-1:0] a_re_i = swap_in ? rd_a[DW-1:0]    : rd_a[2*DW-1:DW];
    wire signed [DW-1:0] a_im_i = swap_in ? rd_a[2*DW-1:DW] : rd_a[DW-1:0];
    wire signed [DW-1:0] b_re_i = swap_in ? rd_b[DW-1:0]    : rd_b[2*DW-1:DW];
    wire signed [DW-1:0] b_im_i = swap_in ? rd_b[2*DW-1:DW] : rd_b[DW-1:0];

    // |v| >= 2^16  <=>  the two top bits differ
    wire in_bad = (stage == 4'd0) && v1 &&
                  ((a_re_i[DW-1] ^ a_re_i[DW-2]) | (a_im_i[DW-1] ^ a_im_i[DW-2]) |
                   (b_re_i[DW-1] ^ b_re_i[DW-2]) | (b_im_i[DW-1] ^ b_im_i[DW-2]));

    reg signed [DW-1:0] a_re_r, a_im_r, b_re_r, b_im_r;
    reg signed [TW-1:0] w_re_r, w_im_r;
    reg                 byp_r;

    always @(posedge clk) begin
        a_re_r <= a_re_i; a_im_r <= a_im_i;
        b_re_r <= b_re_i; b_im_r <= b_im_i;
        w_re_r <= w_re;   w_im_r <= w_im;
        byp_r  <= byp1;
    end

    // ------------------------------------------------------------
    // Cycles t+2, t+3: butterfly
    // ------------------------------------------------------------
    wire signed [DW-1:0] p_re, p_im, q_re, q_im;
    wire                 bf_sat;

    fft2048_bfly #(.DW(DW), .TW(TW)) u_bfly (
        .clk(clk),
        .a_re(a_re_r), .a_im(a_im_r), .b_re(b_re_r), .b_im(b_im_r),
        .w_re(w_re_r), .w_im(w_im_r), .bypass(byp_r), .sh(sh_cur),
        .p_re(p_re), .p_im(p_im), .q_re(q_re), .q_im(q_im), .sat(bf_sat)
    );

    // ------------------------------------------------------------
    // Cycle t+4: write back (P -> idx_a, Q -> idx_b)
    // ------------------------------------------------------------
    wire swap_out = inv_r && (stage == LOG2N - 1);

    wire [2*DW-1:0] P = swap_out ? {p_im, p_re} : {p_re, p_im};
    wire [2*DW-1:0] Q = swap_out ? {q_im, q_re} : {q_re, q_im};

    assign mem_we0    = v4;
    assign mem_we1    = v4;
    assign mem_waddr0 = ba4 ? rb4 : ra4;
    assign mem_waddr1 = ba4 ? ra4 : rb4;
    assign mem_wdata0 = ba4 ? Q : P;
    assign mem_wdata1 = ba4 ? P : Q;

    // One's-complement magnitude (|v| for v >= 0, |v|-1 for v < 0)
    function [DW-2:0] mag1;
        input [DW-1:0] v;
        mag1 = v[DW-1] ? ~v[DW-2:0] : v[DW-2:0];
    endfunction

    wire [DW-2:0] wr_mag = mag1(p_re) | mag1(p_im) | mag1(q_re) | mag1(q_im);

    // Next-stage shift: sh = max(0, p - 14), p = MSB of max_or
    wire [1:0] next_sh = max_or[DW-2] ? 2'd2 :
                         max_or[DW-3] ? 2'd1 : 2'd0;

    wire pipe_empty = ~(v1 | v2 | v3 | v4);

    // ------------------------------------------------------------
    // FSM
    // ------------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state       <= S_IDLE;
            stage       <= 4'd0;
            ctr         <= {RW{1'b0}};
            inv_r       <= 1'b0;
            sel_r       <= 3'd0;
            sh_cur      <= 2'd0;
            exp_acc     <= 5'd0;
            max_or      <= {(DW-1){1'b0}};
            err_in      <= 1'b0;
            err_sat     <= 1'b0;
            fft_busy    <= 1'b0;
            fft_done    <= 1'b0;
            fft_bfp_exp <= 5'd0;
            fft_err     <= 1'b0;
        end else begin
            fft_done <= 1'b0;

            if (in_bad)       err_in  <= 1'b1;
            if (v4 && bf_sat) err_sat <= 1'b1;
            if (v4)           max_or  <= max_or | wr_mag;

            case (state)
                S_IDLE: begin
                    if (fft_start) begin
                        inv_r    <= fft_inverse;
                        sel_r    <= fft_buf_sel;
                        stage    <= 4'd0;
                        ctr      <= {RW{1'b0}};
                        sh_cur   <= 2'd0;          // stage 0: input contract
                        exp_acc  <= 5'd0;
                        max_or   <= {(DW-1){1'b0}};
                        err_in   <= 1'b0;
                        err_sat  <= 1'b0;
                        fft_busy <= 1'b1;
                        state    <= S_RUN;
                    end
                end

                S_RUN: begin
                    ctr <= ctr + 1'b1;
                    if (ctr == {RW{1'b1}})
                        state <= S_DRAIN;
                end

                S_DRAIN: begin
                    if (pipe_empty) begin
                        if (stage == LOG2N - 1) begin
                            state <= S_DONE;
                        end else begin
                            stage   <= stage + 1'b1;
                            ctr     <= {RW{1'b0}};
                            sh_cur  <= next_sh;
                            exp_acc <= exp_acc + next_sh;
                            max_or  <= {(DW-1){1'b0}};
                            state   <= S_RUN;
                        end
                    end
                end

                S_DONE: begin
                    fft_done    <= 1'b1;
                    fft_bfp_exp <= exp_acc;
                    fft_err     <= err_in | err_sat;
                    fft_busy    <= 1'b0;
                    state       <= S_IDLE;
                end
            endcase
        end
    end
endmodule
