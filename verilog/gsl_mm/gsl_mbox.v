// ============================================================
// gsl_mbox — result mailbox for the processor (Milestone C1)
//
//  Serialises gsl_core's SW / MB records, job dispatches, late jobs and
//  engine results (with their lag windows) into a ring of 32-bit words
//  that the processor reads over Avalon-MM. Message format: c1_interface.md
//  section 4. Header: [31:28] type, [27:16] length (words), [15:0] id.
//
//  - Each source has a holding register; a source's events are at least
//    ~20 ms apart (holdoffs, 0.65 ms jobs) while the longest message takes
//    ~1,430 cycles, so an event never arrives while its own source is still
//    pending. If it ever did, it is dropped and counted. Late-job ids have a
//    4-deep queue (several can fire on consecutive cycles).
//  - The lag window is staged during the peak search (lw_we) and copied
//    into the ring when the result arrives.
//  - A message is written completely or dropped (no space): never torn.
//    The write pointer seen by the processor moves only after the last word.
//  - Fixed priority when several sources are pending: late, result, job,
//    SW, MB. Events that happen within the same few microseconds can
//    therefore appear in a different order than they occurred.
// ============================================================
module gsl_mbox #(
    parameter integer AW = 13                  // ring = 2^AW words
)(
    input  wire          clk,
    input  wire          rst_n,
    input  wire          clr,                  // soft reset

    input  wire          sw_v,
    input  wire [15:0]   sw_id,
    input  wire [63:0]   sw_ts,
    input  wire [1:0]    sw_ref,
    input  wire [215:0]  sw_feat,

    input  wire          mb_v,
    input  wire [15:0]   mb_id,
    input  wire [63:0]   mb_ts,
    input  wire [31:0]   mb_idx,
    input  wire [1:0]    mb_ref,
    input  wire          mb_rel,

    input  wire          job_v,
    input  wire          job_type,
    input  wire [15:0]   job_id,
    input  wire [11:0]   job_start,
    input  wire [3:0]    job_men,

    input  wire          late_v,
    input  wire [15:0]   late_id,

    input  wire          res_v,
    input  wire          res_type,
    input  wire [15:0]   res_id,
    input  wire [47:0]   res_lag,
    input  wire [107:0]  res_pk,
    input  wire [107:0]  res_side,
    input  wire [9:0]    res_exp_z,
    input  wire [14:0]   res_exp_r,
    input  wire          res_err,
    input  wire [6:0]    cfg_sw_maxlag,
    input  wire [6:0]    cfg_mb_maxlag,

    input  wire          lw_we,
    input  wire [7:0]    lw_idx,
    input  wire [107:0]  lw_data,

    // processor side
    input  wire [AW-1:0] rd_addr,              // ring read address (1-cycle latency)
    output reg  [31:0]   rd_data,
    input  wire [AW-1:0] rd_ptr,               // processor's read pointer
    output reg  [AW-1:0] wr_ptr,               // published write pointer
    output wire [AW-1:0] level,
    output reg  [31:0]   drop_cnt,
    output reg           drop_sticky
);
    localparam integer N = 1 << AW;
    localparam [2:0] T_SW = 3'd1, T_MB = 3'd2, T_JOB = 3'd3, T_LATE = 3'd4, T_RES = 3'd5;

    // ---------------- ring RAM ----------------
    reg [31:0] ring [0:N-1];
    reg        rwe;
    reg [AW-1:0] rwa;
    reg [31:0] rwd;
    always @(posedge clk) begin
        if (rwe) ring[rwa] <= rwd;
        rd_data <= ring[rd_addr];
    end

    assign level = wr_ptr - rd_ptr;
    wire [AW:0] free = N - 1 - {1'b0, level};

    // ---------------- lag-window staging RAM ----------------
    reg [107:0] lwram [0:255];
    reg [7:0]   lwa;                           // read address (combinational below)
    reg [107:0] lwq;
    always @(posedge clk) begin
        if (lw_we) lwram[lw_idx] <= lw_data;
        lwq <= lwram[lwa];
    end

    // ---------------- holding registers ----------------
    reg         p_sw, p_mb, p_job, p_res;
    reg [15:0]  h_sw_id, h_mb_id, h_job_id, h_res_id;
    reg [63:0]  h_sw_ts, h_mb_ts;
    reg [1:0]   h_sw_ref, h_mb_ref;
    reg [215:0] h_sw_feat;
    reg [31:0]  h_mb_idx;
    reg         h_mb_rel;
    reg [31:0]  h_job_w;
    reg         h_res_type, h_res_err;
    reg [47:0]  h_res_lag;
    reg [107:0] h_res_pk, h_res_side;
    reg [24:0]  h_res_exp;
    reg [6:0]   h_res_L;
    // late-id queue
    reg [15:0]  lq [0:3];
    reg [2:0]   lq_w, lq_r;
    wire        lq_empty = (lq_w == lq_r), lq_full = (lq_w - lq_r == 3'd4);

    // ---------------- packer ----------------
    localparam [2:0] S_IDLE = 3'd0, S_FIX = 3'd1, S_LWRD = 3'd2, S_LWWR = 3'd3, S_DONE = 3'd4;
    reg [2:0]   st;
    reg [2:0]   sel;                           // message type being written
    reg [4:0]   k;                             // fixed-word index
    reg [4:0]   nfix;                          // fixed words in this message
    reg [11:0]  mlen;                          // message length
    reg [AW-1:0] wa;                           // next write address
    reg [7:0]   m, mlast;                      // lag-window row
    reg [2:0]   pp;                            // pair within the row
    reg [15:0]  cur_late;

    // lengths
    wire [6:0]  Lsel   = res_type ? cfg_mb_maxlag : cfg_sw_maxlag;
    wire [11:0] res_len_of_h = 12'd17 + 12'd6 * ({4'd0, h_res_L, 1'b0} + 12'd1);

    function [31:0] sx18;
        input [17:0] v;
        sx18 = {{14{v[17]}}, v};
    endfunction

    // fixed-word multiplexer
    reg [31:0] fw;
    always @* begin
        fw = 32'd0;
        case (sel)
            T_SW: case (k)
                5'd0:  fw = {4'd1, 12'd16, h_sw_id};
                5'd1:  fw = h_sw_ts[31:0];
                5'd2:  fw = h_sw_ts[63:32];
                5'd3:  fw = {30'd0, h_sw_ref};
                5'd4:  fw = sx18(h_sw_feat[54*0+36 +: 18]);
                5'd5:  fw = sx18(h_sw_feat[54*0+9  +: 18]);
                5'd6:  fw = {7'd0, h_sw_feat[54*0 +: 9], 7'd0, h_sw_feat[54*0+27 +: 9]};
                5'd7:  fw = sx18(h_sw_feat[54*1+36 +: 18]);
                5'd8:  fw = sx18(h_sw_feat[54*1+9  +: 18]);
                5'd9:  fw = {7'd0, h_sw_feat[54*1 +: 9], 7'd0, h_sw_feat[54*1+27 +: 9]};
                5'd10: fw = sx18(h_sw_feat[54*2+36 +: 18]);
                5'd11: fw = sx18(h_sw_feat[54*2+9  +: 18]);
                5'd12: fw = {7'd0, h_sw_feat[54*2 +: 9], 7'd0, h_sw_feat[54*2+27 +: 9]};
                5'd13: fw = sx18(h_sw_feat[54*3+36 +: 18]);
                5'd14: fw = sx18(h_sw_feat[54*3+9  +: 18]);
                default: fw = {7'd0, h_sw_feat[54*3 +: 9], 7'd0, h_sw_feat[54*3+27 +: 9]};
            endcase
            T_MB: case (k)
                5'd0:  fw = {4'd2, 12'd5, h_mb_id};
                5'd1:  fw = h_mb_ts[31:0];
                5'd2:  fw = h_mb_ts[63:32];
                5'd3:  fw = h_mb_idx;
                default: fw = {23'd0, h_mb_rel, 6'd0, h_mb_ref};
            endcase
            T_JOB: fw = (k == 5'd0) ? {4'd3, 12'd2, h_job_id} : h_job_w;
            T_LATE: fw = {4'd4, 12'd1, cur_late};
            T_RES: case (k)
                5'd0:  fw = {4'd5, mlen, h_res_id};
                5'd1:  fw = {h_res_type, h_res_err, 7'd0, h_res_L, 16'd0};
                5'd2:  fw = {7'd0, h_res_exp};
                5'd3:  fw = h_res_lag[31:0];
                5'd4:  fw = {16'd0, h_res_lag[47:32]};
                5'd5:  fw = sx18(h_res_pk[18*0 +: 18]);
                5'd6:  fw = sx18(h_res_pk[18*1 +: 18]);
                5'd7:  fw = sx18(h_res_pk[18*2 +: 18]);
                5'd8:  fw = sx18(h_res_pk[18*3 +: 18]);
                5'd9:  fw = sx18(h_res_pk[18*4 +: 18]);
                5'd10: fw = sx18(h_res_pk[18*5 +: 18]);
                5'd11: fw = {14'd0, h_res_side[18*0 +: 18]};
                5'd12: fw = {14'd0, h_res_side[18*1 +: 18]};
                5'd13: fw = {14'd0, h_res_side[18*2 +: 18]};
                5'd14: fw = {14'd0, h_res_side[18*3 +: 18]};
                5'd15: fw = {14'd0, h_res_side[18*4 +: 18]};
                default: fw = {14'd0, h_res_side[18*5 +: 18]};
            endcase
            default: fw = 32'd0;
        endcase
    end

    always @* lwa = m;

    // next message to send (fixed priority) and its length
    reg [2:0]  nsel;
    reg [11:0] nlen;
    always @* begin
        nsel = 3'd0; nlen = 12'd0;
        if      (!lq_empty) begin nsel = T_LATE; nlen = 12'd1;  end
        else if (p_res)     begin nsel = T_RES;  nlen = res_len_of_h; end
        else if (p_job)     begin nsel = T_JOB;  nlen = 12'd2;  end
        else if (p_sw)      begin nsel = T_SW;   nlen = 12'd16; end
        else if (p_mb)      begin nsel = T_MB;   nlen = 12'd5;  end
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            p_sw <= 1'b0; p_mb <= 1'b0; p_job <= 1'b0; p_res <= 1'b0;
            lq_w <= 3'd0; lq_r <= 3'd0;
            st <= S_IDLE; rwe <= 1'b0; wr_ptr <= {AW{1'b0}};
            drop_cnt <= 32'd0; drop_sticky <= 1'b0;
        end else if (clr) begin
            p_sw <= 1'b0; p_mb <= 1'b0; p_job <= 1'b0; p_res <= 1'b0;
            lq_w <= 3'd0; lq_r <= 3'd0;
            st <= S_IDLE; rwe <= 1'b0; wr_ptr <= {AW{1'b0}};
            drop_cnt <= 32'd0; drop_sticky <= 1'b0;
        end else begin
            rwe <= 1'b0;

            // ---- capture events (a busy source drops and counts) ----
            if (sw_v) begin
                if (p_sw) begin drop_cnt <= drop_cnt + 1'b1; drop_sticky <= 1'b1; end
                else begin p_sw <= 1'b1; h_sw_id <= sw_id; h_sw_ts <= sw_ts; h_sw_ref <= sw_ref; h_sw_feat <= sw_feat; end
            end
            if (mb_v) begin
                if (p_mb) begin drop_cnt <= drop_cnt + 1'b1; drop_sticky <= 1'b1; end
                else begin p_mb <= 1'b1; h_mb_id <= mb_id; h_mb_ts <= mb_ts; h_mb_idx <= mb_idx; h_mb_ref <= mb_ref; h_mb_rel <= mb_rel; end
            end
            if (job_v) begin
                if (p_job) begin drop_cnt <= drop_cnt + 1'b1; drop_sticky <= 1'b1; end
                else begin p_job <= 1'b1; h_job_id <= job_id; h_job_w <= {job_type, 11'd0, job_men, 4'd0, job_start}; end
            end
            if (res_v) begin
                if (p_res) begin drop_cnt <= drop_cnt + 1'b1; drop_sticky <= 1'b1; end
                else begin
                    p_res <= 1'b1; h_res_id <= res_id; h_res_type <= res_type; h_res_err <= res_err;
                    h_res_lag <= res_lag; h_res_pk <= res_pk; h_res_side <= res_side;
                    h_res_exp <= {res_exp_r, res_exp_z}; h_res_L <= Lsel;
                end
            end
            if (late_v) begin
                if (lq_full) begin drop_cnt <= drop_cnt + 1'b1; drop_sticky <= 1'b1; end
                else begin lq[lq_w[1:0]] <= late_id; lq_w <= lq_w + 1'b1; end
            end

            // ---- packer ----
            case (st)
                S_IDLE: if (nsel != 3'd0) begin
                    if ({{(AW+1){1'b0}}, nlen} > {12'd0, free}) begin
                        // no space: drop the whole message
                        drop_cnt <= drop_cnt + 1'b1; drop_sticky <= 1'b1;
                        case (nsel)
                            T_LATE: lq_r <= lq_r + 1'b1;
                            T_RES:  p_res <= 1'b0;
                            T_JOB:  p_job <= 1'b0;
                            T_SW:   p_sw  <= 1'b0;
                            default: p_mb <= 1'b0;
                        endcase
                    end else begin
                        sel  <= nsel;
                        mlen <= nlen;
                        k    <= 5'd0;
                        wa   <= wr_ptr;
                        nfix <= (nsel == T_SW) ? 5'd16 : (nsel == T_MB) ? 5'd5 : (nsel == T_JOB) ? 5'd2 :
                                (nsel == T_LATE) ? 5'd1 : 5'd17;
                        if (nsel == T_LATE) cur_late <= lq[lq_r[1:0]];
                        st   <= S_FIX;
                    end
                end

                S_FIX: begin                           // one fixed word per cycle
                    rwe <= 1'b1; rwa <= wa; rwd <= fw;
                    wa  <= wa + 1'b1;
                    k   <= k + 1'b1;
                    if (k == nfix - 1'b1) begin
                        if (sel == T_RES) begin
                            m     <= 8'd0;
                            mlast <= {h_res_L, 1'b0};  // 2L
                            st    <= S_LWRD;
                        end else
                            st <= S_DONE;
                    end
                end

                S_LWRD: begin                          // lwram[m] read issued (lwa = m)
                    pp <= 3'd0;
                    st <= S_LWWR;
                end

                S_LWWR: begin                          // 6 words of row m
                    rwe <= 1'b1; rwa <= wa; rwd <= sx18(lwq[18*pp +: 18]);
                    wa  <= wa + 1'b1;
                    if (pp == 3'd5) begin
                        if (m == mlast) st <= S_DONE;
                        else begin m <= m + 1'b1; st <= S_LWRD; end
                    end else
                        pp <= pp + 1'b1;
                end

                S_DONE: begin                          // publish and release the source
                    wr_ptr <= wa;
                    case (sel)
                        T_LATE: lq_r <= lq_r + 1'b1;
                        T_RES:  p_res <= 1'b0;
                        T_JOB:  p_job <= 1'b0;
                        T_SW:   p_sw  <= 1'b0;
                        default: p_mb <= 1'b0;
                    endcase
                    st <= S_IDLE;
                end

                default: st <= S_IDLE;
            endcase
        end
    end

    integer ii;
    initial for (ii = 0; ii < N; ii = ii + 1) ring[ii] = 32'd0;
endmodule
