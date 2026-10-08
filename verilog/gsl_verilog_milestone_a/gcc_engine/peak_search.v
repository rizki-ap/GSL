`timescale 1ns/1ps

// ============================================================
//  peak_search.v
//  Coarse TDOA peak for the 6 GCC-PHAT pairs, plus the lag window
//  the processor (HPS / Nios) uses for sub-sample refinement.
//
//  Spec: verilog/fpga_blocks.md section 3.5
//
//  Input: the three IFFT outputs, natural order. Each buffer holds
//  two pairs: W0 = r01 + j*r02, W1 = r03 + j*r12, W2 = r13 + j*r23.
//  All three are read at the same address every cycle (different
//  buffers), so one read gives all 6 pairs at one lag.
//
//  Pass 1  (lags -L .. L): per-pair argmax (first maximum wins,
//          same as numpy argmax in gs_det_tdoa.py), and stream every
//          lag to the lag-window output (all 6 pairs in one word).
//  Pass 2  (lags -L .. L): per-pair sidelobe = max |r| over lags
//          more than 2 samples from that pair's peak.
//
//  Peak-to-sidelobe ratio = pk_val / pk_side is computed by the
//  processor (avoids a hardware divider). pk_val <= 0 or a low
//  ratio flags an unreliable pair (e.g. dead microphone).
//
//  Cycles: 2 * (2L + 1) + 4  (L = 100: 406;  L = 20: 86)
// ============================================================
module peak_search #(
    parameter N     = 2048,
    parameter LOG2N = 11,
    parameter DW    = 18
)(
    input  wire                 clk,
    input  wire                 rst_n,

    input  wire                 ps_start,
    input  wire [6:0]           ps_maxlag,     // L, 1 .. 127
    output reg                  ps_busy,
    output reg                  ps_done,

    output wire                 w_re,
    output wire [LOG2N-1:0]     w_raddr,       // same address to W0, W1, W2
    input  wire [2*DW-1:0]      w0_rdata,
    input  wire [2*DW-1:0]      w1_rdata,
    input  wire [2*DW-1:0]      w2_rdata,

    output reg  [6*8-1:0]       pk_lag,        // s8 per pair, pair j at [j*8 +: 8]
    output reg  [6*DW-1:0]      pk_val,        // s18 per pair
    output reg  [6*DW-1:0]      pk_side,       // u18 per pair (max |r| off-peak)

    output reg                  lw_we,
    output reg  [7:0]           lw_idx,        // lag + L, 0 .. 2L
    output reg  [6*DW-1:0]      lw_data        // pair j at [j*DW +: DW]
);
    localparam [1:0] S_IDLE = 2'd0, S_P1 = 2'd1, S_P2 = 2'd2, S_END = 2'd3;

    reg [1:0]         state;
    reg signed [8:0]  m;            // current lag being issued
    reg signed [8:0]  L;
    reg               iss;          // a read is issued this cycle

    assign w_re    = iss;
    assign w_raddr = {{(LOG2N-9){m[8]}}, m};   // m mod N (sign-extend, wraps)

    // read data arrives 1 cycle later
    reg               d_v, d_pass2, d_first;
    reg signed [8:0]  d_m;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) d_v <= 1'b0;
        else        d_v <= iss;
    end
    always @(posedge clk) begin
        d_m     <= m;
        d_pass2 <= (state == S_P2);
        d_first <= (m == -L);
    end

    // the 6 pair values at lag d_m
    wire signed [DW-1:0] v [0:5];
    assign v[0] = w0_rdata[2*DW-1:DW];
    assign v[1] = w0_rdata[DW-1:0];
    assign v[2] = w1_rdata[2*DW-1:DW];
    assign v[3] = w1_rdata[DW-1:0];
    assign v[4] = w2_rdata[2*DW-1:DW];
    assign v[5] = w2_rdata[DW-1:0];

    reg signed [DW-1:0] best [0:5];
    reg signed [8:0]    blag [0:5];
    reg        [DW-1:0] side [0:5];

    function [DW-1:0] absv;                 // |v| as unsigned (fits 2^17)
        input signed [DW-1:0] x;
        absv = x[DW-1] ? -x : x;
    endfunction

    integer j;
    always @(posedge clk) begin
        lw_we <= 1'b0;
        if (d_v && !d_pass2) begin
            // pass 1: argmax, lag-window stream
            for (j = 0; j < 6; j = j + 1)
                if (d_first || v[j] > best[j]) begin
                    best[j] <= v[j];
                    blag[j] <= d_m;
                end
            lw_we   <= 1'b1;
            lw_idx  <= d_m + L;
            lw_data <= {v[5], v[4], v[3], v[2], v[1], v[0]};
        end
        if (d_v && d_pass2) begin
            // pass 2: sidelobe max outside +/-2 of each peak
            for (j = 0; j < 6; j = j + 1) begin
                if (d_first) side[j] <= {DW{1'b0}};
                if ((d_m - blag[j] > 2) || (blag[j] - d_m > 2))
                    if (d_first || absv(v[j]) > side[j])
                        side[j] <= absv(v[j]);
            end
        end
    end

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state   <= S_IDLE;
            m       <= 9'sd0;
            L       <= 9'sd0;
            iss     <= 1'b0;
            ps_busy <= 1'b0;
            ps_done <= 1'b0;
        end else begin
            ps_done <= 1'b0;
            case (state)
                S_IDLE: if (ps_start) begin
                    L       <= {2'b00, ps_maxlag};
                    m       <= -$signed({2'b00, ps_maxlag});
                    iss     <= 1'b1;
                    ps_busy <= 1'b1;
                    state   <= S_P1;
                end
                S_P1: begin
                    if (m == L) begin
                        m     <= -L;
                        state <= S_P2;      // iss stays 1: pass 2 starts next cycle
                    end else
                        m <= m + 9'sd1;
                end
                S_P2: begin
                    if (m == L) begin
                        iss   <= 1'b0;
                        state <= S_END;
                    end else
                        m <= m + 9'sd1;
                end
                S_END: if (!d_v) begin
                    ps_busy <= 1'b0;
                    ps_done <= 1'b1;
                    state   <= S_IDLE;
                end
            endcase
        end
    end

    // register results when pass 2 is complete
    always @(posedge clk)
        if (state == S_END && !d_v)
            for (j = 0; j < 6; j = j + 1) begin
                pk_lag[j*8 +: 8]   <= blag[j][7:0];
                pk_val[j*DW +: DW] <= best[j];
                pk_side[j*DW +: DW]<= side[j];
            end
endmodule
