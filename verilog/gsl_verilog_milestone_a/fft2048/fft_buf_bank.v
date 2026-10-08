`timescale 1ns/1ps

// ============================================================
//  fft_buf_bank.v   (v2: per-buffer client ports)
//
//  NBUF complex buffers (default 5: Z0, Z1, W0, W1, W2), each
//  2048 x 36 bit ({re[17:0], im[17:0]}).
//
//  Each buffer is split into two banks of 1024 words:
//      bank(addr) = ^addr          (XOR of all 11 address bits)
//      row(addr)  = addr[10:1]
//  The two operands of every radix-2 butterfly differ in exactly
//  one address bit, so they always sit in different banks. Each
//  bank is a simple dual-port RAM (1R + 1W), so the engine gets
//  2 reads + 2 writes per cycle with no conflicts.
//
//  Ports
//    engine port : bank-level (engine does the bank mapping).
//                  Owns buffer eng_sel while eng_active = 1.
//    client ports: one per buffer, word-level, 1 read + 1 write
//                  per cycle, 1-cycle read latency. Different
//                  buffers can be accessed in the same cycle, so
//                  xspec_phat can read Z0 and Z1 while writing
//                  W0/W1/W2. Who drives each client port
//                  (gcc_loader, xspec_phat, peak_search, host) is
//                  muxed outside this module (gcc_engine).
//
//  Flattened buses: buffer i uses bits [i*W +: W].
//  A client access to the buffer the engine owns is not allowed;
//  cl_conflict flags it (one cycle later).
// ============================================================
module fft_buf_bank #(
    parameter NBUF = 5,
    parameter DW   = 36,
    parameter AW   = 11
)(
    input  wire                 clk,

    // ---- engine port (bank level) ----
    input  wire                 eng_active,
    input  wire [2:0]           eng_sel,
    input  wire [AW-2:0]        eng_raddr0,
    input  wire [AW-2:0]        eng_raddr1,
    output wire [DW-1:0]        eng_rdata0,
    output wire [DW-1:0]        eng_rdata1,
    input  wire                 eng_we0,
    input  wire                 eng_we1,
    input  wire [AW-2:0]        eng_waddr0,
    input  wire [AW-2:0]        eng_waddr1,
    input  wire [DW-1:0]        eng_wdata0,
    input  wire [DW-1:0]        eng_wdata1,

    // ---- client ports (word level, one per buffer) ----
    input  wire [NBUF-1:0]      cl_re,        // read strobe (conflict check only)
    input  wire [NBUF*AW-1:0]   cl_raddr,
    output wire [NBUF*DW-1:0]   cl_rdata,     // data for cl_raddr of previous cycle
    input  wire [NBUF-1:0]      cl_we,
    input  wire [NBUF*AW-1:0]   cl_waddr,
    input  wire [NBUF*DW-1:0]   cl_wdata,
    output reg                  cl_conflict
);
    localparam RW = AW - 1;

    wire [NBUF*DW-1:0] rd0_flat;
    wire [NBUF*DW-1:0] rd1_flat;
    wire [NBUF-1:0]    own_vec;

    genvar i;
    generate
        for (i = 0; i < NBUF; i = i + 1) begin : g_buf
            wire          own    = eng_active && (eng_sel == i);
            wire [AW-1:0] c_ra   = cl_raddr[i*AW +: AW];
            wire [AW-1:0] c_wa   = cl_waddr[i*AW +: AW];
            wire [DW-1:0] c_wd   = cl_wdata[i*DW +: DW];
            wire          c_wb   = ^c_wa;               // bank of client write
            reg           c_rb_d;                       // bank of client read, delayed

            assign own_vec[i] = own;

            fft_sdp_ram #(.DW(DW), .AW(RW)) u_bank0 (
                .clk   (clk),
                .we    (own ? eng_we0    : (cl_we[i] && !c_wb)),
                .waddr (own ? eng_waddr0 : c_wa[AW-1:1]),
                .wdata (own ? eng_wdata0 : c_wd),
                .raddr (own ? eng_raddr0 : c_ra[AW-1:1]),
                .rdata (rd0_flat[i*DW +: DW])
            );

            fft_sdp_ram #(.DW(DW), .AW(RW)) u_bank1 (
                .clk   (clk),
                .we    (own ? eng_we1    : (cl_we[i] &&  c_wb)),
                .waddr (own ? eng_waddr1 : c_wa[AW-1:1]),
                .wdata (own ? eng_wdata1 : c_wd),
                .raddr (own ? eng_raddr1 : c_ra[AW-1:1]),
                .rdata (rd1_flat[i*DW +: DW])
            );

            always @(posedge clk)
                c_rb_d <= ^c_ra;

            assign cl_rdata[i*DW +: DW] = c_rb_d ? rd1_flat[i*DW +: DW]
                                                 : rd0_flat[i*DW +: DW];
        end
    endgenerate

    // Engine read data (eng_sel is constant for a whole transform)
    assign eng_rdata0 = rd0_flat[eng_sel*DW +: DW];
    assign eng_rdata1 = rd1_flat[eng_sel*DW +: DW];

    always @(posedge clk)
        cl_conflict <= |((cl_re | cl_we) & own_vec);
endmodule
