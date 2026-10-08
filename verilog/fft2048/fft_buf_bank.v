// ============================================================
//  fft_buf_bank.v
//  NBUF complex buffers (default 5: Z0, Z1, W0, W1, W2), each
//  2048 x 36 bit ({re[17:0], im[17:0]}).
//
//  Each buffer is split into two banks of 1024 words:
//      bank(addr) = ^addr          (XOR of all 11 address bits)
//      row(addr)  = addr[10:1]
//  The two operands of every radix-2 butterfly differ in exactly
//  one address bit, so they always sit in different banks. With
//  each bank a simple dual-port RAM (1R + 1W), the engine gets
//  2 reads + 2 writes per cycle with no conflicts.
//
//  Ports
//    engine port : bank-level (engine does the bank mapping),
//                  owns buffer eng_sel while eng_active = 1
//    ext port    : word-level, one access per cycle, for the
//                  producers / consumers (gcc_loader, xspec_phat,
//                  peak_search, testbench). 1-cycle read latency.
//
//  Accessing the buffer the engine currently owns through the
//  ext port is not allowed; ext_conflict flags it.
// ============================================================
`timescale 1ns/1ps

module fft_buf_bank #(
    parameter NBUF = 5,
    parameter DW   = 36,
    parameter AW   = 11
)(
    input  wire              clk,

    // ---- engine port (bank level) ----
    input  wire              eng_active,
    input  wire [2:0]        eng_sel,
    input  wire [AW-2:0]     eng_raddr0,
    input  wire [AW-2:0]     eng_raddr1,
    output wire [DW-1:0]     eng_rdata0,
    output wire [DW-1:0]     eng_rdata1,
    input  wire              eng_we0,
    input  wire              eng_we1,
    input  wire [AW-2:0]     eng_waddr0,
    input  wire [AW-2:0]     eng_waddr1,
    input  wire [DW-1:0]     eng_wdata0,
    input  wire [DW-1:0]     eng_wdata1,

    // ---- external port (word level) ----
    input  wire              ext_en,
    input  wire              ext_we,
    input  wire [2:0]        ext_sel,
    input  wire [AW-1:0]     ext_addr,
    input  wire [DW-1:0]     ext_wdata,
    output wire [DW-1:0]     ext_rdata,
    output reg               ext_conflict
);
    localparam RW = AW - 1;

    wire          ext_bank = ^ext_addr;
    wire [RW-1:0] ext_row  = ext_addr[AW-1:1];

    wire [NBUF*DW-1:0] rd0_flat;
    wire [NBUF*DW-1:0] rd1_flat;

    genvar i;
    generate
        for (i = 0; i < NBUF; i = i + 1) begin : g_buf
            wire own     = eng_active && (eng_sel == i);
            wire ext_hit = ext_en && (ext_sel == i);

            fft_sdp_ram #(.DW(DW), .AW(RW)) u_bank0 (
                .clk   (clk),
                .we    (own ? eng_we0    : (ext_hit && ext_we && !ext_bank)),
                .waddr (own ? eng_waddr0 : ext_row),
                .wdata (own ? eng_wdata0 : ext_wdata),
                .raddr (own ? eng_raddr0 : ext_row),
                .rdata (rd0_flat[i*DW +: DW])
            );

            fft_sdp_ram #(.DW(DW), .AW(RW)) u_bank1 (
                .clk   (clk),
                .we    (own ? eng_we1    : (ext_hit && ext_we &&  ext_bank)),
                .waddr (own ? eng_waddr1 : ext_row),
                .wdata (own ? eng_wdata1 : ext_wdata),
                .raddr (own ? eng_raddr1 : ext_row),
                .rdata (rd1_flat[i*DW +: DW])
            );
        end
    endgenerate

    // Engine read data (eng_sel is constant for a whole transform)
    assign eng_rdata0 = rd0_flat[eng_sel*DW +: DW];
    assign eng_rdata1 = rd1_flat[eng_sel*DW +: DW];

    // External read data: select by the address of the previous cycle
    reg [2:0] ext_sel_d;
    reg       ext_bank_d;
    always @(posedge clk) begin
        ext_sel_d    <= ext_sel;
        ext_bank_d   <= ext_bank;
        ext_conflict <= ext_en && eng_active && (ext_sel == eng_sel);
    end

    assign ext_rdata = ext_bank_d ? rd1_flat[ext_sel_d*DW +: DW]
                                  : rd0_flat[ext_sel_d*DW +: DW];
endmodule
