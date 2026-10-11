`timescale 1ns/1ps

// ============================================================
//  mb_decim.v
//  /5 decimator for the muzzle-blast path: 100 kHz -> 20 kHz.
//
//  Spec: verilog/fpga_blocks.md section 2.6
//
//  63-tap linear-phase low-pass (decim63.hex, Q1.17, DC gain 1):
//  passband 0-4 kHz ripple 0.02 dB, stopband >= 16 kHz below -71 dB.
//      y = sat18( round( sum_{k=0..62} h[k] * x[n-k], 17 ) )
//  An output is produced on every 5th input (input counter = 4 mod 5),
//  once 63 inputs have been seen (full history).
//
//  Timestamps
//    y_ts  = x_ts(newest) - 31  : the input index at the filter centre,
//            i.e. group delay already removed, so MB and SW timestamps
//            compare directly (Delta-t budget 0.5 ms).
//    y_idx = output count (20 kHz index); mb_ring address = y_idx[11:0].
//
//  Timing (area-optimised, v2): one multiplier-accumulator, channels in
//  sequence, one tap per cycle -> 4 x (63 + 2) + 3 = 263 cycles per output
//  (an output is due every 5 frames = 5000 cycles). History: 128 x 72-bit RAM.
//  Same summation order as v1: byte-identical results.
// ============================================================
module mb_decim #(
    parameter TAP_INIT = "decim63.hex"
)(
    input  wire          clk,
    input  wire          rst_n,

    input  wire          x_valid,
    input  wire [71:0]   x_data,         // {ch0..ch3} s18
    input  wire [63:0]   x_ts,

    output reg           y_valid,
    output reg  [71:0]   y_data,
    output reg  [63:0]   y_ts,
    output reg  [31:0]   y_idx           // index of the sample on y_valid
);
    localparam integer NT = 63;
    localparam integer GD = 31;

    function signed [17:0] sat18;
        input signed [47:0] v;
        if (v > 48'sd131071)       sat18 = 18'sd131071;
        else if (v < -48'sd131072) sat18 = -18'sd131072;
        else                       sat18 = v[17:0];
    endfunction

    function signed [47:0] rconv;
        input signed [47:0] v;
        input integer       s;
        reg   signed [48:0] tt;
        begin
            tt = v;
            tt = tt + ((49'sd1 <<< (s-1)) - 49'sd1) + ((v >>> s) & 49'sd1);
            rconv = tt >>> s;
        end
    endfunction

    (* romstyle = "M10K" *) reg [17:0] taps [0:NT-1];
    initial $readmemh(TAP_INIT, taps);

    // history RAM: write newest at wptr. 128 slots for 63 taps: new samples
    // keep arriving while an output is computed (263 cycles; 2.6 frames at
    // 100 clocks per frame), and the 65 spare slots keep them from
    // overwriting taps still to be read, at any frame cadence.
    reg [71:0] hram [0:127];
    reg [6:0]  wptr;
    reg [6:0]  raddr;
    reg [71:0] hq;
    reg [5:0]  taddr;
    reg signed [17:0] hk;
    always @(posedge clk) begin
        if (x_valid) hram[wptr] <= x_data;
        hq <= hram[raddr];
        hk <= taps[taddr];
    end

    reg [2:0]  ph;                  // input counter mod 5
    reg [6:0]  filled;              // inputs seen, saturates at 63
    reg        busy;
    reg [1:0]  ch;                  // channel being filtered
    reg [6:0]  k;                   // tap being issued
    reg [1:0]  pipe;                // read-latency pipeline (valid bits)
    reg [63:0] ts_new;
    reg signed [42:0] acc;
    reg [6:0]  newest;

    wire signed [17:0] hx   = hq[71 - 18*ch -: 18];
    wire signed [42:0] last = acc + hk * hx;          // includes the tap now arriving

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            wptr <= 7'd0; ph <= 3'd0; filled <= 7'd0;
            busy <= 1'b0; ch <= 2'd0; k <= 7'd0; pipe <= 2'b00;
            y_valid <= 1'b0; y_idx <= 32'hFFFF_FFFF;
        end else begin
            y_valid <= 1'b0;
            if (x_valid) begin
                wptr   <= wptr + 1'b1;
                if (filled != NT) filled <= filled + 1'b1;
                ph     <= (ph == 3'd4) ? 3'd0 : ph + 1'b1;
                if (ph == 3'd4 && (filled + 1'b1) >= NT) begin
                    busy   <= 1'b1;
                    ch     <= 2'd0;
                    k      <= 7'd0;
                    newest <= wptr;
                    ts_new <= x_ts;
                    acc    <= 43'sd0;
                end
            end
            if (busy) begin
                // issue tap k of channel ch: history (newest - k), coefficient k
                if (k < NT) begin
                    raddr <= newest - k;
                    taddr <= k[5:0];
                    k     <= k + 1'b1;
                end
                pipe <= {pipe[0], (k < NT)};
                if (pipe[1])
                    acc <= last;
                if (k == NT && pipe == 2'b10) begin       // last tap of this channel arriving
                    y_data[71 - 18*ch -: 18] <= sat18(rconv(last, 17));
                    acc  <= 43'sd0;
                    pipe <= 2'b00;
                    if (ch == 2'd3) begin
                        busy    <= 1'b0;
                        y_valid <= 1'b1;
                        y_idx   <= y_idx + 1'b1;
                        y_ts    <= ts_new - GD;
                    end else begin
                        ch <= ch + 1'b1;
                        k  <= 7'd0;
                    end
                end
            end
        end
    end

    integer jj;
    initial for (jj = 0; jj < 128; jj = jj + 1) hram[jj] = 72'd0;
endmodule
