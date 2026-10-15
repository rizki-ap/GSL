// ============================================================
// gsl_replay_dma — streams a recording from SDRAM into replay_src (C1)
//
//  Frames are 64-bit little-endian words at base + 8*n: the 32-bit word at
//  base + 8n holds ch2/ch3, the word at base + 8n + 4 holds ch0/ch1, so the
//  frame is {word1, word0} = {ch0, ch1, ch2, ch3} as in mb_*.hex.
//
//  Avalon-MM read master, 32-bit, pipelined (readdatavalid), at most
//  MAXOUT reads outstanding. Keeps (FIFO level + frames in flight) below
//  FILL frames. Non-looping: stops after nframes and raises 'done' once
//  every frame has been pushed. Looping: wraps to frame 0, counts passes.
//  'primed' rises when 256 frames are buffered (or all of a short
//  recording), like the Milestone B testbench, so the frame cadence never
//  starts on an empty FIFO.
// ============================================================
module gsl_replay_dma #(
    parameter integer MAXOUT = 16,
    parameter integer FILL   = 1000,
    parameter integer PRIME  = 256
)(
    input  wire          clk,
    input  wire          rst_n,
    input  wire          clr,
    input  wire          enable,           // replay running
    input  wire          loop,
    input  wire [31:0]   base,
    input  wire [31:0]   nframes,

    output reg  [31:0]   avm_address,
    output reg           avm_read,
    input  wire          avm_waitrequest,
    input  wire [31:0]   avm_readdata,
    input  wire          avm_readdatavalid,

    input  wire [10:0]   fifo_level,
    output reg           fifo_we,
    output reg  [63:0]   fifo_wdata,

    output reg  [31:0]   pos,              // frames pushed in this pass
    output reg  [31:0]   loops,            // completed passes
    output reg           done,             // non-looping: all frames pushed
    output reg           primed
);
    reg  [32:0] rq;                         // next word to request (0 .. 2*nframes-1)
    reg  [4:0]  outst;                      // reads outstanding
    reg         half;                       // low word of a frame received
    reg  [31:0] lo;
    reg  [31:0] pushed;                     // frames pushed in total (this run)

    wire [32:0] nwords   = {nframes, 1'b0};
    wire        all_req  = !loop && (rq == nwords);
    wire        accept   = avm_read && !avm_waitrequest;
    wire        inflight_frames_ok = ({21'd0, fifo_level} + {27'd0, outst[4:1]} + 32'd2) < FILL;
    wire        can_req  = enable && (nframes != 32'd0) && !all_req && (outst < MAXOUT - 1) && inflight_frames_ok;
    wire [32:0] rq_next  = (loop && (rq + 33'd1 == nwords)) ? 33'd0 : rq + 33'd1;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            avm_read <= 1'b0; avm_address <= 32'd0;
            rq <= 33'd0; outst <= 5'd0; half <= 1'b0;
            fifo_we <= 1'b0; pos <= 32'd0; loops <= 32'd0; pushed <= 32'd0;
            done <= 1'b0; primed <= 1'b0;
        end else if (clr) begin
            avm_read <= 1'b0;
            rq <= 33'd0; outst <= 5'd0; half <= 1'b0;
            fifo_we <= 1'b0; pos <= 32'd0; loops <= 32'd0; pushed <= 32'd0;
            done <= 1'b0; primed <= 1'b0;
        end else begin
            fifo_we <= 1'b0;

            // ---- request side (Avalon: hold address/read while waitrequest) ----
            if (accept) rq <= rq_next;
            if (!avm_read || accept) begin
                // after an accepted transfer the next request uses the updated index
                if (can_req && !(accept && !loop && rq_next == nwords)) begin
                    avm_read    <= 1'b1;
                    avm_address <= base + {(accept ? rq_next[29:0] : rq[29:0]), 2'b00};
                end else
                    avm_read <= 1'b0;
            end
            outst <= outst + (accept ? 5'd1 : 5'd0) - (avm_readdatavalid ? 5'd1 : 5'd0);

            // ---- response side ----
            if (avm_readdatavalid) begin
                if (!half) begin
                    lo   <= avm_readdata;
                    half <= 1'b1;
                end else begin
                    half       <= 1'b0;
                    fifo_we    <= 1'b1;
                    fifo_wdata <= {avm_readdata, lo};
                    pushed     <= pushed + 1'b1;
                    if (loop && pos + 32'd1 == nframes) begin
                        pos   <= 32'd0;
                        loops <= loops + 1'b1;
                    end else
                        pos <= pos + 1'b1;
                end
            end

            if (!loop && nframes != 32'd0 && pushed == nframes) done <= 1'b1;
            if (fifo_level >= PRIME || done) primed <= 1'b1;
        end
    end
endmodule
