// ============================================================
// tb_gsl_mm — Milestone C1 testbench
//
//  An SDRAM model holds the recording and answers gsl_mm's replay DMA
//  (random wait states, pipelined read latency). A processor model talks
//  to gsl_mm only over the Avalon-MM slave: it writes the 43 configuration
//  registers, starts the replay, polls the mailbox, and rebuilds every
//  message into the same log lines as the Milestone B testbench
//  (S / M / J / X / L / R), so the logs can be compared with Milestone B.
//  The x and y streams are tapped hierarchically to check the replay path.
//
//  Plusargs: +stim=<hex> +nframes=<n> +cfg=<file> +out=<prefix>
//            +slow=<cycles>  extra delay between mailbox word reads
//            +loop=<passes>  loop mode: stop after this many passes
//            +raw=<file>     also write every mailbox word (hex), as the board
//                            bring-up program logs them (tests board_check.py)
// ============================================================
`timescale 1ns/1ps
module tb_gsl_mm;
    parameter FRAME_DIV = 250;
    localparam integer MAXF = 524288;

    reg clk = 1'b0;
    always #5 clk = ~clk;
    reg rst_n = 1'b0;

    // ---------------- DUT ----------------
    reg  [13:0] avs_address = 14'd0;
    reg         avs_read = 1'b0, avs_write = 1'b0;
    reg  [31:0] avs_writedata = 32'd0;
    wire [31:0] avs_readdata;
    wire        avs_readdatavalid;
    wire [31:0] avm_address;
    wire        avm_read;
    reg         avm_waitrequest;
    reg  [31:0] avm_readdata;
    reg         avm_readdatavalid;
    wire        irq;

    gsl_mm #(
        .FRAME_DIV(FRAME_DIV),
        .FD_INIT("../stream/fe_fdelay.hex"), .DEC_INIT("../stream/decim63.hex"),
        .TW_INIT_RE("../fft2048/tw2048_re.hex"), .TW_INIT_IM("../fft2048/tw2048_im.hex"),
        .RSQRT_INIT("../xspec_phat/rsqrt1024.hex")
    ) dut (
        .clk(clk), .rst_n(rst_n),
        .avs_address(avs_address), .avs_read(avs_read), .avs_write(avs_write),
        .avs_writedata(avs_writedata), .avs_readdata(avs_readdata), .avs_readdatavalid(avs_readdatavalid),
        .avm_address(avm_address), .avm_read(avm_read), .avm_waitrequest(avm_waitrequest),
        .avm_readdata(avm_readdata), .avm_readdatavalid(avm_readdatavalid),
        .irq(irq), .adc_valid(1'b0), .adc_data(64'd0)
    );

    // ---------------- SDRAM model ----------------
    reg [63:0] stim [0:MAXF-1];
    localparam integer LAT = 4;
    reg [31:0] lfsr = 32'h1234_5678;
    reg        pv [0:LAT-1];
    reg [31:0] pd [0:LAT-1];
    integer    q;
    always @(posedge clk) lfsr <= {lfsr[30:0], lfsr[31] ^ lfsr[21] ^ lfsr[1] ^ lfsr[0]};
    always @* avm_waitrequest = (lfsr[1:0] == 2'b00);         // ~25% wait states
    always @(posedge clk) begin
        for (q = LAT - 1; q > 0; q = q - 1) begin pv[q] <= pv[q-1]; pd[q] <= pd[q-1]; end
        pv[0] <= avm_read && !avm_waitrequest;
        pd[0] <= avm_address[2] ? stim[avm_address[31:3]][63:32] : stim[avm_address[31:3]][31:0];
        avm_readdatavalid <= pv[LAT-1];
        avm_readdata      <= pd[LAT-1];
    end
    initial for (q = 0; q < LAT; q = q + 1) begin pv[q] = 1'b0; pd[q] = 32'd0; end

    // ---------------- stream taps ----------------
    integer fx, fy, fe, fr;
    function integer s18; input [17:0] v; s18 = $signed(v); endfunction
    always @(posedge clk) if (rst_n) begin
        if (dut.u_core.x_valid)
            $fdisplay(fx, "%0d %0d %0d %0d %0d", $signed(dut.u_core.x_ts),
                      s18(dut.u_core.x_data[71:54]), s18(dut.u_core.x_data[53:36]),
                      s18(dut.u_core.x_data[35:18]), s18(dut.u_core.x_data[17:0]));
        if (dut.u_core.y_valid)
            $fdisplay(fy, "%0d %0d %0d %0d %0d %0d", $signed(dut.u_core.y_ts), dut.u_core.y_idx,
                      s18(dut.u_core.y_data[71:54]), s18(dut.u_core.y_data[53:36]),
                      s18(dut.u_core.y_data[35:18]), s18(dut.u_core.y_data[17:0]));
    end

    // ---------------- bus tasks ----------------
    task bus_write(input [13:0] a, input [31:0] d);
        begin
            @(posedge clk);
            avs_address <= a; avs_writedata <= d; avs_write <= 1'b1;
            @(posedge clk);
            avs_write <= 1'b0;
        end
    endtask

    task bus_read(input [13:0] a, output [31:0] d);
        begin
            @(posedge clk);
            avs_address <= a; avs_read <= 1'b1;
            @(posedge clk);
            avs_read <= 1'b0;
            @(posedge clk);
            while (!avs_readdatavalid) @(posedge clk);
            d = avs_readdata;
        end
    endtask

    // ---------------- processor model ----------------
    reg [8*256:1] stim_file, cfg_file, out_pfx, fname;
    integer nframes, slow, loopn, cf, r, k, v, idle, len, typ, id, j, m, L, nmsg;
    integer n_sw, n_mb, n_job, n_late, n_res;
    reg [31:0] w, lvl, st, rd, msg [0:4095];
    reg [63:0] cycles = 64'd0;
    always @(posedge clk) cycles <= cycles + 1'b1;

    function integer s8; input [7:0] x; s8 = $signed(x); endfunction
    function integer s32; input [31:0] x; s32 = $signed(x); endfunction

    task handle_msg;
        begin
            typ = msg[0][31:28]; len = msg[0][27:16]; id = msg[0][15:0];
            nmsg = nmsg + 1;
            case (typ)
                1: begin
                    n_sw = n_sw + 1;
                    $fwrite(fe, "S %0d %0d %0d", id, $signed({msg[2], msg[1]}), msg[3]);
                    for (j = 0; j < 4; j = j + 1)
                        $fwrite(fe, " %0d %0d %0d %0d", s32(msg[4+3*j]), msg[6+3*j][8:0],
                                s32(msg[5+3*j]), msg[6+3*j][24:16]);
                    $fwrite(fe, "\n");
                end
                2: begin
                    n_mb = n_mb + 1;
                    $fdisplay(fe, "M %0d %0d %0d %0d %0d", id, $signed({msg[2], msg[1]}), msg[3],
                              msg[4][1:0], msg[4][8]);
                end
                3: begin
                    n_job = n_job + 1;
                    $fdisplay(fe, "J %0d %0d %0d %0d", id, msg[1][31], msg[1][11:0], msg[1][19:16]);
                end
                4: begin
                    n_late = n_late + 1;
                    $fdisplay(fe, "X %0d", id);
                end
                5: begin
                    n_res = n_res + 1;
                    L = msg[1][22:16];
                    for (m = 0; m <= 2 * L; m = m + 1) begin
                        $fwrite(fe, "L %0d %0d", id, m);
                        for (j = 0; j < 6; j = j + 1) $fwrite(fe, " %0d", s32(msg[17 + 6*m + j]));
                        $fwrite(fe, "\n");
                    end
                    $fwrite(fe, "R %0d %0d %0d %0d %0d %0d %0d %0d", id, msg[1][31], msg[1][30],
                            msg[2][4:0], msg[2][9:5], msg[2][14:10], msg[2][19:15], msg[2][24:20]);
                    for (j = 0; j < 4; j = j + 1) $fwrite(fe, " %0d", s8(msg[3][8*j +: 8]));
                    for (j = 0; j < 2; j = j + 1) $fwrite(fe, " %0d", s8(msg[4][8*j +: 8]));
                    for (j = 0; j < 6; j = j + 1) $fwrite(fe, " %0d", s32(msg[5+j]));
                    for (j = 0; j < 6; j = j + 1) $fwrite(fe, " %0d", msg[11+j]);
                    $fwrite(fe, "\n");
                end
                default: $display("ERROR: unknown message type %0d", typ);
            endcase
        end
    endtask

    initial begin
        n_sw = 0; n_mb = 0; n_job = 0; n_late = 0; n_res = 0; nmsg = 0;
        if (!$value$plusargs("stim=%s", stim_file)) stim_file = "../stream/mb_s1.hex";
        if (!$value$plusargs("nframes=%d", nframes)) nframes = 0;
        if (!$value$plusargs("cfg=%s", cfg_file)) cfg_file = "../stream/mb_cfg.txt";
        if (!$value$plusargs("out=%s", out_pfx)) out_pfx = "c1";
        if (!$value$plusargs("slow=%d", slow)) slow = 0;
        if (!$value$plusargs("loop=%d", loopn)) loopn = 0;
        fr = 0;
        if ($value$plusargs("raw=%s", fname)) fr = $fopen(fname, "w");
        $readmemh(stim_file, stim, 0, nframes - 1);
        $sformat(fname, "%0s_x.txt",  out_pfx); fx = $fopen(fname, "w");
        $sformat(fname, "%0s_y.txt",  out_pfx); fy = $fopen(fname, "w");
        $sformat(fname, "%0s_ev.txt", out_pfx); fe = $fopen(fname, "w");

        repeat (5) @(posedge clk);
        rst_n <= 1'b1;
        repeat (100) @(posedge clk);                 // mb_det clears its RAM after reset

        bus_read(14'd0, w);
        if (w !== 32'h4753_4C31) $display("ERROR: ID register reads %h", w);
        // configuration
        cf = $fopen(cfg_file, "r");
        for (k = 0; k < 43; k = k + 1) begin
            r = $fscanf(cf, "%d", v);
            bus_write(14'd64 + k, v);
        end
        $fclose(cf);
        for (k = 0; k < 43; k = k + 1) begin          // read back (signed fields sign-extend)
            bus_read(14'd64 + k, w);
        end
        bus_write(14'd16, 32'd0);                     // RP_BASE
        bus_write(14'd17, nframes);                   // RP_NFRAMES
        bus_write(14'd2, (loopn > 0) ? 32'h13 : 32'h03);   // run | replay [| loop]

        rd = 0; idle = 0;
        while (1) begin
            bus_read(14'd14, lvl);                    // MBOX_LEVEL
            if (lvl != 0) begin
                bus_read(14'h2000 + rd[12:0], msg[0]);
                len = msg[0][27:16];
                for (k = 1; k < len; k = k + 1) begin
                    bus_read(14'h2000 + ((rd + k) & 32'h1FFF), msg[k]);
                    if (slow > 0) repeat (slow) @(posedge clk);
                end
                handle_msg;
                if (fr != 0) for (k = 0; k < len; k = k + 1) $fdisplay(fr, "%08h", msg[k]);
                rd = (rd + len) & 32'h1FFF;
                bus_write(14'd13, rd);                // MBOX_RD
                idle = 0;
            end else begin
                bus_read(14'd3, st);                  // STATUS
                if (loopn > 0) begin
                    bus_read(14'd19, w);              // RP_LOOPS
                    if (w >= loopn) bus_write(14'd2, 32'h02);   // stop after N passes
                end
                if ((st[6] || (loopn > 0 && w >= loopn)) && !st[5]) idle = idle + 1;
                else idle = 0;
                if (idle > 2000) begin
                    bus_read(14'd6, w);  $display("  late jobs        : %0d", w);
                    bus_read(14'd7, w);  $display("  pending overflow : %0d", w);
                    bus_read(14'd8, w);  $display("  mailbox drops    : %0d", w);
                    bus_read(14'd3, st); $display("  STATUS           : %h (underflow %0d, res_err %0d)", st, st[0], st[4]);
                    bus_read(14'd19, w); $display("  replay passes    : %0d", w);
                    $display("  messages %0d: SW %0d, MB %0d, jobs %0d, late %0d, results %0d",
                             nmsg, n_sw, n_mb, n_job, n_late, n_res);
                    $display("  E %0d cycles", cycles);
                    $fclose(fx); $fclose(fy); $fclose(fe);
                    if (fr != 0) $fclose(fr);
                    $finish;
                end
                repeat (32) @(posedge clk);
            end
        end
    end
endmodule
