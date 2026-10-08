// ============================================================
//  tb_gcc_engine.v   -- Milestone A testbench
//
//  Data-driven: milestone_a.py prepare writes
//    ma_jobs.txt        line 1: number of jobs
//                       then:   id type start len maxlag lo hi floor
//    ma_ring_<id>.hex   4096 x 72-bit ring image for that job
//  For each job the testbench
//    1. writes the ring image into sw_ring or mb_ring (as fe_cond /
//       mb_decim will), through the ring's write port
//    2. issues the job to gcc_engine and waits for res_valid
//    3. checks: cycle budget (2 ms), res_err, lag-window stream ==
//       correlation buffers, no port conflicts
//    4. dumps results + Z0/Z1 spectra + correlations (debug port) to
//       ma_out_<id>.txt for  milestone_a.py check
//  Finally it runs two jobs back to back to measure throughput.
// ============================================================
`timescale 1ns/1ps

module tb_gcc_engine;

    localparam integer N          = 2048;
    localparam integer JOB_BUDGET = 200000;      // 2 ms @ 100 MHz

    reg clk = 1'b0;
    always #5 clk = ~clk;
    reg rst_n = 1'b0;

    // ------------------------------------------------------------
    // Rings
    // ------------------------------------------------------------
    reg         sw_we = 1'b0, mb_we = 1'b0;
    reg  [11:0] ring_waddr = 12'd0;
    reg  [71:0] ring_wdata = 72'd0;
    wire [11:0] sw_raddr, mb_raddr;
    wire [71:0] sw_rdata, mb_rdata;

    cap_ring u_sw_ring (.clk(clk), .we(sw_we), .waddr(ring_waddr), .wdata(ring_wdata),
                        .raddr(sw_raddr), .rdata(sw_rdata));
    cap_ring u_mb_ring (.clk(clk), .we(mb_we), .waddr(ring_waddr), .wdata(ring_wdata),
                        .raddr(mb_raddr), .rdata(mb_rdata));

    // ------------------------------------------------------------
    // DUT
    // ------------------------------------------------------------
    reg         job_valid = 1'b0, job_type = 1'b0;
    reg  [15:0] job_id = 16'd0;
    reg  [11:0] job_start = 12'd0, job_len = 12'd0;
    reg  [6:0]  job_maxlag = 7'd0;
    reg  [10:0] sw_lo = 11'd0, sw_hi = 11'd1024, mb_lo = 11'd0, mb_hi = 11'd1024;
    reg  [6:0]  sw_floor = 7'd0, mb_floor = 7'd0;
    wire        job_ready, busy;

    wire         res_valid, res_type, res_err;
    wire [15:0]  res_id;
    wire [47:0]  res_lag;
    wire [107:0] res_pk, res_side, lw_data;
    wire [9:0]   res_exp_z;
    wire [14:0]  res_exp_r;
    wire         lw_we;
    wire [7:0]   lw_idx;

    reg          dbg_re = 1'b0;
    reg  [2:0]   dbg_sel = 3'd0;
    reg  [10:0]  dbg_raddr = 11'd0;
    wire [35:0]  dbg_rdata;

    gcc_engine #(
        .TW_INIT_RE("../fft2048/tw2048_re.hex"),
        .TW_INIT_IM("../fft2048/tw2048_im.hex"),
        .RSQRT_INIT("../xspec_phat/rsqrt1024.hex")
    ) dut (
        .clk(clk), .rst_n(rst_n),
        .job_valid(job_valid), .job_ready(job_ready), .job_type(job_type), .job_id(job_id),
        .job_start(job_start), .job_len(job_len), .job_maxlag(job_maxlag),
        .cfg_sw_band_lo(sw_lo), .cfg_sw_band_hi(sw_hi), .cfg_sw_floor(sw_floor),
        .cfg_mb_band_lo(mb_lo), .cfg_mb_band_hi(mb_hi), .cfg_mb_floor(mb_floor),
        .sw_ring_raddr(sw_raddr), .sw_ring_rdata(sw_rdata),
        .mb_ring_raddr(mb_raddr), .mb_ring_rdata(mb_rdata),
        .res_valid(res_valid), .res_type(res_type), .res_id(res_id),
        .res_lag(res_lag), .res_pk(res_pk), .res_side(res_side),
        .res_exp_z(res_exp_z), .res_exp_r(res_exp_r), .res_err(res_err),
        .lw_we(lw_we), .lw_idx(lw_idx), .lw_data(lw_data),
        .dbg_re(dbg_re), .dbg_sel(dbg_sel), .dbg_raddr(dbg_raddr), .dbg_rdata(dbg_rdata),
        .busy(busy)
    );

    // ------------------------------------------------------------
    // Capture the lag-window stream
    // ------------------------------------------------------------
    reg [107:0] lw_mem [0:255];
    integer     lw_count;
    always @(posedge clk)
        if (lw_we) begin
            lw_mem[lw_idx] <= lw_data;
            lw_count = lw_count + 1;
        end

    integer conflicts;
    always @(posedge clk)
        if (dut.cl_conflict) conflicts = conflicts + 1;

    // ------------------------------------------------------------
    // Storage and helpers
    // ------------------------------------------------------------
    reg  [71:0] ring_img [0:4095];
    integer zre [0:2*N-1], zim [0:2*N-1];
    integer rre [0:3*N-1], rim [0:3*N-1];
    integer fails, test_fails, cycles;
    integer j_id, j_type, j_start, j_len, j_maxlag, j_lo, j_hi, j_floor;
    reg [8*32:1] fname;

    function integer s18;
        input [17:0] v;
        s18 = $signed(v);
    endfunction

    function integer s8;
        input [7:0] v;
        s8 = $signed(v);
    endfunction

    task check;
        input cond;
        input [8*64:1] msg;
        begin
            if (!cond) begin
                $display("    FAIL: %0s", msg);
                test_fails = test_fails + 1;
            end
        end
    endtask

    task write_ring;
        input integer typ;
        integer a;
        begin
            for (a = 0; a < 4096; a = a + 1) begin
                @(negedge clk);
                ring_waddr = a;
                ring_wdata = ring_img[a];
                sw_we = (typ == 0);
                mb_we = (typ == 1);
            end
            @(negedge clk);
            sw_we = 1'b0;
            mb_we = 1'b0;
        end
    endtask

    task issue_job;
        begin
            @(negedge clk);
            while (!job_ready) @(negedge clk);
            job_valid  = 1'b1;
            job_id     = j_id;
            job_type   = j_type;
            job_start  = j_start;
            job_len    = j_len;
            job_maxlag = j_maxlag;
            if (j_type == 0) begin sw_lo = j_lo; sw_hi = j_hi; sw_floor = j_floor; end
            else             begin mb_lo = j_lo; mb_hi = j_hi; mb_floor = j_floor; end
            @(negedge clk);
            job_valid = 1'b0;
        end
    endtask

    task wait_result;
        begin
            cycles = 1;
            while (!res_valid) begin
                @(negedge clk);
                cycles = cycles + 1;
                if (cycles > 5 * JOB_BUDGET) begin
                    $display("FATAL: res_valid timeout");
                    $finish;
                end
            end
        end
    endtask

    task read_dbg;
        input integer b;
        input integer to_r;            // 0: into zre/zim slot b, 1: into rre/rim slot b-2
        integer k;
        begin
            for (k = 0; k < N; k = k + 1) begin
                @(negedge clk);
                dbg_re = 1'b1; dbg_sel = b; dbg_raddr = k;
                @(negedge clk);
                dbg_re = 1'b0;
                if (to_r) begin rre[(b-2)*N+k] = s18(dbg_rdata[35:18]); rim[(b-2)*N+k] = s18(dbg_rdata[17:0]); end
                else      begin zre[b*N+k]     = s18(dbg_rdata[35:18]); zim[b*N+k]     = s18(dbg_rdata[17:0]); end
            end
        end
    endtask

    task run_one;
        integer fd, n, j, m, mism, v;
        begin
            test_fails = 0;
            $sformat(fname, "ma_ring_%0d.hex", j_id);
            $readmemh(fname, ring_img);
            write_ring(j_type);
            lw_count = 0;
            issue_job;
            wait_result;

            // read back spectra and correlations through the debug port
            read_dbg(0, 0); read_dbg(1, 0);
            read_dbg(2, 1); read_dbg(3, 1); read_dbg(4, 1);

            // lag-window stream must equal the correlation buffers
            mism = 0;
            for (m = -j_maxlag; m <= j_maxlag; m = m + 1)
                for (j = 0; j < 6; j = j + 1) begin
                    v = (j % 2 == 0) ? rre[(j/2)*N + ((m + N) % N)] : rim[(j/2)*N + ((m + N) % N)];
                    if (s18(lw_mem[m + j_maxlag][j*18 +: 18]) != v) mism = mism + 1;
                end

            $display("  job %0d (%0s): %0d cycles = %0.3f ms, exps Z %0d/%0d R %0d/%0d/%0d, lags %0d %0d %0d %0d %0d %0d",
                     j_id, j_type ? "MB" : "SW", cycles, cycles / 100000.0,
                     res_exp_z[4:0], res_exp_z[9:5], res_exp_r[4:0], res_exp_r[9:5], res_exp_r[14:10],
                     s8(res_lag[7:0]), s8(res_lag[15:8]), s8(res_lag[23:16]),
                     s8(res_lag[31:24]), s8(res_lag[39:32]), s8(res_lag[47:40]));
            check(cycles <= JOB_BUDGET,               "job exceeds 2 ms budget");
            check(res_err == 0,                       "res_err (fft_err or port conflict)");
            check(res_id == j_id,                     "res_id mismatch");
            check(lw_count == 2 * j_maxlag + 1,       "lag-window stream length");
            check(mism == 0,                          "lag-window stream != correlation buffers");

            $sformat(fname, "ma_out_%0d.txt", j_id);
            fd = $fopen(fname, "w");
            $fdisplay(fd, "# id %0d type %0d cycles %0d err %0d", j_id, j_type, cycles, res_err);
            $fdisplay(fd, "# ez0 %0d ez1 %0d er0 %0d er1 %0d er2 %0d",
                      res_exp_z[4:0], res_exp_z[9:5], res_exp_r[4:0], res_exp_r[9:5], res_exp_r[14:10]);
            for (j = 0; j < 6; j = j + 1)
                $fdisplay(fd, "# lag%0d %0d pk%0d %0d side%0d %0d", j, s8(res_lag[j*8 +: 8]),
                          j, s18(res_pk[j*18 +: 18]), j, $unsigned(res_side[j*18 +: 18]));
            for (n = 0; n < N; n = n + 1)
                $fdisplay(fd, "%0d %0d %0d %0d  %0d %0d %0d %0d %0d %0d",
                          zre[n], zim[n], zre[N+n], zim[N+n],
                          rre[n], rim[n], rre[N+n], rim[N+n], rre[2*N+n], rim[2*N+n]);
            $fclose(fd);

            if (test_fails == 0) $display("    -> PASS");
            else                 $display("    -> FAIL (%0d checks)", test_fails);
            fails = fails + test_fails;
        end
    endtask

    // ------------------------------------------------------------
    // Main
    // ------------------------------------------------------------
    integer jf, njobs, i, r, t_first, t_both;

    initial begin
        if ($test$plusargs("vcd")) begin
            if ($test$plusargs("fst")) $dumpfile("tb_gcc_engine.fst");
            else                       $dumpfile("tb_gcc_engine.vcd");
            $dumpvars(0, tb_gcc_engine);
        end
        fails = 0; conflicts = 0;

        jf = $fopen("ma_jobs.txt", "r");
        if (jf == 0) begin
            $display("FATAL: ma_jobs.txt not found -- run: python3 milestone_a.py prepare");
            $finish;
        end
        r = $fscanf(jf, "%d\n", njobs);

        repeat (5) @(negedge clk);
        rst_n = 1'b1;
        repeat (5) @(negedge clk);

        $display("============================================================");
        $display(" Milestone A: gcc_engine  (%0d jobs)", njobs);
        $display("============================================================");

        for (i = 0; i < njobs; i = i + 1) begin
            r = $fscanf(jf, "%d %d %d %d %d %d %d %d\n",
                        j_id, j_type, j_start, j_len, j_maxlag, j_lo, j_hi, j_floor);
            run_one;
        end
        $fclose(jf);

        // Throughput: last job again, twice back to back (ring already loaded)
        $display("");
        $display("  back-to-back: 2 jobs issued together");
        test_fails = 0;
        t_both = 0;
        issue_job;
        wait_result;
        t_first = cycles;
        issue_job;                      // waits for job_ready, then issues
        wait_result;
        t_both = t_first + cycles;
        $display("    total %0d cycles = %0.3f ms for 2 jobs", t_both, t_both / 100000.0);
        check(t_both <= 2 * JOB_BUDGET, "back-to-back exceeds budget");
        fails = fails + test_fails;

        $display("");
        $display("============================================================");
        test_fails = 0;
        check(conflicts == 0, "buffer port conflict");
        fails = fails + test_fails;
        if (fails == 0) $display(" ALL RTL CHECKS PASSED  (now run: python3 milestone_a.py check)");
        else            $display(" %0d CHECK(S) FAILED", fails);
        $display("============================================================");
        $finish;
    end
endmodule
