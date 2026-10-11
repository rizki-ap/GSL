// ============================================================
//  tb_xspec_phat.v
//  System testbench: 4-channel windows -> fft2048 (x2) -> xspec_phat
//  -> fft2048 IFFT (x3) -> coarse TDOA peaks, using the shared
//  fft_buf_bank. The testbench plays gcc_loader / gcc_engine.
//
//  Tests
//    1. Noise, integer delays (0, 7, -12, 25)        SW path, all bins
//    2. Plane-wave N-wave, az 37 el 10, noise A/100  SW path, all bins
//    3. Plane-wave Friedlander, az 200 el 5          MB path (20 kHz)
//    4. Dead mic (ch3 = idle +/-1 LSB), band 20..500, floor 27
//
//  RTL checks here: cycle budgets, fft_err, coarse lags vs truth,
//  band mask, dead-channel suppression, no client-port conflicts.
//  Bit-exactness and sub-sample accuracy are checked by
//      python3 xspec_phat_model.py check tb_xp_*.txt
// ============================================================
`timescale 1ns/1ps

module tb_xspec_phat;

    localparam integer N         = 2048;
    localparam integer LOG2N     = 11;
    localparam integer XP_BUDGET = 6200;
    localparam integer FFT_BUDGET= 11400;
    localparam real    PI        = 3.14159265358979323846;
    localparam real    MIC_L     = 0.30;      // tetrahedron edge, m
    localparam real    C_SND     = 343.0;     // m/s
    localparam real    AMP       = 20000.0;   // signal amplitude (s16 range, as fe_cond)

    reg clk = 1'b0;
    always #5 clk = ~clk;
    reg rst_n = 1'b0;

    // ------------------------------------------------------------
    // fft2048
    // ------------------------------------------------------------
    reg        fft_start = 1'b0, fft_inverse = 1'b0;
    reg  [2:0] fft_buf_sel = 3'd0;
    wire       fft_busy, fft_done, fft_err;
    wire [4:0] fft_bfp_exp;
    wire        mem_active;
    wire [2:0]  mem_sel;
    wire [9:0]  mem_raddr0, mem_raddr1, mem_waddr0, mem_waddr1;
    wire [35:0] mem_rdata0, mem_rdata1, mem_wdata0, mem_wdata1;
    wire        mem_we0, mem_we1;

    fft2048 #(.TW_INIT_RE("../fft2048/tw2048_re.hex"),
              .TW_INIT_IM("../fft2048/tw2048_im.hex")) u_fft (
        .clk(clk), .rst_n(rst_n),
        .fft_start(fft_start), .fft_inverse(fft_inverse), .fft_buf_sel(fft_buf_sel),
        .fft_busy(fft_busy), .fft_done(fft_done), .fft_bfp_exp(fft_bfp_exp), .fft_err(fft_err),
        .mem_active(mem_active), .mem_sel(mem_sel),
        .mem_raddr0(mem_raddr0), .mem_raddr1(mem_raddr1),
        .mem_rdata0(mem_rdata0), .mem_rdata1(mem_rdata1),
        .mem_we0(mem_we0), .mem_we1(mem_we1),
        .mem_waddr0(mem_waddr0), .mem_waddr1(mem_waddr1),
        .mem_wdata0(mem_wdata0), .mem_wdata1(mem_wdata1)
    );

    // ------------------------------------------------------------
    // xspec_phat
    // ------------------------------------------------------------
    reg         xp_start = 1'b0;
    reg  [10:0] band_lo  = 11'd0;
    reg  [10:0] band_hi  = 11'd1024;
    reg  [6:0]  floor_r  = 7'd0;
    reg  [4:0]  ez0_r    = 5'd0, ez1_r = 5'd0;
    wire        xp_busy, xp_done, xp_z_re;
    wire [10:0] xp_z_raddr, xp_w_waddr;
    wire [2:0]  xp_w_we;
    wire [35:0] xp_w_wdata;
    wire [179:0] cl_rdata;

    xspec_phat u_xp (
        .clk(clk), .rst_n(rst_n),
        .xp_start(xp_start), .cfg_band_lo(band_lo), .cfg_band_hi(band_hi),
        .cfg_floor(floor_r), .z0_exp(ez0_r), .z1_exp(ez1_r),
        .xp_busy(xp_busy), .xp_done(xp_done),
        .z_re(xp_z_re), .z_raddr(xp_z_raddr),
        .z0_rdata(cl_rdata[0 +: 36]), .z1_rdata(cl_rdata[36 +: 36]),
        .w_we(xp_w_we), .w_waddr(xp_w_waddr), .w_wdata(xp_w_wdata)
    );

    // ------------------------------------------------------------
    // Host (testbench) access + client-port mux (gcc_engine's job)
    //   buffers 0,1 (Z0,Z1) : read by xspec_phat while it runs
    //   buffers 2..4 (W0..W2): written by xspec_phat while it runs
    // ------------------------------------------------------------
    reg  [4:0]   h_re = 5'd0, h_we = 5'd0;
    reg  [10:0]  h_raddr = 11'd0, h_waddr = 11'd0;
    reg  [179:0] h_wdata = 180'd0;

    wire [4:0]   cl_re, cl_we;
    wire [54:0]  cl_raddr, cl_waddr;
    wire [179:0] cl_wdata;
    wire         cl_conflict;

    genvar g;
    generate
        for (g = 0; g < 5; g = g + 1) begin : g_mux
            if (g < 2) begin : g_z
                assign cl_re[g]               = xp_busy ? xp_z_re    : h_re[g];
                assign cl_raddr[g*11 +: 11]   = xp_busy ? xp_z_raddr : h_raddr;
                assign cl_we[g]               = h_we[g];
                assign cl_waddr[g*11 +: 11]   = h_waddr;
                assign cl_wdata[g*36 +: 36]   = h_wdata[g*36 +: 36];
            end else begin : g_w
                assign cl_re[g]               = h_re[g];
                assign cl_raddr[g*11 +: 11]   = h_raddr;
                assign cl_we[g]               = xp_busy ? xp_w_we[g-2] : h_we[g];
                assign cl_waddr[g*11 +: 11]   = xp_busy ? xp_w_waddr   : h_waddr;
                assign cl_wdata[g*36 +: 36]   = xp_busy ? xp_w_wdata   : h_wdata[g*36 +: 36];
            end
        end
    endgenerate

    fft_buf_bank #(.NBUF(5), .DW(36), .AW(11)) u_bank (
        .clk(clk),
        .eng_active(mem_active), .eng_sel(mem_sel),
        .eng_raddr0(mem_raddr0), .eng_raddr1(mem_raddr1),
        .eng_rdata0(mem_rdata0), .eng_rdata1(mem_rdata1),
        .eng_we0(mem_we0), .eng_we1(mem_we1),
        .eng_waddr0(mem_waddr0), .eng_waddr1(mem_waddr1),
        .eng_wdata0(mem_wdata0), .eng_wdata1(mem_wdata1),
        .cl_re(cl_re), .cl_raddr(cl_raddr), .cl_rdata(cl_rdata),
        .cl_we(cl_we), .cl_waddr(cl_waddr), .cl_wdata(cl_wdata),
        .cl_conflict(cl_conflict)
    );

    // ------------------------------------------------------------
    // Storage
    // ------------------------------------------------------------
    integer x   [0:4*N-1];         // 4 channel windows, zero-padded
    integer zre [0:2*N-1], zim [0:2*N-1];
    integer wre [0:3*N-1], wim [0:3*N-1];
    integer rre [0:3*N-1], rim [0:3*N-1];
    integer rd_re [0:N-1], rd_im [0:N-1];
    real    dly [0:3];             // per-channel arrival delay, samples
    integer coarse [0:5];
    integer ez0, ez1, er0, er1, er2;

    integer fails, test_fails, conflicts, seed;
    integer cycles, exp_val, err_val, xp_cycles;
    integer t_len, t_maxlag, t_lo, t_hi, t_floor;
    real    t_fs;

    // pair (a, b) for pair index j
    function integer pa; input integer j;
        case (j) 0,1,2: pa = 0; 3,4: pa = 1; default: pa = 2; endcase
    endfunction
    function integer pb; input integer j;
        case (j) 0: pb = 1; 1,3: pb = 2; default: pb = 3; endcase
    endfunction

    // ------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------
    function integer bitrev11;
        input integer v;
        integer i;
        begin
            bitrev11 = 0;
            for (i = 0; i < LOG2N; i = i + 1)
                if (v & (1 << i)) bitrev11 = bitrev11 | (1 << (LOG2N - 1 - i));
        end
    endfunction

    function integer qround;
        input real v;
        begin
            if (v >= 0.0) qround =  $rtoi( v + 0.5);
            else          qround = -$rtoi(-v + 0.5);
            if (qround >  32767) qround =  32767;
            if (qround < -32768) qround = -32768;
        end
    endfunction

    function [17:0] s18;
        input integer v;
        s18 = v[17:0];
    endfunction

    function integer from_s18;
        input [17:0] v;
        from_s18 = $signed(v);
    endfunction

    function real sigm;
        input real v;
        begin
            if (v < -60.0)     sigm = 0.0;
            else if (v > 60.0) sigm = 1.0;
            else               sigm = 1.0 / (1.0 + $exp(-v));
        end
    endfunction

    // N-wave: duration 0.3 ms, edges ~5 us
    function real nwave;
        input real t;
        begin
            nwave = (1.0 - 2.0 * t / 0.3e-3) * (sigm(t / 5.0e-6) - sigm((t - 0.3e-3) / 5.0e-6));
        end
    endfunction

    // Friedlander: positive phase 2 ms, onset ~50 us
    function real fried;
        input real t;
        begin
            if (t < -3.0e-3) fried = 0.0;
            else             fried = (1.0 - t / 2.0e-3) * $exp(-t / 2.0e-3) * sigm(t / 50.0e-6);
        end
    endfunction

    // Approx. unit-variance noise (sum of 3 uniforms)
    function real gnoise;
        input integer r1, r2, r3;
        gnoise = ((r1 % 10000) + (r2 % 10000) + (r3 % 10000)) / 10000.0;
    endfunction

    // Plane-wave delays (samples) for the tetrahedral array
    task geometry;
        input real az_deg, el_deg, fs;
        real r_, h_, az, el, ux, uy, uz, a;
        integer i;
        begin
            r_ = MIC_L / $sqrt(3.0);
            h_ = MIC_L * $sqrt(2.0 / 3.0);
            az = az_deg * PI / 180.0;
            el = el_deg * PI / 180.0;
            ux = $cos(el) * $cos(az);
            uy = $cos(el) * $sin(az);
            uz = $sin(el);
            for (i = 0; i < 3; i = i + 1) begin
                a = (90.0 + 120.0 * i) * PI / 180.0;
                dly[i] = -(r_ * $cos(a) * ux + r_ * $sin(a) * uy - h_ / 4.0 * uz) / C_SND * fs;
            end
            dly[3] = -(3.0 * h_ / 4.0 * uz) / C_SND * fs;
        end
    endtask

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

    // ------------------------------------------------------------
    // Engine / bank operations
    // ------------------------------------------------------------
    // gcc_loader's job: Z0 = ch0 + j ch1, Z1 = ch2 + j ch3 at bitrev(n)
    task load_z;
        integer n;
        begin
            for (n = 0; n < N; n = n + 1) begin
                @(negedge clk);
                h_we    = 5'b00011;
                h_waddr = bitrev11(n);
                h_wdata[0  +: 36] = {s18(x[0*N+n]), s18(x[1*N+n])};
                h_wdata[36 +: 36] = {s18(x[2*N+n]), s18(x[3*N+n])};
            end
            @(negedge clk);
            h_we = 5'd0;
        end
    endtask

    task read_buf;
        input integer b;
        input         rev;          // 1: read W[k] from address bitrev(k)
        integer k;
        begin
            for (k = 0; k < N; k = k + 1) begin
                @(negedge clk);
                h_re    = 5'd1 << b;
                h_raddr = rev ? bitrev11(k) : k;
                @(negedge clk);
                h_re = 5'd0;
                rd_re[k] = from_s18(cl_rdata[b*36+18 +: 18]);
                rd_im[k] = from_s18(cl_rdata[b*36    +: 18]);
            end
        end
    endtask

    task run_fft;
        input [2:0] b;
        input       inv;
        begin
            @(negedge clk);
            fft_buf_sel = b; fft_inverse = inv; fft_start = 1'b1;
            @(negedge clk);
            fft_start = 1'b0;
            cycles = 1;
            while (!fft_done) begin
                @(negedge clk);
                cycles = cycles + 1;
            end
            exp_val = fft_bfp_exp;
            err_val = fft_err;
            check(err_val == 0,            "fft_err asserted");
            check(cycles <= FFT_BUDGET,    "fft cycle budget");
        end
    endtask

    task run_xspec;
        begin
            @(negedge clk);
            ez0_r = ez0; ez1_r = ez1;
            band_lo = t_lo; band_hi = t_hi; floor_r = t_floor;
            xp_start = 1'b1;
            @(negedge clk);
            xp_start  = 1'b0;
            xp_cycles = 1;
            while (!xp_done) begin
                @(negedge clk);
                xp_cycles = xp_cycles + 1;
                if (xp_cycles > 4 * XP_BUDGET) begin
                    $display("FATAL: xp_done timeout");
                    $finish;
                end
            end
        end
    endtask

    // ------------------------------------------------------------
    // One complete GCC-PHAT job + checks + dump
    // ------------------------------------------------------------
    task run_job;
        input integer id;
        integer k, p, j, n, best, v, lag, fd, maxabs;
        real truth;
        begin
            load_z;
            run_fft(3'd0, 1'b0); ez0 = exp_val;
            run_fft(3'd1, 1'b0); ez1 = exp_val;
            read_buf(0, 0); for (k = 0; k < N; k = k + 1) begin zre[k]   = rd_re[k]; zim[k]   = rd_im[k]; end
            read_buf(1, 0); for (k = 0; k < N; k = k + 1) begin zre[N+k] = rd_re[k]; zim[N+k] = rd_im[k]; end

            run_xspec;
            $display("  xspec_phat  : %0d cycles (budget %0d)", xp_cycles, XP_BUDGET);
            check(xp_cycles <= XP_BUDGET, "xspec_phat cycle budget");
            for (p = 0; p < 3; p = p + 1) begin
                read_buf(2 + p, 1);
                for (k = 0; k < N; k = k + 1) begin wre[p*N+k] = rd_re[k]; wim[p*N+k] = rd_im[k]; end
            end

            run_fft(3'd2, 1'b1); er0 = exp_val;
            run_fft(3'd3, 1'b1); er1 = exp_val;
            run_fft(3'd4, 1'b1); er2 = exp_val;
            $display("  fft2048     : %0d cycles each, exps Z %0d/%0d, R %0d/%0d/%0d",
                     cycles, ez0, ez1, er0, er1, er2);
            for (p = 0; p < 3; p = p + 1) begin
                read_buf(2 + p, 0);
                for (k = 0; k < N; k = k + 1) begin rre[p*N+k] = rd_re[k]; rim[p*N+k] = rd_im[k]; end
            end

            // coarse peaks (peak_search's job): argmax over lags -L..L
            for (j = 0; j < 6; j = j + 1) begin
                best = -2147483647;
                lag  = 0;
                maxabs = 0;
                for (n = -t_maxlag; n <= t_maxlag; n = n + 1) begin
                    v = (j % 2 == 0) ? rre[(j/2)*N + ((n + N) % N)] : rim[(j/2)*N + ((n + N) % N)];
                    if (v > best) begin best = v; lag = n; end
                    if (v > maxabs) maxabs = v;
                    if (-v > maxabs) maxabs = -v;
                end
                coarse[j] = lag;
                truth = dly[pb(j)] - dly[pa(j)];
                if (id == 4 && pb(j) == 3) begin
                    $display("  pair (%0d,%0d) : dead channel, max|r| = %0d", pa(j), pb(j), maxabs);
                    check(maxabs <= 2, "dead-channel correlation not suppressed");
                end else begin
                    $display("  pair (%0d,%0d) : truth %8.3f  coarse %4d", pa(j), pb(j), truth, lag);
                    check((lag - truth) < 1.0 && (truth - lag) < 1.0, "coarse lag not within 1 sample");
                end
            end

            // band mask: every W bin outside [lo, hi] and its mirror must be 0
            if (t_lo > 0 || t_hi < N/2) begin
                v = 0;
                for (p = 0; p < 3; p = p + 1)
                    for (k = 0; k < N; k = k + 1)
                        if (((k < t_lo) || (k > t_hi)) && ((N - k) % N < t_lo || (N - k) % N > t_hi))
                            if (wre[p*N+k] != 0 || wim[p*N+k] != 0) v = v + 1;
                $display("  band mask   : %0d non-zero bins outside [%0d, %0d]", v, t_lo, t_hi);
                check(v == 0, "band mask leak");
            end

            // dump for the Python bit-true / accuracy check
            case (id)
                1: fd = $fopen("tb_xp_1.txt", "w");
                2: fd = $fopen("tb_xp_2.txt", "w");
                3: fd = $fopen("tb_xp_3.txt", "w");
                default: fd = $fopen("tb_xp_4.txt", "w");
            endcase
            $fdisplay(fd, "# test %0d fs %0d len %0d maxlag %0d lo %0d hi %0d floor %0d",
                      id, $rtoi(t_fs), t_len, t_maxlag, t_lo, t_hi, t_floor);
            $fdisplay(fd, "# d0 %f d1 %f d2 %f d3 %f", dly[0], dly[1], dly[2], dly[3]);
            $fdisplay(fd, "# ez0 %0d ez1 %0d er0 %0d er1 %0d er2 %0d", ez0, ez1, er0, er1, er2);
            for (n = 0; n < N; n = n + 1)
                $fdisplay(fd, "%0d %0d %0d %0d  %0d %0d %0d %0d  %0d %0d %0d %0d %0d %0d  %0d %0d %0d %0d %0d %0d",
                    x[n], x[N+n], x[2*N+n], x[3*N+n],
                    zre[n], zim[n], zre[N+n], zim[N+n],
                    wre[n], wim[n], wre[N+n], wim[N+n], wre[2*N+n], wim[2*N+n],
                    rre[n], rim[n], rre[N+n], rim[N+n], rre[2*N+n], rim[2*N+n]);
            $fclose(fd);
        end
    endtask

    task start_test;
        input integer id;
        input [8*64:1] name;
        begin
            test_fails = 0;
            $display("");
            $display("TEST %0d: %0s", id, name);
        end
    endtask

    task end_test;
        begin
            if (test_fails == 0) $display("  -> PASS");
            else                 $display("  -> FAIL (%0d checks)", test_fails);
            fails = fails + test_fails;
        end
    endtask

    always @(posedge clk)
        if (cl_conflict) conflicts = conflicts + 1;

    // ------------------------------------------------------------
    // Main
    // ------------------------------------------------------------
    integer n, i, di;
    real    t, nz;
    integer noise_buf [0:N+63];

    initial begin
        if ($test$plusargs("vcd")) begin
            if ($test$plusargs("fst")) $dumpfile("tb_xspec_phat.fst");   // with vvp -fst
            else                       $dumpfile("tb_xspec_phat.vcd");
            $dumpvars(0, tb_xspec_phat);
        end
        fails = 0; conflicts = 0; seed = 1234;

        repeat (5) @(negedge clk);
        rst_n = 1'b1;
        repeat (5) @(negedge clk);

        $display("============================================================");
        $display(" xspec_phat system testbench (fft2048 x2 -> xspec_phat -> IFFT x3)");
        $display("============================================================");

        // ---------------- TEST 1 ----------------
        start_test(1, "noise, integer delays (0, 7, -12, 25), SW path, all bins");
        t_fs = 100000.0; t_len = 1800; t_maxlag = 100; t_lo = 0; t_hi = 1024; t_floor = 0;
        dly[0] = 0.0; dly[1] = 7.0; dly[2] = -12.0; dly[3] = 25.0;
        for (n = 0; n < N + 64; n = n + 1) noise_buf[n] = qround(AMP * (($random(seed) % 10000) / 10000.0));
        for (i = 0; i < 4; i = i + 1)
            for (n = 0; n < N; n = n + 1)
                x[i*N+n] = (n < t_len) ? noise_buf[n + 32 - $rtoi(dly[i])] : 0;
        run_job(1);
        end_test;

        // ---------------- TEST 2 ----------------
        start_test(2, "plane-wave N-wave, az 37 el 10, noise A/100, SW path");
        t_fs = 100000.0; t_len = 1800; t_maxlag = 100; t_lo = 0; t_hi = 1024; t_floor = 0;
        geometry(37.0, 10.0, t_fs);
        for (i = 0; i < 4; i = i + 1)
            for (n = 0; n < N; n = n + 1) begin
                t  = (n - 300 - dly[i]) / t_fs;
                nz = gnoise($random(seed), $random(seed), $random(seed));
                x[i*N+n] = (n < t_len) ? qround(AMP * nwave(t) + AMP / 100.0 * nz) : 0;
            end
        run_job(2);
        end_test;

        // ---------------- TEST 3 ----------------
        start_test(3, "plane-wave Friedlander, az 200 el 5, MB path (20 kHz)");
        t_fs = 20000.0; t_len = 1720; t_maxlag = 20; t_lo = 0; t_hi = 1024; t_floor = 0;
        geometry(200.0, 5.0, t_fs);
        for (i = 0; i < 4; i = i + 1)
            for (n = 0; n < N; n = n + 1) begin
                t  = (n - 120 - dly[i]) / t_fs;
                nz = gnoise($random(seed), $random(seed), $random(seed));
                x[i*N+n] = (n < t_len) ? qround(AMP * fried(t) + AMP / 100.0 * nz) : 0;
            end
        run_job(3);
        end_test;

        // ---------------- TEST 4 ----------------
        start_test(4, "dead mic (ch3 idle +/-1 LSB), band 20..500, floor 27");
        t_fs = 100000.0; t_len = 1800; t_maxlag = 100; t_lo = 20; t_hi = 500; t_floor = 27;
        dly[0] = 0.0; dly[1] = 7.0; dly[2] = -12.0; dly[3] = 0.0;
        for (n = 0; n < N + 64; n = n + 1) noise_buf[n] = qround(AMP * (($random(seed) % 10000) / 10000.0));
        for (i = 0; i < 3; i = i + 1)
            for (n = 0; n < N; n = n + 1)
                x[i*N+n] = (n < t_len) ? noise_buf[n + 32 - $rtoi(dly[i])] : 0;
        for (n = 0; n < N; n = n + 1)
            x[3*N+n] = (n < t_len) ? ($random(seed) % 2) : 0;
        run_job(4);
        end_test;

        $display("");
        $display("============================================================");
        test_fails = 0;
        check(conflicts == 0, "client port conflict with engine buffer");
        fails = fails + test_fails;
        if (fails == 0) $display(" ALL TESTS PASSED");
        else            $display(" %0d CHECK(S) FAILED", fails);
        $display("============================================================");
        $finish;
    end

endmodule
