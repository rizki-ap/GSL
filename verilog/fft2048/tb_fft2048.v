// ============================================================
//  tb_fft2048.v
//  Self-checking testbench for fft2048 + fft_buf_bank (v2 client ports).
//
//  Tests (each in a different buffer, to exercise the bank mux):
//    1. DC                 re = 0.25 FS, im = 0                  (buffer 0)
//    2. One frequency      re = 0.40 cos(2*pi*37 n/N)            (buffer 1)
//    3. Two frequencies    re = 0.25 cos(2*pi*100 n/N)
//                             + 0.15 sin(2*pi*333 n/N)           (buffer 2)
//    4. Arbitrary complex  0.3 * complex chirp (bin 20 -> 600)
//                          + uniform complex noise (+/-0.1)     (buffer 3)
//    5. Bonus: IFFT round trip of test 4's spectrum              (buffer 4)
//
//  Every test is checked against a double-precision DFT of the
//  SAME quantised input, so the SNR measures only the engine's
//  fixed-point arithmetic. Tests 1-3 also check bin positions,
//  amplitudes and phases against the analytic expectation.
//
//  Each test writes tb_vec_<n>.txt (input, output, exponent) for
//  the bit-true check:  python3 fft2048_model.py check tb_vec_*.txt
//
//  Run:  make           (Icarus Verilog)
//        make wave      (also dumps tb_fft2048.vcd)
// ============================================================
`timescale 1ns/1ps

module tb_fft2048;

    // ------------------------------------------------------------
    // Parameters and pass criteria
    // ------------------------------------------------------------
    localparam integer N            = 2048;
    localparam integer LOG2N        = 11;
    localparam real    FS           = 131072.0;   // Q1.17 full scale
    localparam real    PI           = 3.14159265358979323846;
    localparam integer CYCLE_BUDGET = 11400;      // per transform
    // Engine arithmetic SNR. 75 dB sits >= 15 dB below the best input
    // SNR the array can deliver (16-bit ADC ~90 dB; acoustic SNR is far
    // lower), so FFT rounding never limits the TDOA estimate.
    // Broadband full-scale input measures ~79-82 dB; tones ~100 dB.
    localparam real    SNR_MIN_DB   = 75.0;
    localparam real    RT_SNR_MIN   = 70.0;       // FFT -> IFFT round trip
    localparam real    AMP_TOL      = 1.0e-3;     // relative, analytic bin checks

    // ------------------------------------------------------------
    // Clock / reset
    // ------------------------------------------------------------
    reg clk = 1'b0;
    always #5 clk = ~clk;                          // 100 MHz
    reg rst_n = 1'b0;

    // ------------------------------------------------------------
    // DUT
    // ------------------------------------------------------------
    reg        fft_start   = 1'b0;
    reg        fft_inverse = 1'b0;
    reg  [2:0] fft_buf_sel = 3'd0;
    wire       fft_busy, fft_done, fft_err;
    wire [4:0] fft_bfp_exp;

    wire        mem_active;
    wire [2:0]  mem_sel;
    wire [9:0]  mem_raddr0, mem_raddr1, mem_waddr0, mem_waddr1;
    wire [35:0] mem_rdata0, mem_rdata1, mem_wdata0, mem_wdata1;
    wire        mem_we0, mem_we1;

    // Host access through the bank's client ports (address/data broadcast,
    // per-buffer strobes)
    reg  [4:0]   h_re    = 5'd0;
    reg  [4:0]   h_we    = 5'd0;
    reg  [10:0]  h_raddr = 11'd0;
    reg  [10:0]  h_waddr = 11'd0;
    reg  [35:0]  h_wdata = 36'd0;
    wire [179:0] cl_rdata;
    wire         cl_conflict;

    fft2048 dut (
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

    fft_buf_bank #(.NBUF(5), .DW(36), .AW(11)) bank (
        .clk(clk),
        .eng_active(mem_active), .eng_sel(mem_sel),
        .eng_raddr0(mem_raddr0), .eng_raddr1(mem_raddr1),
        .eng_rdata0(mem_rdata0), .eng_rdata1(mem_rdata1),
        .eng_we0(mem_we0), .eng_we1(mem_we1),
        .eng_waddr0(mem_waddr0), .eng_waddr1(mem_waddr1),
        .eng_wdata0(mem_wdata0), .eng_wdata1(mem_wdata1),
        .cl_re(h_re), .cl_raddr({5{h_raddr}}), .cl_rdata(cl_rdata),
        .cl_we(h_we), .cl_waddr({5{h_waddr}}), .cl_wdata({5{h_wdata}}),
        .cl_conflict(cl_conflict)
    );

    // ------------------------------------------------------------
    // Storage
    // ------------------------------------------------------------
    integer xin_re  [0:N-1];
    integer xin_im  [0:N-1];
    integer xout_re [0:N-1];
    integer xout_im [0:N-1];
    integer x4_re   [0:N-1];       // copy of test 4 input for the round trip
    integer x4_im   [0:N-1];
    real    ref_re  [0:N-1];
    real    ref_im  [0:N-1];
    real    ctab    [0:N-1];
    real    stab    [0:N-1];
    reg     used    [0:N-1];
    integer top     [0:7];

    integer cycles, exp_val, err_val, exp4;
    integer fails, test_fails, conflicts;
    real    snr_db, max_err_lsb;
    integer seed;

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

    function integer qround;            // round half away from zero
        input real x;
        begin
            if (x >= 0.0) qround =  $rtoi( x + 0.5);
            else          qround = -$rtoi(-x + 0.5);
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

    function real mag;
        input real re, im;
        mag = $sqrt(re * re + im * im);
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

    // ------------------------------------------------------------
    // Buffer access through the bank's external port
    // ------------------------------------------------------------
    // Producer contract: x[n] goes to address bitrev11(n)
    task load_input;
        input [2:0] b;
        integer n;
        begin
            for (n = 0; n < N; n = n + 1) begin
                @(negedge clk);
                h_we    = 5'd1 << b;
                h_waddr = bitrev11(n);
                h_wdata = {s18(xin_re[n]), s18(xin_im[n])};
            end
            @(negedge clk);
            h_we = 5'd0;
        end
    endtask

    // Output is in natural order
    task read_output;
        input [2:0] b;
        integer k;
        begin
            for (k = 0; k < N; k = k + 1) begin
                @(negedge clk);
                h_re    = 5'd1 << b;
                h_raddr = k;
                @(negedge clk);
                h_re = 5'd0;
                xout_re[k] = from_s18(cl_rdata[b*36+18 +: 18]);
                xout_im[k] = from_s18(cl_rdata[b*36    +: 18]);
            end
        end
    endtask

    task run_fft;
        input [2:0] b;
        input       inv;
        begin
            @(negedge clk);
            fft_buf_sel = b;
            fft_inverse = inv;
            fft_start   = 1'b1;
            @(negedge clk);
            fft_start = 1'b0;
            cycles    = 1;
            while (!fft_done) begin
                @(negedge clk);
                cycles = cycles + 1;
                if (cycles > 5 * CYCLE_BUDGET) begin
                    $display("FATAL: fft_done timeout");
                    $finish;
                end
            end
            exp_val = fft_bfp_exp;
            err_val = fft_err;
        end
    endtask

    // ------------------------------------------------------------
    // Reference DFT (double precision) of the quantised input
    //   forward: X[k] = sum x[n] (cos - j sin)(2 pi n k / N)
    //   inverse: X[k] = sum x[n] (cos + j sin)(2 pi n k / N)   (no 1/N)
    // ------------------------------------------------------------
    task compute_ref;
        input inv;
        integer k, n, m;
        real sr, si, c, s;
        begin
            for (k = 0; k < N; k = k + 1) begin
                sr = 0.0;
                si = 0.0;
                for (n = 0; n < N; n = n + 1) begin
                    m  = (n * k) & (N - 1);
                    c  = ctab[m];
                    s  = inv ? -stab[m] : stab[m];
                    sr = sr + xin_re[n] * c + xin_im[n] * s;
                    si = si + xin_im[n] * c - xin_re[n] * s;
                end
                ref_re[k] = sr;
                ref_im[k] = si;
            end
        end
    endtask

    // SNR of (output * 2^exp) against the reference, and the worst
    // single-bin error expressed in output LSBs.
    task evaluate;
        integer k;
        real scale, orr, oi, er, ei, sig, noise, e;
        begin
            scale       = 2.0 ** exp_val;
            sig         = 0.0;
            noise       = 0.0;
            max_err_lsb = 0.0;
            for (k = 0; k < N; k = k + 1) begin
                orr   = xout_re[k] * scale;
                oi    = xout_im[k] * scale;
                er    = orr - ref_re[k];
                ei    = oi  - ref_im[k];
                sig   = sig   + ref_re[k] * ref_re[k] + ref_im[k] * ref_im[k];
                noise = noise + er * er + ei * ei;
                e     = mag(er, ei) / scale;
                if (e > max_err_lsb) max_err_lsb = e;
            end
            snr_db = (noise > 0.0) ? 10.0 * $log10(sig / noise) : 999.0;
        end
    endtask

    // Indices of the 'cnt' largest-magnitude output bins
    task top_bins;
        input integer cnt;
        integer i, k, best;
        real bm, m2;
        begin
            for (k = 0; k < N; k = k + 1) used[k] = 1'b0;
            for (i = 0; i < cnt; i = i + 1) begin
                best = 0;
                bm   = -1.0;
                for (k = 0; k < N; k = k + 1) begin
                    m2 = 1.0 * xout_re[k] * xout_re[k] + 1.0 * xout_im[k] * xout_im[k];
                    if (!used[k] && m2 > bm) begin
                        bm   = m2;
                        best = k;
                    end
                end
                used[best] = 1'b1;
                top[i]     = best;
            end
        end
    endtask

    function in_top;
        input integer bin, cnt;
        integer i;
        begin
            in_top = 1'b0;
            for (i = 0; i < cnt; i = i + 1)
                if (top[i] == bin) in_top = 1'b1;
        end
    endfunction

    // Output bin value in true (unscaled) Q1.17 units
    function real out_re_true;
        input integer k;
        out_re_true = xout_re[k] * (2.0 ** exp_val);
    endfunction

    function real out_im_true;
        input integer k;
        out_im_true = xout_im[k] * (2.0 ** exp_val);
    endfunction

    function real rel_err;
        input real got, want;
        rel_err = (got - want) / want;
    endfunction

    // Dump input/output for the Python bit-true model
    task dump_vectors;
        input integer id;
        input inv;
        integer fd, n;
        begin
            case (id)
                1: fd = $fopen("tb_vec_1.txt", "w");
                2: fd = $fopen("tb_vec_2.txt", "w");
                3: fd = $fopen("tb_vec_3.txt", "w");
                4: fd = $fopen("tb_vec_4.txt", "w");
                default: fd = $fopen("tb_vec_5.txt", "w");
            endcase
            $fdisplay(fd, "# test %0d inverse %0d exp %0d err %0d", id, inv, exp_val, err_val);
            for (n = 0; n < N; n = n + 1)
                $fdisplay(fd, "%0d %0d %0d %0d", xin_re[n], xin_im[n], xout_re[n], xout_im[n]);
            $fclose(fd);
        end
    endtask

    // Common steps after a transform: read back, reference, report
    task finish_test;
        input integer id;
        input [2:0]   b;
        input         inv;
        begin
            read_output(b);
            compute_ref(inv);
            evaluate;
            dump_vectors(id, inv);
            $display("  cycles      : %0d (budget %0d)", cycles, CYCLE_BUDGET);
            $display("  bfp_exp     : %0d", exp_val);
            $display("  SNR         : %0.1f dB (min %0.1f)", snr_db, SNR_MIN_DB);
            $display("  max bin err : %0.2f LSB", max_err_lsb);
            check(cycles <= CYCLE_BUDGET, "cycle budget exceeded");
            check(err_val == 0,           "fft_err asserted");
            check(snr_db >= SNR_MIN_DB,   "SNR below minimum");
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
        input integer id;
        begin
            if (test_fails == 0) $display("  -> PASS");
            else                 $display("  -> FAIL (%0d checks)", test_fails);
            fails = fails + test_fails;
        end
    endtask

    // ------------------------------------------------------------
    // Monitor: ext port must never touch the engine's buffer
    // ------------------------------------------------------------
    always @(posedge clk)
        if (cl_conflict) conflicts = conflicts + 1;

    // ------------------------------------------------------------
    // Main
    // ------------------------------------------------------------
    integer n, k;
    integer k1, k2a, k2b;
    real    amp, ampb, ph, f0, f1, want, got;

    initial begin
        if ($test$plusargs("vcd")) begin
            if ($test$plusargs("fst")) $dumpfile("tb_fft2048.fst");   // with vvp -fst
            else                       $dumpfile("tb_fft2048.vcd");
            $dumpvars(0, tb_fft2048);
        end

        fails     = 0;
        conflicts = 0;
        seed      = 20261008;

        for (n = 0; n < N; n = n + 1) begin
            ctab[n] = $cos(2.0 * PI * n / N);
            stab[n] = $sin(2.0 * PI * n / N);
        end

        repeat (5) @(negedge clk);
        rst_n = 1'b1;
        repeat (5) @(negedge clk);

        $display("============================================================");
        $display(" fft2048 testbench   N=%0d  clk=100 MHz", N);
        $display("============================================================");

        // --------------------------------------------------------
        // TEST 1: DC
        // --------------------------------------------------------
        start_test(1, "DC (re = 0.25 FS)");
        amp = 0.25;
        for (n = 0; n < N; n = n + 1) begin
            xin_re[n] = qround(amp * FS);
            xin_im[n] = 0;
        end
        load_input(3'd0);
        run_fft(3'd0, 1'b0);
        finish_test(1, 3'd0, 1'b0);
        top_bins(1);
        want = N * amp * FS;
        $display("  peak bin    : %0d", top[0]);
        $display("  X[0]        = %0.1f %+0.1fj (expect %0.1f)", out_re_true(0), out_im_true(0), want);
        check(top[0] == 0,                                   "peak not at bin 0");
        check(within_tol(rel_err(out_re_true(0), want)),        "X[0] amplitude");
        end_test(1);

        // --------------------------------------------------------
        // TEST 2: one frequency
        // --------------------------------------------------------
        k1  = 37;
        amp = 0.40;
        start_test(2, "one frequency (0.40 cos, bin 37)");
        for (n = 0; n < N; n = n + 1) begin
            xin_re[n] = qround(amp * FS * $cos(2.0 * PI * k1 * n / N));
            xin_im[n] = 0;
        end
        load_input(3'd1);
        run_fft(3'd1, 1'b0);
        finish_test(2, 3'd1, 1'b0);
        top_bins(2);
        want = N * amp * FS / 2.0;
        $display("  peak bins   : %0d, %0d", top[0], top[1]);
        $display("  X[%0d]       = %0.1f %+0.1fj (expect %0.1f)", k1,
                 out_re_true(k1), out_im_true(k1), want);
        check(in_top(k1, 2) && in_top(N - k1, 2),             "peaks not at +/- bin 37");
        check(within_tol(rel_err(out_re_true(k1),     want)),    "X[37] amplitude");
        check(within_tol(rel_err(out_re_true(N - k1), want)),    "X[N-37] amplitude");
        end_test(2);

        // --------------------------------------------------------
        // TEST 3: two frequencies (cos + sin -> checks phase too)
        // --------------------------------------------------------
        k2a  = 100;
        k2b  = 333;
        amp  = 0.25;
        ampb = 0.15;
        start_test(3, "two frequencies (0.25 cos bin 100 + 0.15 sin bin 333)");
        for (n = 0; n < N; n = n + 1) begin
            xin_re[n] = qround(FS * (amp  * $cos(2.0 * PI * k2a * n / N) +
                                     ampb * $sin(2.0 * PI * k2b * n / N)));
            xin_im[n] = 0;
        end
        load_input(3'd2);
        run_fft(3'd2, 1'b0);
        finish_test(3, 3'd2, 1'b0);
        top_bins(4);
        $display("  peak bins   : %0d, %0d, %0d, %0d", top[0], top[1], top[2], top[3]);
        $display("  X[%0d]      = %0.1f %+0.1fj (expect %0.1f real)", k2a,
                 out_re_true(k2a), out_im_true(k2a), N * amp * FS / 2.0);
        $display("  X[%0d]      = %0.1f %+0.1fj (expect %0.1fj)", k2b,
                 out_re_true(k2b), out_im_true(k2b), -N * ampb * FS / 2.0);
        check(in_top(k2a, 4) && in_top(N - k2a, 4) &&
              in_top(k2b, 4) && in_top(N - k2b, 4),             "peaks not at bins 100/333");
        // cos -> +N*A/2 real at +k and -k
        check(within_tol(rel_err(out_re_true(k2a),      N * amp * FS / 2.0)), "X[100] amplitude");
        check(within_tol(rel_err(out_re_true(N - k2a),  N * amp * FS / 2.0)), "X[N-100] amplitude");
        // sin -> -j*N*B/2 at +k, +j*N*B/2 at -k
        check(within_tol(rel_err(out_im_true(k2b),     -N * ampb * FS / 2.0)), "X[333] phase/amplitude");
        check(within_tol(rel_err(out_im_true(N - k2b),  N * ampb * FS / 2.0)), "X[N-333] phase/amplitude");
        end_test(3);

        // --------------------------------------------------------
        // TEST 4: arbitrary complex signal
        //   0.3 * exp(j*phi(n)), linear chirp bin 20 -> 600
        //   + uniform complex noise in +/- 0.1
        // --------------------------------------------------------
        amp = 0.30;
        f0  = 20.0;
        f1  = 600.0;
        start_test(4, "arbitrary complex (chirp 20->600 + noise)");
        for (n = 0; n < N; n = n + 1) begin
            ph        = 2.0 * PI * (f0 * n + 0.5 * (f1 - f0) * n * n / N) / N;
            xin_re[n] = qround(FS * (amp * $cos(ph) + 0.1 * ($random(seed) % 10000) / 10000.0));
            xin_im[n] = qround(FS * (amp * $sin(ph) + 0.1 * ($random(seed) % 10000) / 10000.0));
            x4_re[n]  = xin_re[n];
            x4_im[n]  = xin_im[n];
        end
        load_input(3'd3);
        run_fft(3'd3, 1'b0);
        finish_test(4, 3'd3, 1'b0);
        exp4 = exp_val;
        end_test(4);

        // --------------------------------------------------------
        // TEST 5 (bonus): IFFT round trip of test 4's spectrum.
        //   The spectrum is halved first so it meets the input
        //   contract (|x| < 0.5 FS), as xspec_phat would.
        //   Expected: IFFT(FFT(x)) = N * x
        // --------------------------------------------------------
        start_test(5, "IFFT round trip of test 4");
        for (n = 0; n < N; n = n + 1) begin
            xin_re[n] = qround(xout_re[n] / 2.0);
            xin_im[n] = qround(xout_im[n] / 2.0);
        end
        load_input(3'd4);
        run_fft(3'd4, 1'b1);
        finish_test(5, 3'd4, 1'b1);
        begin : roundtrip
            real s_sig, s_err, xr, xi, sc;
            s_sig = 0.0;
            s_err = 0.0;
            sc    = (2.0 ** (exp_val + exp4 + 1)) / N;
            for (n = 0; n < N; n = n + 1) begin
                xr    = xout_re[n] * sc - x4_re[n];
                xi    = xout_im[n] * sc - x4_im[n];
                s_sig = s_sig + 1.0 * x4_re[n] * x4_re[n] + 1.0 * x4_im[n] * x4_im[n];
                s_err = s_err + xr * xr + xi * xi;
            end
            $display("  round trip  : %0.1f dB vs original x (min %0.1f)",
                     10.0 * $log10(s_sig / s_err), RT_SNR_MIN);
            check(10.0 * $log10(s_sig / s_err) >= RT_SNR_MIN, "round-trip SNR");
        end
        end_test(5);

        // --------------------------------------------------------
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

    // |x| <= AMP_TOL helper (named to read like a check)
    function within_tol;
        input real x;
        within_tol = (x < 0.0 ? -x : x) <= AMP_TOL;
    endfunction

endmodule
