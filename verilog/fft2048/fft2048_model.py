#!/usr/bin/env python3
"""
fft2048_model.py
Bit-true Python model of verilog fft2048 (radix-2 DIT, s18 Q1.17,
block-floating-point, convergent rounding, IFFT via re/im swap).

It reproduces the RTL exactly, so it can stand in for the hardware
when you evaluate fixed-point TDOA error in the Python pipeline
(gs_det_tdoa.py), per timing_budget.md section 4.3.

Usage
    python3 fft2048_model.py check tb_vec_*.txt   # RTL vs model, must be bit-exact
    python3 fft2048_model.py selftest             # model vs numpy, random input

Library
    from fft2048_model import fft2048
    out_re, out_im, exp, err = fft2048(x_re, x_im, inverse=False)
    # true result = (out_re + 1j*out_im) * 2**exp   (Q1.17 integer units)
"""
import sys

from gen_twiddle2048 import twiddles

N = 2048
LOG2N = 11
FRAC = 17
DW = 18
MAXV = (1 << (DW - 1)) - 1
MINV = -(1 << (DW - 1))

_TW = twiddles(N)


def bitrev(v, bits=LOG2N):
    r = 0
    for i in range(bits):
        if v & (1 << i):
            r |= 1 << (bits - 1 - i)
    return r


def rnd_sat(x, s):
    """Convergent rounding of integer x by s bits, then saturate to s18.
    Mirrors fft2048_bfly.v rnd_sat()."""
    t = x + ((1 << (s - 1)) - 1) + ((x >> s) & 1)
    y = t >> s
    if y > MAXV:
        return MAXV, True
    if y < MINV:
        return MINV, True
    return y, False


def mag1(v):
    """One's-complement magnitude, 17 bits (mirrors fft2048.v mag1())."""
    return (v if v >= 0 else ~v) & ((1 << (DW - 1)) - 1)


def fft2048(x_re, x_im, inverse=False):
    """Bit-true transform of natural-order integer input (s18).
    Returns (out_re, out_im, bfp_exp, err) with natural-order output."""
    assert len(x_re) == N and len(x_im) == N
    # Producer contract: x[n] stored at address bitrev(n)
    mem_re = [0] * N
    mem_im = [0] * N
    for n in range(N):
        mem_re[bitrev(n)] = int(x_re[n])
        mem_im[bitrev(n)] = int(x_im[n])

    err = False
    # Input contract |v| < 2^16
    for v in mem_re + mem_im:
        if not (-(1 << 16) <= v < (1 << 16)):
            err = True
            break

    if inverse:                      # stage-0 read swap
        mem_re, mem_im = mem_im, mem_re

    sh = 0
    exp = 0
    for s in range(LOG2N):
        half = 1 << s
        max_or = 0
        S = FRAC + sh
        for c in range(N // 2):
            ia = ((c >> s) << (s + 1)) | (c & (half - 1))
            ib = ia | half
            k = (c & (half - 1)) << (LOG2N - 1 - s)
            ar, ai = mem_re[ia], mem_im[ia]
            br, bi = mem_re[ib], mem_im[ib]
            if k == 0:
                mr, mi = br << FRAC, bi << FRAC
            else:
                wr, wi = _TW[k]
                mr = br * wr - bi * wi
                mi = br * wi + bi * wr
            ax, ay = ar << FRAC, ai << FRAC
            pr, s1 = rnd_sat(ax + mr, S)
            pi, s2 = rnd_sat(ay + mi, S)
            qr, s3 = rnd_sat(ax - mr, S)
            qi, s4 = rnd_sat(ay - mi, S)
            err |= s1 or s2 or s3 or s4
            mem_re[ia], mem_im[ia] = pr, pi
            mem_re[ib], mem_im[ib] = qr, qi
            max_or |= mag1(pr) | mag1(pi) | mag1(qr) | mag1(qi)
        if s < LOG2N - 1:
            sh = 2 if max_or & (1 << 16) else (1 if max_or & (1 << 15) else 0)
            exp += sh

    if inverse:                      # last-stage write swap
        mem_re, mem_im = mem_im, mem_re
    return mem_re, mem_im, exp, err


def check(files):
    ok = True
    for fn in files:
        with open(fn) as f:
            hdr = f.readline().split()
            inv = int(hdr[hdr.index('inverse') + 1])
            exp_rtl = int(hdr[hdr.index('exp') + 1])
            rows = [list(map(int, line.split())) for line in f if line.strip()]
        xr = [r[0] for r in rows]
        xi = [r[1] for r in rows]
        orr = [r[2] for r in rows]
        oi = [r[3] for r in rows]
        mr, mi, exp, _ = fft2048(xr, xi, inverse=bool(inv))
        mism = sum(1 for k in range(N) if mr[k] != orr[k] or mi[k] != oi[k])
        good = (mism == 0 and exp == exp_rtl)
        ok &= good
        print(f"{fn}: inverse={inv} exp rtl={exp_rtl} model={exp}  "
              f"mismatching bins={mism}  -> {'BIT-EXACT' if good else 'MISMATCH'}")
    print("ALL BIT-EXACT" if ok else "MISMATCHES FOUND")
    return ok


def selftest():
    import numpy as np
    rng = np.random.default_rng(1)
    for inv in (False, True):
        x = rng.uniform(-0.45, 0.45, N) + 1j * rng.uniform(-0.45, 0.45, N)
        xr = np.round(x.real * (1 << FRAC)).astype(int)
        xi = np.round(x.imag * (1 << FRAC)).astype(int)
        orr, oi, exp, err = fft2048(list(xr), list(xi), inv)
        out = (np.array(orr) + 1j * np.array(oi)) * 2.0 ** exp
        xin = xr + 1j * xi
        ref = np.fft.ifft(xin) * N if inv else np.fft.fft(xin)
        snr = 10 * np.log10(np.sum(abs(ref) ** 2) / np.sum(abs(out - ref) ** 2))
        print(f"{'IFFT' if inv else 'FFT '}: exp={exp} err={err} SNR={snr:.1f} dB")


if __name__ == '__main__':
    if len(sys.argv) >= 2 and sys.argv[1] == 'check':
        sys.exit(0 if check(sys.argv[2:]) else 1)
    elif len(sys.argv) >= 2 and sys.argv[1] == 'selftest':
        selftest()
    else:
        print(__doc__)
