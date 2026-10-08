#!/usr/bin/env python3
"""
gen_twiddle2048.py
Generate twiddle ROM init files for fft2048.

    W_N^k = cos(2*pi*k/N) - j*sin(2*pi*k/N),   k = 0 .. N/2-1

Outputs (one 5-digit hex word per line, 18-bit two's complement, Q1.17):
    tw2048_re.hex :  round( cos(2*pi*k/N) * 2^17)
    tw2048_im.hex :  round(-sin(2*pi*k/N) * 2^17)

Notes
  * k = 0 (W = 1.0) is not representable in Q1.17 (max +0.99999). The engine
    bypasses the multiplier for k = 0, so the ROM value there is unused.
  * k = N/4 gives -sin = -1.0 = -131072, which IS representable, so the
    W = -j twiddle is exact.
  * Rounding is round-half-up (floor(x + 0.5)), not Python's banker's round(),
    so the bit-true model and the ROM always agree.

Usage:  python3 gen_twiddle2048.py [outdir]
"""
import math
import os
import sys

N = 2048
FRAC = 17
BITS = 18


def q(x):
    v = math.floor(x * (1 << FRAC) + 0.5)
    return max(-(1 << (BITS - 1)), min((1 << (BITS - 1)) - 1, v))


def twiddles(n=N):
    """List of (w_re, w_im) integers for k = 0 .. n/2-1."""
    return [(q(math.cos(2 * math.pi * k / n)), q(-math.sin(2 * math.pi * k / n)))
            for k in range(n // 2)]


def hex18(v):
    return format(v & ((1 << BITS) - 1), '05X')


def main(outdir='.'):
    tw = twiddles()
    with open(os.path.join(outdir, 'tw2048_re.hex'), 'w') as f:
        f.write('\n'.join(hex18(r) for r, _ in tw) + '\n')
    with open(os.path.join(outdir, 'tw2048_im.hex'), 'w') as f:
        f.write('\n'.join(hex18(i) for _, i in tw) + '\n')
    print(f"Generated tw2048_re.hex / tw2048_im.hex ({len(tw)} entries, Q1.17)")
    for k in (0, 1, N // 8, N // 4, N // 2 - 1):
        print(f"  k={k:4d}: re={tw[k][0]:+7d}  im={tw[k][1]:+7d}")


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '.')
