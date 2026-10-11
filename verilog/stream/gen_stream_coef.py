#!/usr/bin/env python3
"""
gen_stream_coef.py
Coefficient ROMs for the Milestone B stream path.

  fe_fdelay.hex  fe_cond fractional-delay trim: 9-tap windowed sinc,
                 phases f = p/128, p = 0..128 (negative trims use the
                 mirrored taps). Q1.16 in s18 (1.0 = 65536).
                 Address = p*9 + tap  (1161 entries).
                 p = 0 is an exact unit impulse at tap 4 (transparent),
                 p = 128 an exact one-sample shift.
  decim63.hex    mb_decim: 63-tap linear-phase low-pass (Kaiser beta 6,
                 fc 7 kHz @ 100 kHz). Q1.17 in s18, DC gain exactly 1.
                 Passband 0-4 kHz ripple 0.02 dB, stopband >= 16 kHz
                 below -71 dB (spec: 0.1 dB / -60 dB).

Usage: python3 gen_stream_coef.py [outdir]
"""
import math
import os
import sys

import numpy as np

FD_TAPS, FD_CENTER, FD_PHASES = 9, 4, 128
DEC_TAPS = 63


def hex18(v):
    return format(int(v) & 0x3FFFF, '05X')


def fdelay_table():
    rows = []
    for p in range(FD_PHASES + 1):
        f = p / FD_PHASES
        k = np.arange(FD_TAPS) - FD_CENTER - f
        w = 0.5 * (1 + np.cos(np.pi * k / (FD_CENTER + 1.5)))
        h = np.sinc(k) * w
        h = h / h.sum()
        q = np.round(h * 65536).astype(int)
        q[int(np.argmax(np.abs(q)))] += 65536 - int(q.sum())     # exact DC gain
        rows.append([int(v) for v in q])
    return rows


def decim_taps():
    from scipy.signal import firwin
    h = firwin(DEC_TAPS, 7000, window=('kaiser', 6.0), fs=100000)
    q = np.round(h / h.sum() * (1 << 17)).astype(int)
    q[DEC_TAPS // 2] += (1 << 17) - int(q.sum())
    return [int(v) for v in q]


def main(outdir='.'):
    fd = fdelay_table()
    with open(os.path.join(outdir, 'fe_fdelay.hex'), 'w') as f:
        for row in fd:
            for v in row:
                f.write(hex18(v) + '\n')
    dec = decim_taps()
    with open(os.path.join(outdir, 'decim63.hex'), 'w') as f:
        for v in dec:
            f.write(hex18(v) + '\n')
    print(f"fe_fdelay.hex: {len(fd)} phases x {FD_TAPS} taps; p=0 -> {fd[0]}")
    print(f"decim63.hex  : {len(dec)} taps, sum {sum(dec)}")


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '.')
