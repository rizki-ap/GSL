#!/usr/bin/env python3
"""
xspec_phat_model.py
Bit-true model of verilog xspec_phat, plus the whole fixed-point GCC-PHAT
chain (pack -> fft2048 x2 -> xspec_phat -> fft2048 IFFT x3 -> peaks) and
float references, so the FPGA path can be scored against gs_det_tdoa.py.

Usage
    python3 xspec_phat_model.py check tb_xp_*.txt
        1. Z0/Z1   : RTL vs fft2048 model           (must be bit-exact)
        2. W0..W2  : RTL vs xspec_phat model        (must be bit-exact)
        3. R0..R2  : RTL vs fft2048 IFFT model      (must be bit-exact)
        4. TDOA accuracy of the RTL correlations, refined like the
           HPS / Nios would, against float GCC-PHAT and the truth.

Library
    W = xspec_phat(Z0, Z1, lo=0, hi=1024)      # Z: list of (re, im), natural order
    res = gcc_chain([x0, x1, x2, x3], max_lag)  # full fixed-point chain
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'fft2048'))
sys.path.insert(0, HERE)

from fft2048_model import fft2048, N          # noqa: E402
from gen_rsqrt_rom import rsqrt_rom           # noqa: E402

HALF = N // 2
OMAX = 32767
ROM = rsqrt_rom()
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
MASK36 = (1 << 36) - 1


# ----------------------------------------------------------------------
# Bit-true xspec_phat
# ----------------------------------------------------------------------
def half_rnd(s):
    return (s + 1) >> 1


def unpack(zk, zm):
    zr, zi = zk
    mr, mi = zm
    xa = (half_rnd(zr + mr), half_rnd(zi - mi))
    xb = (half_rnd(zi + mi), half_rnd(mr - zr))
    return xa, xb


def _mag1(v):
    return (v if v >= 0 else ~v) & MASK36


def _rnd_clamp(x, s):
    t = x + ((1 << (s - 1)) - 1) + ((x >> s) & 1)
    t >>= s
    return max(-OMAX, min(OMAX, t))


def phat(xa, xb, esum=0, floor=0):
    """G = Xb * conj(Xa), normalised to ~32767. Mirrors xspec_phat_unit.v."""
    ar, ai = xa
    br, bi = xb
    gr = br * ar + bi * ai
    gi = bi * ar - br * ai
    mor = _mag1(gr) | _mag1(gi)
    m = mor.bit_length() - 1 if mor else 0
    if (gr == 0 and gi == 0) or (m + esum < floor):
        return 0, 0
    if m >= 16:
        nr, ni = gr >> (m - 16), gi >> (m - 16)
    else:
        nr, ni = gr << (16 - m), gi << (16 - m)
    q = nr * nr + ni * ni
    if q >= (1 << 35):
        q = (1 << 35) - 1
    e17 = (q >> 34) & 1
    idx = ((q >> 26) if e17 else (q >> 24)) & 0x3FF
    y = ROM[idx]
    return _rnd_clamp(nr * y, 18 + e17), _rnd_clamp(ni * y, 18 + e17)


def xspec_phat(Z0, Z1, lo=0, hi=HALF, e0=0, e1=0, floor=0):
    """Return W0, W1, W2 as lists of (re, im) in NATURAL order
    (the RTL stores them at bit-reversed addresses)."""
    W = [[(0, 0)] * N for _ in range(3)]
    for k in range(HALF + 1):
        km = (N - k) % N
        x0, x1 = unpack(Z0[k], Z0[km])
        x2, x3 = unpack(Z1[k], Z1[km])
        X = [x0, x1, x2, x3]
        if lo <= k <= hi:
            es = [e0 + e0, e0 + e1, e0 + e1, e0 + e1, e0 + e1, e1 + e1]
            G = [phat(X[a], X[b], es[j], floor) for j, (a, b) in enumerate(PAIRS)]
        else:
            G = [(0, 0)] * 6
        for p in range(3):
            (gar, gai), (gbr, gbi) = G[2 * p], G[2 * p + 1]
            W[p][k] = (gar - gbi, gai + gbr)
            W[p][km] = (gar + gbi, gbr - gai)      # written last, like the RTL
    return W


# ----------------------------------------------------------------------
# Full fixed-point chain
# ----------------------------------------------------------------------
def coarse_peak(r, L):
    """Argmax of real correlation r (natural order) over lags -L..L."""
    best, lag = None, 0
    for n in list(range(N - L, N)) + list(range(0, L + 1)):
        if best is None or r[n] > best:
            best, lag = r[n], (n if n <= L else n - N)
    return lag


def lag_window(r, L):
    """Correlation samples for lags -L..L (what the HPS / Nios receives)."""
    return [r[(m + N) % N] for m in range(-L, L + 1)]


def refine(win, L, up=16, taps=16):
    """Windowed-sinc refinement around the coarse peak, as planned for the
    HPS / Nios: evaluate on a 1/up grid within +/-1 sample. Returns lag in
    samples (float)."""
    lags = range(-L, L + 1)
    c = max(lags, key=lambda m: win[m + L])
    best_t, best_v = float(c), -1e300
    for i in range(-up, up + 1):
        t = c + i / up
        v = 0.0
        for m in range(max(-L, c - taps), min(L, c + taps) + 1):
            d = t - m
            w = 0.5 * (1 + math.cos(math.pi * d / (taps + 1)))      # Hann
            s = 1.0 if d == 0 else math.sin(math.pi * d) / (math.pi * d)
            v += win[m + L] * s * w
        if v > best_v:
            best_t, best_v = t, v
    return best_t


def gcc_chain(x, L, lo=0, hi=HALF, floor=0):
    """x: 4 channel windows (ints, natural, zero-padded to N)."""
    z0r, z0i, e0, err0 = fft2048(x[0], x[1])
    z1r, z1i, e1, err1 = fft2048(x[2], x[3])
    Z0 = list(zip(z0r, z0i))
    Z1 = list(zip(z1r, z1i))
    W = xspec_phat(Z0, Z1, lo, hi, e0, e1, floor)
    R, er = [], []
    for p in range(3):
        rr, ri, e, err = fft2048([w[0] for w in W[p]], [w[1] for w in W[p]], inverse=True)
        R.append((rr, ri))
        er.append(e)
    corr = []
    for p in range(3):
        corr += [R[p][0], R[p][1]]
    return dict(Z0=Z0, Z1=Z1, exp_z=(e0, e1), W=W, R=R, exp_r=er, corr=corr,
                err=err0 or err1)


# ----------------------------------------------------------------------
# Float references
# ----------------------------------------------------------------------
def gcc_float_corr(xa, xb, n, lo=0, hi=None):
    """Float GCC-PHAT correlation at integer lags (numpy), length n,
    with the same band mask as xspec_phat."""
    import numpy as np
    Xa = np.fft.rfft(xa, n=n)
    Xb = np.fft.rfft(xb, n=n)
    G = Xb * np.conj(Xa)
    G /= (np.abs(G) + 1e-12)
    k = np.arange(len(G))
    G[(k < lo) | (k > (n // 2 if hi is None else hi))] = 0
    return np.fft.irfft(G, n=n)


def gcc_phat_repo(x1, x2, fs, interp_factor=16, max_lag_s=0.001):
    """Verbatim algorithm of gs_det_tdoa.gcc_phat (returns seconds)."""
    import numpy as np
    n = len(x1)
    X1 = np.fft.rfft(x1, n=n)
    X2 = np.fft.rfft(x2, n=n)
    G = X2 * np.conj(X1)
    G /= (np.abs(G) + 1e-12)
    nf = n * interp_factor
    cc = np.fft.fftshift(np.fft.irfft(G, n=nf))
    lags = np.arange(-nf // 2, nf // 2) / (fs * interp_factor)
    mask = np.abs(lags) <= max_lag_s
    idx = np.argmax(cc[mask])
    return float(lags[mask][idx])


# ----------------------------------------------------------------------
# Checker for tb_xspec_phat dumps
# ----------------------------------------------------------------------
def _read_dump(fn):
    meta = {}
    rows = []
    with open(fn) as f:
        for line in f:
            if line.startswith('#'):
                tok = line[1:].split()
                for i in range(0, len(tok) - 1, 2):
                    try:
                        meta[tok[i]] = float(tok[i + 1])
                    except ValueError:
                        meta[tok[i]] = tok[i + 1]
            elif line.strip():
                rows.append(list(map(int, line.split())))
    return meta, rows


def check(files, verbose=True):
    import numpy as np
    ok = True
    for fn in files:
        meta, rows = _read_dump(fn)
        L = int(meta['maxlag'])
        lo, hi = int(meta['lo']), int(meta['hi'])
        floor = int(meta.get('floor', 0))
        fs = meta['fs']
        wl = int(meta['len'])
        d = [meta[f'd{i}'] for i in range(4)]
        x = [[r[i] for r in rows] for i in range(4)]
        col = lambda c: [r[c] for r in rows]           # noqa: E731
        rtl_Z0 = list(zip(col(4), col(5)))
        rtl_Z1 = list(zip(col(6), col(7)))
        rtl_W = [list(zip(col(8 + 2 * p), col(9 + 2 * p))) for p in range(3)]
        rtl_R = [(col(14 + 2 * p), col(15 + 2 * p)) for p in range(3)]
        rtl_ez = (int(meta['ez0']), int(meta['ez1']))
        rtl_er = [int(meta[f'er{p}']) for p in range(3)]

        print(f"=== {fn}: {meta.get('name', '')}")
        # 1-3. bit-exact stages
        m = gcc_chain(x, L, lo, hi, floor)
        zok = (m['Z0'] == rtl_Z0 and m['Z1'] == rtl_Z1 and m['exp_z'] == rtl_ez)
        W_from_rtlZ = xspec_phat(rtl_Z0, rtl_Z1, lo, hi, rtl_ez[0], rtl_ez[1], floor)
        wbad = sum(1 for p in range(3) for k in range(N) if W_from_rtlZ[p][k] != rtl_W[p][k])
        rok = all(list(m['R'][p][0]) == rtl_R[p][0] and list(m['R'][p][1]) == rtl_R[p][1]
                  for p in range(3)) and m['exp_r'] == rtl_er
        print(f"  Z0/Z1 (fft2048)     : {'BIT-EXACT' if zok else 'MISMATCH'}  exps {rtl_ez}")
        print(f"  W0..W2 (xspec_phat) : {'BIT-EXACT' if wbad == 0 else f'MISMATCH ({wbad} bins)'}")
        print(f"  R0..R2 (IFFT)       : {'BIT-EXACT' if rok else 'MISMATCH'}  exps {tuple(rtl_er)}")
        ok &= zok and wbad == 0 and rok

        # 4. TDOA accuracy
        corr = []
        for p in range(3):
            corr += [rtl_R[p][0], rtl_R[p][1]]
        xa = [np.array(x[i][:wl], dtype=float) for i in range(4)]
        dead = [max(abs(v) for v in x[i][:wl]) <= 2 for i in range(4)]   # idle / unplugged
        if verbose:
            print("  pair   truth   coarse  fpga16x  float16x  repo16x | fixpt err  (samples; 1 sample ="
                  f" {1e6 / fs:.0f} us)")
        worst_fix = 0.0
        for j, (a, b) in enumerate(PAIRS):
            truth = d[b] - d[a]
            if dead[a] or dead[b]:
                zero = max(abs(v) for v in corr[j]) <= 2
                print(f"  ({a},{b})  dead channel -> correlation suppressed (max|r| <= 2): {zero}")
                ok &= zero
                continue
            win = lag_window(corr[j], L)
            coarse = coarse_peak(corr[j], L)
            t_fpga = refine(win, L)
            fl = gcc_float_corr(np.pad(xa[a], (0, N - wl)), np.pad(xa[b], (0, N - wl)), N, lo, hi)
            t_float = refine([fl[(m + N) % N] for m in range(-L, L + 1)], L)
            t_repo = gcc_phat_repo(xa[a], xa[b], fs, 16, L / fs) * fs
            # fixed-point error on a fine grid (separates it from the 1/16 grid)
            f_fpga = refine(win, L, up=256)
            f_float = refine([fl[(m + N) % N] for m in range(-L, L + 1)], L, up=256)
            fix = f_fpga - f_float
            worst_fix = max(worst_fix, abs(fix))
            if verbose:
                print(f"  ({a},{b}) {truth:+8.3f} {coarse:+5d} {t_fpga:+9.4f} {t_float:+9.4f}"
                      f" {t_repo:+9.4f} | {fix:+.4f} ({fix * 1e6 / fs:+.3f} us)")
        print(f"  worst fixed-point TDOA error: {worst_fix * 1e6 / fs:.3f} us (budget 1 us)")
        ok &= worst_fix * 1e6 / fs <= 1.0
    print("ALL CHECKS PASSED" if ok else "CHECKS FAILED")
    return ok


if __name__ == '__main__':
    if len(sys.argv) >= 2 and sys.argv[1] == 'check':
        sys.exit(0 if check(sys.argv[2:]) else 1)
    print(__doc__)
