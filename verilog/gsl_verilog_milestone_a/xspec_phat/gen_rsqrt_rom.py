#!/usr/bin/env python3
"""
gen_rsqrt_rom.py
Reciprocal-square-root ROM for xspec_phat (PHAT normalisation).

The PHAT unit normalises |G|^2 into q = m * 2^(2e), m in [1, 4),
and looks up y = 1/sqrt(m) with a 10-bit index idx = floor(m * 256)
(only idx = 256 .. 1023 are used).

    rom[idx] = round( 2^17 * (32767/32768) / sqrt((idx + 0.5) / 256) )

  * Midpoint sampling -> max relative error ~0.1 % (magnitude only;
    PHAT phase comes from G itself and is not affected).
  * The 32767/32768 factor makes a unit-magnitude output land at
    32767 (0.25 FS) instead of 32768, so the packed IFFT input
    W = G_a + j*G_b stays below 0.5 FS (fft2048 input contract).
  * All values < 2^17, so they fit a positive s18.

Output: rsqrt1024.hex (1024 lines, 5 hex digits).
Usage:  python3 gen_rsqrt_rom.py [outdir]
"""
import math
import os
import sys

SIZE = 1024
SCALE = (1 << 17) * 32767.0 / 32768.0


def rsqrt_rom():
    rom = [0] * SIZE
    for idx in range(256, SIZE):
        m = (idx + 0.5) / 256.0
        rom[idx] = int(math.floor(SCALE / math.sqrt(m) + 0.5))
    return rom


def main(outdir='.'):
    rom = rsqrt_rom()
    with open(os.path.join(outdir, 'rsqrt1024.hex'), 'w') as f:
        f.write('\n'.join(format(v, '05X') for v in rom) + '\n')
    print(f"Generated rsqrt1024.hex: idx 256..1023 -> {rom[256]} .. {rom[1023]}")


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else '.')
