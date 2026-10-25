# UMA-8 v2 mic geometry (center mic = origin)

Source: ODAS config for miniDSP UMA-8 (introlab/odas, config/odaslive/minidsp.cfg). The miniDSP v2 manual and product brief give no mic coordinates; the manual's mic diagram is missing from the PDF.

Layout: 1 center mic + 6 mics on a circle, r = 43 mm, 60° apart (regular hexagon). Board is 90 mm in diameter. Mics are bottom-port and fire through the board, so all z = 0 (planar).

| Mic | x (m) | y (m) | z (m) | angle from +y |
|---|---|---|---|---|
| 0 (center) | 0.0000 | 0.0000 | 0 | – |
| 1 | 0.0000 | +0.0430 | 0 | 0° |
| 2 | +0.0372 | +0.0215 | 0 | 60° |
| 3 | +0.0372 | -0.0215 | 0 | 120° |
| 4 | 0.0000 | -0.0430 | 0 | 180° |
| 5 | -0.0372 | -0.0215 | 0 | 240° |
| 6 | -0.0372 | +0.0215 | 0 | 300° |

Exact values: 37.24 mm = 43·cos30°, 21.5 mm = 43·sin30°.

Caveats:
- Numbering, channel order and rotation are ODAS's convention. ODAS's xmos.cfg uses the same 43 mm hexagon rotated 30° with a different order. Confirm the order on the real board by tapping each mic in RAW-firmware mode (8 ch; ch8 is the spare PDM port).
- Measure opposite outer mics with calipers to confirm 86 mm.

Derived for TDOA: the largest baseline is 86 mm, so the maximum TDOA is about 251 µs (about 12 samples at 48 kHz). Adjacent spacing is 43 mm, so spatial aliasing starts above about 4 kHz.
