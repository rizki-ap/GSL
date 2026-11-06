# UMA-8 v2 mic geometry (center mic = origin)

**Radius in use: r = 42.5 mm.** ODAS gives 43 mm, rounded to the nearest mm, so the true value is somewhere in 42.5–43.5 mm. A ruler measurement on the real board gave 42–42.5 mm. 42.5 mm is consistent with both. Confirm with calipers: opposite ring mics should be about 85 mm apart.

Source for layout and numbering: ODAS config for miniDSP UMA-8 (introlab/odas, config/odaslive/minidsp.cfg). The miniDSP v2 manual and product brief give no mic coordinates; the manual's mic diagram is missing from the PDF.

Layout: 1 center mic + 6 mics on a circle, 60° apart (regular hexagon). Board is 90 mm in diameter. Mics are bottom-port and fire through the board, so all z = 0 (planar).

| Mic | x (m) | y (m) | z (m) | angle from +y |
|---|---|---|---|---|
| 0 (center) | 0.0000 | 0.0000 | 0 | – |
| 1 | 0.0000 | +0.0425 | 0 | 0° |
| 2 | +0.0368 | +0.02125 | 0 | 60° |
| 3 | +0.0368 | -0.02125 | 0 | 120° |
| 4 | 0.0000 | -0.0425 | 0 | 180° |
| 5 | -0.0368 | -0.02125 | 0 | 240° |
| 6 | -0.0368 | +0.02125 | 0 | 300° |

Exact values: 36.806 mm = 42.5·cos30°, 21.25 mm = 42.5·sin30°.

## Configuration used in GSL (decision, Oct–Nov 2026)

`-mic_config=UMA8` in the GSL Python pipeline uses 4 mics: #0 (centre), #1, #3 and #5.

- The three ring mics form an equilateral triangle with 73.6 mm sides; the aperture is 73.6 mm.
- Coordinates: 0,0,0; 0,0.0425,0; 0.036806,-0.02125,0; -0.036806,-0.02125,0
- `channel_map = 0,1,3,5` for 8-channel RAW recordings.
- The array is planar. Azimuth is fully determined. Elevation sign comes from `planar_source_side`, and elevation is poorly conditioned near the horizon.

## Caveats

- Numbering, channel order and rotation follow ODAS's convention. ODAS's xmos.cfg uses the same hexagon rotated 30° with a different order. Confirm the order on the real board by tapping each mic in RAW-firmware mode (8 ch; ch8 is the spare PDM port).
- A uniform radius error does not change azimuth. It scales the in-plane slowness, so it mainly biases elevation (about 1% radius error gives roughly 8° of elevation error near the horizon).

## Derived for TDOA

- 4-mic configuration: the largest baseline is 73.6 mm, so the maximum TDOA is about 215 µs.
- All 7 mics: the largest baseline is 85 mm, so the maximum TDOA is about 248 µs.
- Adjacent ring spacing is 42.5 mm, so spatial aliasing starts above about 4 kHz.

Reference : https://www.minidsp.com/images/documents/UMA-8%20v2%20User%20manual.pdf
