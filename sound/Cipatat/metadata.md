# Cipatat gunshot dataset – metadata

Reference (ground-truth) values for the Cipatat test-range recordings: the theoretical values of the reference shots, as supplied with the recordings.

## Recording format

All files in this folder share the same format:

| Property | Value |
|---|---|
| Channels | 4 |
| Sample rate | 51 200 Hz |
| Sample format | 24-bit PCM (`PCM_24`) |
| Naming | `shot-NN_audio_4ch_cut.wav`, cut from the original recordings |
| Cut window | The shot window stored with each original recording (since Nov 2026; see *Original recordings*). Mostly 5 s with the SW 1.4–2.5 s in. |

## Rule: development set (Nov 2026)

**Only shots of 12.7 mm calibre and below are used for detection-algorithm development** (choosing filters, thresholds and features). The larger calibres are kept in the folder but are not used to tune the detector.

| Set | Shots | Calibres | Distance |
|---|---|---|---|
| **Development set** | 01, 02, 03, 04, 05, 06, 15, 16, 17 (+ 06a, an SW-only window of 06) | 5.56 mm to 12.7 mm | 88–764 m |
| Excluded | 07, 08, 09, 10, 11, 12, 13, 14 | 14.5 mm to 120 mm, RPG-7 | 350–1470 m |

- **Effective size:** **9 recordings** (01–06, 15, 16, 17), where shot-17 holds 5 single shots at 88 m. shot-06a is the same recording as shot-06, cut so that it ends before the MB (an SW-only test). The analyses below were made on 8 shots (01–06, 15, 16), before shot-17 was replaced.
- **In the analysis script:** `python/ana_sw_mb_range.py ..\sound\Cipatat -max_caliber_mm 12.7` analyses only the development set.

## Columns

| Column | Meaning |
|---|---|
| Shot | Reference shot ID |
| WAV file | Matching file in this folder |
| Weapon | Weapon fired |
| Caliber | Ammunition caliber / designation |
| Type | Single shot or burst |
| Rounds | Number of rounds fired |
| Azimuth (°) | Azimuth of the shooter |
| Elevation (°) | Elevation angle of the shooter |
| Distance (m) | Shooter-to-array distance |
| DCPA (m) | Distance at Closest Point of Approach: how close the bullet trajectory passes to the array (miss distance) |

## Reference values

| Shot | WAV file | Weapon | Caliber | Type | Rounds | Azimuth (°) | Elevation (°) | Distance (m) | DCPA (m) | Duration (s) |
|---|---|---|---|---|---:|---:|---:|---:|---:|---:|
| Shot_01 | shot-01_audio_4ch_cut.wav | M16 | 5.56×45 mm | Single | 1 | 269 | 0 | 306 | 20 | 5.000 |
| Shot_02 | shot-02_audio_4ch_cut.wav | AK47 | 7.62×39 mm | Single | 1 | 85 | −1 | 306 | 20 | 5.000 |
| Shot_03 | shot-03_audio_4ch_cut.wav | Dragunov | 7.62×54R mm | Single | 1 | 355 | 6 | 500 | 10 | 5.000 |
| Shot_04 | shot-04_audio_4ch_cut.wav | MINIMI | 7.62×51 mm | Burst | 5 | 1 | 1 | 304 | 30 | 5.000 |
| Shot_05 | shot-05_audio_4ch_cut.wav | AX338 / Accuracy | 8.6×70 mm | Single | 1 | 0 | 7 | 517 | 35 | 5.000 |
| Shot_06 | shot-06_audio_4ch_cut.wav | STEYR | 12.7×99 mm | Single | 1 | 11 | 2 | 764 | 30 | 5.000 |
| Shot_07 | shot-07_audio_4ch_cut.wav | Cannon | 14.5×114 mm | Single | 1 | 11 | 2 | 752 | 15 | 5.000 |
| Shot_08 | shot-08_audio_4ch_cut.wav | Cannon | 20 mm | Single | 1 | 0 | 0 | 1302 | 52 | 5.000 |
| Shot_09 | shot-09_audio_4ch_cut.wav | Cannon | 25×137 mm | Single | 1 | 356 | 1 | 1032 | 50 | 5.000 |
| Shot_10 | shot-10_audio_4ch_cut.wav | TP-T-GDOTS | 30×173 mm | Single | 1 | 1 | 0 | 1204 | 73 | 5.000 |
| Shot_11 | shot-11_audio_4ch_cut.wav | Cannon | T40 (fin-stabilized dart and target-practice rounds) | Single | 1 | 1 | 0 | 1470 | 16 | 6.000 |
| Shot_12 | shot-12_audio_4ch_cut.wav | Cannon | 105 mm HEP shell | Single | 1 | 355 | 1 | 1356 | 67 | 5.000 |
| Shot_13 | shot-13_audio_4ch_cut.wav | RPG-7 | 73 mm PG-7V | Single | 1 | 356 | 0 | 350 | – | 10.000 |
| Shot_14 | shot-14_audio_4ch_cut.wav | Cannon | 120 mm BSCC | Single | 1 | 2 | 1 | 1270 | 123 | 5.000 |
| Shot_15 | shot-15_audio_4ch_cut.wav | AK47 | 7.62×39 mm | Single | 1 | 0 | −10 | 300 | 20 | 5.000 |
| Shot_16 | shot-16_audio_4ch_cut.wav | M16 | 5.56×45 mm | Single | 1 | 175 | 0 | 304 | 24 | 5.000 |
| Multiple_Shot | shot-17_audio_4ch_cut.wav | AK47 | 7.62×39 mm | Single | 5 | 28 | 0 | 88 | 4 | 38.000 |

## Notes and open points

- **Shot_13 (RPG-7):** no DCPA given in the source (`----`).
- **Multiple_Shot → shot-17:** confirmed from the original recording, which holds five single shots about 7 s apart (see *Original recordings*).
- **shot-06a_audio_4ch_cut.wav:** has no entry in the source table. It is the same recording as shot-06, cut 0–11 s: the SW is at 10.04 s and the cut ends before the MB, so it is an SW-only case (detection without localization).
- **shot-18_audio_4ch_cut.wav:** not in the source table. It is a 34.6 s cut of the same recording as shot-17 (6.42–41.04 s of the original), so it duplicates shot-17.
- **Cuts replaced (Nov 2026):** all `shot-NN` WAVs were replaced by the shot windows stored with the original recordings. Where old and new cuts overlap the samples are identical, so the analyses below stay valid; only the event times inside the files moved (e.g. shot-03's SW is now at 2.43 s instead of 7.01 s). The shot-03 noise analysis used the earlier 15.6 s cut.
- **Unknowns to confirm:**
  - The reference direction and sign convention of *Azimuth*. The PILAR-V geometry fit below is consistent with clockwise from north.
  - ~~The array geometry and channel order~~ — now known: see "Array geometry – `PILAR-V`" below.
- **Source:** the reference spreadsheet supplied with the recordings. Numbers were normalized to a decimal point and no thousands separator.

## SW/MB range check – development set (Nov 2026)

Rerun of the analysis under the development rule above, with `python/ana_sw_mb_range.py -max_caliber_mm 12.7`. That script uses the same SW/MB classifier as the gs_det pipeline (`python/gs_det_events.py`, Stage 1 `method = classify`). Full table: `sw_mb_range_results.md` (also `.csv` and `sw_mb_range_overview.png`).

**Method**

- **Classification:** every event is classified as shockwave (SW) or muzzle blast (MB) by its own frequency content and time/duration:
  - **SW:** strong (≥ 25 dB above noise) and impulsive (crest ≥ 14 dB) in 1–4 kHz, with a near-instant front (10–90 % rise ≤ 0.08 ms).
  - **MB:** clear in **150–1000 Hz**, not dominated by high frequencies, and standing out from what came before it (≥ 11 dB above noise and ≥ 6 dB above the preceding 50–250 ms).
- **No order required:** a SW and an MB are each optional, and the MB does not have to follow a SW. When a SW is present, its echo tail is not searched for the MB.
- **Range:** solved from the SW→MB delay Dt with a drag-aware bullet model, v(s) = v0·e^(−s/λ), using the reference DCPA as the miss distance and c = 347 m/s (about 27 °C).

**Muzzle-blast band, re-chosen on the development set:**

- **Sweep:** 68 single bands were tested on the 8 effective shots.
- **Choice: 150–1000 Hz.** The weakest muzzle blast still stands +8.2 dB above any other sound after the SW echo tail (median +14.9 dB). The weakest is shot-03, at 14 dB above noise.
- **Before, with all calibres: 100–1000 Hz.** Its weakest muzzle blast stands only +5.9 dB clear, because wind noise below 150 Hz hurts shot-03.
- **Why the low edge could move up:** the cannons needed 100–150 Hz, but without them the low edge can rise. The best bands all start at 150–300 Hz and end at 600–1200 Hz; the top few differ by only about 0.4 dB.

**Shockwave settings, re-tuned on the development set:**

| Feature | Real shockwaves | Muzzle blasts |
|---|---|---|
| 1–4 kHz level above noise | ≥ 37 dB | (not a discriminator) |
| Impulsiveness (crest) | ≥ 16.6 dB | ≤ 11.9 dB where the front is ≤ 0.08 ms (shot-04's muzzle blast) |
| 10–90 % rise of the front | ≤ 0.059 ms | (some are also sharp) |

- **Crest minimum raised from 12 to 14 dB.** This leaves at least 2 dB of margin on both sides. With 12 dB, shot-02's muzzle blast (crest 11.6 dB) was only 0.4 dB from being called a shockwave.
- **"Short low-frequency body" alternative removed** (`sw_max_lf_width_ms = 0`). Only the 105/120 mm cannons needed it; every development shockwave passes on its sharp front.
- **Band kept at 1–4 kHz.** 1.5–6 kHz and 2–8 kHz raise the shockwave crest by about 1.5 dB, but they mistime shot-04's muzzle blast.
- **Effect:** detection results and ranges are unchanged, and fewer echoes are labelled as shockwaves (15 events instead of about 25).

**Result:** SW and MB found in every development shot.

| Shot | Calibre | Computed / reference (m) | Error |
|---|---|---|---:|
| 01 M16 | 5.56 mm | 303 / 306 | −0.9 % |
| 02 AK47 | 7.62×39 | 324 / 306 | +5.8 % |
| 03 Dragunov | 7.62×54R | 448 / 500 | −10.3 % |
| 04 MINIMI (burst) | 7.62×51 | 341 / 304 | +12.1 % |
| 05 AX338 | 8.6×70 | 487 / 517 | −5.8 % |
| 06 STEYR | 12.7×99 | 744 / 764 | −2.6 % |
| 15 AK47 | 7.62×39 | 297 / 300 | −1.1 % |
| 16 M16 | 5.56 mm | 298 / 304 | −1.9 % |

6 of the 8 shots are within ±5.8 %. The ranges are within 1 m of the earlier 100–1000 Hz run.

**Detector tests on cut recordings (development shots only):**

- **SW removed (8 cuts):** all found as MB only, at the right time. The exception is shot-04, where the cut still contained a later round of the burst, so SW+MB was correct.
- **MB removed (8 cuts):** all 8 found as SW only, with no false MB.

**Notes**

- **Bias:** most errors are negative. Either the air was warmer than 27 °C, or the bullets were slower than the table assumes. Use the measured temperature of the test day if it is known.
- **shot-03 (Dragunov):** fits about 770 m/s against 500 m, but the acceptance window stored with the original recording points to about 431 m (see *Original recordings* below). Against 431 m the computed range is within 2 %. It has the weakest muzzle blast in the set, because its background noise in the MB band is 15–46 dB louder than in the other shots, with engine-like tones inside the MB band (see *Background noise – shot-03* below).
- **shot-04 (MINIMI burst):** the 5 rounds are about 75 ms apart, and all 5 pairs give the same delay.
- **shot-17 (Multiple_Shot):** the WAV in this folder does not fit 88 m with DCPA 4 m, and it turns out not to come from the Multiple_Shot recording at all. The real recording fits well (+8 %); see *Original recordings* below.

**Earlier analysis with all calibres (before the rule):** 12 of the 15 usable shots were within ±5.8 %. Findings on the excluded shots:

- **shot-14 (120 BSCC):** fits a projectile of about 1180 m/s, not 1700 m/s.
- **shot-11 (T40):** fits about 1300 m/s, not 1500 m/s.
- **shot-13 (RPG-7):** subsonic, so it was correctly found as MB only.

## Array geometry – `PILAR-V` (Nov 2026)

The recordings come from PILAR-V acoustic sensors ("Trackball"). At least 12 different units were used (see *Original recordings* below), but every recording header lists the same geometry and sensitivity, so the geometry below is the nominal design. In the gs_det / gs_gen configs it is `mic_config = PILAR-V`; details are in `devices/mic_pilar-v.md`.

| WAV channel | Sensor channel | x (m) | y (m) | z (m) | Position |
|---|---|---:|---:|---:|---|
| 0 | 1 | 0.2055 | 0.000 | 0.000 | Ring, 0° |
| 1 | 2 | −0.1028 | 0.178 | 0.000 | Ring, 120° |
| 2 | 3 | −0.1028 | −0.178 | 0.000 | Ring, 240° |
| 3 | 4 | 0.000 | 0.000 | 0.2412 | Top |

- **Shape:** a ring of 3 mics (radius 0.2055 m, 0.356 m apart) plus 1 mic 0.2412 m above the centre (0.317 m from each ring mic).
- **Maximum time difference:** about 1.04 ms.
- **Sensitivity:** 8.8 mV/Pa per mic.
- **Orientation, from the fit below:** x = north (heading 1.1°), y = west, z = up. Compass azimuth = −atan2(y, x) of the direction vector.
- **Header field `param_float_0x44 = 31.3`:** the same in every recording, from every unit and year, so it is **not** the air temperature (an earlier version of this file assumed it was). The air temperature is a separate per-recording value; see *Original recordings* below.

## Direction and range check – development set (Nov 2026)

Uses the SW/MB times from the detector (settings above), the `PILAR-V` geometry and 8 shots: 01–06, 15, 16.

### Time differences between channels

- **GCC-PHAT (the pipeline's Stage 2) fails on these files:** it returns about 0 µs for almost every pair. Background noise is common to all channels on several shots (correlation 0.8–0.9), and PHAT whitening locks onto that zero-delay noise.
- **Plain cross-correlation is used here instead:** it is applied to the band-passed event (SW 1–8 kHz, 4.5 ms window; MB 150–1000 Hz, 15 ms window). The delays it gives (up to about ±1 ms) agree with per-channel onset times to within about 20–50 µs.

### Azimuth and elevation of the shooter (from the muzzle blast)

The array geometry is not fitted; only the array heading (1.1°) is, as one value for all shots. c = 349.8 m/s (at the time, `param_float_0x44` was taken to be the temperature; see above).

| Shot | Reference az / el (°) | Estimated az / el (°) | Error az | Error el |
|---|---|---|---:|---:|
| 01 | 269 / 0 | 272.8 / −1.8 | +3.8° | −1.8° |
| 02 | 85 / −1 | 84.4 / 0.8 | −0.6° | +1.8° |
| 03 | 355 / 6 | 355.7 / 6.0 | +0.7° | 0.0° |
| 04 | 1 / 1 | 358.7 / 2.4 | −2.3° | +1.4° |
| 05 | 0 / 7 | 1.7 / 7.7 | +1.7° | +0.7° |
| 06 | 11 / 2 | 9.8 / 0.1 | −1.2° | −1.9° |
| 15 | 0 / −10 | 355.9 / −5.9 | −4.1° | +4.1° |
| 16 | 175 / 0 | 177.0 / 1.2 | +2.0° | +1.2° |

**RMS error: azimuth 2.4°, elevation 2.0°.** Without the heading correction (x exactly north), the azimuth RMS is 2.6°. Five of the 8 shots come from about 0° azimuth.

### Range

| Method | RMS error | Mean | Worst |
|---|---:|---:|---:|
| Drag-aware model, SW→MB delay + **reference** DCPA, c = 347 m/s (`ana_sw_mb_range.py`) | 6.4 % | −0.6 % | 12.1 % (shot-04) |
| Acoustic only, PILAR constant-speed formula d = c·Δt / (1 − m·s) (m = MB direction, s = SW direction) | 34 % | +30 % | 50 % (shot-06) |

- **The constant-speed formula overestimates the range:** the bullets slow down, as predicted in `python/drag_model.md`.
- **A range that uses only acoustic measurements needs a drag-aware solver.**
- **With c = 349.8 m/s (31.3 °C, now known not to be a temperature)**, the drag-aware errors are: 01 +0.6 %, 02 +8.3 %, 03 −8.5 %, 04 +14.0 %, 06 −1.1 %, 15 +1.2 %, 16 −0.5 %.

## Background noise – shot-03 (Nov 2026)

Measured on 13.7 s of noise-only signal: the whole file except 6.9–8.8 s (SW, MB and their echoes). Levels are in dBFS, because the recorder's full-scale pressure is unknown. Figure: `shot-03_noise_characteristics.png` in this folder.

**Much louder than the rest of the set.** The shockwave peak is normal (−25.5 dBFS), but the noise is not:

| Shot | Noise 10–150 Hz | Noise 150–1000 Hz | Noise 1–4 kHz | SW peak above 1–4 kHz noise | MB SNR (detector) |
|---|---:|---:|---:|---:|---:|
| 03 | −63.4 | −72.0 | −79.3 | 54 dB | 13.6 dB |
| Other development shots | −105 to −67 | −118 to −87 | −119 to −95 | 75–97 dB | 21–54 dB |

**Spectrum.** The noise falls steadily with frequency.

| Band | Share of noise power | Spectral slope |
|---|---:|---|
| Below 20 Hz | 71 % | – |
| 20–150 Hz | 20 % | f^−1.9 |
| 150–1000 Hz (MB band) | 5 % | f^−1.5 |
| 1–10 kHz | 1.2 % | f^−2.4 |
| 10–25.6 kHz | 0.02 % | steep roll-off |

**Three sources make up the noise:**

- **Wind, below about 50 Hz.**
  - It carries most of the power and comes in gusts: peaks up to 12 dB above the median level, and a kurtosis of 5 in 10–150 Hz.
  - The microphones barely agree on it (coherence < 0.2), as expected for wind turbulence at each microphone rather than a sound wave.
  - The second half of the file is 2.4 dB noisier than the first.
- **Machinery tones, 80–850 Hz.** These are steady and come from fixed directions, below the horizon. They look like engines or generators.
  - **A ~83.6 Hz series** (83.6, 167, 251 Hz) arrives from azimuth ≈ 75°, elevation ≈ −40°.
  - **A ~412 Hz line, with its harmonic near 820 Hz,** arrives from azimuth ≈ 22°, elevation ≈ −25°, very consistently in every 1 s window.
  - The mics agree strongly on these tones: coherence 0.6–0.87 at 80–300 Hz. They fall **inside the MB band**, so they cause shot-03's weak muzzle blast. The shooter is at 355°, +6°.
  - A weak trace of the 412 Hz line is also in shot-02.
- **Insects, 3.8 and 4.15–4.2 kHz.**
  - Narrow tones that pulse about 3 times a second. Only present at 0–6 s and from 14 s on.
  - They come from azimuth ≈ 115°, and the mics agree on them (coherence up to 0.87 at 4.2 kHz).
  - They fall inside the SW band (1–4 kHz) but are 50 dB below the shockwave, so they do no harm here.
- **No mains hum:** there are no 50 Hz lines.

**What this means:**

- **Denoising before detection:** remove the steady tones in the MB band with notch filters, an adaptive line canceller or spectral subtraction. This is the *Denoising* block of the PILAR-V chain.
- **Spatial filtering:** because the tones are coherent and come from fixed directions, a null steered towards 22° and 75° would also suppress them.
- **Simulation:** the generator's default `noise_slope = −2.4` matches 1–10 kHz, but this noise is flatter in the MB band (f^−1.5). It is also partly the same at every microphone and contains tones, while the generator adds independent random noise to each microphone.

## Sound levels – estimated calibration (Nov 2026)

The WAV files store values relative to full scale (FS). The sensor header gives the microphone sensitivity (8.8 mV/Pa) but not the recorder's full-scale voltage, so FS was **estimated from physics**: the measured shockwave peaks were compared with the Whitham N-wave theory for each bullet at its reference miss distance (DCPA).

**Method.** Peak overpressure p = 0.53·p0·(M²−1)^(1/8)·d / (b^(3/4)·l^(1/4)) and N-wave duration T = 1.82·M·b^(1/4)·d / (c·(M²−1)^(3/8)·l^(1/4)).

- d, l: bullet diameter and length; b: DCPA; M: bullet Mach number where the shockwave was emitted (drag model of `ana_sw_mb_range.py`).
- p0 = 98 kPa and c = 347 m/s.
- Peaks are measured after a 100 Hz high-pass and averaged over the 4 channels.

| Shot | Predicted SW peak | Predicted / measured N-wave duration | Measured SW peak | Pa per FS |
|---|---:|---:|---:|---:|
| 01 | 92 Pa | 215 / 195 µs | 0.066 FS | 1396 |
| 02 | 110 Pa | 301 / 283 µs | 0.078 FS | 1404 |
| 03 | 181 Pa | 238 / 273 µs | 0.046 FS | 3961 |
| 04 | 87 Pa | 310 / 195 µs | 0.105 FS | 820 |
| 05 | 80 Pa | 326 / 332 µs | 0.034 FS | 2362 |
| 06 | 121 Pa | 429 / 410 µs | 0.078 FS | 1557 |
| 15 | 110 Pa | 300 / 293 µs | 0.067 FS | 1659 |
| 16 | 80 Pa | 225 / 195 µs | 0.053 FS | 1524 |

**Result: 1 FS ≈ 1520 Pa (157.6 dB SPL peak), uncertain by about ±3 dB.**

- **Five shots agree within ±1 dB:** 01, 02, 06, 15 and 16. Their measured N-wave durations also match the theory within 10 %, which confirms the bullet and DCPA values.
- **Three shots disagree, but their durations show why:**
  - **shot-03:** the N-wave is longer than predicted, so the miss distance was probably larger than 10 m.
  - **shot-04:** the N-wave is much shorter than predicted, so the miss distance was probably much smaller than 30 m.
  - **shot-05:** the duration matches, so the reason is unclear: possibly a different bullet.
- **Implied full-scale voltage:** 1520 Pa × 8.8 mV/Pa ≈ 13.4 V. A ±10 V recorder input would give 1136 Pa (155.1 dB), 2.5 dB from this estimate. **The recorder's real input range would confirm or correct this.**

**Levels with this calibration** (peaks are averages over the 4 channels; MB after a 20 Hz high-pass):

| Shot | SW peak | MB peak | MB extrapolated to 1 m* | Background noise (20 Hz–20 kHz RMS) |
|---|---:|---:|---:|---:|
| 01 M16 | 134 dB (101 Pa) | 98 dB | 148 dB | 90 dB |
| 02 AK47 | 136 dB (119 Pa) | 109 dB | 159 dB | 77 dB |
| 03 Dragunov | 131 dB (68 Pa) | 105 dB | 159 dB | 93 dB |
| 04 MINIMI | 138 dB (160 Pa) | 105 dB | 154 dB | 82 dB |
| 05 AX338 | 128 dB (52 Pa) | 85 dB | 139 dB | 52 dB |
| 06 STEYR | 135 dB (118 Pa) | 92 dB | 149 dB | 57 dB |
| 15 AK47 | 134 dB (102 Pa) | 103 dB | 152 dB | 60 dB |
| 16 M16 | 132 dB (82 Pa) | 98 dB | 147 dB | 89 dB |

\*Spherical spreading only (+20·log10 of the distance), with no correction for ground or air absorption, so these are lower bounds.

- **Plausibility check:** rifle muzzle blasts are typically 155–165 dB at 1 m. The values here are similar or a few dB lower, as expected without the ground and absorption correction. The low values for shot-05 and shot-06 may come from muzzle brakes, which send the blast sideways.
- **shot-03 noise:** 94 dB in 10–150 Hz (wind), 86 dB in 150–1000 Hz (machine tones), 78 dB in 1–4 kHz. A level of 86 dB in the MB band suggests a running engine or generator quite close to the sensor.

## Original recordings (Nov 2026)

The WAV files in this folder are cut from longer original recordings. The cut samples are identical to the originals. The originals also store, per recording, the sensor unit, a clock time, the shot window and an acceptance window around the reference values.

| Shot | Sensor serial | Recorder clock (UTC) | Air temp. (°C)* | Original length (s) | Shot window (s) |
|---|---|---|---:|---:|---|
| 01 | 10345 | 2018-10-18 13:27 | 26 | 14 | 7–12 |
| 02 | 10100 | 2014-11-20 15:32 | 9 | 20 | 10–15 |
| 03 | (blank) | – | 20 | 43 | 32–37 |
| 04 | 00000 | – | 14 | 10 | 5–10 |
| 05 | 10118 | 2015-09-08 14:05 | 27 | 19 | 8–13 |
| 06 (and 06a) | 10100 | 2015-09-09 13:47 | 35 | 18 | 8–13 (06a: 0–11) |
| 07 | 10118 | 2015-09-09 13:23 | 35 | 15 | 2–7 |
| 08 | 1407004 | 2021-11-18 13:04 | 11 | 15 | 8–13 |
| 09 | 10204 | – | 27 | 11 | 3–8 |
| 10 | 10328 | 2021-06-02 11:56 | 27 | 24 | 18–23 |
| 11 | 10499 | 2019-07-22 13:01 | 38 | 23 | 14–20 |
| 12 | 10302 | 2021-06-01 08:39 | 22 | 42 | 27–32 |
| 13 | 10293 | 2018-01-23 09:41 | 2 | 181 | 64–74 |
| 14 | 10319 | 2021-06-09 09:43 | 30 | 23 | 15–20 |
| 15 | 10100 | – | 33 | 30 | 10–15 |
| 16 | 00000 | 2018-10-18 14:32 | 24 | 18 | 10–15 |
| 17 (Multiple_Shot) | 10468 | 2024-07-06 12:00 | 26 | 46 | 5–43 |
| Simultaneous_Shot (no WAV here) | 10468 | 2025-02-23 10:44 | 22 | 8 | 2–9 |

\*A per-recording value stored with every 1 s audio block. It is constant within a recording and follows the season (2 °C in January, 35–38 °C in July–September), so it is taken to be the air temperature in °C. This is an interpretation, not documented.

**What this changes:**

- **Not one campaign.** The set was compiled from at least 12 sensor units over 2014–2025, so conditions (site, weather, unit) differ from shot to shot. Shots 01 and 16 were recorded on the same day.
- **Sample-clock drift.** The block time-stamps advance by 51 200 ± 15 samples per block in some recordings (02, 04, 15). The audio itself is continuous, so this only matters for absolute timing.
- **Speed of sound per shot.** Using c = 331.3 + 0.606·T from each recording's temperature instead of one value (347 m/s), the drag-aware range errors of the development set become:

| Shot | T (°C) | c (m/s) | Error, c = 347 m/s | Error, per-shot c |
|---|---:|---:|---:|---:|
| 01 | 26 | 347.1 | −0.9 % | −0.9 % |
| 02 | 9 | 336.8 | +5.7 % | −2.7 % |
| 03 | 20 | 343.4 | −10.3 % | −12.6 % (+1.4 % against 431 m) |
| 04 | 14 | 339.8 | +12.1 % | +7.6 % |
| 05 | 27 | 347.7 | −5.8 % | −5.4 % |
| 06 | 35 | 352.5 | −2.7 % | +0.5 % |
| 15 | 33 | 351.3 | −1.1 % | +2.4 % |
| 16 | 24 | 345.8 | −1.9 % | −2.5 % |
| **RMS** | | | **6.4 %** | **5.8 %**, or **3.7 %** with shot-03 against 431 m |

**Acceptance windows.** For most shots the stored window is the reference ±10 % in range, ±2° in azimuth and ±3° in elevation. Two shots differ:

- **shot-03:** 388–474 m, which is ±10 % around **431 m**, not the 500 m in the reference table.
- **shot-12:** 1325–1387 m (±31 m, about ±2.3 %).

**shot-17 (Multiple_Shot), redone from the original recording.** The original holds **5 single shots, about 7 s apart** (6.92, 13.97, 21.05, 28.38 and 35.92 s), as the reference row says (single, 5 rounds). Each shot has a strong SW followed by an MB 122–126 ms later. The previous `shot-17_audio_4ch_cut.wav` did not match any part of this recording (correlation ≤ 0.36) and has been replaced by the 5–43 s window of the original (shots at 1.92, 8.97, 16.05, 23.38 and 30.92 s in the new file).

| Shot in recording | SW→MB delay | Range (drag model, DCPA 4 m, c = 347.1 m/s) |
|---|---:|---:|
| 1 | 122.4 ms | 93.6 m |
| 2 | 123.1 ms | 94.1 m |
| 3 | 123.5 ms | 94.4 m |
| 4 | 125.5 ms | 96.0 m |
| 5 | 124.5 ms | 95.2 m |

- **Result:** mean 94.7 m against the reference 88 m (+7.6 %). The expected delay at 88 m is 114.9 ms.
- **Detector:** at 88 m the MB arrives inside the SW echo period, which the classifier skips by default (`mb_in_tail = false`). So the MB here was taken as the strongest low-frequency event (MB band clearly above 1–4 kHz) 50–300 ms after the SW. This is the same limitation seen at short range in the Cadre dataset.

**Simultaneous_Shot.** An extra 8 s recording with no WAV in this folder and no reference values. It contains two SWs 34 ms apart (5.126 and 5.160 s), a weaker third at 5.186 s, and an MB-band event at 5.588 s.

### Range check on the new cuts

`ana_sw_mb_range.py -max_caliber_mm 12.7` on the replaced cuts gives the same delays as before for 01, 02, 03, 05, 06, 15 and 16. Two differences:

- **shot-04 (MINIMI burst):** the first round's MB now has a high-frequency excess of 11.1 dB, just over the `mb_max_hf_excess_db = 10` limit (8.8 dB in the old cut; the longer lead-in gives a lower high-frequency noise estimate). The detector then pairs the SW with the second round's MB (+75 ms), giving 403 m (+32 %). With the limit at 12 dB, all 8 shots give their earlier delays again.
- **shot-06a and shot-17:** SW only. shot-06a ends before the MB by design. In shot-17 the MB falls inside the SW echo period (see above).


## Ground reflection (Nov 2026)

Measured with `python/ana_ground_reflection.py`: the clean twin's direct N-wave (`../Cipatat_clean`) is subtracted from the recording, and what follows 0.35–4.5 ms after the direct wave is matched against the same N-wave shape. A channel counts as a detection when the match correlation is at least 0.6.

| Shot | Detected on | Delay after the direct SW (ms) | Reflected / direct | Array height (m) | Grazing angle |
|---|---|---|---:|---:|---:|
| 17 (5 shots) | ch 1, 3, 4 | ring 2.5–3.0, top 0.45–0.5 later | 0.26–0.28 | 1.29–1.34 | 20–23° |
| 03 | all 4 | 2.76–2.92 | 0.33–0.54 | 5.4 | 5° |
| 01, 02, 04 | 1–2 channels | 0.7–1.4 | 0.18–0.32 | – | – |
| 05, 06, 15, 16 | none clear | – | – | – | – |

- **Array height without knowing the trajectory:** for a flat ground the delay is 2 (h + z) sin φ / c, so the top mic (z = 0.2412 m) hears the reflection later than the ring. The ratio gives h = 0.2412 · d_ring / (d_top − d_ring).
- **shot-17 confirms the image-source model:** with the array 1.2–1.3 m above ground and the bullet passing 4.3 m to the side and about 1 m above the array, the generator reproduces both the direct SW time differences (8 µs RMS) and the reflection delays of every channel (within about 0.05 ms). The direct+reflected waveform then matches the recording slightly better than the free-field one (r 0.76 → 0.78–0.83 over 5 ms).
- **Why most shots show no separate reflection:** when the shock front meets the ground at a shallow angle (bullet passing level with the array), the delay 2h sin φ / c falls below about 0.4 ms. The reflection then overlaps the direct N-wave and cannot be separated by this method.
- **shot-03:** the delay points to a sensor about 5 m above the reflecting ground (or another reflecting surface), on a different site from shot-17.
- **Weaker than earlier estimates.** An earlier note gave reflections of 0.34–1.0× the SW amplitude; those were raw peaks after the SW, which include other echoes and noise. Matching the N-wave shape gives 0.2–0.3 on the clear cases.
- **Generator settings from this:** `[ground] array_height_m = 1.3`, `reflection_coeff_sw = 0.33` (reflected/direct 0.27 divided by the model's spreading factor (b/b_reflected)^0.75 ≈ 0.8). The blast reflection cannot be measured here (it merges with the direct blast); `reflection_coeff_mb = 0.9` is a typical grass value.
