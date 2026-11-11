# Cipatat clean dataset – metadata

Clean twins of the Cipatat development recordings (`../Cipatat`). Each file holds **only the shockwave (SW) and the muzzle blast (MB)**: no background noise, no reverberation, no ground reflection, no echoes. They are synthesised with the generator's Stage 1 (`python/gs_gen_clean_signal.py`) and fitted to the matching real recording. Made with `python/util_gen_cipatat_clean.py` (Nov 2026).

Real recordings cannot simply be denoised into clean ones: noise and reflections overlap the events in time and frequency. So these files are model signals placed where the real events are.

## Format

| Property | Value |
|---|---|
| Files | `shot-NN_audio_4ch_clean.wav` + `shot-NN_audio_4ch_clean.json` (fit, measured and model values, ground truth) |
| Channels, rate, format | 4 ch, 51 200 Hz, `PCM_24`: the same as the real cuts |
| Length | the same as the matching real cut, sample for sample |
| Scale | Cipatat full-scale units: 1 FS = 1520 Pa (the estimate in `../Cipatat/metadata.md`) |
| Synthesis | at 512 kHz, then low-pass filtered and decimated, like the anti-aliasing filter of a real ADC |
| Also | `clean_summary.csv` (one line per event), `clean_vs_real_overview.png` |

## Shots

| Clean file | From | Events |
|---|---|---|
| shot-01, 02, 03, 05, 06, 15, 16 | the cut with the same number | 1 SW + 1 MB |
| shot-04 | MINIMI burst | 5 rounds, each SW + MB |
| shot-06a | SW-only window of shot-06 | 1 SW. The MB would arrive after the end of the file, as in the real cut |
| shot-17 | Multiple_Shot | 5 single shots, each SW + MB |

- **Not made:** 07–14 (14.5 mm to 120 mm cannons and RPG-7). They are outside the development rule and the generator has no model for these projectiles.

## How each file is made

1. **Scene** from the reference table: calibre (generator bullet library, exponential drag), azimuth, elevation, distance, miss distance, and the air temperature of the recording (speed of sound and air absorption at 50 % humidity). PILAR-V geometry, x = north, y = west, z = up.
2. **Bullet side:** the trajectory passes the array on the left or the right at the miss distance. The side whose SW time differences match the recording is used: it fits within 26–151 µs RMS, the other side is off by 540–1320 µs.
3. **Level and duration fit:** for every event, SW peak and N-wave duration, MB peak and positive-phase duration are measured on the recording and on the clean signal with the same estimators. The generator's `amp_scale` and `dur_scale` are adjusted until they agree. Each round of a burst is fitted on its own.
4. **Timing:** each event is placed at its arrival on each channel in the recording: first one shift for the whole array, then a per-channel correction of at most ±0.2 ms, both by cross-correlation. So the **time differences between channels are those of the recording**, not of the reference scene. A channel without a clear optimum keeps the array shift. With `--align array` the model's time differences are kept instead.

## Fit quality

Pearson r is the correlation of the clean and the real signal on aligned windows (SW 1 ms after a 100 Hz high-pass, MB 12 ms after a 50 Hz high-pass), averaged over the 4 channels. Real peaks and durations are measured on the recording; model values are the generator's before the fit. Peaks are averages over the channels.

| Shot | SW peak, real / model (Pa) | N-wave duration, real / model (µs) | SW r | MB peak, real / model (Pa) | MB positive phase, real / model (ms) | MB r | SW→MB delay, real / model (ms) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 01 | 106 / 93 | 248 / 210 | 0.93 | 2.2 / 35 | 2.1 / 0.53 | 0.72 | 443 / 450 |
| 02 | 132 / 115 | 339 / 307 | 0.84 | 5.5 / 35 | 0.7 / 0.53 | 0.52 | 344 / 354 |
| 03 | 74 / 189 | 325 / 240 | 0.89 | 3.2 / 21 | 1.4 / 0.53 | 0.31 | 607 / 680 |
| 04 (5 rounds) | 175–191 / 91 | 255–273 / 318 | 0.82–0.84 | 3.0–3.5 / 35 | 3.3–3.5 / 0.53 | 0.67–0.71 | 436–443 / 408 |
| 05 | 54 / 83 | 391 / 325 | 0.84 | 0.22 / 20 | 6.0 / 0.53 | 0.16 | 709 / 748 |
| 06 (and 06a) | 124 / 124 | 458 / 424 | 0.87 | 0.67 / 14 | 3.4 / 0.53 | 0.87 | 1094 / 1096 |
| 15 | 103 / 115 | 334 / 297 | 0.87 | 2.7 / 35 | 2.5 / 0.53 | 0.81 | 318 / 316 |
| 16 | 88 / 81 | 249 / 220 | 0.94 | 1.8 / 35 | 3.1 / 0.53 | 0.77 | 427 / 440 |
| 17 (5 shots) | 325–396 / 417 | 221–229 / 193 | 0.89–0.91 | 16–46 / 120 | 1.1–1.8 / 0.53 | 0.80–0.84 | 121–124 / 115 |

- **Shockwave: good match, r = 0.82–0.94.** The model peaks are within about ±20 % of the real ones, except the three shots that `../Cipatat/metadata.md` already flags: shot-03 (model 2.5× too strong; miss distance probably larger), shot-04 (2× too weak; miss distance probably smaller) and shot-05 (1.5× too strong).
- **Durations here are measured on the band-limited, upsampled signal.** They are 40–60 µs longer than those in `../Cipatat/metadata.md`, which used a different estimator. Real and clean are measured the same way, so the fit is not affected.
- **Muzzle blast: the levels come from the recordings, not from the model.** At 88–764 m the generator's blast model (Kinney-Graham) is 3–90 times too strong and too short; it is used far outside its valid range. Its waveform shape (Friedlander) is kept.
- **Muzzle blast: weaker match, r = 0.16–0.87.** The real far-field blasts have a much deeper negative phase than a Friedlander pulse (see the figure), and they carry a ground reflection. shot-03 (machinery noise) and shot-05 (0.22 Pa blast, near the noise) fit worst, because their blast barely stands out of the noise.
- **Timing:** the per-channel SW corrections are mostly below 50 µs. The top microphone (channel 4) arrives 70–160 µs earlier than the scene predicts in shots 02, 04, 06 and 17. This suggests that the bullets passed somewhat above the array, while the scene assumes a trajectory turned horizontally by the miss distance.

## What is left out, on purpose

- **Background noise:** wind, machinery tones, insects (see *Background noise – shot-03* in `../Cipatat/metadata.md`).
- **Reverberation, echoes and the ground-reflected SW:** in the recordings that reflection follows 0.6–2.9 ms after the SW, at 0.34–1.0× its amplitude.
- **Microphone and recorder response:** apart from the anti-aliasing low-pass. The real N-waves have a smaller, more rounded rear shock than the clean ones; part of that is the microphone response.

## Uses

- **Clean reference for the detector:** run the gs_det pipeline on clean and real files to separate errors caused by noise and reflections from errors of the method.
- **Controlled noise:** add noise at known SNR, or real noise cut from the recordings, to the clean files.
- **Denoising tests:** a known clean target for each real recording. Only approximate: the clean waveform is a model, not the true noise-free signal.

## Reproduce

```
cd python
python util_gen_cipatat_clean.py                 # ../sound/Cipatat -> ../sound/Cipatat_clean
python util_gen_cipatat_clean.py --align array   # model time differences instead of the recording's
python util_gen_cipatat_clean.py --no-fit        # pure model levels and durations
```
