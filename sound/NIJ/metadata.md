# NIJ gunshot dataset (2014-DN-BX-K034) – metadata

Rifle shots recorded by Routh & Maher (Montana State University) at 500 kHz with 12 microphones on a semicircle around the gun. Paper: *Recording anechoic gunshot waveforms of several firearms at 500 kilohertz sampling rate*, POMA 26, 030001 (2016); NIJ grant 2014-DN-BX-K034.

These files are 70 ms cuts around each shot, converted from the original `.mat` files. The originals (4–15 s per shot, plus the .22, .45 and 12-gauge shots) and the converter `nij_mat_to_wav.py` are kept outside the repo, in the Dataset-online `NIJ_Grant` folder.

## Files

| Files | Firearm | Ammunition | Shots | Shockwave |
|---|---|---|---:|---|
| `record_MB_AR15r_shot01–10.wav` | Stag Arms AR-15 | 5.56 NATO XM855, 62 gr (~945 m/s) | 10 | yes, mics 1–3 |
| `record_MB_308r_shot01–10.wav` | Surgeon / AI rifle | .308 Winchester (7.62×51) | 10 | yes, mics 1–3 |

Each WAV has a `.json` with the firearm, mic azimuths and sensitivities, polarity and the cut position in the original file.

## Format

| Property | Value |
|---|---|
| Channels | 12, one per microphone, in mic order 1–12 |
| Sample rate | 500 000 Hz |
| Sample format | 32-bit float |
| Length | 70 ms: 10 ms before the first arrival, 60 ms after |
| Values | ADC counts / 32768 (±1.0 = ±5 V at the ADC); peaks reach 0.12–0.16 |
| Polarity | corrected: the stored data were inverted and have been multiplied by −1, so compression is positive |

## Geometry

- **Microphones:** 12 × G.R.A.S. 40DP (1/8 in, 0.68–0.91 mV/Pa, flat 6.5 Hz–70 kHz, up to 175 dB), individually calibrated (sensitivities in the `.json`).
- **Layout:** semicircle of radius 3 m centred on the shooting position. Mic 1 is on the line of fire (0°), then every 16.4° to mic 12 at 180°, directly behind the shooter.
- **Height:** shooter and microphones 3 m above the ground, so the ground reflection arrives about 10 ms after the direct sound (inside the cut).
- **Not an array for localization:** the microphones surround the gun. The gun was not exactly at the centre: the muzzle-blast arrival times spread by about 1–1.5 ms across the mics, in the same pattern for every shot.

## What is in the signals

- **Shockwave on mics 1–3 (0°, 16°, 33°) for both rifles; none from mic 4 (49°) on.**
  - The N-wave comes before the blast: AR-15 by 2.5, 1.8 and 0.6 ms on mics 1–3.
  - Duration about 96–105 µs (AR-15) and 128–135 µs (.308), sampled with about 50 samples per N-wave.
  - Whitham's duration formula puts mic 1 about 1 m from the trajectory.
  - Keep 500 kHz for shockwave work: below about 200 kHz the N-wave shrinks to a few samples and loses its peak.
- **Muzzle blast on all 12 mics, strongly directional:** the peak drops by about 29 dB (AR-15) from the front (mic 1) to the back (mic 12).
- **Level, approximate:** the paper gives 7 000–8 700 Pa on-axis for the AR-15 at 3 m, about 1.7 Pa per ADC count on mic 1 (≈ 172 dB). Inline attenuator settings are not documented, so treat this as approximate.
- **Mic 12 (behind the shooter)** is weak and can carry isolated noise spikes.

## Use in GSL

- **Generator:** measured muzzle-blast waveform, duration and directivity per firearm (the generator's blast is omnidirectional Kinney-Graham); near-field N-waves (b ≈ 1 m), sharply resolved.
- **Classifier:** clean rifle muzzle-blast and shockwave templates.
- **Not for:** localization tests (the mics surround the gun) or long-range behaviour (3 m only, no propagation).

## Reproduce the cuts

From the folder with the original `.mat` files:

```
python nij_mat_to_wav.py "MB_AR15r_shot*.mat" "MB_308r_shot*.mat" --cut
```
