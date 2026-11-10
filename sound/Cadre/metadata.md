# Cadre gunshot dataset – metadata

A sample of the Gunshot Audio Forensics Dataset:

- **Source:** Cadre Research, released June 2018, https://cadreforensics.com/audio/. See `AudioDatafileDescription.pdf` in this folder.
- **Funding:** NIJ Grant 2016-DN-BX-0183, "Development of Computational Methods for the Audio Analysis of Gunshots".
- **Full dataset:** recorded in summer 2017 in rural Arizona. It covers 20 firearms, 20 recording positions and 4 recording devices, about 10,000 shots in total.

## Naming

`XX_YYYY_SZZ.wav`. These files also start with a weapon prefix (`M16_`).

| Part | Meaning |
|---|---|
| XX | Recording device: ZM = Zoom H4N, IP = iPhone 7, SA = Samsung Galaxy S7 |
| YYYY | Experiment ID: geometry number + letter (A/B) |
| SZZ | Shot number |

Angles are measured from the line of fire, with 0° pointing towards the target.

## Files in this folder

| Files | Weapon | Device | Experiment | Geometry | Shots |
|---|---|---|---|---|---:|
| M16_IP_017A_S01–S06 | M16 (5.56×45 mm) | iPhone 7 | 017A | −30°, 20 m | 6 |
| M16_IP_017B_S01–S06 | M16 (5.56×45 mm) | iPhone 7 | 017B | −30°, 40 m | 6 |

- **Format:** mono, 44.1 kHz, 16-bit PCM, about 2 s per file.
- **No direction finding:** each file is a single microphone.
- **Clipping:** 3 files have one sample at full scale (017A_S01, 017A_S05, 017B_S02).

## Analysis (Nov 2026)

### Shockwave and muzzle blast

At −30° the recorder is ahead of the muzzle, inside the bullet's cone of sound. **Every shot has both a shockwave and a muzzle blast.**

- **Shockwave:** the first sharp spike, about 0.3 ms long, with the same shape at both distances.
- **Muzzle blast:** at 20 m it shows best in 500–2000 Hz, and at 40 m in 2–8 kHz.

| Geometry | Measured SW→MB delay (6 shots) | Expected delay* | Range from delay | Error |
|---|---|---|---:|---:|
| 017A: −30°, 20 m | 10.5 ± 0.3 ms | 12.0–12.6 ms | 17.6 m | −12 % |
| 017B: −30°, 40 m | 20.3 ± 0.6 ms | 23.7–25.1 ms | 34.1 m | −15 % |

\*From the drag model in `python/ana_sw_mb_range.py`, for M855 at v0 = 940 m/s, λ = 950 m, with miss distance = D·sin 30° and c = 343–352 m/s (Arizona summer). The range column uses c = 352 m/s.

- **The pairing is right:** the delay doubles when the distance doubles (ratio 1.93), as expected for a fixed angle.
- **Systematic short range:** the delays are 13–16 % shorter than predicted. Possible causes:
  - **Angle:** an angle of about 33° instead of 30° would explain it on its own.
  - **Bullet speed:** a bullet of about 900 m/s instead of 940 m/s.
  - **Timing:** the muzzle-blast onset definition could be slightly off.

### Recording artifact

Every file has a small, sharp event about 23 ms before the shockwave, with the level rising between it and the shockwave.

- It sits at −22.95 ms in 9 of the 12 files and −21.6 to −21.8 ms in the other 3, at both distances. That is about 1024 samples, the frame size of the AAC compression iPhones use.
- This points to compression **pre-echo**, not an acoustic event.
- If so, the files are decoded AAC, so the timing and shape of very sharp sounds are not fully reliable.

### Detector

The gs_det classifier (`python/gs_det_events.py`) fails on these files: in 12 shots it found a shockwave in 7 and a plausible muzzle blast in none.

- **Muzzle blast hidden by design:** at 20–40 m the muzzle blast arrives inside the shockwave's echo period, which the classifier does not search by default. `mb_in_tail = true` is meant for this case but still needs testing.
- **iPhone processing:** compression and automatic gain control change the sound, so it no longer looks like the real-world recordings the classifier was tuned on.

### Next step

The full dataset has 20 positions per weapon (0.5 m to 150 m, 0° to 180°) and also includes the Zoom H4N, a recorder that does not compress the audio. Recordings from it would allow proper tests of range estimation and of the classifier across many distances and angles.
