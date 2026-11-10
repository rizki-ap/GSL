# Zenodo gunshot dataset – metadata

A subset of the "Edge-collected gunshot audio" dataset:

- **Data:** https://zenodo.org/records/7004819
- **Paper:** Kabealo et al., *Data in Brief* 48 (2023) 109091, https://doi.org/10.1016/j.dib.2023.109091. A copy of the paper is in this folder.
- **Metadata:** `gunshot-audio-all-metadata.csv` (semicolon-separated) holds the dataset's own metadata for all files.

## Recording setup

| Property | Value |
|---|---|
| Location | Eastern Florida State College firing range, Cocoa, Florida, USA (outdoor) |
| Range size | North half of a range divided by a brick wall: about 85 ft wide × 75 ft deep (≈ 26 × 23 m) |
| Shooter | Certified firearms instructor, standing, firing at a target downrange |
| Array | MiniDSP UMA-8 USB mic array V2, on a Raspberry Pi Zero W |
| Files | One mono WAV per microphone, `<name>_ch0.wav` … `_ch6.wav` (7 mics), plus `<name>_mean.wav` |
| Format | 44.1 kHz, 32-bit float, 3.75 s per file |
| Gunshot time | Labelled at 1.75 s in every file |

The recorders were placed around the shooter (beside and behind), not downrange along the line of fire.

## Recordings in this folder

| Name | Source file (uuid_version) | Weapon | Caliber | Device | Distance* | Direction | Shots in file |
|---|---|---|---|---|---:|---|---:|
| glock_a283 | a283ea70-…-ba4e8507455a_v0 | Glock 17 | 9 mm | P17 | 50 | NORTH-EAST | 1 |
| SmithWesson38cal_47V0 | 47ec5352-…-b13901c7eb4d_v0 | Smith & Wesson | .38 | P19 | 30 | SOUTH-WEST | 1 |
| remington12g_0aV1 | 0a07b229-…-c94cbc7b1487_v1 | Remington 870 | 12 gauge | P16 | 40 | NORTH-EAST | 1 |
| ruger_04_1s | 0495b86d-…-0cfbc160cf3d_v1 | Ruger AR-556 | .223 | P17 | 50 | NORTH-EAST | 1 |
| ruger_ar556_223cal_97V0 | 972e5199-…-ec7c2903f2ee_v0 | Ruger AR-556 | .223 | P21 | 90 | SOUTH-EAST | 1 |

\* The CSV gives no unit. It is most likely **feet**: the range is only about 85 × 75 ft.

## Analysis (Nov 2026)

### Shockwave and muzzle blast

None of the five recordings has a clearly separate shockwave and muzzle blast. Each has one main bang at 1.749 s, which is the muzzle blast.

| Name | Bullet speed vs sound | Events seen | Shockwave? |
|---|---|---|---|
| glock_a283 | ~350–370 m/s, at or just above the speed of sound | Bang at 1.749 s, echo at 1.822 s | No |
| SmithWesson38cal_47V0 | ~270 m/s, slower than sound | One long bang, 1.745–1.80 s | No (none expected) |
| remington12g_0aV1 | Pellets ~400 m/s, only slightly faster than sound | Bang at 1.75 s, small later events | No clear one |
| ruger_04_1s | ~900 m/s, well above the speed of sound | Bang at 1.749 s, echo at 1.822 s | No |
| ruger_ar556_223cal_97V0 | ~900 m/s | One bang at 1.749 s | Possibly, merged with the blast |

- **Why the rifle files have no separate crack:** a bullet's shockwave only reaches a microphone ahead of the muzzle, inside the bullet's cone of sound. These recorders were beside and behind the shooter.
- **The 1.822 s event is an echo:** it sits 72 ms after the bang, the same in the Glock and Ruger files, both from recorder P17. The Glock is barely supersonic, so this is a reflection, probably off the brick wall that splits the range.
- **ruger_ar556_223cal_97V0:** its first pressure drop is much sharper than in the other rifle file (about 0.07 ms against 0.4 ms). That fits a shockwave arriving almost together with the blast, but at 44.1 kHz the two cannot be separated.
- **Polarity:** all blasts start with a negative swing, so the recordings appear to have inverted polarity.

### Data quality

| Name | Spread in arrival time across the 7 channels | Comment |
|---|---:|---|
| ruger_ar556_223cal_97V0 | 0.07 ms | Good. The UMA-8's maximum possible is 0.25 ms (85 mm aperture) |
| glock_a283 | 0.52 ms | More than the array allows |
| ruger_04_1s | 0.73 ms | More than the array allows |
| SmithWesson38cal_47V0 | 3.1 ms | Channels not time-aligned with each other; unusable for direction |
| remington12g_0aV1 | 4.4 ms | Same, and 18 samples are clipped (\|x\| ≥ 1.0) |

Only **ruger_ar556_223cal_97V0** is synchronised well enough for direction finding.

### Detector

The gs_det classifier (`python/gs_det_events.py`) labels the bang in the Glock, S&W and Ruger files as a **shockwave**. Heard from 10–30 m, a muzzle blast is as sharp and high-pitched as a crack heard from far away. The classifier was tuned on Cipatat (300–1470 m), so close-range recordings are outside what it handles.

### Use

- **Good for:** muzzle-blast detection, and bearing with ruger_ar556_223cal_97V0.
- **Not suitable for:** shockwave work, or range from the shockwave-to-blast delay.
