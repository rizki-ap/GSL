# Cipatat gunshot dataset – metadata

Reference (ground-truth) values for the Cipatat test-range recordings, converted from `Resultats.xlsx` (theoretical values of the reference shots).

## Recording format

All files in this folder share the same format:

| Property | Value |
|---|---|
| Channels | 4 |
| Sample rate | 51 200 Hz |
| Sample format | 24-bit PCM (`PCM_24`) |
| Naming | `shot-NN_audio_4ch_cut.wav`, cut from the original `Shot_NN.tkb` recordings |

## Columns

| Column | Original (FR) | Meaning |
|---|---|---|
| Shot | Nom du fichier | Source recording (`.tkb`) |
| WAV file | – | Matching file in this folder |
| Weapon | Arme | Weapon fired |
| Caliber | Calibre (mm) | Ammunition caliber / designation |
| Type | Shot type | Single shot or burst |
| Rounds | Nb Cp (nombre de coups) | Number of rounds fired |
| Azimuth (°) | Azm | Azimuth of the shooter |
| Elevation (°) | Sit (site) | Elevation (site angle) of the shooter |
| Distance (m) | Dst | Shooter-to-array distance |
| DCPA (m) | DCPA | Distance at Closest Point of Approach: how close the bullet trajectory passes to the array (miss distance) |

## Reference values

| Shot | WAV file | Weapon | Caliber | Type | Rounds | Azimuth (°) | Elevation (°) | Distance (m) | DCPA (m) | Duration (s) |
|---|---|---|---|---|---:|---:|---:|---:|---:|---:|
| Shot_01 | shot-01_audio_4ch_cut.wav | M16 | 5.56×45 mm | Single | 1 | 269 | 0 | 306 | 20 | 2.003 |
| Shot_02 | shot-02_audio_4ch_cut.wav | AK47 | 7.62×39 mm | Single | 1 | 85 | −1 | 306 | 20 | 5.502 |
| Shot_03 | shot-03_audio_4ch_cut.wav | Dragunov | 7.62×54R mm | Single | 1 | 355 | 6 | 500 | 10 | 15.577 |
| Shot_04 | shot-04_audio_4ch_cut.wav | MINIMI | 7.62×51 mm | Burst | 5 | 1 | 1 | 304 | 30 | 2.957 |
| Shot_05 | shot-05_audio_4ch_cut.wav | AX338 / Accuracy | 8.6×70 mm | Single | 1 | 0 | 7 | 517 | 35 | 5.506 |
| Shot_06 | shot-06_audio_4ch_cut.wav | STEYR | 12.7×99 mm | Single | 1 | 11 | 2 | 764 | 30 | 5.504 |
| Shot_07 | shot-07_audio_4ch_cut.wav | Cannon | 14.5×114 mm | Single | 1 | 11 | 2 | 752 | 15 | 5.507 |
| Shot_08 | shot-08_audio_4ch_cut.wav | Cannon | 20 mm | Single | 1 | 0 | 0 | 1302 | 52 | 5.505 |
| Shot_09 | shot-09_audio_4ch_cut.wav | Cannon | 25×137 mm | Single | 1 | 356 | 1 | 1032 | 50 | 5.505 |
| Shot_10 | shot-10_audio_4ch_cut.wav | TP-T-GDOTS | 30×173 mm | Single | 1 | 1 | 0 | 1204 | 73 | 4.016 |
| Shot_11 | shot-11_audio_4ch_cut.wav | Cannon | T40 (fin-stabilized dart and target-practice rounds) | Single | 1 | 1 | 0 | 1470 | 16 | 5.503 |
| Shot_12 | shot-12_audio_4ch_cut.wav | Cannon | 105 mm HEP shell | Single | 1 | 355 | 1 | 1356 | 67 | 5.511 |
| Shot_13 | shot-13_audio_4ch_cut.wav | RPG-7 | 73 mm PG-7V | Single | 1 | 356 | 0 | 350 | – | 9.486 |
| Shot_14 | shot-14_audio_4ch_cut.wav | Cannon | 120 mm BSCC | Single | 1 | 2 | 1 | 1270 | 123 | 6.663 |
| Shot_15 | shot-15_audio_4ch_cut.wav | AK47 | 7.62×39 mm | Single | 1 | 0 | −10 | 300 | 20 | 5.502 |
| Shot_16 | shot-16_audio_4ch_cut.wav | M16 | 5.56×45 mm | Single | 1 | 175 | 0 | 304 | 24 | 5.503 |
| Multiple_Shot | shot-17_audio_4ch_cut.wav | AK47 | 7.62×39 mm | Single | 5 | 28 | 0 | 88 | 4 | 3.375 |

## Notes and open points

- **Shot_13 (RPG-7):** no DCPA given in the source (`----`).
- **Multiple_Shot → shot-17:** the mapping is assumed (it is the only remaining file). The source lists the type as *Single* with 5 rounds, i.e. five single shots.
- **shot-06a_audio_4ch_cut.wav:** has no entry in the source table, and it is byte-identical to `shot-06_audio_4ch_cut.wav` (a duplicate).
- **Long files:** shot-03 (15.6 s), shot-13 (9.5 s) and shot-14 (6.7 s) contain peaks spread over several seconds. Check whether they hold more than one event before running the detector on them.
- **Unknowns to confirm:**
  - The reference direction and sign convention of *Azimuth* (e.g. clockwise from north or from the array's axis).
  - The array geometry and channel order used for these recordings.
- **Source:** `Resultats.xlsx`. Numbers with French formatting were normalized: decimal comma → point (5,56 → 5.56), and thousands separator removed (1 032 → 1032).
