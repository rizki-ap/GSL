# Zenodo gunshot dataset – metadata

A curated subset of the "Edge-collected gunshot audio" dataset, kept as a **muzzle-blast-only detection and classification test set** (Nov 2026; it replaces the earlier 5 per-channel recordings).

- **Data:** https://zenodo.org/records/7004819 (the full set, 1 993 WAVs, is kept outside the repo in the Dataset-online `Zenodo` folder)
- **Paper:** Kabealo et al., *Data in Brief* 48 (2023) 109091, https://doi.org/10.1016/j.dib.2023.109091 (copy in this folder)
- **Dataset metadata:** `gunshot-audio-all-metadata.csv` (semicolon-separated, all files of the full set)

## Files

| | |
|---|---|
| Recordings | 40: 10 per firearm, 229 labelled shots |
| Naming | `zen_<firearm>_<recorder>_<distance>ft_<direction>_s<session>_<uuid8>.wav` |
| Format | mono, 44.1 kHz, 32-bit float, 3.1–11.8 s; the UMA-8 array's **mean channel**, copied unchanged |
| Labels | `labels.csv`: firearm, calibre, recorder, distance, direction, session, time stamp, number of shots, **shot times (s)**, weakest labelled shot (dB above the median level), clipped samples, source file |

| Firearm (`<firearm>`) | Calibre | Recordings | Shots |
|---|---|---:|---:|
| Ruger AR-556 (`ar15-223`) | .223 | 10 | 90 |
| Glock 17 (`glock17-9mm`) | 9 mm | 10 | 59 |
| Smith & Wesson (`sw38`) | .38 | 10 | 40 |
| Remington 870 (`rem870-12ga`) | 12 gauge | 10 | 40 |

**Selection:** recordings made with the UMA-8 array (the dataset's main recorder), one per combination of recorder, distance and session where possible, preferring files with many shots. Distances 30–90 ft (9–27 m, the CSV gives no unit; feet fits the 85 × 75 ft range), both sessions (2021-11-15 and 2021-11-22), all six recorder positions. No file is clipped.

## Recording setup

- **Site:** Eastern Florida State College firing range, Cocoa, Florida (outdoor, north half of a range divided by a brick wall, about 26 × 23 m).
- **Recorders:** Raspberry Pi Zero W with a MiniDSP UMA-8 v2 array, placed beside and behind the shooter (the `direction` column: quadrant of the recorder relative to the shooter), not downrange.
- **Echo:** a strong reflection, probably off the brick wall, follows the blast (72 ms at recorder P17).

## Review for the GSL pipeline (Nov 2026)

| Use | Verdict | Why |
|---|---|---|
| Muzzle-blast detection, false alarms, shot counting | **Yes** | Many clean blasts at 9–27 m with labelled times, four firearm types, several shots per file |
| Firearm classification (rifle / pistol / shotgun) | **Yes** | Same site and recorders for all four classes |
| Shockwave detection, SW→MB range, calibre | **No** | No shockwave: the recorders were beside and behind the shooter, outside the bullet's Mach cone; the pistols and shotgun are subsonic or barely supersonic |
| Direction finding (TDOA) | **No** | See below |
| Generator calibration | Weak | Uncalibrated level, 44.1 kHz, close range only |

**Why the multichannel files are not used.** The 7 UMA-8 channels of a recording are stored as separate files, and their relative timing is not physical:

- In most recordings the channels carry fixed offsets of 0.4–4 ms against each other (the same for every shot in the file), more than the 0.25 ms the 85 mm array allows.
- In the few that look aligned (spread 0.04–0.13 ms), including the recording previously kept here (972e5199), the delays fit a plane wave only at an apparent speed of 1 000–2 500 m/s instead of 343 m/s: the source would have to be almost overhead, not 9–27 m away on the ground.
- There is no array orientation and the recorder position is only a quadrant, so a bearing could not be checked anyway.

**Detector results (Nov 2026).** `python/gs_det_events.py` with `mb_max_hf_excess_db = 12`:

- **Blasts are found as blasts.** 226 of the 229 labelled shots come out as MB candidates; only 1 is taken for a shockwave. (The earlier per-channel files had been labelled as shockwaves: heard from 9–27 m, a blast is almost as sharp and high-pitched as a distant shockwave.) Onset error: median −0.2 ms against the labels.
- **Several shots per file** (`detect_shots`, mono, no shockwave): 218 of 229 shots found, plus 160 extra events. These are mostly echoes off the brick wall, 25–200 ms after each shot. Two blasts closer than 90 ms cannot both be shots (project limit: bursts of at most 600 rpm), so the weaker is dropped; that removes most echoes (930 extras before this rule) and 3 labelled shots. A later echo cannot be told from a round with one channel, so the detector sets `echo_risk` on an event that follows a ≥ 3 dB stronger one within 300 ms: 106 of the 160 extras are flagged, none of the true shots.
- **What an array adds:** with a synchronised array, echoes are rejected by direction (`direction_gate_deg`). This cannot be tested here (see above); it was tested on generator bursts instead.

**Labels.** The shot times come from the dataset. A few labelled shots are weak (17–30 dB above the median level, e.g. the second shot of `zen_rem870-12ga_P21_90ft_SE_s1_1fab478a.wav`); check those before using them as positives.
