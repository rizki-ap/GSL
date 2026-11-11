# util_signal_similarity.py – measured vs synthesized gunshots

Implements the evaluation of `publication/Paper-Generator.md`, Section 4.3, for one or more pairs of recordings (measured, synthesized).

## Metrics

| Paper section | Metric | How it is computed |
|---|---|---|
| 4.3.1 | Pearson r, SNR_synth | on each channel's event window after removing the best-fit lag |
| 4.3.1 | Peak error ε_Δp | peak \|x\| of the aligned windows; amplitudes are never normalized |
| 4.3.1 | Duration error ε_T | T measured from each signal: N-wave T_N (front shock to rear shock), blast T_mb (positive phase, Friedlander T) |
| 4.3.2 | LSD | Tukey-window event PSD (4096-point zero padding), averaged over channels and recordings |
| 4.3.2 | ε_f_pk (blast) | peak of the smoothed averaged PSD; the Friedlander value 1/(2π T_mb) is also reported |
| 4.3.2 | ε_f_null (N-wave) | first null searched at 0.6–1.6 × 1.43/T_N, with T_N measured from the same signal |
| 4.3.2 | Coherence | ensemble MSC pooled over all segments of all windows (channels × recordings); reported **only if K ≥ 8** segments, because with few segments it is biased towards 1 (with K = 1 it is exactly 1 even for unrelated noise) |
| 4.3.3 | TDOA MAE / RMSE / Bias | all pairs i < j, plus the N−1 independent pairs (reference channel 1); plain cross-correlation (default) or GCC-PHAT, sinc-upsampled ×16 |
| 4.3.4 | Shock-to-blast delay error | per channel, (t_mb − t_sw) synthesized minus measured |

**TDOA reference:** with `--synth-json` (the generator's JSON), measured TDOAs are compared with the model labels, as the paper defines. Without it, the synthesized TDOAs are compared with the measured ones.

## Defaults

| Setting | N-wave | Blast |
|---|---|---|
| Window | 1.5 ms, from 0.25 ms before the onset | 20 ms, from 2 ms before the onset |
| Band | 200 Hz – min(20 kHz, 0.45 fs) | 50 Hz – 4 kHz (paper Table 2) |

- **N-wave band:** wider than the paper's 8 kHz, because the N-wave null 1.43/T_N lies near 5–7 kHz for T_N = 0.2–0.3 ms.
- **High-pass:** both signals are high-passed at the band's lower edge before windowing (wind, recoil and drift are not modelled); `--no-highpass` turns this off.
- **Event times:** found with `gs_det_events.detect_sw_mb` (measured) or read from the generator JSON (synthesized). They can also be given by hand (`--meas-sw-ms` …). The detector does not find the synthetic events reliably, so use `--synth-json` for generator output.
- **Polarity:** a mismatch is counted when the first major excursions differ in sign or the best alignment is inverted.

## Examples

```
# one Cipatat shot against a generator run (scene rebuilt from the reference values)
python util_signal_similarity.py ..\sound\Cipatat\shot-01_audio_4ch_cut.wav work\sim01_clean.wav ^
       --synth-json work\sim01.json --meas-scale 1520 --plot work\sim_plots --json work\sim01_eval.json

# pool several recordings for the paper tables
python util_signal_similarity.py --pair m1.wav s1.wav --pair m2.wav s2.wav ^
       --synth-json s1.json --synth-json s2.json --meas-scale 1520 --json table.json
```

- **`--meas-scale 1520`** converts PILAR-V full scale to Pa (the estimate in `sound/Cipatat/metadata.md`); the generator's `*_clean.wav` is already in Pa. Without a common scale the peak error and SNR_synth are meaningless.
- **Library use:** `from util_signal_similarity import load_wav, evaluate, read_generator_json`.

## First result (Nov 2026): Cipatat shot-01 vs generator

Scene: M16 5.56 mm, 306 m, azimuth 269°, miss distance 20 m, 26 °C, PILAR-V, library muzzle velocity and drag.

| Event | Pearson r | Peak error | T meas / synth | Spectrum | TDOA, measured vs model |
|---|---:|---:|---|---|---:|
| N-wave | 0.82 | 28 % | 0.234 / 0.205 ms | f_null 5.9 / 7.0 kHz, coherence > 0.8 over 4.5 kHz | MAE 34 µs |
| Blast | −0.14 | 920 % | 1.6 / 0.66 ms | f_pk 88 / 275 Hz | MAE 7 µs |

- **Shockwave:** reproduced well in shape and timing.
- **Muzzle blast at 300 m:** far too strong (about 20 Pa against about 2 Pa) and too short in the generator. Its Kinney-Graham model is extrapolated far beyond its valid range at this distance; the generator warns about this.
- **Bullet side:** the shockwave TDOAs also showed which side of the array the bullet passed: 34 µs MAE on one side against 810 µs on the other.
