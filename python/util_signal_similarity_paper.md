Created the paper-focused evaluation version:

[Download `util_signal_similarity_paper.py`](sandbox:/mnt/data/util_signal_similarity_paper.py)

It is **584 lines** and has been syntax-checked successfully.

### What this version implements

**§4.3.1 — Waveform morphology**

* Pearson \(r\)
* `SNR_synth`
* peak-amplitude error \(\epsilon_{\Delta p}\)
* duration error \(\epsilon_T\)
* best-fit lag before morphology comparison
* no independent peak normalization

**§4.3.2 — Frequency-domain similarity**

* Tukey-window PSD
* 4096-point zero-padded FFT
* 50 Hz–8 kHz evaluation band
* averaged PSD
* MD-defined LSD
* spectral peak-frequency error
* N-wave first non-trivial null-frequency error
* pooled ensemble magnitude-squared coherence:

  $$
  \bar C(f)=
  \frac{|\sum_kS_{xy,k}|^2}
       {(\sum_kS_{xx,k})(\sum_kS_{yy,k})}
  $$
* \(C>0.8\) bandwidth

### Example

```bash
python util_signal_similarity_paper.py \
    measured.wav \
    synthetic.wav \
    --event-type nwave \
    --event-start-ms 2 \
    --event-duration-ms 10 \
    --json result.json
```

For muzzle blast:

```bash
python util_signal_similarity_paper.py \
    measured.wav \
    synthetic.wav \
    --event-type blast \
    --event-start-ms 5 \
    --event-duration-ms 20
```

You can also specify different event windows for measured and synthetic signals:

```bash
--measured-start-ms ...
--synth-start-ms ...
--measured-duration-ms ...
--synth-duration-ms ...
```

### One important qualification

The MD itself still leaves several experimental parameters as `[ ]` or unspecified—particularly the **event windows**, Tukey parameter, coherence segment length, and spectral smoothing method. I therefore exposed these as explicit parameters rather than silently inventing them.

Also, the `paper_evaluate_many()` API is included specifically so the final experiment can pool **multiple sensors and recordings**, which is what the MD's ensemble definitions require.

The next step I would recommend is to **integrate this directly into `GSL/python/` and modify `project.md` so every reported Table 2/3 number can be reproduced with one command**.
