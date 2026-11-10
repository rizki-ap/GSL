#!/usr/bin/env python3
"""
gs_det_events.py
================
DETECTION library: array-level shockwave (SW) / muzzle blast (MB) event
classifier, used by gs_det_signal_prepare.py when [detection] method = classify
(and by ana_sw_mb_range.py). Detection side only -- imports nothing from the
generation block.

Every event in the recording is classified by its own frequency content and
time/duration. A SW and an MB are each optional: a recording gives SW+MB,
SW only, MB only (subsonic weapon, slowed bullet, shot fired away from the
array) or nothing. The MB does NOT have to follow a SW.

1. Candidate events: peaks of the 1-4 kHz envelope and of the MB-band envelopes
   (default 100-1000 Hz and 300-1000 Hz) rising >= 10 dB above their noise floor.
2. Features of each event:
     crest     impulsiveness of the 1-4 kHz envelope (peak vs its 20 ms surroundings)
     rise      10-90 % rise time of the pressure front (median over channels)
     LF width  width of the 50-800 Hz envelope at half height
     HF-LF     1-4 kHz level above noise minus MB-band level above noise
     prom      MB-band level vs the median 50-250 ms before the event
3. SW = 1-4 kHz >= sw_min_hf_snr_db above noise, crest >= sw_min_crest_db, and a
   near-instant front (rise <= sw_max_rise_ms) or a short low-frequency body
   (LF width <= sw_max_lf_width_ms). Of several SWs (burst, echoes): the earliest
   within 6 dB of the strongest.
4. MB candidates = not SW, MB band >= mb_min_lf_snr_db above noise, and
   HF-LF <= mb_max_hf_excess_db. Score = prom - 0.3 * HF-LF (HF-LF clipped to
   +-15 dB). The MB is the earliest candidate within 3 dB of the best score.
5. Only if a SW was found: its first 25 ms belong to the SW, and its echo tail
   (until the 50 Hz-4 kHz level stays within 10 dB of the noise floor for 30 ms)
   is not searched for the MB -- low-frequency echoes of large calibres look
   exactly like an MB. mb_in_tail = True accepts a clearly low-frequency event
   inside the tail when none is found after it (very close shots).
6. Onsets: first sample where the envelope reaches 30 % of the event's peak
   (SW: 0.1 ms envelope of the SW band; MB: causal-filter envelope, see
   causal_envelope()).

The default thresholds were tuned on real recordings (Cipatat test range, 4 ch,
51.2 kHz: 17 shots from 5.56 mm to 120 mm, 300-1470 m). They are NOT tuned on
the generator's output.
"""

import numpy as np
from scipy.signal import butter, hilbert, sosfilt, sosfiltfilt

DEFAULT_PARAMS = dict(
    sw_band_hz=(1000.0, 4000.0),
    mb_bands_hz=((100.0, 1000.0), (300.0, 1000.0)),
    sw_min_hf_snr_db=25.0,
    sw_min_crest_db=12.0,
    sw_max_rise_ms=0.08,
    sw_max_lf_width_ms=3.0,
    mb_min_lf_snr_db=15.0,
    mb_max_hf_excess_db=10.0,
    mb_in_tail=False,
)


# --------------------------------------------------------------------------- helpers
def bandpass(x, fs, lo, hi, order=4):
    hi = min(hi, 0.45 * fs)
    sos = butter(order, [lo, hi], btype="bandpass", fs=fs, output="sos")
    return sosfiltfilt(sos, x, axis=0)


def envelope(x, fs, smooth_s):
    """Hilbert envelope, averaged over channels, smoothed with a moving average."""
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    e = np.abs(hilbert(x, axis=0)).mean(axis=1)
    w = max(1, int(smooth_s * fs))
    return np.convolve(e, np.ones(w) / w, mode="same")


def causal_envelope(x, fs, band, smooth_s):
    """Envelope of a CAUSAL band-pass (sosfilt) with a trailing moving average. Used for
    onsets: the zero-phase filter rings BEFORE the event (tens of ms for a 100 Hz high-pass),
    which put muzzle-blast onsets ~20 ms early on clean signals; a causal filter cannot
    respond before the event (it is ~1 ms late instead)."""
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    lo, hi = band
    sos = butter(4, [lo, min(hi, 0.45 * fs)], btype="bandpass", fs=fs, output="sos")
    e = np.abs(hilbert(sosfilt(sos, x, axis=0), axis=0)).mean(axis=1)
    w = max(1, int(smooth_s * fs))
    return np.convolve(e, np.ones(w) / w, mode="full")[:len(e)]


def peak_list(env, fs, noise, i0, i1, k=8, excl_s=0.03):
    """Up to k strongest peaks of env[i0:i1] (at least excl_s apart): [(index, dB above noise)]."""
    seg = env[i0:i1].copy()
    out = []
    for _ in range(k):
        if seg.size == 0 or seg.max() <= 0:
            break
        p = int(np.argmax(seg))
        out.append((i0 + p, 20 * np.log10(seg[p] / noise)))
        seg[max(0, p - int(excl_s * fs)):p + int(excl_s * fs)] = 0
    return out


def onset(env, peak_idx, fs, back_s, frac=0.3):
    """First sample (searching back from the peak, at most back_s) at frac of the peak."""
    a = max(0, peak_idx - int(back_s * fs))
    seg = env[a:peak_idx + 1]
    below = np.flatnonzero(seg < frac * env[peak_idx])
    return a + int(below[-1]) + 1 if below.size else a


def onset_sustained(env, peak_idx, fs, back_s, frac=0.3, hold_s=0.0005):
    """First sample in [peak - back_s, peak] from which env stays >= frac * peak for hold_s.
    More robust than onset() on jagged (causal, lightly smoothed) envelopes, where the last
    dip below frac can sit well after the true start of a long event."""
    a = max(0, peak_idx - int(back_s * fs))
    above = env[a:peak_idx + 1] >= frac * env[peak_idx]
    w = max(1, int(hold_s * fs))
    run = np.convolve(above.astype(int), np.ones(w, dtype=int), mode="full")[w - 1:]  # run[i] = sum above[i:i+w]
    ok = np.flatnonzero(run[:len(above)] >= min(w, len(above)))
    return a + int(ok[0]) if ok.size else peak_idx


def _width50(e, pk, fs, max_s):
    thr, lim = 0.5 * e[pk], int(max_s * fs)
    a = b = pk
    while a > 0 and e[a] >= thr and pk - a < lim:
        a -= 1
    while b < len(e) - 1 and e[b] >= thr and b - pk < lim:
        b += 1
    return (b - a) / fs * 1e3


def _rise_ms(xbb, ph, fs):
    """10-90 % rise of the signed broadband pressure at the event, median over channels."""
    out = []
    for ch in range(xbb.shape[1]):
        s = xbb[:, ch]
        w0, w1 = max(0, ph - int(0.001 * fs)), ph + int(0.001 * fs)
        p = w0 + int(np.argmax(np.abs(s[w0:w1])))
        v = s[p]
        seg = np.sign(v) * s[max(0, p - int(0.005 * fs)):p + 1]
        lo = np.flatnonzero(seg < 0.1 * abs(v))
        i10 = int(lo[-1]) if lo.size else 0
        i90 = i10 + int(np.flatnonzero(seg[i10:] >= 0.9 * abs(v))[0])
        out.append((i90 - i10) / fs * 1e3)
    return float(np.median(out))


class EventFeatures:
    """Per-recording envelopes; .at(idx) gives the features of the event near sample idx."""

    def __init__(self, x, fs, params):
        self.fs = fs
        self.xbb = bandpass(x, fs, 20, 15000)
        lo, hi = params["sw_band_hz"]
        self.hf = envelope(bandpass(x, fs, lo, hi), fs, 1e-4)         # sharp HF envelope
        self.lfw = envelope(bandpass(x, fs, 50, 800), fs, 1e-3)        # for LF width
        k5 = np.ones(int(5e-3 * fs)) / int(5e-3 * fs)
        self.hf5 = np.convolve(self.hf, k5, mode="same")
        self.n_hf5 = np.median(self.hf5) + 1e-30
        self.mb_bands = [tuple(b) for b in params["mb_bands_hz"]]
        self.lf5s = [envelope(bandpass(x, fs, lo, hi), fs, 5e-3) for lo, hi in self.mb_bands]
        self.n_lf5s = [np.median(e) + 1e-30 for e in self.lf5s]

    def at(self, idx):
        fs = self.fs
        a, b = max(0, idx - int(0.01 * fs)), min(len(self.hf), idx + int(0.02 * fs))
        ph = a + int(np.argmax(self.hf[a:b]))
        pw = a + int(np.argmax(self.lfw[a:b]))
        win = self.hf[max(0, ph - int(0.005 * fs)):ph + int(0.015 * fs)]
        snr_h = 20 * np.log10(self.hf5[a:b].max() / self.n_hf5)
        best = None
        for k, (env, n) in enumerate(zip(self.lf5s, self.n_lf5s)):
            pl_k = a + int(np.argmax(env[a:b]))
            snr_k = 20 * np.log10(env[pl_k] / n)
            if best is None or snr_k > best[2]:
                best = (k, pl_k, snr_k)
        kb, pl, snr_l = best
        env = self.lf5s[kb]
        bg = np.median(env[max(0, pl - int(0.25 * fs)):max(1, pl - int(0.05 * fs))]) + 1e-30
        return dict(idx_hf=ph, idx_lf=pl, band_k=kb, snr_hf=snr_h, snr_lf=snr_l,
                    crest=20 * np.log10(self.hf[ph] / (win.mean() + 1e-30)),
                    rise_ms=_rise_ms(self.xbb, ph, fs),
                    lf_width_ms=_width50(self.lfw, pw, fs, 0.05),
                    hf_excess=snr_h - snr_l,
                    prominence=20 * np.log10(env[pl] / bg))


def is_sw_like(f, p):
    """Shockwave = strong, impulsive SW-band event with a near-instant front or a short LF body."""
    return (f["snr_hf"] >= p["sw_min_hf_snr_db"] and f["crest"] >= p["sw_min_crest_db"]
            and (f["rise_ms"] <= p["sw_max_rise_ms"] or f["lf_width_ms"] <= p["sw_max_lf_width_ms"]))


# --------------------------------------------------------------------------- main entry
def detect_sw_mb(x, fs, params=None, verbose=False):
    """Array-level SW/MB detection on x (n_samples, n_channels).

    Returns dict: t_sw, t_mb (onset times in s, or None), sw_idx / mb_idx (peak samples),
    n_sw (SWs in the burst train), t_sw_tail_end (end of the SW echo tail), mb_band (Hz),
    sw_features / mb_features, mb_hidden_possible (True if a SW tail was the only place
    left to look), event_class ('SW+MB', 'SW only', 'MB only', 'none'), events (all
    classified events with their features)."""
    p = dict(DEFAULT_PARAMS)
    if params:
        p.update({k: v for k, v in params.items() if v is not None})
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    F = EventFeatures(x, fs, p)

    cand = set()
    for env, noise, excl in [(F.hf5, F.n_hf5, 0.005)] + [(e, n, 0.01) for e, n in zip(F.lf5s, F.n_lf5s)]:
        # zero-phase filters ring at the file ends: skip the first 5 ms and the last 60 ms
        for i, s in peak_list(env, fs, noise, int(0.005 * fs), len(env) - int(0.06 * fs), k=25, excl_s=excl):
            if s >= 10.0:
                cand.add(i)
    events = []
    for i in sorted(cand):
        f = F.at(i)
        if any(abs(f["idx_hf"] - e["idx_hf"]) < int(0.003 * fs) for e in events):
            continue                                   # same event found in several bands
        f["t"] = f["idx_hf"] / fs
        f["cls"] = ("SW" if is_sw_like(f, p) else
                    "MB?" if f["snr_lf"] >= p["mb_min_lf_snr_db"] and f["hf_excess"] <= p["mb_max_hf_excess_db"]
                    else "-")
        events.append(f)

    res = dict(t_sw=None, t_mb=None, sw_idx=None, mb_idx=None, n_sw=0, t_sw_tail_end=None,
               mb_band=None, sw_features=None, mb_features=None, mb_hidden_possible=False,
               events=events, params=p)
    sws = [e for e in events if e["cls"] == "SW"]
    coda = None
    sw = None
    if sws:
        best = max(e["snr_hf"] for e in sws)
        train = sorted((e for e in sws if e["snr_hf"] >= best - 6.0), key=lambda e: e["t"])
        sw = train[0]
        res.update(t_sw=onset(F.hf, sw["idx_hf"], fs, 0.002) / fs, sw_idx=sw["idx_hf"],
                   n_sw=len(train), sw_features=sw)
        env_bb = envelope(bandpass(x, fs, 50, 4000), fs, 0.01)
        quiet = env_bb < np.median(env_bb) * 10 ** (10 / 20)
        cs = np.concatenate([[0], np.cumsum(quiet)])
        w, k = int(0.03 * fs), train[-1]["idx_hf"] + int(0.01 * fs)
        while k + w < len(env_bb) and cs[k + w] - cs[k] < w:
            k += int(0.002 * fs)
        res["t_sw_tail_end"] = k / fs
        coda = (sw["idx_hf"], k)
        for e in events:                               # the SW's own low-frequency body
            if e["cls"] == "MB?" and any(-0.005 <= e["t"] - s["t"] <= 0.025 for s in sws):
                e["cls"] = "SW-tail"

    mbs = [e for e in events if e["cls"] == "MB?"]
    for e in mbs:
        e["score"] = e["prominence"] - 0.3 * float(np.clip(e["hf_excess"], -15, 15))
    pool = mbs
    if coda:
        outside = [e for e in mbs if not (coda[0] <= e["idx_lf"] < coda[1])]
        inside = ([e for e in mbs if e["hf_excess"] <= 0.0 and e["score"] >= 15.0]
                  if p["mb_in_tail"] else [])
        pool = outside or inside
        res["mb_hidden_possible"] = not outside
    if pool:
        top = max(e["score"] for e in pool)
        mb = min((e for e in pool if e["score"] >= top - 3.0), key=lambda e: e["t"])
        mb["cls"] = "MB"
        back = 0.05
        if sw is not None and mb["idx_lf"] > sw["idx_hf"]:   # keep the onset out of the SW
            back = max(0.002, min(back, (mb["idx_lf"] - sw["idx_hf"]) / fs - 0.025))
        band = F.mb_bands[mb["band_k"]]
        ec = causal_envelope(x, fs, band, 1e-3)
        a, b = max(0, mb["idx_lf"] - int(0.01 * fs)), min(len(ec), mb["idx_lf"] + int(0.02 * fs))
        pk = a + int(np.argmax(ec[a:b]))
        res.update(t_mb=onset_sustained(ec, pk, fs, back + (pk - mb["idx_lf"]) / fs) / fs, mb_idx=pk,
                   mb_band=F.mb_bands[mb["band_k"]], mb_features=mb)

    res["event_class"] = ("SW+MB" if res["t_sw"] is not None and res["t_mb"] is not None else
                          "SW only" if res["t_sw"] is not None else
                          "MB only" if res["t_mb"] is not None else "none")
    if verbose:
        for e in events:
            if e["cls"] in ("-", "SW-tail"):
                continue
            print(f"    {e['t']:8.4f} s  {e['cls']:3s}  HF {e['snr_hf']:5.1f} dB  LF {e['snr_lf']:5.1f} dB  "
                  f"crest {e['crest']:5.1f}  rise {e['rise_ms']:6.3f} ms  LFwidth {e['lf_width_ms']:5.1f} ms  "
                  f"HF-LF {e['hf_excess']:+6.1f}  prom {e['prominence']:5.1f}"
                  + (f"  score {e['score']:5.1f}" if "score" in e else ""))
    return res


def channel_onset(sig, fs, center_idx, half_win_s, smooth_s, band, back_s=None, frac=0.3, causal=False):
    """Onset of one channel near an array-level event: band-pass, envelope, local peak within
    +-half_win_s of center_idx (the array-level peak), then searching back at most back_s
    (default half_win_s) for the first sample at frac of that peak. Returns sample index."""
    back_s = half_win_s if back_s is None else back_s
    a = max(0, center_idx - int(half_win_s * fs))
    b = min(len(sig), center_idx + int(half_win_s * fs))
    pad = int((back_s + 0.05) * fs)
    lo_i, hi_i = max(0, a - pad), min(len(sig), b + int(0.05 * fs))
    e = (causal_envelope(sig[lo_i:hi_i], fs, band, smooth_s) if causal else
         envelope(bandpass(sig[lo_i:hi_i], fs, *band), fs, smooth_s))
    pk = (a - lo_i) + int(np.argmax(e[a - lo_i:b - lo_i]))
    if causal:
        return lo_i + onset_sustained(e, pk, fs, back_s, frac)
    return lo_i + onset(e, pk, fs, back_s, frac)
