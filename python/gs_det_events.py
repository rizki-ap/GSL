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

1. Candidate events: peaks of the 1-4 kHz envelope and of the MB-band envelope(s)
   (default one band, 150-1000 Hz: best single band on the Cipatat development set,
   shots of 12.7 mm and below) rising >= 10 dB
   above their noise floor.
2. Features of each event:
     crest     impulsiveness of the 1-4 kHz envelope (peak vs its 20 ms surroundings)
     rise      10-90 % rise time of the pressure front (median over channels)
     LF width  width of the 50-800 Hz envelope at half height
     HF-LF     1-4 kHz level above noise minus MB-band level above noise
     prom      MB-band level vs the median 50-250 ms before the event
3. SW = 1-4 kHz >= sw_min_hf_snr_db above noise, crest >= sw_min_crest_db, and a
   near-instant front (rise <= sw_max_rise_ms). Optional alternative to the front
   test: a short low-frequency body (LF width <= sw_max_lf_width_ms; 0 = off, the
   default -- it was only needed for 105/120 mm cannons, outside the 12.7 mm rule). Of several SWs (burst, echoes): the earliest
   within 6 dB of the strongest.
4. MB candidates = not SW, MB band >= mb_min_lf_snr_db above noise,
   HF-LF <= mb_max_hf_excess_db, and prom >= mb_min_prominence_db. Score = prom - 0.3 * HF-LF (HF-LF clipped to
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
51.2 kHz). Development rule (Nov 2026): only Cipatat shots of 12.7 mm and below
are used (5.56 mm to 12.7 mm, 300-764 m); see sound/Cipatat/metadata.md. They are
NOT tuned on the generator's output.
"""

import numpy as np
from scipy.signal import butter, hilbert, sosfilt, sosfiltfilt

DEFAULT_PARAMS = dict(
    sw_band_hz=(1000.0, 4000.0),
    mb_bands_hz=((150.0, 1000.0),),
    sw_min_hf_snr_db=25.0,
    sw_min_crest_db=14.0,
    sw_max_rise_ms=0.10,
    sw_max_lf_width_ms=0.0,         # 0 = rule off (it was only needed for 105/120 mm cannons)
    mb_min_lf_snr_db=11.0,
    mb_max_hf_excess_db=12.0,
    mb_min_prominence_db=6.0,
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


def onset_sustained(env, peak_idx, fs, back_s, frac=0.3, hold_s=0.0005, floor=0.0):
    """First sample in [peak - back_s, peak] from which env stays above
    floor + frac * (peak - floor) for hold_s. floor = the envelope's noise level, so a weak
    event (peak only a few times the noise) is not started early inside the noise.
    More robust than onset() on jagged (causal, lightly smoothed) envelopes, where the last
    dip below the threshold can sit well after the true start of a long event."""
    a = max(0, peak_idx - int(back_s * fs))
    thr = floor + frac * max(env[peak_idx] - floor, 0.0)
    above = env[a:peak_idx + 1] >= thr
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

    def at(self, idx, lf_origin=False):
        """Features of the event near sample idx. lf_origin: idx is a peak of a muzzle-blast-band
        envelope; the event is then anchored on that low-frequency peak (its HF peak is searched
        only from 8 ms before to 3 ms after it), so that a strong HF event nearby (a shockwave or
        its echo in a burst) does not take its place."""
        fs = self.fs
        if lf_origin:
            k = int(np.argmax([e[idx] / n for e, n in zip(self.lf5s, self.n_lf5s)]))
            env = self.lf5s[k]
            a0, b0 = max(0, idx - int(0.01 * fs)), min(len(env), idx + int(0.01 * fs))
            pl0 = a0 + int(np.argmax(env[a0:b0]))
            a, b = max(0, pl0 - int(0.008 * fs)), min(len(self.hf), pl0 + int(0.003 * fs))
        else:
            a, b = max(0, idx - int(0.01 * fs)), min(len(self.hf), idx + int(0.02 * fs))
        ph = a + int(np.argmax(self.hf[a:b]))
        pw = a + int(np.argmax(self.lfw[a:b]))
        win = self.hf[max(0, ph - int(0.005 * fs)):ph + int(0.015 * fs)]
        snr_h = 20 * np.log10(self.hf5[a:b].max() / self.n_hf5)
        best = None
        la, lb = (max(0, ph - int(0.003 * fs)), min(len(self.hf), ph + int(0.012 * fs))) if lf_origin else (a, b)
        for k, (env, n) in enumerate(zip(self.lf5s, self.n_lf5s)):
            pl_k = la + int(np.argmax(env[la:lb]))
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
    width_ok = p["sw_max_lf_width_ms"] > 0 and f["lf_width_ms"] <= p["sw_max_lf_width_ms"]
    return (f["snr_hf"] >= p["sw_min_hf_snr_db"] and f["crest"] >= p["sw_min_crest_db"]
            and (f["rise_ms"] <= p["sw_max_rise_ms"] or width_ok))


def classify_events(x, fs, p, k=25):
    """Candidate events of x (n_samples, n_channels) with their features and a first class:
    'SW', 'MB?' (muzzle-blast candidate) or '-'. Returns (EventFeatures, events sorted by time).
    k = maximum number of candidate peaks taken from each envelope."""
    F = EventFeatures(x, fs, p)
    cand = {}
    for kind, (env, noise, excl) in [("hf", (F.hf5, F.n_hf5, 0.005))] + [("lf", (e, n, 0.01)) for e, n in zip(F.lf5s, F.n_lf5s)]:
        # zero-phase filters ring at the file ends: skip the first 5 ms and the last 60 ms
        for i, s_db in peak_list(env, fs, noise, int(0.005 * fs), len(env) - int(0.06 * fs), k=k, excl_s=excl):
            if s_db >= 10.0:
                cand.setdefault(i, kind)
    events = []
    for i in sorted(cand):
        lf_origin = cand[i] == "lf"
        f = F.at(i, lf_origin)
        for _ in range(3 if not lf_origin else 0):     # re-centre the feature window on the peak
            if f["idx_hf"] == i:                       # (a window cut just before the peak would
                break                                  # report a point on its rising edge)
            i = f["idx_hf"]
            f = F.at(i)
        dup = [e for e in events if abs(f["idx_hf"] - e["idx_hf"]) < int(0.003 * fs)]
        if dup:                                        # same event found in several bands
            continue
        f["t"] = f["idx_hf"] / fs
        f["cls"] = ("SW" if is_sw_like(f, p) else
                    "MB?" if (f["snr_lf"] >= p["mb_min_lf_snr_db"] and f["hf_excess"] <= p["mb_max_hf_excess_db"]
                              and f["prominence"] >= p["mb_min_prominence_db"])
                    else "-")
        events.append(f)
    return F, sorted(events, key=lambda e: e["t"])


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
    F, events = classify_events(x, fs, p)

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
        res.update(t_mb=onset_sustained(ec, pk, fs, back + (pk - mb["idx_lf"]) / fs,
                                        floor=float(np.median(ec))) / fs, mb_idx=pk,
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


def nwave_duration(sig, fs, k_onset, up=8, hp_hz=100.0, max_tn_s=0.0008):
    """N-wave duration T_N (s) of one channel, sub-sample: front shock (steepest rise) to rear
    shock (steepest rise after the minimum). The signal around the onset sample k_onset is
    high-passed (hp_hz, removes offsets and wind) and upsampled x up, so T_N is not limited to
    whole samples (19.5 us at 51.2 kHz). The rear shock is searched within max_tn_s of the
    front: a ground reflection or echo follows later (0.6-3 ms in Cipatat). NaN if unclear."""
    from scipy.signal import resample_poly
    a, b = k_onset - int(0.005 * fs), k_onset + int(0.005 * fs)
    if a < 0 or b > len(sig):
        return float("nan")
    seg = sosfiltfilt(butter(2, hp_hz, "highpass", fs=fs, output="sos"), np.asarray(sig[a:b], float))
    i0 = k_onset - a
    w = resample_poly(seg[i0 - int(0.00025 * fs): i0 + int(0.0015 * fs)], up, 1)
    fsu = fs * up
    m = np.abs(w)
    if m.max() <= 0:
        return float("nan")
    w = w if w[int(np.argmax(m > 0.5 * m.max()))] >= 0 else -w      # first major excursion > 0
    d = np.diff(w)
    if len(d) < 10:
        return float("nan")
    dd = d[: max(5, int(0.6 * len(d)))]
    i_front = int(np.argmax(dd > 0.5 * dd.max()))           # EARLIEST steep rise (the rear is as steep)
    i_front += int(np.argmax(d[i_front:i_front + up // 2 + 1]))
    i_end = min(len(w), i_front + max(4, int(max_tn_s * fsu)))
    i_min = i_front + int(np.argmin(w[i_front:i_end]))
    if i_min >= len(d) - 1:
        return float("nan")
    i_rear = i_min + int(np.argmax(d[i_min:i_min + max(3, int(0.00015 * fsu))]))
    return float((i_rear - i_front) / fsu) if i_rear > i_front else float("nan")


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
        return lo_i + onset_sustained(e, pk, fs, back_s, frac, floor=float(np.percentile(e, 20)))
    return lo_i + onset(e, pk, fs, back_s, frac)


# --------------------------------------------------------------------------- several shots
SHOT_DEFAULTS = dict(
    direction_gate_deg=20.0,   # an event further than this from the shooter direction is an echo
    pair_tol_ms=8.0,           # SW -> MB delay tolerance between rounds (at least) ...
    pair_tol_frac=0.02,        # ... or this fraction of the delay
    max_dt_s=5.0,              # longest SW -> MB delay considered (2000 m: at most ~3.5 s)
    min_dt_s=0.030,            # shortest SW -> MB delay (a SW's own LF body follows within 25 ms)
    sw_train_db=12.0,          # SWs within this of the strongest form the burst train
    mb_tail_max_angle_deg=75.0,  # one SW, MB in its echo tail: accepted if SW and MB come from
                               # directions at most this far apart (Mach cone; needs an array)
    tail_hf_allowance_db=3.0,  # bursts: HF limit raised by this inside another round's SW tail
    min_shot_interval_s=0.090, # MB-only shots: of two blasts closer than this, the weaker is an
                               # echo (project limit 600 rpm = 100 ms, minus 10 % cadence spread)
)


def _xcorr_lag(a, b, fs, max_lag_s, up=8):
    """Delay (s) of b relative to a (positive: b later), positive peak of the cross-correlation."""
    from scipy.signal import resample_poly
    ua, ub = resample_poly(a, up, 1), resample_poly(b, up, 1)
    cc = np.correlate(ub, ua, mode="full")
    mid = len(ua) - 1
    m = int(max_lag_s * fs * up)
    w = cc[max(0, mid - m): mid + m + 1]
    return (int(np.argmax(w)) - min(m, mid)) / (fs * up)


def event_direction(x, fs, t0, band, mic_pos, c=343.0, pre=0.005, post=0.030):
    """Direction TO the source of the event starting near t0 (s): band-pass, cross-correlate
    every channel with channel 0 in [t0 - pre, t0 + post], least-squares plane wave.
    Returns (unit vector or None, rms residual in s, planar flag). For a planar array only the
    horizontal components are meaningful (vertical set to 0)."""
    mic_pos = np.asarray(mic_pos, dtype=float)
    n = x.shape[1]
    if mic_pos is None or n < 3 or len(mic_pos) != n:
        return None, None, None
    aperture = max(np.linalg.norm(mic_pos[i] - mic_pos[j]) for i in range(n) for j in range(n))
    max_lag = 1.2 * aperture / c
    i0, i1 = max(0, int((t0 - pre - max_lag) * fs)), min(len(x), int((t0 + post + max_lag) * fs))
    if i1 - i0 < 16:
        return None, None, None
    y = bandpass(x[i0:i1], fs, *band)
    k0, k1 = int(max_lag * fs), len(y) - int(max_lag * fs)
    seg = y[k0:k1] * np.hanning(k1 - k0)[:, None]
    tau = np.array([0.0] + [_xcorr_lag(seg[:, 0], seg[:, i], fs, max_lag) for i in range(1, n)])
    A = mic_pos[1:] - mic_pos[0]
    rhs = -c * tau[1:]
    planar = np.ptp(mic_pos[:, 2]) < 1e-6
    if planar:
        d2, *_ = np.linalg.lstsq(A[:, :2], rhs, rcond=None)
        d = np.array([d2[0], d2[1], 0.0])
    else:
        d, *_ = np.linalg.lstsq(A, rhs, rcond=None)
    res = float(np.sqrt(np.mean((A @ d - rhs) ** 2)) / c)
    nrm = np.linalg.norm(d)
    return (d / nrm if nrm > 0 else None), res, bool(planar)


def _angle_deg(u, v, planar=False):
    if u is None or v is None:
        return None
    if planar:
        u, v = np.array([u[0], u[1], 0.0]), np.array([v[0], v[1], 0.0])
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu == 0 or nv == 0:
        return None
    return float(np.degrees(np.arccos(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))))


def detect_shots(x, fs, params=None, mic_pos=None, c=343.0, verbose=False):
    """Several shots in one recording (bursts, repeated shots). Returns a dict:
      shots      list, in time order, of dicts like detect_sw_mb's result for ONE shot
                 (t_sw, t_mb, sw_idx, mb_idx, mb_band, event_class, sw_features, mb_features)
                 plus 'paired' (MB matched to this SW by the common delay), 'dt_s',
                 'mb_direction' / 'sw_direction' (unit vectors to the source, if an array is given)
      rejected   MB candidates left out, each with the reason: 'direction' (echo from another
                 direction), 'unpaired' (supersonic shots found, but no SW at the common delay),
                 'sw_body' (low-frequency part of a SW), 'min_shot_interval' (MB only: weaker of two
                 blasts closer than min_shot_interval_s)
      dt_common_s, shooter_direction, echo_rejection (how echoes were rejected)
    How rounds are told from echoes. Project limit: bursts of at most 600 rpm (rounds >= 100 ms
    apart), but echoes also arrive 25-200 ms after a shot, so a hold-off alone is not enough:
      1. Direction, when there is no SW (needs a synchronised array, mic_pos): rounds come from the shooter, echoes
         from reflectors. Each MB candidate's direction is measured; the shooter direction is the
         one with the most level-weighted support; candidates more than direction_gate_deg away
         are echoes.
      2. Supersonic shots: every round has its own SW and the SW -> MB delay is nearly the same
         for all rounds. The delay matched by the most SWs is taken (the earliest such delay: an
         echo of the MB train forms a second regular train, always later); each SW is paired with
         the MB at that delay (+- pair_tol); MB candidates without a SW partner are not shots.
         Directions are then only checked (direction_inconsistent), not used to drop a pair.
      3. No SW: of two MB candidates closer than min_shot_interval_s (90 ms) only the stronger is
         kept. Without an array, the later echoes cannot be told from rounds: every remaining
         candidate is reported, flagged 'echo_risk' if it follows a >= 3 dB stronger one within
         300 ms. Supersonic shots use no hold-off: SW pairing already separates the rounds.
    A delay is preferred whose MBs lie outside the SWs' echo tails; a single SW whose only
    candidates are inside its tail stays 'SW only' (unless mb_in_tail, or an array shows the SW
    and the MB within mb_tail_max_angle_deg of each other: shots at 50-100 m), but a delay
    inside the tails that several SWs share (e.g. shots at 90 m) is accepted."""
    p = dict(DEFAULT_PARAMS)
    p.update(SHOT_DEFAULTS)
    if params:
        p.update({k: v for k, v in params.items() if v is not None})
    x = np.asarray(x, dtype=float)
    if x.ndim == 1:
        x = x[:, None]
    dur = len(x) / fs
    F, events = classify_events(x, fs, p, k=max(25, int(12 * dur)))
    use_dir = mic_pos is not None and x.shape[1] >= 3 and len(mic_pos) == x.shape[1]
    gate = p["direction_gate_deg"]

    # ---- shockwaves: the burst train
    sws = [e for e in events if e["cls"] == "SW"]
    if sws:
        best = max(e["snr_hf"] for e in sws)
        sws = [e for e in sws if e["snr_hf"] >= best - p["sw_train_db"]]
    for e in sws:
        e["t_on"] = onset(F.hf, e["idx_hf"], fs, 0.002) / fs
        if use_dir:
            e["dir"], e["dir_res"], e["planar"] = event_direction(x, fs, e["t_on"], (1000.0, 8000.0),
                                                                  mic_pos, c, pre=0.0015, post=0.003)

    # ---- muzzle-blast candidates and their onsets
    rejected = []
    cands = []
    # In a burst a round's MB can fall 5-60 ms after another round's SW, where that SW's
    # reverberation adds high frequencies: there the HF limit is raised by tail_hf_allowance_db.
    # Pairing then decides; unpaired candidates are dropped.
    def in_sw_tail(e):
        return any(0.005 <= e["t"] - sw["t"] <= 0.060 for sw in sws)

    def is_cand(e):
        if e["cls"] == "MB?":
            return True
        return (e["cls"] == "-" and len(sws) >= 2 and in_sw_tail(e)
                and e["snr_lf"] >= p["mb_min_lf_snr_db"] and e["prominence"] >= p["mb_min_prominence_db"]
                and e["hf_excess"] <= p["mb_max_hf_excess_db"] + p.get("tail_hf_allowance_db", 0.0))

    for e in events:
        if not is_cand(e):
            continue
        # a SW's own low-frequency body follows it within 25 ms; min_dt_s (30 ms) keeps it from
        # being paired with that SW, while a round's MB that happens to fall just after another
        # round's SW (bursts) stays a candidate
        e["near_sw"] = any(-0.005 <= e["t"] - s["t"] <= 0.025 for s in sws)
        band = F.mb_bands[e["band_k"]]
        ec = causal_envelope(x, fs, band, 1e-3)
        a, b = max(0, e["idx_lf"] - int(0.01 * fs)), min(len(ec), e["idx_lf"] + int(0.02 * fs))
        pk = a + int(np.argmax(ec[a:b]))
        e["mb_idx"], e["mb_band"] = pk, band
        back = 0.01 + (pk - e["idx_lf"]) / fs
        prev_sw = [sw["idx_hf"] for sw in sws if sw["idx_hf"] < pk]
        if prev_sw:                                    # keep the onset search out of a preceding SW
            back = max(0.001, min(back, (pk - max(prev_sw)) / fs - 0.003))
        e["t_on"] = onset_sustained(ec, pk, fs, back, floor=float(np.median(ec))) / fs
        cands.append(e)

    for e in cands:
        e["score"] = e["prominence"] - 0.3 * float(np.clip(e["hf_excess"], -15, 15))
        if use_dir:
            e["dir"], e["dir_res"], e["planar"] = event_direction(x, fs, e["t_on"], e["mb_band"], mic_pos, c)
    planar = bool(cands[0].get("planar")) if cands else False

    def within(u, v):
        ang = _angle_deg(u, v, planar)
        return ang is not None and ang <= gate

    def mean_dir(items):
        v = np.sum([10 ** (f["snr_lf"] / 20) * f["dir"] for f in items], axis=0)
        return v / np.linalg.norm(v)

    shooter_dir, shots, dt_common = None, [], None
    echo_rejection = "none (no array direction)"
    if sws:
        # ---- supersonic: pair every SW with the MB at the common SW -> MB delay.
        # SW echo tails (until the 50 Hz-4 kHz level stays within 10 dB of the noise for 30 ms):
        # an MB candidate inside one is less likely (low-frequency echoes of the SW)
        env_bb = envelope(bandpass(x, fs, 50, 4000), fs, 0.01)
        quiet = env_bb < np.median(env_bb) * 10 ** (10 / 20)
        cs = np.concatenate([[0], np.cumsum(quiet)])
        w30 = int(0.03 * fs)
        tails = []
        for sw in sws:
            k = sw["idx_hf"] + int(0.01 * fs)
            while k + w30 < len(env_bb) and cs[k + w30] - cs[k] < w30:
                k += int(0.002 * fs)
            tails.append((sw["idx_hf"], k))
        for m in cands:
            m["in_sw_tail"] = any(a0 <= m["idx_lf"] < a1 for a0, a1 in tails)
        best = None
        for s0 in sws:
            for m0 in cands:
                d = m0["t_on"] - s0["t_on"]
                if not p["min_dt_s"] <= d <= p["max_dt_s"]:
                    continue
                tol = max(p["pair_tol_ms"] * 1e-3, p["pair_tol_frac"] * d)
                matched = []
                for s in sws:
                    near = [m for m in cands if abs((m["t_on"] - s["t_on"]) - d) <= tol]
                    if near:
                        matched.append(min(near, key=lambda m: abs((m["t_on"] - s["t_on"]) - d)))
                # Several SWs: most SWs paired, then the EARLIEST delay (an echo of the MB train forms a second, equally regular train,
                # always later). One SW: as detect_sw_mb -- outside the tail, best blast score.
                if len(sws) > 1:     # (in a burst the SW tails cover the rounds' MBs too)
                    key = (len(matched), -d)
                else:
                    key = (len(matched), -sum(m["in_sw_tail"] for m in matched),
                           round(sum(m["score"] for m in matched) / len(matched) / 3.0), -d)
                if best is None or key > best[0]:
                    best = (key, d, matched)
        used = set()
        if best is not None and best[0][0] < 2 and best[0][1] == -best[0][0] and not p["mb_in_tail"]:
            # One SW and its only candidate MBs inside the SW echo tail (shots closer than
            # ~100 m, or low-frequency echoes of the SW). With an array the pair is checked: the
            # SW comes from at most mb_tail_max_angle_deg away from the shooter (Mach-cone
            # angle; 49-57 deg in Cipatat), while a SW echo can come from anywhere.
            ok = False
            if use_dir and best[2]:
                ang = _angle_deg(sws[0].get("dir"), best[2][0].get("dir"), False)
                ok = ang is not None and ang <= p.get("mb_tail_max_angle_deg", 0.0)
            if not ok:
                best = None     # SW only
        if best is not None:
            d = best[1]
            tol = max(p["pair_tol_ms"] * 1e-3, p["pair_tol_frac"] * d)
            hits = []
            for s in sws:
                near = [m for m in cands if id(m) not in used and abs((m["t_on"] - s["t_on"]) - d) <= tol]
                mb = min(near, key=lambda m: abs((m["t_on"] - s["t_on"]) - d)) if near else None
                if mb is not None:
                    used.add(id(mb))
                    hits.append(mb["t_on"] - s["t_on"])
                shots.append([s, mb])
            dt_common = float(np.median(hits)) if hits else None
        else:
            shots = [[s, None] for s in sws]
        for m in cands:
            if id(m) not in used:
                rejected.append((m, "sw_body" if m.get("near_sw") else "unpaired"))
        echo_rejection = "SW pairing (common SW -> MB delay)"
        paired = [sm[1] for sm in shots if sm[1] is not None and sm[1].get("dir") is not None]
        if use_dir and paired:
            # The pairing is the stronger evidence: in a fast burst an MB window can overlap
            # another round's SW, which bends its measured direction. A paired MB far from the
            # other rounds' direction is only flagged, not dropped.
            centre = max(paired, key=lambda e: sum(10 ** (f["snr_lf"] / 20) for f in paired
                                                   if within(e["dir"], f["dir"])))
            shooter_dir = mean_dir([f for f in paired if within(centre["dir"], f["dir"])])
            for sm in shots:
                m = sm[1]
                if m is None or m.get("dir") is None:
                    continue
                m["angle_to_shooter_deg"] = _angle_deg(shooter_dir, m["dir"], planar)
                m["direction_inconsistent"] = not within(shooter_dir, m["dir"])
    else:
        # ---- no SW: direction gate, the shooter direction being the one with the most
        # level-weighted support (rounds share it, echoes scatter)
        good = [e for e in cands if e.get("dir") is not None]
        if use_dir and good:
            centre = max(good, key=lambda e: sum(10 ** (f["snr_lf"] / 20) for f in good
                                                 if within(e["dir"], f["dir"])))
            shooter_dir = mean_dir([f for f in good if within(centre["dir"], f["dir"])])
            keep = []
            for e in cands:
                e["angle_to_shooter_deg"] = _angle_deg(shooter_dir, e.get("dir"), planar)
                if e["angle_to_shooter_deg"] is None or e["angle_to_shooter_deg"] > gate:
                    rejected.append((e, "direction"))
                else:
                    keep.append(e)
            cands = keep
            echo_rejection = f"direction (gate {gate:g} deg)"
        # Project limit: bursts of at most 600 rpm, so two real blasts are >= 100 ms apart. Of two
        # candidates closer than min_shot_interval_s the stronger is kept, the other is an echo.
        gap = p.get("min_shot_interval_s", 0.0)
        if gap > 0:
            kept = []
            for e in sorted(cands, key=lambda e: -e["snr_lf"]):
                if any(abs(e["t"] - k["t"]) < gap for k in kept):
                    rejected.append((e, "min_shot_interval"))
                else:
                    kept.append(e)
            cands = sorted(kept, key=lambda e: e["t"])
        for m in cands:
            shots.append([None, m])
        if not use_dir:
            for i, m in enumerate(cands):
                m["echo_risk"] = any(0 < m["t"] - o["t"] <= 0.3 and o["snr_lf"] - m["snr_lf"] >= 3.0
                                     for o in cands[:i])

    out = []
    for s, m in sorted(shots, key=lambda sm: (sm[0] or sm[1])["t_on"]):
        if s is None and m is None:
            continue
        r = dict(t_sw=s["t_on"] if s else None, sw_idx=s["idx_hf"] if s else None,
                 t_mb=m["t_on"] if m else None, mb_idx=m["mb_idx"] if m else None,
                 mb_band=m["mb_band"] if m else None, sw_features=s, mb_features=m,
                 paired=bool(s and m), dt_s=(m["t_on"] - s["t_on"]) if (s and m) else None,
                 sw_direction=(s.get("dir").tolist() if s and s.get("dir") is not None else None),
                 mb_direction=(m.get("dir").tolist() if m and m.get("dir") is not None else None),
                 angle_to_shooter_deg=(m.get("angle_to_shooter_deg") if m else None),
                 direction_inconsistent=bool(m.get("direction_inconsistent")) if m else False,
                 echo_risk=bool(m.get("echo_risk")) if m else False)
        r["event_class"] = "SW+MB" if s and m else ("SW only" if s else "MB only")
        out.append(r)
    if verbose:
        for i, r in enumerate(out):
            f = lambda v: "--" if v is None else f"{v:.4f}"
            print(f"    shot {i + 1}: {r['event_class']:7s} SW {f(r['t_sw'])}  MB {f(r['t_mb'])}"
                  + (f"  dt {r['dt_s'] * 1e3:.1f} ms" if r["dt_s"] else "")
                  + (f"  {r['angle_to_shooter_deg']:.1f} deg" if r["angle_to_shooter_deg"] is not None else ""))
        for e, why in rejected:
            print(f"    rejected {e['t']:.4f} s ({why})")
    return dict(shots=out, rejected=[dict(t=e["t"], reason=why, snr_lf=e["snr_lf"],
                                          angle_to_shooter_deg=e.get("angle_to_shooter_deg"))
                                     for e, why in rejected],
                dt_common_s=dt_common,
                shooter_direction=(shooter_dir.tolist() if shooter_dir is not None else None),
                echo_rejection=echo_rejection, n_events=len(events), params=p)
