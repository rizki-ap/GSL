#!/usr/bin/env python3
"""
gs_gen_add_noise.py
===================
STAGE 2 of 3. Reads the noiseless per-mic signal written by gs_gen_clean_signal.py and adds the
ENVIRONMENT: noise, interfering sounds, echoes and reverberation. Nothing hardware-specific:
microphone response, analog gain and ADC quantization are Stage 3 (gs_gen_apply_adc.py).
Everything fixed by the scene geometry (direct shockwave and blast, ground reflection, air
absorption) is Stage 1.

Run Stage 1 once, then this as often as you like with different [noise] settings.

Usage
-----
  python gs_gen_add_noise.py gunshot_sim_clean.wav                     # gen_config.ini
  python gs_gen_add_noise.py gunshot_sim_clean.wav my_config.ini
  python gs_gen_add_noise.py gunshot_sim_clean.wav my_config.ini -o my_shot
  python gs_gen_add_noise.py gunshot_sim_clean.wav -layers colored,wind,tones -snr 20

Noise models ([noise] model)
----------------------------
  simple     white noise of noise_floor_pa on each mic (unchanged)
  realistic  room reverb (rt60) + independent colored noise (noise_rms_pa, noise_slope) on each
             mic (unchanged; kept for old configurations)
  layers     any combination of the layers below ([noise] layers = colored, wind, ...), each
             configured in its own [noise.<layer>] section

Layers (model = layers)
-----------------------
  Convolutive (applied to the Stage-1 signal):
    room      diffuse reverberation tail: rt60, drr_db (direct-to-reverberant ratio)
    echoes    a few discrete outdoor echoes (treeline, buildings): random (n, delay_ms range,
              gain range, lowpass_hz) or given (list = delay_ms, gain, azimuth, elevation; ...).
              Each echo arrives from its own direction (its own TDOAs).
  Additive:
    colored   broadband noise, piecewise power-law spectrum (breaks_hz, slopes), part of it
              identical on all mics (common_fraction, as in real recordings)
    wind      gusty low-frequency noise, different at each mic (corner_hz, slope, gust_depth,
              gust_rate_hz)
    tones     harmonic series from fixed directions (machinery):
              sources = f0_hz, n_harmonics, level_pa, azimuth_deg, elevation_deg; ...
    chirps    pulsed narrow-band tones from fixed directions (insects):
              sources = freq_hz, level_pa, azimuth_deg, elevation_deg, pulse_rate_hz, duty; ...
    diffuse   noise arriving from all directions (n_directions plane waves): the physically
              correct inter-mic coherence for the array
    real      noise from real recordings (files, scale to Pa); loud impulsive parts are cut out
              (exclude_events)
    impulses  interfering impulsive sounds (rate_per_s, level_pa range, kinds = blast, click,
              nwave), from random directions, kept guard_s away from the true events; their
              times are listed in the .json for false-alarm tests

  Directions: azimuth counter-clockwise from +x of the array frame, elevation up from the
  horizontal (the generator's convention). For PILAR-V (x = north, y = west) the compass bearing
  is -azimuth.

Level ([noise] level_mode)
--------------------------
  absolute  every additive layer at the level given in its section (Pa)
  snr       the additive layers (except impulses) are scaled together so that
            20 log10(reference peak / total noise rms) = snr_db; snr_ref = shockwave,
            muzzle_blast (Stage-1 ground-truth peak, mean over mics) or peak (max |signal|)

Output
------
  <output>_noisy.wav   N-channel, 32-bit FLOAT PCM, Pa units (not yet mic-filtered or quantized)
  <output>_noisy.json  Stage-1 ground truth + noise details (per-layer settings, realized random
                       values, rms per layer, achieved SNR)
The channel count follows the input .wav; -mic_config is optional (default: the preset recorded
in the Stage-1 .json). Reads only [noise], [noise.*], [atmosphere] (speed of sound, if the
Stage-1 .json lacks it) and [output].
"""

import argparse
import json
from datetime import datetime, timezone

import numpy as np
from scipy.io import wavfile
from scipy.signal import butter, fftconvolve, resample_poly, sosfiltfilt

import gs_gen_physic as gp

LAYERS = ("room", "echoes", "colored", "wind", "tones", "chirps", "diffuse", "real", "impulses")
CONVOLUTIVE = ("room", "echoes")
SCALED_BY_SNR = ("colored", "wind", "tones", "chirps", "diffuse", "real")

# Defaults of each [noise.<layer>] section. Levels roughly follow the Cipatat test range
# (sound/Cipatat/metadata.md, background noise of shot-03, scaled down to a quieter day).
LAYER_DEFAULTS = {
    "room":     dict(rt60="1.25", drr_db="10.0"),
    "echoes":   dict(n="3", delay_ms="20, 400", gain="0.05, 0.25", lowpass_hz="3000",
                     elevation_deg="0, 20", list=""),
    "colored":  dict(rms_pa="0.02", breaks_hz="1000", slopes="-1.5, -2.4", fmin_hz="20",
                     common_fraction="0.3"),
    "wind":     dict(rms_pa="0.1", corner_hz="20", slope="-2.0", gust_depth="0.6", gust_rate_hz="0.3"),
    "tones":    dict(sources="83.6, 3, 0.02, -75, -40; 412, 2, 0.01, -22, -25", drift_pct="0.1"),
    "chirps":   dict(sources="3800, 0.003, -115, 0, 3.0, 0.3; 4180, 0.003, -115, 0, 3.0, 0.3"),
    "diffuse":  dict(rms_pa="0.01", slope="-2.0", fmin_hz="20", n_directions="64"),
    "real":     dict(files="", scale="1.0", rms_pa="0", exclude_events="true"),
    "impulses": dict(rate_per_s="0.2", level_pa="0.5, 5", kinds="blast, click, nwave", guard_s="0.05",
                     elevation_deg="-10, 30"),
}


# =============================================================================================
# I/O
# =============================================================================================
def load_clean_wav(wav_path, json_path=None):
    """Load Stage 1's float32 clean .wav (+ its .json, looked up by convention if not given)."""
    fs, data = wavfile.read(wav_path)
    if data.dtype != np.float32:
        print(f"WARNING: expected float32 clean .wav, got {data.dtype}. "
              f"Proceeding, but check this is really a Stage-1 output.")
    if data.ndim == 1:
        data = data[:, None]
    signals = [data[:, i].astype(np.float64) for i in range(data.shape[1])]

    if json_path is None:
        json_path = wav_path.rsplit(".wav", 1)[0] + ".json"
        if not json_path.endswith("_clean.json") and wav_path.endswith("_clean.wav"):
            json_path = wav_path[:-len("_clean.wav")] + "_clean.json"
    clean_meta = None
    try:
        with open(json_path) as f:
            clean_meta = json.load(f)
    except FileNotFoundError:
        print(f"WARNING: no companion metadata file found at '{json_path}' -- "
              f"output .json will be missing ground-truth physics info.")
    return signals, fs, clean_meta


def write_noisy_wav(path, signals, fs):
    """Float32 PCM, still in Pa units (lossless hand-off to Stage 3)."""
    n = min(len(s) for s in signals)
    wavfile.write(path, fs, np.stack([s[:n] for s in signals], axis=1).astype(np.float32))


# =============================================================================================
# small helpers
# =============================================================================================
def _floats(s, n=None):
    v = [float(p) for p in str(s).replace(";", ",").split(",") if p.strip()]
    if n is not None and len(v) != n:
        raise ValueError(f"expected {n} numbers, got '{s}'")
    return v


def _sources(s, n):
    """'a, b, c; d, e, f' -> list of n-tuples."""
    return [_floats(row, n) for row in str(s).split(";") if row.strip()]


def _layer_params(cfg, name):
    p = dict(LAYER_DEFAULTS[name])
    p.update({k.lower(): v for k, v in cfg["noise"].get("layer_params", {}).get(name, {}).items()})
    unknown = set(p) - set(LAYER_DEFAULTS[name])
    if unknown:
        raise ValueError(f"[noise.{name}]: unknown key(s) {sorted(unknown)} -- valid: {sorted(LAYER_DEFAULTS[name])}")
    return p


def unit_vector(az_deg, el_deg):
    az, el = np.radians(az_deg), np.radians(el_deg)
    return np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])


def plane_wave_delays(mic_pos, u, c):
    """Arrival delay (s) at each mic, relative to the array origin, of a plane wave coming FROM
    direction u: mics closer to the source hear it earlier."""
    return -(np.asarray(mic_pos) @ u) / c


def random_direction(rng, el_range=(-90.0, 90.0)):
    az = rng.uniform(0.0, 360.0)
    lo, hi = np.sin(np.radians(el_range[0])), np.sin(np.radians(el_range[1]))
    el = np.degrees(np.arcsin(rng.uniform(lo, hi)))
    return float(az), float(el)


def shaped_spectrum(n, fs, rng, shape):
    """Random complex spectrum (rfft bins) with amplitude sqrt(PSD) = shape(f)."""
    f = np.fft.rfftfreq(n, 1.0 / fs)
    X = (rng.standard_normal(len(f)) + 1j * rng.standard_normal(len(f))) / np.sqrt(2)
    return f, X * shape(f)


def piecewise_psd_amp(breaks, slopes, fmin):
    """Amplitude (sqrt of PSD) of a continuous piecewise power law: PSD ~ f^slopes[k] between
    consecutive breaks (Hz); zero below fmin."""
    bounds = [max(float(fmin), 1e-3)] + [float(b) for b in breaks] + [np.inf]

    def shape(f):
        f = np.asarray(f, dtype=float)
        a = np.zeros_like(f)
        base = 0.0                                   # log PSD at the start of the segment
        for k, sl in enumerate(slopes):
            lo, hi = bounds[k], bounds[k + 1]
            seg = (f >= lo) & (f < hi)
            a[seg] = np.exp(0.5 * (base + sl * np.log(f[seg] / lo)))
            if np.isfinite(hi):
                base += sl * np.log(hi / lo)
        return a
    return shape


def _rms(x):
    return float(np.sqrt(np.mean(np.square(x)))) if np.size(x) else 0.0


def _scale_to_rms(x, rms):
    r = _rms(x)
    return x * (rms / r) if r > 0 else x


# =============================================================================================
# layers -- each returns (array [n, n_mics] in Pa, info dict)
# =============================================================================================
def layer_colored(n, fs, mic_pos, c, rng, p):
    m = len(mic_pos)
    breaks, slopes = _floats(p["breaks_hz"]), _floats(p["slopes"])
    if len(slopes) != len(breaks) + 1:
        raise ValueError("[noise.colored]: slopes needs one value more than breaks_hz")
    shape = piecewise_psd_amp(breaks, slopes, float(p["fmin_hz"]))
    cf = float(p["common_fraction"])
    if not 0.0 <= cf <= 1.0:
        raise ValueError("[noise.colored]: common_fraction must be in [0, 1]")
    _, Xc = shaped_spectrum(n, fs, rng, shape)
    common = np.fft.irfft(Xc, n)
    out = np.zeros((n, m))
    for i in range(m):
        _, Xi = shaped_spectrum(n, fs, rng, shape)
        out[:, i] = np.sqrt(cf) * _scale_to_rms(common, 1.0) + np.sqrt(1 - cf) * _scale_to_rms(np.fft.irfft(Xi, n), 1.0)
    out = _scale_to_rms(out, float(p["rms_pa"]))
    return out, dict(breaks_hz=breaks, slopes=slopes, common_fraction=cf)


def layer_wind(n, fs, mic_pos, c, rng, p):
    m = len(mic_pos)
    fc, slope = float(p["corner_hz"]), float(p["slope"])
    shape = lambda f: np.where(f >= 1.0, 1.0 / np.sqrt(1.0 + (np.maximum(f, 1e-6) / fc) ** (-slope)), 0.0)
    # gusts: one slowly varying envelope for the whole array (gusts are larger than the array),
    # independent turbulence at each microphone
    rate, depth = float(p["gust_rate_hz"]), float(p["gust_depth"])
    g = rng.standard_normal(n)
    if rate > 0:
        sos = butter(2, min(rate, 0.45 * fs), btype="low", fs=fs, output="sos")
        g = sosfiltfilt(sos, g)
    g = (g - g.mean()) / (g.std() + 1e-12)
    env = np.maximum(0.1, 1.0 + depth * g)
    out = np.zeros((n, m))
    for i in range(m):
        _, X = shaped_spectrum(n, fs, rng, shape)
        out[:, i] = np.fft.irfft(X, n) * env
    out = _scale_to_rms(out, float(p["rms_pa"]))
    return out, dict(corner_hz=fc, slope=slope, gust_depth=depth, gust_rate_hz=rate)


def layer_tones(n, fs, mic_pos, c, rng, p):
    t = np.arange(n) / fs
    drift = float(p["drift_pct"]) / 100.0
    out = np.zeros((n, len(mic_pos)))
    info = []
    for f0, nh, level, az, el in _sources(p["sources"], 5):
        tau = plane_wave_delays(mic_pos, unit_vector(az, el), c)
        # slow random frequency wander (engine speed), the same for all harmonics
        w = np.cumsum(rng.standard_normal(n)) / np.sqrt(fs)
        w = sosfiltfilt(butter(1, 0.5, btype="low", fs=fs, output="sos"), w)
        w = w / (np.abs(w).max() + 1e-12)
        finst = f0 * (1.0 + drift * w)
        phase = 2 * np.pi * np.cumsum(finst) / fs + rng.uniform(0, 2 * np.pi)
        for i in range(len(mic_pos)):
            ph = np.interp(t - tau[i], t, phase)
            for k in range(1, int(nh) + 1):
                out[:, i] += level * np.sqrt(2) / k * np.cos(k * ph + k * 0.7)
        info.append(dict(f0_hz=f0, n_harmonics=int(nh), level_pa_rms_fundamental=level,
                         azimuth_deg=az, elevation_deg=el))
    return out, dict(sources=info, drift_pct=drift * 100)


def layer_chirps(n, fs, mic_pos, c, rng, p):
    t = np.arange(n) / fs
    out = np.zeros((n, len(mic_pos)))
    info = []
    for f, level, az, el, rate, duty in _sources(p["sources"], 6):
        tau = plane_wave_delays(mic_pos, unit_vector(az, el), c)
        # pulse train with jittered period, raised-cosine gates
        gate = np.zeros(n)
        period = 1.0 / rate
        tk = rng.uniform(0, period)
        L = max(4, int(duty * period * fs))
        win = np.hanning(L)
        while tk < n / fs:
            i0 = int(tk * fs)
            seg = win[: max(0, min(L, n - i0))]
            gate[i0:i0 + len(seg)] = np.maximum(gate[i0:i0 + len(seg)], seg)
            tk += period * rng.uniform(0.8, 1.2)
        ph0 = rng.uniform(0, 2 * np.pi)
        for i in range(len(mic_pos)):
            g = np.interp(t - tau[i], t, gate, left=0.0, right=0.0)
            out[:, i] += level * np.sqrt(2) * g * np.cos(2 * np.pi * f * (t - tau[i]) + ph0)
        info.append(dict(freq_hz=f, level_pa_rms_in_pulse=level, azimuth_deg=az, elevation_deg=el,
                         pulse_rate_hz=rate, duty=duty))
    return out, dict(sources=info)


def layer_diffuse(n, fs, mic_pos, c, rng, p):
    m, nd = len(mic_pos), int(p["n_directions"])
    shape = piecewise_psd_amp([], [float(p["slope"])], float(p["fmin_hz"]))
    f = np.fft.rfftfreq(n, 1.0 / fs)
    Y = np.zeros((m, len(f)), dtype=complex)
    for _ in range(nd):
        u = unit_vector(*random_direction(rng))
        _, X = shaped_spectrum(n, fs, rng, shape)
        tau = plane_wave_delays(mic_pos, u, c)
        for i in range(m):
            Y[i] += X * np.exp(-2j * np.pi * f * tau[i])
    out = np.stack([np.fft.irfft(Y[i], n) for i in range(m)], axis=1)
    out = _scale_to_rms(out, float(p["rms_pa"]))
    return out, dict(n_directions=nd, slope=float(p["slope"]))


def _quiet_segments(x, fs, exclude):
    """(start, end) sample ranges without loud impulsive events: 20 ms frames whose energy is
    more than 10 dB above the median frame are cut out, from 50 ms before to 1 s after."""
    if not exclude:
        return [(0, len(x))]
    hp = sosfiltfilt(butter(2, 100.0, btype="high", fs=fs, output="sos"), x, axis=0)
    L = int(0.02 * fs)
    nfr = len(x) // L
    e = np.array([np.mean(hp[k * L:(k + 1) * L] ** 2) for k in range(nfr)])
    bad = np.zeros(len(x), bool)
    for k in np.where(e > 10 * np.median(e))[0]:
        bad[max(0, k * L - int(0.05 * fs)):min(len(x), (k + 1) * L + int(1.0 * fs))] = True
    segs, i = [], 0
    while i < len(x):
        if bad[i]:
            i += 1
            continue
        j = i
        while j < len(x) and not bad[j]:
            j += 1
        if j - i >= int(0.25 * fs):
            segs.append((i, j))
        i = j
    return segs


def layer_real(n, fs, mic_pos, c, rng, p):
    import soundfile as sf
    m = len(mic_pos)
    files = [s.strip() for s in p["files"].split(",") if s.strip()]
    if not files:
        raise ValueError("[noise.real]: give files = noise1.wav, noise2.wav, ...")
    scale, exclude = float(p["scale"]), str(p["exclude_events"]).lower() in ("1", "true", "yes", "on")
    pieces = []
    for fpath in files:
        x, fsf = sf.read(fpath, dtype="float64", always_2d=True)
        if x.shape[1] < m:
            raise ValueError(f"[noise.real]: {fpath} has {x.shape[1]} channels, need {m}")
        x = x[:, :m] * scale
        if fsf != fs:
            from math import gcd
            g = gcd(int(fs), int(fsf))
            x = resample_poly(x, int(fs) // g, int(fsf) // g, axis=0)
        for a, b in _quiet_segments(x, fs, exclude):
            pieces.append((fpath, a, b, x[a:b]))
    if not pieces:
        raise ValueError("[noise.real]: no quiet segment of at least 0.25 s in the given files")
    out, used, k, fade = np.zeros((0, m)), [], 0, int(0.01 * fs)
    order = rng.permutation(len(pieces))
    while len(out) < n:
        fpath, a, b, seg = pieces[order[k % len(order)]]
        k += 1
        start = int(rng.integers(0, max(1, len(seg) - n))) if len(seg) > n else 0
        seg = seg[start:start + n].copy()
        if len(out) >= fade and len(seg) > 2 * fade:          # crossfade the joint
            r = np.linspace(0, 1, fade)[:, None]
            out[-fade:] = out[-fade:] * (1 - r) + seg[:fade] * r
            seg = seg[fade:]
        out = np.vstack([out, seg])
        used.append(dict(file=fpath, start_s=(a + start) / fs, duration_s=len(seg) / fs))
    out = out[:n] - out[:n].mean(axis=0)
    if float(p["rms_pa"]) > 0:
        out = _scale_to_rms(out, float(p["rms_pa"]))
    return out, dict(segments=used, scale=scale, exclude_events=exclude)


def _impulse_waveform(kind, fs, rng, level):
    if kind == "blast":
        tp = rng.uniform(0.002, 0.010)
        tt = np.arange(int(6 * tp * fs)) / fs
        w = (1 - tt / tp) * np.exp(-tt / tp)
        par = dict(t_pos_s=tp)
    elif kind == "click":
        f, tau = rng.uniform(1000, 5000), rng.uniform(0.0002, 0.001)
        tt = np.arange(int(8 * tau * fs)) / fs
        w = np.sin(2 * np.pi * f * tt) * np.exp(-tt / tau)
        par = dict(freq_hz=f, decay_s=tau)
    elif kind == "nwave":
        T = rng.uniform(0.0001, 0.0004)
        tt = np.arange(int(T * fs) + 2) / fs
        w = np.where(tt < T, 1 - 2 * tt / T, 0.0)
        par = dict(T_N_s=T)
    else:
        raise ValueError(f"[noise.impulses]: unknown kind '{kind}' -- use blast, click, nwave")
    w = w / (np.abs(w).max() + 1e-12) * level
    return w, par


def layer_impulses(n, fs, mic_pos, c, rng, p, event_times):
    m = len(mic_pos)
    rate, guard = float(p["rate_per_s"]), float(p["guard_s"])
    lo, hi = _floats(p["level_pa"], 2)
    kinds = [k.strip() for k in p["kinds"].split(",") if k.strip()]
    el_rng = tuple(_floats(p["elevation_deg"], 2))
    out = np.zeros((n, m))
    info, t = [], 0.0
    while rate > 0:
        t += rng.exponential(1.0 / rate)
        if t >= n / fs:
            break
        if any(abs(t - te) < guard for te in event_times):
            continue
        kind = kinds[int(rng.integers(len(kinds)))]
        level = float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
        w, par = _impulse_waveform(kind, fs, rng, level)
        az, el = random_direction(rng, el_rng)
        tau = plane_wave_delays(mic_pos, unit_vector(az, el), c)
        pad = int(0.003 * fs)
        snippet = np.concatenate([np.zeros(pad), w, np.zeros(pad)])
        i0 = int(t * fs) - pad
        for i in range(m):
            y = gp.fractional_delay(snippet, tau[i] * fs + (t * fs - int(t * fs)))
            a, b = max(0, i0), min(n, i0 + len(y))
            if b > a:
                out[a:b, i] += y[a - i0:b - i0]
        info.append(dict(t_s=round(t, 6), kind=kind, peak_pa=level, azimuth_deg=az, elevation_deg=el, **par))
    return out, dict(events=info, rate_per_s=rate)


def apply_room(dry, fs, rng, p):
    rt60, drr = float(p["rt60"]), float(p["drr_db"])
    out = np.zeros_like(dry)
    for i in range(dry.shape[1]):
        ir = gp.make_room_ir(rt60, fs, seed=int(rng.integers(2 ** 31)))
        tail = ir[1:]
        e_tail = np.sum(tail ** 2)
        if e_tail > 0:
            ir[1:] = tail * np.sqrt(10 ** (-drr / 10.0) / e_tail)   # ir[0] = 1 (direct)
        out[:, i] = fftconvolve(dry[:, i], ir, mode="full")[:len(dry)]
    return out, dict(rt60=rt60, drr_db=drr)


def apply_echoes(dry, fs, mic_pos, c, rng, p):
    """Discrete echoes of the whole Stage-1 signal: echo k = gain x low-passed channel-1 signal,
    delayed by delay_ms, arriving from its own direction (per-mic plane-wave delays)."""
    n, m = dry.shape
    if p["list"].strip():
        echoes = [dict(delay_ms=d, gain=g, azimuth_deg=az, elevation_deg=el)
                  for d, g, az, el in _sources(p["list"], 4)]
    else:
        d_lo, d_hi = _floats(p["delay_ms"], 2)
        g_lo, g_hi = _floats(p["gain"], 2)
        el_rng = tuple(_floats(p["elevation_deg"], 2))
        echoes = []
        for _ in range(int(p["n"])):
            az, el = random_direction(rng, el_rng)
            echoes.append(dict(delay_ms=float(rng.uniform(d_lo, d_hi)), gain=float(rng.uniform(g_lo, g_hi)),
                               azimuth_deg=az, elevation_deg=el))
    src = dry[:, 0]
    lp = float(p["lowpass_hz"])
    if 0 < lp < 0.45 * fs:
        src = sosfiltfilt(butter(2, lp, btype="low", fs=fs, output="sos"), src)
    out = np.zeros_like(dry)
    for e in echoes:
        tau = plane_wave_delays(mic_pos, unit_vector(e["azimuth_deg"], e["elevation_deg"]), c)
        for i in range(m):
            out[:, i] += e["gain"] * gp.fractional_delay(src, (e["delay_ms"] * 1e-3 + tau[i] - tau[0]) * fs)
    return out, dict(echoes=echoes, lowpass_hz=lp, source="channel 1 of the Stage-1 signal")


# =============================================================================================
# noise models
# =============================================================================================
def _event_times(gt):
    """Arrival times (s, relative to the start of the signal) of the true events."""
    if not gt:
        return []
    t0, out = gt.get("t_master_start_s", 0.0), []
    for pm in gt.get("per_mic", []):
        for key in ("shockwave", "muzzle_blast", "shockwave_reflected", "muzzle_blast_reflected"):
            t = (pm.get(key) or {}).get("t_arrive_s")
            if t is not None:
                out.append(t - t0)
    return out


def _reference_peak(gt, dry, ref):
    if ref == "peak" or not gt:
        return float(np.abs(dry).max()), "peak"
    key = ref
    vals = [pm[key]["peak_dP_pa"] for pm in gt.get("per_mic", [])
            if pm.get(key, {}).get("generated") and pm[key].get("peak_dP_pa") is not None]
    if not vals:
        print(f"WARNING: no {ref} in the Stage-1 signal -- SNR referenced to the signal peak instead")
        return float(np.abs(dry).max()), "peak"
    return float(np.mean(vals)), ref


def apply_layers(ideal_signals, fs, cfg, gt, mic_pos, c):
    noi, seed = cfg["noise"], cfg["output"]["seed"]
    layers = noi["layers"]
    bad = [x for x in layers if x not in LAYERS]
    if bad:
        raise ValueError(f"[noise] layers: unknown {bad} -- valid: {', '.join(LAYERS)}")
    dry = np.stack(ideal_signals, axis=1)
    n, m = dry.shape
    if len(mic_pos) != m:
        raise ValueError(f"{len(mic_pos)} mic positions for {m} channels")
    rngs = {name: np.random.default_rng([seed, k]) for k, name in enumerate(LAYERS)}
    info = {}

    wet = dry.copy()
    for name in CONVOLUTIVE:
        if name in layers:
            p = _layer_params(cfg, name)
            if name == "room":
                wet, inf = apply_room(dry, fs, rngs[name], p)
                inf["rms_pa"] = _rms(wet - dry)
            else:
                y, inf = apply_echoes(dry, fs, mic_pos, c, rngs[name], p)
                wet = wet + y
                inf["rms_pa"] = _rms(y)
            info[name] = inf

    additive, scaled = {}, np.zeros((n, m))
    for name in layers:
        if name in CONVOLUTIVE:
            continue
        p = _layer_params(cfg, name)
        if name == "impulses":
            y, inf = layer_impulses(n, fs, mic_pos, c, rngs[name], p, _event_times(gt))
        else:
            y, inf = globals()[f"layer_{name}"](n, fs, mic_pos, c, rngs[name], p)
            scaled += y
        additive[name] = y
        inf["rms_pa"] = _rms(y)
        info[name] = inf

    ref_peak, ref_used = _reference_peak(gt, dry, noi["snr_ref"])
    gain = 1.0
    if noi["level_mode"] == "snr":
        r = _rms(scaled)
        if r > 0:
            gain = ref_peak / r / 10 ** (noi["snr_db"] / 20.0)
            for name in SCALED_BY_SNR:
                if name in additive:
                    additive[name] *= gain
                    info[name]["rms_pa"] *= gain
        else:
            print("WARNING: level_mode = snr, but no noise layer to scale")
    total = wet + sum(additive.values()) if additive else wet
    noise_only = sum((additive[k] for k in additive if k in SCALED_BY_SNR), np.zeros((n, m)))
    r = _rms(noise_only)
    meta = dict(model="layers", layers=layers, level_mode=noi["level_mode"], snr_gain=gain,
                noise_rms_pa=r, noise_rms_pa_per_mic=[_rms(noise_only[:, i]) for i in range(m)],
                snr_reference=ref_used, snr_reference_peak_pa=ref_peak,
                snr_achieved_db=(20 * np.log10(ref_peak / r) if r > 0 else None),
                speed_of_sound_mps=c, layer_details=info)
    return [total[:, i] for i in range(m)], meta


def apply_noise(ideal_signals, t_master_len, fs, cfg, gt=None, mic_pos=None, c=343.0):
    """Environmental noise/reverb ONLY. Output stays in Pa units, float -- no mic response,
    no calibration, no quantization (Stage 3)."""
    noi, out = cfg["noise"], cfg["output"]
    n_mics = len(ideal_signals)
    seed = out["seed"]

    if noi["model"] == "simple":
        rng_adc = np.random.default_rng(seed)
        pa_signals = [ideal_signals[i] + noi["noise_floor_pa"] * rng_adc.standard_normal(t_master_len)
                      for i in range(n_mics)]
        return pa_signals, dict(model="simple", noise_floor_pa=noi["noise_floor_pa"])
    if noi["model"] == "realistic":
        pa_signals = []
        for i in range(n_mics):
            dry = ideal_signals[i]
            rir = gp.make_room_ir(noi["rt60"], fs, seed=1000 + seed + i)
            wet = fftconvolve(dry, rir, mode="full")[:len(dry)]
            if np.abs(dry).max() > 0:
                wet *= np.abs(dry).max() / (np.abs(wet).max() + 1e-12)
            noise = gp.make_colored_noise(len(dry), fs, noi["noise_slope"], noi["noise_rms_pa"],
                                          seed=2000 + seed + i)
            pa_signals.append(wet + noise)
        return pa_signals, dict(model="realistic", rt60=noi["rt60"], noise_rms_pa=noi["noise_rms_pa"],
                                noise_slope=noi["noise_slope"])
    if noi["model"] == "layers":
        return apply_layers(ideal_signals, fs, cfg, gt, mic_pos, c)
    raise ValueError(f"Unknown noise.model '{noi['model']}' -- use 'simple', 'realistic' or 'layers'.")


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"not JSON serialisable: {type(o)}")


def main():
    parser = argparse.ArgumentParser(
        description="Stage 2: add environmental noise, echoes and reverberation to a Stage-1 clean signal.")
    parser.add_argument("clean_wav", help="Path to the *_clean.wav from gs_gen_clean_signal.py")
    parser.add_argument("config", nargs="?", default="gen_config.ini",
                        help="Path to .ini config file (default: gen_config.ini)")
    parser.add_argument("-o", "--output", default=None,
                        help="Output basename (overrides output.basename in config); "
                             "produces <output>_noisy.wav and <output>_noisy.json")
    parser.add_argument("--clean-json", default=None,
                        help="Explicit path to the Stage-1 .json (auto-detected by default)")
    parser.add_argument("-layers", "--layers", default=None,
                        help="comma-separated layers (sets model = layers), e.g. colored,wind,tones")
    parser.add_argument("-snr", "--snr", type=float, default=None,
                        help="target SNR in dB (sets level_mode = snr)")
    parser.add_argument("-seed", "--seed", type=int, default=None, help="random seed (overrides output.seed)")
    gp.add_mic_config_argument(parser)
    args = parser.parse_args()

    ideal_signals, fs, clean_meta = load_clean_wav(args.clean_wav, args.clean_json)
    cfg = gp.load_config(args.config, args.mic_config or gp.upstream_mic_config(clean_meta))
    if args.output:
        cfg["output"]["basename"] = args.output
    if args.layers:
        cfg["noise"]["model"] = "layers"
        cfg["noise"]["layers"] = [x.strip().lower() for x in args.layers.split(",") if x.strip()]
    if args.snr is not None:
        cfg["noise"]["level_mode"], cfg["noise"]["snr_db"] = "snr", args.snr
    if args.seed is not None:
        cfg["output"]["seed"] = args.seed

    gt = clean_meta.get("ground_truth") if clean_meta else None
    mic_pos = np.asarray(gt["mic_pos"]) if gt and "mic_pos" in gt else gp.build_mic_positions(cfg)
    c = float(gt["speed_of_sound_mps"]) if gt and "speed_of_sound_mps" in gt \
        else gp.speed_of_sound(cfg["atmosphere"]["temp_c"])
    if len(mic_pos) != len(ideal_signals):
        mic_pos = gp.build_mic_positions(cfg)

    desc = cfg["noise"]["model"] + (f" ({', '.join(cfg['noise']['layers'])})" if cfg["noise"]["model"] == "layers" else "")
    print(f"[Stage 2] Adding noise model '{desc}' to {args.clean_wav} "
          f"({len(ideal_signals)} channels, {fs} Hz)")
    pa_signals, noise_meta = apply_noise(ideal_signals, len(ideal_signals[0]), fs, cfg, gt, mic_pos, c)

    wav_path = f"{cfg['output']['basename']}_noisy.wav"
    json_path = f"{cfg['output']['basename']}_noisy.json"
    write_noisy_wav(wav_path, pa_signals, fs)
    meta = dict(
        generated_at_utc=datetime.now(timezone.utc).isoformat(),
        stage="noisy_signal",
        wav_file=wav_path,
        wav_format="float32 PCM, unnormalized Pa units (mic/ADC not yet applied)",
        source_clean_wav=args.clean_wav,
        sample_rate_hz=fs,
        n_channels=len(pa_signals),
        n_samples=len(pa_signals[0]),
        duration_s=len(pa_signals[0]) / fs,
        config_used=cfg,
        noise_processing=noise_meta,
        clean_signal_ground_truth=gt,
    )
    with open(json_path, "w") as f:
        json.dump(meta, f, indent=2, default=_json_default)

    print(f"[Stage 2] Wrote {wav_path}  ({len(pa_signals)} channels, {fs} Hz, "
          f"{len(pa_signals[0]) / fs:.3f} s, still Pa units)")
    if noise_meta.get("model") == "layers":
        for name, d in noise_meta["layer_details"].items():
            extra = f", {len(d['events'])} impulses" if name == "impulses" else ""
            print(f"[Stage 2]   {name:9s} rms {d.get('rms_pa', 0):.4g} Pa{extra}")
        if noise_meta["snr_achieved_db"] is not None:
            print(f"[Stage 2]   SNR {noise_meta['snr_achieved_db']:.1f} dB "
                  f"({noise_meta['snr_reference']} peak {noise_meta['snr_reference_peak_pa']:.3g} Pa "
                  f"/ noise rms {noise_meta['noise_rms_pa']:.3g} Pa)")
    print(f"[Stage 2] Wrote {json_path}")
    print(f"[Stage 2] Next: python gs_gen_apply_adc.py {wav_path} {args.config}")


if __name__ == "__main__":
    main()
