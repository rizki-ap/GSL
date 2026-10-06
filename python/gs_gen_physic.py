"""
gs_gen_physic.py  --  clean, consolidated physics + config library
===================================================================
Save this file as  gs_gen_physic.py  (it is self-contained: it does NOT import
itself or any other physics module).

Used by the pipeline
  gs_gen_clean_signal.py / gs_gen_clean_signal_claude.py -> <name>_clean.wav/.json
  gs_gen_add_noise.py                                    -> <name>_noisy.wav/.json
  gs_gen_apply_adc.py                                    -> <name>_adc.wav/.json

Contents
  PART 1  the ORIGINAL library (same names, signatures and return values), with
          a few bug fixes -- see "BEHAVIOUR FIXES".
  PART 2  the extended model (wind, temperature-based c, Whitham-type N-wave,
          Hopkinson-Cranz/Kinney-Graham blast, causal spectral absorption,
          arbitrary arrays/trajectories). Used by gs_gen_clean_signal_claude.py.

BEHAVIOUR FIXES (all switched off together by LEGACY_COMPAT = True)
  1. shockwave_arrival / nwave_duration: Mach <= 1 now raises ValueError
     (before: silently returned NaN).
  2. absorption: water-vapour concentration now includes the pressure factor
     (identical at 101.325 kPa, wrong before at other pressures).
  3. N-wave absorption path: the shock travels R = M b / sqrt(M^2-1), not b.
  4. fractional_delay: zero-padded (before: circular, so the end of the buffer
     wrapped to the start).
  5. quantize: full_scale now scales the step, and the top code is 2^(bits-1)-1
     (before: step ignored full_scale; +full scale produced code 2^(bits-1)).
  6. make_colored_noise: DC bin removed (before: kept, and with slope -2.4 the
     "noise" was mostly a DC offset).
  Run an unmodified stage script in legacy mode with  GS_GEN_LEGACY_COMPAT=1.
  Items 1 and 3 (raise on Mach <= 1, shock-ray path) have no legacy switch for
  the extended functions; nwave_waveform(legacy=True) restores item 3.
  Unchanged on purpose: make_tetrahedron (rounded coordinates, identical to the
  copy in gs_det_signal_prepare used by the detector; pass exact=True for the
  regular tetrahedron) and mic_bandpass (zero-phase; pass causal=True for a
  causal filter).

NOT FIXED, only documented (modelling choices, not bugs)
  * nwave_duration / nwave_peak_pressure (legacy model) are empirical
    reference-point laws (b^1/2 and b^-1/2) that were fitted to your hardware
    chain; they are NOT the classical Whitham scaling (b^1/4, b^-3/4).
  * friedlander_* (legacy model) is a reference-point law (1/r, r^1/3); it is not
    a Hopkinson-Cranz fit with an equivalent charge W.
  * make_room_ir: the diffuse tail is ~40 dB stronger than the direct path
    (direct-to-reverberant ratio about -40 dB at rt60 = 1.25 s).
  * make_colored_noise: `rms` is the rms BEFORE the microphone band-pass.

VERIFY BEFORE RELYING ON ABSOLUTE LEVELS
  * WHITHAM_K_P = 0.53 and WHITHAM_K_T = 1.82 are the classical coefficients as
    remembered from Whitham (1952) / Maher (2006).
  * The Kinney-Graham expressions are written from memory (free-air burst).
  * The 9 mm entries have NO fitted muzzle-blast constants (None).
Not meant to be run directly.
"""

import configparser
import os
import sys

import numpy as np
from scipy.fft import next_fast_len
from scipy.signal import butter, sosfilt, sosfiltfilt

# ===========================================================================
# Global switch for the behaviour fixes (see header)
# ===========================================================================
# Environment override so the unmodified stage scripts can be run in legacy mode:
#   GS_GEN_LEGACY_COMPAT=1 python gs_gen_clean_signal.py ...
LEGACY_COMPAT = os.environ.get("GS_GEN_LEGACY_COMPAT", "0") == "1"


def set_legacy_compat(flag):
    """True: reproduce the pre-fix behaviour of items 2-6 above exactly."""
    global LEGACY_COMPAT
    LEGACY_COMPAT = bool(flag)


def _legacy(flag=None):
    return LEGACY_COMPAT if flag is None else bool(flag)


# ===========================================================================
# PART 1 -- original library
# ===========================================================================
# Physical constants -------------------------------------------------------
C = 343.0              # speed of sound (m/s), about 20 C (ISA sea level is 340.3)
C_LEGACY = C
P_ATM = 101325.0       # atmospheric pressure (Pa)

# Bullet library. Original keys: L, dP0_sw, b0_sw, P_REF_MB, R_REF_MB, T_POS_REF.
# Added keys: d (diameter, m) and v0 (nominal muzzle velocity, m/s; informational).
# None = not calibrated (must be supplied in [muzzle_blast] or use kinney_graham).
BULLET_LIBRARY = {
    '7.62_NATO': dict(L=0.028, dP0_sw=7.5, b0_sw=50.0,
                      P_REF_MB=200.0, R_REF_MB=10.0, T_POS_REF=0.003,
                      d=7.82e-3, v0=850.0),
    '5.56_NATO': dict(L=0.023, dP0_sw=7.5, b0_sw=50.0,
                      P_REF_MB=200.0, R_REF_MB=10.0, T_POS_REF=0.003 / 4.0,
                      d=5.56e-3, v0=940.0),
    '9mm_Parabellum': dict(L=0.0155, dP0_sw=None, b0_sw=None,
                           P_REF_MB=None, R_REF_MB=None, T_POS_REF=None,
                           d=9.02e-3, v0=370.0),
    'Glock_17_9mm': dict(L=0.0155, dP0_sw=None, b0_sw=None,      # 124 gr FMJ
                         P_REF_MB=None, R_REF_MB=None, T_POS_REF=None,
                         d=9.02e-3, v0=375.0),
}
BULLET_LIBRARY_CLAUDE = BULLET_LIBRARY     # alias used by gs_gen_clean_signal_claude.py


# Config loading -- shared schema, each stage only reads the sections it needs
DEFAULTS = {
    "geometry":    {"l_array": "0.30"},
    "adc":         {"fs": "100000", "bit_depth": "16",
                    "mic_lo": "40.0", "mic_hi": "20000.0"},
    "calibration": {"real_noise_rms": "0.000595", "real_rt60": "1.32",
                    "real_peak_norm": "0.098", "real_noise_slope": "-2.46"},
    "analog_chain": {"enabled": "false", "mic_sensitivity_mv_per_pa": "22.4",
                     "preamp_gain_db": "40.0", "adc_vref_peak_v": "2.5"},
    "adc_multi_channel": {"type": "simultaneous", "conversion_time_us": "2.0",
                          "gain_mismatch_pct": "0.0", "offset_mismatch_mv": "0.0",
                          "mismatch_seed": "42"},
    "bullet":      {"caliber": "5.56_NATO", "mach": "2.5"},
    "trajectory":  {"range": "200.0", "y_miss": "10.0"},
    "atmosphere":  {"enabled": "true", "temp_c": "20.0", "humidity_pct": "50.0",
                    "pressure_kpa": "101.325"},
    "timing":      {"pre_roll": "0.010", "post_roll": "0.150"},
    "noise":       {"model": "realistic", "noise_floor_pa": "0.005",
                    "rt60": "1.25", "noise_rms_pa": "0.075", "noise_slope": "-2.4"},
    "output":      {"basename": "gunshot_sim", "seed": "0"},
}


def _validate_config(cfg):
    """Fail early, with a clear message, on values that would otherwise give NaN/inf."""
    def need(cond, msg):
        if not cond:
            raise ValueError(f"config: {msg}")
    need(cfg["geometry"]["l_array"] > 0, "geometry.l_array must be > 0")
    a = cfg["adc"]
    need(a["fs"] > 0 and a["bit_depth"] >= 2, "adc.fs must be > 0 and adc.bit_depth >= 2")
    need(0 < a["mic_lo"] < a["mic_hi"], "adc.mic_lo must be > 0 and < adc.mic_hi")
    at = cfg["atmosphere"]
    need(0.0 <= at["humidity_pct"] <= 100.0, "atmosphere.humidity_pct must be in [0, 100]")
    need(at["pressure_kpa"] > 0 and at["temp_c"] > -273.15, "atmosphere pressure/temperature invalid")
    need(cfg["bullet"]["mach"] >= 0, "bullet.mach must be >= 0")
    need(cfg["timing"]["pre_roll"] >= 0 and cfg["timing"]["post_roll"] >= 0, "timing roll values must be >= 0")
    need(cfg["noise"]["rt60"] > 0, "noise.rt60 must be > 0")
    need(cfg["adc_multi_channel"]["type"] in ("simultaneous", "multiplexed"),
         "adc_multi_channel.type must be 'simultaneous' or 'multiplexed'")


def load_config(path):
    """Load an .ini config file, falling back to DEFAULTS for anything missing
    (including a missing file entirely). Returns a nested dict of RESOLVED
    values (native types) -- used for computation AND echoed back verbatim into
    the output .json for reproducibility.

    Per-channel hardware: optional [analog_chain_ch0]..[analog_chain_ch3]
    sections override mic_sensitivity_mv_per_pa / preamp_gain_db for that
    channel; anything not overridden falls back to [analog_chain]. Result is in
    cfg["analog_chain_per_channel"], a list of 4 dicts (the noise/ADC stages
    assume 4 channels)."""
    cp = configparser.ConfigParser()
    cp.read_dict(DEFAULTS)
    if path is not None:
        found = cp.read(path)
        if not found:
            print(f"WARNING: config file '{path}' not found -- using built-in defaults.",
                  file=sys.stderr)

    n_channels = 4
    nominal = dict(mic_sensitivity_mv_per_pa=cp.getfloat("analog_chain", "mic_sensitivity_mv_per_pa"),
                   preamp_gain_db=cp.getfloat("analog_chain", "preamp_gain_db"))
    per_channel = []
    for i in range(n_channels):
        section = f"analog_chain_ch{i}"
        ch = dict(nominal)
        if cp.has_section(section):
            if cp.has_option(section, "mic_sensitivity_mv_per_pa"):
                ch["mic_sensitivity_mv_per_pa"] = cp.getfloat(section, "mic_sensitivity_mv_per_pa")
            if cp.has_option(section, "preamp_gain_db"):
                ch["preamp_gain_db"] = cp.getfloat(section, "preamp_gain_db")
        per_channel.append(ch)

    cfg = {
        "geometry": dict(l_array=cp.getfloat("geometry", "l_array")),
        "adc": dict(fs=cp.getint("adc", "fs"),
                    bit_depth=cp.getint("adc", "bit_depth"),
                    mic_lo=cp.getfloat("adc", "mic_lo"),
                    mic_hi=cp.getfloat("adc", "mic_hi")),
        "calibration": dict(real_noise_rms=cp.getfloat("calibration", "real_noise_rms"),
                            real_rt60=cp.getfloat("calibration", "real_rt60"),
                            real_peak_norm=cp.getfloat("calibration", "real_peak_norm"),
                            real_noise_slope=cp.getfloat("calibration", "real_noise_slope")),
        "analog_chain": dict(enabled=cp.getboolean("analog_chain", "enabled"),
                             mic_sensitivity_mv_per_pa=nominal["mic_sensitivity_mv_per_pa"],
                             preamp_gain_db=nominal["preamp_gain_db"],
                             adc_vref_peak_v=cp.getfloat("analog_chain", "adc_vref_peak_v")),
        "analog_chain_per_channel": per_channel,
        "adc_multi_channel": dict(
            type=cp.get("adc_multi_channel", "type"),
            conversion_time_us=cp.getfloat("adc_multi_channel", "conversion_time_us"),
            gain_mismatch_pct=cp.getfloat("adc_multi_channel", "gain_mismatch_pct"),
            offset_mismatch_mv=cp.getfloat("adc_multi_channel", "offset_mismatch_mv"),
            mismatch_seed=cp.getint("adc_multi_channel", "mismatch_seed")),
        "bullet": dict(caliber=cp.get("bullet", "caliber"),
                       mach=cp.getfloat("bullet", "mach")),
        "trajectory": dict(range=cp.getfloat("trajectory", "range"),
                           y_miss=cp.getfloat("trajectory", "y_miss")),
        "atmosphere": dict(enabled=cp.getboolean("atmosphere", "enabled"),
                           temp_c=cp.getfloat("atmosphere", "temp_c"),
                           humidity_pct=cp.getfloat("atmosphere", "humidity_pct"),
                           pressure_kpa=cp.getfloat("atmosphere", "pressure_kpa")),
        "timing": dict(pre_roll=cp.getfloat("timing", "pre_roll"),
                       post_roll=cp.getfloat("timing", "post_roll")),
        "noise": dict(model=cp.get("noise", "model"),
                      noise_floor_pa=cp.getfloat("noise", "noise_floor_pa"),
                      rt60=cp.getfloat("noise", "rt60"),
                      noise_rms_pa=cp.getfloat("noise", "noise_rms_pa"),
                      noise_slope=cp.getfloat("noise", "noise_slope")),
        "output": dict(basename=cp.get("output", "basename"),
                       seed=cp.getint("output", "seed")),
    }
    _validate_config(cfg)
    return cfg


# Section 1 -- Sensor geometry ----------------------------------------------
def make_tetrahedron(edge_length, exact=False):
    """Four microphone positions, edge ~ edge_length (m).

    exact=False (default) uses the rounded coordinates that gs_det_signal_prepare
    (the detector) also uses, so generator and detector agree; the six edges then
    differ by up to ~0.15 mm at 0.30 m. exact=True gives a regular tetrahedron."""
    if exact:
        raw = np.array([[0.0, 0.0, 1.0],
                        [0.0, 2.0 * np.sqrt(2.0) / 3.0, -1.0 / 3.0],
                        [-np.sqrt(6.0) / 3.0, -np.sqrt(2.0) / 3.0, -1.0 / 3.0],
                        [np.sqrt(6.0) / 3.0, -np.sqrt(2.0) / 3.0, -1.0 / 3.0]], dtype=float)
    else:
        raw = np.array([[0.000, 0.000, 1.000], [0.000, 0.943, -0.333],
                        [-0.816, -0.471, -0.333], [0.816, -0.471, -0.333]], dtype=float)
    scale = edge_length / np.linalg.norm(raw[0] - raw[1])
    return raw * scale


# Section 2 -- Mic sensitivity / ADC (used only by the noise / ADC stages) -----
def mic_bandpass(x, fs, lo, hi, order=4, causal=False):
    """Band-pass microphone response. causal=False (default, original): zero-phase
    sosfiltfilt (squared magnitude, no delay, pre-ringing). causal=True: sosfilt."""
    hi_eff = min(hi, fs / 2 * 0.99)
    if not 0 < lo < hi_eff:
        raise ValueError(f"mic_bandpass needs 0 < lo < hi (got lo={lo}, hi={hi_eff} at fs={fs})")
    sos = butter(order, [lo, hi_eff], btype='band', fs=fs, output='sos')
    return sosfilt(sos, x) if causal else sosfiltfilt(sos, x)


def quantize(x, bits, full_scale=1.0, legacy=None):
    """Round to a signed `bits`-bit grid spanning +-full_scale; clip to the code range.

    Fixed: the step is full_scale/2^(bits-1) and the largest code is 2^(bits-1)-1.
    legacy=True (or LEGACY_COMPAT): the old behaviour (step 1/2^(bits-1) whatever
    full_scale is, and +full_scale maps to code 2^(bits-1))."""
    levels = 2 ** (bits - 1)
    x = np.asarray(x, dtype=float)
    if _legacy(legacy):
        return np.round(np.clip(x, -full_scale, full_scale) * levels) / levels
    code = np.clip(np.round(x / full_scale * levels), -levels, levels - 1)
    return code / levels * full_scale


# Analog gain chain (mic sensitivity -> preamp gain -> ADC reference):
#   v_mic = S_mic . p ; v_out = 10^(G_dB/20) . v_mic ; code = v_out / V_ref_peak
#   => PA_TO_NORM = S_mic . 10^(G_dB/20) / V_ref_peak
def dbv_per_pa_to_v_per_pa(dbv_per_pa):
    """Convert a mic sensitivity in dBV/Pa (datasheet convention, e.g. -38 dBV/Pa)
    to linear V/Pa."""
    return 10 ** (dbv_per_pa / 20.0)


def gain_chain_scale(mic_sensitivity_v_per_pa, preamp_gain_db, adc_vref_peak_v):
    """PA_TO_NORM scale factor from explicit hardware specs: mic sensitivity (V/Pa),
    preamp gain (dB), and the ADC peak input voltage that maps to code +-1.0."""
    return mic_sensitivity_v_per_pa * 10 ** (preamp_gain_db / 20.0) / adc_vref_peak_v


# Shared-ADC, multi-channel effects ------------------------------------------
#   * channel gain/offset mismatch: small (~0.1-0.5 %).
#   * timing skew if the ADC is MULTIPLEXED: channel i is delayed by
#     i * conversion_time, which directly corrupts the microsecond-scale timing
#     a TDOA array measures.
def fractional_delay(x, delay_samples, zero_pad=None):
    """Delay a signal by a (possibly non-integer) number of samples with an
    FFT-domain linear phase shift (exact for band-limited signals).

    Fixed: the signal is zero-padded, so nothing wraps from the end of the buffer
    to the start (zero_pad=False / LEGACY_COMPAT: the old circular behaviour).
    Samples shifted past the end (or start) are discarded, like a real delay line."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    if _legacy(None if zero_pad is None else (not zero_pad)):
        X = np.fft.rfft(x)
        freqs = np.fft.rfftfreq(n)
        return np.fft.irfft(X * np.exp(-2j * np.pi * freqs * delay_samples), n=n)
    pad = int(np.ceil(abs(delay_samples))) + 256
    nfft = next_fast_len(n + 2 * pad)
    xp = np.zeros(nfft)
    xp[pad:pad + n] = x
    freqs = np.fft.rfftfreq(nfft)
    y = np.fft.irfft(np.fft.rfft(xp) * np.exp(-2j * np.pi * freqs * delay_samples), n=nfft)
    return y[pad:pad + n]


# Section 5.1 -- Shockwave (Mach cone) arrival ------------------------------------
def shockwave_arrival(sensor_pos, bullet_origin, v_hat, M, c):
    """Closed-form Mach-cone arrival for a straight, constant-speed bullet (v_b = M c).
    Returns (t_arrive, t_emit, a, b); raises ValueError for M <= 1 or when the
    sensor is outside the region reached by the Mach cone."""
    if not np.isfinite(M) or M <= 1.0:
        raise ValueError(f"shockwave_arrival needs Mach > 1 (got {M})")
    sensor_pos = np.asarray(sensor_pos, dtype=float)
    bullet_origin = np.asarray(bullet_origin, dtype=float)
    v_hat = np.asarray(v_hat, dtype=float)
    beta = np.sqrt(M**2 - 1.0)
    r = sensor_pos - bullet_origin
    a = float(np.dot(r, v_hat))
    b = float(np.linalg.norm(r - a * v_hat))
    if a - b / beta < 0:
        raise ValueError(
            f"Sensor outside Mach cone (a={a:.1f} < b/beta={b/beta:.1f}). "
            "Increase trajectory.range or reduce trajectory.y_miss.")
    return (a + b * beta) / (M * c), (a - b / beta) / (M * c), a, b   # t_arrive, t_emit, a, b


# Section 5.3 -- N-wave, LEGACY reference-point model ------------------------------
# NOTE: empirical laws (duration ~ b^1/2, amplitude ~ b^-1/2) fitted to the project's
# hardware chain. They are not the classical Whitham scaling (see nwave_params).
def nwave_duration(b, M, c, L):
    if not np.isfinite(M) or M <= 1.0:
        raise ValueError(f"nwave_duration needs Mach > 1 (got {M})")
    return (2.0 / (M * c)) * np.sqrt(b * L / np.sqrt(M**2 - 1.0))


def nwave_peak_pressure(b, b_ref, dP_ref, freq_hz=None, distance_m=None, atmosphere=None):
    dP = dP_ref * np.sqrt(b_ref / b)
    if atmosphere is not None and freq_hz is not None and distance_m is not None:
        dP, _ = apply_atmospheric_attenuation(dP, freq_hz, distance_m, atmosphere)
    return dP


def nwave_waveform(t, t_arrive, b, M, c, L, b_ref, dP_ref, atmosphere=None, legacy=None):
    """Linear N-wave (+dP at the leading shock, 0 at T_N/2, -dP at T_N).

    atmosphere: None (no absorption) or an atmosphere-params dict. The peak is
    attenuated at the characteristic frequency 1/T_N over the SHOCK RAY length
    R = M b / sqrt(M^2-1) (fixed; the old code used the perpendicular distance b,
    which under-states the path by up to ~80 % near Mach 1.2)."""
    T_N = nwave_duration(b, M, c, L)
    f_char = 1.0 / T_N
    path = b if _legacy(legacy) else M * b / np.sqrt(M**2 - 1.0)
    dP = nwave_peak_pressure(b, b_ref, dP_ref, freq_hz=f_char, distance_m=path, atmosphere=atmosphere)
    p = np.zeros(len(t))
    mask = (t >= t_arrive) & (t < t_arrive + T_N)
    tau = t[mask] - t_arrive
    p[mask] = dP * (1.0 - 2.0 * tau / T_N)
    return p, T_N, dP


# Section 5.6 -- Muzzle blast, LEGACY reference-point model ---------------------------
# NOTE: dP = P_ref (R_ref/r), t_pos = T_ref (r/R_ref)^(1/3) is a reference-point scaling,
# not a Hopkinson-Cranz fit with an equivalent charge (see blast_params).
def friedlander_params(r, P_REF_MB, R_REF_MB, T_POS_REF, atmosphere=None):
    if r <= 0:
        raise ValueError("friedlander_params needs r > 0")
    dP = P_REF_MB * (R_REF_MB / r)
    t_pos = T_POS_REF * (r / R_REF_MB) ** (1.0 / 3.0)
    if atmosphere is not None:
        dP, _ = apply_atmospheric_attenuation(dP, 1.0 / t_pos, r, atmosphere)
    return dP, t_pos


def friedlander_waveform(t, t_arrive, r, P_REF_MB, R_REF_MB, T_POS_REF, atmosphere=None):
    """The blast propagates spherically over the full range r, so the geometric (1/r)
    and the atmospheric attenuation use the same distance r."""
    dP, t_pos = friedlander_params(r, P_REF_MB, R_REF_MB, T_POS_REF, atmosphere=atmosphere)
    p = np.zeros(len(t))
    mask = t >= t_arrive
    tau = t[mask] - t_arrive
    p[mask] = dP * (1.0 - tau / t_pos) * np.exp(-tau / t_pos)
    return p, t_pos, dP


# Atmospheric absorption (ISO 9613-1) --------------------------------------------------
# Reference values reproduced at 4000 Hz, 20 C, 101.325 kPa: 109.8 dB/km at 10 % RH
# and 23.1 dB/km at 70 % RH (published: ~109 and ~23).
DEFAULT_ATMOSPHERE = dict(temp_c=20.0, humidity_pct=50.0, pressure_kpa=101.325)


def atmospheric_absorption_coefficient(freq_hz, temp_c=20.0, humidity_pct=50.0,
                                       pressure_kpa=101.325, legacy=None):
    """ISO 9613-1 absorption coefficient (dB/m) for a pure tone at freq_hz.

    Fixed: h = h_r (p_sat/p_r) / (p_a/p_r) as in the standard; the old code left out
    the pressure factor (same result only at p_a = 101.325 kPa)."""
    if not 0.0 <= humidity_pct <= 100.0 or pressure_kpa <= 0 or temp_c <= -273.15:
        raise ValueError("invalid atmosphere (humidity 0-100 %, pressure > 0, temperature > -273.15 C)")
    T = temp_c + 273.15       # K
    T0 = 293.15               # reference temperature, 20 C
    T01 = 273.16              # triple-point temperature
    Pr = 101.325              # reference pressure, kPa
    Pa = pressure_kpa
    f = np.asarray(freq_hz, dtype=float)

    psat_over_pr = 10 ** (-6.8346 * (T01 / T) ** 1.261 + 4.6151)
    h = humidity_pct * psat_over_pr                    # molar concentration of water vapour, %
    if not _legacy(legacy):
        h = h * (Pr / Pa)

    frO = (Pa / Pr) * (24 + 4.04e4 * h * (0.02 + h) / (0.391 + h))
    frN = (Pa / Pr) * (T / T0) ** (-0.5) * (9 + 280 * h * np.exp(-4.170 * ((T / T0) ** (-1.0 / 3.0) - 1)))

    term1 = 1.84e-11 * (Pr / Pa) * (T / T0) ** 0.5
    term2 = (T / T0) ** (-2.5) * (
        0.01275 * np.exp(-2239.1 / T) / (frO + f**2 / frO) +
        0.1068 * np.exp(-3352.0 / T) / (frN + f**2 / frN))
    return 8.686 * f**2 * (term1 + term2)              # dB/m


def atmospheric_attenuation_db(freq_hz, distance_m, temp_c=20.0, humidity_pct=50.0, pressure_kpa=101.325):
    return atmospheric_absorption_coefficient(freq_hz, temp_c, humidity_pct, pressure_kpa) * distance_m


def apply_atmospheric_attenuation(dP, freq_hz, distance_m, atmosphere):
    """Reduce a peak-pressure value dP by ISO 9613-1 absorption at freq_hz over distance_m.
    atmosphere: dict with temp_c/humidity_pct/pressure_kpa (missing keys fall back to
    DEFAULT_ATMOSPHERE). Returns (dP_attenuated, attenuation_db)."""
    params = dict(DEFAULT_ATMOSPHERE)
    params.update(atmosphere)
    atten_db = atmospheric_attenuation_db(freq_hz, distance_m, **params)
    return dP * 10 ** (-atten_db / 20.0), float(atten_db)


# Section 6.1 -- Noise / reverberation (used only by gs_gen_add_noise.py) -------------
def make_room_ir(rt60, fs, duration=None, seed=None):
    """Synthetic room impulse response: exponentially decaying Gaussian noise (-60 dB in
    amplitude after rt60) with ir[0] = 1. NOTE: the diffuse tail carries ~9000x the energy
    of the direct sample (direct-to-reverberant ratio about -40 dB at rt60 = 1.25 s)."""
    if rt60 <= 0:
        raise ValueError("rt60 must be > 0")
    rng = np.random.default_rng(seed)
    if duration is None:
        duration = rt60 * 1.2
    n = int(duration * fs)
    t = np.arange(n) / fs
    ir = 10 ** (-3.0 * t / rt60) * rng.standard_normal(n)
    ir[0] = 1.0
    return ir


def make_colored_noise(n, fs, slope, rms, seed=None, remove_dc=None):
    """Noise with power spectrum ~ f^slope, scaled to the given rms (rms is the value
    BEFORE any later band-pass). Fixed: the DC bin is removed (remove_dc=False /
    LEGACY_COMPAT keeps it, as before)."""
    rng = np.random.default_rng(seed)
    X = np.fft.rfft(rng.standard_normal(n))
    freqs = np.fft.rfftfreq(n, d=1 / fs)
    freqs[0] = freqs[1]
    X = X * (freqs ** (slope / 2.0))
    if not _legacy(None if remove_dc is None else (not remove_dc)):
        X[0] = 0.0
    colored = np.fft.irfft(X, n=n)
    return colored / (np.sqrt(np.mean(colored**2)) + 1e-12) * rms


# ===========================================================================
# PART 2 -- extended model (wind, temperature-based c, Whitham N-wave, ...)
# ===========================================================================
WHITHAM_K_P = 0.53      # dP = K_P p0 (M^2-1)^(1/8) d / (L^(1/4) b^(3/4))    [VERIFY]
WHITHAM_K_T = 1.82      # T  = K_T M d (b/L)^(1/4) / (c (M^2-1)^(3/8))       [VERIFY]
B_MIN = 0.05            # m, lower bound on b in the amplitude/duration laws (avoids b -> 0)

nwave_duration_legacy = nwave_duration     # aliases kept for gs_gen_clean_signal_claude.py
load_config_legacy = load_config


def speed_of_sound(temp_c):
    """c = 331.3 sqrt(T/273.15), T in kelvin (m/s)."""
    return 331.3 * np.sqrt((temp_c + 273.15) / 273.15)


def _vec(s, n=3):
    parts = [p for p in str(s).replace(";", ",").split(",") if p.strip() != ""]
    v = [float(p) for p in parts]
    if len(v) != n:
        raise ValueError(f"expected {n} comma-separated numbers, got '{s}'")
    return v


EXT_DEFAULTS = {
    "geometry":     {"array_type": "tetrahedron", "mic_positions": "", "mic_positions_file": ""},
    "bullet":       {"muzzle_velocity_mps": "0", "diameter_mm": "0", "length_mm": "0"},
    "trajectory":   {"mode": "legacy", "shooter_pos": "", "direction": "",
                     "azimuth_deg": "0", "elevation_deg": "0"},
    "propagation":  {"speed_of_sound": "temperature", "wind_mps": "0,0,0",
                     "absorption_mode": "spectral", "absorption_causal": "true"},
    "shockwave":    {"model": "whitham", "mach_threshold": "1.1", "outside_cone": "zero",
                     "amp_scale": "1.0", "dur_scale": "1.0"},
    "muzzle_blast": {"model": "reference", "amp_scale": "1.0", "dur_scale": "1.0",
                     "w_g_tnt": "2.0", "p_ref_pa": "", "r_ref_m": "", "t_pos_ref_s": ""},
}


def load_config_ext(path):
    """Original resolved config (all legacy sections) plus cfg['ext'] with the extended
    options. The legacy loader ignores the extra sections/keys, so the same .ini file also
    works with the noise / ADC stages."""
    cfg = load_config(path)
    cp = configparser.ConfigParser()
    cp.read_dict(EXT_DEFAULTS)
    if path is not None:
        cp.read(path)

    def opt_float(sec, key):
        v = cp.get(sec, key).strip()
        return float(v) if v != "" else None

    arr_type = cp.get("geometry", "array_type").strip().lower()
    mic_positions = None
    if arr_type == "custom":
        inline = cp.get("geometry", "mic_positions").strip()
        fpath = cp.get("geometry", "mic_positions_file").strip()
        if inline:
            mic_positions = [_vec(row) for row in inline.split(";") if row.strip()]
        elif fpath:
            mic_positions = np.loadtxt(fpath, delimiter=",", comments="#", ndmin=2).tolist()
        else:
            raise ValueError("geometry.array_type = custom needs mic_positions or mic_positions_file")

    traj_mode = cp.get("trajectory", "mode").strip().lower()
    shooter_pos = _vec(cp.get("trajectory", "shooter_pos")) if cp.get("trajectory", "shooter_pos").strip() else None
    direction = _vec(cp.get("trajectory", "direction")) if cp.get("trajectory", "direction").strip() else None

    cfg["ext"] = dict(
        geometry=dict(array_type=arr_type, mic_positions=mic_positions),
        bullet=dict(muzzle_velocity_mps=cp.getfloat("bullet", "muzzle_velocity_mps"),
                    diameter_mm=cp.getfloat("bullet", "diameter_mm"),
                    length_mm=cp.getfloat("bullet", "length_mm")),
        trajectory=dict(mode=traj_mode, shooter_pos=shooter_pos, direction=direction,
                        azimuth_deg=cp.getfloat("trajectory", "azimuth_deg"),
                        elevation_deg=cp.getfloat("trajectory", "elevation_deg")),
        propagation=dict(speed_of_sound=cp.get("propagation", "speed_of_sound").strip().lower(),
                         wind_mps=_vec(cp.get("propagation", "wind_mps")),
                         absorption_mode=cp.get("propagation", "absorption_mode").strip().lower(),
                         absorption_causal=cp.getboolean("propagation", "absorption_causal")),
        shockwave=dict(model=cp.get("shockwave", "model").strip().lower(),
                       mach_threshold=cp.getfloat("shockwave", "mach_threshold"),
                       outside_cone=cp.get("shockwave", "outside_cone").strip().lower(),
                       amp_scale=cp.getfloat("shockwave", "amp_scale"),
                       dur_scale=cp.getfloat("shockwave", "dur_scale")),
        muzzle_blast=dict(model=cp.get("muzzle_blast", "model").strip().lower(),
                          amp_scale=cp.getfloat("muzzle_blast", "amp_scale"),
                          dur_scale=cp.getfloat("muzzle_blast", "dur_scale"),
                          w_g_tnt=cp.getfloat("muzzle_blast", "w_g_tnt"),
                          p_ref_pa=opt_float("muzzle_blast", "p_ref_pa"),
                          r_ref_m=opt_float("muzzle_blast", "r_ref_m"),
                          t_pos_ref_s=opt_float("muzzle_blast", "t_pos_ref_s")),
    )
    _validate_ext(cfg["ext"])
    return cfg


load_config_claude = load_config_ext        # name used by gs_gen_clean_signal_claude.py


def _validate_ext(ext):
    def ok(name, v, allowed):
        if v not in allowed:
            raise ValueError(f"{name} = '{v}' not in {allowed}")
    ok("geometry.array_type", ext["geometry"]["array_type"], ("tetrahedron", "custom"))
    ok("trajectory.mode", ext["trajectory"]["mode"], ("legacy", "explicit"))
    ok("propagation.speed_of_sound", ext["propagation"]["speed_of_sound"], ("temperature", "fixed"))
    ok("propagation.absorption_mode", ext["propagation"]["absorption_mode"],
       ("spectral", "characteristic", "off"))
    ok("shockwave.model", ext["shockwave"]["model"], ("whitham", "legacy"))
    ok("shockwave.outside_cone", ext["shockwave"]["outside_cone"], ("zero", "error"))
    ok("muzzle_blast.model", ext["muzzle_blast"]["model"], ("reference", "kinney_graham"))


# Scene construction ---------------------------------------------------------------------
def build_mic_positions(cfg):
    g = cfg["ext"]["geometry"]
    if g["array_type"] == "custom":
        pos = np.asarray(g["mic_positions"], dtype=float)
        if pos.ndim != 2 or pos.shape[1] != 3 or len(pos) < 2:
            raise ValueError("mic positions must be an (N, 3) array with N >= 2")
        return pos
    return make_tetrahedron(cfg["geometry"]["l_array"])


def build_trajectory(cfg):
    """Returns (bullet_origin = shooter position, unit direction)."""
    t = cfg["ext"]["trajectory"]
    if t["mode"] == "legacy":
        tr = cfg["trajectory"]
        return np.array([-tr["range"], tr["y_miss"], 0.0]), np.array([1.0, 0.0, 0.0])
    if t["shooter_pos"] is None:
        raise ValueError("trajectory.mode = explicit needs shooter_pos")
    origin = np.asarray(t["shooter_pos"], dtype=float)
    if t["direction"] is not None:
        v = np.asarray(t["direction"], dtype=float)
    else:
        az, el = np.radians(t["azimuth_deg"]), np.radians(t["elevation_deg"])
        v = np.array([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)])
    n = np.linalg.norm(v)
    if n == 0:
        raise ValueError("trajectory direction has zero length")
    return origin, v / n


def resolve_bullet(cfg, c):
    """Library entry (copy) with config overrides, and the Mach number."""
    name = cfg["bullet"]["caliber"]
    if name not in BULLET_LIBRARY:
        raise ValueError(f"Unknown bullet.caliber '{name}' -- valid: {list(BULLET_LIBRARY)}")
    bl = dict(BULLET_LIBRARY[name])
    eb = cfg["ext"]["bullet"]
    if eb["diameter_mm"] > 0:
        bl["d"] = eb["diameter_mm"] * 1e-3
    if eb["length_mm"] > 0:
        bl["L"] = eb["length_mm"] * 1e-3
    if eb["muzzle_velocity_mps"] > 0:
        M = eb["muzzle_velocity_mps"] / c
    elif cfg["bullet"]["mach"] > 0:
        M = cfg["bullet"]["mach"]
    else:
        M = bl["v0"] / c
    return bl, float(M)


# Shockwave geometry (closed form) ---------------------------------------------------------
def shock_geometry(sensor_pos, bullet_origin, v_hat, M, c, wind=None):
    """Mach-cone arrival for a straight, constant-speed bullet (v_b = M c).

      q = x_i - x_s,  a = q.u,  b = |q - a u|,  beta = sqrt(M^2-1)
      emission coordinate  s = a - b/beta          (valid if s >= 0)
      acoustic path        R = M b / beta
      arrival delay        tau = (a + b beta) / (M c)       (still air)

    wind (3-vector, m/s) is applied to the acoustic leg only (first-order):
      tau ~= s/(M c) + R / (c + wind . r_hat).

    Returns a dict (keys a, b, s, beta, in_cone, R, R_shock, x_emit, t_emit, t_arrive,
    c_eff). If the sensor is outside the region reached by the Mach cone (s < 0),
    in_cone is False and the time fields are None. Raises ValueError for M <= 1,
    non-finite M, or a non-positive effective sound speed."""
    if not np.isfinite(M) or M <= 1.0:
        raise ValueError(f"shock_geometry needs Mach > 1 (got {M})")
    sensor_pos = np.asarray(sensor_pos, dtype=float)
    bullet_origin = np.asarray(bullet_origin, dtype=float)
    v_hat = np.asarray(v_hat, dtype=float)
    beta = float(np.sqrt(M ** 2 - 1.0))
    q = sensor_pos - bullet_origin
    a = float(np.dot(q, v_hat))
    b = float(np.linalg.norm(q - a * v_hat))
    s = a - b / beta
    rec = dict(a=a, b=b, s=s, beta=beta, in_cone=bool(s >= 0.0),
               R=None, R_shock=None, x_emit=None, t_emit=None, t_arrive=None, c_eff=None)
    if s < 0.0:
        return rec
    R = M * b / beta
    x_emit = bullet_origin + s * v_hat
    t_emit = s / (M * c)
    c_eff = float(c)
    if wind is not None and R > 0:
        c_eff = float(c) + float(np.dot(wind, (sensor_pos - x_emit) / R))
        if c_eff <= 0.0:
            raise ValueError(f"effective sound speed must be positive (got {c_eff:.4g} m/s)")
    rec.update(R=R, R_shock=R, x_emit=x_emit, t_emit=t_emit, c_eff=c_eff,
               t_arrive=t_emit + R / c_eff)
    return rec


def muzzle_arrival(sensor_pos, shooter_pos, c, wind=None):
    q = np.asarray(sensor_pos, dtype=float) - np.asarray(shooter_pos, dtype=float)
    r = float(np.linalg.norm(q))
    c_eff = float(c)
    if wind is not None and r > 0:
        c_eff = float(c) + float(np.dot(wind, q / r))
        if c_eff <= 0.0:
            raise ValueError(f"effective sound speed must be positive (got {c_eff:.4g} m/s)")
    return r / c_eff, r, c_eff


# N-wave amplitude and duration (extended) ------------------------------------------------------
def nwave_params(b, M, c, bullet, model="whitham", amp_scale=1.0, dur_scale=1.0, p_atm=P_ATM):
    """Peak overpressure (Pa) and total duration T_N (s) at perpendicular distance b.

    whitham: classical weak-shock scaling (b^-3/4, b^1/4; depends on d, L, M)
        dP  = K_P p_atm (M^2-1)^(1/8) d / (L^(1/4) b^(3/4))
        T_N = K_T M d (b/L)^(1/4) / (c (M^2-1)^(3/8))
    legacy: the original reference-point model (b^-1/2, b^1/2).
    b is clamped to B_MIN (5 cm) to avoid the b -> 0 singularity."""
    if not np.isfinite(M) or M <= 1.0:
        raise ValueError(f"N-wave model requires M > 1; got M={M:.6g}")
    be = max(float(b), B_MIN)
    k = M * M - 1.0
    if model == "whitham":
        dP = WHITHAM_K_P * p_atm * k ** (1 / 8) * bullet["d"] / (bullet["L"] ** 0.25 * be ** 0.75)
        T = WHITHAM_K_T * M * bullet["d"] * (be / bullet["L"]) ** 0.25 / (c * k ** (3 / 8))
    elif model == "legacy":
        if bullet.get("dP0_sw") is None:
            raise ValueError("legacy N-wave constants are not defined for this bullet")
        dP = bullet["dP0_sw"] * np.sqrt(bullet["b0_sw"] / be)
        T = nwave_duration(be, M, c, bullet["L"])
    else:
        raise ValueError(model)
    return float(dP * amp_scale), float(T * dur_scale)


def nwave_waveform_claude(t, t_arrive, dP, T_N):
    """Linear N-wave: +dP at the leading shock (t_arrive), 0 at T_N/2, -dP at T_N."""
    p = np.zeros(len(t))
    mask = (t >= t_arrive) & (t < t_arrive + T_N)
    p[mask] = dP * (1.0 - 2.0 * (t[mask] - t_arrive) / T_N)
    return p


# Muzzle blast (extended) ---------------------------------------------------------------------------
def kinney_graham(Z):
    """Free-air burst fits (Kinney & Graham 1985) [VERIFY].
    Z in m/kg^(1/3). Returns (p_so / p_atm, t_d / W^(1/3) in s/kg^(1/3))."""
    Z = float(Z)
    ps = 808.0 * (1.0 + (Z / 4.5) ** 2) / np.sqrt(
        (1.0 + (Z / 0.048) ** 2) * (1.0 + (Z / 0.32) ** 2) * (1.0 + (Z / 1.35) ** 2))
    td_ms = 980.0 * (1.0 + (Z / 0.54) ** 10) / (
        (1.0 + (Z / 0.02) ** 3) * (1.0 + (Z / 0.74) ** 6) * np.sqrt(1.0 + (Z / 6.9) ** 2))
    return float(ps), float(td_ms * 1e-3)


def blast_params(r, bullet, model="reference", amp_scale=1.0, dur_scale=1.0,
                 p_atm=P_ATM, w_kg=None, overrides=None):
    """Peak overpressure (Pa) and positive-phase duration (s) at muzzle range r."""
    if r <= 0:
        raise ValueError("muzzle range must be > 0")
    if model == "reference":
        ov = overrides or {}
        P_REF = ov.get("p_ref_pa") if ov.get("p_ref_pa") is not None else bullet["P_REF_MB"]
        R_REF = ov.get("r_ref_m") if ov.get("r_ref_m") is not None else bullet["R_REF_MB"]
        T_REF = ov.get("t_pos_ref_s") if ov.get("t_pos_ref_s") is not None else bullet["T_POS_REF"]
        if P_REF is None or R_REF is None or T_REF is None:
            raise ValueError("muzzle-blast reference constants are missing for this bullet: set "
                             "[muzzle_blast] p_ref_pa, r_ref_m, t_pos_ref_s, or use model = kinney_graham")
        dP = P_REF * (R_REF / r)
        t_pos = T_REF * (r / R_REF) ** (1.0 / 3.0)
    elif model == "kinney_graham":
        if w_kg is None or w_kg <= 0:
            raise ValueError("kinney_graham needs a positive charge W (muzzle_blast.w_g_tnt)")
        Z = r / w_kg ** (1.0 / 3.0)
        ps, td = kinney_graham(Z)
        dP, t_pos = ps * p_atm, td * w_kg ** (1.0 / 3.0)
    else:
        raise ValueError(model)
    return float(dP * amp_scale), float(t_pos * dur_scale)


def friedlander_waveform_claude(t, t_arrive, dP, t_pos):
    p = np.zeros(len(t))
    mask = t >= t_arrive
    tau = t[mask] - t_arrive
    p[mask] = dP * (1.0 - tau / t_pos) * np.exp(-tau / t_pos)
    return p


# Atmospheric absorption (extended) -----------------------------------------------------------------
def characteristic_attenuation(dP, T_char, distance_m, atmosphere):
    """Original behaviour: one attenuation value at f = 1/T_char applied to the peak."""
    return apply_atmospheric_attenuation(dP, 1.0 / T_char, distance_m, atmosphere)


def spectral_absorption(x, fs, distance_m, atmosphere, causal=True):
    """Apply ISO 9613-1 absorption |H(f)| = 10^(-alpha(f) d / 20) to one event.

    causal = True builds the minimum-phase filter with that magnitude (folded cepstrum),
    so nothing appears before the geometric arrival time; False applies the magnitude only
    (zero phase, small pre-ringing)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    nfft = 1 << int(np.ceil(np.log2(max(2 * n, 16))))
    f = np.fft.rfftfreq(nfft, 1.0 / fs)
    params = dict(DEFAULT_ATMOSPHERE)
    params.update(atmosphere)
    mag = 10.0 ** (-atmospheric_attenuation_db(f, distance_m, **params) / 20.0)
    if causal:
        logm = np.log(np.maximum(mag, 1e-30))
        full = np.concatenate([logm, logm[-2:0:-1]])
        cep = np.fft.ifft(full).real
        fold = np.zeros(nfft)
        fold[0] = cep[0]
        fold[1:nfft // 2] = 2.0 * cep[1:nfft // 2]
        fold[nfft // 2] = cep[nfft // 2]
        H = np.exp(np.fft.fft(fold))[: nfft // 2 + 1]
    else:
        H = mag
    return np.fft.irfft(np.fft.rfft(x, nfft) * H, nfft)[:n]
