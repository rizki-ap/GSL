"""
gs_gen_physic.py -- physics, geometry and configuration library of the gunshot signal generator
================================================================================================

Used by
  gs_gen_clean_signal.py  Stage 1  noiseless multichannel pressure signal  -> <name>_clean.wav / .json
  gs_gen_add_noise.py     Stage 2  reverberation and noise                  -> <name>_noisy.wav / .json
  gs_gen_apply_adc.py     Stage 3  microphone, analog chain, ADC            -> <name>.wav / .json
Not meant to be run directly.

What is modelled
  Ballistic shockwave   Mach-cone emission point and arrival time (closed form), N-wave
                        amplitude and duration from Whitham's weak-shock scaling.
  Muzzle blast          Friedlander pulse; amplitude and duration from Hopkinson-Cranz
                        scaled distance with the Kinney-Graham free-air fits.
  Propagation           speed of sound from temperature, first-order uniform wind,
                        ISO 9613-1 absorption applied as a causal frequency-domain filter.
  Geometry              any number of sensors, explicit shooter position and bullet direction.

Assumptions and validity (see the physics audit for the numbers)
  * straight bullet path at constant speed -- the speed loss of a real bullet shifts arrival
    times and TDOAs noticeably beyond a few tens of metres
  * free field: no ground reflection, no muzzle directivity, no refraction or turbulence
  * wind: uniform; applied to the acoustic path only
  * the N-wave result is a far-field result (the wave must be fully formed)
  * the Kinney-Graham fits are published for scaled distance Z <= 40 m/kg^(1/3); beyond that
    the blast is extrapolated (a warning is issued once). W is an empirical surrogate for a
    gun, not a real charge

To calibrate against recordings: [shockwave] amp_scale, dur_scale; [muzzle_blast] w_g_tnt,
amp_scale, dur_scale.
To verify against the sources: WHITHAM_K_P (amplitude coefficient) and the Kinney-Graham
duration expression.
"""

import configparser
import sys
import warnings

import numpy as np
from scipy.fft import next_fast_len
from scipy.signal import butter, sosfiltfilt

# ===========================================================================
# Constants
# ===========================================================================
P_ATM = 101325.0        # reference atmospheric pressure (Pa)

# Whitham weak-shock coefficients for the N-wave of a projectile
#   dP  = WHITHAM_K_P * p_atm * (M^2-1)^(1/8) * d / (L^(1/4) * b^(3/4))
#   T_N = WHITHAM_K_T * M * d * (b/L)^(1/4) / (c * (M^2-1)^(3/8))
WHITHAM_K_P = 0.53
WHITHAM_K_T = 1.82
B_MIN = 0.05            # m, lower bound on b in the N-wave laws (avoids the b -> 0 singularity)
KG_Z_MAX = 40.0         # m/kg^(1/3), upper end of the published Kinney-Graham range

# Projectiles: d = diameter (m), L = length (m), v0 = nominal muzzle velocity (m/s).
# Nominal values; override with [bullet] diameter_mm, length_mm, muzzle_velocity_mps.
BULLET_LIBRARY = {
    "5.56_NATO":      dict(d=5.56e-3, L=0.023,  v0=940.0),
    "7.62_NATO":      dict(d=7.82e-3, L=0.028,  v0=850.0),
    "9mm_Parabellum": dict(d=9.02e-3, L=0.0155, v0=370.0),
}


# ===========================================================================
# Configuration
# ===========================================================================
DEFAULTS = {
    "geometry":    {"l_array": "0.30", "array_type": "tetrahedron",
                    "mic_positions": "", "mic_positions_file": ""},
    "adc":         {"fs": "100000", "bit_depth": "16",
                    "mic_lo": "40.0", "mic_hi": "20000.0"},
    "calibration": {"real_noise_rms": "0.000595", "real_rt60": "1.32",
                    "real_peak_norm": "0.098", "real_noise_slope": "-2.46"},
    "analog_chain": {"enabled": "false", "mic_sensitivity_mv_per_pa": "22.4",
                     "preamp_gain_db": "40.0", "adc_vref_peak_v": "2.5"},
    "adc_multi_channel": {"type": "simultaneous", "conversion_time_us": "2.0",
                          "gain_mismatch_pct": "0.0", "offset_mismatch_mv": "0.0",
                          "mismatch_seed": "42"},
    "bullet":      {"caliber": "5.56_NATO", "mach": "2.5", "muzzle_velocity_mps": "0",
                    "diameter_mm": "0", "length_mm": "0"},
    "trajectory":  {"mode": "range_miss", "range": "200.0", "y_miss": "10.0",
                    "shooter_pos": "", "direction": "", "azimuth_deg": "0", "elevation_deg": "0"},
    "atmosphere":  {"enabled": "true", "temp_c": "20.0", "humidity_pct": "50.0",
                    "pressure_kpa": "101.325"},
    "propagation": {"wind_mps": "0,0,0"},
    "shockwave":   {"mach_threshold": "1.1", "amp_scale": "1.0", "dur_scale": "1.0"},
    "muzzle_blast": {"w_g_tnt": "2.0", "amp_scale": "1.0", "dur_scale": "1.0"},
    "timing":      {"pre_roll": "0.010", "post_roll": "0.150"},
    "noise":       {"model": "realistic", "noise_floor_pa": "0.005",
                    "rt60": "1.25", "noise_rms_pa": "0.075", "noise_slope": "-2.4"},
    "output":      {"basename": "gunshot_sim", "seed": "0"},
}


def _vec(s, n=3):
    """'x, y, z' (or 'x; y; z') -> list of n floats."""
    parts = [p for p in str(s).replace(";", ",").split(",") if p.strip() != ""]
    v = [float(p) for p in parts]
    if len(v) != n:
        raise ValueError(f"expected {n} comma-separated numbers, got '{s}'")
    return v


def _validate_config(cfg):
    """Fail early, with a clear message, on values that would otherwise give NaN or inf."""
    def need(cond, msg):
        if not cond:
            raise ValueError(f"config: {msg}")
    g, a, at = cfg["geometry"], cfg["adc"], cfg["atmosphere"]
    need(g["l_array"] > 0, "geometry.l_array must be > 0")
    need(g["array_type"] in ("tetrahedron", "custom"), "geometry.array_type must be 'tetrahedron' or 'custom'")
    if g["array_type"] == "custom":
        need(g["mic_positions"] is not None and len(g["mic_positions"]) >= 2,
             "geometry.array_type = custom needs mic_positions or mic_positions_file with >= 2 sensors")
    need(a["fs"] > 0 and a["bit_depth"] >= 2, "adc.fs must be > 0 and adc.bit_depth >= 2")
    need(0 < a["mic_lo"] < a["mic_hi"], "adc.mic_lo must be > 0 and < adc.mic_hi")
    need(0.0 <= at["humidity_pct"] <= 100.0, "atmosphere.humidity_pct must be in [0, 100]")
    need(at["pressure_kpa"] > 0 and at["temp_c"] > -273.15, "atmosphere pressure/temperature invalid")
    b = cfg["bullet"]
    need(b["mach"] >= 0 and b["muzzle_velocity_mps"] >= 0, "bullet.mach and muzzle_velocity_mps must be >= 0")
    need(b["diameter_mm"] >= 0 and b["length_mm"] >= 0, "bullet.diameter_mm and length_mm must be >= 0")
    t = cfg["trajectory"]
    need(t["mode"] in ("range_miss", "explicit"), "trajectory.mode must be 'range_miss' or 'explicit'")
    if t["mode"] == "explicit":
        need(t["shooter_pos"] is not None, "trajectory.mode = explicit needs shooter_pos")
    need(cfg["timing"]["pre_roll"] >= 0 and cfg["timing"]["post_roll"] >= 0, "timing roll values must be >= 0")
    need(cfg["noise"]["rt60"] > 0, "noise.rt60 must be > 0")
    need(cfg["adc_multi_channel"]["type"] in ("simultaneous", "multiplexed"),
         "adc_multi_channel.type must be 'simultaneous' or 'multiplexed'")
    sw, mb = cfg["shockwave"], cfg["muzzle_blast"]
    need(sw["mach_threshold"] >= 1.0, "shockwave.mach_threshold must be >= 1")
    need(sw["amp_scale"] > 0 and sw["dur_scale"] > 0, "shockwave amp_scale and dur_scale must be > 0")
    need(mb["w_g_tnt"] > 0 and mb["amp_scale"] > 0 and mb["dur_scale"] > 0,
         "muzzle_blast w_g_tnt, amp_scale and dur_scale must be > 0")


def load_config(path):
    """Load an .ini file, falling back to DEFAULTS for anything missing (including a missing
    file). Returns a nested dict of resolved values (native types, lists instead of arrays),
    used for the computation AND echoed into every output .json for reproducibility.

    Per-channel hardware: optional [analog_chain_ch0]..[analog_chain_ch3] sections override
    mic_sensitivity_mv_per_pa / preamp_gain_db for that channel (anything not overridden uses
    [analog_chain]); result in cfg["analog_chain_per_channel"], a list of 4 dicts (the noise and
    ADC stages assume 4 channels)."""
    cp = configparser.ConfigParser()
    cp.read_dict(DEFAULTS)
    if path is not None:
        if not cp.read(path):
            print(f"WARNING: config file '{path}' not found -- using built-in defaults.", file=sys.stderr)

    nominal = dict(mic_sensitivity_mv_per_pa=cp.getfloat("analog_chain", "mic_sensitivity_mv_per_pa"),
                   preamp_gain_db=cp.getfloat("analog_chain", "preamp_gain_db"))
    per_channel = []
    for i in range(4):
        section, ch = f"analog_chain_ch{i}", dict(nominal)
        if cp.has_section(section):
            for key in ("mic_sensitivity_mv_per_pa", "preamp_gain_db"):
                if cp.has_option(section, key):
                    ch[key] = cp.getfloat(section, key)
        per_channel.append(ch)

    def opt_vec(section, key):
        s = cp.get(section, key).strip()
        return _vec(s) if s else None

    array_type = cp.get("geometry", "array_type").strip().lower()
    mic_positions = None
    if array_type == "custom":
        inline = cp.get("geometry", "mic_positions").strip()
        fpath = cp.get("geometry", "mic_positions_file").strip()
        if inline:
            mic_positions = [_vec(row) for row in inline.split(";") if row.strip()]
        elif fpath:
            mic_positions = np.loadtxt(fpath, delimiter=",", comments="#", ndmin=2).tolist()

    cfg = {
        "geometry": dict(l_array=cp.getfloat("geometry", "l_array"), array_type=array_type,
                         mic_positions=mic_positions),
        "adc": dict(fs=cp.getint("adc", "fs"), bit_depth=cp.getint("adc", "bit_depth"),
                    mic_lo=cp.getfloat("adc", "mic_lo"), mic_hi=cp.getfloat("adc", "mic_hi")),
        "calibration": dict(real_noise_rms=cp.getfloat("calibration", "real_noise_rms"),
                            real_rt60=cp.getfloat("calibration", "real_rt60"),
                            real_peak_norm=cp.getfloat("calibration", "real_peak_norm"),
                            real_noise_slope=cp.getfloat("calibration", "real_noise_slope")),
        "analog_chain": dict(enabled=cp.getboolean("analog_chain", "enabled"),
                             mic_sensitivity_mv_per_pa=nominal["mic_sensitivity_mv_per_pa"],
                             preamp_gain_db=nominal["preamp_gain_db"],
                             adc_vref_peak_v=cp.getfloat("analog_chain", "adc_vref_peak_v")),
        "analog_chain_per_channel": per_channel,
        "adc_multi_channel": dict(type=cp.get("adc_multi_channel", "type"),
                                  conversion_time_us=cp.getfloat("adc_multi_channel", "conversion_time_us"),
                                  gain_mismatch_pct=cp.getfloat("adc_multi_channel", "gain_mismatch_pct"),
                                  offset_mismatch_mv=cp.getfloat("adc_multi_channel", "offset_mismatch_mv"),
                                  mismatch_seed=cp.getint("adc_multi_channel", "mismatch_seed")),
        "bullet": dict(caliber=cp.get("bullet", "caliber"), mach=cp.getfloat("bullet", "mach"),
                       muzzle_velocity_mps=cp.getfloat("bullet", "muzzle_velocity_mps"),
                       diameter_mm=cp.getfloat("bullet", "diameter_mm"),
                       length_mm=cp.getfloat("bullet", "length_mm")),
        "trajectory": dict(mode=cp.get("trajectory", "mode").strip().lower(),
                           range=cp.getfloat("trajectory", "range"), y_miss=cp.getfloat("trajectory", "y_miss"),
                           shooter_pos=opt_vec("trajectory", "shooter_pos"),
                           direction=opt_vec("trajectory", "direction"),
                           azimuth_deg=cp.getfloat("trajectory", "azimuth_deg"),
                           elevation_deg=cp.getfloat("trajectory", "elevation_deg")),
        "atmosphere": dict(enabled=cp.getboolean("atmosphere", "enabled"),
                           temp_c=cp.getfloat("atmosphere", "temp_c"),
                           humidity_pct=cp.getfloat("atmosphere", "humidity_pct"),
                           pressure_kpa=cp.getfloat("atmosphere", "pressure_kpa")),
        "propagation": dict(wind_mps=_vec(cp.get("propagation", "wind_mps"))),
        "shockwave": dict(mach_threshold=cp.getfloat("shockwave", "mach_threshold"),
                          amp_scale=cp.getfloat("shockwave", "amp_scale"),
                          dur_scale=cp.getfloat("shockwave", "dur_scale")),
        "muzzle_blast": dict(w_g_tnt=cp.getfloat("muzzle_blast", "w_g_tnt"),
                             amp_scale=cp.getfloat("muzzle_blast", "amp_scale"),
                             dur_scale=cp.getfloat("muzzle_blast", "dur_scale")),
        "timing": dict(pre_roll=cp.getfloat("timing", "pre_roll"), post_roll=cp.getfloat("timing", "post_roll")),
        "noise": dict(model=cp.get("noise", "model"), noise_floor_pa=cp.getfloat("noise", "noise_floor_pa"),
                      rt60=cp.getfloat("noise", "rt60"), noise_rms_pa=cp.getfloat("noise", "noise_rms_pa"),
                      noise_slope=cp.getfloat("noise", "noise_slope")),
        "output": dict(basename=cp.get("output", "basename"), seed=cp.getint("output", "seed")),
    }
    _validate_config(cfg)
    return cfg


# ===========================================================================
# Scene construction
# ===========================================================================
def make_tetrahedron(edge_length):
    """Four microphone positions (m) of a tetrahedral array with edge ~ edge_length.

    The coordinates are the rounded ones also used by the detector (gs_det_signal_prepare), so
    generator and detector see the same array; the six edges differ by up to ~0.15 mm at 0.30 m."""
    raw = np.array([[0.000, 0.000, 1.000], [0.000, 0.943, -0.333],
                    [-0.816, -0.471, -0.333], [0.816, -0.471, -0.333]], dtype=float)
    return raw * (edge_length / np.linalg.norm(raw[0] - raw[1]))


def build_mic_positions(cfg):
    """(N, 3) sensor positions from the config."""
    g = cfg["geometry"]
    if g["array_type"] == "custom":
        pos = np.asarray(g["mic_positions"], dtype=float)
        if pos.ndim != 2 or pos.shape[1] != 3 or len(pos) < 2:
            raise ValueError("mic positions must be an (N, 3) array with N >= 2")
        return pos
    return make_tetrahedron(g["l_array"])


def build_trajectory(cfg):
    """(shooter position = bullet origin, unit direction of the bullet).

    range_miss: shooter at (-range, y_miss, 0), bullet along +x, array centroid at the origin.
    explicit:   shooter_pos plus direction (vector) or azimuth/elevation (ENU, azimuth 0 = +x)."""
    t = cfg["trajectory"]
    if t["mode"] == "range_miss":
        return np.array([-t["range"], t["y_miss"], 0.0]), np.array([1.0, 0.0, 0.0])
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
    """(bullet parameters dict, Mach number). Speed priority: muzzle_velocity_mps, then mach,
    then the library's nominal v0. diameter_mm and length_mm override the library."""
    name = cfg["bullet"]["caliber"]
    if name not in BULLET_LIBRARY:
        raise ValueError(f"Unknown bullet.caliber '{name}' -- valid: {list(BULLET_LIBRARY)}")
    bl = dict(BULLET_LIBRARY[name])
    b = cfg["bullet"]
    if b["diameter_mm"] > 0:
        bl["d"] = b["diameter_mm"] * 1e-3
    if b["length_mm"] > 0:
        bl["L"] = b["length_mm"] * 1e-3
    if b["muzzle_velocity_mps"] > 0:
        M = b["muzzle_velocity_mps"] / c
    elif b["mach"] > 0:
        M = b["mach"]
    else:
        M = bl["v0"] / c
    return bl, float(M)


# ===========================================================================
# Propagation and arrival times
# ===========================================================================
def speed_of_sound(temp_c):
    """c = 331.3 sqrt(T / 273.15), T in kelvin (m/s). Dry air; humid air is ~0.2 % faster."""
    return 331.3 * np.sqrt((temp_c + 273.15) / 273.15)


def shock_geometry(sensor_pos, bullet_origin, v_hat, M, c, wind=None):
    """Mach-cone arrival for a straight, constant-speed bullet (v_b = M c).

      q = x_i - x_s,  a = q.u,  b = |q - a u|,  beta = sqrt(M^2 - 1)
      emission coordinate  s = a - b/beta          (valid if s >= 0)
      acoustic path        R = M b / beta
      arrival delay        tau = (a + b beta) / (M c)         (still air)

    wind (3-vector, m/s) acts on the acoustic leg only (first order):
      tau ~= s/(M c) + R / (c + wind . r_hat).

    Returns a dict (a, b, s, beta, in_cone, R, x_emit, t_emit, t_arrive, c_eff). If the sensor is
    outside the region reached by the Mach cone (s < 0), in_cone is False and the time fields are
    None. Raises ValueError for M <= 1 or a non-positive effective sound speed."""
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
               R=None, x_emit=None, t_emit=None, t_arrive=None, c_eff=None)
    if s < 0.0:
        return rec
    R = M * b / beta
    x_emit = bullet_origin + s * v_hat
    t_emit = s / (M * c)
    c_eff = float(c)
    if wind is not None and R > 0:
        c_eff += float(np.dot(wind, (sensor_pos - x_emit) / R))
        if c_eff <= 0.0:
            raise ValueError(f"effective sound speed must be positive (got {c_eff:.4g} m/s)")
    rec.update(R=R, x_emit=x_emit, t_emit=t_emit, c_eff=c_eff, t_arrive=t_emit + R / c_eff)
    return rec


def muzzle_arrival(sensor_pos, shooter_pos, c, wind=None):
    """(arrival delay, range, effective sound speed) of the muzzle blast at a sensor."""
    q = np.asarray(sensor_pos, dtype=float) - np.asarray(shooter_pos, dtype=float)
    r = float(np.linalg.norm(q))
    c_eff = float(c)
    if wind is not None and r > 0:
        c_eff += float(np.dot(wind, q / r))
        if c_eff <= 0.0:
            raise ValueError(f"effective sound speed must be positive (got {c_eff:.4g} m/s)")
    return r / c_eff, r, c_eff


# ===========================================================================
# Shockwave: N-wave amplitude, duration, waveform
# ===========================================================================
def nwave_params(b, M, c, bullet, amp_scale=1.0, dur_scale=1.0, p_atm=P_ATM):
    """Peak overpressure (Pa) and total duration T_N (s) of the N-wave at perpendicular distance b
    from the trajectory (Whitham weak-shock scaling: dP ~ b^-3/4, T_N ~ b^1/4).

      dP  = K_P p_atm (M^2-1)^(1/8) d / (L^(1/4) b^(3/4))
      T_N = K_T M d (b/L)^(1/4) / (c (M^2-1)^(3/8))

    amp_scale and dur_scale are calibration factors. b is limited below by B_MIN."""
    if not np.isfinite(M) or M <= 1.0:
        raise ValueError(f"N-wave model requires M > 1; got M={M:.6g}")
    be = max(float(b), B_MIN)
    k = M * M - 1.0
    dP = WHITHAM_K_P * p_atm * k ** (1 / 8) * bullet["d"] / (bullet["L"] ** 0.25 * be ** 0.75)
    T = WHITHAM_K_T * M * bullet["d"] * (be / bullet["L"]) ** 0.25 / (c * k ** (3 / 8))
    return float(dP * amp_scale), float(T * dur_scale)


def nwave_waveform(t, t_arrive, dP, T_N):
    """Linear N-wave: +dP at the leading shock (t_arrive), 0 at T_N/2, -dP at T_N."""
    p = np.zeros(len(t))
    mask = (t >= t_arrive) & (t < t_arrive + T_N)
    p[mask] = dP * (1.0 - 2.0 * (t[mask] - t_arrive) / T_N)
    return p


# ===========================================================================
# Muzzle blast: Hopkinson-Cranz scaling with the Kinney-Graham fits
# ===========================================================================
_kg_warned = False


def kinney_graham(Z):
    """Free-air burst fits (Kinney & Graham 1985); Z in m/kg^(1/3).
    Returns (peak overpressure / p_atm, positive-phase duration in s per kg^(1/3))."""
    Z = float(Z)
    ps = 808.0 * (1.0 + (Z / 4.5) ** 2) / np.sqrt(
        (1.0 + (Z / 0.048) ** 2) * (1.0 + (Z / 0.32) ** 2) * (1.0 + (Z / 1.35) ** 2))
    td_ms = 980.0 * (1.0 + (Z / 0.54) ** 10) / (
        (1.0 + (Z / 0.02) ** 3) * (1.0 + (Z / 0.74) ** 6) * np.sqrt(1.0 + (Z / 6.9) ** 2))
    return float(ps), float(td_ms * 1e-3)


def blast_params(r, w_kg, amp_scale=1.0, dur_scale=1.0, p_atm=P_ATM):
    """Peak overpressure (Pa) and positive-phase duration (s) of the muzzle blast at range r (m),
    for an equivalent charge w_kg (kg): Z = r / W^(1/3); dP = p_atm f(Z); t_pos = W^(1/3) g(Z).
    Warns once if Z exceeds the published range of the fits (KG_Z_MAX)."""
    global _kg_warned
    if r <= 0 or w_kg <= 0:
        raise ValueError("blast_params needs r > 0 and w_kg > 0")
    Z = r / w_kg ** (1.0 / 3.0)
    if Z > KG_Z_MAX and not _kg_warned:
        _kg_warned = True
        warnings.warn(f"muzzle blast: scaled distance Z = {Z:.0f} m/kg^(1/3) is beyond the published "
                      f"range of the Kinney-Graham fits (Z <= {KG_Z_MAX:.0f}); the blast is extrapolated "
                      "(the duration stops growing). Calibrate with muzzle_blast amp_scale/dur_scale.",
                      RuntimeWarning, stacklevel=2)
    ps, td = kinney_graham(Z)
    return float(ps * p_atm * amp_scale), float(td * w_kg ** (1.0 / 3.0) * dur_scale)


def friedlander_waveform(t, t_arrive, dP, t_pos):
    """Friedlander pulse dP (1 - tau/t_pos) exp(-tau/t_pos) for tau = t - t_arrive >= 0."""
    p = np.zeros(len(t))
    mask = t >= t_arrive
    tau = t[mask] - t_arrive
    p[mask] = dP * (1.0 - tau / t_pos) * np.exp(-tau / t_pos)
    return p


# ===========================================================================
# Atmospheric absorption (ISO 9613-1)
# ===========================================================================
DEFAULT_ATMOSPHERE = dict(temp_c=20.0, humidity_pct=50.0, pressure_kpa=101.325)


def atmospheric_absorption_coefficient(freq_hz, temp_c=20.0, humidity_pct=50.0, pressure_kpa=101.325):
    """ISO 9613-1 absorption coefficient (dB/m) for a pure tone at freq_hz."""
    if not 0.0 <= humidity_pct <= 100.0 or pressure_kpa <= 0 or temp_c <= -273.15:
        raise ValueError("invalid atmosphere (humidity 0-100 %, pressure > 0, temperature > -273.15 C)")
    T = temp_c + 273.15       # K
    T0 = 293.15               # reference temperature, 20 C
    T01 = 273.16              # triple-point temperature
    Pr = 101.325              # reference pressure, kPa
    Pa = pressure_kpa
    f = np.asarray(freq_hz, dtype=float)

    psat_over_pr = 10 ** (-6.8346 * (T01 / T) ** 1.261 + 4.6151)
    h = humidity_pct * psat_over_pr * (Pr / Pa)        # molar concentration of water vapour, %

    frO = (Pa / Pr) * (24 + 4.04e4 * h * (0.02 + h) / (0.391 + h))
    frN = (Pa / Pr) * (T / T0) ** (-0.5) * (9 + 280 * h * np.exp(-4.170 * ((T / T0) ** (-1.0 / 3.0) - 1)))

    term1 = 1.84e-11 * (Pr / Pa) * (T / T0) ** 0.5
    term2 = (T / T0) ** (-2.5) * (
        0.01275 * np.exp(-2239.1 / T) / (frO + f ** 2 / frO) +
        0.1068 * np.exp(-3352.0 / T) / (frN + f ** 2 / frN))
    return 8.686 * f ** 2 * (term1 + term2)            # dB/m


def atmospheric_attenuation_db(freq_hz, distance_m, temp_c=20.0, humidity_pct=50.0, pressure_kpa=101.325):
    return atmospheric_absorption_coefficient(freq_hz, temp_c, humidity_pct, pressure_kpa) * distance_m


def spectral_absorption(x, fs, distance_m, atmosphere):
    """Apply ISO 9613-1 absorption, |H(f)| = 10^(-alpha(f) d / 20), to one event signal.

    The filter is minimum-phase (folded cepstrum), so nothing appears before the geometric
    arrival time. atmosphere: dict with temp_c / humidity_pct / pressure_kpa (missing keys fall
    back to DEFAULT_ATMOSPHERE)."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    nfft = 1 << int(np.ceil(np.log2(max(2 * n, 16))))
    f = np.fft.rfftfreq(nfft, 1.0 / fs)
    params = dict(DEFAULT_ATMOSPHERE)
    params.update(atmosphere)
    mag = 10.0 ** (-atmospheric_attenuation_db(f, distance_m, **params) / 20.0)
    logm = np.log(np.maximum(mag, 1e-30))
    full = np.concatenate([logm, logm[-2:0:-1]])
    cep = np.fft.ifft(full).real
    fold = np.zeros(nfft)
    fold[0] = cep[0]
    fold[1:nfft // 2] = 2.0 * cep[1:nfft // 2]
    fold[nfft // 2] = cep[nfft // 2]
    H = np.exp(np.fft.fft(fold))[: nfft // 2 + 1]
    return np.fft.irfft(np.fft.rfft(x, nfft) * H, nfft)[:n]


# ===========================================================================
# Noise and reverberation (Stage 2)
# ===========================================================================
def make_room_ir(rt60, fs, duration=None, seed=None):
    """Synthetic room impulse response: exponentially decaying Gaussian noise (-60 dB in
    amplitude after rt60) with ir[0] = 1. The diffuse tail carries about 9000 times the energy of
    the direct sample at rt60 = 1.25 s (direct-to-reverberant ratio about -40 dB)."""
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


def make_colored_noise(n, fs, slope, rms, seed=None):
    """Noise with power spectrum ~ f^slope and zero mean, scaled to the given rms (the rms before
    any later band-pass)."""
    rng = np.random.default_rng(seed)
    X = np.fft.rfft(rng.standard_normal(n))
    freqs = np.fft.rfftfreq(n, d=1 / fs)
    freqs[0] = freqs[1]
    X = X * (freqs ** (slope / 2.0))
    X[0] = 0.0
    colored = np.fft.irfft(X, n=n)
    return colored / (np.sqrt(np.mean(colored ** 2)) + 1e-12) * rms


# ===========================================================================
# Microphone, analog chain, ADC (Stage 3)
# ===========================================================================
def mic_bandpass(x, fs, lo, hi, order=4):
    """Band-pass microphone response, zero-phase (sosfiltfilt): the magnitude response is squared
    and there is no delay."""
    hi_eff = min(hi, fs / 2 * 0.99)
    if not 0 < lo < hi_eff:
        raise ValueError(f"mic_bandpass needs 0 < lo < hi (got lo={lo}, hi={hi_eff} at fs={fs})")
    sos = butter(order, [lo, hi_eff], btype="band", fs=fs, output="sos")
    return sosfiltfilt(sos, x)


def quantize(x, bits, full_scale=1.0):
    """Round to a signed `bits`-bit grid spanning +-full_scale and clip to the code range
    (step full_scale / 2^(bits-1), largest code 2^(bits-1) - 1)."""
    levels = 2 ** (bits - 1)
    code = np.clip(np.round(np.asarray(x, dtype=float) / full_scale * levels), -levels, levels - 1)
    return code / levels * full_scale


def dbv_per_pa_to_v_per_pa(dbv_per_pa):
    """Convert a mic sensitivity in dBV/Pa (datasheet convention, e.g. -38 dBV/Pa) to V/Pa."""
    return 10 ** (dbv_per_pa / 20.0)


def gain_chain_scale(mic_sensitivity_v_per_pa, preamp_gain_db, adc_vref_peak_v):
    """Pa -> normalized code scale: mic sensitivity (V/Pa) x preamp gain / ADC peak input voltage
    that maps to code +-1.0."""
    return mic_sensitivity_v_per_pa * 10 ** (preamp_gain_db / 20.0) / adc_vref_peak_v


def fractional_delay(x, delay_samples):
    """Delay a signal by a (possibly non-integer) number of samples with an FFT-domain linear
    phase shift (exact for band-limited signals). The signal is zero-padded, so nothing wraps
    around; samples shifted past either end are discarded."""
    x = np.asarray(x, dtype=float)
    n = len(x)
    pad = int(np.ceil(abs(delay_samples))) + 256
    nfft = next_fast_len(n + 2 * pad)
    xp = np.zeros(nfft)
    xp[pad:pad + n] = x
    freqs = np.fft.rfftfreq(nfft)
    y = np.fft.irfft(np.fft.rfft(xp) * np.exp(-2j * np.pi * freqs * delay_samples), n=nfft)
    return y[pad:pad + n]
