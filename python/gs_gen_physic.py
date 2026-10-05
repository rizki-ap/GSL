"""
gs_gen_physic_claude.py
=======================
GPT-corrected physics library for the gunshot signal generator (a superset of
gs_gen_physic.py -- the original file is NOT modified and must sit in the
same folder; unchanged helpers are re-exported from it).

What this adds relative to gs_gen_physic.py
-------------------------------------------
1. Shock geometry returned as a full record (emission point, ray length,
   Mach-cone validity) instead of raising when a sensor is outside the cone.
   Optional first-order wind on the acoustic leg.
2. Speed of sound from temperature (331.3*sqrt(T/273.15)) instead of the
   fixed 343 m/s (switchable back: propagation.speed_of_sound = fixed).
3. N-wave amplitude / duration from the classical Whitham scaling
   (dP ~ b^-3/4, T_N ~ b^1/4, both depending on bullet diameter d, bullet
   length L and Mach number M) with calibration scale factors. The original
   reference-point model is kept as model = legacy.
4. Muzzle blast: original reference-point scaling (model = reference) or
   Hopkinson-Cranz scaled-distance with the Kinney-Graham free-air fits
   (model = kinney_graham, needs the equivalent charge W).
5. Frequency-dependent ISO 9613-1 absorption applied as a causal
   (minimum-phase) filter per event and path length (absorption_mode =
   spectral). The original single-frequency scalar attenuation is kept as
   absorption_mode = characteristic.
6. Shock activation threshold M_th (default 1.1).
7. Arbitrary sensor arrays (tetrahedron, inline list, or CSV file) and an
   explicit shooter position + bullet direction.
8. Ground-truth TDOA matrices for all sensor pairs and ordering checks.

VERIFY BEFORE RELYING ON ABSOLUTE LEVELS
----------------------------------------
* WHITHAM_K_P = 0.53 and WHITHAM_K_T = 1.82 are the classical coefficients
  as remembered from Whitham (1952) / Maher (2006). Check them against the
  sources and then fit nwave amp_scale / dur_scale to your recordings.
* The Kinney-Graham expressions are written from memory (Kinney & Graham,
  Explosive Shocks in Air, 1985). Check them before use. They describe a
  free-air spherical burst; a ground-level burst is usually modelled with a
  larger effective charge. For a gun, W is only an empirical surrogate.
* The 9 mm entry has NO fitted muzzle-blast constants (None). Supply them in
  [muzzle_blast] (p_ref_pa, r_ref_m, t_pos_ref_s) or use kinney_graham.
"""

import configparser
import sys

import numpy as np

import gs_gen_physic as _gp
from gs_gen_physic import (                                   # re-exported, unchanged
    C as C_LEGACY, P_ATM, make_tetrahedron, fractional_delay,
    atmospheric_absorption_coefficient, atmospheric_attenuation_db,
    apply_atmospheric_attenuation, DEFAULT_ATMOSPHERE,
    make_room_ir, make_colored_noise, mic_bandpass, quantize,
    nwave_duration as nwave_duration_legacy,
)

load_config_legacy = _gp.load_config

__version__ = "1.1-gpt-final"

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
WHITHAM_K_P = 0.53      # dP = K_P p0 (M^2-1)^(1/8) d / (L^(1/4) b^(3/4))   [VERIFY]
WHITHAM_K_T = 1.82      # T  = K_T M d (b/L)^(1/4) / (c (M^2-1)^(3/8))      [VERIFY]
B_MIN = 0.05            # m, lower bound on b in amplitude/duration laws (avoids b->0 singularity)

BULLET_LIBRARY_CLAUDE = {
    # d: bullet diameter (m); L: bullet length (m); v0: nominal muzzle velocity (m/s, informational)
    # dP0_sw/b0_sw/P_REF_MB/R_REF_MB/T_POS_REF: original reference-point constants (legacy models)
    '7.62_NATO': dict(d=7.82e-3, L=0.028, v0=850.0, dP0_sw=7.5, b0_sw=50.0,
                      P_REF_MB=200.0, R_REF_MB=10.0, T_POS_REF=0.003),
    '5.56_NATO': dict(d=5.56e-3, L=0.023, v0=940.0, dP0_sw=7.5, b0_sw=50.0,
                      P_REF_MB=200.0, R_REF_MB=10.0, T_POS_REF=0.003 / 4.0),
    # 9 mm: NOT calibrated. Blast constants must be supplied (None = missing).
    '9mm_Parabellum': dict(d=9.02e-3, L=0.0155, v0=370.0, dP0_sw=None, b0_sw=None,
                           P_REF_MB=None, R_REF_MB=None, T_POS_REF=None),
}


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
def speed_of_sound(temp_c):
    """c = 331.3 sqrt(T/273.15), T in kelvin (m/s)."""
    return 331.3 * np.sqrt((temp_c + 273.15) / 273.15)


def _vec(s, n=3):
    parts = [p for p in str(s).replace(";", ",").split(",") if p.strip() != ""]
    v = [float(p) for p in parts]
    if len(v) != n:
        raise ValueError(f"expected {n} comma-separated numbers, got '{s}'")
    return v


# ---------------------------------------------------------------------------
# Configuration (legacy sections via the original loader + extensions)
# ---------------------------------------------------------------------------
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


def load_config_claude(path):
    """Original resolved config (all legacy sections, unchanged) plus cfg['ext']
    with the new options. Unknown sections are ignored by the legacy loader, so
    the same .ini file also works with the unmodified noise / ADC stages."""
    cfg = load_config_legacy(path)
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


def _validate_ext(ext):
    ok = lambda name, v, allowed: (v in allowed) or (_ for _ in ()).throw(
        ValueError(f"{name} = '{v}' not in {allowed}"))
    ok("geometry.array_type", ext["geometry"]["array_type"], ("tetrahedron", "custom"))
    ok("trajectory.mode", ext["trajectory"]["mode"], ("legacy", "explicit"))
    ok("propagation.speed_of_sound", ext["propagation"]["speed_of_sound"], ("temperature", "fixed"))
    ok("propagation.absorption_mode", ext["propagation"]["absorption_mode"],
       ("spectral", "characteristic", "off"))
    ok("shockwave.model", ext["shockwave"]["model"], ("whitham", "legacy"))
    ok("shockwave.outside_cone", ext["shockwave"]["outside_cone"], ("zero", "error"))
    ok("muzzle_blast.model", ext["muzzle_blast"]["model"], ("reference", "kinney_graham"))


# ---------------------------------------------------------------------------
# Scene construction
# ---------------------------------------------------------------------------
def build_mic_positions(cfg):
    g = cfg["ext"]["geometry"]
    if g["array_type"] == "custom":
        pos = np.asarray(g["mic_positions"], dtype=float)
        if pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError("mic positions must be an (N, 3) array")
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
    if name not in BULLET_LIBRARY_CLAUDE:
        raise ValueError(f"Unknown bullet.caliber '{name}' -- valid: {list(BULLET_LIBRARY_CLAUDE)}")
    bl = dict(BULLET_LIBRARY_CLAUDE[name])
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


# ---------------------------------------------------------------------------
# Shockwave geometry (closed form) -- see paper Sec. 3.2
# ---------------------------------------------------------------------------

def _validate_supersonic_mach(M: float, mach_threshold: float = 1.1) -> None:
    """Raise a clear error when a ballistic shock model is not applicable."""
    if not np.isfinite(M):
        raise ValueError(f"Mach number must be finite, got {M!r}")
    if M <= 1.0:
        raise ValueError(
            f"Supersonic shock model requires M > 1; got M={M:.6g}."
        )
    if mach_threshold is not None and M < float(mach_threshold):
        raise ValueError(
            f"Mach number M={M:.6g} is below the configured shock-model "
            f"threshold {mach_threshold:.6g}."
        )


def _effective_sound_speed(c: float, wind: np.ndarray | None,
                           direction: np.ndarray) -> float:
    """Return acoustic propagation speed along a ray, including optional wind."""
    if wind is None:
        return float(c)
    c_eff = float(c) + float(np.dot(wind, direction))
    if c_eff <= 0.0:
        raise ValueError(
            f"Effective sound speed must be positive; got {c_eff:.6g} m/s."
        )
    return c_eff


def shock_geometry(
    origin: np.ndarray,
    direction: np.ndarray,
    mic_pos: np.ndarray,
    M: float,
    c: float,
    wind: np.ndarray | None = None,
    mach_threshold: float = 1.1,
) -> dict:
    """
    Return the retarded shock-wave geometry and arrival timing.

    The cone condition is evaluated before any square-root of M^2-1.
    This prevents the previous M<=1 crash and makes the configured Mach
    threshold operational.
    """
    _validate_supersonic_mach(M, mach_threshold)

    origin = np.asarray(origin, dtype=float)
    direction = np.asarray(direction, dtype=float)
    mic_pos = np.asarray(mic_pos, dtype=float)

    dn = np.linalg.norm(direction)
    if dn == 0.0:
        raise ValueError("Trajectory direction must be non-zero.")
    vhat = direction / dn

    r = mic_pos - origin
    a = float(np.dot(r, vhat))
    bvec = r - a * vhat
    b = float(np.linalg.norm(bvec))

    beta = float(np.sqrt(M * M - 1.0))
    s = a - b / beta
    in_cone = bool(s >= 0.0)

    result = {
        "a": a,
        "b": b,
        "beta": beta,
        "s": s,
        "in_cone": in_cone,
        "R_shock": float(M * b / beta),
        "x_emit": None,
        "t_emit": None,
        "t_arrive": None,
        "c_eff": None,
    }

    if not in_cone:
        return result

    x_emit = origin + s * vhat
    t_emit = s / (M * float(c))

    # Shock propagation direction from emission point to microphone.
    ray = mic_pos - x_emit
    ray_norm = float(np.linalg.norm(ray))
    if ray_norm == 0.0:
        # Degenerate but well-defined limiting case.
        rhat = vhat
    else:
        rhat = ray / ray_norm

    c_eff = _effective_sound_speed(float(c), wind, rhat)
    t_arrive = t_emit + ray_norm / c_eff

    result.update(
        x_emit=x_emit,
        t_emit=float(t_emit),
        t_arrive=float(t_arrive),
        c_eff=float(c_eff),
        R_shock=float(ray_norm),
    )
    return result

def muzzle_arrival(sensor_pos, shooter_pos, c, wind=None):
    q = np.asarray(sensor_pos, dtype=float) - shooter_pos
    r = float(np.linalg.norm(q))
    c_eff = c
    if wind is not None and r > 0:
        c_eff = c + float(np.dot(wind, q / r))
    return r / c_eff, r, c_eff


# ---------------------------------------------------------------------------
# N-wave (shockwave) amplitude and duration
# ---------------------------------------------------------------------------
def nwave_params(b, M, c, bullet, model="whitham", amp_scale=1.0, dur_scale=1.0, p_atm=P_ATM, b_min: float = B_MIN):
    """Peak overpressure (Pa) and total duration T_N (s) at perpendicular distance b.

    whitham: classical weak-shock scaling (b^-3/4, b^1/4; depends on d, L, M)
        dP  = K_P p_atm (M^2-1)^(1/8) d / (L^(1/4) b^(3/4))
        T_N = K_T M d (b/L)^(1/4) / (c (M^2-1)^(3/8))
    legacy: the original reference-point model (b^-1/2, b^1/2).
    """
    be = float(b)
    if be < float(b_min):
        raise ValueError(
            f"Perpendicular distance b={be:.6g} m is below b_min="
            f"{b_min:.6g} m. Set b_min=0 to disable this guard."
        )
    if be <= 0.0:
        raise ValueError("Perpendicular distance b must be > 0.")
    if M <= 1.0:
        raise ValueError(f"N-wave model requires M > 1; got M={M:.6g}.")
    k = M * M - 1.0
    if model == "whitham":
        dP = WHITHAM_K_P * p_atm * k ** (1 / 8) * bullet["d"] / (bullet["L"] ** 0.25 * be ** 0.75)
        T = WHITHAM_K_T * M * bullet["d"] * (be / bullet["L"]) ** 0.25 / (c * k ** (3 / 8))
    elif model == "legacy":
        if bullet.get("dP0_sw") is None:
            raise ValueError("legacy N-wave constants are not defined for this bullet")
        dP = bullet["dP0_sw"] * np.sqrt(bullet["b0_sw"] / be)
        T = nwave_duration_legacy(be, M, c, bullet["L"])
    else:
        raise ValueError(model)
    return float(dP * amp_scale), float(T * dur_scale)


def nwave_waveform_claude(t, t_arrive, dP, T_N):
    """Linear N-wave: +dP at the leading shock (t_arrive), 0 at T_N/2, -dP at T_N."""
    p = np.zeros(len(t))
    mask = (t >= t_arrive) & (t < t_arrive + T_N)
    p[mask] = dP * (1.0 - 2.0 * (t[mask] - t_arrive) / T_N)
    return p


# ---------------------------------------------------------------------------
# Muzzle blast amplitude and duration
# ---------------------------------------------------------------------------
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


# ---------------------------------------------------------------------------
# Atmospheric absorption
# ---------------------------------------------------------------------------
def characteristic_attenuation(dP, T_char, distance_m, atmosphere):
    """Original behaviour: one attenuation value at f = 1/T_char applied to the peak."""
    return apply_atmospheric_attenuation(dP, 1.0 / T_char, distance_m, atmosphere)


def spectral_absorption(x, fs, distance_m, atmosphere, causal=True):
    """Apply ISO 9613-1 absorption |H(f)| = 10^(-alpha(f) d / 20) to one event.

    causal = True builds the minimum-phase filter with that magnitude (folded
    cepstrum), so nothing appears before the geometric arrival time; False
    applies the magnitude only (zero phase, small pre-ringing).
    """
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
