#!/usr/bin/env python3
"""
ana_ground_reflection.py
========================
Measures the ground-reflected shockwave in the Cipatat recordings, using the clean twins of
sound/Cipatat_clean (util_gen_cipatat_clean.py) as templates of the direct N-wave.

Per event and channel:
  1. residual = recording - clean direct wave (both high-passed at 300 Hz)
  2. the direct-wave template is cross-correlated with the residual 0.35-4.5 ms after the
     direct arrival (x 8 upsampling); the best match gives the reflection delay, its
     least-squares gain relative to the direct wave and its correlation r
  3. a channel counts as a valid detection if r >= 0.6 and the delay is not at the lower limit

Array height, independent of the trajectory: for a flat ground the reflection delay is
2 (h + z) sin(phi) / c, so the top microphone (z = 0.2412 m) is later than the ring (z = 0):
    h = 0.2412 * d_ring / (d_top - d_ring),   sin(phi) = c d_ring / (2 h)
(phi: grazing angle of the shock front on the ground).

Usage
  python ana_ground_reflection.py                    # ../sound/Cipatat + ../sound/Cipatat_clean
  python ana_ground_reflection.py --csv out.csv
"""
import argparse, csv, json, os
import numpy as np
import soundfile as sf
from scipy.signal import butter, sosfiltfilt, resample_poly

HERE = os.path.dirname(os.path.abspath(__file__))
FS, UP, Z_TOP = 51200, 8, 0.2412
LAG_MIN, LAG_MAX, R_MIN = 0.00035, 0.0045, 0.6


def hpf(x, fc=300.0):
    return sosfiltfilt(butter(4, fc, btype="high", fs=FS, output="sos"), x, axis=0)


def measure(real_wav, clean_wav, clean_json):
    x, _ = sf.read(real_wav)
    c, _ = sf.read(clean_wav)
    m = json.load(open(clean_json, encoding="utf-8"))
    xf, cf = hpf(x), hpf(c)
    c_snd = m["scene"]["speed_of_sound_mps"]
    out = []
    for ev in m["events"]:
        per = []
        for ch, on in enumerate(ev["shockwave"]["t_arrive_s"]):
            i0, L, J = int((on - 0.0003) * FS), int(0.0012 * FS), int(0.005 * FS)
            tpl = resample_poly(cf[i0:i0 + L, ch], UP, 1)
            res = resample_poly(xf[i0:i0 + J + L, ch] - cf[i0:i0 + J + L, ch], UP, 1)
            cc = np.correlate(res, tpl, mode="valid")
            lags = np.arange(len(cc)) / (FS * UP)
            ok = (lags > LAG_MIN) & (lags < LAG_MAX)
            k = int(np.argmax(np.where(ok, cc, -np.inf)))
            seg = res[k:k + len(tpl)]
            g = float(np.dot(seg, tpl) / np.dot(tpl, tpl))
            r = float(np.corrcoef(seg, tpl)[0, 1])
            valid = bool(r >= R_MIN and lags[k] > LAG_MIN + 2 / (FS * UP))
            per.append(dict(delay_ms=float(lags[k] * 1e3), gain=g, r=r, valid=valid))
        ring = [p["delay_ms"] for p in per[:3] if p["valid"]]
        h = phi = None
        if len(ring) >= 2 and per[3]["valid"] and per[3]["delay_ms"] > np.mean(ring) + 0.02:
            d_r, d_t = np.mean(ring) * 1e-3, per[3]["delay_ms"] * 1e-3
            h = Z_TOP * d_r / (d_t - d_r)
            s = c_snd * d_r / (2 * h)
            phi = float(np.degrees(np.arcsin(min(1.0, s))))
        out.append(dict(event=ev["index"], channels=per, array_height_m=h, grazing_angle_deg=phi,
                        gain_mean=float(np.mean([p["gain"] for p in per if p["valid"]])) if any(p["valid"] for p in per) else None))
    return out


def main():
    ap = argparse.ArgumentParser(description="Ground-reflected shockwave in the Cipatat recordings.")
    ap.add_argument("--real", default=os.path.join(HERE, "..", "sound", "Cipatat"))
    ap.add_argument("--clean", default=os.path.join(HERE, "..", "sound", "Cipatat_clean"))
    ap.add_argument("--csv", default=None)
    a = ap.parse_args()
    rows = []
    for sh in ("01", "02", "03", "04", "05", "06", "15", "16", "17"):
        res = measure(os.path.join(a.real, f"shot-{sh}_audio_4ch_cut.wav"),
                      os.path.join(a.clean, f"shot-{sh}_audio_4ch_clean.wav"),
                      os.path.join(a.clean, f"shot-{sh}_audio_4ch_clean.json"))
        for e in res:
            ch = e["channels"]
            row = dict(shot=sh, event=e["event"],
                       **{f"delay_ms_ch{i+1}": round(p["delay_ms"], 3) if p["valid"] else "" for i, p in enumerate(ch)},
                       **{f"gain_ch{i+1}": round(p["gain"], 2) if p["valid"] else "" for i, p in enumerate(ch)},
                       n_valid=sum(p["valid"] for p in ch),
                       array_height_m=round(e["array_height_m"], 2) if e["array_height_m"] else "",
                       grazing_angle_deg=round(e["grazing_angle_deg"], 1) if e["grazing_angle_deg"] else "",
                       gain_mean=round(e["gain_mean"], 2) if e["gain_mean"] is not None else "")
            rows.append(row)
            print(row)
    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
