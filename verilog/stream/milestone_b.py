#!/usr/bin/env python3
"""
milestone_b.py -- Milestone B: continuous detection chain (gsl_core)

  python3 milestone_b.py prepare [--python-dir <GSL>/python]
      * coefficient ROMs (gen_stream_coef.py)
      * 10 single-shot scenarios from YOUR generator (gs_gen_clean_signal ->
        gs_gen_add_noise -> gs_gen_apply_adc), simple noise model, 50 ms
        lead-in: 200 m, 350-430 m, 600 m, 50 m (minimum range), 1500 m
        (maximum range) at 0.05 / 0.2 / 0.5 Pa noise, and three with bullet
        drag (5.56 / 7.62 at 1000 m -- 5.56 subsonic, blast only -- and 12.7
        at 1500 m). Needs the drag model in python/gs_gen_physic.py.
      * 4 MISMATCH scenarios at 200 m: the generator's world differs from what
        the detector assumes (mic positions +/-1 mm, shockwave and blast shapes,
        40 C air) while the detector configuration stays untouched
      * an 8-shot burst at 900 rpm (6667 samples apart) built from the s1 clean
        signal and passed through your noise and ADC stages
      Writes mb_<scn>.hex (frames), mb_<scn>_truth.json, mb_cfg.txt, mb_list.txt

  python3 milestone_b.py check
      For every scenario (after the simulation wrote mb_<scn>_*.txt):
        1. bit-true models of fe_cond, mb_decim, sw_det, mb_det, evt_sched:
           x stream, y stream, SW / MB records, jobs  -> must match exactly
        2. engine results for every job vs the Milestone A chain model on the
           model's own windows -> must match exactly (lags, peaks, sidelobes,
           exponents, lag windows)
        3. detection vs ground truth: every shot gives exactly one SW and one
           MB event (one MB only if the bullet is subsonic), Delta-t (MB - SW)
           error <= max(0.5 ms, 0.5 % of the gap)
        4. no replay underflow, no late jobs, no overflow, no engine errors
      Informational: onset errors, TDOA error per pair, MB bearing error,
      and (single shots) the onsets found by python/gs_det_signal_prepare.py.

Configuration order in mb_cfg.txt (43 integers, read by tb_gsl_core.v):
  fe_gain[4] fe_offset[4] fe_dtrim[4] fe_dcblk_en
  sw_thr_mult2 sw_noise_min sw_holdoff sw_dt_min sw_feat_win sw_sym
  mb_thr_mult mb_settle_mult mb_dwell mb_holdoff mb_noise_min mb_settle_timeout
  sw_pre sw_len sw_maxlag mb_pre mb_len mb_maxlag
  sw_band_lo sw_band_hi sw_floor mb_band_lo mb_band_hi mb_floor
  mb_mask_pre mb_mask_post sw_mask_pre sw_mask_post sw_tn_max sw_pk_mult
"""
import argparse
import json
import math
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'gcc_engine'))
sys.path.insert(0, os.path.join(HERE, '..', 'xspec_phat'))
sys.path.insert(0, os.path.join(HERE, '..', 'fft2048'))
sys.path.insert(0, HERE)

from xspec_phat_model import gcc_chain, refine, PAIRS, N          # noqa: E402
from milestone_a import (run_generator, find_python_dir, bearing,  # noqa: E402
                         ang_diff, read_ring)                       # noqa: F401
from milestone_a import detector_array                              # noqa: E402
import gen_stream_coef                                              # noqa: E402

# Separation rule: generator parameters (gen_config.ini, gs_gen_physic, ground-truth JSON)
# are used ONLY to synthesise the recordings and to SCORE results. Detector computations
# (FPGA configuration below, bearing) never read them: the bearing uses det_config.ini
# through the detector's own code.

CFG = dict(
    fe_gain=[16384] * 4, fe_offset=[0] * 4, fe_dtrim=[0] * 4, fe_dcblk_en=1,
    sw_thr_mult2=100, sw_noise_min=16, sw_holdoff=2000, sw_dt_min=1000, sw_feat_win=256,
    sw_sym=32,
    mb_thr_mult=20, mb_settle_mult=9, mb_dwell=100, mb_holdoff=400, mb_noise_min=16,
    mb_settle_timeout=4000,
    sw_pre=300, sw_len=1800, sw_maxlag=100, mb_pre=120, mb_len=1720, mb_maxlag=20,
    sw_band_lo=20, sw_band_hi=512, sw_floor=0, mb_band_lo=10, mb_band_hi=205, mb_floor=0,
    mb_mask_pre=20, mb_mask_post=60, sw_mask_pre=100, sw_mask_post=2000, sw_tn_max=60, sw_pk_mult=9)

BURST_SHOTS, BURST_GAP = 8, 6667            # 900 rpm
SCENARIOS = [
    dict(name='s1', desc='single shot, 200 m, miss 10 m, noise 0.005 Pa',
         patch={'noise': {'model': 'simple'}, 'timing': {'pre_roll': '0.05'}}),
    dict(name='s2', desc='single shot, shooter (350, 250, 20) m, noise 0.05 Pa',
         patch={'noise': {'model': 'simple', 'noise_floor_pa': '0.05'},
                'trajectory': {'mode': 'explicit', 'shooter_pos': '350,250,20',
                               'direction': '-350,-235,-20'},
                'timing': {'pre_roll': '0.05'}}),
    dict(name='s3', desc='single shot, 600 m, miss 25 m, noise 0.02 Pa',
         patch={'noise': {'model': 'simple', 'noise_floor_pa': '0.02'},
                'trajectory': {'range': '600.0', 'y_miss': '25.0'},
                'timing': {'pre_roll': '0.05'}}),
    dict(name='r50', desc='single shot, 50 m (minimum range), miss 10 m, noise 0.005 Pa',
         patch={'noise': {'model': 'simple'}, 'timing': {'pre_roll': '0.05'},
                'trajectory': {'range': '50.0', 'y_miss': '10.0'}}),
] + [
    dict(name=f'L1500_n{n[2:]}', desc=f'single shot, 1500 m (maximum range), miss 25 m, noise {n} Pa',
         patch={'noise': {'model': 'simple', 'noise_floor_pa': n}, 'timing': {'pre_roll': '0.05'},
                'trajectory': {'range': '1500.0', 'y_miss': '25.0'}})
    for n in ('0.05', '0.20', '0.50')
] + [
    dict(name=f'd{cal[:4].replace(".", "")}_{rng}', desc=f'drag: {cal}, {rng} m, miss 25 m, noise 0.05 Pa',
         patch={'noise': {'model': 'simple', 'noise_floor_pa': '0.05'}, 'timing': {'pre_roll': '0.05'},
                'bullet': {'caliber': cal, 'mach': '0', 'drag_model': 'exponential'},
                'trajectory': {'range': f'{rng}.0', 'y_miss': '25.0'}})
    for cal, rng in (('5.56_NATO', 1000), ('7.62_NATO', 1000), ('12.7_NATO', 1500))
] + [
    # MISMATCH scenarios: the generator's world deliberately differs from what the detector
    # assumes (det_config.ini, FPGA configuration). The detector is NOT re-tuned for them.
    dict(name=name, desc=f'mismatch: {desc}, 200 m',
         patch=dict({'noise': {'model': 'simple', 'noise_floor_pa': '0.05'}, 'timing': {'pre_roll': '0.05'}},
                    **extra))
    for name, desc, extra in (
        ('m_tol', 'mic positions off by up to 1 mm (generator only)',
         {'geometry': {'array_type': 'custom', 'mic_positions': '@perturb_1mm'}}),
        ('m_sw', 'shockwave 1.8x longer, 0.5x amplitude', {'shockwave': {'dur_scale': '1.8', 'amp_scale': '0.5'}}),
        ('m_mb', 'muzzle blast 0.4x duration (sharper)', {'muzzle_blast': {'dur_scale': '0.4'}}),
        ('m_hot', 'air at 40 C (c = 355 m/s); detector assumes 343', {'atmosphere': {'temp_c': '40.0'}}))
]


def resolve_patch(sc, pydir):
    """SYNTHESIS side: expand placeholders in a scenario patch (uses the generator)."""
    g = sc['patch'].get('geometry', {})
    if g.get('mic_positions') == '@perturb_1mm':
        sys.path.insert(0, pydir)
        import gs_gen_physic as gp
        nom = np.array(gp.build_mic_positions(gp.load_config(os.path.join(pydir, 'gen_config.ini'))))
        pert = nom + np.random.default_rng(11).uniform(-1e-3, 1e-3, nom.shape)
        g['mic_positions'] = '; '.join(','.join(f'{v:.6f}' for v in row) for row in pert)
    return sc


def cfg_list(c):
    return (c['fe_gain'] + c['fe_offset'] + c['fe_dtrim'] + [c['fe_dcblk_en']] +
            [c[k] for k in ('sw_thr_mult2', 'sw_noise_min', 'sw_holdoff', 'sw_dt_min',
                            'sw_feat_win', 'sw_sym', 'mb_thr_mult', 'mb_settle_mult', 'mb_dwell',
                            'mb_holdoff', 'mb_noise_min', 'mb_settle_timeout', 'sw_pre',
                            'sw_len', 'sw_maxlag', 'mb_pre', 'mb_len', 'mb_maxlag',
                            'sw_band_lo', 'sw_band_hi', 'sw_floor', 'mb_band_lo',
                            'mb_band_hi', 'mb_floor', 'mb_mask_pre', 'mb_mask_post',
                            'sw_mask_pre', 'sw_mask_post', 'sw_tn_max', 'sw_pk_mult')])


# ======================================================================
# prepare
# ======================================================================
def write_stim(fn, x):
    with open(fn, 'w') as f:
        for row in x:
            v = 0
            for c in row:
                v = (v << 16) | (int(c) & 0xFFFF)
            f.write(format(v, '016X') + '\n')


def truth_from_json(g, fs, offsets=(0,)):
    t0 = g['t_master_start_s']
    P = np.array(g['mic_pos'])
    sv = np.array(g['shooter_pos']) - P.mean(axis=0)
    shots = []
    for off in offsets:
        sw = [(p['shockwave']['t_arrive_s'] - t0) * fs + off for p in g['per_mic']] \
            if all(p['shockwave'].get('t_arrive_s') is not None for p in g['per_mic']) else None
        mb = [(p['muzzle_blast']['t_arrive_s'] - t0) * fs + off for p in g['per_mic']]
        shots.append(dict(sw=sw, mb=mb))
    return dict(shots=shots,                      # scoring only: never fed to the detector
                dir=[math.degrees(math.atan2(sv[1], sv[0])),
                     math.degrees(math.asin(sv[2] / np.linalg.norm(sv)))])


def prepare(pydir):
    from scipy.io import wavfile
    gen_stream_coef.main(HERE)
    if pydir is None:
        pydir = find_python_dir()
    if pydir is None or not os.path.isfile(os.path.join(os.path.abspath(pydir), 'gs_gen_clean_signal.py')):
        sys.exit("ERROR: generator (python/gs_gen_clean_signal.py) not found; use --python-dir")
    pydir = os.path.abspath(pydir)
    print(f"Generator: {pydir}")
    work = os.path.join(HERE, 'gen_work')
    names = []
    for sc in SCENARIOS:
        sc = resolve_patch({**sc, 'patch': {k: dict(v) for k, v in sc['patch'].items()}}, pydir)
        wav, tj = run_generator(pydir, work, sc)
        fs, x = wavfile.read(wav)
        g = json.load(open(tj))['ground_truth']
        write_stim(os.path.join(HERE, f"mb_{sc['name']}.hex"), x)
        tr = truth_from_json(g, fs)
        tr.update(name=sc['name'], desc=sc['desc'], nframes=len(x), wav=wav)
        json.dump(tr, open(os.path.join(HERE, f"mb_{sc['name']}_truth.json"), 'w'))
        names.append((sc['name'], len(x)))
        print(f"  {sc['name']}: {len(x)} frames ({len(x) / fs:.3f} s)  {sc['desc']}")

    # ---- 8-shot burst at 900 rpm from the s1 clean signal ----
    fs, clean = wavfile.read(os.path.join(work, 's1_clean.wav'))
    L0 = len(clean)
    burst = np.zeros((L0 + (BURST_SHOTS - 1) * BURST_GAP, clean.shape[1]), dtype=np.float32)
    for k in range(BURST_SHOTS):
        burst[k * BURST_GAP:k * BURST_GAP + L0] += clean
    wavfile.write(os.path.join(work, 'burst_clean.wav'), fs, burst)
    shutil.copy(os.path.join(work, 's1_clean.json'), os.path.join(work, 'burst_clean.json'))
    ini = os.path.join(work, 's1.ini')
    for step in (['gs_gen_add_noise.py', 'burst_clean.wav', ini, '-o', 'burst'],
                 ['gs_gen_apply_adc.py', 'burst_noisy.wav', ini, '-o', 'burst']):
        r = subprocess.run([sys.executable, os.path.join(pydir, step[0])] + step[1:],
                           cwd=work, capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout, r.stderr)
            sys.exit(f"generator step {step[0]} failed for burst")
    fs, x = wavfile.read(os.path.join(work, 'burst.wav'))
    g = json.load(open(os.path.join(work, 's1_clean.json')))['ground_truth']
    write_stim(os.path.join(HERE, 'mb_burst.hex'), x)
    tr = truth_from_json(g, fs, [k * BURST_GAP for k in range(BURST_SHOTS)])
    tr.update(name='burst', desc=f'{BURST_SHOTS} shots at 900 rpm (s1 geometry)',
              nframes=len(x), wav=os.path.join(work, 'burst.wav'))
    json.dump(tr, open(os.path.join(HERE, 'mb_burst_truth.json'), 'w'))
    names.append(('burst', len(x)))
    print(f"  burst: {len(x)} frames ({len(x) / fs:.3f} s)  {tr['desc']}")

    with open(os.path.join(HERE, 'mb_cfg.txt'), 'w') as f:
        f.write('\n'.join(str(v) for v in cfg_list(CFG)) + '\n')
    with open(os.path.join(HERE, 'mb_list.txt'), 'w') as f:
        for n, nf in names:
            f.write(f"{n} {nf}\n")
    print("Wrote mb_cfg.txt, mb_list.txt and stimulus / truth files")


# ======================================================================
# bit-true models
# ======================================================================
def rconv(v, s):
    """Convergent rounding of integer array / int v by s bits."""
    v = np.asarray(v, dtype=object) if isinstance(v, list) else v
    return (v + ((1 << (s - 1)) - 1) + ((v >> s) & 1)) >> s


def sat18(v):
    return np.clip(v, -131072, 131071)


def read_hex18(fn):
    out = []
    for line in open(fn):
        v = int(line.strip(), 16)
        out.append(v - (1 << 18) if v & 0x20000 else v)
    return out


def fe_cond_model(src, cfg):
    """src: [n, 4] int16 frames. Returns x [n, 4] (x_ts = n - 4)."""
    tab = read_hex18(os.path.join(HERE, 'fe_fdelay.hex'))
    n = len(src)
    x = np.zeros((n, 4), dtype=np.int64)
    for c in range(4):
        a = src[:, c].astype(np.int64) + cfg['fe_offset'][c]
        g = sat18(rconv(a * cfg['fe_gain'][c], 14))
        d = cfg['fe_dtrim'][c]
        p = min(abs(d), 128)
        coef = [tab[p * 9 + (8 - k if d < 0 else k)] for k in range(9)]
        acc = np.convolve(g, np.array(coef, dtype=np.int64))[:n]
        y = sat18(rconv(acc, 16))
        if cfg['fe_dcblk_en']:
            z = np.zeros(n, dtype=np.int64)
            dacc, yp = 0, 0
            yl = y.tolist()
            for i in range(n):
                dacc = dacc + ((yl[i] - yp) << 13) - (dacc >> 13)
                yp = yl[i]
                t = (dacc + 4095 + ((dacc >> 13) & 1)) >> 13
                z[i] = max(-131072, min(131071, t))
            x[:, c] = z
        else:
            x[:, c] = y
    return x


def mb_decim_model(x):
    """x: [n, 4] (x index i has x_ts = i - 4). Returns y [m, 4], y_ts, src index."""
    h = np.array(read_hex18(os.path.join(HERE, 'decim63.hex')), dtype=np.int64)
    n = len(x)
    cnt = np.arange(1, n + 1)
    sel = np.where((cnt % 5 == 0) & (cnt >= 63))[0]
    y = np.zeros((len(sel), 4), dtype=np.int64)
    for c in range(4):
        acc = np.convolve(x[:, c], h)[:n]
        y[:, c] = sat18(rconv(acc[sel], 17))
    return y, (sel - 4) - 31, sel


def sw_det_model(x, cfg):
    """Returns records [(t_emit_i, ts, ref, feats)], blank[n], confirmations
    [(i, ts)], qual[n] (qualifying after sample n). Mirrors sw_det.v."""
    m2, nmin = cfg['sw_thr_mult2'], cfg['sw_noise_min']
    fw, sym = cfg['sw_feat_win'], cfg['sw_sym']
    ap = [False] * 4
    nv = False
    thr = [None] * 4
    bsum = [0] * 4
    bcnt, clean, hc, bk, fa, fi = 0, True, 0, 0, False, 0
    qa, qcnt, qpos, qneg, qpo, qno = False, 0, 0, 0, 0, 0
    pp, po, pn, no = [0] * 4, [0] * 4, [0] * 4, [0] * 4
    recs, blank, confs = [], np.zeros(len(x), dtype=bool), []
    qual = np.zeros(len(x), dtype=bool)
    rts, rref = 0, 0
    for i, row in enumerate(x.tolist()):
        rise = [False] * 4
        for c in range(4):
            v = row[c]
            sq = v * v
            abv = nv and sq > thr[c]
            rise[c] = abv and not ap[c]
            ap[c] = abv
            bsum[c] += sq
            if qa and c == rref:
                if v > qpos:
                    qpos, qpo = v, qcnt
                if v < qneg:
                    qneg, qno = v, qcnt
            if fa:
                if v > pp[c]:
                    pp[c], po[c] = v, fi
                if v < pn[c]:
                    pn[c], no[c] = v, fi
        hact, qa_old = hc != 0, qa
        trig = nv and not hact and not qa and any(rise)
        if trig:
            rts, rref = i - 4, rise.index(True)
            qa, qcnt, qpos, qneg, qpo, qno = True, 1, row[rref], row[rref], 0, 0
            fa, fi = True, 1
            pp, pn, po, no = list(row), list(row), [0] * 4, [0] * 4
        else:
            if hc:
                hc -= 1
            if bk:
                bk -= 1
            reject = False
            if qa:
                if qcnt == 99:
                    qa = False
                    P, Q = max(qpos, 0), max(-qneg, 0)
                    mx = max(P, Q)
                    if (mx and 64 * min(P, Q) >= sym * mx and abs(qno - qpo) <= cfg['sw_tn_max']
                            and mx * mx >= cfg['sw_pk_mult'] * thr[rref]):
                        hc, bk = cfg['sw_holdoff'], cfg['sw_dt_min']
                        confs.append((i, rts))
                    else:
                        fa, reject = False, True
                else:
                    qcnt += 1
            if fa and not reject:
                if fi == fw - 1:
                    fa = False
                    recs.append((i, rts, rref, [(pp[c], po[c], pn[c], no[c]) for c in range(4)]))
                else:
                    fi += 1
        if bcnt == 511:
            bcnt = 0
            accept = clean and not trig and not hact and not qa_old
            clean = True
            if accept:
                nms = [b >> 9 for b in bsum]
                nv = True
                thr = [m2 * max(v, nmin) for v in nms]
            bsum = [0] * 4
        else:
            bcnt += 1
            if trig or hact or qa_old:
                clean = False
        blank[i] = bk != 0
        qual[i] = qa
    return recs, blank, confs, qual


def mb_det_model(y, y_ts, src_i, blank, cfg, qual=None, confs=()):
    """Returns records [(t_emit_i, ts, idx, ref, reliable)]. Events triggered
    while sw_det is qualifying are held, then dropped (shockwave confirmed) or
    released (rejected) -- mirrors mb_det.v."""
    conf_at = set(c[0] for c in confs)
    if qual is None:
        qual = np.zeros(len(blank), dtype=bool)
    pend, last_x = None, -1
    tm, sm, nmin = cfg['mb_thr_mult'], cfg['mb_settle_mult'], cfg['mb_noise_min']
    dl = [[0] * 16 for _ in range(4)]
    E = [0] * 4
    ap = [False] * 4
    nv = False
    thrE = thrS = thrI = [None] * 4
    bsum = [0] * 4
    bcnt, clean, hc, qc, age, arm_l = 0, True, 0, 0, 0, False
    recs = []
    yl = y.tolist()
    for m, row in enumerate(yl):
        # resolve a held event on x samples up to this y sample's newest x
        if pend is not None:
            for j in range(last_x + 1, src_i[m] + 1):
                if j in conf_at:
                    pend, hc = None, 0
                    break
                if not qual[j]:
                    recs.append((j,) + pend)
                    pend = None
                    break
        last_x = src_i[m]
        blk = bool(blank[src_i[m]])
        ql = bool(qual[src_i[m]])
        rise, quiet = [False] * 4, [False] * 4
        for c in range(4):
            sq = row[c] * row[c]
            En = E[c] + sq - dl[c][15]
            dl[c] = [sq] + dl[c][:15]
            E[c] = En
            abv = nv and En > thrE[c]
            rise[c] = abv and not ap[c]
            ap[c] = abv
            quiet[c] = nv and En < thrS[c]
            bsum[c] += sq
        hact = hc != 0
        settled = qc >= cfg['mb_dwell']
        tmo = age >= cfg['mb_settle_timeout']
        armed = nv and not blk and not hact and (arm_l or settled or tmo)
        trig = armed and any(rise)
        rel = arm_l or settled
        qc = min(qc + 1, 65535) if all(quiet) else 0
        age = 0 if blk else min(age + 1, 65535)
        arm_l = False if (trig or blk or not nv) else (arm_l or settled)
        if trig:
            hc = cfg['mb_holdoff']
            ref = rise.index(True)
            back = 0
            for a in range(15, -1, -1):
                if dl[ref][a] > thrI[ref]:
                    back = a
                    break
            rec = (int(y_ts[m]) - 5 * back, m - back, ref, int(rel))
            if ql:
                pend = rec
            else:
                recs.append((src_i[m],) + rec)
        elif hc:
            hc -= 1
        if bcnt == 127:
            bcnt = 0
            accept = clean and not trig and not hact and not blk
            clean = True
            if accept:
                nms = [b >> 7 for b in bsum]
                nv = True
                ne = [max(v, nmin) for v in nms]
                thrE = [tm * v * 16 for v in ne]
                thrS = [sm * v * 16 for v in ne]
                thrI = [tm * v for v in ne]
            bsum = [0] * 4
        else:
            bcnt += 1
            if trig or hact or blk:
                clean = False
    return recs


def sched_model(sw_recs, mb_recs, confs, src_i, cfg):
    """Event ids in arrival order (SW before MB in the same frame); jobs per id.
    MB jobs get shockwave masks from the last 4 SW confirmations made up to
    the frame in which the MB window is released (mirrors evt_sched.v)."""
    ev = [(r[0], 0, r) for r in sw_recs] + [(r[0], 1, r) for r in mb_recs]
    ev.sort(key=lambda e: (e[0], e[1]))   # same frame: SW record before MB
    out = []
    for k, (_, typ, r) in enumerate(ev):
        if typ == 0:
            first = r[1] - cfg['sw_pre']
            start = first % 4096
            rel_i = first + cfg['sw_len'] - 1 + 4          # x index of the last window sample
            hist = [m[1] for m in mb_recs if m[0] < rel_i][-4:]
            w_hi = first + cfg['sw_len'] - 1
            masks = []
            for h in hist:
                lo, hi = h - cfg['sw_mask_pre'], h + cfg['sw_mask_post']
                if h > r[1] and lo <= w_hi:
                    masks.append((lo, hi))
            out.append(dict(id=k, type=0, ts=r[1], ref=r[2], feat=r[3], start=start,
                            win_first=first, ts0=first, masks=masks))
        else:
            first = r[2] - cfg['mb_pre']
            start = first % 4096
            ts0 = r[1] - 5 * cfg['mb_pre']
            last_m = first + cfg['mb_len'] - 1
            rel_i = src_i[last_m] if last_m < len(src_i) else 1 << 62
            hist = [c[1] for c in confs if c[0] <= rel_i][-4:]
            w_lo, w_hi = ts0, ts0 + 5 * (cfg['mb_len'] - 1)
            masks = []
            for h in hist:
                lo, hi = h - 5 * cfg['mb_mask_pre'], h + 5 * cfg['mb_mask_post']
                if hi >= w_lo and lo <= w_hi:
                    masks.append((lo, hi))
            out.append(dict(id=k, type=1, ts=r[1], idx=r[2], ref=r[3], rel=r[4], start=start,
                            win_first=first, ts0=ts0, masks=masks))
    return out


# ======================================================================
# check
# ======================================================================
def load_stim(fn, n):
    out = np.zeros((n, 4), dtype=np.int64)
    for i, line in enumerate(open(fn)):
        if i >= n:
            break
        v = int(line.strip(), 16)
        for c in range(4):
            w = (v >> (48 - 16 * c)) & 0xFFFF
            out[i, c] = w - 65536 if w & 0x8000 else w
    return out


def parse_ev(fn):
    S, M, J, X, R, Lw, E = {}, {}, {}, [], {}, {}, None
    for line in open(fn):
        t = line.split()
        if not t:
            continue
        v = list(map(int, t[1:]))
        if t[0] == 'S':
            S[v[0]] = dict(ts=v[1], ref=v[2], feat=[tuple(v[3 + 4 * c:7 + 4 * c]) for c in range(4)])
        elif t[0] == 'M':
            M[v[0]] = dict(ts=v[1], idx=v[2], ref=v[3], rel=v[4])
        elif t[0] == 'J':
            J[v[0]] = dict(type=v[1], start=v[2], men=v[3] if len(v) > 3 else 0)
        elif t[0] == 'X':
            X.append(v[0])
        elif t[0] == 'L':
            Lw.setdefault(v[0], {})[v[1]] = v[2:8]
        elif t[0] == 'R':
            R[v[0]] = dict(type=v[1], err=v[2], ez=(v[3], v[4]), er=v[5:8],
                           lag=v[8:14], pk=v[14:20], side=v[20:26])
        elif t[0] == 'E':
            E = dict(underflow=v[0], late=v[1], ovf=v[2], frames=v[3], cycles=v[4])
    return S, M, J, X, R, Lw, E


def repo_detector(pydir, wav):
    """Run python/gs_det_signal_prepare.py; return (min t_sw, min t_mb) in samples."""
    try:
        work = os.path.join(HERE, 'gen_work')
        r = subprocess.run([sys.executable, os.path.join(pydir, 'gs_det_signal_prepare.py'), wav,
                            os.path.join(pydir, 'det_config.ini'), '-o', 'repo_det'],
                           cwd=work, capture_output=True, text=True, timeout=300)
        d = json.load(open(os.path.join(work, 'repo_det_prepared.json')))
        chans = d.get('channels') or d.get('per_channel') or []
        tsw = [c['t_sw'] for c in chans if c.get('t_sw') is not None]
        tmb = [c['t_mb'] for c in chans if c.get('t_mb') is not None]
        return (min(tsw) * 1e5 if tsw else None, min(tmb) * 1e5 if tmb else None)
    except Exception:
        return None


def check(only=None):
    names = [l.split() for l in open(os.path.join(HERE, 'mb_list.txt')) if l.strip()]
    if only:
        names = [n for n in names if n[0] in only]
    pydir = find_python_dir()
    all_ok = True
    summary = []
    for name, nfr in names:
        nfr = int(nfr)
        tr = json.load(open(os.path.join(HERE, f'mb_{name}_truth.json')))
        print(f"\n=== {name}: {tr['desc']}  ({nfr} frames, {nfr / 1e5:.3f} s)")
        src = load_stim(os.path.join(HERE, f'mb_{name}.hex'), nfr)
        ok = True

        # ---- 1. streams ----
        x = fe_cond_model(src, CFG)
        rx = np.loadtxt(os.path.join(HERE, f'mb_{name}_x.txt'), dtype=np.int64, ndmin=2)
        x_ok = len(rx) == len(x) and np.array_equal(rx[:, 0], np.arange(len(x)) - 4) \
            and np.array_equal(rx[:, 1:], x)
        y, y_ts, src_i = mb_decim_model(x)
        ry = np.loadtxt(os.path.join(HERE, f'mb_{name}_y.txt'), dtype=np.int64, ndmin=2)
        y_ok = len(ry) == len(y) and np.array_equal(ry[:, 0], y_ts) \
            and np.array_equal(ry[:, 1], np.arange(len(y))) and np.array_equal(ry[:, 2:], y)
        print(f"  fe_cond  x stream ({len(x)} samples): {'BIT-EXACT' if x_ok else 'MISMATCH'}")
        print(f"  mb_decim y stream ({len(y)} samples): {'BIT-EXACT' if y_ok else 'MISMATCH'}")
        ok &= x_ok and y_ok

        # ---- 1b. detectors + scheduler ----
        sw_recs, blank, confs, qual = sw_det_model(x, CFG)
        mb_recs = mb_det_model(y, y_ts, src_i, blank, CFG, qual, confs)
        evs = sched_model(sw_recs, mb_recs, confs, src_i, CFG)
        S, M, J, X, R, Lw, E = parse_ev(os.path.join(HERE, f'mb_{name}_ev.txt'))
        rec_ok = True
        for e in evs:
            if e['type'] == 0:
                r = S.get(e['id'])
                rec_ok &= r is not None and (r['ts'], r['ref'], r['feat']) == \
                    (e['ts'], e['ref'], [tuple(f) for f in e['feat']])
            else:
                r = M.get(e['id'])
                rec_ok &= r is not None and (r['ts'], r['idx'], r['ref'], r['rel']) == \
                    (e['ts'], e['idx'], e['ref'], e['rel'])
            rec_ok &= J.get(e['id'], {}).get('start') == e['start'] and J[e['id']]['type'] == e['type'] \
                and bin(J[e['id']]['men']).count('1') == len(e['masks'])
        rec_ok &= len(S) == sum(e['type'] == 0 for e in evs) and len(M) == sum(e['type'] == 1 for e in evs)
        print(f"  sw_det / mb_det / evt_sched: {len(S)} SW + {len(M)} MB records, {len(J)} jobs: "
              f"{'BIT-EXACT' if rec_ok else 'MISMATCH'}")
        ok &= rec_ok

        # ---- 2. engine results on the model's windows ----
        res_ok = True
        for e in evs:
            if e['type'] == 0:
                i0 = e['win_first'] + 4
                win = [list(x[i0:i0 + CFG['sw_len'], c]) + [0] * (N - CFG['sw_len']) for c in range(4)]
                for n in range(CFG['sw_len']):
                    tn = e['ts0'] + n
                    if any(a <= tn <= b for a, b in e['masks']):
                        for c in range(4):
                            win[c][n] = 0
                L, lo, hi, fl = CFG['sw_maxlag'], CFG['sw_band_lo'], CFG['sw_band_hi'], CFG['sw_floor']
            else:
                i0 = e['win_first']
                win = [list(y[i0:i0 + CFG['mb_len'], c]) + [0] * (N - CFG['mb_len']) for c in range(4)]
                for n in range(CFG['mb_len']):
                    tn = e['ts0'] + 5 * n
                    if any(a <= tn <= b for a, b in e['masks']):
                        for c in range(4):
                            win[c][n] = 0
                L, lo, hi, fl = CFG['mb_maxlag'], CFG['mb_band_lo'], CFG['mb_band_hi'], CFG['mb_floor']
            win = [[int(v) for v in w] for w in win]
            m = gcc_chain(win, L, lo, hi, fl)
            r = R.get(e['id'])
            if r is None:
                res_ok = False
                continue
            lags, pks, sides = [], [], []
            for j in range(6):
                cr = m['corr'][j]
                lw = [cr[(mm + N) % N] for mm in range(-L, L + 1)]
                lag = next(mm for mm in range(-L, L + 1) if lw[mm + L] == max(lw))
                side = max([abs(lw[mm + L]) for mm in range(-L, L + 1) if abs(mm - lag) > 2] + [0])
                lags.append(lag)
                pks.append(lw[lag + L])
                sides.append(side)
                lw_rtl = [Lw.get(e['id'], {}).get(k, [None] * 6)[j] for k in range(2 * L + 1)]
                res_ok &= lw_rtl == lw
            res_ok &= (r['lag'], r['pk'], r['side'], tuple(r['ez']), list(r['er'])) == \
                (lags, pks, sides, tuple(m['exp_z']), list(m['exp_r'])) and r['err'] == 0
            e['lw'] = {j: [Lw[e['id']][k][j] for k in range(2 * L + 1)] for j in range(6)} \
                if e['id'] in Lw and len(Lw[e['id']]) == 2 * L + 1 else None
            e['L'] = L
        print(f"  gcc_engine results ({len(R)} jobs): {'BIT-EXACT' if res_ok else 'MISMATCH'}")
        ok &= res_ok

        # ---- 3. health ----
        health = E is not None and E['underflow'] == 0 and E['late'] == 0 and E['ovf'] == 0 and not X
        print(f"  underflow {E['underflow']}, late {E['late']}, overflow {E['ovf']}  "
              f"({E['cycles'] / 1e8:.3f} s simulated at 100 MHz)")
        ok &= health

        # ---- 4. detection vs truth ----
        det_ok = True
        sws = sorted([e for e in evs if e['type'] == 0], key=lambda e: e['ts'])
        mbs = sorted([e for e in evs if e['type'] == 1], key=lambda e: e['ts'])
        dts = []
        print("  shot   SW onset err   MB onset err   Delta-t err   worst TDOA err (SW / MB)   MB bearing err (az / el)")
        for k, sh in enumerate(tr['shots']):
            tmb = min(sh['mb'])
            b = [e for e in mbs if abs(e['ts'] - tmb) < 500]
            if sh['sw'] is None:
                # subsonic bullet: no shockwave exists; expect a blast-only event
                if len(b) != 1 or b[0].get('lw') is None:
                    print(f"  {k}: no shockwave (subsonic); MB events {len(b)}  -> MISSED / EXTRA")
                    det_ok = False
                    continue
                b = b[0]
                tau = [refine(b['lw'][j], b['L']) / 2e4 for j in range(6)]
                w = max(abs(tau[j] - (sh['mb'][pb] - sh['mb'][pa]) / 1e5) * 1e6 for j, (pa, pb) in enumerate(PAIRS))
                az, el = bearing(tau, *detector_array())
                print(f"  {k}:   no shockwave   {(b['ts'] - tmb) * 10:+7.1f} us        n/a         "
                      f"  n/a / {w:5.1f} us           {ang_diff(az, tr['dir'][0]):+.2f} / {el - tr['dir'][1]:+.2f} deg"
                      f"   (blast-only: bearing, no range)")
                continue
            tsw = min(sh['sw'])
            s = [e for e in sws if abs(e['ts'] - tsw) < 500]
            if len(s) != 1 or len(b) != 1 or s[0].get('lw') is None or b[0].get('lw') is None:
                print(f"  {k}: SW events {len(s)}, MB events {len(b)}  -> MISSED / EXTRA")
                det_ok = False
                continue
            s, b = s[0], b[0]
            dterr = (b['ts'] - s['ts']) - (tmb - tsw)
            dts.append(dterr)
            worst = []
            for ev_, arr, fs in ((s, sh['sw'], 1e5), (b, sh['mb'], 2e4)):
                w = 0.0
                tau = []
                for j, (pa, pb) in enumerate(PAIRS):
                    t = refine(ev_['lw'][j], ev_['L']) / fs
                    tau.append(t)
                    w = max(w, abs(t - (arr[pb] - arr[pa]) / 1e5) * 1e6)
                worst.append(w)
                ev_['tau'] = tau
            az, el = bearing(b['tau'], *detector_array())
            print(f"  {k}:   {(s['ts'] - tsw) * 10:+7.1f} us     {(b['ts'] - tmb) * 10:+7.1f} us    "
                  f"{dterr * 10:+7.1f} us      {worst[0]:5.1f} / {worst[1]:5.1f} us"
                  f"           {ang_diff(az, tr['dir'][0]):+.2f} / {el - tr['dir'][1]:+.2f} deg")
        extra = len(sws) + len(mbs) - sum(2 if sh['sw'] is not None else 1 for sh in tr['shots'])
        if extra > 0:
            print(f"  {extra} extra events (false triggers)")
            det_ok = False
        elif extra < 0:
            print(f"  {-extra} events missing")
            det_ok = False
        nmask = sum(len(e['masks']) for e in evs)
        if nmask:
            print(f"  masks applied: {sum(len(e['masks']) for e in evs if e['type'] == 1)} SW-in-MB, "
                  f"{sum(len(e['masks']) for e in evs if e['type'] == 0)} MB-in-SW")
        # Delta-t budget: 500 us (timing_budget.md 4.4, set by the 50 m case), or
        # 0.5 % of the true gap at long range, where the range error it causes
        # (d_range/range = d_dt/dt) stays ~60x inside the +/-30 % spec.
        gaps = [min(sh['mb']) - min(sh['sw']) for sh in tr['shots'] if sh['sw'] is not None] or [0]
        dt_budget = max(50.0, 0.005 * max(gaps))
        dt_ok = (not dts and all(sh['sw'] is None for sh in tr['shots'])) or \
            (bool(dts) and max(abs(d) for d in dts) <= dt_budget)
        if dts:
            print(f"  worst Delta-t error {max(abs(d) for d in dts) * 10:.1f} us "
                  f"(budget {dt_budget * 10:.0f} us = max(500 us, 0.5% of {max(gaps) / 100:.0f} ms))")
        ok &= det_ok and dt_ok

        if name != 'burst' and pydir and tr['shots'][0]['sw'] is not None:
            rd = repo_detector(pydir, tr['wav'])
            if rd and rd[0] is not None:
                tsw, tmb = min(tr['shots'][0]['sw']), min(tr['shots'][0]['mb'])
                print(f"  python/gs_det_signal_prepare.py: SW onset err {(rd[0] - tsw) * 10:+.1f} us, "
                      f"MB onset err {((rd[1] or float('nan')) - tmb) * 10:+.1f} us (informational)")
        print(f"  -> {'PASS' if ok else 'FAIL'}")
        summary.append((name, ok))
        all_ok &= ok

    print("\n  (bearing uses the detector's own geometry and speed of sound from"
          " python/det_config.ini; generator ground truth is used for scoring only)")
    print("\n============================================================")
    for n, o in summary:
        print(f"  {n:6s} {'PASS' if o else 'FAIL'}")
    print(f" Milestone B: {'PASSED' if all_ok else 'FAILED'}")
    print("============================================================")
    return all_ok


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['prepare', 'check'])
    ap.add_argument('--python-dir', default=None)
    ap.add_argument('--only', nargs='*', default=None, help='check only these scenarios')
    ap.add_argument('--set', nargs='*', default=[], metavar='KEY=VALUE',
                    help='override configuration values (also written to mb_cfg.txt by prepare)')
    a = ap.parse_args()
    for kv in a.set:
        k, v = kv.split('=')
        if k not in CFG:
            sys.exit(f"unknown configuration key {k}")
        CFG[k] = int(v)
    if a.cmd == 'prepare':
        prepare(a.python_dir)
    else:
        sys.exit(0 if check(a.only) else 1)
