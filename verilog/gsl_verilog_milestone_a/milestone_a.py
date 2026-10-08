#!/usr/bin/env python3
"""
milestone_a.py -- Milestone A: GCC-PHAT core (gcc_engine) verification

  python3 milestone_a.py prepare   [--python-dir ../../python]
      Builds the job set for tb_gcc_engine.v:
        * 4 synthetic jobs (integer-delay noise, plane-wave N-wave,
          plane-wave Friedlander on the MB path, dead mic + band + floor)
        * realistic jobs from YOUR generator (gs_gen_clean_signal ->
          gs_gen_add_noise -> gs_gen_apply_adc), 3 scenarios, one SW and
          one MB job each, windowed around the generator's ground-truth
          arrival times. The MB path uses a Python stand-in for mb_decim
          (63-tap linear-phase FIR, /5, group delay removed).
      Writes ma_jobs.txt, ma_ring_<id>.hex, ma_prep.json.
      Ring start addresses are spread over the ring (incl. wrap-around).

  python3 milestone_a.py check
      Reads ma_out_<id>.txt from the testbench and reports, per job:
        1. loader + fft2048 + xspec_phat + IFFT + peak_search vs the
           bit-true model (must be bit-exact, incl. exponents, lags,
           peak values, sidelobes)
        2. fixed-point TDOA error: RTL correlations vs float GCC-PHAT on
           the same windows and band, both refined with the planned
           windowed-sinc method  (budget 1 us)
        3. TDOA error vs ground truth, and vs the gs_det_tdoa.gcc_phat
           algorithm (informational)
        4. MB jobs: shooter bearing (azimuth / elevation) from the 6 TDOAs
           vs ground truth (informational -- this is what PILAR's
           +/-2 deg azimuth spec is about)
"""
import argparse
import configparser
import json
import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'xspec_phat'))
sys.path.insert(0, os.path.join(HERE, '..', 'fft2048'))

from xspec_phat_model import (gcc_chain, coarse_peak, lag_window, refine,     # noqa: E402
                              gcc_float_corr, gcc_phat_repo, PAIRS, N)

RING = 4096
AMP = 20000.0
SW_CFG = dict(pre=300, len=1800, maxlag=100)          # det_config [tdoa], 100 kHz
MB_CFG = dict(pre=120, len=1720, maxlag=20)           # 6 + 80 ms @ 20 kHz
SW_BAND = (20, 512)                                   # ~1-25 kHz  (recommended)
MB_BAND = (10, 205)                                   # ~0.1-2 kHz (recommended)


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
def tetra(L=0.30):
    R, h = L / math.sqrt(3), L * math.sqrt(2 / 3)
    P = [[R * math.cos(math.radians(a)), R * math.sin(math.radians(a)), -h / 4] for a in (90, 210, 330)]
    P.append([0.0, 0.0, 3 * h / 4])
    return np.array(P)


def unit(az_deg, el_deg):
    az, el = math.radians(az_deg), math.radians(el_deg)
    return np.array([math.cos(el) * math.cos(az), math.cos(el) * math.sin(az), math.sin(el)])


def plane_delays(P, u, c, fs):
    """Arrival offset (samples) per mic for a plane wave coming FROM direction u."""
    return list(-(P @ u) / c * fs)


def sat16(v):
    return int(max(-32768, min(32767, round(v))))


def write_ring(fn, win, start):
    """win: 4 x len windows (ints). Place sample n at ring[(start+n) % 4096];
    fill the rest with a recognisable pattern (must never be read)."""
    ring = [[0x1AAAA, 0x15555, 0x1AAAA, 0x15555] for _ in range(RING)]   # junk
    for n in range(len(win[0])):
        ring[(start + n) % RING] = [win[i][n] for i in range(4)]
    with open(fn, 'w') as f:
        for w in ring:
            v = 0
            for c in w:                                  # {ch0, ch1, ch2, ch3}, 18 bits each
                v = (v << 18) | (c & 0x3FFFF)
            f.write(format(v, '018X') + '\n')
    return ring


def read_ring(fn):
    ring = []
    with open(fn) as f:
        for line in f:
            v = int(line.strip(), 16)
            ch = [(v >> s) & 0x3FFFF for s in (54, 36, 18, 0)]
            ring.append([c - (1 << 18) if c & 0x20000 else c for c in ch])
    return ring


def windows_from_ring(ring, start, length):
    """What gcc_loader does: 4 channels from the same ring address, zero-padded."""
    w = [[0] * N for _ in range(4)]
    for n in range(length):
        word = ring[(start + n) % RING]
        for i in range(4):
            w[i][n] = word[i]
    return w


def bearing(tdoa_s, P, c):
    """Least-squares plane-wave direction (towards the source) from 6 TDOAs.
    tau_ab = t_b - t_a = -((p_b - p_a) . u) / c."""
    D = np.array([P[b] - P[a] for a, b in PAIRS])
    u = np.linalg.lstsq(D, -c * np.asarray(tdoa_s), rcond=None)[0]
    u /= np.linalg.norm(u)
    return math.degrees(math.atan2(u[1], u[0])), math.degrees(math.asin(max(-1, min(1, u[2]))))


def ang_diff(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


# ----------------------------------------------------------------------
# synthetic jobs
# ----------------------------------------------------------------------
def synthetic_jobs(rng):
    jobs = []
    P = tetra()

    # 1. noise, integer delays
    d = [0, 7, -12, 25]
    nz = rng.uniform(-1, 1, SW_CFG['len'] + 64) * AMP
    win = [[sat16(nz[n + 32 - d[i]]) for n in range(SW_CFG['len'])] for i in range(4)]
    jobs.append(dict(name='synthetic noise, integer delays', type=0, start=1000, win=win,
                     d=[float(v) for v in d], fs=100000.0, band=(0, 1024), floor=0, **SW_CFG))

    # 2. plane-wave N-wave (SW), ring wrap-around
    fs = 100000.0
    dd = plane_delays(P, unit(37, 10), 343.0, fs)
    sg = lambda x: 0.5 * (1 + np.tanh(x / 2))                       # noqa: E731
    def nw(t): return (1 - 2 * t / 0.3e-3) * (sg(t / 5e-6) - sg((t - 0.3e-3) / 5e-6))
    n = np.arange(SW_CFG['len'])
    win = [[sat16(v) for v in AMP * nw((n - 300 - dd[i]) / fs) + rng.normal(0, AMP / 100, n.size)]
           for i in range(4)]
    jobs.append(dict(name='synthetic plane-wave N-wave az 37 el 10 (ring wrap)', type=0, start=4000,
                     win=win, d=[300 + v for v in dd], fs=fs, band=(0, 1024), floor=0,
                     P=P.tolist(), c=343.0, true_dir=[37.0, 10.0], **SW_CFG))

    # 3. plane-wave Friedlander (MB path, 20 kHz)
    fs = 20000.0
    dd = plane_delays(P, unit(200, 5), 343.0, fs)
    def fr(t): return np.where(t < -3e-3, 0.0, (1 - t / 2e-3) * np.exp(-np.clip(t, -3e-3, None) / 2e-3) * sg(t / 50e-6))
    n = np.arange(MB_CFG['len'])
    win = [[sat16(v) for v in AMP * fr((n - 120 - dd[i]) / fs) + rng.normal(0, AMP / 100, n.size)]
           for i in range(4)]
    jobs.append(dict(name='synthetic plane-wave Friedlander az 200 el 5 (MB path)', type=1, start=2500,
                     win=win, d=[120 + v for v in dd], fs=fs, band=(0, 1024), floor=0,
                     P=P.tolist(), c=343.0, true_dir=[200.0, 5.0], **MB_CFG))

    # 4. dead mic + band mask + floor
    d = [0, 7, -12, 0]
    nz = rng.uniform(-1, 1, SW_CFG['len'] + 64) * AMP
    win = [[sat16(nz[n + 32 - d[i]]) for n in range(SW_CFG['len'])] for i in range(3)]
    win.append([int(v) for v in rng.integers(-1, 2, SW_CFG['len'])])
    jobs.append(dict(name='synthetic dead mic ch3, band 20-500, floor 27', type=0, start=0, win=win,
                     d=[float(v) for v in d], fs=100000.0, band=(20, 500), floor=27, dead=[3],
                     **SW_CFG))
    return jobs


# ----------------------------------------------------------------------
# realistic jobs from the repo generator
# ----------------------------------------------------------------------
SCENARIOS = [
    dict(name='s1', desc='default: 200 m, miss 10 m, simple noise 0.005 Pa',
         patch={'noise': {'model': 'simple'}}),
    dict(name='s2', desc='explicit: shooter (350, 250, 20) m, simple noise 0.05 Pa',
         patch={'noise': {'model': 'simple', 'noise_floor_pa': '0.05'},
                'trajectory': {'mode': 'explicit', 'shooter_pos': '350,250,20',
                               'direction': '-350,-235,-20'}}),
    dict(name='s3', desc='600 m, miss 25 m, simple noise 0.02 Pa',
         patch={'noise': {'model': 'simple', 'noise_floor_pa': '0.02'},
                'trajectory': {'range': '600.0', 'y_miss': '25.0'}}),
]
STARTS = {'s1': (3000, 3950), 's2': (50, 2000), 's3': (4090, 700)}


def run_generator(pydir, workdir, sc):
    os.makedirs(workdir, exist_ok=True)
    cp = configparser.ConfigParser(inline_comment_prefixes=('#', ';'))
    cp.read(os.path.join(pydir, 'gen_config.ini'))
    for sec, kv in sc['patch'].items():
        for k, v in kv.items():
            cp[sec][k] = v
    cp['output']['basename'] = sc['name']
    ini = os.path.join(workdir, sc['name'] + '.ini')
    with open(ini, 'w') as f:
        cp.write(f)
    b = sc['name']
    steps = [['gs_gen_clean_signal.py', ini],
             ['gs_gen_add_noise.py', b + '_clean.wav', ini],
             ['gs_gen_apply_adc.py', b + '_noisy.wav', ini]]
    for s in steps:
        r = subprocess.run([sys.executable, os.path.join(pydir, s[0])] + s[1:], cwd=workdir,
                           capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stdout, r.stderr)
            raise RuntimeError(f"generator step {s[0]} failed for {b}")
    return os.path.join(workdir, b + '.wav'), os.path.join(workdir, b + '_clean.json')


def realistic_jobs(pydir, workdir):
    from scipy.io import wavfile
    from scipy.signal import firwin, lfilter
    jobs = []
    h = firwin(63, 7000, window=('kaiser', 6.0), fs=100000)      # mb_decim stand-in
    for sc in SCENARIOS:
        wav, truth = run_generator(pydir, workdir, sc)
        fs, x = wavfile.read(wav)
        x = x.astype(np.int64).T                                   # int16 ADC codes -> ring
        g = json.load(open(truth))['ground_truth']
        t0, c = g['t_master_start_s'], g['speed_of_sound_mps']
        P = np.array(g['mic_pos'])
        cen = P.mean(axis=0)
        sv = np.array(g['shooter_pos']) - cen
        sh_az = math.degrees(math.atan2(sv[1], sv[0]))
        sh_el = math.degrees(math.asin(sv[2] / np.linalg.norm(sv)))
        st_sw, st_mb = STARTS[sc['name']]

        # SW job, 100 kHz
        ta = [(p['shockwave']['t_arrive_s'] - t0) * fs for p in g['per_mic']]
        base = int(round(min(ta))) - SW_CFG['pre']
        win = [[int(v) for v in x[i, base:base + SW_CFG['len']]] for i in range(4)]
        jobs.append(dict(name=f"{sc['name']} SW: {sc['desc']}", type=0, start=st_sw, win=win,
                         d=[v - base for v in ta], fs=float(fs), band=SW_BAND, floor=0,
                         P=P.tolist(), c=c,
                         repo_win=[list(map(int, x[i, base:base + SW_CFG['len']])) for i in range(4)],
                         repo_fs=float(fs), repo_d=[v - base for v in ta], **SW_CFG))

        # MB job, decimated to 20 kHz (mb_decim stand-in)
        y = lfilter(h, 1.0, x.astype(float), axis=1)[:, 31::5]
        fsd = fs / 5
        ta = [(p['muzzle_blast']['t_arrive_s'] - t0) * fsd for p in g['per_mic']]
        base = int(round(min(ta))) - MB_CFG['pre']
        win = [[sat16(v) for v in y[i, base:base + MB_CFG['len']]] for i in range(4)]
        # full-rate window covering the same span, for the repo-algorithm comparison
        ta_full = [(p['muzzle_blast']['t_arrive_s'] - t0) * fs for p in g['per_mic']]
        bf = base * 5
        jobs.append(dict(name=f"{sc['name']} MB: {sc['desc']}", type=1, start=st_mb, win=win,
                         d=[v - base for v in ta], fs=fsd, band=MB_BAND, floor=0,
                         P=P.tolist(), c=c, true_dir=[sh_az, sh_el],
                         repo_win=[list(map(int, x[i, bf:bf + 5 * MB_CFG['len']])) for i in range(4)],
                         repo_fs=float(fs), repo_d=[v - bf for v in ta_full], **MB_CFG))
    return jobs


# ----------------------------------------------------------------------
# prepare
# ----------------------------------------------------------------------
def prepare(pydir):
    rng = np.random.default_rng(20261008)
    jobs = synthetic_jobs(rng)
    if os.path.isdir(pydir):
        jobs += realistic_jobs(pydir, os.path.join(HERE, 'gen_work'))
    else:
        print(f"NOTE: {pydir} not found -- realistic jobs skipped (use --python-dir)")
    meta = []
    with open(os.path.join(HERE, 'ma_jobs.txt'), 'w') as jf:
        jf.write(f"{len(jobs)}\n")
        for k, jb in enumerate(jobs, start=1):
            write_ring(os.path.join(HERE, f'ma_ring_{k}.hex'), jb['win'], jb['start'])
            jf.write(f"{k} {jb['type']} {jb['start']} {jb['len']} {jb['maxlag']} "
                     f"{jb['band'][0]} {jb['band'][1]} {jb['floor']}\n")
            m = {kk: v for kk, v in jb.items() if kk != 'win'}
            m['id'] = k
            meta.append(m)
            print(f"  job {k}: {jb['name']}  (ring start {jb['start']})")
    json.dump(meta, open(os.path.join(HERE, 'ma_prep.json'), 'w'))
    print(f"Wrote ma_jobs.txt, ma_ring_1..{len(jobs)}.hex, ma_prep.json")


# ----------------------------------------------------------------------
# check
# ----------------------------------------------------------------------
def read_out(fn):
    meta, rows = {}, []
    for line in open(fn):
        if line.startswith('#'):
            t = line[1:].split()
            for i in range(0, len(t) - 1, 2):
                meta[t[i]] = int(t[i + 1])
        elif line.strip():
            rows.append(list(map(int, line.split())))
    return meta, rows


def peak_model(r, L):
    lag = coarse_peak(r, L)
    pk = r[lag % N]
    side = max([abs(r[m % N]) for m in range(-L, L + 1) if abs(m - lag) > 2] + [0])
    return lag, pk, side


def check():
    jobs = json.load(open(os.path.join(HERE, 'ma_prep.json')))
    all_ok = True
    worst_fix_all = 0.0
    for jb in jobs:
        k = jb['id']
        L, fs, lo, hi, fl = jb['maxlag'], jb['fs'], jb['band'][0], jb['band'][1], jb['floor']
        us = 1e6 / fs
        meta, rows = read_out(os.path.join(HERE, f'ma_out_{k}.txt'))
        print(f"\n=== job {k}: {jb['name']}")
        print(f"  {meta['cycles']} cycles = {meta['cycles'] / 1e5:.3f} ms   "
              f"band {lo}-{hi} floor {fl}   ring start {jb['start']}")

        # 1. bit-true chain from the ring image (exercises loader addressing)
        ring = read_ring(os.path.join(HERE, f'ma_ring_{k}.hex'))
        x = windows_from_ring(ring, jb['start'], jb['len'])
        m = gcc_chain(x, L, lo, hi, fl)
        col = lambda c: [r[c] for r in rows]                                # noqa: E731
        z_ok = (list(zip(col(0), col(1))) == m['Z0'] and list(zip(col(2), col(3))) == m['Z1']
                and (meta['ez0'], meta['ez1']) == m['exp_z'])
        r_ok = all(col(4 + 2 * p) == list(m['R'][p][0]) and col(5 + 2 * p) == list(m['R'][p][1])
                   for p in range(3)) and [meta[f'er{p}'] for p in range(3)] == m['exp_r']
        pk_ok = True
        for j in range(6):
            lag, pk, side = peak_model(m['corr'][j], L)
            pk_ok &= (meta[f'lag{j}'], meta[f'pk{j}'], meta[f'side{j}']) == (lag, pk, side)
        ok = z_ok and r_ok and pk_ok and meta['err'] == 0
        print(f"  bit-exact: loader+FFT (Z) {'OK' if z_ok else 'MISMATCH'}, "
              f"xspec+IFFT (R) {'OK' if r_ok else 'MISMATCH'}, "
              f"peak_search {'OK' if pk_ok else 'MISMATCH'}")

        # 2-3. accuracy
        corr = m['corr']
        dead = jb.get('dead', [])
        d = jb['d']
        xa = [np.array(x[i][:jb['len']], float) for i in range(4)]
        print("  pair     truth    FPGA   float | fixpt err  | FPGA-truth   PSR" +
              ("   repo-truth" if 'repo_win' in jb else ""))
        tdoa_fpga = []
        worst_fix = 0.0
        for j, (a, b) in enumerate(PAIRS):
            lag, pk, side = meta[f'lag{j}'], meta[f'pk{j}'], meta[f'side{j}']
            psr = pk / side if side > 0 else float('inf')
            if a in dead or b in dead:
                print(f"  ({a},{b})   dead channel: peak {pk}, sidelobe {side} -> flagged")
                ok &= pk <= 2
                tdoa_fpga.append(None)
                continue
            truth = (d[b] - d[a]) * us
            win = lag_window(corr[j], L)
            t_f = refine(win, L) * us
            fl_c = gcc_float_corr(np.pad(xa[a], (0, N - jb['len'])), np.pad(xa[b], (0, N - jb['len'])),
                                  N, lo, hi)
            fwin = [fl_c[(mm + N) % N] for mm in range(-L, L + 1)]
            t_fl = refine(fwin, L) * us
            fix = (refine(win, L, up=256) - refine(fwin, L, up=256)) * us
            worst_fix = max(worst_fix, abs(fix))
            tdoa_fpga.append(t_f * 1e-6)
            line = (f"  ({a},{b}) {truth:+8.1f} {t_f:+7.1f} {t_fl:+7.1f} | {fix:+6.3f} us  | "
                    f"{t_f - truth:+7.1f} us {psr:6.1f}")
            if 'repo_win' in jb:
                rw = [np.array(v, float) for v in jb['repo_win']]
                t_r = gcc_phat_repo(rw[a], rw[b], jb['repo_fs'], 16, L / fs) * 1e6
                tr_r = (jb['repo_d'][b] - jb['repo_d'][a]) * 1e6 / jb['repo_fs']
                line += f"   {t_r - tr_r:+8.1f} us"
            print(line)
        print(f"  worst fixed-point TDOA error: {worst_fix:.3f} us (budget 1 us)")
        worst_fix_all = max(worst_fix_all, worst_fix)
        ok &= worst_fix <= 1.0 and meta['cycles'] <= 200000

        # 4. bearing (jobs with a known source direction)
        if 'true_dir' in jb and None not in tdoa_fpga:
            P = np.array(jb['P'])
            az, el = bearing(tdoa_fpga, P, jb['c'])
            taz, tel = jb['true_dir']
            print(f"  bearing from FPGA TDOAs: az {az:7.2f} el {el:6.2f}   truth az {taz:7.2f} el {tel:6.2f}"
                  f"   error az {ang_diff(az, taz):+.2f} el {el - tel:+.2f} deg")
        print(f"  -> {'PASS' if ok else 'FAIL'}")
        all_ok &= ok

    print("\n============================================================")
    print(f" Milestone A: {'PASSED' if all_ok else 'FAILED'}  "
          f"(bit-exact, fixed-point error <= 1 us [worst {worst_fix_all:.3f}], <= 2 ms per job)")
    print("============================================================")
    return all_ok


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['prepare', 'check'])
    ap.add_argument('--python-dir', default=os.path.join(HERE, '..', '..', 'python'))
    a = ap.parse_args()
    if a.cmd == 'prepare':
        prepare(a.python_dir)
    else:
        sys.exit(0 if check() else 1)
