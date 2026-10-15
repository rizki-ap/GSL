#!/usr/bin/env python3
"""Milestone C1 check: gsl_mm (registers, mailbox, replay DMA) against Milestone B.

Every Milestone B scenario is replayed through gsl_mm: the recording comes from
an SDRAM model through the replay DMA, and a processor model reads all results
over the Avalon-MM slave. Pass criteria:

  1. x and y streams byte-identical to Milestone B (replay path delivers every
     frame in order),
  2. every record, job, late job, result and lag window identical to the
     Milestone B log (compared as sets of lines: the mailbox may order events
     that happen within microseconds differently),
  3. no mailbox drops, late jobs, overflows or underflows,
  4. overflow test: a slow processor loses whole messages (counted), never
     parts of one,
  5. loop test: 3 passes of one recording give the same events every pass.

Needs the Milestone B outputs: run `make` (or `make prepare sim`) in ../stream first.

Usage: python3 c1_check.py [--only name ...] [--skip-stress]
"""
import os, sys, subprocess, re, argparse

HERE = os.path.dirname(os.path.abspath(__file__))
STR = os.path.normpath(os.path.join(HERE, '..', 'stream'))
SIM = os.path.join(HERE, 'obj_dir', 'Vtb_gsl_mm')


def scenarios():
    lst = os.path.join(STR, 'mb_list.txt')
    if not os.path.exists(lst):
        sys.exit("../stream/mb_list.txt not found: run `make` (or `make prepare sim`) in ../stream first")
    return [(l.split()[0], int(l.split()[1])) for l in open(lst) if l.strip()]


def run(name, n, out, extra=()):
    cmd = [SIM, f'+stim={STR}/mb_{name}.hex', f'+nframes={n}', f'+cfg={STR}/mb_cfg.txt', f'+out={out}', *extra]
    r = subprocess.run(cmd, cwd=HERE, capture_output=True, text=True)
    txt = r.stdout
    def num(label):
        m = re.search(label + r'\s*:\s*(\d+)', txt)
        return int(m.group(1)) if m else -1
    st = re.search(r'STATUS\s*:\s*([0-9a-fA-F]+)', txt)
    return dict(late=num('late jobs'), ovf=num('pending overflow'), drops=num('mailbox drops'),
                status=int(st.group(1), 16) if st else -1, passes=num('replay passes'), log=txt)


def ev_lines(path):
    return [l.rstrip('\n') for l in open(path) if l.strip() and not l.startswith('E ')]


def check_scenario(name, n):
    out = f'c1_{name}'
    ref = {k: os.path.join(STR, f'mb_{name}_{k}.txt') for k in ('x', 'y', 'ev')}
    if not all(os.path.exists(p) for p in ref.values()):
        return False, "Milestone B output missing: run `make sim` in ../stream"
    r = run(name, n, os.path.join(HERE, out))
    probs = []
    for k in ('x', 'y'):
        if open(os.path.join(HERE, f'{out}_{k}.txt'), 'rb').read() != open(ref[k], 'rb').read():
            probs.append(f'{k} stream differs')
    a, b = ev_lines(os.path.join(HERE, f'{out}_ev.txt')), ev_lines(ref['ev'])
    if sorted(a) != sorted(b):
        probs.append(f'events differ ({len(a)} vs {len(b)} lines)')
    order = 'same order' if a == b else 'same content, some events reordered'
    if r['drops'] != 0 or r['late'] != 0 or r['ovf'] != 0:
        probs.append(f"drops {r['drops']} late {r['late']} overflow {r['ovf']}")
    if r['status'] < 0 or r['status'] & 0x11:
        probs.append(f"STATUS {r['status']:#x} (underflow or result error)")
    msg = f"{len(b)} log lines identical ({order}); x, y streams identical" if not probs else '; '.join(probs)
    return not probs, msg


def overflow_test():
    """A processor 8000 cycles slower per word than the hardware: the 8192-word ring
    fills during the burst. Whole messages may be dropped (counted); every message
    that does arrive must be complete and correct."""
    name = 'burst'
    n = dict(scenarios())[name]
    out = os.path.join(HERE, 'c1_ovf')
    r = run(name, n, out, ['+slow=8000'])
    got = ev_lines(out + '_ev.txt')
    ref = set(ev_lines(os.path.join(STR, f'mb_{name}_ev.txt')))
    foreign = [l for l in got if l not in ref]
    # every received result must carry its complete lag window
    torn = 0
    for l in got:
        if l.startswith('R '):
            rid, rtype = l.split()[1], int(l.split()[2])
            L = 20 if rtype else 100
            nl = sum(1 for x in got if x.startswith(f'L {rid} '))
            if nl != 2 * L + 1:
                torn += 1
    ok = r['drops'] > 0 and not foreign and torn == 0
    nres = sum(1 for l in got if l.startswith('R '))
    return ok, (f"{r['drops']} messages dropped and counted; {nres} of 16 results arrived, "
                f"all complete and identical to Milestone B" if ok else
                f"drops {r['drops']}, foreign lines {len(foreign)}, torn results {torn}")


def loop_test():
    """Loop mode, 3 passes of one recording: same events in every pass."""
    name = 's1'
    n = dict(scenarios())[name]
    out = os.path.join(HERE, 'c1_loop')
    r = run(name, n, out, ['+loop=3'])
    got = ev_lines(out + '_ev.txt')
    S = [l.split() for l in got if l.startswith('S ')]
    M = [l.split() for l in got if l.startswith('M ')]
    R = [l.split() for l in got if l.startswith('R ')]
    passes = r['passes']
    detail = []
    ok = passes >= 3 and len(S) == len(M) == passes and len(R) == 2 * passes
    if ok:
        for p in range(1, passes):
            dts = int(S[p][2]) - int(S[0][2]) - p * n           # onset vs pass 0, in samples
            dtm = int(M[p][2]) - int(M[0][2]) - p * n
            lag_sw = R[2 * p][8:14] == R[0][8:14]
            lag_mb = R[2 * p + 1][8:14] == R[1][8:14]
            detail.append(f"pass {p + 1}: SW onset {dts:+d}, MB onset {dtm:+d} samples, "
                          f"coarse lags {'equal' if lag_sw and lag_mb else 'differ'}")
            if abs(dts) > 1 or abs(dtm) > 5 or not (lag_sw and lag_mb):
                ok = False
    bad = r['drops'] or r['late'] or r['ovf'] or (r['status'] & 0x11)
    ok = ok and not bad
    return ok, f"{passes} passes, {len(S)} SW / {len(M)} MB events, {len(R)} results; " + '; '.join(detail)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--only', nargs='*')
    ap.add_argument('--skip-stress', action='store_true')
    a = ap.parse_args()
    if not os.path.exists(SIM):
        sys.exit("build first: make build")
    allok = True
    print("============================================================")
    print(" Milestone C1: gsl_mm (registers, mailbox, replay DMA)")
    print("============================================================")
    for name, n in scenarios():
        if a.only and name not in a.only:
            continue
        ok, msg = check_scenario(name, n)
        allok &= ok
        print(f"  {name:10s} {'PASS' if ok else 'FAIL'}  {msg}")
    if not a.skip_stress and not a.only:
        for label, fn in (('overflow', overflow_test), ('loop', loop_test)):
            ok, msg = fn()
            allok &= ok
            print(f"  {label:10s} {'PASS' if ok else 'FAIL'}  {msg}")
    print("============================================================")
    print(f" Milestone C1: {'PASSED' if allok else 'FAILED'}")
    print("============================================================")
    sys.exit(0 if allok else 1)


if __name__ == '__main__':
    main()
