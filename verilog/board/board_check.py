#!/usr/bin/env python3
"""GSL board check (Milestone C3).

  python3 board_check.py prepare [name ...]   recordings for SDRAM: ../stream/mb_<name>.hex
                                              -> rec/mb_<name>.bin (all scenarios if none given)
  python3 board_check.py check <name> [log]   compare a board log (default log/<name>.bin,
                                              saved by gsl_wait_dump) with Milestone B
  python3 board_check.py selftest <raw.txt> <name>
                                              test this script on the C1 testbench's raw
                                              mailbox words (+raw=...), no board needed

Pass criterion for `check`: every record, job, result and lag window the board
produced is identical to the Milestone B simulation (../stream/mb_<name>_ev.txt),
with no late jobs, overflows, mailbox drops, underflows or result errors.
"""
import os, sys, struct

HERE = os.path.dirname(os.path.abspath(__file__))
STR = os.path.normpath(os.path.join(HERE, '..', 'stream'))
LOG_MAGIC = 0x47534C4C


def s(v, bits):
    v &= (1 << bits) - 1
    return v - (1 << bits) if v >> (bits - 1) else v


def prepare(names):
    lst = os.path.join(STR, 'mb_list.txt')
    if not names:
        if not os.path.exists(lst):
            sys.exit("../stream/mb_list.txt not found: run `make prepare` in ../stream")
        names = [l.split()[0] for l in open(lst) if l.strip()]
    os.makedirs(os.path.join(HERE, 'rec'), exist_ok=True)
    for n in names:
        src = os.path.join(STR, f'mb_{n}.hex')
        frames = [int(l, 16) for l in open(src) if l.strip()]
        dst = os.path.join(HERE, 'rec', f'mb_{n}.bin')
        with open(dst, 'wb') as f:
            f.write(b''.join(struct.pack('<Q', v) for v in frames))   # low word first
        print(f"  {n:10s} {len(frames):7d} frames -> rec/mb_{n}.bin ({8 * len(frames) / 1e6:.1f} MB)")


def messages_to_lines(words):
    """Rebuild the Milestone B log lines (S / M / J / X / L / R) from mailbox words."""
    out, i, bad = [], 0, 0
    while i < len(words):
        h = words[i]
        typ, ln, mid = h >> 28, (h >> 16) & 0xFFF, h & 0xFFFF
        if ln == 0 or i + ln > len(words):
            bad += 1
            break
        m = words[i:i + ln]
        if typ == 1:
            t = s(m[1] | (m[2] << 32), 64)
            line = f"S {mid} {t} {m[3]}"
            for c in range(4):
                line += f" {s(m[4 + 3 * c], 32)} {m[6 + 3 * c] & 0x1FF} {s(m[5 + 3 * c], 32)} {(m[6 + 3 * c] >> 16) & 0x1FF}"
            out.append(line)
        elif typ == 2:
            out.append(f"M {mid} {s(m[1] | (m[2] << 32), 64)} {m[3]} {m[4] & 3} {(m[4] >> 8) & 1}")
        elif typ == 3:
            out.append(f"J {mid} {m[1] >> 31} {m[1] & 0xFFF} {(m[1] >> 16) & 0xF}")
        elif typ == 4:
            out.append(f"X {mid}")
        elif typ == 5:
            L = (m[1] >> 16) & 0x7F
            for k in range(2 * L + 1):
                out.append(f"L {mid} {k} " + ' '.join(str(s(m[17 + 6 * k + j], 32)) for j in range(6)))
            e = m[2]
            line = f"R {mid} {m[1] >> 31} {(m[1] >> 30) & 1} {e & 31} {(e >> 5) & 31} {(e >> 10) & 31} {(e >> 15) & 31} {(e >> 20) & 31}"
            lags = [s(m[3] >> (8 * j), 8) for j in range(4)] + [s(m[4] >> (8 * j), 8) for j in range(2)]
            line += ''.join(f" {v}" for v in lags)
            line += ''.join(f" {s(m[5 + j], 32)}" for j in range(6))
            line += ''.join(f" {m[11 + j]}" for j in range(6))
            out.append(line)
        else:
            bad += 1
        i += ln
    return out, bad


def compare(name, words, hdr=None):
    ref_path = os.path.join(STR, f'mb_{name}_ev.txt')
    if not os.path.exists(ref_path):
        sys.exit(f"{ref_path} not found: run `make` (or `make prepare sim`) in ../stream")
    ref = [l.rstrip('\n') for l in open(ref_path) if l.strip() and not l.startswith('E ')]
    got, bad = messages_to_lines(words)
    ok = sorted(got) == sorted(ref) and bad == 0
    print(f"  {len(words)} words, {len(got)} log lines ({len(ref)} expected)"
          + (", same order" if got == ref else ""))
    if hdr:
        st, late, ovf, drops, ev, res, pos, nfr = hdr
        print(f"  board: SW {ev & 0xFFFF}, MB {ev >> 16}, results {res}; frames {pos}/{nfr}; "
              f"late {late}, overflow {ovf}, mailbox drops {drops}; STATUS 0x{st:02X}")
        if late or ovf or drops or (st & 0x11) or pos != nfr:
            ok = False
            print("  FAIL: errors reported by the board (or not every frame replayed)")
    if not ok:
        missing = sorted(set(ref) - set(got))[:3]
        extra = sorted(set(got) - set(ref))[:3]
        for l in missing: print("  missing:", l[:100])
        for l in extra:   print("  unexpected:", l[:100])
        if bad: print(f"  {bad} malformed message(s)")
    print(f" {name}: {'PASS - identical to simulation' if ok else 'FAIL'}")
    return ok


def check(name, path):
    raw = open(path, 'rb').read()
    w = struct.unpack(f'<{len(raw) // 4}I', raw[:len(raw) // 4 * 4])
    if w[0] != LOG_MAGIC:
        sys.exit(f"{path}: not a GSL log (magic 0x{w[0]:08X})")
    nw = w[1]
    words = list(w[16:16 + nw])
    if len(words) != nw:
        sys.exit(f"{path}: truncated ({len(words)} of {nw} words)")
    return compare(name, words, hdr=(w[2], w[3], w[4], w[5], w[6], w[7], w[8], w[9]))


def main():
    a = sys.argv[1:]
    if not a:
        sys.exit(__doc__)
    if a[0] == 'prepare':
        prepare(a[1:])
    elif a[0] == 'check' and len(a) >= 2:
        sys.exit(0 if check(a[1], a[2] if len(a) > 2 else os.path.join(HERE, 'log', f'{a[1]}.bin')) else 1)
    elif a[0] == 'selftest' and len(a) == 3:
        words = [int(l, 16) for l in open(a[1]) if l.strip()]
        sys.exit(0 if compare(a[2], words) else 1)
    else:
        sys.exit(__doc__)


if __name__ == '__main__':
    main()
