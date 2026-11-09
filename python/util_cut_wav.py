#!/usr/bin/env python3
"""
util_cut_wav.py
===============
Cuts a multichannel .wav from PRE seconds before the first peak to POST
seconds after the last peak (default 1 s / 1 s). Useful for trimming long
recordings down to the part that holds the event(s) before plotting or
running the gs_det_* pipeline.

A "peak" is any sample whose |amplitude| reaches the threshold on any of the
detection channels. The threshold is either:
  - relative to the loudest sample in the file (default: -20 dB, i.e. 10 % of
    the maximum), or
  - absolute in dBFS with -abs (integer PCM files only, e.g. -th -30 -abs).

The cut keeps the original sample format (16-bit, 24-bit, 32-bit, float) and
all channels. If PRE/POST reach past the file edges, the cut stops at the edge.

Usage
-----
  python util_cut_wav.py rec.wav                      # -> rec_cut.wav, 1 s before/after
  python util_cut_wav.py rec.wav -pre 0.5 -post 2
  python util_cut_wav.py rec.wav -th -12              # only peaks within 12 dB of the loudest
  python util_cut_wav.py rec.wav -th -30 -abs         # absolute threshold, -30 dBFS
  python util_cut_wav.py rec.wav -chan 0 1 3 5        # detect peaks on these channels only
  python util_cut_wav.py rec.wav -o short.wav
"""

import argparse
import os

import numpy as np
import soundfile as sf

INT_FULL_SCALE = {"PCM_U8": 2.0**7, "PCM_S8": 2.0**7, "PCM_16": 2.0**15,
                  "PCM_24": 2.0**31, "PCM_32": 2.0**31}   # 24-bit is read as int32


def non_negative(value):
    v = float(value)
    if v < 0:
        raise argparse.ArgumentTypeError("must be >= 0")
    return v


def main():
    p = argparse.ArgumentParser(
        description="Cut a wav from PRE s before the first peak to POST s after the last peak.",
        add_help=False)
    p.add_argument("-h", "-help", "--help", action="help", help="show this help message and exit")
    p.add_argument("wav_path", help="input .wav file")
    p.add_argument("-o", default=None, metavar="OUT",
                   help="output .wav (default: <input>_cut.wav)")
    p.add_argument("-pre", type=non_negative, default=0.5, metavar="SEC",
                   help="seconds kept before the first peak (default: 1)")
    p.add_argument("-post", type=non_negative, default=5.0, metavar="SEC",
                   help="seconds kept after the last peak (default: 1)")
    p.add_argument("-th", type=float, default=-20.0, metavar="DB",
                   help="peak threshold in dB, relative to the loudest sample "
                        "(default: -20); with -abs it is in dBFS")
    p.add_argument("-abs", action="store_true",
                   help="treat -th as an absolute level in dBFS instead of relative")
    p.add_argument("-chan", type=int, nargs="+", default=None, metavar="CH",
                   help="channels (0-based) used to find peaks (default: all); "
                        "all channels are always kept in the output")
    args = p.parse_args()

    if args.th > 0:
        p.error("-th must be <= 0 dB")

    info = sf.info(args.wav_path)
    is_int = info.subtype in INT_FULL_SCALE
    data, fs = sf.read(args.wav_path, dtype="int32" if is_int else "float64", always_2d=True)
    n, n_ch = data.shape
    print(f"{args.wav_path}: {n_ch} ch @ {fs} Hz, {info.subtype}, {n / fs:.3f} s")

    chans = list(range(n_ch)) if args.chan is None else list(dict.fromkeys(args.chan))
    bad = [c for c in chans if not 0 <= c < n_ch]
    if bad:
        p.error(f"-chan {bad} out of range: file has channels 0..{n_ch - 1}")

    level = np.abs(data[:, chans].astype(np.float64)).max(axis=1)
    loudest = level.max()
    if loudest == 0:
        raise SystemExit("File is silent -- nothing to cut.")

    if args.abs:
        if not is_int:
            p.error("-abs needs an integer PCM wav (float files have no fixed full scale)")
        thr = 10 ** (args.th / 20) * INT_FULL_SCALE[info.subtype]
        thr_txt = f"{args.th:g} dBFS"
    else:
        thr = 10 ** (args.th / 20) * loudest
        thr_txt = f"{args.th:g} dB re loudest sample"

    hits = np.flatnonzero(level >= thr)
    if hits.size == 0:
        full_scale = INT_FULL_SCALE.get(info.subtype)
        peak_txt = f"{20 * np.log10(loudest / full_scale):.1f} dBFS" if full_scale else f"{loudest:.4g}"
        raise SystemExit(f"No sample reaches the threshold ({thr_txt}); loudest is {peak_txt}.")

    first, last = int(hits[0]), int(hits[-1])
    lo = max(0, first - int(round(args.pre * fs)))
    hi = min(n, last + int(round(args.post * fs)) + 1)

    print(f"Threshold {thr_txt} on ch {chans}: first peak at {first / fs:.4f} s, "
          f"last peak at {last / fs:.4f} s ({hits.size} samples above threshold)")
    if first - lo < round(args.pre * fs):
        print(f"  note: only {(first - lo) / fs:.3f} s available before the first peak")
    if hi - 1 - last < round(args.post * fs):
        print(f"  note: only {(hi - 1 - last) / fs:.3f} s available after the last peak")

    out = args.o or os.path.splitext(args.wav_path)[0] + "_cut.wav"
    sf.write(out, data[lo:hi], fs, subtype=info.subtype)
    print(f"Saved {out}: {lo / fs:.4f} s .. {hi / fs:.4f} s of the input "
          f"({(hi - lo) / fs:.3f} s, {n_ch} ch, {info.subtype})")


if __name__ == "__main__":
    main()
