#!/usr/bin/env python3
"""
util_trigger_record_uma8.py
===========================
Listens to the UMA-8 continuously and saves a short WAV clip only when the
signal goes above a threshold: PRE seconds before the peak and POST seconds
after it (default 0.5 s / 0.5 s). Keeps listening for more events until
Ctrl+C or until -n events have been saved.

How it works
------------
- Audio is captured continuously into a rolling buffer, so the PRE part is
  already in memory when the trigger fires.
- Trigger: the first sample on any monitored channel with |x| >= threshold.
- The clip is centred on the actual peak, searched within -look seconds
  after the trigger (default 5 ms), not on the first threshold crossing.
- After a clip is saved, detection resumes right after the clip ends.

Output modes
------------
- Default (multiplex): all clips go into ONE wav, uma8_events_<start time>.wav,
  separated by -gap seconds of silence (default 0 s). A matching .csv lists
  each event's wall-clock time and where its clip and peak sit in the wav.
  The file is finalized after every event, so it stays readable even if the
  script is killed.
- -split: each clip is saved as its own wav, uma8_<event time>.wav.

Usage
-----
  python util_trigger_record_uma8.py                    # -12 dBFS, 0.5 s / 0.5 s, 24-bit, one wav
  python util_trigger_record_uma8.py -split             # one wav per event
  python util_trigger_record_uma8.py -gap 0             # multiplex with no silence between clips
  python util_trigger_record_uma8.py -th -20            # more sensitive
  python util_trigger_record_uma8.py -pre 0.3 -post 1.0
  python util_trigger_record_uma8.py -n 5 -o shots      # stop after 5 events, save in .\\shots
  python util_trigger_record_uma8.py -mon 0 1 3 5       # trigger only on these channels
  python util_trigger_record_uma8.py -bits 16
  python util_trigger_record_uma8.py -max_wait 60       # give up after 60 s without a trigger

It stops on Ctrl+C, after -n events, or when no trigger has happened for
-max_wait seconds (default 120 s, counted from start and reset after each
event; -max_wait 0 waits forever).

While running it prints the live peak level (dBFS) once per second, which
helps pick a threshold: set -th a few dB above the background level.
"""

import argparse
import os
import queue
import sys
import time
from datetime import datetime

import numpy as np
import sounddevice as sd
import soundfile as sf

MIN_CHANNELS = 6


def find_uma8(min_channels):
    """Find the UMA-8 input device, trying WASAPI first (MME often caps at 2 channels)."""
    hostapis = sd.query_hostapis()
    order = sorted(range(len(hostapis)),
                   key=lambda i: 0 if "WASAPI" in hostapis[i]["name"] else 1)
    for api in order:
        for idx, dev in enumerate(sd.query_devices()):
            name = dev["name"].lower()
            if (dev["hostapi"] == api
                    and dev["max_input_channels"] >= min_channels
                    and any(k in name for k in ("minidsp", "uma", "vocalfusion"))):
                return idx, dev
    return None, None


def positive_float(value):
    v = float(value)
    if v <= 0:
        raise argparse.ArgumentTypeError("must be > 0")
    return v


def parse_args():
    p = argparse.ArgumentParser(
        description="Save UMA-8 clips around peaks above a threshold.", add_help=False)
    p.add_argument("-h", "-help", "--help", action="help",
                   help="show this help message and exit")
    p.add_argument("-th", type=float, default=-12.0, metavar="DBFS",
                   help="trigger threshold in dBFS, e.g. -12 (default: -12)")
    p.add_argument("-pre", type=positive_float, default=0.5, metavar="SEC",
                   help="seconds kept before the peak (default: 0.5)")
    p.add_argument("-post", type=positive_float, default=0.5, metavar="SEC",
                   help="seconds kept after the peak (default: 0.5)")
    p.add_argument("-look", type=positive_float, default=0.005, metavar="SEC",
                   help="window after the trigger to search for the true peak (default: 0.005)")
    p.add_argument("-mon", type=int, nargs="+", default=None, metavar="CH",
                   help="channel indices used for triggering (default: all)")
    p.add_argument("-n", type=int, default=0,
                   help="stop after N saved events (default: 0 = run until Ctrl+C)")
    p.add_argument("-o", default=".", metavar="DIR",
                   help="output folder (default: current folder)")
    p.add_argument("-bits", type=int, choices=(16, 24), default=24,
                   help="WAV sample resolution in bits (default: 24)")
    p.add_argument("-max_wait", type=float, default=120.0, metavar="SEC",
                   help="stop if no trigger for this many seconds, counted from start "
                        "and reset after each event (default: 120, 0 = no limit)")
    p.add_argument("-split", action="store_true",
                   help="save each event in its own wav (default: all events in one wav)")
    p.add_argument("-gap", type=float, default=0, metavar="SEC",
                   help="silence between events in the single wav (default: 0)")
    args = p.parse_args()
    if args.th >= 0:
        p.error("-th must be below 0 dBFS")
    if args.max_wait < 0:
        p.error("-max_wait must be >= 0")
    if args.gap < 0:
        p.error("-gap must be >= 0")
    return args


def main():
    args = parse_args()

    device, info = find_uma8(MIN_CHANNELS)
    if device is None:
        print(f"No UMA-8 input with {MIN_CHANNELS}+ channels found. Available devices:\n")
        print(sd.query_devices())
        raise SystemExit("Check the firmware and the Windows recording format.")

    hostapi = sd.query_hostapis(info["hostapi"])["name"]
    # WASAPI shared mode only accepts the device's own channel count and rate.
    channels = info["max_input_channels"]
    fs = int(info["default_samplerate"]) if "WASAPI" in hostapi else 48000

    mon = list(range(channels)) if args.mon is None else args.mon
    bad = [c for c in mon if not 0 <= c < channels]
    if bad:
        raise SystemExit(f"-mon channels {bad} out of range 0..{channels - 1}")

    # 24-bit is captured as int32 (sounddevice can't do packed 24-bit) and saved as PCM_24.
    rec_dtype, wav_subtype, full_scale = {16: ("int16", "PCM_16", 2.0**15),
                                          24: ("int32", "PCM_24", 2.0**31)}[args.bits]
    thr = 10 ** (args.th / 20) * full_scale

    pre_n = int(round(args.pre * fs))
    post_n = int(round(args.post * fs))
    look_n = max(1, int(round(args.look * fs)))
    keep_n = pre_n + look_n          # history kept while waiting for a trigger

    os.makedirs(args.o, exist_ok=True)

    print(f"Using [{device}] {info['name']} ({hostapi}, {channels} ch @ {fs} Hz, {args.bits}-bit)")
    print(f"Trigger: {args.th:.1f} dBFS on ch {mon}; clip = {args.pre}s before + "
          f"{args.post}s after peak")

    # Single-wav (multiplex) mode: one open file plus a CSV index of events.
    mux_wav = mux_csv = None
    mux_pos = 0           # samples written to the single wav so far
    gap = np.zeros((int(round(args.gap * fs)), channels), dtype=rec_dtype)
    if args.split:
        print(f"Saving one wav per event in {os.path.abspath(args.o)}")
    else:
        stem = os.path.join(args.o, "uma8_events_" + datetime.now().strftime("%Y%m%d_%H%M%S"))
        mux_path = stem + ".wav"
        mux_wav = sf.SoundFile(mux_path, "w", samplerate=fs, channels=channels,
                               subtype=wav_subtype)
        mux_csv = open(stem + ".csv", "w", newline="")
        mux_csv.write("event,wall_time,clip_start_s,clip_end_s,peak_s,peak_dbfs\n")
        print(f"Saving all events in {os.path.abspath(mux_path)} (index: {stem}.csv)")
    wait_txt = f"max wait {args.max_wait:g}s per trigger" if args.max_wait else "no time limit"
    print(f"Listening... ({wait_txt}; Ctrl+C to stop)\n")

    q = queue.Queue()
    overflows = [0]

    def callback(indata, frames, time_info, status):
        if status.input_overflow:
            overflows[0] += 1
        q.put(indata.copy())

    buf = np.empty((0, channels), dtype=rec_dtype)
    searched = 0          # rows of buf already checked for a trigger
    trig = None           # row index of the trigger in buf
    trig_time = None
    saved = 0
    level_max = 0.0
    last_print = time.monotonic()

    try:
        with sd.InputStream(device=device, channels=channels, samplerate=fs,
                            dtype=rec_dtype, callback=callback):
            wait_start = time.monotonic()
            timed_out = False
            while args.n == 0 or saved < args.n:
                if (args.max_wait and trig is None
                        and time.monotonic() - wait_start >= args.max_wait):
                    timed_out = True
                    break
                try:
                    block = q.get(timeout=0.5)
                except queue.Empty:
                    continue
                buf = np.concatenate((buf, block))

                blk_abs = np.abs(block[:, mon].astype(np.float64)).max(axis=1)
                level_max = max(level_max, blk_abs.max(initial=0.0))
                now = time.monotonic()
                if now - last_print >= 1.0:
                    db = 20 * np.log10(level_max / full_scale + 1e-12)
                    print(f"\r  peak level {db:7.1f} dBFS   events saved: {saved}   ",
                          end="", flush=True)
                    level_max, last_print = 0.0, now

                if trig is None:
                    seg = np.abs(buf[searched:, mon].astype(np.float64)).max(axis=1)
                    hits = np.flatnonzero(seg >= thr)
                    if hits.size:
                        trig = searched + int(hits[0])
                        trig_time = datetime.now()
                    else:
                        searched = len(buf)
                        if len(buf) > keep_n:          # drop old history
                            cut = len(buf) - keep_n
                            buf = buf[cut:]
                            searched -= cut
                        continue

                if len(buf) < trig + look_n + post_n:
                    continue                          # still collecting the POST part

                win = np.abs(buf[trig:trig + look_n, mon].astype(np.float64)).max(axis=1)
                peak = trig + int(np.argmax(win))
                lo, hi = peak - pre_n, peak + post_n
                if lo < 0:
                    print(f"\n  note: only {peak / fs:.3f}s available before the peak")
                    lo = 0
                clip = buf[lo:hi]
                peak_db = 20 * np.log10(np.abs(buf[peak, mon].astype(np.float64)).max()
                                        / full_scale + 1e-12)

                saved += 1
                if args.split:
                    name = "uma8_" + trig_time.strftime("%Y%m%d_%H%M%S_%f")[:-3] + ".wav"
                    path = os.path.join(args.o, name)
                    sf.write(path, clip, fs, subtype=wav_subtype)
                    print(f"\n  [{saved}] peak {peak_db:.1f} dBFS -> {path} "
                          f"({len(clip) / fs:.3f}s, {channels} ch)")
                else:
                    if saved > 1 and len(gap):
                        mux_wav.write(gap)
                        mux_pos += len(gap)
                    start = mux_pos
                    mux_wav.write(clip)
                    mux_wav.flush()               # keep the wav header valid after each event
                    mux_pos += len(clip)
                    peak_s = (start + peak - lo) / fs
                    mux_csv.write(f"{saved},{trig_time.isoformat(timespec='milliseconds')},"
                                  f"{start / fs:.6f},{mux_pos / fs:.6f},{peak_s:.6f},"
                                  f"{peak_db:.2f}\n")
                    mux_csv.flush()
                    print(f"\n  [{saved}] peak {peak_db:.1f} dBFS -> event at "
                          f"{peak_s:.3f}s in {os.path.basename(mux_path)}")

                buf = buf[hi:]                        # resume detection after the clip
                searched = 0
                trig = None
                wait_start = time.monotonic()         # reset the max-wait timer
    except KeyboardInterrupt:
        timed_out = False
    finally:
        if mux_wav is not None:
            mux_wav.close()
            mux_csv.close()
            if saved == 0:                        # don't leave empty files behind
                os.remove(mux_path)
                os.remove(stem + ".csv")

    print(f"\n\nStopped. {saved} event(s) saved.")
    if timed_out:
        print(f"No trigger for {args.max_wait:g}s (-max_wait), so recording ended.")
    if overflows[0]:
        print(f"Warning: {overflows[0]} input overflow(s) -- some audio may have been dropped.")


if __name__ == "__main__":
    main()
