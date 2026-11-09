import argparse

import sounddevice as sd
import soundfile as sf

CHANNELS = 6           # minimum channels the device must expose
SAMPLE_RATE = 48000
DURATION = 4          # default record duration (s), override with -t
OUTPUT = "sample.wav"


def positive_float(value):
    t = float(value)
    if t <= 0:
        raise argparse.ArgumentTypeError("duration must be > 0")
    return t


parser = argparse.ArgumentParser(
    description="Record all UMA-8 channels to a WAV file.", add_help=False)
parser.add_argument("-h", "-help", "--help", action="help",
                    help="show this help message and exit")
parser.add_argument("-t", type=positive_float, default=DURATION, metavar="SEC",
                    help=f"record duration in seconds (default: {DURATION})")
parser.add_argument("-bits", type=int, choices=(16, 24), default=24,
                    help="WAV sample resolution in bits (default: 24)")
args = parser.parse_args()
DURATION = args.t

# sounddevice can't capture packed 24-bit into numpy, so 24-bit is captured
# as int32 (top 24 bits used) and written as PCM_24.
REC_DTYPE, WAV_SUBTYPE = {16: ("int16", "PCM_16"),
                          24: ("int32", "PCM_24")}[args.bits]


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


device, info = find_uma8(CHANNELS)
if device is None:
    print("No UMA-8 input with 6+ channels found. Available devices:\n")
    print(sd.query_devices())
    raise SystemExit("Check the firmware and the Windows recording format.")

print(f"Using [{device}] {info['name']} "
      f"({sd.query_hostapis(info['hostapi'])['name']}, "
      f"{info['max_input_channels']} ch max)")

# WASAPI shared mode only accepts the device's own channel count and rate,
# so open the stream with exactly what Windows reports (8 ch on RAW firmware).
# All channels are saved; the GSL pipeline picks mics via channel_map (0,1,3,5).
rec_channels = info["max_input_channels"]
if "WASAPI" in sd.query_hostapis(info["hostapi"])["name"]:
    SAMPLE_RATE = int(info["default_samplerate"])

print(f"Recording {DURATION}s, {rec_channels} ch @ {SAMPLE_RATE} Hz, {args.bits}-bit...")

audio = sd.rec(int(DURATION * SAMPLE_RATE), samplerate=SAMPLE_RATE,
               channels=rec_channels, dtype=REC_DTYPE, device=device)
sd.wait()

sf.write(OUTPUT, audio, SAMPLE_RATE, subtype=WAV_SUBTYPE)
print(f"Saved {OUTPUT} with shape {audio.shape} (samples, channels)")
