import sounddevice as sd
import soundfile as sf

CHANNELS = 6
SAMPLE_RATE = 48000
DURATION = 10          # seconds
OUTPUT = "sample.wav"


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
print(f"Recording {DURATION}s, {CHANNELS} ch @ {SAMPLE_RATE} Hz...")

audio = sd.rec(int(DURATION * SAMPLE_RATE), samplerate=SAMPLE_RATE,
               channels=CHANNELS, dtype="int16", device=device)
sd.wait()

sf.write(OUTPUT, audio, SAMPLE_RATE, subtype="PCM_16")
print(f"Saved {OUTPUT} with shape {audio.shape} (samples, channels)")
