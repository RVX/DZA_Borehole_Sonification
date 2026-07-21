"""DZA01 seismic data toolkit — fetch, plot, and sonify in one place.

Copyright (C) 2026 Victor Mazon
Licensed under the GNU General Public License v3.0 (see LICENSE).

A single entry point for the DZA01 (KIT / KB network) seismic station:

  1. FETCH   -> download + process a waveform window, save as .mseed
  2. PLOT    -> render a styled spectrogram + waveform plot (.png)
  3. SONIFY  -> speed the waveform up into an audible .wav
  4. PLAY    -> play a .wav straight from the terminal

Each run handles a single listening piece (see --duration). Use --listen-minutes
to target an exact sonification length: the script calculates the raw recording
length needed (when fetching) or the speed-up factor needed (when reusing an
existing file) to hit that length precisely.

Folder layout (created automatically next to this script):
  datasets/plot/               waveform + spectrogram plots (.png)
  datasets/mseed/               raw processed waveform data (.mseed)
  datasets/sonifications/       sonified audio (.wav)

Examples
--------
  py DZA01.py                                      # fetch + plot + sonify (15 min piece)
  py DZA01.py fetch plot                           # only fetch new data and plot it
  py DZA01.py sonify play --speed-up 300           # sonify latest file and play it
  py DZA01.py --list                               # list saved .mseed files
  py DZA01.py --pick 2 plot sonify play            # use file #2 from --list
  py DZA01.py --latency 2 --max-latency 60         # near real-time, backs off up to 1h if needed
  py DZA01.py --listen-minutes 15 --speed-up 20    # fetch exactly enough data (5h) for a 15 min listen
  py DZA01.py --listen-minutes 15                  # same, but at the 100x audible floor (fetches 25h)
  py DZA01.py --pick 0 sonify --listen-minutes 15  # re-stretch an existing file to a 15 min listen
"""
import os
import sys
import glob
import shutil
import argparse
import subprocess

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.signal import spectrogram
from obspy.clients.fdsn import Client
from obspy import UTCDateTime, read
from scipy.io import wavfile

# ---------------------------------------------------------------------------
# Station / network configuration
# ---------------------------------------------------------------------------
NETWORK = "KB"
STATION = "DZA1*"
CHANNEL = "HH*"
FDSN_URL = "http://ws.gpi.kit.edu"

# ---------------------------------------------------------------------------
# Folder layout
# ---------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATASETS_DIR = os.path.join(BASE_DIR, "datasets")
PLOT_DIR = os.path.join(DATASETS_DIR, "plot")
MSEED_DIR = os.path.join(DATASETS_DIR, "mseed")
SONIFY_DIR = os.path.join(DATASETS_DIR, "sonifications")

for _dir in (PLOT_DIR, MSEED_DIR, SONIFY_DIR):
    os.makedirs(_dir, exist_ok=True)

# pyTREMOR-inspired dark theme
BG_COLOR = "#0b0c10"
FG_COLOR = "0.85"
GRID_COLOR = "0.3"
COMPONENT_COLORS = {
    "Z": "#00e5ff",  # cyan
    "N": "#ff9100",  # orange
    "1": "#ff9100",
    "E": "#ff1744",  # red/pink
    "2": "#ff1744",
}

DEFAULT_SPEED_UP = 200   # audio playback speed multiplier
DEFAULT_DURATION = 15.0  # default listening piece length, minutes

# Below this speed-up, seismic ground motion (bandpassed to 0.5-10 Hz) is
# barely audible as anything more than a faint rumble. Used as the implicit
# floor for --listen-minutes when the user hasn't set --speed-up explicitly.
MIN_LISTEN_SPEED_UP = 100.0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(
        description="Fetch, plot and/or sonify DZA01 seismic waveform data.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "actions", nargs="*",
        default=None,
        help="What to do: any combination of fetch, plot, sonify, play, all "
             "(default: fetch plot sonify)",
    )

    fetch_group = parser.add_argument_group("fetch options")
    fetch_group.add_argument(
        "--duration", type=float, default=DEFAULT_DURATION,
        help=f"Length of the data window to fetch, in minutes (default: {DEFAULT_DURATION:.0f})",
    )
    fetch_group.add_argument(
        "--latency", type=float, default=2.0,
        help="How many minutes back from now to end the window (smaller = closer "
             "to real time, but the server may not have the data yet). Default: 2",
    )
    fetch_group.add_argument(
        "--max-latency", type=float, default=60.0,
        help="Ceiling (minutes) the script will back off to if --latency is too "
             "aggressive and no data is available yet. FDSN data lag can vary "
             "(e.g. around 10 min, more in CET winter than CEST summer), so this "
             "is kept generous by default. Default: 60",
    )
    fetch_group.add_argument(
        "--latency-step", type=float, default=2.0,
        help="How many minutes to add to the latency each retry when backing off. "
             "Default: 2",
    )

    sonify_group = parser.add_argument_group("sonify options")
    sonify_group.add_argument(
        "--speed-up", type=float, default=None,
        help=f"Playback speed multiplier for sonification (default: {DEFAULT_SPEED_UP}, or "
             f"{MIN_LISTEN_SPEED_UP:.0f}x if --listen-minutes is set and no --speed-up is "
             f"given, since ground motion is barely audible below that). A 15 min recording "
             f"at 200x becomes ~4.5s of audio -- lower this for longer audio (e.g. 20x turns "
             f"15 min into ~45s).",
    )
    sonify_group.add_argument(
        "--channel", default=None, metavar="SEED_ID_SUBSTRING",
        help="Which trace to sonify when multiple stations/channels are present, matched "
             "as a substring against the trace id (e.g. 'DZA11' or 'HHZ' or 'DZA13.00.HH1'). "
             "Default: the first trace in the file.",
    )
    sonify_group.add_argument(
        "--listen-minutes", type=float, default=None, metavar="MINUTES",
        help="Target length of the sonified audio, in minutes. If fetching new data, this "
             "calculates and overrides --duration (raw minutes needed = MINUTES x --speed-up, "
             f"or MINUTES x {MIN_LISTEN_SPEED_UP:.0f} if --speed-up isn't given, since that's "
             "the floor for ground motion to be clearly audible). If reusing an existing file, "
             "this instead overrides --speed-up to stretch/compress that recording to exactly "
             "this listening length (warns if the result falls below the audible floor).",
    )

    select_group = parser.add_argument_group("file selection (skip fetching)")
    select_group.add_argument(
        "--list", action="store_true",
        help="List saved .mseed files (with index numbers) and exit",
    )
    select_group.add_argument(
        "--pick", type=int, default=None, metavar="INDEX",
        help="Use the .mseed file at INDEX from --list instead of fetching new data",
    )
    select_group.add_argument(
        "--file", default=None,
        help="Use a specific .mseed file path instead of fetching new data",
    )

    return parser.parse_args()


VALID_ACTIONS = {"fetch", "plot", "sonify", "play", "all"}


# ---------------------------------------------------------------------------
# Fetch
# ---------------------------------------------------------------------------
def fetch_once(client, latency_min, duration_min):
    """Request inventory + waveforms for a given latency/duration (both in minutes)."""
    endtime = UTCDateTime() - latency_min * 60
    starttime = endtime - duration_min * 60
    inv = client.get_stations(
        network=NETWORK, station=STATION, starttime=starttime, endtime=endtime,
        level="response",
    )
    st = client.get_waveforms(NETWORK, STATION, "*", CHANNEL, starttime, endtime)
    if len(st) == 0:
        raise ValueError("no traces returned (data not available yet)")
    return st, inv, starttime, endtime


def describe_stream(st, prefix="info"):
    """Print a short summary of every trace in a Stream (station/channel/duration)."""
    print(f"[{prefix}] {len(st)} trace(s) in this file "
          f"(station pattern '{STATION}' can match multiple sub-stations/channels):")
    for tr in st:
        duration_s = tr.stats.npts / tr.stats.sampling_rate
        print(f"  - {tr.id}  |  {tr.stats.sampling_rate:.1f} Hz  |  {duration_s:.1f}s  |  {tr.stats.npts} samples")


def fetch_closest_to_real_time(client, latency_min, duration_min, max_latency_min, latency_step_min):
    """Try to get data as close to 'now' as possible, backing off only if needed."""
    latency = latency_min
    while True:
        try:
            st, inv, starttime, endtime = fetch_once(client, latency, duration_min)
            print(f"[fetch] Got data with {latency:.1f} min latency (duration {duration_min:.1f} min)")
            return st, inv, starttime, endtime
        except Exception as exc:
            if latency >= max_latency_min:
                raise RuntimeError(
                    f"Could not fetch data even at max latency of {max_latency_min:.1f} min: {exc}"
                ) from exc
            next_latency = min(latency + latency_step_min, max_latency_min)
            print(f"[fetch] No data at {latency:.1f} min latency ({exc}); retrying at {next_latency:.1f} min...")
            latency = next_latency


def do_fetch(args):
    """Fetch, process and save one waveform window. Returns the saved .mseed path."""
    client = Client(FDSN_URL)
    st, inv, starttime, endtime = fetch_closest_to_real_time(
        client, args.latency, args.duration, args.max_latency, args.latency_step
    )

    st.remove_response(inventory=inv)
    st.detrend("demean")
    st.detrend("linear")
    st.filter("bandpass", freqmin=0.5, freqmax=10.0, corners=4)
    st.taper(0.125)

    timestamp = starttime.strftime("%Y-%m-%d-%H-%M")
    mseed_path = os.path.join(MSEED_DIR, f"DZA1_{timestamp}.mseed")
    st.write(mseed_path, format="MSEED")
    print(f"[fetch] Saved waveform data to {mseed_path}")
    describe_stream(st, prefix="fetch")
    return mseed_path


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def component_color(channel_code):
    comp = channel_code[-1].upper() if channel_code else "Z"
    return COMPONENT_COLORS.get(comp, "#76ff03")


def rms_envelope(data, sampling_rate, window_s=2.0):
    """Rolling RMS envelope, used as a translucent fill behind each waveform."""
    window = max(1, int(window_s * sampling_rate))
    squared = data.astype(np.float64) ** 2
    kernel = np.ones(window) / window
    return np.sqrt(np.convolve(squared, kernel, mode="same"))


def pick_primary_for_spectrogram(st):
    """Prefer a vertical (Z) component for the top spectrogram panel."""
    for tr in st:
        if tr.stats.channel.upper().endswith("Z"):
            return tr
    return st[0]


def style_axes(ax):
    ax.set_facecolor(BG_COLOR)
    ax.tick_params(colors=FG_COLOR, labelsize=8)
    for spine in ax.spines.values():
        spine.set_color(GRID_COLOR)


def do_plot(mseed_path):
    st = read(mseed_path)
    base = os.path.splitext(os.path.basename(mseed_path))[0]
    plot_path = os.path.join(PLOT_DIR, f"{base}.png")

    n_traces = len(st)
    primary = pick_primary_for_spectrogram(st)
    sr = primary.stats.sampling_rate

    fig = plt.figure(figsize=(12, 3 + 1.3 * n_traces), facecolor=BG_COLOR)
    gs = fig.add_gridspec(1 + n_traces, 1, height_ratios=[2.2] + [1] * n_traces, hspace=0.15)

    # -- top panel: spectrogram of the primary (vertical) channel --
    spec_ax = fig.add_subplot(gs[0])
    nperseg = int(max(32, min(sr * 4, len(primary.data))))
    noverlap = int(nperseg * 0.9)
    f, t, Sxx = spectrogram(primary.data.astype(np.float64), fs=sr, nperseg=nperseg, noverlap=noverlap)
    Sxx_db = 10 * np.log10(Sxx + 1e-30)
    vmin, vmax = np.percentile(Sxx_db, [5, 99.5])
    mesh = spec_ax.pcolormesh(t, f, Sxx_db, cmap="inferno", shading="auto", vmin=vmin, vmax=vmax)
    spec_ax.set_ylim(0, min(12, sr / 2))
    spec_ax.set_ylabel("Frequency (Hz)", color=FG_COLOR)
    spec_ax.set_title(f"Spectrogram - {primary.id}", color="white", fontsize=10, loc="left")
    spec_ax.set_xticklabels([])
    style_axes(spec_ax)
    cbar = fig.colorbar(mesh, ax=spec_ax, pad=0.01, fraction=0.02)
    cbar.set_label("Power (dB)", color=FG_COLOR)
    cbar.ax.yaxis.set_tick_params(color=FG_COLOR)
    plt.setp(cbar.ax.get_yticklabels(), color=FG_COLOR)

    # -- one panel per trace: waveform + RMS envelope --
    max_points = 200_000  # cap rendered points per trace for speed/memory on long recordings
    for i, tr in enumerate(st):
        ax = fig.add_subplot(gs[i + 1])
        times = tr.times()
        data = tr.data
        step = max(1, len(data) // max_points)
        if step > 1:
            times = times[::step]
            data = data[::step]
        color = component_color(tr.stats.channel)
        env = rms_envelope(data, tr.stats.sampling_rate / step)
        ax.fill_between(times, -env, env, color=color, alpha=0.25, linewidth=0)
        ax.plot(times, data, color=color, linewidth=0.5)
        ax.set_xlim(times[0], times[-1])
        ax.text(
            0.005, 0.85, tr.id, transform=ax.transAxes, color="white", fontsize=8,
            family="monospace", va="top",
            bbox=dict(boxstyle="round", facecolor="#1c1f26", edgecolor=GRID_COLOR, alpha=0.85),
        )
        style_axes(ax)
        if i < n_traces - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("Time (s)", color=FG_COLOR)

    title = f"{NETWORK}.{STATION}  |  {primary.stats.starttime} - {primary.stats.endtime} UTC"
    fig.suptitle(title, color="white", fontsize=11, y=0.995)
    fig.subplots_adjust(left=0.07, right=0.93, top=0.94, bottom=0.05, hspace=0.15)
    fig.savefig(plot_path, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)

    print(f"[plot] Saved plot to {plot_path} ({n_traces} traces shown)")
    return plot_path


# ---------------------------------------------------------------------------
# Sonify
# ---------------------------------------------------------------------------
def pick_trace(st, channel_filter=None):
    """Pick one trace from a Stream, warning about the others being ignored."""
    if channel_filter:
        matches = [tr for tr in st if channel_filter.lower() in tr.id.lower()]
        if not matches:
            available = ", ".join(tr.id for tr in st)
            raise SystemExit(
                f"--channel '{channel_filter}' matched no trace. Available: {available}"
            )
        tr = matches[0]
    else:
        tr = st[0]

    if len(st) > 1:
        others = [t.id for t in st if t.id != tr.id]
        print(f"[sonify] {len(st)} channels available; using '{tr.id}'. "
              f"Ignoring: {', '.join(others)}. Use --channel to pick a different one.")
    return tr


def do_sonify(mseed_path, speed_up_factor=DEFAULT_SPEED_UP, channel_filter=None):
    st = read(mseed_path)
    tr = pick_trace(st, channel_filter)

    # normalize to -1..1 then scale to 16-bit PCM range
    data = tr.data.astype(np.float64)
    data -= data.mean()
    peak = np.max(np.abs(data))
    if peak > 0:
        data /= peak
    audio = (data * 32767).astype(np.int16)

    # speeding up = play the same samples at a higher sample rate
    wav_sample_rate = int(tr.stats.sampling_rate * speed_up_factor)
    input_duration_s = tr.stats.npts / tr.stats.sampling_rate
    output_duration_s = input_duration_s / speed_up_factor

    base = os.path.splitext(os.path.basename(mseed_path))[0]
    wav_path = os.path.join(SONIFY_DIR, f"{base}_{int(speed_up_factor)}x.wav")
    wavfile.write(wav_path, wav_sample_rate, audio)
    print(f"[sonify] Saved audio to {wav_path}")
    print(f"[sonify]   trace: {tr.id}  |  input: {input_duration_s:.1f}s  ->  "
          f"output: {output_duration_s:.1f}s at {speed_up_factor:.0f}x speed "
          f"(wav sample rate {wav_sample_rate} Hz)")
    return wav_path


# ---------------------------------------------------------------------------
# Play
# ---------------------------------------------------------------------------
def do_play(wav_path):
    print(f"[play] Playing {wav_path}")
    if sys.platform.startswith("win"):
        import winsound
        winsound.PlaySound(wav_path, winsound.SND_FILENAME)
    elif sys.platform == "darwin":
        subprocess.run(["afplay", wav_path], check=False)
    else:
        for player in ("paplay", "aplay", "ffplay"):
            if shutil.which(player):
                cmd = [player, wav_path] if player != "ffplay" else [player, "-nodisp", "-autoexit", wav_path]
                subprocess.run(cmd, check=False)
                return
        print(f"[play] No audio player found. Open manually: {wav_path}")


# ---------------------------------------------------------------------------
# File selection helpers
# ---------------------------------------------------------------------------
def list_mseed_files():
    return sorted(glob.glob(os.path.join(MSEED_DIR, "*.mseed")))


def print_file_list():
    files = list_mseed_files()
    if not files:
        print(f"No .mseed files found in {MSEED_DIR}.")
        return
    print("Saved .mseed files:")
    for i, f in enumerate(files):
        print(f"  [{i}] {os.path.basename(f)}")


def find_wav_for(mseed_path, speed_up_factor):
    base = os.path.splitext(os.path.basename(mseed_path))[0]
    candidate = os.path.join(SONIFY_DIR, f"{base}_{int(speed_up_factor)}x.wav")
    if os.path.exists(candidate):
        return candidate
    # fall back to any sonification of this same recording, any speed
    matches = sorted(glob.glob(os.path.join(SONIFY_DIR, f"{base}_*x.wav")))
    return matches[-1] if matches else None


def stream_duration_seconds(st):
    tr = st[0]
    return tr.stats.npts / tr.stats.sampling_rate


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    args = parse_args()

    if args.actions:
        invalid = set(args.actions) - VALID_ACTIONS
        if invalid:
            raise SystemExit(
                f"Invalid action(s): {', '.join(sorted(invalid))}. "
                f"Choose from: {', '.join(sorted(VALID_ACTIONS))}"
            )

    if args.duration <= 0:
        raise SystemExit(f"--duration must be > 0 (got {args.duration})")
    if args.latency < 0:
        raise SystemExit(f"--latency must be >= 0 (got {args.latency})")
    if args.max_latency < args.latency:
        raise SystemExit("--max-latency must be >= --latency")
    if args.speed_up is not None and args.speed_up <= 0:
        raise SystemExit(f"--speed-up must be > 0 (got {args.speed_up})")
    if args.listen_minutes is not None and args.listen_minutes <= 0:
        raise SystemExit(f"--listen-minutes must be > 0 (got {args.listen_minutes})")

    if args.list:
        print_file_list()
        return

    actions = set(args.actions) if args.actions else {"fetch", "plot", "sonify"}
    if "all" in actions:
        actions = {"fetch", "plot", "sonify"}

    # resolve which .mseed file to use
    mseed_path = None
    if args.file:
        mseed_path = args.file
        actions.discard("fetch")
    elif args.pick is not None:
        files = list_mseed_files()
        if not files or not (0 <= args.pick < len(files)):
            print_file_list()
            raise SystemExit(f"Invalid --pick index: {args.pick}")
        mseed_path = files[args.pick]
        actions.discard("fetch")

    fetch_needed = "fetch" in actions
    speed_up_explicit = args.speed_up is not None

    if args.listen_minutes and fetch_needed:
        # Without an explicit --speed-up, use MIN_LISTEN_SPEED_UP (100x) as a
        # sensible floor rather than silently multiplying by the plain
        # DEFAULT_SPEED_UP (200x). Below ~100x, seismic ground motion is
        # barely audible as more than a faint rumble, so a "listen" piece
        # needs at least that much compression to actually be worth hearing.
        # Pass --speed-up explicitly to request a different amount of raw
        # data on purpose (e.g. --listen-minutes 15 --speed-up 20 -> 5h fetch).
        effective_speed_up = args.speed_up if speed_up_explicit else MIN_LISTEN_SPEED_UP
        required_duration = args.listen_minutes * effective_speed_up
        hours = required_duration / 60
        print(
            f"[info] To get a {args.listen_minutes:.1f} min sonification at "
            f"{effective_speed_up:.0f}x speed-up, fetching {required_duration:.1f} min "
            f"({hours:.2f} h) of raw data (overrides --duration)."
        )
        if hours > 6:
            print(
                f"[info] That's a large request ({hours:.1f} h). Lower --speed-up "
                f"to fetch less data, if desired."
            )
        args.duration = required_duration

    if fetch_needed:
        mseed_path = do_fetch(args)
    elif mseed_path is None:
        files = list_mseed_files()
        if not files:
            raise SystemExit(
                "No .mseed files available and 'fetch' not requested. "
                "Run with 'fetch' or use --file/--pick."
            )
        mseed_path = files[-1]
        print(f"[info] Using latest saved file: {mseed_path}")

    if ("plot" in actions or "sonify" in actions or "play" in actions) and not fetch_needed:
        describe_stream(read(mseed_path), prefix="info")

    if args.listen_minutes and not fetch_needed:
        input_duration_s = stream_duration_seconds(read(mseed_path))
        target_s = args.listen_minutes * 60
        args.speed_up = input_duration_s / target_s
        print(
            f"[info] Existing recording is {input_duration_s / 60:.1f} min long; using "
            f"speed-up {args.speed_up:.2f}x to produce a {args.listen_minutes:.1f} min "
            f"sonification (overrides --speed-up)."
        )
        if args.speed_up < MIN_LISTEN_SPEED_UP:
            needed_min = target_s * MIN_LISTEN_SPEED_UP / 60
            print(
                f"[warn] {args.speed_up:.2f}x is below the ~{MIN_LISTEN_SPEED_UP:.0f}x "
                f"floor for audible ground motion; this piece will likely sound like a "
                f"faint rumble. Fetch a longer recording (>= {needed_min:.0f} min) to "
                f"reach {MIN_LISTEN_SPEED_UP:.0f}x for the same {args.listen_minutes:.1f} "
                f"min listen length."
            )
    elif args.listen_minutes and fetch_needed:
        # We just fetched exactly enough data for this listen length; speed-up
        # is whatever we assumed above (MIN_LISTEN_SPEED_UP unless the user
        # set --speed-up explicitly).
        args.speed_up = args.speed_up if speed_up_explicit else MIN_LISTEN_SPEED_UP
    elif args.speed_up is None:
        args.speed_up = DEFAULT_SPEED_UP

    if "plot" in actions:
        do_plot(mseed_path)

    wav_path = None
    if "sonify" in actions:
        wav_path = do_sonify(mseed_path, args.speed_up, args.channel)

    if "play" in actions:
        if wav_path is None:
            wav_path = find_wav_for(mseed_path, args.speed_up)
        if wav_path is None:
            wav_path = do_sonify(mseed_path, args.speed_up, args.channel)
        do_play(wav_path)


if __name__ == "__main__":
    main()
