"""sonifimeteor.py -- sonify the multi-network meteoroid-event recording.

Copyright (C) 2026 Victor Mazon
Licensed under the GNU General Public License v3.0 (see LICENSE).

Data source: Dario Eickhoff (TU Dresden), shared 2026-09 via
https://datashare.tu-dresden.de/s/qpYLB3ZssrSjQ5W. The .mseed file has
instrument response already removed, and each trace is individually cut
tight around the same meteoroid's seismic arrival at that station
(~1 min before, ~4 min after) -- so absolute start times differ per
station (closer stations see the arrival sooner) and there's no single
"primary" channel the way there is for a DZA single-station fetch.

Runs standalone, independent of DZA01.py -- point it at any similarly-cut
meteor .mseed file with --file. Output goes to
datasets/sonifications_sonifimeteor/ (a separate folder from DZA01.py's
datasets/sonifications/, so a meteor run and a normal DZA01.py run can be
triggered side by side without colliding on filenames).

Examples
--------
  py sonifimeteor.py                              # sonify every trace in the default file
  py sonifimeteor.py list                         # print available station/channel IDs and exit
  py sonifimeteor.py --station DEP52 sonify        # only the KB.DEP52 station (3 components)
  py sonifimeteor.py --speed-up 40 sonify          # a punchier, more compressed 40x version
  py sonifimeteor.py --file other_event.mseed sonify
"""
import os
import sys
import argparse

import numpy as np
from obspy import read
from scipy.io import wavfile

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
METEOR_INPUT_DIR = os.path.join(BASE_DIR, "4_sonic")
DEFAULT_METEOR_MSEED = os.path.join(
    METEOR_INPUT_DIR, "4_sonic_koblenz_waveforms_response_removed.mseed"
)
DATASETS_DIR = os.path.join(BASE_DIR, "datasets")
SONIFIMETEOR_DIR = os.path.join(DATASETS_DIR, "sonifications_sonifimeteor")
os.makedirs(SONIFIMETEOR_DIR, exist_ok=True)

# The recording is only ~5-6 minutes per station, already cut tight around the
# arrival -- a much smaller speed-up than DZA01.py's 24h-recording default (200x)
# is enough to make it audible: 20x turns a 5 min clip into ~15s and keeps the
# impulsive arrival + coda clearly shaped, rather than compressing it into an
# unrecognizable blip.
DEFAULT_SPEED_UP = 20.0
# Matches the band used in check.py (Dario's own quick-look script) to make the
# arrival stand out against background noise.
DEFAULT_FREQMIN = 1.0
DEFAULT_FREQMAX = 10.0
# Percentage-based tapering eats a disproportionate share of a short clip (see
# DZA01.py's _process_and_save for the full story) -- these traces are only
# ~300s long, so the taper is capped to a small fixed length rather than a
# percentage, to avoid clicking at the edges without chewing into the coda.
DEFAULT_TAPER_MAX_LENGTH_S = 2.0

VALID_ACTIONS = {"sonify", "list"}


def load_stream(mseed_path):
    """Read a .mseed file into a Stream, with clear user-facing errors instead of a
    raw traceback, and drop any zero-length traces."""
    if not os.path.exists(mseed_path):
        raise SystemExit(f"File not found: {mseed_path}")
    try:
        st = read(mseed_path)
    except Exception as exc:
        raise SystemExit(f"Could not read {mseed_path} as MSEED data: {exc}") from exc

    empty_ids = [tr.id for tr in st if tr.stats.npts == 0]
    if empty_ids:
        print(f"[warn] Dropping {len(empty_ids)} zero-length trace(s): {', '.join(empty_ids)}",
              file=sys.stderr)
        st.traces = [tr for tr in st if tr.stats.npts > 0]
    if len(st) == 0:
        raise SystemExit(f"{mseed_path} contains no usable (non-empty) traces.")
    return st


def describe_meteor_stream(st):
    print(f"[sonifimeteor] {len(st)} trace(s) in this file:")
    for tr in sorted(st, key=lambda t: t.id):
        duration_s = tr.stats.npts / tr.stats.sampling_rate
        print(f"  - {tr.id:20s} {tr.stats.sampling_rate:6.1f} Hz  {duration_s:6.1f}s  "
              f"{tr.stats.starttime} - {tr.stats.endtime}")


def _process_trace(tr, freqmin, freqmax, taper_max_length):
    """Detrend, bandpass and lightly taper a copy of the trace. Response removal is
    NOT applied here -- the source file is already response-removed."""
    tr = tr.copy()
    tr.detrend("demean")
    tr.detrend("linear")
    tr.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4, zerophase=False)
    tr.taper(0.05, max_length=taper_max_length)
    return tr


def _sonify_one_trace(tr, out_dir, base_name, speed_up_factor):
    data = tr.data.astype(np.float64)
    data -= data.mean()
    peak = np.max(np.abs(data))
    if peak == 0:
        print(f"[warn] {tr.id}: trace is all-zero after processing; writing silent audio.",
              file=sys.stderr)
    else:
        data /= peak
    audio = (data * 32767).astype(np.int16)

    wav_sample_rate = int(tr.stats.sampling_rate * speed_up_factor)
    input_duration_s = tr.stats.npts / tr.stats.sampling_rate
    output_duration_s = input_duration_s / speed_up_factor

    tag = tr.id.replace(".", "-")
    wav_path = os.path.join(out_dir, f"{base_name}_{tag}_{int(speed_up_factor)}x.wav")
    wavfile.write(wav_path, wav_sample_rate, audio)
    print(f"[sonifimeteor]   {tr.id}: {input_duration_s:.1f}s -> {output_duration_s:.2f}s "
          f"at {speed_up_factor:.0f}x (wav sample rate {wav_sample_rate} Hz) -> {wav_path}")
    return wav_path


def do_sonify(st, mseed_path, station_filter, speed_up_factor, freqmin, freqmax, taper_max_length):
    if station_filter and station_filter.strip().lower() != "all":
        traces = [tr for tr in st if station_filter.lower() in tr.id.lower()]
        if not traces:
            available = ", ".join(sorted(tr.id for tr in st))
            raise SystemExit(f"--station '{station_filter}' matched no trace. Available: {available}")
    else:
        traces = list(st)

    base_name = os.path.splitext(os.path.basename(mseed_path))[0]
    print(f"[sonifimeteor] Sonifying {len(traces)}/{len(st)} trace(s) from "
          f"{os.path.basename(mseed_path)} -> {SONIFIMETEOR_DIR}")

    wav_paths = []
    failures = []
    for tr in traces:
        try:
            processed = _process_trace(tr, freqmin, freqmax, taper_max_length)
            wav_paths.append(_sonify_one_trace(processed, SONIFIMETEOR_DIR, base_name, speed_up_factor))
        except Exception as exc:
            print(f"[warn] Skipping {tr.id}: {exc}", file=sys.stderr)
            failures.append(tr.id)

    print(f"[sonifimeteor] Done: {len(wav_paths)}/{len(traces)} trace(s) sonified successfully.")
    if failures:
        print(f"[sonifimeteor] {len(failures)} trace(s) failed: {', '.join(failures)}", file=sys.stderr)
    return wav_paths


def parse_args():
    parser = argparse.ArgumentParser(
        description="Sonify a meteoroid-event seismic recording (multi-network, per-station "
                     "arrival-aligned cut).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "actions", nargs="*", default=None,
        help="'sonify' (default) and/or 'list' (print available station/channel IDs)",
    )
    parser.add_argument(
        "--file", default=DEFAULT_METEOR_MSEED,
        help=f"Path to the meteor .mseed file (default: {DEFAULT_METEOR_MSEED})",
    )
    parser.add_argument(
        "--station", default=None, metavar="SEED_ID_SUBSTRING",
        help="Only sonify traces whose id contains this substring (e.g. 'DEP52' or 'KB' or "
             "'HHZ'). Default: sonify every trace in the file -- this is a whole-network event "
             "capture, so 'all stations' is the natural default (unlike DZA01.py, which "
             "defaults to one primary channel).",
    )
    parser.add_argument(
        "--speed-up", type=float, default=DEFAULT_SPEED_UP,
        help=f"Playback speed multiplier (default: {DEFAULT_SPEED_UP:.0f}).",
    )
    parser.add_argument(
        "--freqmin", type=float, default=DEFAULT_FREQMIN, metavar="HZ",
        help=f"Lower bandpass corner in Hz (default: {DEFAULT_FREQMIN:g}).",
    )
    parser.add_argument(
        "--freqmax", type=float, default=DEFAULT_FREQMAX, metavar="HZ",
        help=f"Upper bandpass corner in Hz (default: {DEFAULT_FREQMAX:g}).",
    )
    parser.add_argument(
        "--taper-max-length", type=float, default=DEFAULT_TAPER_MAX_LENGTH_S, metavar="SECONDS",
        help=f"Max taper length per side, in seconds (default: {DEFAULT_TAPER_MAX_LENGTH_S:g}).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    actions = {a.lower() for a in args.actions} if args.actions else {"sonify"}
    invalid = actions - VALID_ACTIONS
    if invalid:
        raise SystemExit(
            f"Invalid action(s): {', '.join(sorted(invalid))}. "
            f"Choose from: {', '.join(sorted(VALID_ACTIONS))}"
        )
    if args.speed_up <= 0:
        raise SystemExit(f"--speed-up must be > 0 (got {args.speed_up})")
    if args.freqmin <= 0:
        raise SystemExit(f"--freqmin must be > 0 (got {args.freqmin})")
    if args.freqmax <= args.freqmin:
        raise SystemExit(f"--freqmax ({args.freqmax}) must be > --freqmin ({args.freqmin})")
    if args.taper_max_length < 0:
        raise SystemExit(f"--taper-max-length must be >= 0 (got {args.taper_max_length})")

    st = load_stream(args.file)

    if "list" in actions:
        describe_meteor_stream(st)
        if actions == {"list"}:
            return

    if "sonify" in actions:
        do_sonify(st, args.file, args.station, args.speed_up, args.freqmin, args.freqmax,
                  args.taper_max_length)


if __name__ == "__main__":
    main()
