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
  py sonifimeteor.py merge                         # one composite track: all stations on a
                                                    # shared real-world timeline, panned L/R by
                                                    # geographic (east-west) position, so the
                                                    # listener hears the arrival sweep across the
                                                    # network the way the reference meteor.mp4
                                                    # animation shows it visually
"""
import json
import os
import sys
import argparse

import numpy as np
from obspy import read
from scipy import signal
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

# 'merge' resamples every station onto one shared clock -- 100 Hz covers the large majority of
# traces in the reference file natively (only a couple of 200 Hz / 20 Hz stations need resampling).
MERGE_TARGET_RATE_HZ = 100.0
# Public FDSN/EIDA providers tried (in order) to resolve station lat/lon for stereo panning.
# Not every network here is publicly archived -- e.g. KB is a private deployment specific to this
# study -- so this is best-effort; unresolved stations are simply centered rather than failing.
MERGE_COORD_PROVIDERS = ["ORFEUS", "EIDA", "RESIF", "GFZ", "BGR", "RASPISHAKE"]
STATION_COORDS_CACHE = os.path.join(METEOR_INPUT_DIR, "station_coordinates_cache.json")

VALID_ACTIONS = {"sonify", "list", "merge"}


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


def _pick_representative_traces(traces):
    """One trace per station (prefer the vertical/Z component) so a 3-component station
    isn't triple-counted -- and so there's a single, unambiguous point to pan in the mix."""
    by_station = {}
    for tr in traces:
        by_station.setdefault((tr.stats.network, tr.stats.station), []).append(tr)
    reps = []
    for trs in by_station.values():
        vertical = [t for t in trs if t.stats.channel.upper().endswith("Z")]
        reps.append(vertical[0] if vertical else sorted(trs, key=lambda t: t.stats.channel)[0])
    return reps


def _fetch_station_coords(pairs):
    """Best-effort lat/lon lookup for (network, station) pairs via public FDSN/EIDA providers,
    cached to disk since the same ~33 stations are looked up on every run. Some networks here
    (e.g. KB, a private deployment for this study) simply have no public station metadata --
    that's expected, not an error, and those stations just end up centered instead of panned."""
    cache = {}
    if os.path.exists(STATION_COORDS_CACHE):
        with open(STATION_COORDS_CACHE, "r", encoding="utf-8") as f:
            raw = json.load(f)
        cache = {tuple(k.split(".", 1)): (tuple(v) if v else None) for k, v in raw.items()}

    missing_by_net = {}
    for net, sta in pairs:
        if (net, sta) not in cache:
            missing_by_net.setdefault(net, []).append(sta)

    if missing_by_net:
        from obspy.clients.fdsn import Client

        any_success = False
        for net, stas in missing_by_net.items():
            resolved = {}
            for provider in MERGE_COORD_PROVIDERS:
                try:
                    client = Client(provider, timeout=15)
                    inv = client.get_stations(network=net, station=",".join(sorted(set(stas))),
                                               level="station")
                    for network_obj in inv:
                        for sta_obj in network_obj:
                            resolved[sta_obj.code] = (sta_obj.latitude, sta_obj.longitude)
                    any_success = True
                    break
                except Exception:
                    continue
            for sta in stas:
                cache[(net, sta)] = resolved.get(sta)

        # Only persist "not found" for networks we actually managed to talk to a server about --
        # if every provider failed for every network, it's likely we're offline right now, and we
        # don't want to permanently cache false negatives just because of a transient outage.
        if any_success:
            serializable = {f"{n}.{s}": (list(v) if v else None) for (n, s), v in cache.items()}
            with open(STATION_COORDS_CACHE, "w", encoding="utf-8") as f:
                json.dump(serializable, f, indent=2, sort_keys=True)

    return cache


def do_merge(st, mseed_path, station_filter, speed_up_factor, freqmin, freqmax, taper_max_length,
             use_pan):
    """Time-align every selected station's trace on one real-world clock (using each trace's own
    absolute UTC start time, which already differs station-to-station by true travel time), and
    optionally pan each station left/right by geographic (east-west) position -- so playback
    recreates the wavefront sweeping across the network, the way meteor.mp4 shows it visually."""
    if station_filter and station_filter.strip().lower() != "all":
        traces = [tr for tr in st if station_filter.lower() in tr.id.lower()]
        if not traces:
            available = ", ".join(sorted(tr.id for tr in st))
            raise SystemExit(f"--station '{station_filter}' matched no trace. Available: {available}")
    else:
        traces = list(st)

    reps = _pick_representative_traces(traces)
    if len(reps) < 2:
        raise SystemExit(f"merge needs at least 2 stations; only matched {len(reps)}. "
                          "Try a broader --station filter or omit --station.")

    processed = []
    for tr in reps:
        try:
            processed.append(_process_trace(tr, freqmin, freqmax, taper_max_length))
        except Exception as exc:
            print(f"[warn] Skipping {tr.id} from merge: {exc}", file=sys.stderr)
    if len(processed) < 2:
        raise SystemExit("Fewer than 2 traces survived processing; cannot build a merge.")

    pans = {tr.id: 0.0 for tr in processed}
    if use_pan:
        pairs = sorted({(tr.stats.network, tr.stats.station) for tr in processed})
        print(f"[sonifimeteor] Looking up coordinates for {len(pairs)} station(s) "
              f"(cache: {os.path.relpath(STATION_COORDS_CACHE, BASE_DIR)})...")
        coords = _fetch_station_coords(pairs)
        resolved = {tr.id: coords.get((tr.stats.network, tr.stats.station)) for tr in processed}
        unresolved = sorted({f"{tr.stats.network}.{tr.stats.station}"
                              for tr in processed if not resolved[tr.id]})
        if unresolved:
            print(f"[warn] No public coordinates for {len(unresolved)} station(s) -- centered "
                  f"in the stereo field instead of panned: {', '.join(unresolved)}", file=sys.stderr)

        lons = [c[1] for c in resolved.values() if c]
        if len(lons) >= 2:
            lon_min, lon_max = min(lons), max(lons)
            span = (lon_max - lon_min) or 1.0
            for tr in processed:
                c = resolved[tr.id]
                if c:
                    pans[tr.id] = ((c[1] - lon_min) / span) * 2.0 - 1.0  # -1 west .. +1 east
        else:
            print("[warn] Not enough resolved coordinates to pan -- merge will be centered/mono.",
                  file=sys.stderr)
            use_pan = False

    for tr in processed:
        if abs(tr.stats.sampling_rate - MERGE_TARGET_RATE_HZ) > 1e-6:
            # obspy's own Trace.resample() calls the now-removed ndarray.newbyteorder() on
            # numpy>=2.0, so resample via scipy directly instead.
            n_new = int(round(tr.stats.npts * MERGE_TARGET_RATE_HZ / tr.stats.sampling_rate))
            tr.data = signal.resample(tr.data, n_new).astype(np.float64)
            tr.stats.sampling_rate = MERGE_TARGET_RATE_HZ

    ref_start = min(tr.stats.starttime for tr in processed)
    ref_end = max(tr.stats.endtime for tr in processed)
    n_samples = int(round((ref_end - ref_start) * MERGE_TARGET_RATE_HZ)) + 1

    global_peak = max(np.max(np.abs(tr.data)) for tr in processed if tr.data.size)
    if global_peak == 0:
        raise SystemExit("All selected traces are silent after processing; nothing to merge.")

    left = np.zeros(n_samples, dtype=np.float64)
    right = np.zeros(n_samples, dtype=np.float64)

    print(f"[sonifimeteor] Merging {len(processed)} station(s) onto one "
          f"{(n_samples / MERGE_TARGET_RATE_HZ):.1f}s real-time timeline (t=0 at {ref_start})"
          f"{' with geographic panning' if use_pan else ' (centered, no panning)'}:")
    for tr in sorted(processed, key=lambda t: t.stats.starttime):
        data = tr.data.astype(np.float64) / global_peak
        offset = max(0, int(round((tr.stats.starttime - ref_start) * MERGE_TARGET_RATE_HZ)))
        end = min(offset + data.size, n_samples)
        n = end - offset
        if n <= 0:
            continue
        pan = pans.get(tr.id, 0.0)
        theta = (pan + 1.0) * (np.pi / 4.0)  # equal-power panning: -1->left, 0->center, +1->right
        left[offset:end] += data[:n] * np.cos(theta)
        right[offset:end] += data[:n] * np.sin(theta)
        print(f"  {tr.id:20s} pan={pan:+.2f}  arrives at t={offset / MERGE_TARGET_RATE_HZ:6.1f}s")

    peak = max(np.max(np.abs(left)), np.max(np.abs(right)))
    if peak > 0:
        left /= peak
        right /= peak
    audio = (np.column_stack([left, right]) * 32767).astype(np.int16)

    wav_sample_rate = int(MERGE_TARGET_RATE_HZ * speed_up_factor)
    base_name = os.path.splitext(os.path.basename(mseed_path))[0]
    tag = "merged-panned" if use_pan else "merged-centered"
    wav_path = os.path.join(SONIFIMETEOR_DIR,
                             f"{base_name}_{tag}_{len(processed)}stations_{int(speed_up_factor)}x.wav")
    wavfile.write(wav_path, wav_sample_rate, audio)
    output_duration_s = n_samples / MERGE_TARGET_RATE_HZ / speed_up_factor
    print(f"[sonifimeteor] Wrote merged track: {n_samples / MERGE_TARGET_RATE_HZ:.1f}s real time "
          f"-> {output_duration_s:.2f}s at {speed_up_factor:.0f}x -> {wav_path}")
    return wav_path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Sonify a meteoroid-event seismic recording (multi-network, per-station "
                     "arrival-aligned cut).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "actions", nargs="*", default=None,
        help="'sonify' (default), 'list' (print available station/channel IDs), and/or 'merge' "
             "(one composite track with every station time-aligned and panned by geography)",
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
    parser.add_argument(
        "--no-pan", action="store_true",
        help="For 'merge': skip the FDSN/EIDA coordinate lookup and stereo panning, producing a "
             "centered (still time-aligned) track. Useful offline or if you don't want the "
             "network calls.",
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

    if "merge" in actions:
        do_merge(st, args.file, args.station, args.speed_up, args.freqmin, args.freqmax,
                 args.taper_max_length, use_pan=not args.no_pan)


if __name__ == "__main__":
    main()
