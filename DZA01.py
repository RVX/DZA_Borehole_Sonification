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
  datasets/plot/                waveform + spectrogram plots (.png)
  datasets/mseed/               raw processed waveform data (.mseed) + fetch metadata (.json)
  datasets/sonifications/       sonified audio (.wav)
  datasets/maps/                station location maps (.png)

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
  py DZA01.py --hours-back 24 --speed-up 100       # last 24h of data, sonified at 100x (pytremor-style)
  py DZA01.py --pick 0 sonify --channel all        # sonify every channel in the file, one .wav each
  py DZA01.py --days-back 10 --speed-up 100        # last 10 days as 10 daily files, all plotted+sonified
  py DZA01.py map                                  # print station coordinates/depths and save a map plot
  py DZA01.py --sites 1,3 fetch plot sonify        # fetch site 1 (DZA11/13) and site 3 (DZA31/33) together
  py DZA01.py --freqmin 0.0001 --freqmax 0.01 --speed-up 20000 --days-back 10 fetch sonify
                                                    # ultra-low-frequency band (<10 mHz), see README
"""
import os
import sys
import glob
import json
import math
import colorsys
import shutil
import argparse
import subprocess
from datetime import datetime, timedelta

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.signal import spectrogram
from obspy.clients.fdsn import Client
from obspy import UTCDateTime, read
from scipy.io import wavfile

# ---------------------------------------------------------------------------
# Station / network configuration
# ---------------------------------------------------------------------------
NETWORK = "KB"
CHANNEL = "HH*"
FDSN_URL = "http://ws.gpi.kit.edu"

# DZA station naming: DZA<site><1|3>, e.g. DZA11, DZA13, DZA31, DZA33...
#   suffix 1 -> surface sensor (0 m depth)
#   suffix 3 -> borehole sensor (~240 m depth, +/- a few meters, per station metadata)
# Confirmed by the station scientist (KIT/GPI), July 2026. Additional site numbers may
# be added later -- if so, add them to STATION_SITES below rather than hardcoding.
BOREHOLE_DEPTH_M = 240.0

# lat/lon/elevation are hardcoded (rather than always queried live) so plotting a
# station-location context map never needs a network call. Confirmed via `map`
# action against the live FDSN service, July 2026. `None` means not yet installed/
# no confirmed coordinates.
STATION_SITES = {
    1: dict(
        surface="DZA11", borehole="DZA13", status="active",
        lat=51.32358, lon=14.24694, elevation_m=134.6,
        note="Streaming continuously since project start.",
    ),
    3: dict(
        surface="DZA31", borehole="DZA33", status="active",
        lat=51.24889, lon=14.15348, elevation_m=211.5,
        note="Streaming, but may be less stable than site 1: solar-powered in a "
             "forest region with limited light.",
    ),
    4: dict(
        surface="DZA41", borehole="DZA43", status="planned",
        lat=None, lon=None, elevation_m=None,
        note="Installation planned for ~mid-to-end September 2026, pending drilling; "
             "no data expected before then.",
    ),
    6: dict(
        surface="DZA61", borehole="DZA63", status="planned",
        lat=None, lon=None, elevation_m=None,
        note="Installation planned for ~mid-to-end September 2026, pending drilling; "
             "no data expected before then.",
    ),
}
ACTIVE_SITE_NUMBERS = sorted(n for n, info in STATION_SITES.items() if info["status"] == "active")

# Broadband sensor models per surface station and their eigenperiod (the number in the
# model name). DZA11 uses a shorter-eigenperiod sensor than every other surface site;
# see the "Sensor differences" section in README.md for what that implies.
DEFAULT_SURFACE_SENSOR = "Trillium Horizon 120s"
DEFAULT_EIGENPERIOD_S = 120.0
SURFACE_SENSOR_OVERRIDES = {"DZA11": "Trillium Compact 20s"}
EIGENPERIOD_S_OVERRIDES = {"DZA11": 20.0}


def sensor_for_station(station_code):
    """Return (model_name, eigenperiod_s) for a surface station code."""
    model = SURFACE_SENSOR_OVERRIDES.get(station_code, DEFAULT_SURFACE_SENSOR)
    eigenperiod = EIGENPERIOD_S_OVERRIDES.get(station_code, DEFAULT_EIGENPERIOD_S)
    return model, eigenperiod


def _build_station_info_by_code():
    """Precompute a station-code -> metadata lookup (site number, depth, sensor,
    status) for every known station, so plot/map code never has to re-derive this
    from the station code string itself."""
    lookup = {}
    for site_num, info in STATION_SITES.items():
        sensor_model, eigenperiod = sensor_for_station(info["surface"])
        lookup[info["surface"]] = dict(
            site=site_num, depth_label="surface", depth_m=0.0,
            sensor_model=sensor_model, eigenperiod_s=eigenperiod,
            status=info["status"], note=info["note"],
        )
        lookup[info["borehole"]] = dict(
            site=site_num, depth_label="borehole", depth_m=BOREHOLE_DEPTH_M,
            sensor_model=None, eigenperiod_s=None,
            status=info["status"], note=info["note"],
        )
    return lookup


STATION_INFO_BY_CODE = _build_station_info_by_code()


def build_station_pattern(sites):
    """Turn a list of site numbers (e.g. [1, 3]) into an FDSN comma-separated station
    code list (e.g. 'DZA11,DZA13,DZA31,DZA33'), covering both the surface and
    borehole sensor at each site."""
    codes = []
    for site in sites:
        info = STATION_SITES.get(site)
        if info is None:
            raise SystemExit(
                f"Unknown site number {site}. Known sites: {sorted(STATION_SITES)} "
                f"(active: {ACTIVE_SITE_NUMBERS})"
            )
        codes.append(info["surface"])
        codes.append(info["borehole"])
    return ",".join(codes)


def build_file_prefix(sites):
    """Filename prefix for saved .mseed/.png/.wav files. Kept as the historical 'DZA1'
    for the original single-site (site 1) default, and made explicit for any other
    site selection so multi-site fetches never collide/overwrite each other."""
    if sites == [1]:
        return "DZA1"
    return "DZA_sites-" + "-".join(str(s) for s in sites)


def parse_sites_arg(value, fallback):
    """Parse the --sites CLI value ('all', 'active', or comma-separated site numbers)
    into a sorted list of site numbers. `fallback` is used when --sites wasn't given."""
    if value is None:
        return fallback
    value = value.strip().lower()
    if value == "all":
        return sorted(STATION_SITES)
    if value == "active":
        return list(ACTIVE_SITE_NUMBERS)
    try:
        return sorted({int(v.strip()) for v in value.split(",") if v.strip()})
    except ValueError:
        raise SystemExit(
            f"--sites must be 'all', 'active', or comma-separated site numbers (got {value!r})"
        )


def print_site_advisories(sites):
    """Print any known operational caveats for the selected sites (stability, install
    status) so they're never a silent surprise."""
    for site in sites:
        info = STATION_SITES.get(site)
        if not info:
            continue
        if info["status"] != "active":
            print(f"[info] Site {site} ({info['surface']}/{info['borehole']}) is not yet "
                  f"active: {info['note']}")
        elif info["note"]:
            print(f"[info] Site {site} ({info['surface']}/{info['borehole']}): {info['note']}")


# ---------------------------------------------------------------------------
# Folder layout
# ---------------------------------------------------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATASETS_DIR = os.path.join(BASE_DIR, "datasets")
PLOT_DIR = os.path.join(DATASETS_DIR, "plot")
MSEED_DIR = os.path.join(DATASETS_DIR, "mseed")
SONIFY_DIR = os.path.join(DATASETS_DIR, "sonifications")
MAP_DIR = os.path.join(DATASETS_DIR, "maps")
ASSETS_DIR = os.path.join(BASE_DIR, "assets")
GERMANY_OUTLINE_PATH = os.path.join(ASSETS_DIR, "germany_outline.json")

for _dir in (PLOT_DIR, MSEED_DIR, SONIFY_DIR, MAP_DIR):
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
# Germany civil time (CET/CEST), for showing local time alongside UTC on plots.
# Implemented manually (EU DST rule: last Sunday of March 01:00 UTC to last
# Sunday of October 01:00 UTC) rather than via the stdlib `zoneinfo` module, so
# no `tzdata` package is required on Windows (zoneinfo has no bundled IANA
# database there) -- keeps the "no extra dependencies" guarantee intact.
# ---------------------------------------------------------------------------
def _last_sunday(year, month):
    """The last Sunday of a given (year, month), as a naive datetime at 00:00."""
    if month == 12:
        next_month_first = datetime(year + 1, 1, 1)
    else:
        next_month_first = datetime(year, month + 1, 1)
    day = next_month_first - timedelta(days=1)
    while day.weekday() != 6:  # Monday=0 ... Sunday=6
        day -= timedelta(days=1)
    return day


def _cet_offset_hours(utc_dt):
    """German civil time UTC offset (2 for CEST/summer, 1 for CET/winter) for a
    naive-UTC datetime, per the EU DST rule."""
    year = utc_dt.year
    dst_start = _last_sunday(year, 3).replace(hour=1)
    dst_end = _last_sunday(year, 10).replace(hour=1)
    return 2 if dst_start <= utc_dt < dst_end else 1


def format_local_time(utc_dt):
    """Return (formatted local time string, tz abbreviation) for a naive-UTC
    datetime, in German civil time (CET/CEST)."""
    offset_h = _cet_offset_hours(utc_dt)
    local_dt = utc_dt + timedelta(hours=offset_h)
    tz_label = "CEST" if offset_h == 2 else "CET"
    return local_dt.strftime("%Y-%m-%d %H:%M"), tz_label


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
        help="What to do: any combination of fetch, plot, sonify, play, map, all "
             "(default: fetch plot sonify)",
    )

    station_group = parser.add_argument_group("station selection")
    station_group.add_argument(
        "--sites", default=None, metavar="SITES",
        help="Which DZA site number(s) to use, comma-separated (e.g. '1' or '1,3'), or "
             "'active' for all currently-streaming sites, or 'all' for every known site "
             "including ones not yet installed. Each site has a surface sensor (DZAx1) and "
             f"a borehole sensor at ~{BOREHOLE_DEPTH_M:.0f} m depth (DZAx3), both fetched "
             f"together. Known sites: {sorted(STATION_SITES)} (active: {ACTIVE_SITE_NUMBERS}). "
             "Default: site 1 (DZA11/DZA13) for fetch/plot/sonify/play; every known site for "
             "'map'. See README for per-site notes (stability, install status, sensor model).",
    )

    fetch_group = parser.add_argument_group("fetch options")
    fetch_group.add_argument(
        "--duration", type=float, default=DEFAULT_DURATION,
        help=f"Length of the data window to fetch, in minutes (default: {DEFAULT_DURATION:.0f})",
    )
    fetch_group.add_argument(
        "--hours-back", type=float, default=None, metavar="HOURS",
        help="Convenience alternative to --duration, in hours instead of minutes (pytremor-style "
             "lookback window): fetch the last HOURS of data, ending --latency minutes ago. "
             "e.g. --hours-back 24 for 'the last 24h from now'. Overrides --duration.",
    )
    fetch_group.add_argument(
        "--days-back", type=int, default=None, metavar="DAYS",
        help="Batch mode: fetch DAYS separate ~24h windows, stepping back one full day at a "
             "time from now (day 1 ends --latency minutes ago; day 2 is the 24h before that; "
             "etc). Each day is saved, plotted and sonified independently -- useful for "
             "building a multi-day dataset for review, e.g. --days-back 10 for the last 10 "
             "days as 10 separate daily files. Overrides --duration/--hours-back/--listen-minutes; "
             "'play' is skipped in this mode (too many files to play automatically).",
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
    fetch_group.add_argument(
        "--freqmin", type=float, default=0.5, metavar="HZ",
        help="Lower bandpass corner in Hz (default: 0.5). Lower this for long-period work, "
             "e.g. ocean-generated hum/microseism studies or the sub-10 mHz band researched "
             "at BFO (see README's 'Ultra-low-frequency band' section) -- e.g. --freqmin 0.0001. "
             "Below ~0.0083 Hz (the 120s sensors' corner), instrument-response removal becomes "
             "unreliable; the script will print a caution but not block the request.",
    )
    fetch_group.add_argument(
        "--freqmax", type=float, default=10.0, metavar="HZ",
        help="Upper bandpass corner in Hz (default: 10.0).",
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
             "Pass 'all' to sonify every channel in the file, each saved as its own .wav. "
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


VALID_ACTIONS = {"fetch", "plot", "sonify", "play", "map", "all"}


# ---------------------------------------------------------------------------
# Fetch
# ---------------------------------------------------------------------------
def fetch_once(client, latency_min, duration_min, station_pattern):
    """Request inventory + waveforms for a given latency/duration (both in minutes)."""
    endtime = UTCDateTime() - latency_min * 60
    starttime = endtime - duration_min * 60
    inv = client.get_stations(
        network=NETWORK, station=station_pattern, starttime=starttime, endtime=endtime,
        level="response",
    )
    st = client.get_waveforms(NETWORK, station_pattern, "*", CHANNEL, starttime, endtime)
    if len(st) == 0:
        raise ValueError("no traces returned (data not available yet)")
    return st, inv, starttime, endtime


def describe_stream(st, prefix="info", station_pattern=None):
    """Print a short summary of every trace in a Stream (station/channel/duration)."""
    pattern_note = f"station pattern '{station_pattern}'" if station_pattern else "the requested station(s)"
    print(f"[{prefix}] {len(st)} trace(s) in this file ({pattern_note} can match "
          f"multiple sub-stations/channels):")
    for tr in st:
        duration_s = tr.stats.npts / tr.stats.sampling_rate
        print(f"  - {tr.id}  |  {tr.stats.sampling_rate:.1f} Hz  |  {duration_s:.1f}s  |  {tr.stats.npts} samples")


def fetch_closest_to_real_time(client, latency_min, duration_min, max_latency_min, latency_step_min, station_pattern):
    """Try to get data as close to 'now' as possible, backing off only if needed."""
    latency = latency_min
    while True:
        try:
            st, inv, starttime, endtime = fetch_once(client, latency, duration_min, station_pattern)
            print(
                f"[fetch] Got data from {starttime} to {endtime} UTC "
                f"({latency:.1f} min latency margin from now, {duration_min:.1f} min duration). "
                f"Adjust the margin with --latency."
            )
            return st, inv, starttime, endtime
        except Exception as exc:
            if latency >= max_latency_min:
                raise RuntimeError(
                    f"Could not fetch data even at max latency of {max_latency_min:.1f} min: {exc}"
                ) from exc
            next_latency = min(latency + latency_step_min, max_latency_min)
            print(f"[fetch] No data at {latency:.1f} min latency ({exc}); retrying at {next_latency:.1f} min...")
            latency = next_latency


def do_fetch(args, station_pattern, file_prefix):
    """Fetch, process and save one waveform window. Returns the saved .mseed path."""
    client = Client(FDSN_URL)
    st, inv, starttime, endtime = fetch_closest_to_real_time(
        client, args.latency, args.duration, args.max_latency, args.latency_step, station_pattern
    )
    return _process_and_save(
        st, inv, starttime, prefix="fetch", file_prefix=file_prefix,
        freqmin=args.freqmin, freqmax=args.freqmax, station_pattern=station_pattern,
    )


def _process_and_save(st, inv, starttime, prefix="fetch", file_prefix="DZA1",
                       freqmin=0.5, freqmax=10.0, station_pattern=None):
    """Apply the standard processing chain and save as .mseed. Returns the saved path."""
    st.remove_response(inventory=inv)
    st.detrend("demean")
    st.detrend("linear")
    st.filter("bandpass", freqmin=freqmin, freqmax=freqmax, corners=4)
    st.taper(0.125)

    timestamp = starttime.strftime("%Y-%m-%d-%H-%M")
    mseed_path = os.path.join(MSEED_DIR, f"{file_prefix}_{timestamp}.mseed")
    st.write(mseed_path, format="MSEED")
    print(f"[{prefix}] Saved waveform data to {mseed_path}")
    describe_stream(st, prefix=prefix, station_pattern=station_pattern)

    # small .json sidecar recording the fetch/filter parameters, so `plot` can later
    # show the exact bandpass used even when re-run separately from `fetch`
    channel_orientations = {}
    for net in inv:
        for sta in net:
            for cha in sta:
                key = f"{net.code}.{sta.code}.{cha.location_code}.{cha.code}"
                channel_orientations[key] = dict(azimuth=cha.azimuth, dip=cha.dip)

    metadata = dict(
        network=NETWORK, station_pattern=station_pattern, channel=CHANNEL,
        freqmin=freqmin, freqmax=freqmax,
        starttime=str(st[0].stats.starttime), endtime=str(st[0].stats.endtime),
        channel_orientations=channel_orientations,
    )
    metadata_path = os.path.splitext(mseed_path)[0] + ".json"
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)

    return mseed_path


def do_fetch_batch(days_back, latency_min, max_latency_min, latency_step_min,
                    station_pattern, file_prefix, freqmin, freqmax):
    """Batch-fetch `days_back` separate ~24h windows, going back day by day from now
    (day 1 = the most recent 24h ending `latency_min` minutes ago, day 2 = the 24h
    before that, etc). Each day is fetched, processed and saved independently -- a
    day with no/partial data is skipped with a warning rather than aborting the
    whole batch. Returns a chronologically sorted list of saved .mseed paths.
    """
    client = Client(FDSN_URL)
    day_minutes = 24 * 60
    mseed_paths = []
    print(f"[fetch] Batch mode: {days_back} day(s) x 24h, starting {latency_min:.1f} min "
          f"before now and stepping back one day at a time.")
    for day_index in range(days_back):
        day_latency = latency_min + day_index * day_minutes
        try:
            st, inv, starttime, endtime = fetch_closest_to_real_time(
                client, day_latency, day_minutes, day_latency + max_latency_min, latency_step_min,
                station_pattern,
            )
        except Exception as exc:
            print(f"[fetch] Day {day_index + 1}/{days_back}: no data available ({exc}) -- skipping.")
            continue
        mseed_path = _process_and_save(
            st, inv, starttime, prefix="fetch", file_prefix=file_prefix,
            freqmin=freqmin, freqmax=freqmax, station_pattern=station_pattern,
        )
        print(f"[fetch] Day {day_index + 1}/{days_back} done -> {mseed_path}")
        mseed_paths.append(mseed_path)

    mseed_paths.sort()  # filenames are timestamp-based, so this is chronological order
    print(f"[fetch] Batch complete: {len(mseed_paths)}/{days_back} day(s) fetched successfully.")
    return mseed_paths


# ---------------------------------------------------------------------------
# Map (station locations)
# ---------------------------------------------------------------------------
def haversine_km(lat1, lon1, lat2, lon2):
    """Great-circle distance between two lat/lon points, in km."""
    r = 6371.0
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dphi = math.radians(lat2 - lat1)
    dlambda = math.radians(lon2 - lon1)
    a = math.sin(dphi / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlambda / 2) ** 2
    return 2 * r * math.asin(math.sqrt(a))


_germany_outline_cache = None


def load_germany_outline():
    """Load the bundled low-resolution Germany border outline (a simple list of
    [lon, lat] points, see assets/germany_outline.json) for use as visual context
    behind the station markers. Cached after first read. Returns an empty list
    (and the map is drawn without a country outline) if the asset is missing."""
    global _germany_outline_cache
    if _germany_outline_cache is not None:
        return _germany_outline_cache
    try:
        with open(GERMANY_OUTLINE_PATH) as f:
            _germany_outline_cache = json.load(f)["coordinates"]
    except (OSError, ValueError, KeyError):
        _germany_outline_cache = []
    return _germany_outline_cache


def site_points_from_registry(site_nums):
    """Build the site_points structure (site_num -> lat/lon/elevation/stations) from
    the hardcoded STATION_SITES registry -- no network call needed, so this can be
    used from do_plot without slowing every plot down with a metadata query. Sites
    without confirmed coordinates yet (e.g. planned installs) are omitted."""
    points = {}
    for site_num in site_nums:
        info = STATION_SITES.get(site_num)
        if not info or info.get("lat") is None:
            continue
        stations = []
        for code in (info["surface"], info["borehole"]):
            s_info = STATION_INFO_BY_CODE.get(code)
            if not s_info:
                continue
            sensor_note = f", {s_info['sensor_model']}" if s_info["depth_label"] == "surface" else ""
            stations.append((code, s_info["depth_label"], s_info["depth_m"], sensor_note))
        points[site_num] = dict(
            lat=info["lat"], lon=info["lon"], elevation=info["elevation_m"], stations=stations,
        )
    return points


def draw_station_map(ax, site_points, germany_outline=None, title=None, marker_size=200):
    """Draw a station-location map onto `ax`: the Germany border outline (if
    available) for country-level context, plus one marker per site (surface and
    borehole sensors are colocated, so a single marker covers both)."""
    ax.set_facecolor(BG_COLOR)
    site_nums = sorted(site_points)

    if germany_outline:
        poly_lons = [p[0] for p in germany_outline]
        poly_lats = [p[1] for p in germany_outline]
        ax.fill(poly_lons, poly_lats, facecolor="#1c2530", edgecolor="#4a5a6a",
                 linewidth=1.0, zorder=1)
        mean_lat = sum(poly_lats) / len(poly_lats)
    elif site_points:
        mean_lat = sum(p["lat"] for p in site_points.values()) / len(site_points)
    else:
        mean_lat = 51.0
    ax.set_aspect(1.0 / max(math.cos(math.radians(mean_lat)), 1e-6))

    cmap = plt.get_cmap("tab10")
    if germany_outline:
        # Whole-country context: nearby sites (a few km apart) can end up almost on
        # top of each other at this zoom level, so full per-marker text boxes would
        # overlap. Use short numbered markers instead, plus one combined legend box
        # with the full station details, anchored in a fixed corner.
        legend_lines = []
        for i, site_num in enumerate(site_nums):
            pt = site_points[site_num]
            color = cmap(i % 10)
            ax.scatter(pt["lon"], pt["lat"], s=marker_size * 1.6, color=color,
                       edgecolor="white", linewidth=2.0, zorder=3)
            ax.annotate(str(site_num), xy=(pt["lon"], pt["lat"]), xytext=(0, 0),
                        textcoords="offset points", color="white", fontsize=8,
                        fontweight="bold", ha="center", va="center", zorder=4)
            legend_lines.append(f"Site {site_num}:")
            for code, label, depth, note in pt["stations"]:
                depth_part = f", ~{depth:.0f} m" if label == "borehole" else ""
                legend_lines.append(f"  {code} ({label}{depth_part}{note})")
        if legend_lines:
            ax.text(
                0.02, 0.02, "\n".join(legend_lines), transform=ax.transAxes,
                color="white", fontsize=7, family="monospace", va="bottom", ha="left",
                zorder=5,
                bbox=dict(boxstyle="round", facecolor="#1c1f26", edgecolor=GRID_COLOR, alpha=0.92),
            )
    else:
        for i, site_num in enumerate(site_nums):
            pt = site_points[site_num]
            color = cmap(i % 10)
            ax.scatter(pt["lon"], pt["lat"], s=marker_size, color=color, edgecolor="white",
                       linewidth=1.2, zorder=3)
            station_lines = []
            for code, label, depth, note in pt["stations"]:
                depth_part = f", ~{depth:.0f} m" if label == "borehole" else ""
                station_lines.append(f"{code} ({label}{depth_part}{note})")
            label_text = f"Site {site_num}\n" + "\n".join(station_lines)
            ax.annotate(
                label_text, xy=(pt["lon"], pt["lat"]), xytext=(10, 7), textcoords="offset points",
                color="white", fontsize=7.5, family="monospace",
                bbox=dict(boxstyle="round", facecolor="#1c1f26", edgecolor=GRID_COLOR, alpha=0.9),
                zorder=4,
            )

    if germany_outline:
        # Whole-country context: zoom to the outline itself, not just the site
        # cluster, so it's obvious *where in Germany* the sites are.
        ax.set_xlim(min(poly_lons) - 0.3, max(poly_lons) + 0.3)
        ax.set_ylim(min(poly_lats) - 0.3, max(poly_lats) + 0.3)
    elif site_points:
        lons = [p["lon"] for p in site_points.values()]
        lats = [p["lat"] for p in site_points.values()]
        lon_pad = max((max(lons) - min(lons)) * 0.4, 0.05)
        lat_pad = max((max(lats) - min(lats)) * 0.4, 0.05)
        ax.set_xlim(min(lons) - lon_pad, max(lons) + lon_pad)
        ax.set_ylim(min(lats) - lat_pad, max(lats) + lat_pad)
    else:
        ax.text(0.5, 0.5, "No location data available", color=FG_COLOR, fontsize=9,
                ha="center", va="center", transform=ax.transAxes)

    ax.set_xlabel("Longitude (\u00b0E)", color=FG_COLOR, fontsize=8)
    ax.set_ylabel("Latitude (\u00b0N)", color=FG_COLOR, fontsize=8)
    if title:
        ax.set_title(title, color="white", fontsize=10)
    ax.grid(True, color=GRID_COLOR, linewidth=0.4, alpha=0.5)
    style_axes(ax)


def draw_depth_scale(ax, site_points):
    """Draw a simple vertical depth cross-section: one column per site, marking the
    surface sensor (0 m) and borehole sensor (~240 m) so the vertical separation
    between the two sensor types is obvious at a glance. Uses depth metadata only
    (no coordinates needed), so this works even for sites without map data."""
    ax.set_facecolor(BG_COLOR)
    site_nums = sorted(site_points)
    cmap = plt.get_cmap("tab10")

    all_depths = [depth for pt in site_points.values() for _, _, depth, _ in pt["stations"]]
    max_depth = max(all_depths) if all_depths else BOREHOLE_DEPTH_M

    ax.axhline(0, color=FG_COLOR, linewidth=1.0, alpha=0.8, zorder=1)
    ax.text(-0.5, 0, " ground level", color=FG_COLOR, fontsize=7, ha="left", va="bottom")

    for i, site_num in enumerate(site_nums):
        pt = site_points[site_num]
        color = cmap(i % 10)
        depths = sorted(depth for _, _, depth, _ in pt["stations"])
        if len(depths) > 1:
            ax.plot([i, i], [depths[0], depths[-1]], color=color, linewidth=1.2,
                     linestyle=":", zorder=2)
        for code, label, depth, _note in pt["stations"]:
            marker = "^" if label == "surface" else "v"
            ax.scatter(i, depth, s=90, color=color, edgecolor="white", linewidth=1.0,
                       marker=marker, zorder=3)
            ax.annotate(f"{code}\n~{depth:.0f} m", xy=(i, depth), xytext=(7, 0),
                        textcoords="offset points", color="white", fontsize=6.5,
                        family="monospace", va="center")

    if site_nums:
        ax.set_xlim(-0.6, len(site_nums) - 0.4 + 0.6)
        ax.set_xticks(range(len(site_nums)))
        tick_labels = []
        for site_num in site_nums:
            label = f"Site {site_num}"
            # Annotate each site with the distance to its nearest neighbor among the
            # other sites shown (using the hardcoded STATION_SITES coordinates, so
            # this works here even though this panel otherwise only needs depths).
            nearest = None
            for other in site_nums:
                if other == site_num:
                    continue
                a, b = STATION_SITES.get(site_num), STATION_SITES.get(other)
                if not a or not b or a.get("lat") is None or b.get("lat") is None:
                    continue
                dist_km = haversine_km(a["lat"], a["lon"], b["lat"], b["lon"])
                if nearest is None or dist_km < nearest[0]:
                    nearest = (dist_km, other)
            if nearest:
                label += f"\n({nearest[0]:.1f} km to Site {nearest[1]})"
            tick_labels.append(label)
        ax.set_xticklabels(tick_labels, color=FG_COLOR, fontsize=7.5)
    else:
        ax.set_xlim(-1, 1)
        ax.set_xticks([])
        ax.text(0.5, 0.5, "No depth data available", color=FG_COLOR, fontsize=9,
                ha="center", va="center", transform=ax.transAxes)

    ax.set_ylim(max_depth * 1.15, -max_depth * 0.08)  # 0 near top, deeper goes down
    ax.set_ylabel("Depth (m)", color=FG_COLOR, fontsize=8)
    ax.set_title("Sensor depth", color="white", fontsize=10)
    ax.grid(True, axis="y", color=GRID_COLOR, linewidth=0.4, alpha=0.5)
    style_axes(ax)


def do_map(sites):
    """Print coordinates/depth/status for the given sites' stations and save a
    self-contained matplotlib station-map plot (no extra dependencies needed --
    this does not use ObsPy's Inventory.plot/cartopy). A wide, fixed time range is
    used for the metadata query so planned-but-not-yet-installed stations don't
    cause an error -- they simply won't appear in the result.
    """
    station_pattern = build_station_pattern(sites)
    client = Client(FDSN_URL)
    print(f"[map] Querying station metadata for: {station_pattern}")
    inv = client.get_stations(
        network=NETWORK, station=station_pattern,
        starttime=UTCDateTime(2020, 1, 1), endtime=UTCDateTime(2035, 1, 1),
        level="channel",
    )

    # site_num -> dict(lat, lon, elevation, stations=[(code, depth_label, depth_m, sensor_note)])
    site_points = {}
    found_any = False
    for net in inv:
        for sta in net:
            info = STATION_INFO_BY_CODE.get(sta.code)
            depth_label = info["depth_label"] if info else ("borehole" if sta.code.endswith("3") else "surface")
            depth_m = info["depth_m"] if info else (BOREHOLE_DEPTH_M if depth_label == "borehole" else 0.0)
            sensor_note = f", {info['sensor_model']}" if info and depth_label == "surface" else ""
            found_any = True
            print(
                f"  - {sta.code}  ({depth_label}, ~{depth_m:.0f} m{sensor_note})  "
                f"lat={sta.latitude:.5f}  lon={sta.longitude:.5f}  elevation={sta.elevation:.1f} m"
            )
            if info:
                print(f"      status: {info['status']}  |  {info['note']}")
                site_points.setdefault(info["site"], dict(
                    lat=sta.latitude, lon=sta.longitude, elevation=sta.elevation, stations=[],
                ))["stations"].append((sta.code, depth_label, depth_m, sensor_note))

    if not found_any:
        print(f"[map] No station metadata found for {station_pattern}.")
        return inv, None

    site_nums = sorted(site_points)
    if len(site_nums) > 1:
        print("[map] Inter-site distances:")
        for a_idx in range(len(site_nums)):
            for b_idx in range(a_idx + 1, len(site_nums)):
                a, b = site_nums[a_idx], site_nums[b_idx]
                dist_km = haversine_km(
                    site_points[a]["lat"], site_points[a]["lon"],
                    site_points[b]["lat"], site_points[b]["lon"],
                )
                print(f"    site {a} <-> site {b}: {dist_km:.2f} km")

    # -- self-contained matplotlib map (no cartopy/basemap dependency), with the
    #    bundled Germany outline for country-level context, plus a depth panel --
    germany_outline = load_germany_outline()
    fig, (map_ax, depth_ax) = plt.subplots(
        1, 2, figsize=(12, 6.5), facecolor=BG_COLOR, gridspec_kw={"width_ratios": [2.1, 1]},
    )
    draw_station_map(map_ax, site_points, germany_outline=germany_outline,
                      title=f"DZA network station map ({station_pattern})", marker_size=220)
    draw_depth_scale(depth_ax, site_points)

    map_path = os.path.join(MAP_DIR, "station_map.png")
    fig.tight_layout()
    fig.savefig(map_path, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f"[map] Saved station map to {map_path}")
    return inv, map_path


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
def component_color(channel_code, depth_label=None):
    """Base color for a channel's component (Z/N/E), shaded by sensor depth: surface
    sensors get a brighter/lighter version, borehole sensors a darker/muted one, so
    depth is visually encoded in the waveform color too, not just the text label."""
    comp = channel_code[-1].upper() if channel_code else "Z"
    base = COMPONENT_COLORS.get(comp, "#76ff03")
    return _shade_for_depth(base, depth_label)


def _shade_for_depth(hex_color, depth_label):
    """Lighten (surface) or darken (borehole) a hex color's lightness (HLS), keeping
    hue/saturation intact. `depth_label` of None leaves the color unchanged."""
    r, g, b = mcolors.to_rgb(hex_color)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    if depth_label == "borehole":
        l = max(0.0, l * 0.55)
    elif depth_label == "surface":
        l = min(0.95, l * 1.15 + 0.05)
    r2, g2, b2 = colorsys.hls_to_rgb(h, l, s)
    return (r2, g2, b2)


def rms_envelope(data, sampling_rate, window_s=2.0):
    """Rolling RMS envelope, used as a translucent fill behind each waveform."""
    # np.convolve(..., mode="same") returns an array of length max(len(data),
    # len(kernel)) -- if the requested window is longer than a very short trace
    # (e.g. a few seconds of data from a telemetry-latency-truncated channel),
    # the result would come back *longer* than `data`, breaking the caller's
    # assumption that times/data/env all share the same length. Clamp the
    # window so the kernel is never longer than the data being convolved.
    if len(data) == 0:
        return np.array([], dtype=np.float64)
    window = max(1, min(int(window_s * sampling_rate), len(data)))
    squared = data.astype(np.float64) ** 2
    kernel = np.ones(window) / window
    return np.sqrt(np.convolve(squared, kernel, mode="same"))


def pick_primary_for_spectrogram(st):
    """Prefer a vertical (Z) component for the top spectrogram panel -- and among
    multiple Z channels (e.g. several sites/sensors in one plot), prefer whichever
    has the most complete data, so the spectrogram is as fully populated as
    possible rather than picking a channel that happens to have a telemetry gap."""
    z_traces = [tr for tr in st if tr.stats.channel.upper().endswith("Z")]
    candidates = z_traces or list(st)
    return max(candidates, key=lambda tr: tr.stats.npts / tr.stats.sampling_rate)



def style_axes(ax):
    ax.set_facecolor(BG_COLOR)
    ax.tick_params(colors=FG_COLOR, labelsize=8)
    for spine in ax.spines.values():
        spine.set_color(GRID_COLOR)


def load_fetch_metadata(mseed_path):
    """Load the .json sidecar written alongside a .mseed by _process_and_save, if any.
    Returns a dict with at least 'freqmin'/'freqmax' (falling back to the tool's
    historical defaults for older files saved before this sidecar existed)."""
    metadata_path = os.path.splitext(mseed_path)[0] + ".json"
    defaults = {"freqmin": 0.5, "freqmax": 10.0, "channel_orientations": {}}
    if not os.path.exists(metadata_path):
        return defaults
    try:
        with open(metadata_path) as f:
            data = json.load(f)
    except (OSError, ValueError):
        return defaults
    defaults.update(data)
    return defaults


def _pick_time_axis_unit(duration_s):
    """Pick a human-friendly unit for the waveform x-axis based on the recording's
    actual duration, so a 10-hour fetch reads as '0-10' on a 'Time (hours)' axis
    instead of '0-36000' raw seconds -- matching how long the sonification of the
    same recording actually lasts (before any --speed-up is applied)."""
    if duration_s >= 2 * 3600:
        return 3600.0, "Time (hours)", "hours"
    if duration_s >= 180:
        return 60.0, "Time (minutes)", "minutes"
    return 1.0, "Time (seconds)", "seconds"


def _depth_points_for_stations(station_codes):
    """Build a site_num -> {'stations': [...]} structure purely for the depth-scale
    panel, using whichever stations are actually present -- no coordinates needed,
    so this works even for the (unlikely) case of a site with data but no confirmed
    lat/lon yet."""
    points = {}
    for code in station_codes:
        info = STATION_INFO_BY_CODE.get(code)
        if not info:
            continue
        sensor_note = f", {info['sensor_model']}" if info["depth_label"] == "surface" else ""
        points.setdefault(info["site"], dict(stations=[]))["stations"].append(
            (code, info["depth_label"], info["depth_m"], sensor_note)
        )
    return points


def do_plot(mseed_path):
    st = read(mseed_path)
    base = os.path.splitext(os.path.basename(mseed_path))[0]
    plot_path = os.path.join(PLOT_DIR, f"{base}.png")
    metadata = load_fetch_metadata(mseed_path)
    freqmin, freqmax = metadata["freqmin"], metadata["freqmax"]

    n_traces = len(st)
    primary = pick_primary_for_spectrogram(st)
    sr = primary.stats.sampling_rate
    stations_shown = sorted({tr.stats.station for tr in st})
    site_nums_shown = sorted({
        STATION_INFO_BY_CODE[code]["site"] for code in stations_shown if code in STATION_INFO_BY_CODE
    })
    map_site_points = site_points_from_registry(site_nums_shown)
    depth_site_points = _depth_points_for_stations(stations_shown)
    germany_outline = load_germany_outline()
    channel_orientations = metadata.get("channel_orientations", {})

    # Consistent site -> color mapping shared across the map, depth panel, and each
    # waveform panel's label accent, so a given site is visually the same color
    # everywhere in the figure (same ordering/index as draw_station_map/draw_depth_scale).
    tab10 = plt.get_cmap("tab10")
    site_color_map = {site_num: tab10(i % 10) for i, site_num in enumerate(site_nums_shown)}

    # Every trace is plotted on one shared, absolute time axis (rather than each
    # panel auto-zooming to its own data range) so that a channel with less data
    # than the others -- a real near-real-time telemetry gap between sensors, not
    # a bug -- shows up as visible blank space in the *right place*, instead of
    # being silently stretched to fill the panel and hiding the gap. This is also
    # what was cutting the spectrogram panel short: its (single) primary trace
    # sometimes has less data than other channels, but its axis was previously
    # stretched to the full shared duration.
    t0 = min(tr.stats.starttime for tr in st)
    duration_s = max(
        (tr.stats.starttime - t0) + tr.stats.npts / tr.stats.sampling_rate for tr in st
    )
    time_divisor, time_label, time_unit_word = _pick_time_axis_unit(duration_s)

    fig = plt.figure(figsize=(17, max(6.5, 3 + 1.3 * n_traces)), facecolor=BG_COLOR)
    outer = fig.add_gridspec(1, 2, width_ratios=[2.5, 1], wspace=0.24)
    left_gs = outer[0].subgridspec(1 + n_traces, 1, height_ratios=[2.2] + [1] * n_traces, hspace=0.15)
    right_gs = outer[1].subgridspec(2, 1, height_ratios=[1.15, 1], hspace=0.45)

    # -- right column: station location map (Germany context) + sensor depth scale --
    map_ax = fig.add_subplot(right_gs[0])
    depth_ax = fig.add_subplot(right_gs[1])
    draw_station_map(map_ax, map_site_points, germany_outline=germany_outline,
                      title="Station location", marker_size=140)
    draw_depth_scale(depth_ax, depth_site_points)

    # -- top-left panel: spectrogram of the primary (vertical) channel --
    spec_ax = fig.add_subplot(left_gs[0])
    nperseg = int(max(32, min(sr * 4, len(primary.data))))
    noverlap = int(nperseg * 0.9)
    f, t, Sxx = spectrogram(primary.data.astype(np.float64), fs=sr, nperseg=nperseg, noverlap=noverlap)
    Sxx_db = 10 * np.log10(Sxx + 1e-30)
    vmin, vmax = np.percentile(Sxx_db, [5, 99.5])
    primary_offset_s = primary.stats.starttime - t0
    mesh = spec_ax.pcolormesh((t + primary_offset_s) / time_divisor, f, Sxx_db, cmap="inferno",
                               shading="auto", vmin=vmin, vmax=vmax)
    y_top = min(sr / 2, max(12.0, freqmax * 1.5))
    spec_ax.set_ylim(0, y_top)
    spec_ax.set_xlim(0, duration_s / time_divisor)
    spec_ax.set_ylabel("Frequency (Hz)", color=FG_COLOR)
    spec_title = f"Spectrogram - {primary.id}"
    primary_duration_s = primary.stats.npts / primary.stats.sampling_rate
    primary_coverage_pct = 100 * primary_duration_s / duration_s
    if primary_coverage_pct < 99:
        spec_title += f"  ({primary_coverage_pct:.0f}% of window -- see gap note below)"
    spec_ax.set_title(spec_title, color="white", fontsize=10, loc="left")
    spec_ax.set_xticklabels([])
    style_axes(spec_ax)
    cbar = fig.colorbar(mesh, ax=spec_ax, pad=0.01, fraction=0.02)
    cbar.set_label("Power (dB)", color=FG_COLOR)
    cbar.ax.yaxis.set_tick_params(color=FG_COLOR)
    plt.setp(cbar.ax.get_yticklabels(), color=FG_COLOR)

    # -- mark the applied bandpass on the spectrogram, so the visible band is never
    #    ambiguous at a glance --
    for edge, label in ((freqmin, "freqmin"), (freqmax, "freqmax")):
        if 0 < edge < y_top:
            spec_ax.axhline(edge, color="white", linewidth=0.8, linestyle="--", alpha=0.6)
            spec_ax.text(
                0.995, edge, f"{label}={edge:g} Hz", transform=spec_ax.get_yaxis_transform(),
                color="white", fontsize=7, va="bottom", ha="right",
                bbox=dict(boxstyle="round,pad=0.15", facecolor=BG_COLOR, edgecolor="none", alpha=0.7),
            )

    # -- mark the dominant frequency band (highest average power across the whole
    #    window) so "what's the loudest band" is answered at a glance --
    mean_power_per_freq = Sxx_db.mean(axis=1)
    peak_freq_hz = f[np.argmax(mean_power_per_freq)]
    if 0 < peak_freq_hz < y_top:
        spec_ax.axhline(peak_freq_hz, color="#39ff14", linewidth=1.0, linestyle="-", alpha=0.8)
        spec_ax.text(
            0.005, peak_freq_hz, f"peak~{peak_freq_hz:.2g} Hz", transform=spec_ax.get_yaxis_transform(),
            color="#39ff14", fontsize=7, va="bottom", ha="left",
            bbox=dict(boxstyle="round,pad=0.15", facecolor=BG_COLOR, edgecolor="none", alpha=0.7),
        )

    # -- one panel per trace: waveform + RMS envelope, labeled with site metadata --
    max_points = 200_000  # cap rendered points per trace for speed/memory on long recordings
    for i, tr in enumerate(st):
        ax = fig.add_subplot(left_gs[i + 1])
        offset_s = tr.stats.starttime - t0
        trace_duration_s = tr.stats.npts / tr.stats.sampling_rate
        full_data = tr.data.astype(np.float64)
        peak_amp = np.max(np.abs(full_data))
        rms_amp = np.sqrt(np.mean(full_data ** 2))
        times = tr.times() + offset_s
        data = tr.data
        step = max(1, len(data) // max_points)
        if step > 1:
            times = times[::step]
            data = data[::step]

        info = STATION_INFO_BY_CODE.get(tr.stats.station)
        depth_label = info["depth_label"] if info else None
        color = component_color(tr.stats.channel, depth_label)
        env = rms_envelope(data, tr.stats.sampling_rate / step)
        times_scaled = times / time_divisor
        ax.fill_between(times_scaled, -env, env, color=color, alpha=0.25, linewidth=0)
        ax.plot(times_scaled, data, color=color, linewidth=0.5)
        # shared absolute axis across every panel (see note above `t0` for why) --
        # a trace with less data than the window simply stops early, showing the
        # gap instead of hiding it.
        ax.set_xlim(0, duration_s / time_divisor)

        if info and info["depth_label"] == "surface":
            site_line = f"surface, {info['sensor_model']} ({info['eigenperiod_s']:.0f}s eigenperiod)"
        elif info:
            site_line = f"borehole, ~{info['depth_m']:.0f} m depth"
        else:
            site_line = None

        label_lines = [f"{tr.id}  |  {tr.stats.sampling_rate:.0f} Hz"]
        if site_line:
            label_lines.append(site_line)

        orientation = channel_orientations.get(tr.id)
        if orientation and orientation.get("azimuth") is not None and orientation.get("dip") is not None:
            label_lines.append(f"az {orientation['azimuth']:.0f}\u00b0, dip {orientation['dip']:.0f}\u00b0")

        label_lines.append(f"peak {peak_amp:.2e} m/s | rms {rms_amp:.2e} m/s")

        coverage_pct = 100 * trace_duration_s / duration_s
        if coverage_pct < 99:
            end_gap_s = duration_s - (offset_s + trace_duration_s)
            gap_bits = []
            if offset_s > max(1.0, duration_s * 0.01):
                gap_bits.append(f"starts {offset_s:.0f}s late")
            if end_gap_s > max(1.0, duration_s * 0.01):
                gap_bits.append(f"ends {end_gap_s:.0f}s early")
            gap_note = ", ".join(gap_bits) or "partial data"
            label_lines.append(f"{coverage_pct:.0f}% of window ({gap_note})")
        label_text = "\n".join(label_lines)

        accent_color = site_color_map.get(info["site"]) if info else GRID_COLOR
        ax.text(
            0.005, 0.92, label_text, transform=ax.transAxes, color="white", fontsize=8,
            family="monospace", va="top",
            bbox=dict(boxstyle="round", facecolor="#1c1f26", edgecolor=accent_color, linewidth=1.6, alpha=0.88),
        )
        style_axes(ax)
        if i < n_traces - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel(time_label, color=FG_COLOR)

    start_local_str, tz_label = format_local_time(primary.stats.starttime.datetime)
    end_local_str, _ = format_local_time(primary.stats.endtime.datetime)
    title = (
        f"{NETWORK}.{','.join(stations_shown)}  |  {primary.stats.starttime} - {primary.stats.endtime} UTC  "
        f"({start_local_str} - {end_local_str} {tz_label})  |  bandpass {freqmin:g}-{freqmax:g} Hz"
    )
    fig.suptitle(title, color="white", fontsize=11, y=0.995)
    fig.subplots_adjust(left=0.045, right=0.97, top=0.94, bottom=0.05, hspace=0.15)
    fig.savefig(plot_path, dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)

    print(f"[plot] Saved plot to {plot_path} ({n_traces} traces shown, {freqmin:g}-{freqmax:g} Hz bandpass, "
          f"{duration_s / time_divisor:.1f} {time_unit_word})")
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
              f"Ignoring: {', '.join(others)}. Use --channel to pick a different one, "
              f"or --channel all to sonify every one of them.")
    return tr


def _sonify_one_trace(tr, mseed_path, speed_up_factor, channel_tag=None):
    """Sonify a single trace to a .wav file. channel_tag, if given, is embedded in the
    filename so multiple channels from the same recording don't collide (used by
    --channel all)."""
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
    if channel_tag:
        wav_path = os.path.join(SONIFY_DIR, f"{base}_{channel_tag}_{int(speed_up_factor)}x.wav")
    else:
        wav_path = os.path.join(SONIFY_DIR, f"{base}_{int(speed_up_factor)}x.wav")
    wavfile.write(wav_path, wav_sample_rate, audio)
    print(f"[sonify] Saved audio to {wav_path}")
    print(f"[sonify]   trace: {tr.id}  |  input: {input_duration_s:.1f}s  ->  "
          f"output: {output_duration_s:.1f}s at {speed_up_factor:.0f}x speed "
          f"(wav sample rate {wav_sample_rate} Hz)")
    return wav_path


def do_sonify(mseed_path, speed_up_factor=DEFAULT_SPEED_UP, channel_filter=None):
    st = read(mseed_path)

    if channel_filter and channel_filter.strip().lower() == "all":
        print(f"[sonify] --channel all: sonifying all {len(st)} channel(s) into separate .wav files.")
        wav_paths = []
        # A recording with a mid-stream data gap (not just a start/end truncation)
        # is read back as *multiple* Trace objects sharing the same seed ID (ObsPy
        # does not auto-merge disjoint segments). Without disambiguation, every
        # segment after the first would silently overwrite the previous file under
        # the same channel_tag -- losing real data rather than erroring, so each
        # repeat of a tag gets a "_segN" suffix.
        seen_tag_counts = {}
        for tr in st:
            channel_tag = tr.id.replace(".", "-")
            seen_tag_counts[channel_tag] = seen_tag_counts.get(channel_tag, 0) + 1
            occurrence = seen_tag_counts[channel_tag]
            if occurrence > 1:
                channel_tag = f"{channel_tag}_seg{occurrence}"
            wav_paths.append(_sonify_one_trace(tr, mseed_path, speed_up_factor, channel_tag))
        return wav_paths

    tr = pick_trace(st, channel_filter)
    return _sonify_one_trace(tr, mseed_path, speed_up_factor)


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


def find_wav_for(mseed_path, speed_up_factor, channel_filter=None):
    base = os.path.splitext(os.path.basename(mseed_path))[0]

    if channel_filter and channel_filter.strip().lower() == "all":
        matches = sorted(glob.glob(os.path.join(SONIFY_DIR, f"{base}_*_{int(speed_up_factor)}x.wav")))
        return matches if matches else None

    candidate = os.path.join(SONIFY_DIR, f"{base}_{int(speed_up_factor)}x.wav")
    if os.path.exists(candidate):
        return candidate
    # fall back to any single-channel sonification of this same recording, any speed
    # (the [0-9] guard avoids matching the per-channel files from --channel all)
    matches = sorted(glob.glob(os.path.join(SONIFY_DIR, f"{base}_[0-9]*x.wav")))
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
    if args.hours_back is not None and args.hours_back <= 0:
        raise SystemExit(f"--hours-back must be > 0 (got {args.hours_back})")
    if args.latency < 0:
        raise SystemExit(f"--latency must be >= 0 (got {args.latency})")
    if args.max_latency < args.latency:
        raise SystemExit("--max-latency must be >= --latency")
    if args.speed_up is not None and args.speed_up <= 0:
        raise SystemExit(f"--speed-up must be > 0 (got {args.speed_up})")
    if args.listen_minutes is not None and args.listen_minutes <= 0:
        raise SystemExit(f"--listen-minutes must be > 0 (got {args.listen_minutes})")
    if args.days_back is not None and args.days_back <= 0:
        raise SystemExit(f"--days-back must be > 0 (got {args.days_back})")
    if args.freqmin <= 0:
        raise SystemExit(f"--freqmin must be > 0 (got {args.freqmin})")
    if args.freqmax <= args.freqmin:
        raise SystemExit(f"--freqmax ({args.freqmax}) must be > --freqmin ({args.freqmin})")

    if args.list:
        print_file_list()
        return

    actions = set(args.actions) if args.actions else {"fetch", "plot", "sonify"}
    if "all" in actions:
        actions = {"fetch", "plot", "sonify"}

    if "map" in actions:
        map_sites = parse_sites_arg(args.sites, fallback=sorted(STATION_SITES))
        do_map(map_sites)
        actions.discard("map")
        if not actions:
            return

    fetch_sites = parse_sites_arg(args.sites, fallback=[1])
    station_pattern = build_station_pattern(fetch_sites)
    file_prefix = build_file_prefix(fetch_sites)

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

    if fetch_needed:
        print_site_advisories(fetch_sites)
        if args.freqmin < 1.0 / DEFAULT_EIGENPERIOD_S:
            print(
                f"[caution] --freqmin {args.freqmin:g} Hz is below the ~"
                f"{1.0 / DEFAULT_EIGENPERIOD_S:.4f} Hz corner of the 120s surface sensors "
                f"(and the 0.05 Hz corner of DZA11's 20s sensor). Instrument-response removal "
                f"amplifies noise heavily below a sensor's corner frequency, so results this low "
                f"should be treated as exploratory/artistic rather than scientifically validated "
                f"without further review. See README's 'Ultra-low-frequency band' section."
            )

    if args.days_back is not None and fetch_needed:
        if args.listen_minutes:
            print("[info] --days-back overrides --listen-minutes; ignoring --listen-minutes.")
        effective_speed_up = args.speed_up if speed_up_explicit else DEFAULT_SPEED_UP
        mseed_paths = do_fetch_batch(
            args.days_back, args.latency, args.max_latency, args.latency_step,
            station_pattern, file_prefix, args.freqmin, args.freqmax,
        )
        if not mseed_paths:
            raise SystemExit("Batch fetch produced no files (no data available for any requested day).")
        for i, path in enumerate(mseed_paths):
            print(f"[batch] Processing file {i + 1}/{len(mseed_paths)}: {os.path.basename(path)}")
            if "plot" in actions:
                do_plot(path)
            if "sonify" in actions:
                do_sonify(path, effective_speed_up, args.channel)
        if "play" in actions:
            print("[batch] 'play' is skipped in --days-back batch mode (too many files to auto-play).")
        print(f"[batch] Done: {len(mseed_paths)} day(s) fetched and processed.")
        return

    if args.hours_back is not None and fetch_needed:
        args.duration = args.hours_back * 60
        print(
            f"[info] --hours-back {args.hours_back:.1f}h -> fetching the last "
            f"{args.duration:.0f} min of data (overrides --duration)."
        )

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
        mseed_path = do_fetch(args, station_pattern, file_prefix)
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
        describe_stream(read(mseed_path), prefix="info", station_pattern=station_pattern)

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
            wav_path = find_wav_for(mseed_path, args.speed_up, args.channel)
        if wav_path is None:
            wav_path = do_sonify(mseed_path, args.speed_up, args.channel)
        for p in (wav_path if isinstance(wav_path, list) else [wav_path]):
            do_play(p)


if __name__ == "__main__":
    main()
