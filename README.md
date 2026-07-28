# DZA01 Seismic Sonification Toolkit

A single-file Python CLI that pulls real seismic waveform data from the
**DZA seismic station network (`KB` network, KIT/GPI, Germany)**, plots
it as a spectrogram + waveform figure, and **sonifies** it: it speeds up
the recording so ground motion normally far below the range of human
hearing (0.5–10 Hz by default, configurable) becomes audible sound.

The DZA network consists of several **sites**, each pairing a surface
broadband sensor with a borehole sensor at ~240 m depth (see
[Station network](#station-network)). The tool defaults to the original
site (`DZA11`/`DZA13`) but any site or combination can be selected with
`--sites`, and a new `map` action prints/plots station locations.

It is plain, auditable scientific tooling, not a black box: standard FDSN
retrieval via [ObsPy](https://docs.obspy.org/), standard instrument-response
removal and bandpass filtering, and a transparent "resample the timebase"
audification method — no proprietary DSP, no hidden steps.

**Author:** Victor Mazon, July 2026
**Scientific Support:** Mike Lindner
**License:** GNU General Public License v3.0 (see [LICENSE](LICENSE))

---

## Pipeline overview

**Operationally** — one command chains four independent stages, each of
which can also be run on its own against a previously saved file:

```mermaid
flowchart LR
    A["FDSN web service<br/>ws.gpi.kit.edu<br/>KB.DZA11 + KB.DZA13 (default site), HH*<br/>6 traces @ 100 Hz"] --> B["fetch<br/>remove_response -&gt; detrend -&gt; bandpass freqmin-freqmax -&gt; taper"]
    B --> C[(".mseed")]
    C --> D["plot<br/>spectrogram (Z) + per-channel waveform panels"]
    C --> E["sonify<br/>pick 1 or all channels<br/>resample timebase x speed-up"]
    D --> F([".png"])
    E --> G([".wav"])
    G --> H(["play"])
```

(Diagram shows the default single-site case; `--sites` can add more
stations to the same pipeline, and `--freqmin`/`--freqmax` override the
default 0.5–10 Hz bandpass, e.g. for the [ULF band](#ultra-low-frequency-ulf-research-band).)

**Conceptually** — sonification here is nothing more than a frequency
shift: the same samples are kept, only the clock they're played back at
changes, moving inaudible ground motion up into the audible range:

```mermaid
flowchart LR
    A["Ground motion<br/>0.5-10 Hz<br/>inaudible / felt, not heard"] -->|"x speed-up<br/>e.g. 100x-200x"| B["Playback audio<br/>50 Hz-2 kHz+<br/>audible"]
```

**Example output** — the `plot` action for a real 24h `--days-back` fetch:
a spectrogram of the vertical channel on top (with the applied bandpass
marked directly on it) and one labeled waveform panel per trace below
(all 6 available channels, each tagged with its depth and sensor model),
with the x-axis automatically scaled to the recording's actual duration
(here, 0-24 hours — matching how long the sonification of the same data
lasts before any `--speed-up`). The right-hand side shows exactly where
in Germany the site is and how deep each sensor sits below ground, so the
whole figure is self-explanatory without cross-referencing the README:

![Example spectrogram, waveform, station map, and sensor-depth plot of a 24h DZA01 recording](docs/images/example_plot.png)

---

## Table of contents

- [What it does](#what-it-does)
- [Station network](#station-network)
- [Sensor differences and recommended filter band](#sensor-differences-and-recommended-filter-band)
- [Ultra-low-frequency (ULF) research band](#ultra-low-frequency-ulf-research-band)
- [How it works](#how-it-works)
  - [1. Fetch](#1-fetch)
  - [2. Plot](#2-plot)
  - [3. Sonify](#3-sonify)
  - [4. Play](#4-play)
  - [5. Map](#5-map)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Command-line reference](#command-line-reference)
- [Recipes](#recipes)
- [Sharing a full sonification batch with a collaborator](#sharing-a-full-sonification-batch-with-a-collaborator)
- [Understanding `--listen-minutes`](#understanding---listen-minutes)
- [Understanding `--hours-back` and `--days-back`](#understanding---hours-back-and---days-back)
- [Understanding the multi-channel output](#understanding-the-multi-channel-output)
- [Output folder layout](#output-folder-layout)
- [Data source and limitations](#data-source-and-limitations)
- [Known limitations / things to be aware of](#known-limitations--things-to-be-aware-of)
- [Contributing](#contributing)
- [Citation](#citation)
- [License](#license)

---

## What it does

The script (`DZA01.py`) is a single entry point with four actions that can
be combined freely on the command line:

| Action   | What it does                                                              |
|----------|----------------------------------------------------------------------------|
| `fetch`  | Downloads a window of raw waveform data from the FDSN web service, removes the instrument response, filters it, and saves it as a `.mseed` file. |
| `plot`   | Renders a dark, pyTREMOR-inspired figure: a spectrogram of the vertical channel on top, and one colored waveform-with-envelope panel per channel below. Saved as `.png`. |
| `sonify` | Speeds up one channel's waveform so it becomes audible, and writes it out as a `.wav` file. |
| `play`   | Plays a `.wav` file straight from the terminal (Windows/macOS/Linux). |
| `map`    | Prints coordinates/depth/status for the selected site(s) and saves a station map plot (self-contained Matplotlib rendering, no `cartopy`/internet required). Does not fetch waveform data. |

By default, running the script with no arguments does `fetch plot sonify`
on a fresh 15-minute window of data for the default site (`DZA11`/`DZA13`).
Use `--sites` to select a different site or combination of sites (see
[Station network](#station-network)).

## Station network

The DZA network is a set of numbered **sites**; each site pairs a surface
broadband sensor (`DZAx1`, at the surface, ~0 m depth) with a borehole
sensor (`DZAx3`, at ~240 m depth ± a few meters) at the same location.
`--sites` selects one or more site numbers (comma-separated), or the
keywords `active` (every currently-streaming site) / `all` (every known
site, including ones not yet installed — they're simply omitted from the
result rather than causing an error).

| Site | Stations | Status | Location (lat, lon, elevation) | Notes |
|---|---|---|---|---|
| 1 | `DZA11` (surface), `DZA13` (borehole) | **active** | 51.32358, 14.24694, 134.6 m | Streaming continuously since project start. Default site if `--sites` is not given. |
| 3 | `DZA31` (surface), `DZA33` (borehole) | **active** | 51.24889, 14.15348, 211.5 m | Streaming, but may be less stable than site 1: solar-powered in a forest region with limited light. |
| 4 | `DZA41` (surface), `DZA43` (borehole) | planned | — | Expected to be installed around mid-to-end of September 2026 (hard to pinpoint exactly). |
| 6 | `DZA61` (surface), `DZA63` (borehole) | planned | — | Expected to be installed around mid-to-end of September 2026 (hard to pinpoint exactly). |

Run `python DZA01.py map` (or `--sites 1,3 map`, etc.) to print this same
information live from the FDSN service (plus the great-circle distance
between sites) and save a plotted map + depth cross-section of the selected
sites to `datasets/maps/station_map.png`. Both are rendered directly with
Matplotlib (no extra dependencies, no `cartopy`/internet-tile requirement) —
the map draws a simplified Germany outline (bundled as
[assets/germany_outline.json](assets/germany_outline.json), so it works
fully offline) with one marker per site, and the depth panel shows the
surface vs. borehole sensor depth for each site side by side. Site markers
use bold, thick-outlined dots (sized for visibility even when two sites are
only ~10 km apart at whole-country zoom) and are colored consistently
across the map, the depth panel, and every waveform panel's label border in
`plot`, so a given site is the same color everywhere. The same two
panels are also embedded automatically on the right-hand side of every
`plot` output (see [Example output](#example-output) above), scoped to only
the site(s) actually present in that recording — no separate `map` run
needed just to know where the data came from.

![Example station map for the two currently-active DZA sites](docs/images/station_map.png)

## Sensor differences and recommended filter band

Not all surface sensors in the network are the same model, which affects
how reliable low-frequency (long-period) signals are at each station:

- **`DZA11`** uses a **Trillium Compact 20s** (20-second eigenperiod).
- **All other surface stations** (`DZA31`, `DZA41`, `DZA61`, ...) use a
  **Trillium Horizon 120s** (120-second eigenperiod).

In practice: amplitudes of signals with periods **below 20s** (i.e. above
0.05 Hz) are less reliable at `DZA11` than at the 120s stations for
small/weak signals — though `DZA11` is perfectly usable for most
applications, including ocean microseism and strong teleseismic events.
It may be less reliable for weak atmospheric-noise-band signals. **Signals
with periods above 20s are reliable at all stations.** The default
bandpass (`--freqmin 0.5 --freqmax 10.0`, i.e. periods of 0.1–2s) is well
within the reliable range of every station in the network regardless of
sensor model.

## Ultra-low-frequency (ULF) research band

The borehole stations are also used for ongoing research into the
**sub-10 mHz band** (periods longer than ~100s) at the
[Black Forest Observatory (BFO)](https://www.bfo-fr.de/), including a new
PhD project — see the related
[publication](https://publikationen.bibliothek.kit.edu/1000140172) for
background. This is a fundamentally different regime from normal
seismology: it overlaps with tidal, atmospheric-pressure, and
gravitational-wave-adjacent signal research.

`--freqmin`/`--freqmax` expose this band to the tool, e.g.:

```bash
# ULF band (<10 mHz), sped up 20000x, over a 10-day window
python DZA01.py --days-back 10 --freqmin 0.0001 --freqmax 0.01 --speed-up 20000 fetch plot sonify
```

**Caution:** instrument-response removal becomes scientifically
unreliable below a sensor's corner frequency — roughly **0.0083 Hz**
(1/120s) for the Trillium Horizon 120s stations, or **0.05 Hz** (1/20s)
for `DZA11`'s Trillium Compact 20s. The script prints a warning (not a
hard block) when `--freqmin` is set below that threshold, since deliberate
exploration below it may still be scientifically or artistically useful —
just treat results in that regime as exploratory rather than validated
without further review from a domain expert.

## How it works

### 1. Fetch

`do_fetch()` uses `obspy.clients.fdsn.Client` to query the KIT/GPI FDSN
web service (`http://ws.gpi.kit.edu`) for exactly one time window:

- **Network:** `KB`
- **Station:** determined by `--sites` (default: site 1, i.e. `DZA11,DZA13`).
  Each selected site contributes its surface and borehole stations, each
  with 3 component channels — so a single-site fetch returns **6
  traces**, and each additional site adds up to 6 more, all sampled at
  100 Hz (see
  [Understanding the multi-channel output](#understanding-the-multi-channel-output)
  and [Station network](#station-network)).
- **Channel:** `HH*` — high-broadband, high-gain seismometer channels.
- **Window length:** `--duration` minutes (default 15), or `--hours-back`
  hours, or a full day per file in `--days-back` batch mode — see
  [Understanding `--hours-back` and `--days-back`](#understanding---hours-back-and---days-back).

FDSN servers typically lag "now" by a few minutes before the latest data
is available. Instead of guessing a fixed delay, the script requests as
close to "now" as possible (`--latency`, default 2 minutes) and, if the
server has no data yet, retries with progressively larger latency
(`--latency-step`, default 2 minutes) up to a ceiling (`--max-latency`,
default 60 minutes). Every fetch prints the exact UTC start/end time and
latency margin actually used, so the window is never opaque:

```
[fetch] Got data from 2026-07-22T15:00:02Z to 2026-07-22T15:01:02Z UTC
(10.0 min latency margin from now, 1.0 min duration). Adjust the margin with --latency.
```

Once data is retrieved, the standard ObsPy processing chain is applied:

1. `remove_response()` — deconvolve the instrument response using the
   station inventory (metadata) fetched alongside the waveform, converting
   raw counts into physical ground-motion units.
2. `detrend("demean")` + `detrend("linear")` — remove DC offset and linear
   drift.
3. `filter("bandpass", freqmin=0.5, freqmax=10.0, corners=4)` — a 4-pole
   Butterworth bandpass, `0.5–10 Hz` by default (overridable with
   `--freqmin`/`--freqmax`, e.g. for the
   [ULF research band](#ultra-low-frequency-ulf-research-band)). The
   default excludes very-long-period drift and high-frequency instrument
   noise for typical seismological work.
4. `taper(0.125)` — a cosine taper on both ends to avoid edge artifacts
   from the filter.

The result is saved as MiniSEED (`.mseed`), the standard seismological
waveform exchange format, under `datasets/mseed/`.

### 2. Plot

`do_plot()` renders a single dark-themed figure per `.mseed` file, designed to be
read "at a glance" for scientific work — every relevant piece of context the
tool knows is shown directly on the figure rather than hidden in the terminal:

- **Title:** network, station codes, exact UTC time range plus the
  equivalent **German civil local time (CET/CEST)**, and the bandpass
  actually applied (`freqmin`–`freqmax`), read back from the `.json` metadata
  sidecar written next to the `.mseed` file at fetch time (falls back to the
  tool's 0.5–10 Hz default for older files saved before this existed). Local
  time is computed with a manual EU-DST rule (last Sunday of March/October) —
  no `zoneinfo`/`tzdata` dependency needed.
- **Top panel:** a spectrogram (`scipy.signal.spectrogram`, `inferno`
  colormap) of the vertical (`*Z`) channel with the most complete data
  available (see below), in dB, clipped to the 5th–99.5th percentile for
  contrast, with **dashed reference lines at `freqmin` and `freqmax`**, plus
  a **solid green line marking the dominant frequency** (highest
  average power across the whole window), so both the applied passband and
  the actual loudest band are visible at a glance.
- **One panel per trace below:** the raw waveform in a component-specific
  color (cyan for Z, orange for N/1, red/pink for E/2), **shaded brighter
  for surface sensors and darker for borehole sensors** (same hue, HLS
  lightness adjusted via `colorsys`) so depth is visually encoded in the
  waveform color itself, not just the text label. A translucent 2-second
  rolling RMS envelope fill sits behind it, and each label's border is
  colored to match that trace's site marker on the map/depth panels (see
  below), so a viewer can trace "this waveform ↔ this map dot ↔ this depth
  column" by color alone. The label itself shows the trace ID, **sampling
  rate**, **depth (surface/borehole) and, for surface sensors, the sensor
  model and eigenperiod**, the **channel orientation (azimuth/dip)** read
  from the station inventory at fetch time (when available), and the
  **peak amplitude and RMS level** for that trace (in m/s, computed from
  the full-resolution data) — so which physical sensor produced each panel,
  how it's oriented, and how strong the signal was is always clear without
  cross-referencing the README.
- **One shared, absolute time axis for every panel** (spectrogram +
  all waveforms), instead of each panel auto-zooming to its own data. Near-
  real-time fetches routinely have channels with slightly different data
  latency/availability; rather than silently stretching a shorter trace to
  fill its panel (which is what previously made the spectrogram look
  "unfinished" — its channel simply had less data than the others), any
  trace covering less than 99% of the full window now visibly stops early/
  late in the *right place*, and its label gains a note like
  `63% of window (ends 67s early)` so the gap is impossible to miss or
  misread as a plotting bug.
- **Time axis scaled to the actual recording length:** seconds, minutes,
  or hours are picked automatically (e.g. a 10h fetch shows a 0–10
  "Time (hours)" axis) so the plot always reads at the same real-world
  duration as the sonification made from the same file (before any
  `--speed-up`).
- **Right-hand side:** the same station-location map and sensor-depth
  cross-section produced by the `map` action (see [below](#5-map)), scoped
  automatically to just the site(s) present in this particular recording —
  so the figure never needs a separate `map` run to know *where* the data
  came from or *how deep* each sensor sits. When more than one site is
  shown, the depth panel's x-axis labels also note the **distance to the
  nearest other active site** (great-circle, via the same `haversine_km`
  helper used by `map`), e.g. `Site 1 (10.5 km to Site 3)`.

For long recordings, waveform panels are downsampled (capped at ~200,000
rendered points per trace) purely for plotting speed/memory — this only
affects the picture, not the sonified audio or the saved `.mseed` data.

Saved as `.png` under `datasets/plot/`.

### 3. Sonify

`do_sonify()` normalizes trace data to the -1..1 range and writes it out
as 16-bit PCM `.wav` audio. By default it picks **one** trace (see
[Understanding the multi-channel output](#understanding-the-multi-channel-output));
pass `--channel all` to sonify **every** trace in the file instead, each
written to its own `.wav` (`{name}_{channel}_{speed}x.wav`).

The "speeding up" is deliberately simple and transparent: the same
samples are kept, but the **declared sample rate of the `.wav` file** is
multiplied by the speed-up factor. A seismometer sampling at 100 Hz
written out at `100 Hz * 200 = 20000 Hz` plays back 200x faster, shifting
all frequency content up by exactly that factor. This is a **linear
time/frequency stretch**, not a pitch-preserving time-stretch — by
design, since the goal is to move inaudible low-frequency ground motion
(< 10 Hz) into the human hearing range, not to disguise the speed change.

```
output_duration = input_duration / speed_up_factor
```

### 4. Play

`do_play()` plays a `.wav` file using the OS's native player:
`winsound` on Windows (Python standard library, no install needed),
`afplay` on macOS, and `paplay`/`aplay`/`ffplay` (whichever is found) on
Linux.

### 5. Map

`do_map()` queries station metadata (not waveform data) for the selected
sites over a wide, fixed time range, so planned-but-not-yet-installed
stations simply don't appear rather than causing an error. For each
station found, it prints the code, depth (surface/borehole), sensor model
(for surface stations), coordinates, elevation, and site status/notes; if
more than one site is selected, it also prints the great-circle distance
between every pair of sites. It then renders two panels with Matplotlib and
saves them to `datasets/maps/station_map.png`:

- **Station location map:** a simplified Germany national border (bundled
  as [assets/germany_outline.json](assets/germany_outline.json), a small
  local asset — no `cartopy`, no internet map tiles, works fully offline)
  with one colored, numbered marker per site (surface and borehole sensors
  share the same coordinates, so they're combined into a single point),
  plus a legend box with the full station details for each site.
- **Sensor depth cross-section:** one column per site showing the surface
  sensor at 0 m and the borehole sensor at ~240 m, connected by a dotted
  line, so the vertical separation between the two sensor types is obvious
  at a glance.

This uses only the dependencies already in `requirements.txt`; no `cartopy`
or internet access is required at plot time. See
[Station network](#station-network) for an example of the printed output.

## Installation

Requires Python 3.9+.

```bash
pip install -r requirements.txt
```

Dependencies: [ObsPy](https://docs.obspy.org/) (FDSN client + seismological
processing), NumPy, SciPy (spectrogram + WAV I/O), Matplotlib (plotting,
including the `map` action's station-location plot — no extra/optional
dependencies needed).

## Quick start

```bash
# Fetch a fresh 15-minute window, plot it, and sonify it (all defaults)
python DZA01.py

# List every .mseed file you've saved so far
python DZA01.py --list

# Re-plot and re-sonify a previously saved file (index 2 from --list)
python DZA01.py --pick 2 plot sonify

# Sonify and immediately play
python DZA01.py sonify play
```

## Command-line reference

```
python DZA01.py [actions ...] [options]
```

**Actions** (positional, any combination, default `fetch plot sonify`):

- `fetch` — download a new window of data
- `plot` — render the styled plot
- `sonify` — produce the `.wav`
- `play` — play a `.wav`
- `map` — print/plot station locations for the selected site(s) (no waveform data fetched)
- `all` — shorthand for `fetch plot sonify`

**Station selection:**

| Option | Default | Description |
|---|---|---|
| `--sites SITES` | site `1` (`active` sites for `map`) | Comma-separated site number(s) (e.g. `1` or `1,3`), or `active` (all currently-streaming sites), or `all` (every known site, including not-yet-installed ones). See [Station network](#station-network). |

**Fetch options:**

| Option | Default | Description |
|---|---|---|
| `--duration MIN` | `15` | Length of the window to fetch, in minutes. |
| `--hours-back HOURS` | — | Convenience alternative to `--duration`, in hours. Overrides `--duration`. |
| `--days-back N` | — | Batch mode: fetch `N` separate ~24h windows, one full day at a time, stepping back from now. Overrides `--duration`/`--hours-back`/`--listen-minutes`; `play` is skipped. See [below](#understanding---hours-back-and---days-back). |
| `--latency MIN` | `2` | How far back from "now" to end the window. Smaller = closer to real time, but the server may not have the data yet. |
| `--max-latency MIN` | `60` | Ceiling the script backs off to if `--latency` is too aggressive. |
| `--latency-step MIN` | `2` | How much to increase latency per retry when backing off. |
| `--freqmin HZ` | `0.5` | Lower bandpass corner. Lower this for long-period work, e.g. the [ULF research band](#ultra-low-frequency-ulf-research-band) (`--freqmin 0.0001`). |
| `--freqmax HZ` | `10.0` | Upper bandpass corner. |

**Sonify options:**

| Option | Default | Description |
|---|---|---|
| `--speed-up N` | `200` | Playback speed multiplier. 15 min of data at 200x ≈ 4.5 s of audio. |
| `--channel SUBSTRING` \| `all` | first trace | Which of the (typically 6) channels to sonify, matched as a case-insensitive substring against the trace ID, e.g. `DZA11`, `HHZ`, `DZA13.00.HH1`. Pass `all` to sonify every channel into its own `.wav`. |
| `--listen-minutes MIN` | — | Target sonification length; see [below](#understanding---listen-minutes). |

**File selection (skip fetching):**

| Option | Description |
|---|---|
| `--list` | List saved `.mseed` files with index numbers, then exit. |
| `--pick INDEX` | Use the file at `INDEX` from `--list` instead of fetching. |
| `--file PATH` | Use a specific `.mseed` file path instead of fetching. |

Run `python DZA01.py -h` at any time for the full built-in help text.

## Recipes

```bash
# Near-real-time: try to grab data as close to "now" as possible
python DZA01.py --latency 0.5 --max-latency 60

# Sonify a different channel from the ones already fetched
python DZA01.py --pick 0 sonify --channel HH2

# A slower, longer sonification (20x instead of the default 200x)
python DZA01.py sonify --speed-up 20

# A 15-minute listening piece at the 100x audible floor (fetches 25h of raw data)
python DZA01.py --listen-minutes 15 fetch plot sonify play

# Grab the last 24 hours in one window instead of computing --duration in minutes
python DZA01.py --hours-back 24 --speed-up 100 fetch plot sonify

# Sonify every channel in a saved file, not just the default one
python DZA01.py --pick 0 sonify --channel all

# Build a 10-day dataset of 24h files (fetch + plot + sonify each day, no play)
python DZA01.py --days-back 10 --speed-up 100 fetch plot sonify

# Print/plot station locations for the two active sites
python DZA01.py --sites 1,3 map

# Fetch two sites together instead of just the default one
python DZA01.py --sites 1,3 fetch plot sonify

# Ultra-low-frequency (ULF) band (<10 mHz) over a 10-day window, sped up a lot
python DZA01.py --days-back 10 --freqmin 0.0001 --freqmax 0.01 --speed-up 20000 fetch plot sonify

# Full "every possible sound source" batch: both active sites, every channel,
# one combined plot + one .wav per channel -- see docs/sonification_batch_2026-07-27.md
# for a worked example index of exactly this command's output, ready to share
# with a collaborating scientist.
python DZA01.py --sites 1,3 --hours-back 24 --speed-up 100 --channel all fetch plot sonify
```

## Sharing a full sonification batch with a collaborator

[`docs/sonification_batch_2026-07-27.md`](docs/sonification_batch_2026-07-27.md)
is a worked example of preparing every currently possible "sound source" (all
channels, both active sites, one 24h window) as a single deliverable: the
combined plot, the raw `.mseed`, and one `.wav` per channel, with a table
explaining what each file is (site, depth, sensor model, orientation,
duration, and any data-coverage caveats) plus a few suggested comparisons
(surface vs. borehole, sensor model differences, component orientation).
Regenerate it any time with the command shown in that file, swapping
`--sites 1,3` for `--sites all` once more sites are installed.


## Understanding `--listen-minutes`

`--speed-up` alone tells you *how much* to compress a recording, but not
*how long* the result will be — that depends on how much raw data you
fetched. `--listen-minutes` flips the question around: **"I want the
final audio to last exactly N minutes — figure out the rest for me."**
It behaves differently depending on whether you're fetching new data or
reusing a saved file:

- **When fetching** (`fetch` is one of the actions): the amount of raw
  data to download is computed as `--listen-minutes * --speed-up`. If you
  don't also pass `--speed-up` explicitly, it defaults to **100x** — the
  `MIN_LISTEN_SPEED_UP` floor, below which bandpassed (0.5–10 Hz) ground
  motion is barely more than a faint, inaudible rumble. So
  `--listen-minutes 15` alone fetches **25 hours** of raw data by
  default, in order to actually be worth listening to. Pass `--speed-up`
  explicitly to request a *different* amount of raw data on purpose, e.g.
  `--listen-minutes 15 --speed-up 20` fetches only 5 hours and compresses
  it into a 15-minute piece (quieter/less audible, but faster to fetch).
  The script prints the raw duration it's about to fetch, and warns if
  it's a large request (> 6 hours), since long fetches are slower and
  produce larger files.
- **When reusing an existing file** (`--pick`/`--file`, no `fetch`):
  `--speed-up` is instead **computed for you** as
  `(existing recording length) / --listen-minutes`, so a saved 25-hour
  file becomes a 15-minute sonification at 100x. If the saved file is too
  short to reach the 100x floor for the requested listen length (e.g. a
  15-minute file can only reach ~1x for a 15-minute listen), the script
  prints a warning and tells you how much raw data you'd need to fetch to
  reach 100x instead.


## Understanding `--hours-back` and `--days-back`

Both are alternatives to specifying `--duration` in minutes, aimed at
pytremor-style "lookback window" workflows:

- **`--hours-back HOURS`** fetches a single window of that many hours,
  ending `--latency` minutes before now. It overrides `--duration` but
  behaves exactly like a normal `fetch` otherwise — one `.mseed` file,
  one optional `plot`, one optional `sonify`/`play`.
- **`--days-back N`** is batch mode: it runs `fetch` **N times**, each
  time for a distinct ~24h window one full day further back than the
  last (day 0 = most recent day before latency, day 1 = the day before
  that, and so on), producing **N separate `.mseed` files**. If `plot`
  and/or `sonify` are also requested, each of the N files is processed
  in turn. A failed day (e.g. no data available that far back) is
  skipped with a warning rather than aborting the whole batch. `play` is
  always skipped in this mode — it exists to build a dataset, not to
  listen interactively. `--days-back` overrides `--duration`,
  `--hours-back`, and `--listen-minutes` if any are also given.

```bash
# 10 days of 24h files, plotted and sonified at 100x — a batch dataset
python DZA01.py --days-back 10 --speed-up 100 fetch plot sonify
```

## Understanding the multi-channel output

Each selected site contributes **two** physical sub-stations — a surface
station (e.g. `DZA11`) and a borehole station (e.g. `DZA13`) — each
recording **three** components (a vertical `Z` and two horizontals, named
`N`/`E` or `1`/`2` depending on the sub-station). That means **every site
adds 6 traces**, all at 100 Hz. The default (`--sites 1`, i.e. `DZA1*`)
returns 6 traces; `--sites 1,3` returns 12; and so on:

```
KB.DZA11.00.HHZ   KB.DZA11.00.HHN   KB.DZA11.00.HHE
KB.DZA13.00.HHZ   KB.DZA13.00.HH1   KB.DZA13.00.HH2
```

`plot` shows **every** trace in the file. `sonify` only ever turns **one**
trace into audio at a time (audio is inherently single-channel here) — by
default the first trace in the file, or whichever one matches
`--channel`. The script always prints which channel it picked and which
ones it's ignoring, so this is never silent/surprising:

```
[sonify] 6 channels available; using 'KB.DZA11.00.HHZ'. Ignoring: KB.DZA11.00.HHN,
KB.DZA11.00.HHE, KB.DZA13.00.HHZ, KB.DZA13.00.HH1, KB.DZA13.00.HH2. Use --channel to
pick a different one.
```

## Output folder layout

Created automatically next to the script on first run:

```
datasets/
├── mseed/           raw, processed waveform data (.mseed) + fetch metadata (.json)
├── plot/             spectrogram + waveform plots (.png)
├── sonifications/    sonified audio (.wav)
└── maps/             station map plots from the `map` action (.png)
```

Saved data is not tracked in git (see `.gitignore`) — only the folder
structure is, via `.gitkeep` placeholders. Re-run `fetch` to regenerate
data locally; the FDSN archive is the source of truth.

## Data source and limitations

- Data is retrieved live from the KIT/GPI FDSN web service
  (`http://ws.gpi.kit.edu`), a public academic seismic data service. Note
  that this endpoint is plain HTTP (not HTTPS) — that is the server's own
  configuration, not something this script controls, so treat retrieved
  data as authentic-but-unverified in transit (there is no cryptographic
  integrity check on the wire). This is standard practice for many
  academic FDSN endpoints and is not a concern for casual/artistic use,
  but should be kept in mind for any use where data provenance matters.
- Availability and latency of "real-time" data depend entirely on the
  station and network operator; this script cannot guarantee data exists
  for any given time window (see `--max-latency` backoff behavior).
- This is a small utility script, not a package — there are no automated
  tests. It has been manually exercised for the CLI paths described in
  this README (short and multi-hour fetches, single/multi-channel
  sonification, plotting, latency backoff, and `--listen-minutes` in both
  directions).

## Known limitations / things to be aware of

- Sonification always converts to mono 16-bit PCM; if the underlying
  ground motion is silent/near-zero for a whole window, the resulting
  audio will be effectively silent too (this is physically accurate, not
  a bug).
- Very large `--listen-minutes`/`--speed-up` combinations can trigger
  multi-hour FDSN downloads; the script warns above 6 hours of raw data
  but does not hard-block larger requests, since that may be a deliberate
  choice for some listening pieces.
- `--speed-up`, `--duration`, and `--listen-minutes` are validated to be
  positive; `--latency`/`--max-latency` are validated to be non-negative
  and consistent with each other. Malformed FDSN responses or network
  errors surface as Python exceptions with the underlying error message.

## Contributing

Issues and pull requests are welcome — this is meant to be a small,
readable reference implementation for turning seismic waveform data into
sound, not a polished framework. Feel free to fork it, adapt the station/
network/channel configuration at the top of `DZA01.py` for a different
seismic station, or extend the sonification approach (e.g. pitch-preserving
time-stretch, multi-channel stereo mixes, etc.).

## Citation

If this tool is useful in an academic or creative context, a citation or
acknowledgment along these lines is appreciated:

> Mazon, V. (2026). *DZA01 Seismic Sonification Toolkit* [Software].

## License

Licensed under the **GNU General Public License v3.0** — see
[LICENSE](LICENSE) for the full text. In short: you are free to use,
study, modify, and redistribute this software (including commercially),
provided derivative works are also released under the GPL-3.0 and you
preserve the copyright/license notices.

Copyright (C) 2026 Victor Mazon.
