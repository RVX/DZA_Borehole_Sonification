# DZA network — full sonification batch (2026-07-27/28, 24h, both active sites)

Prepared for the collaborating scientist (KIT/GPI) to share every currently
possible "sound source" from the DZA network in one batch: all channels at
both currently-active sites, sonified from the same 24h recording so they can
be compared directly against one another.

## Source recording

| | |
|---|---|
| Command | `python DZA01.py --sites 1,3 --hours-back 24 --speed-up 100 --channel all fetch plot sonify` |
| Time window (UTC) | 2026-07-27T08:30:51 – 2026-07-28T08:30:51 |
| Time window (local, Germany) | 2026-07-27 10:30 – 2026-07-28 10:30 CEST |
| Bandpass applied | 0.5–10 Hz |
| Raw waveform data | `datasets/mseed/DZA_sites-1-3_2026-07-27-08-30.mseed` (+ `.json` metadata sidecar) |
| Combined plot (spectrogram + all 15 waveform panels + map + depth) | `datasets/plot/DZA_sites-1-3_2026-07-27-08-30.png` |
| Speed-up | 100x (a 24h recording becomes ~14.4 min of audio) |

Open the `.png` first — every waveform panel is labeled with the exact same
metadata as the table below (sensor depth/model, channel orientation,
peak/RMS amplitude, and any data-coverage gap), color-matched to the site
markers on the station map and depth cross-section on the right.

## Sites covered

| Site | Surface sensor | Borehole sensor | Distance to other site |
|---|---|---|---|
| 1 | DZA11 — Trillium Compact 20s | DZA13 — ~240 m depth | 10.5 km to Site 3 |
| 3 | DZA31 — Trillium Horizon 120s | DZA33 — ~240 m depth | 10.5 km to Site 1 |

Sites 4 and 6 are not included — not yet installed (planned ~September 2026).

## Per-channel .wav files

All files are in `datasets/sonifications/`, 16-bit PCM WAV, sample rate =
trace sampling rate (100 Hz) × 100 = 10,000 Hz.

A channel that has a gap partway through the 24h window (not just a
start/end truncation, but a genuine mid-recording drop-out) is split by
ObsPy into more than one file, suffixed `_segN`, in chronological order —
this is normal for near-real-time telemetry and does not indicate a
problem with the sensor.

| # | File (suffix after the shared prefix) | Trace ID | Depth | Orientation (az/dip) | Input duration | Output duration | Coverage note |
|---|---|---|---|---|---|---|---|
| 1 | `KB-DZA11-00-HHZ_100x.wav` | KB.DZA11.00.HHZ | surface | 0°/-90° (vertical) | 86279 s | 862.8 s | full |
| 2 | `KB-DZA11-00-HHN_100x.wav` | KB.DZA11.00.HHN | surface | 0°/0° (north) | 86265 s | 862.6 s | full |
| 3 | `KB-DZA11-00-HHE_100x.wav` | KB.DZA11.00.HHE | surface | 90°/0° (east) | 86261 s | 862.6 s | full |
| 4 | `KB-DZA13-00-HHZ_100x.wav` | KB.DZA13.00.HHZ | borehole ~240 m | 0°/-90° (vertical) | 86400 s | 864.0 s | full |
| 5 | `KB-DZA13-00-HH1_100x.wav` | KB.DZA13.00.HH1 | borehole ~240 m | 0°/0° | 86400 s | 864.0 s | full |
| 6 | `KB-DZA13-00-HH2_100x.wav` | KB.DZA13.00.HH2 (segment 1) | borehole ~240 m | 90°/0° | 3604 s | 36.0 s | early segment only, before gap |
| 7 | `KB-DZA13-00-HH2_seg2_100x.wav` | KB.DZA13.00.HH2 (segment 2) | borehole ~240 m | 90°/0° | 82800 s | 828.0 s | remainder after gap |
| 8 | `KB-DZA31-00-HHZ_100x.wav` | KB.DZA31.00.HHZ (segment 1) | surface | 0°/-90° (vertical) | 3897 s | 39.0 s | early segment only, before gap |
| 9 | `KB-DZA31-00-HHZ_seg2_100x.wav` | KB.DZA31.00.HHZ (segment 2) | surface | 0°/-90° (vertical) | 81895 s | 818.9 s | remainder after gap |
| 10 | `KB-DZA31-00-HHN_100x.wav` | KB.DZA31.00.HHN | surface | 0°/0° (north) | 85788 s | 857.9 s | full |
| 11 | `KB-DZA31-00-HHE_100x.wav` | KB.DZA31.00.HHE | surface | 90°/0° (east) | 85790 s | 857.9 s | full |
| 12 | `KB-DZA33-00-HHZ_100x.wav` | KB.DZA33.00.HHZ | borehole ~240 m | 0°/-90° (vertical) | 85439 s | 854.4 s | 99% (ends early) |
| 13 | `KB-DZA33-00-HH1_100x.wav` | KB.DZA33.00.HH1 | borehole ~240 m | 0°/0° | 85439 s | 854.4 s | 99% (ends early) |
| 14 | `KB-DZA33-00-HH2_100x.wav` | KB.DZA33.00.HH2 (segment 1) | borehole ~240 m | 90°/0° | 21394 s | 213.9 s | early segment only, before gap |
| 15 | `KB-DZA33-00-HH2_seg2_100x.wav` | KB.DZA33.00.HH2 (segment 2) | borehole ~240 m | 90°/0° | 64048 s | 640.5 s | remainder after gap |

(Every filename above is prefixed with `DZA_sites-1-3_2026-07-27-08-30_`.)

## Listening suggestions

- **Surface vs. borehole, same site:** compare #1 (`DZA11 Z`, surface) with
  #4 (`DZA13 Z`, borehole) — the borehole sensor should sound noticeably
  quieter/cleaner (less wind/traffic noise coupling at depth).
- **Sensor model difference:** DZA11 (20s eigenperiod) vs. DZA31 (120s
  eigenperiod) are the two surface sensor models in the network — compare
  #1 and #8/#9 for how the different sensor response shapes the sound.
  Amplitudes are absolute (m/s), so relative loudness between files is
  meaningful, not just an artifact of normalization.
- **Component/orientation:** Z (vertical) vs. N/E (horizontal) panels on the
  same station can be compared directly since they're the same time window.

## Regenerating or extending this batch

```
python DZA01.py --sites 1,3 --hours-back 24 --speed-up 100 --channel all fetch plot sonify
```

Swap `--sites 1,3` for `--sites all` once sites 4/6 are installed to include
them automatically. Change `--speed-up` to compare bands at a different
compression ratio (e.g. `--speed-up 20000` for the ultra-low-frequency band,
combined with `--freqmin`/`--freqmax`, see the main README).
