# DZA network — full sonification batch (2026-07-27/28, 24h, both active sites)

Prepared for the collaborating scientist (KIT/GPI) to share every currently
possible "sound source" from the DZA network in one batch: all channels at
both currently-active sites, sonified from the same 24h recording so they can
be compared directly against one another.

This is the current/clean batch — the `datasets/sonifications/` folder was
wiped of older runs and one-off test files first, so this is the only run in
it and there's no ambiguity about which files go together.

## Source recording

| | |
|---|---|
| Command | `python DZA01.py --sites 1,3 --hours-back 24 --speed-up 100 --channel all fetch plot sonify` |
| Time window (UTC) | 2026-07-27T10:00:07.16 – 2026-07-28T10:00:07.16 |
| Time window (local, Germany) | 2026-07-27 12:00 – 2026-07-28 12:00 CEST |
| Bandpass applied | 0.5–10 Hz |
| Raw waveform data | `datasets/mseed/DZA_sites-1-3_2026-07-27-10-00.mseed` (+ `.json` metadata sidecar) |
| Combined plot (spectrogram + all 14 waveform panels + map + depth) | `datasets/plot/DZA_sites-1-3_2026-07-27-10-00.png` |
| Speed-up | 100x (a 24h recording becomes ~14.4 min of audio) |

Open the `.png` first — every waveform panel is labeled (in a single line, or
two where there's a channel orientation and/or a data-coverage gap to note)
with the exact same metadata as the table below: sensor depth/model, channel
orientation, peak/RMS amplitude, and any data-coverage gap, color-matched to
the site markers on the station map and depth cross-section on the right.

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
| 1 | `KB-DZA11-00-HHZ_100x.wav` | KB.DZA11.00.HHZ | surface | 0°/-90° (vertical) | 86364.1 s | 863.6 s | full |
| 2 | `KB-DZA11-00-HHN_100x.wav` | KB.DZA11.00.HHN | surface | 0°/0° (north) | 86346.7 s | 863.5 s | full |
| 3 | `KB-DZA11-00-HHE_100x.wav` | KB.DZA11.00.HHE | surface | 90°/0° (east) | 86342.9 s | 863.4 s | full |
| 4 | `KB-DZA13-00-HHZ_100x.wav` | KB.DZA13.00.HHZ | borehole ~240 m | 0°/-90° (vertical) | 86400.0 s | 864.0 s | full |
| 5 | `KB-DZA13-00-HH1_100x.wav` | KB.DZA13.00.HH1 | borehole ~240 m | 0°/0° | 86400.0 s | 864.0 s | full |
| 6 | `KB-DZA13-00-HH2_100x.wav` | KB.DZA13.00.HH2 | borehole ~240 m | 90°/0° | 86400.0 s | 864.0 s | full |
| 7 | `KB-DZA31-00-HHZ_100x.wav` | KB.DZA31.00.HHZ | surface | 0°/-90° (vertical) | 85382.1 s | 853.8 s | 99% (ends 1018s early) |
| 8 | `KB-DZA31-00-HHN_100x.wav` | KB.DZA31.00.HHN | surface | 0°/0° (north) | 85381.0 s | 853.8 s | 99% (ends 1019s early) |
| 9 | `KB-DZA31-00-HHE_100x.wav` | KB.DZA31.00.HHE (segment 1) | surface | 90°/0° (east) | 81876.0 s | 818.8 s | 95% (early segment, ends 4524s before gap) |
| 10 | `KB-DZA31-00-HHE_seg2_100x.wav` | KB.DZA31.00.HHE (segment 2) | surface | 90°/0° (east) | 3506.6 s | 35.1 s | 4% (short remainder after gap, starts 81874s late) |
| 11 | `KB-DZA33-00-HHZ_100x.wav` | KB.DZA33.00.HHZ | borehole ~240 m | 0°/-90° (vertical) | 83585.6 s | 835.9 s | 97% (ends 2814s early) |
| 12 | `KB-DZA33-00-HH1_100x.wav` | KB.DZA33.00.HH1 | borehole ~240 m | 0°/0° | 83584.6 s | 835.8 s | 97% (ends 2815s early) |
| 13 | `KB-DZA33-00-HH2_100x.wav` | KB.DZA33.00.HH2 (segment 1) | borehole ~240 m | 90°/0° | 16037.9 s | 160.4 s | 19% (early segment only, before gap) |
| 14 | `KB-DZA33-00-HH2_seg2_100x.wav` | KB.DZA33.00.HH2 (segment 2) | borehole ~240 m | 90°/0° | 67550.2 s | 675.5 s | 78% (remainder after gap) |

(Every filename above is prefixed with `DZA_sites-1-3_2026-07-27-10-00_`.)

## Listening suggestions

- **Surface vs. borehole, same site:** compare #1 (`DZA11 Z`, surface) with
  #4 (`DZA13 Z`, borehole) — the borehole sensor should sound noticeably
  quieter/cleaner (less wind/traffic noise coupling at depth).
- **Sensor model difference:** DZA11 (20s eigenperiod) vs. DZA31 (120s
  eigenperiod) are the two surface sensor models in the network — compare
  #1 and #9/#10 for how the different sensor response shapes the sound.
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
