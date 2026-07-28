# DZA network — full sonification batch (2026-07-27/28, 24h, both active sites)

Prepared for the collaborating scientist (KIT/GPI) to share every currently
possible "sound source" from the DZA network in one batch: all channels at
both currently-active sites, sonified from the same 24h recording so they can
be compared directly against one another.

This is the current/clean batch — the `datasets/sonifications/` folder was
wiped of older runs and one-off test files first, so this is the only run in
it and there's no ambiguity about which files go together. It was also
re-fetched after fixing a bug where the short segments before/after a
mid-recording gap (see below) spent up to 25% of their own duration fading
in/out, making them sound weak/broken compared to the full-length channels --
the taper is now capped to a fixed 20 real-world seconds per side, regardless
of segment length.

## Source recording

| | |
|---|---|
| Command | `python DZA01.py --sites 1,3 --hours-back 24 --speed-up 100 --channel all fetch plot sonify` |
| Time window (UTC) | 2026-07-27T12:42:39.11 – 2026-07-28T12:42:39.11 |
| Time window (local, Germany) | 2026-07-27 14:42 – 2026-07-28 14:42 CEST |
| Bandpass applied | 0.5–10 Hz |
| Raw waveform data | `datasets/mseed/DZA_sites-1-3_2026-07-27-12-42.mseed` (+ `.json` metadata sidecar) |
| Combined plot (spectrogram + all 15 waveform panels + map + depth) | `datasets/plot/DZA_sites-1-3_2026-07-27-12-42.png` |
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
| 1 | `KB-DZA11-00-HHZ_100x.wav` | KB.DZA11.00.HHZ | surface | 0°/-90° (vertical) | 86277.8 s | 862.8 s | full |
| 2 | `KB-DZA11-00-HHN_100x.wav` | KB.DZA11.00.HHN | surface | 0°/0° (north) | 86262.9 s | 862.6 s | full |
| 3 | `KB-DZA11-00-HHE_100x.wav` | KB.DZA11.00.HHE | surface | 90°/0° (east) | 86260.1 s | 862.6 s | full |
| 4 | `KB-DZA13-00-HHZ_100x.wav` | KB.DZA13.00.HHZ | borehole ~240 m | 0°/-90° (vertical) | 86400.0 s | 864.0 s | full |
| 5 | `KB-DZA13-00-HH1_100x.wav` | KB.DZA13.00.HH1 | borehole ~240 m | 0°/0° | 86400.0 s | 864.0 s | full |
| 6 | `KB-DZA13-00-HH2_100x.wav` | KB.DZA13.00.HH2 | borehole ~240 m | 90°/0° | 86400.0 s | 864.0 s | full |
| 7 | `KB-DZA31-00-HHZ_100x.wav` | KB.DZA31.00.HHZ | surface | 0°/-90° (vertical) | 84440.1 s | 844.4 s | 98% (ends 1960s early) |
| 8 | `KB-DZA31-00-HHN_100x.wav` | KB.DZA31.00.HHN | surface | 0°/0° (north) | 84438.3 s | 844.4 s | 98% (ends 1962s early) |
| 9 | `KB-DZA31-00-HHE_100x.wav` | KB.DZA31.00.HHE (segment 1) | surface | 90°/0° (east) | 72124.0 s | 721.2 s | 83% (first segment, ends 14276s before first gap) |
| 10 | `KB-DZA31-00-HHE_seg2_100x.wav` | KB.DZA31.00.HHE (segment 2) | surface | 90°/0° (east) | 8661.1 s | 86.6 s | 10% (between two gaps, starts 72122s late) |
| 11 | `KB-DZA31-00-HHE_seg3_100x.wav` | KB.DZA31.00.HHE (segment 3) | surface | 90°/0° (east) | 3658.9 s | 36.6 s | 4% (short remainder after second gap, starts 80781s late) |
| 12 | `KB-DZA33-00-HHZ_100x.wav` | KB.DZA33.00.HHZ | borehole ~240 m | 0°/-90° (vertical) | 84909.2 s | 849.1 s | 98% (ends 1491s early) |
| 13 | `KB-DZA33-00-HH1_100x.wav` | KB.DZA33.00.HH1 | borehole ~240 m | 0°/0° | 84910.2 s | 849.1 s | 98% (ends 1490s early) |
| 14 | `KB-DZA33-00-HH2_100x.wav` | KB.DZA33.00.HH2 (segment 1) | borehole ~240 m | 90°/0° | 6285.9 s | 62.9 s | 7% (early segment only, before gap) |
| 15 | `KB-DZA33-00-HH2_seg2_100x.wav` | KB.DZA33.00.HH2 (segment 2) | borehole ~240 m | 90°/0° | 78625.7 s | 786.3 s | 91% (remainder after gap) |

(Every filename above is prefixed with `DZA_sites-1-3_2026-07-27-12-42_`.)

Note: DZA31 HHE has *two* gaps in this particular 24h window (three
segments) -- a slightly less stable stretch of telemetry from the
solar-powered Site 3 station, not a processing artifact. Segments #10/#11 are
short but should now sound like a clean excerpt from start to finish (no
audible fade eating into the content) thanks to the taper fix above.

## Listening suggestions

- **Surface vs. borehole, same site:** compare #1 (`DZA11 Z`, surface) with
  #4 (`DZA13 Z`, borehole) — the borehole sensor should sound noticeably
  quieter/cleaner (less wind/traffic noise coupling at depth).
- **Sensor model difference:** DZA11 (20s eigenperiod) vs. DZA31 (120s
  eigenperiod) are the two surface sensor models in the network — compare
  #1 and #9/#10/#11 for how the different sensor response shapes the sound.
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
