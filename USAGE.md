# SeismoFK — Usage Manual

A step-by-step walkthrough of **SeismoFK**, the desktop tool for **frequency–wavenumber (FK) array analysis** of infrasound and seismic array data.

This manual covers the everyday workflow in the GUI, the three tools reachable from the bottom of the main window (XML Creator, Event Database, Plot Spectrogram), the standalone command-line scripts, and the output files SeismoFK produces. For background on what FK array analysis is and for installation instructions, see [`README.md`](README.md).

---

## Contents

1. [Introduction](#1-introduction)
2. [Getting started](#2-getting-started)
3. [Step-by-step FK analysis workflow](#3-step-by-step-fk-analysis-workflow)
4. [Plot Spectrogram](#4-plot-spectrogram)
5. [XML Creator / Editor](#5-xml-creator--editor)
6. [Event Database](#6-event-database)
7. [Command-line tools](#7-command-line-tools)
8. [Output files](#8-output-files)
9. [Troubleshooting and tips](#9-troubleshooting-and-tips)

---

## 1. Introduction

When an infrasound or seismic signal crosses an **array** of closely spaced sensors, it reaches each element at slightly different times. **FK (frequency–wavenumber) array analysis** measures those time delays to estimate two key properties of the incoming wavefield:

- **Back-azimuth** — the direction the signal arrives *from*.
- **Apparent (trace) velocity / slowness** — how fast the wavefront sweeps across the array.

SeismoFK runs FK analysis in sliding time windows and, for each window, reports the **semblance** (a 0–1 coherence measure), the **Fisher ratio**, the **FK power**, the **back-azimuth**, and the **apparent velocity**. It also forms a **delay-and-sum beam** steered to the dominant detected direction, and compares the measured back-azimuth against the back-azimuth expected from a user-supplied source location.

This manual assumes SeismoFK is already installed (see the README's *Installation* section).

---

## 2. Getting started

### 2.1 Launching the application

From the project directory, with your Python environment active:

```bash
python Infra_Analysis.py
```

This opens the main window, titled **SeismoFK — Infrasound FK Array Analysis**.

![Main analysis window](screenshots/01_main_window.png)

The main window is organised top-to-bottom into:

1. **File Selection** — choose a MiniSEED waveform file and a station inventory (XML).
2. **Waveform Preview** — a multi-trace plot of the loaded data; click on it to pick the analysis start time.
3. **Analysis Parameters** — the FK frequency band, window settings, start time, duration, event name, expected source location, and the optional origin-time / celerity controls.
4. **Status / Progress / Run** — a status line, a progress bar, and the green **Run FK Analysis** button.
5. **Tools row** — three buttons: **XML Creator / Editor**, **Event Database**, and **Plot Spectrogram**.

### 2.2 Preparing station metadata (XML inventories)

FK analysis needs **station metadata** — the coordinates and instrument response of every sensor in the array — supplied as **StationXML** files. As noted in the README's *Station metadata* section, the repository does **not** ship station XML files; you supply your own. There are three ways to obtain them:

- Run `convert_ims.py` against an FDSN StationXML to generate per-station files (see [Section 7](#7-command-line-tools)).
- Build station arrays interactively with the **XML Creator** (see [Section 5](#5-xml-creator--editor)).
- Use the **+ Add XML** button in the main window to import an existing StationXML file.

SeismoFK discovers inventories from the user data directory and, when running
from a checkout, two directories next to `Infra_Analysis.py`:

- **`~/.seismofk/XML/`** — inventories imported with **+ Add XML** or saved by the XML Creator. Override the root with `SEISMOFK_DATA_DIR`.

- **`XML_IM/`** — intended for IMS / multi-station array inventories. If it contains any `.xml` files, a single **"★ All IMS Stations (XML_IM/)"** entry appears first in the *Inventory* dropdown; selecting it merges *every* `.xml` file in that directory into one combined inventory. Each individual file is also listed separately, tagged `[IMS]`.
- **`XML/`** — for individual / custom station files. Each `.xml` file is listed as its own entry.

Click **⟳ Refresh** at any time to rescan these directories — for example, after creating a new XML with the XML Creator.
Select the inventory for the waveform's array explicitly. Use **All IMS
Stations** only when a waveform contains channels from multiple IMS arrays.
If an inventory has several metadata epochs for one channel, SeismoFK uses the
waveform trace time for coordinate lookup; resolve overlapping, conflicting
epochs in the source StationXML.

---

## 3. Step-by-step FK analysis workflow

This is the core path through SeismoFK: load waveforms → select an inventory → set parameters → preview and pick a start time → run the analysis → review results.

Use **Check data readiness** after loading a waveform. It reports missing
metadata, unsuitable array geometry, invalid frequency bands, gaps in common
time coverage, and channels that cannot provide pressure calibration. Resolve
blocking issues before interpreting FK or PMCC results.

### Step 1 — Load MiniSEED waveform data

In the **File Selection** panel, click **Browse** next to *MiniSEED*. Select **one or more** MiniSEED files (`.mseed`, `.msd`, `.ms`). If you select several files, their traces are read and **merged automatically** into a single stream, and the field shows how many files were merged. Gaps are filled and masked samples are replaced with zeros so processing has a continuous record.

Once data is loaded, the **Waveform Preview** populates with one normalised panel per trace, labelled with `network.station.channel`, the sampling rate, and the trace start time.

### Step 2 — Select a station inventory

In the **Inventory** dropdown, choose the StationXML entry that matches the array in your waveform data:

- For an IMS array, pick the relevant `[IMS]` entry, or **★ All IMS Stations (XML_IM/)** to merge everything in `XML_IM/`.
- For a custom array, pick its file from the `XML/` group.

If the inventory you need is not listed, click **+ Add XML** to copy a StationXML file into `~/.seismofk/XML/`, or build one with the XML Creator and then click **⟳ Refresh**. Set `SEISMOFK_DATA_DIR` to change the user data location. Checkout `XML/` and `XML_IM/` files remain discoverable.

> The inventory must contain coordinates and instrument response for every sensor present in the loaded waveform. See [Section 9](#9-troubleshooting-and-tips) for what SeismoFK requires of the response.

### Step 3 — Preview the waveforms and pick an analysis start time

The **Waveform Preview** is interactive:

- **Left-click anywhere on the preview plot** to set the FK analysis **start time**. A red dashed line marks the pick, the *Pick:* label shows the absolute time, and the **Start Time** parameter field below updates to match.
- Click **✖ Clear Pick** to remove the pick and reset the start time.
- Click **⤢ Open in Window** to open the larger **Waveform Viewer** dialog. There you can apply a band-pass filter (set *Low* / *High* Hz and click **Apply Filter**, or **Reset** to undo), left-click to pick a start time on the bigger plot, then click **✔ Confirm** to send the pick back to the main window.

> **Note:** picking a start time is a convenience for setting the *Start Time* parameter. You can also set *Start Time* manually (Step 4).

### Step 4 — Set the analysis parameters

In the **Analysis Parameters** panel:

**Row 1 — FK and window settings**

| Field | Meaning | Default |
|---|---|---|
| **Min Freq (Hz)** | Lower corner of the band-pass filter / FK band. | 0.5 |
| **Max Freq (Hz)** | Upper corner of the band-pass filter / FK band. | 6.0 |
| **Window (s)** | Length of each sliding FK analysis window. | 20.0 |
| **Window step** | Fraction of window length between FK windows (0.01–0.99); `0.1` means a 10% step and 90% overlap. | 0.1 |
| **Semb. Threshold** | Semblance value at/above which a window counts as a *detection* in the results plots. | 0.30 |

**Row 2 — event and geometry**

| Field | Meaning | Default |
|---|---|---|
| **Start Time** | Analysis start time (`yyyy-MM-dd HH:mm:ss`). Set by a preview pick, or typed manually. | current time |
| **Duration** | Total length of data analysed, in seconds, from the start time. | 900 |
| **Event Name** | Label used for output filenames and the database record. | `Event` |
| **Lat** / **Lon** | Expected **source** latitude / longitude — used to compute the *expected back-azimuth* shown for comparison in the results. | 0 / 0 |

**Row 3 — optional event physics**

- **Origin Time** — tick the checkbox and set a known event origin time to enable this.
- **Celerity (m/s)** — tick the checkbox and set the infrasound propagation speed (typical range ~220–340 m/s).

When **both** Origin Time and Celerity are set, SeismoFK overlays the **expected infrasound arrival time** on the beam-waveform panel of the results figure, using the source-to-array distance and the celerity. These two controls are optional; leave them unchecked for a plain FK run.

**Row 4 — uncertainty (v1.2.1)**

- **Estimate uncertainty (bootstrap)** — tick this to bootstrap the detected back-azimuth and apparent velocity over the detection windows. The results title gains a line such as `Bootstrap 95% CI (104 windows, blocks of 10): baz 84.2° [80.5, 86.5]  app. vel 348 m/s [334.1, 365.1] (limited by slowness grid)`, and a green cross-hair spanning those intervals is overlaid on the Back-Azimuth Detection Map. That example is a real I48 run: the interval excludes the expected 87.2° because the ~65° background windows also pass semblance 0.3, which is the mixture effect described below.
  - The intervals are **percentile** intervals from a **moving-block** bootstrap. Overlapping FK windows are correlated (a 0.1 step means 90% overlap), so blocks of window/step consecutive windows are resampled.
  - They are never narrower than the FK slowness-grid resolution (half the cell diagonal, 0.16/√2 s/km by default); `(limited by slowness grid)` marks when that floor sets the width. To go below it, use a finer slowness grid.
  - The bootstrap mixes every detected window. If several arrivals from different directions pass the semblance threshold, the interval describes their mixture, not one arrival; the PMCC families give per-arrival intervals.

  Leave it unchecked for the standard run.

### Step 5 — Run the FK analysis

Click the green **▶ Run FK Analysis** button. Processing runs in a **background thread**, so the window stays responsive, and the status line and progress bar update as it proceeds through these stages:

1. Load the (pre-merged) waveform stream.
2. Load the selected inventory (merging all XMLs if a directory entry was chosen).
3. Merge traces and clean masked samples.
4. **Remove the instrument response** (full ObsPy response removal; if that fails, SeismoFK falls back to dividing by the scalar instrument sensitivity).
5. **Run FK array processing** — band-pass filtering, then `array_processing` over the sliding windows.
6. Compute the **delay-and-sum beam** steered to the median detected direction.

If any input is missing (no MiniSEED file, no inventory selected) SeismoFK warns you before starting. If processing fails, an error dialog reports the cause.

### Step 6 — Interpret the results

When processing finishes, the **SeismoFK — Analysis Results** window opens with a six-panel figure:

**Left column (shared time axis):**

1. **Fisher** — Fisher ratio per window. A dashed line marks the threshold corresponding to your semblance threshold; an annotation reports how many windows were detected.
2. **Back-Az (°)** — back-azimuth per window. A dashed blue line marks the **expected** back-azimuth (from your source Lat/Lon), with a shaded ±20° acceptance band.
3. **App. Vel. (m/s)** — apparent velocity per window, with a shaded 300–380 m/s reference band (a typical infrasound range).
4. **Beam waveform** — the delay-and-sum beam (pressure in Pa). If Origin Time and Celerity were set, a dashed orange line marks the expected arrival.

**Right column:**

5. **Back-Azimuth Detection Map** — a polar plot (North up, clockwise) showing detections as points (angle = back-azimuth, radius = apparent velocity, colour = semblance), a semblance-weighted rose histogram, and the expected back-azimuth line.
6. **Array Geometry** — sensor positions as east/north offsets (km) from the array centroid.

In all panels, **detections** (windows at/above the semblance threshold) are drawn larger and colour-coded by semblance on a *plasma* scale; non-detections are shown small and grey. A coherent signal typically shows clustered detections at a consistent back-azimuth near the expected line.

From the results window you can:

- **💾 Save Figure (300 DPI)** — export the figure as PNG, PDF, or SVG.
- **🗄 Save to Database** — open the *Save Event* dialog (see [Section 6](#6-event-database)) to classify and archive the run.
- **✖ Close** — close the results window.

SeismoFK also **auto-saves** the results figure as `FK_<event>.jpg` in the working directory each time results are shown (see [Section 8](#8-output-files)).

### Step 7 — Advanced array methods (v1.2.1)

Below the **▶ Run FK Analysis** button, an **“Advanced (v1.2.1)”** row provides three
extra methods that open in their own windows. They reuse the **loaded MiniSEED**
and the **selected inventory** (and, where relevant, the picked **Start Time** and
**Window** length), so just load data, pick a window, and click. None of them
require running the normal FK analysis first, and none change its output.

| Button | What it does | Inputs used |
|---|---|---|
| **📡 Array Response** | Plots the theoretical **Array Response Function (ARF)** for the current sensor geometry. Use it to judge slowness resolution (main-lobe width) and spatial **aliasing** (secondary peaks) *before* trusting a slowness estimate. The white contour marks the half-power resolution limit. The dialog has an **Export figure** action. | Inventory geometry, Min/Max Freq. |
| **🎯 Slowness Map (Capon/MUSIC)** | High-resolution adaptive-beamforming **slowness map** for one window. Pick the method (**Capon / MVDR**, **MUSIC**, or **Bartlett**) in the dialog and click **↻ Recompute**; the peak (white star) gives the back-azimuth, apparent velocity and a sharpness measure (printed below the plot). Reference velocity circles and the expected-azimuth line are overlaid. Use **Export figure** to save the current map. | Picked Start Time + Window, Min/Max Freq; MUSIC source count. |
| **🔬 PMCC Detector** | Runs a **PMCC-style** detector from the picked start time over **Span** (the main **Duration**, and at least 20 analysis windows), or over the **Whole record**. Four panels share the full time axis. From top to bottom: the band-passed waveform of the first sensor (in Pa when StationXML sensitivities are available); the mean correlation of **every** window and band, with rejected windows in grey so arrivals can be judged against the noise floor and the correlation threshold; back-azimuth; and trace velocity. Accepted pixels are grouped into **families**, one per coherent arrival. Pixels join a family when they are adjacent in time and frequency and agree within **Family Δ azimuth** (10° by default) and 15% in trace velocity; a family needs at least **Family min pixels** (5). Each family's time span is shaded and labelled (for example `F1 87° 346 m/s`). Family pixels are coloured by centre frequency, and isolated pixels are drawn as hollow circles. Tune **Bands**, **Window**, **Min correlation**, **Max closure**, the velocity range and the family settings, then click **Run PMCC**. The summary lists each family's time span, back-azimuth and velocity with **95% confidence intervals**, and pixel count, plus why the other windows were rejected. A family whose direction changes significantly between its early and late pixels, or between its low and high bands, is marked **⚠ mixed**, meaning it probably contains two sources. Grouping uses single linkage, so two nearby sources can chain into one family through intermediate pixels; the flag reduces this risk but doesn't guarantee separation. The intervals come from a statistical error model checked on synthetic plane waves; they haven't been validated against ground-truth field events and exclude model errors such as wavefront curvature or wind shear. **Export families…** saves the family table as CSV, and **Export figure** saves the panels. The hop is a quarter of the window, so each arrival yields enough pixels to form a family. Needs **≥ 3 sensors**. | Start Time + Duration, Min/Max Freq. |

> These dialogs estimate **direction** from inter-sensor time delays, so they work
> on the raw loaded waveforms (band-pass filtered internally) and do **not** require
> instrument-response removal. For continuous, long-term batch processing of the
> same methods, use the command-line driver in [Section 7.4](#74-seismofk_clipy--long-term-batch-array-analysis).

PMCC aligns sensors to their common time interval and sampling grid before
cross-correlation. If the result shows *mixed directions* or a low direction
agreement, inspect the waveform and station coordinates before interpreting
the summary as a single arrival direction.

---

## 4. Plot Spectrogram

The **Spectrogram** button opens a multi-panel PSD spectrogram viewer for the loaded waveforms.

**Prerequisite:** load a MiniSEED file. StationXML is optional: when valid pressure sensitivity is available for every trace, the plot uses calibrated pressure; otherwise it clearly labels the display as raw counts. Do not compare raw-count PSD values with pressure PSD values.

**Time range:** after a preview pick, the spectrogram runs from the picked start through the main window's **Duration** setting (900 s by default). **Window length** controls FK and PMCC calculations; it does not limit the spectrogram. Without a pick, the whole loaded stream is shown. The spectrogram window displays the actual data range, which may be shorter near the end of a recording. If the picked range contains no samples, SeismoFK warns and uses the full stream.

The **Spectrograms** window shows one panel per trace, with a common colour scale. The default segment length and overlap adapt to the loaded window so short selections produce multiple time slices. Pressure plots use dB re 20 µPa²/Hz; raw-count plots use dB re 1 count²/Hz. It contains:

- **Spectrogram Parameters** control row — edit and re-apply:
  - **Bandpass min / max (Hz)** — band-pass filter corners applied before the spectrogram (seeded from the main window's Min/Max Freq).
  - **Freq. max (Hz)** — upper frequency limit shown on the panels.
  - **nperseg** — spectrogram segment length in samples.
  - **noverlap** — segment overlap in samples (must satisfy `0 ≤ noverlap < nperseg`).
  - **Smooth** — Gaussian averaging across nearby PSD bins (default 1.0). Set to 0 to inspect unsmoothed bins. This changes the displayed PSD, so use 0 when comparing exact bin values.
- **↻ Recompute** — re-runs the spectrogram with the current control-row values and redraws in place. Invalid entries (e.g. min ≥ max, non-numeric) raise a clear warning.
- A matplotlib **navigation toolbar** for pan/zoom.
- **💾 Save Figure (300 DPI)** — export the spectrogram figure as PNG, PDF, or SVG.
- **Save plotted data** — export per-channel frequency bins, UTC time bins, raw PSD, displayed PSD, and filter/smoothing settings to a NumPy `.npz` file. Open it with `numpy.load(path)` and parse `metadata_json` for the processing settings.
- **✖ Close**.

If a trace cannot be processed, the window reports the reason below the controls. When none can be plotted, the canvas shows an explicit no-data state.

---

## 5. XML Creator / Editor

The **🛠 XML Creator / Editor** button opens the **StationXML Creator** — a tool for building or editing the station inventories that FK analysis requires, without hand-writing XML. It can also be launched standalone:

```bash
python xml_creator.py
```

![XML Creator / Editor](screenshots/02_xml_creator.png)

The creator window has three sections:

### 5.1 Network

Set the **network code**, **source**, and the network **start / end dates** that apply to the inventory as a whole.

### 5.2 Stations

The **Stations** table holds one row per sensor element, with columns *Station Code*, *Latitude*, *Longitude*, *Elevation (m)*, and *Site Name*. Use the buttons to manage rows:

- **+ Add Station** — append a new station row.
- **⧉ Duplicate** — copy the selected station row.
- **− Remove** — delete the selected station row.
- **📂 Import CSV** — bulk-load stations from a CSV file. Each row should provide *code, latitude, longitude, elevation* and, optionally, a *site name*; a header row is auto-detected and skipped.

### 5.3 Channels

Select a station first, then edit its channels in the **Channels** table (*Ch. Code*, *Location*, *Sample Rate (Hz)*, *Sensitivity*, *Ref. Freq (Hz)*, *Input Units*, *Output Units*).

- **Preset dropdown + Apply Preset to New Channel** — start a channel from a preset. Presets include *Custom*, generic infrasound configurations, HLW Array configurations, IMS/CTBTO configurations, and seismic velocity channels (HHZ / BHZ). Each preset fills in a sensible channel code, sample rate, sensitivity, reference frequency, and input/output units.
- **Copy Channels → All Stations** — apply the current station's channel set to every station (handy for arrays where all elements share the same instrument).
- **+ Add Channel**, **⧉ Duplicate**, **− Remove** — manage individual channel rows.

> Infrasound channels typically use input units **PA** and output units **COUNTS**, with the sensitivity expressed in counts per pascal. The instrument **Sensitivity** value is what SeismoFK uses if full response removal is unavailable — set it correctly for your sensor.

### 5.4 Validate, preview, and save

- **🔍 Preview XML** — show the generated StationXML text in a dialog.
- **✔ Validate** — check the entered values (codes present, numeric lat/lon/elevation, numeric sample rate and sensitivity, at least one channel per station). If the basic checks pass, the XML is additionally parsed with ObsPy to confirm it loads correctly.
- **💾 Save XML** — write the StationXML file. The default save location is `~/.seismofk/XML/`, so the new inventory is auto-discovered the next time you click **⟳ Refresh** in the main window.
- **📂 Load XML** — open an existing StationXML and populate the editor for editing.

---

## 6. Event Database

The **🗄 Event Database** button opens a browser onto SeismoFK's local SQLite database of archived analyses (`~/.seismofk/fk_events.db` by default). An existing checkout `fk_events.db` remains in use for backward compatibility.

![Event database](screenshots/03_event_database.png)

### 6.1 Saving an event

After an FK run, click **🗄 Save to Database** in the results window. The **Save Event to Database** dialog shows an analysis summary and lets you:

- Pick an **Event Classification** from a fixed vocabulary: *Unknown, Explosion / Blast, Mining, Volcanic, Earthquake, Meteor / Bolide, Aircraft / Sonic Boom, Ocean / Microbaroms, Industrial, Noise / Artifact*.
- Add an optional free-text **Note**.

Click **💾 Save**. SeismoFK stores the analysis parameters, the FK results summary (median back-azimuth, median velocity, expected back-azimuth, detection counts), the source location and optional event physics, file references, and a **PNG snapshot of the results figure** — all in one database row. The database (`fk_events.db`) is created automatically on the first save.

### 6.2 Browsing the database

The **Event Database** window lists every archived event in a table (ID, saved time, event name, array, classification, median back-azimuth, median velocity, detections, frequency band, note). From here you can:

- **Filter by classification** — use the dropdown at the top to show only one event type, or *All*.
- **⟳ Refresh** — reload the table.
- **✏ Edit Classification / Note** — change the classification or note of the selected event.
- **🖼 View Figure** — display the stored analysis figure for the selected event; the viewer also offers **💾 Save as PNG** to export it.
- **🗑 Delete Selected** — permanently delete the selected event (with a confirmation prompt).
- **📤 Export to CSV** — dump the entire events table to a CSV file.

---

## 7. Command-line tools

These standalone scripts support the GUI workflow and can be run directly from a terminal.

### 7.1 `convert_ims.py` — split an FDSN StationXML into per-station files

Converts a combined IMS station file (`all_IMS_sts.xml`, FDSN StationXML) into one XML file per IMS station, grouped by sensor element, written into `XML_IM/`.

```bash
python convert_ims.py [src_file] [-o OUT_DIR]
```

| Argument | Meaning | Default |
|---|---|---|
| `src_file` (positional, optional) | Source IMS station XML file. | `all_IMS_sts.xml` next to the script |
| `-o`, `--out OUT_DIR` | Output directory for the per-station XML files. | `XML_IM/` next to the script |

Examples:

```bash
# Use the defaults (all_IMS_sts.xml → XML_IM/)
python convert_ims.py

# Explicit source file and output directory
python convert_ims.py /path/to/all_IMS_sts.xml -o XML_IM
```

The script processes all `BDF` channels, groups them by sensor element into individual `<Station>` entries, and prints a summary per array. Files written to `XML_IM/` are then auto-discovered by the main window.

### 7.2 `export_stations.py` — alternative IMS station splitter

An alternative splitter that also exports individual stations from `all_IMS_sts.xml` into `XML_IM/`, with each unique channel group becoming its own station.

```bash
python export_stations.py
```

This script takes **no command-line arguments**: it reads `all_IMS_sts.xml` from next to the script and writes per-station files to `XML_IM/`. It fails with a clear message if the source file is missing.

Use `convert_ims.py` for new conversions because it accepts source and output
paths; `export_stations.py` remains available for its fixed local-file workflow.

### 7.3 `plot_spectrogram_window.py` — standalone spectrogram exporter

Generates spectrogram and filtered-waveform PNG figures for a time window of an array dataset, headless (no GUI). Requires **SciPy** (`pip install scipy`).

```bash
python plot_spectrogram_window.py -i INPUT_FILE -x INVENTORY_FILE [options]
```

| Flag | Meaning | Default |
|---|---|---|
| `-i`, `--input-file` | **Required.** Waveform file (any format ObsPy can read, e.g. MiniSEED). | — |
| `-x`, `--inventory-file` | **Required.** StationXML inventory providing instrument response. | — |
| `-o`, `--output-dir` | Directory for the exported PNGs. | `spectrogram_exports` |
| `--window-start` | UTC start of the analysis window (ISO 8601). | `2026-04-01T23:40:00` |
| `--window-end` | UTC end of the analysis window (ISO 8601). | `2026-04-02T00:20:00` |
| `--dpi` | Resolution of the exported figures. | `300` |
| `--reference-pressure` | Reference pressure (Pa) for the dB-scaled PSD. | `2e-05` (20 µPa) |
| `--freq-max` | Maximum frequency (Hz) shown on the spectrogram. | `5.0` |
| `--nperseg` | Spectrogram segment length (samples). | `512` |
| `--noverlap` | Spectrogram segment overlap (samples). | `460` |
| `--filter-freqmin` | Lower band-pass corner (Hz). | `0.5` |
| `--filter-freqmax` | Upper band-pass corner (Hz). | `6.0` |

Run `python plot_spectrogram_window.py --help` to see the flags with their live defaults. Example:

```bash
python plot_spectrogram_window.py \
  -i data/array.mseed \
  -x XML_IM/I57US.xml \
  --window-start 2026-04-01T23:40:00 \
  --window-end   2026-04-02T00:20:00 \
  -o spectrogram_exports
```

This writes a combined spectrogram PNG and a filtered-waveform PNG into the output directory, with filenames stamped by the time window.

> The in-window **Plot Spectrogram** feature ([Section 4](#4-plot-spectrogram)) and this CLI share the same core spectrogram routine, so they produce consistent results.

### 7.4 `seismofk_cli.py` — long-term batch array analysis

**(New in v1.2.1.)** A headless driver for analysing **long, continuous** data
sets. It streams an arbitrary time range in fixed chunks, runs the chosen
array-analysis method on each chunk, and appends the detections to a Parquet
archive with a JSON metadata sidecar and a master `index.csv`. This is the
recommended path for monitoring campaigns spanning days, months or years.

```bash
python seismofk_cli.py --mseed MSEED [--inventory INVENTORY] \
    --start ISO --end ISO --method {bartlett,capon,music,pmcc} [options]
```

| Flag | Meaning | Default |
|---|---|---|
| `--mseed` | **Required.** MiniSEED file, glob pattern, or directory. | — |
| `--inventory` | StationXML file or directory of `.xml` files. | GUI inventory folder (`~/.seismofk/XML`, or `$SEISMOFK_DATA_DIR/XML`) |
| `--start`, `--end` | **Required.** Analysis time range (ISO 8601). | — |
| `--method` | `fk` runs conventional FK exactly as the GUI does: semblance, Fisher and beam, with the GUI's six-panel figure. `bartlett` / `capon` / `music` run a high-resolution FK scan, and `pmcc` runs the PMCC detector. | `capon` |
| `--event-lat`, `--event-lon` | Presumed source location (`fk`): draws the expected back-azimuth in figures. | none |
| `--semb-thresh` | `fk`: a window with semblance at or above this value is a detection (same as the GUI), and `--bootstrap` uses these windows. | `0.3` |
| `--chunk` | Processing chunk length (s) — one detection file per chunk, and the unit of memory use (only one chunk is held at a time). | `3600` |
| `--pad` | Seconds read on each side of a chunk to absorb response/filter edge transients. | one window |
| `--max-gap` | Skip a chunk when more than this fraction of its samples is missing. Gaps are zero-filled, and zeros are coherent across sensors, so they can produce false detections. | `0.2` |
| `--fmin`, `--fmax` | Analysis band (Hz). | `0.5`, `4.0` |
| `--win-length` | Sliding window length (s). | `30` |
| `--overlap` | FK window overlap fraction (0–1). | `0.5` |
| `--smax`, `--ds` | Slowness grid half-extent and step (s/km). | `4.0`; `0.16` for `fk` (as the GUI), `0.1` otherwise |
| `--n-sources` | MUSIC signal-subspace dimension. | `1` |
| `--pmcc-bands` | Number of log-spaced PMCC bands. | `6` |
| `--step` | PMCC window hop (s). A quarter window gives each arrival enough pixels to form families. | `win-length/4` |
| Uncertainty | Every PMCC pixel carries `baz_std`, `vel_std`, the arrival-time error `delay_err`, and its slowness covariance (`sx`, `sy`, `s_cxx`, `s_cxy`, `s_cyy`). Every family carries `baz_ci95`, `vel_ci95` (95% half-widths), a 95% slowness ellipse (`ell_major`, `ell_minor` in s/km, `ell_azimuth`), `chi2_p` (fit to a single plane wave), and `mixed` / `mixed_dbaz` (direction change within the family). `index.csv` counts `n_mixed_families`. | — |
| `--family-min-pixels`, `--family-baz-tol`, `--family-vel-tol` | PMCC family grouping: minimum linked pixels, maximum back-azimuth difference (°), and maximum relative velocity difference between linked pixels. Families are written to `OUT_DIR/families/<chunk run_id>` and counted in `index.csv` (`n_families`, `families_path`). An arrival that straddles a chunk boundary can be split into two families. | `5`, `10`, `0.15` |
| `--consistency`, `--corr-min` | PMCC closure-residual (s) and min correlation gates. | `0.05`, `0.5` |
| `--min-vel`, `--max-vel` | Trace-velocity acceptance gate (m/s, PMCC), and the velocity axis of saved figures. | `200`, `600` |
| `--min-coherence` | FK methods: archive only windows whose peak/mean slowness power is at least this value (always ≥ 1; ~3 isolates coherent arrivals). `0` keeps every window. `--bootstrap` uses the same windows. | `0` |
| `--bootstrap` | Add baz/vel bootstrap uncertainty to each `index.csv` row (FK methods): the mean, standard deviation and 95% percentile interval (`baz_ci95_lo/hi`, `vel_ci95_lo/hi`), from a moving-block bootstrap (`boot_block` windows per block) floored at the slowness-grid resolution (`boot_grid_limited`). | off |
| `--resp-output` / `--no-response` | Response output (`VEL`/`DISP`/`ACC`/`DEF`) for seismic sensors, or skip removal. Pressure (Pa) responses are applied as-is, so infrasound data stays in Pa for any choice. | `VEL` |
| `--check-only` | Run the same data-readiness checks as the GUI (**Check data readiness**) on the first chunk and exit: status 0 if ready, 2 if blocked. | off |
| `--skip-checks` | Process even if the readiness checks report blocking issues (collinear array, band above Nyquist, missing coordinates). Without it the run stops with status 2 before writing anything. | off |
| `--csv` | Force CSV detection tables instead of Parquet. | off |
| `--method noise` options | `--noise-window` sets the RMS window in seconds (for example 10, 60 for 1 min, or 300 for 5 min; default 60). `--noise-broadband` turns off the band-pass, which otherwise runs from `--fmin` to `--fmax`. `--noise-warn-db` / `--noise-fault-db` (3 / 6 dB) set the sensor-offset flags, and `--noise-min-coverage` (0.9) is the minimum fraction of real data per window. Levels are in dB re 20 µPa when calibrated to Pa. The per-window table goes to `OUT_DIR/noise/`, and whole-run statistics to `<run_id>_noise_stats.csv` plus `<run_id>_noise_hourly.csv` (hour-of-day medians). Each `index.csv` row has `array_L50_db`, `array_Leq_db` and `noise_flags`. | `60`, off, `3`, `6`, `0.9` |
| `--save-figures` | Save one figure per chunk under `OUT_DIR/figures/`. For `fk` this is the same six-panel figure as the GUI results window: Fisher, back-azimuth, velocity, beam, polar detection map and array geometry. For the other methods it is the band-passed waveform, then correlation (PMCC, with rejected windows in grey) or peak/mean power (Capon/MUSIC/Bartlett), back-azimuth and trace velocity over the whole chunk. Also saves a whole-run summary `<run_id>_summary.<fmt>`. Works without a display. | off |
| `--figure-format`, `--figure-dpi` | Figure file type (`png`, `pdf`, `svg`) and resolution. | `png`, `150` |
| `--array-name`, `--out-dir`, `--verbose` | Array label, output directory, progress detail. | `ARRAY`, `results`, off |
| `--version` | Print the SeismoFK version and exit. | — |

Every chunk gets an `index.csv` row, including chunks without detections, so
coverage gaps are visible. `status` is `ok`, `no_detections`, `no_data`,
`gap_skipped` or `failed`; failed rows include an `error` message. Rows also
record `n_channels`, `gap_fraction`, `units` (for example `Pa`, `m/s`, or
`counts`), `calibration` (`full_response`, `scalar_sensitivity`, or `none`),
and `figure_path` when figures are saved. `median_baz` is a circular median, so
arrivals from near north do not average to south. The JSON sidecar stores the
readiness findings.

Each chunk is analysed one window past its end, and every window is assigned to
the chunk where it starts. An arrival that straddles a chunk boundary is
therefore neither lost nor counted twice.

Run `python seismofk_cli.py --help` for the full list. Example — one month of an
IMS array with Capon FK and bootstrap uncertainty:

```bash
python seismofk_cli.py \
  --mseed "data/2025-05/*.mseed" --inventory XML_IM/ \
  --start 2025-05-01T00:00:00 --end 2025-06-01T00:00:00 \
  --method capon --fmin 0.5 --fmax 4.0 --win-length 30 --chunk 3600 \
  --array-name I31 --out-dir results/ --bootstrap
```

**Memory and very long runs (months / years).** The driver does **not** load the
whole dataset into memory. On start-up it scans the waveform files *headers only*
to build a time index, then for each chunk it reads **only** the files overlapping
that chunk (plus `--pad` seconds), processes them, writes the detections, and
frees the data before moving on. Peak memory therefore scales with `--chunk`
(one chunk at a time), **not** with the total duration — so a one-year run uses
about the same memory as a one-hour run. Because each chunk is written
immediately and the master `index.csv` is appended per chunk, a run that is
interrupted keeps all completed chunks. Repeating the command creates a new
run ID and preserves earlier files; a run with failed chunks exits nonzero.
(For data already organised as an SDS
archive, an `obspy.clients.filesystem.sds` client is an alternative back-end; the
header-index approach here works for any flat collection of files.)

**Output archive layout** (under `--out-dir`):

```
results/
├── detections/<run_id>.parquet   # one detection table per chunk (CSV if pyarrow absent)
├── families/<run_id>.parquet  # PMCC families (one row per arrival)
├── <run_id>.json                 # run metadata + provenance (parameters, channels, version)
└── index.csv                     # master catalogue, one row per chunk
```

Query the whole archive later with pandas:

```python
import pandas as pd
df = pd.read_parquet("results/detections")     # all chunks at once
df = df[df.app_vel.between(280, 450)]           # cheap columnar filter
```

> **Why Parquet + JSON?** Detection rows are tabular, typed and append-heavy —
> a columnar format (Parquet) stores them ~5–10× smaller than CSV, preserves
> dtypes, and scans fast across years of data. JSON is kept only for the small,
> nested per-run metadata. `index.csv` stays CSV because it is tiny and meant to
> be browsed by eye. See [README_v1.2.1.md](README_v1.2.1.md) for the full rationale.

---

## 8. Output files

Analysis CSV and automatic figure files are written relative to the directory
where SeismoFK is run. User data is kept in `~/.seismofk/` by default.

| File | Produced by | Contents |
|---|---|---|
| `Output_<event>.csv` | Each FK run | Per-window FK results: time, semblance, FK power, Fisher ratio, back-azimuth, apparent velocity, slowness. |
| `FK_<event>.jpg` | Each FK run (auto-saved) | The six-panel results figure, saved automatically as a 300-DPI JPEG. |
| `~/.seismofk/fk_events.db` | First time you *Save to Database* | SQLite database of archived, classified events, including a PNG snapshot of each results figure. An existing checkout database remains in use. |
| Exported figures | *Save Figure* in the results / spectrogram windows | PNG / PDF / SVG at 300 DPI, saved to a location you choose. |
| `spectrogram_data.npz` | *Save plotted data* in the spectrogram window | Per-channel time/frequency bins, raw and displayed PSD, and processing metadata. |
| Spectrogram PNGs | `plot_spectrogram_window.py` (CLI) | Combined spectrogram and filtered-waveform images, in the chosen output directory. |
| `<NETWORK>_array.xml` | XML Creator *Save XML* | A StationXML inventory, saved to `~/.seismofk/XML/` by default. |
| Per-station XMLs | `convert_ims.py` / `export_stations.py` | One StationXML file per IMS station, written to `XML_IM/`. |
| `results/` archive | `seismofk_cli.py` (v1.2.1) | Long-term detection archive: `detections/<run_id>.parquet` (or `.csv`), per-run `<run_id>.json` metadata sidecar, and master `index.csv`. |

For `Output_<event>.csv` and `FK_<event>.jpg`, characters outside letters,
digits, hyphens, and underscores in *Event Name* become underscores. The JPEG
is saved automatically when results are shown.

---

## 9. Troubleshooting and tips

**No inventory listed on startup.**
Add a StationXML file with **+ Add XML** (or build one with the XML Creator), then click **⟳ Refresh**. See [Section 2.2](#22-preparing-station-metadata-xml-inventories).

**"Please select an inventory file" when running.**
No entry is selected in the *Inventory* dropdown. Choose one, or add an XML with **+ Add XML**.

**"Stream is empty after trimming."**
The *Start Time* (plus *Duration*) falls outside the time span of the loaded waveform. Check the trace start/end times shown in the Waveform Preview panel labels, and pick a start time inside that window.

**Response removal fails / falls back to scalar sensitivity.**
SeismoFK first attempts full ObsPy instrument-response removal. If the inventory lacks full response information, FK processing falls back to scalar **instrument sensitivity**. The GUI spectrogram uses pressure when every channel has valid **PA → COUNTS** sensitivity; otherwise it displays clearly labeled raw-count PSD. Set *Sensitivity*, *Input Units*, and *Output Units* in the XML Creator for calibrated pressure plots.

**Traces at different sample rates.**
If the loaded traces have mixed sampling rates, SeismoFK resamples them all to a common rate before FK processing — no action needed, but be aware the analysis runs at the lower rate.

**Spectrogram: "noverlap must satisfy 0 ≤ noverlap < nperseg".**
In the spectrogram window's control row, *noverlap* must be smaller than *nperseg*. Reduce *noverlap* or increase *nperseg* and click **↻ Recompute**.

**No detections in the results.**
If few or no windows clear the semblance threshold, try widening the frequency band, adjusting the analysis window length, lowering the *Semb. Threshold*, or confirming the *Start Time* and *Duration* actually bracket the signal.

**Picking vs. typing the start time.**
Clicking the Waveform Preview is the quickest way to set *Start Time*, but the *Start Time* field can always be edited directly if you know the exact time.

The ±20° azimuth band and 300–380 m/s velocity band are visual reference guides,
not statistical confidence limits. Interpret detections using the array geometry,
signal quality, and the analysis settings for your data.

---

*SeismoFK — © 2024–2026 Islam Hamama, National Research Institute of Astronomy and Geophysics (NRIAG), Egypt · islam.hamama@nriag.sci.eg*
