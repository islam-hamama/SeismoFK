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

SeismoFK discovers inventories automatically from two directories next to `Infra_Analysis.py`:

- **`XML_IM/`** — intended for IMS / multi-station array inventories. If it contains any `.xml` files, a single **"★ All IMS Stations (XML_IM/)"** entry appears first in the *Inventory* dropdown; selecting it merges *every* `.xml` file in that directory into one combined inventory. Each individual file is also listed separately, tagged `[IMS]`.
- **`XML/`** — for individual / custom station files. Each `.xml` file is listed as its own entry.

Click **⟳ Refresh** at any time to rescan these directories — for example, after creating a new XML with the XML Creator.

---

## 3. Step-by-step FK analysis workflow

This is the core path through SeismoFK: load waveforms → select an inventory → set parameters → preview and pick a start time → run the analysis → review results.

### Step 1 — Load MiniSEED waveform data

In the **File Selection** panel, click **Browse** next to *MiniSEED*. Select **one or more** MiniSEED files (`.mseed`, `.msd`, `.ms`). If you select several files, their traces are read and **merged automatically** into a single stream, and the field shows how many files were merged. Gaps are filled and masked samples are replaced with zeros so processing has a continuous record.

Once data is loaded, the **Waveform Preview** populates with one normalised panel per trace, labelled with `network.station.channel`, the sampling rate, and the trace start time.

### Step 2 — Select a station inventory

In the **Inventory** dropdown, choose the StationXML entry that matches the array in your waveform data:

- For an IMS array, pick the relevant `[IMS]` entry, or **★ All IMS Stations (XML_IM/)** to merge everything in `XML_IM/`.
- For a custom array, pick its file from the `XML/` group.

If the inventory you need is not listed, click **+ Add XML** to copy a StationXML file into `XML/`, or build one with the XML Creator and then click **⟳ Refresh**.

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
| **Overlap** | Fractional window step / overlap (between 0.01 and 0.99). | 0.1 |
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
- **Celerity (m/s)** — tick the checkbox and set the infrasound propagation speed (typical range ~300–360 m/s).

When **both** Origin Time and Celerity are set, SeismoFK overlays the **expected infrasound arrival time** on the beam-waveform panel of the results figure, using the source-to-array distance and the celerity. These two controls are optional; leave them unchecked for a plain FK run.

> *For review:* please confirm the intended meaning of **Overlap** — it is passed to ObsPy's `array_processing` as `win_frac` (the window step as a fraction of window length).

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

> *For review:* SeismoFK also **auto-saves** the results figure as a JPEG named `FK_<event>.jpg` in the working directory each time results are shown (see [Section 8](#8-output-files)).

---

## 4. Plot Spectrogram

The **📊 Plot Spectrogram** button (tools row) opens a multi-panel PSD spectrogram viewer for the loaded data — useful for inspecting the frequency content of a signal alongside or instead of FK analysis.

**Prerequisites:** a MiniSEED file must be loaded and an inventory must be selected (the spectrogram needs the instrument response to convert counts to pressure).

**Time window:** the spectrogram is computed for **the picked time segment** if you have picked a start time in the preview — specifically the window from the pick to *pick + Window (s)*. If you have **not** picked a time, the **whole loaded stream** is used instead. If the picked window happens to contain no samples, SeismoFK warns and falls back to the full stream.

The **Spectrograms** window shows one spectrogram panel per trace (magma colour map, PSD in dB referenced to 20 µPa). It contains:

- **Spectrogram Parameters** control row — edit and re-apply:
  - **Bandpass min / max (Hz)** — band-pass filter corners applied before the spectrogram (seeded from the main window's Min/Max Freq).
  - **Freq. max (Hz)** — upper frequency limit shown on the panels.
  - **nperseg** — spectrogram segment length in samples.
  - **noverlap** — segment overlap in samples (must satisfy `0 ≤ noverlap < nperseg`).
- **↻ Recompute** — re-runs the spectrogram with the current control-row values and redraws in place. Invalid entries (e.g. min ≥ max, non-numeric) raise a clear warning.
- A matplotlib **navigation toolbar** for pan/zoom.
- **💾 Save Figure (300 DPI)** — export the spectrogram figure as PNG, PDF, or SVG.
- **✖ Close**.

If a particular trace cannot be processed (for example, missing or unusable response), that trace is skipped and the window reports which ones were omitted; the remaining traces are still plotted.

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
- **💾 Save XML** — write the StationXML file. The default save location is the `XML/` directory next to `Infra_Analysis.py`, so the new inventory is auto-discovered the next time you click **⟳ Refresh** in the main window.
- **📂 Load XML** — open an existing StationXML and populate the editor for editing.

---

## 6. Event Database

The **🗄 Event Database** button opens a browser onto SeismoFK's local SQLite database of archived analyses (`fk_events.db`).

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

> *For review:* `convert_ims.py` and `export_stations.py` overlap in purpose (both split `all_IMS_sts.xml` into `XML_IM/`). Please confirm which one is the recommended path for users, since only `convert_ims.py` exposes configurable input/output paths.

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

---

## 8. Output files

SeismoFK writes files relative to the directory it is run from (the project directory).

| File | Produced by | Contents |
|---|---|---|
| `Output_<event>.csv` | Each FK run | Per-window FK results: time, semblance, FK power, Fisher ratio, back-azimuth, apparent velocity, slowness. |
| `FK_<event>.jpg` | Each FK run (auto-saved) | The six-panel results figure, saved automatically as a 300-DPI JPEG. |
| `fk_events.db` | First time you *Save to Database* | SQLite database of archived, classified events, including a PNG snapshot of each results figure. |
| Exported figures | *Save Figure* in the results / spectrogram windows | PNG / PDF / SVG at 300 DPI, saved to a location you choose. |
| Spectrogram PNGs | `plot_spectrogram_window.py` (CLI) | Combined spectrogram and filtered-waveform images, in the chosen output directory. |
| `<NETWORK>_array.xml` | XML Creator *Save XML* | A StationXML inventory, saved to `XML/` by default. |
| Per-station XMLs | `convert_ims.py` / `export_stations.py` | One StationXML file per IMS station, written to `XML_IM/`. |

> *For review:* please confirm the exact event name used in `Output_<event>.csv` and `FK_<event>.jpg` — SeismoFK sanitises the *Event Name* (replacing non-alphanumeric characters with underscores) before building these filenames. Also confirm whether the auto-saved `FK_<event>.jpg` is intended behaviour or a debugging convenience.

---

## 9. Troubleshooting and tips

**"No XML or XML_IM directory found" on startup.**
Neither `XML/` nor `XML_IM/` exists next to `Infra_Analysis.py`. Create at least one, add a StationXML file (or build one with the XML Creator), and click **⟳ Refresh**. See [Section 2.2](#22-preparing-station-metadata-xml-inventories).

**"Please select an inventory file" when running.**
No entry is selected in the *Inventory* dropdown. Choose one, or add an XML with **+ Add XML**.

**"Stream is empty after trimming."**
The *Start Time* (plus *Duration*) falls outside the time span of the loaded waveform. Check the trace start/end times shown in the Waveform Preview panel labels, and pick a start time inside that window.

**Response removal fails / falls back to scalar sensitivity.**
SeismoFK first attempts full ObsPy instrument-response removal. If the inventory lacks full response information, it falls back to dividing each trace by the scalar **instrument sensitivity**. For this fallback — and for the spectrogram tools, which depend on it — every channel must carry a sensitivity value. The spectrogram routine specifically expects infrasound sensitivity in **PA → COUNTS** units; channels with other or missing units are skipped. When building inventories with the XML Creator, set the *Sensitivity*, *Input Units*, and *Output Units* fields correctly for your sensor.

**Traces at different sample rates.**
If the loaded traces have mixed sampling rates, SeismoFK resamples them all to a common rate before FK processing — no action needed, but be aware the analysis runs at the lower rate.

**Spectrogram: "noverlap must satisfy 0 ≤ noverlap < nperseg".**
In the spectrogram window's control row, *noverlap* must be smaller than *nperseg*. Reduce *noverlap* or increase *nperseg* and click **↻ Recompute**.

**No detections in the results.**
If few or no windows clear the semblance threshold, try widening the frequency band, adjusting the analysis window length, lowering the *Semb. Threshold*, or confirming the *Start Time* and *Duration* actually bracket the signal.

**Picking vs. typing the start time.**
Clicking the Waveform Preview is the quickest way to set *Start Time*, but the *Start Time* field can always be edited directly if you know the exact time.

> *For review:* the scientific interpretation in [Section 3, Step 6](#step-6--interpret-the-results) — the meaning of the reference bands, the ±20° acceptance band, and what constitutes a confident detection — should be checked against your intended methodology.

---

*SeismoFK — © 2024–2025 Islam Hamama, National Research Institute of Astronomy and Geophysics (NRIAG), Egypt · islam.hamama@nriag.sci.eg*
