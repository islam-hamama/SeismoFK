#!/usr/bin/env python3
"""
seismofk_cli.py — command-line batch driver for SeismoFK v1.2.1.

Copyright (c) 2024-2026 Islam Hamama
Contact: islam.hamama@nriag.sci.eg
Licensed under the MIT License.

Designed for *long-term* infrasound/seismic array analysis: it streams through
an arbitrarily long time range in fixed chunks, runs the chosen array-analysis
method on each chunk, and appends the detections to a Parquet archive with a
JSON metadata sidecar and a master index.csv (see results_io.py).

Examples
--------
# One month of data, Capon high-resolution FK, hourly chunks:
python seismofk_cli.py \
    --mseed "data/2025-05/*.mseed" --inventory XML/ \
    --start 2025-05-01T00:00:00 --end 2025-06-01T00:00:00 \
    --method capon --fmin 0.5 --fmax 4.0 \
    --win-length 30 --overlap 0.5 --chunk 3600 \
    --array-name I31 --out-dir results/ --bootstrap

# PMCC detector with 6 log-spaced bands:
python seismofk_cli.py --mseed data/ --inventory XML/ \
    --start 2025-05-01 --end 2025-05-02 \
    --method pmcc --fmin 0.1 --fmax 8.0 --pmcc-bands 6 \
    --win-length 30 --array-name I31 --out-dir results/

Methods: fk (conventional FK exactly as the GUI: semblance, Fisher, beam),
         bartlett | capon | music (high-resolution FK scan), pmcc, and
         noise (RMS noise levels in dB re 20 µPa per window and sensor).
"""

import argparse
import collections
import glob
import os
import re
import sys
import traceback
from datetime import datetime, timezone
from uuid import uuid4

import numpy as np
import pandas as pd
from obspy import read, read_inventory, UTCDateTime
from obspy.core.inventory import Inventory

import fk_analysis as fk
import noise_levels as nl
import pmcc as pmcc_mod
import result_plots
import results_io
from app_paths import inventory_dir
from data_readiness import assess_stream


# ─────────────────────────────────────────────────────────────────────────────
#  Data loading / preprocessing  (mirrors the GUI ProcessThread pipeline)
# ─────────────────────────────────────────────────────────────────────────────

def _list_files(spec):
    """Resolve a file / glob / directory spec into a sorted list of paths."""
    if os.path.isdir(spec):
        return sorted(glob.glob(os.path.join(spec, "*")))
    return sorted(glob.glob(spec)) or [spec]


def build_file_index(spec, verbose=False):
    """
    Scan the waveform files **headers only** and return a time index, so that
    no sample data is held in memory.  Returns

        index        : list of (starttime, endtime, path) — one per trace
        channel_ids  : sorted list of unique trace ids (NET.STA.LOC.CHA)

    `read(headonly=True)` reads only MiniSEED record headers, so indexing a
    year of files is cheap and bounded — the data itself is loaded later, one
    chunk at a time, by :func:`read_window`.
    """
    files = _list_files(spec)
    index, ids = [], set()
    for f in files:
        try:
            st = read(f, headonly=True)
        except Exception:
            continue
        for tr in st:
            index.append((tr.stats.starttime, tr.stats.endtime, f))
            ids.add(tr.id)
    if not index:
        raise FileNotFoundError(f"No readable waveform data from '{spec}'")
    if verbose:
        print(f"[INFO] Indexed {len(index)} trace segments across "
              f"{len(set(p for _, _, p in index))} file(s).")
    return index, sorted(ids)


def read_window(index, t0, t1, pad=0.0):
    """
    Read only the data for [t0 - pad, t1 + pad] by selecting the indexed files
    that overlap the window.  For MiniSEED, `read(starttime=, endtime=)` loads
    only the overlapping records, so memory stays bounded to one window.

    Returns ``(stream, gap_fraction)``: a merged Stream with gaps filled with 0,
    and the fraction of [t0, t1] samples (over all channels read) that were
    missing before filling.  Returns ``(None, 1.0)`` if no data overlap.
    Zero-filled gaps are coherent across sensors, so a large gap fraction can
    produce false detections — see ``--max-gap``.
    """
    lo, hi = t0 - pad, t1 + pad
    paths = sorted({p for (s, e, p) in index if e >= lo and s <= hi})
    if not paths:
        return None, 1.0
    st = None
    for p in paths:
        try:
            seg = read(p, starttime=lo, endtime=hi)
            st = seg if st is None else st + seg
        except Exception:
            pass
    if st is None or len(st) == 0:
        return None, 1.0
    st.merge()                                  # gaps stay masked for counting
    gap_fraction = _gap_fraction(st, t0, t1)
    for tr in st:
        if isinstance(tr.data, np.ma.MaskedArray):
            tr.data = tr.data.filled(0)
    return st, gap_fraction


def _gap_fraction(st, t0, t1):
    """Fraction of samples in [t0, t1] that are missing across *st*."""
    missing = total = 0
    for tr in st:
        expected = int(round((t1 - t0) * tr.stats.sampling_rate))
        present = int(np.ma.count(tr.slice(t0, t1).data))
        total += expected
        missing += max(expected - present, 0)
    return missing / total if total else 1.0


def load_inv(spec):
    """Read inventory from a file or merge every .xml in a directory."""
    if os.path.isdir(spec):
        inv = Inventory()
        loaded = 0
        for xf in sorted(glob.glob(os.path.join(spec, "*.xml"))):
            try:
                inv += read_inventory(xf)
                loaded += 1
            except Exception:
                pass
        if not loaded:
            raise ValueError(f"No readable StationXML files in '{spec}'")
        return inv
    return read_inventory(spec)


def preprocess(st, inv, output="VEL", pre_filt=(0.1, 0.5, 9.0, 10.0),
               water_level=60, verbose=True):
    """
    Merge, fill gaps, and remove instrument response (scalar fallback).

    Returns ``(stream, units, calibration)``: the per-trace physical units
    (e.g. 'Pa' for infrasound, even with output='VEL') and either
    'full_response' or 'scalar_sensitivity'.  Raises if any trace has no
    usable sensitivity, rather than mixing calibrated and raw data.
    """
    st.merge(fill_value=0)
    for tr in st:
        if isinstance(tr.data, np.ma.MaskedArray):
            tr.data = tr.data.filled(0)
    raw_data = [tr.data.copy() for tr in st]
    try:
        fk.require_response_stages(inv, st)
        st.remove_response(inventory=inv, output=output,
                           pre_filt=pre_filt, water_level=water_level)
        return st, [fk.calibrated_units(inv, tr, output) for tr in st], \
            "full_response"
    except Exception as e:
        if verbose:
            print(f"[WARN] response removal failed ({e}); scalar sensitivity.")
        missing = []
        for tr, data in zip(st, raw_data):
            tr.data = data
            try:
                resp = inv.get_response(tr.id, datetime=tr.stats.starttime)
                sens = resp.instrument_sensitivity.value
                if not sens:
                    raise ValueError("zero sensitivity")
                tr.data = tr.data.astype(float) / sens
            except Exception:
                missing.append(tr.id)
        if missing:
            raise ValueError("No usable instrument sensitivity for: "
                             + ", ".join(missing)) from e
    return st, [fk.calibrated_units(inv, tr, "DEF") for tr in st], \
        "scalar_sensitivity"


def summarise_units(units):
    """One label for a chunk: the shared unit, or 'mixed' if traces differ."""
    distinct = sorted(set(units))
    return distinct[0] if len(distinct) == 1 else "mixed"


def preflight(index, inv, args, t_start, t_end, pad):
    """
    Run the GUI's data-readiness checks on the first chunk that has data.

    Returns the list of findings (empty if no chunk in range has data).
    """
    cstart = t_start
    while cstart < t_end:
        cend = min(cstart + args.chunk, t_end)
        st, _ = read_window(index, cstart, cend, pad=pad)
        if st is not None:
            return assess_stream(st, inv, fmin=args.fmin, fmax=args.fmax,
                                 start=cstart, duration=cend - cstart)
        cstart = cend
    return []


# ─────────────────────────────────────────────────────────────────────────────
#  Per-chunk analysis
# ─────────────────────────────────────────────────────────────────────────────

def analyse_chunk(st, args, cstart, cend, diagnostics=None):
    """Run the selected method on one trimmed chunk → detection DataFrame.

    *diagnostics* (PMCC only) receives rejection counts and every candidate
    window, used for the signal-versus-noise figure.
    """
    if args.method == "pmcc":
        bands = pmcc_mod.log_bands(args.fmin, args.fmax, args.pmcc_bands)
        return pmcc_mod.pmcc(
            st, bands, window_sec=args.win_length, step_sec=args.step,
            consistency_max=args.consistency, corr_min=args.corr_min,
            min_app_vel=args.min_vel, max_app_vel=args.max_vel,
            stime=cstart, etime=cend, verbose=args.verbose,
            diagnostics=diagnostics)
    # bartlett | capon | music
    return fk.fk_scan(
        st, args.fmin, args.fmax, win_length=args.win_length,
        overlap=args.overlap, method=args.method,
        smax=args.smax, ds=args.ds, stime=cstart, etime=cend,
        n_sources=args.n_sources, verbose=args.verbose)


def analyse_fk_chunk(st, inv, args, cstart, aend):
    """
    Conventional FK on one chunk with the GUI pipeline
    (:func:`fk_analysis.fk_array`).  Returns ``(windows, result)``: a
    DataFrame with one row per window (time = window centre) and the full
    result dict used by the GUI-style figure.
    """
    result = fk.fk_array(
        st, inv, args.fmin, args.fmax, args.win_length,
        1.0 - args.overlap,               # fk_array takes the step fraction
        cstart, float(aend - cstart), args.array_name,
        args.event_lat, args.event_lon,
        max_slowness_skm=args.smax, slowness_step_skm=args.ds,
        abs_semb_thresh=args.semb_thresh, save_csv=False,
        verbose=args.verbose)
    centre = pd.to_timedelta(args.win_length / 2.0, unit="s")
    windows = pd.DataFrame({
        "time": pd.DatetimeIndex(result["time"]) + centre,
        "baz": result["bazi"], "app_vel": result["app_vel"],
        "slowness": result["slowness"], "semblance": result["semblance"],
        "fisher": result["fisher"]})
    return windows, result


def noise_chunk(st, inv, args, cstart, cend, aend, last):
    """
    Noise levels for one chunk: StationXML-sensitivity calibration (exact in
    the sensor passband, and no response pre-filter to cut the band), then
    RMS per ``--noise-window`` on the absolute time grid.  Windows are kept
    when they *start* in [cstart, cend) — the look-ahead read lets the last
    one finish — so chunking never drops or duplicates a window.

    Returns (levels, units, calibration).
    """
    if args.no_response:
        units, calibration = "counts", "none"
    else:
        units = fk.scalar_calibrate(st, inv)
        calibration = "scalar_sensitivity" if units != "counts" else "none"
    fmin, fmax = (None, None) if args.noise_broadband else (args.fmin, args.fmax)
    levels = nl.noise_levels(st, args.noise_window, fmin, fmax, t0=cstart,
                             t1=aend, units=units,
                             min_coverage=args.noise_min_coverage)
    starts = pd.to_datetime(levels["time"])
    keep = starts < pd.Timestamp(cend.datetime) if not last else starts == starts
    return levels[keep.to_numpy()].reset_index(drop=True), units, calibration


def rows_in_chunk(df, cstart, cend, win_length, last):
    """
    Keep the windows that belong to [cstart, cend).

    Each chunk is analysed one window past its end so a window straddling
    the boundary is still computed; its centre decides which chunk keeps it
    (centre < cend + win/2, i.e. it starts before cend).  The last chunk
    keeps everything up to the requested end.
    """
    if df is None or df.empty:
        return df
    times = pd.to_datetime(df["time"])
    upper = cend if last else cend + win_length / 2.0
    keep = (times >= pd.Timestamp(cstart.datetime)) & (
        (times <= pd.Timestamp(upper.datetime)) if last
        else (times < pd.Timestamp(upper.datetime)))
    return df[keep.to_numpy()].reset_index(drop=True)


def save_chunk_figure(path, df, st, args, t0, t1, units, windows=None,
                      title=None, families=None):
    """Write the per-chunk result figure to *path*.  For ``fk`` *df* is the
    fk_array result dict; otherwise the detection DataFrame."""
    from matplotlib.figure import Figure

    if args.method == "noise":
        fig = Figure(figsize=(12, 8), constrained_layout=True)
        stats = nl.noise_statistics(df, warn_db=args.noise_warn_db,
                                    fault_db=args.noise_fault_db, units=units)
        result_plots.noise_figure(
            fig, df, stats, units=units, band_label=noise_band_label(args),
            t0=t0, t1=t1, warn_db=args.noise_warn_db,
            fault_db=args.noise_fault_db, title=title)
        fig.savefig(path, dpi=args.figure_dpi)
        return
    if args.method == "fk":
        # Same six-panel figure as the GUI results window.
        fig = Figure(figsize=(15, 11))
        result_plots.fk_results_figure(
            fig, df, {"bootstrap": args.bootstrap,
                      "overlap": 1.0 - args.overlap,   # step fraction
                      "slowness_step": args.ds,
                      "event_lat": args.event_lat,
                      "event_lon": args.event_lon},
            event_name=title, fmin=args.fmin, fmax=args.fmax,
            semb_thresh=args.semb_thresh)
        fig.savefig(path, dpi=args.figure_dpi, bbox_inches="tight")
        return
    fig = Figure(figsize=(11, 8.5), constrained_layout=True)
    if args.method == "pmcc":
        result_plots.pmcc_figure(
            fig, df, windows, st, fmin=args.fmin, fmax=args.fmax, t0=t0, t1=t1,
            corr_min=args.corr_min, min_vel=args.min_vel,
            max_vel=args.max_vel, units=units, title=title,
            families=families, window_sec=args.win_length)
    else:
        result_plots.fk_scan_figure(
            fig, df, st, fmin=args.fmin, fmax=args.fmax, t0=t0, t1=t1,
            method=args.method, units=units,
            coherence_min=args.min_coherence or None, min_vel=args.min_vel,
            max_vel=args.max_vel, title=title)
    fig.savefig(path, dpi=args.figure_dpi)


def noise_band_label(args):
    return ("broadband" if args.noise_broadband
            else f"{args.fmin:g}–{args.fmax:g} Hz")


def save_summary_figure(path, parts, args, t_start, t_end, families=None):
    """Write the whole-run overview figure from the collected detections."""
    from matplotlib.figure import Figure

    df = pd.concat(parts, ignore_index=True) if parts else None
    fams = pd.concat(families, ignore_index=True) if families else None
    fig = Figure(figsize=(11, 6), constrained_layout=True)
    result_plots.run_summary_figure(
        fig, df, families=fams, method=args.method, t0=t_start, t1=t_end,
        fmin=args.fmin,
        fmax=args.fmax,
        min_vel=args.min_vel, max_vel=args.max_vel,
        title=(f"{args.array_name}  ·  {args.method.upper()}  ·  "
               f"{args.fmin:g}–{args.fmax:g} Hz  ·  "
               f"{t_start.strftime('%Y-%m-%d %H:%M')} – "
               f"{t_end.strftime('%Y-%m-%d %H:%M')} UTC"))
    fig.savefig(path, dpi=args.figure_dpi)


def make_run_tag(array_name, t_start, t_end, method):
    """Return a unique file-safe prefix for one CLI invocation."""
    invocation = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    return (f"{array_name}_{t_start.strftime('%Y%m%dT%H%M%S')}_"
            f"{t_end.strftime('%Y%m%dT%H%M%S')}_{method}_"
            f"{invocation}_{uuid4().hex[:8]}")


# ─────────────────────────────────────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────────────────────────────────────

_EXAMPLES = """\
examples:
  # One month of data, Capon high-resolution FK, hourly chunks, with bootstrap CI
  python seismofk_cli.py \\
      --mseed "data/2025-05/*.mseed" --inventory XML_IM/ \\
      --start 2025-05-01T00:00:00 --end 2025-06-01T00:00:00 \\
      --method capon --fmin 0.5 --fmax 4.0 \\
      --win-length 30 --overlap 0.5 --chunk 3600 \\
      --array-name I31 --out-dir results/ --bootstrap

  # PMCC detector (IMS-style) with 6 log-spaced bands
  python seismofk_cli.py --mseed data/ --inventory XML_IM/ \\
      --start 2025-05-01 --end 2025-05-02 \\
      --method pmcc --fmin 0.1 --fmax 8.0 --pmcc-bands 6 \\
      --win-length 30 --array-name I31 --out-dir results/

  # MUSIC, two sources, 10-minute chunks, force CSV output (no pyarrow)
  python seismofk_cli.py --mseed array.mseed --inventory inv.xml \\
      --start 2025-05-01T03:00:00 --end 2025-05-01T04:00:00 \\
      --method music --n-sources 2 --chunk 600 \\
      --array-name HLW --out-dir results/ --csv

  # One day, conventional Bartlett FK, skip response removal (data already in Pa)
  python seismofk_cli.py --mseed data/ --inventory XML_IM/ \\
      --start 2025-05-01 --end 2025-05-02 --method bartlett \\
      --fmin 0.5 --fmax 4.0 --no-response --array-name I31 --out-dir results/

Query the archive afterwards:
  python -c "import pandas as pd; print(pd.read_parquet('results/detections').describe())"
"""


# Show argument defaults *and* keep the examples block formatted verbatim.
class _HelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                     argparse.RawDescriptionHelpFormatter):
    pass


def build_parser():
    p = argparse.ArgumentParser(
        prog="seismofk_cli",
        description="SeismoFK v1.2.1 long-term batch array analysis.",
        epilog=_EXAMPLES,
        formatter_class=_HelpFormatter)
    p.add_argument("--version", action="version",
                   version=f"%(prog)s {fk.__version__}")

    io = p.add_argument_group("data I/O")
    io.add_argument("--mseed", required=True,
                    help="MiniSEED file, glob pattern, or directory.")
    io.add_argument("--inventory", default=None,
                    help="StationXML file or directory of .xml files "
                         "(default: the GUI's inventory folder, "
                         "$SEISMOFK_DATA_DIR/XML or ~/.seismofk/XML).")
    io.add_argument("--out-dir", default="results",
                    help="Output directory for the detection archive.")
    io.add_argument("--array-name", default="ARRAY",
                    help="Array label used in run IDs / index (letters, digits, ._-).")

    tm = p.add_argument_group("time range")
    tm.add_argument("--start", required=True, help="Analysis start (ISO 8601).")
    tm.add_argument("--end", required=True, help="Analysis end (ISO 8601).")
    tm.add_argument("--chunk", type=float, default=3600.0,
                    help="Processing chunk length in seconds (also bounds peak "
                         "memory — only one chunk is held at a time).")
    tm.add_argument("--pad", type=float, default=None,
                    help="Seconds of data read on each side of a chunk to absorb "
                         "response/filter edge transients (default: one window).")
    tm.add_argument("--max-gap", type=float, default=0.2,
                    help="Skip chunks whose fraction of missing (zero-filled) "
                         "samples exceeds this value (0-1).")

    an = p.add_argument_group("analysis")
    an.add_argument("--method", default="capon",
                    choices=["fk", "bartlett", "capon", "music", "pmcc",
                             "noise"],
                    help="fk = conventional FK exactly as in the GUI "
                         "(semblance, Fisher, beam; GUI-style figures); "
                         "bartlett/capon/music = high-resolution FK scan; "
                         "pmcc = PMCC detector; noise = RMS noise levels "
                         "(dB re 20 µPa) per window and sensor.")
    an.add_argument("--fmin", type=float, default=0.5, help="Band min (Hz).")
    an.add_argument("--fmax", type=float, default=4.0, help="Band max (Hz).")
    an.add_argument("--win-length", type=float, default=30.0,
                    help="Sliding window length (s).")
    an.add_argument("--overlap", type=float, default=0.5,
                    help="FK window overlap fraction (0-1).")
    an.add_argument("--step", type=float, default=None,
                    help="PMCC window hop (s); default win-length/4, which "
                         "gives arrivals enough pixels to form families.")
    an.add_argument("--smax", type=float, default=4.0,
                    help="Slowness grid half-extent (s/km).")
    an.add_argument("--ds", type=float, default=None,
                    help="Slowness grid step (s/km); default 0.16 for fk "
                         "(as the GUI) and 0.1 otherwise.")
    an.add_argument("--event-lat", type=float, default=None,
                    help="Presumed source latitude (fk): draws the expected "
                         "back-azimuth in figures.")
    an.add_argument("--event-lon", type=float, default=None,
                    help="Presumed source longitude (fk).")
    an.add_argument("--n-sources", type=int, default=1,
                    help="MUSIC signal-subspace dimension.")
    an.add_argument("--min-vel", type=float, default=200.0,
                    help="Min accepted trace velocity (m/s, PMCC gate); "
                         "also the velocity axis of saved figures.")
    an.add_argument("--max-vel", type=float, default=600.0,
                    help="Max accepted trace velocity (m/s, PMCC gate); "
                         "also the velocity axis of saved figures.")
    an.add_argument("--pmcc-bands", type=int, default=6,
                    help="Number of log-spaced PMCC bands.")
    an.add_argument("--consistency", type=float, default=0.05,
                    help="PMCC max closure residual (s).")
    an.add_argument("--corr-min", type=float, default=0.5,
                    help="PMCC min mean correlation.")
    an.add_argument("--family-min-pixels", type=int, default=5,
                    help="PMCC: minimum linked pixels for a family (one "
                         "detection).")
    an.add_argument("--family-baz-tol", type=float, default=10.0,
                    help="PMCC: max back-azimuth difference (deg) between "
                         "linked pixels.")
    an.add_argument("--family-vel-tol", type=float, default=0.15,
                    help="PMCC: max relative trace-velocity difference "
                         "between linked pixels.")

    pp = p.add_argument_group("preprocessing")
    pp.add_argument("--resp-output", default="VEL",
                    choices=["VEL", "DISP", "ACC", "DEF"],
                    help="ObsPy response output for seismic sensors. Pressure "
                         "(Pa) responses are applied as-is, so infrasound "
                         "data stays in Pa for any choice.")
    pp.add_argument("--no-response", action="store_true",
                    help="Skip instrument-response removal.")

    ck = p.add_argument_group("data readiness")
    ck.add_argument("--check-only", action="store_true",
                    help="Run the readiness checks on the first chunk and exit "
                         "(status 2 if any blocking issue is found).")
    ck.add_argument("--skip-checks", action="store_true",
                    help="Process even if the readiness checks report "
                         "blocking issues.")

    nz = p.add_argument_group("noise levels (--method noise)")
    nz.add_argument("--noise-window", type=float, default=60.0,
                    help="RMS window (s): e.g. 10, 60 (1 min), 300 (5 min).")
    nz.add_argument("--noise-broadband", action="store_true",
                    help="No band-pass (default: --fmin to --fmax).")
    nz.add_argument("--noise-warn-db", type=float, default=3.0,
                    help="Warn when a sensor's median offset from the other "
                         "sensors reaches this (dB).")
    nz.add_argument("--noise-fault-db", type=float, default=6.0,
                    help="Fault at this offset (dB; 6 dB = factor 2).")
    nz.add_argument("--noise-min-coverage", type=float, default=0.9,
                    help="Minimum fraction of real (non-gap) samples for a "
                         "window's level.")

    fg = p.add_argument_group("figures")
    fg.add_argument("--save-figures", action="store_true",
                    help="Save a figure per chunk (band-passed waveform, "
                         "correlation or peak/mean power, back-azimuth, trace "
                         "velocity) under OUT_DIR/figures/, plus a whole-run "
                         "summary figure.")
    fg.add_argument("--figure-format", default="png",
                    choices=["png", "pdf", "svg"])
    fg.add_argument("--figure-dpi", type=int, default=150)

    ex = p.add_argument_group("extras")
    ex.add_argument("--min-coherence", type=float, default=0.0,
                    help="FK methods: archive only windows whose peak/mean "
                         "slowness power is at least this value (>= 1; 0 "
                         "keeps every window). Also selects the windows "
                         "used by --bootstrap.")
    ex.add_argument("--bootstrap", action="store_true",
                    help="Bootstrap baz/vel uncertainty per chunk (FK methods).")
    ex.add_argument("--n-boot", type=int, default=1000)
    ex.add_argument("--semb-thresh", type=float, default=0.3,
                    help="fk: semblance (0-1) at or above which a window is a "
                         "detection, as in the GUI; also selects the windows "
                         "used by --bootstrap.")
    ex.add_argument("--csv", action="store_true",
                    help="Force CSV output instead of Parquet.")
    ex.add_argument("--verbose", action="store_true")
    return p


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.array_name):
        parser.error("--array-name must use letters, digits, dots, underscores or hyphens")
    if args.chunk <= 0 or args.win_length <= 0:
        parser.error("--chunk and --win-length must be positive")
    if args.pad is not None and args.pad < 0:
        parser.error("--pad must be non-negative")
    if not (0 < args.fmin < args.fmax):
        parser.error("frequency band must satisfy 0 < --fmin < --fmax")
    if not (0 <= args.overlap < 1):
        parser.error("--overlap must be in [0, 1)")
    if args.step is not None and args.step <= 0:
        parser.error("--step must be positive")
    if args.step is None:
        args.step = args.win_length / 4.0
    if (args.family_min_pixels < 1 or not (0 < args.family_baz_tol <= 180)
            or args.family_vel_tol <= 0):
        parser.error("family options: --family-min-pixels >= 1, "
                     "0 < --family-baz-tol <= 180, --family-vel-tol > 0")
    if args.smax <= 0 or (args.ds is not None and args.ds <= 0) \
            or args.n_sources < 1:
        parser.error("--smax, --ds, and --n-sources must be positive")
    if args.pmcc_bands < 1 or not (0 <= args.corr_min <= 1):
        parser.error("--pmcc-bands must be positive and --corr-min in [0, 1]")
    if args.consistency < 0 or not (0 < args.min_vel < args.max_vel):
        parser.error("PMCC requires non-negative consistency and 0 < min-vel < max-vel")
    if args.n_boot < 1 or args.min_coherence < 0:
        parser.error("--n-boot must be positive and --min-coherence non-negative")
    if not (0 <= args.semb_thresh < 1):
        parser.error("--semb-thresh must be in [0, 1)")
    if (args.event_lat is None) != (args.event_lon is None):
        parser.error("give both --event-lat and --event-lon, or neither")
    if args.ds is None:
        args.ds = 0.16 if args.method == "fk" else 0.1
    if args.figure_dpi < 20:
        parser.error("--figure-dpi must be at least 20")
    if not (0 <= args.max_gap <= 1):
        parser.error("--max-gap must be in [0, 1]")
    if args.method == "noise":
        if not (0 < args.noise_window <= args.chunk):
            parser.error("--noise-window must be positive and <= --chunk")
        if not (0 < args.noise_warn_db < args.noise_fault_db):
            parser.error("need 0 < --noise-warn-db < --noise-fault-db")
        if not (0 < args.noise_min_coverage <= 1):
            parser.error("--noise-min-coverage must be in (0, 1]")
    if args.inventory is None:
        args.inventory = str(inventory_dir())
        if not os.path.isdir(args.inventory):
            parser.error(f"--inventory not given and {args.inventory} does not exist")
    t_start = UTCDateTime(args.start)
    t_end = UTCDateTime(args.end)
    if t_end <= t_start:
        sys.exit("[ERROR] --end must be after --start")

    print(f"[INFO] Indexing waveform files from {args.mseed} ...")
    index, channel_ids = build_file_index(args.mseed, verbose=args.verbose)
    print(f"[INFO] Loading inventory from {args.inventory} ...")
    inv = load_inv(args.inventory)

    data_start = min(s for s, _, _ in index)
    data_end = max(e for _, e, _ in index)
    n_sta = len(channel_ids)
    if args.method == "music" and args.n_sources >= n_sta:
        parser.error("--n-sources must be smaller than the channel count")
    # Edge padding (s) so per-chunk response removal / filtering transients fall
    # outside the analysed window.  Default: one analysis window.
    pad = args.pad if args.pad is not None else float(args.win_length)
    print(f"[INFO] {n_sta} channels; data spans {data_start} → {data_end}")

    # Same readiness checks as the GUI, on the first chunk that has data.
    findings = preflight(index, inv, args, t_start, t_end, pad)
    geometry_only = ("Collinear array geometry", "Too few sensors")
    blocking = [f for f in findings if f.severity == "error"
                and not (args.method == "noise"
                         and (f.title in geometry_only
                              or f.title.startswith("Missing coordinates")))]
    for f in findings:
        if f.severity != "info" or args.check_only:
            print(f"[{f.severity.upper()}] {f.title}: {f.detail}")
    if args.check_only:
        sys.exit(2 if blocking else 0)
    if blocking and not args.skip_checks:
        print(f"[ERROR] {len(blocking)} blocking readiness issue(s); fix them "
              "or pass --skip-checks.", file=sys.stderr)
        sys.exit(2)

    run_tag = make_run_tag(args.array_name, t_start, t_end, args.method)

    # Run-level JSON metadata sidecar (written once for the whole invocation).
    metadata = {
        "software": "SeismoFK", "version": fk.__version__,
        "array": args.array_name, "method": args.method,
        "start": str(t_start), "end": str(t_end), "chunk_sec": args.chunk,
        "freq_min": args.fmin, "freq_max": args.fmax,
        "win_length": args.win_length, "overlap": args.overlap,
        "smax": args.smax, "ds": args.ds, "pad_sec": pad,
        "n_channels": n_sta,
        "channels": channel_ids,
        "data_start": str(data_start), "data_end": str(data_end),
        "response_removed": (not args.no_response),
        "resp_output": args.resp_output,
        "max_gap": args.max_gap,
        "min_coherence": args.min_coherence,
        "semb_thresh": args.semb_thresh,
        "noise_window": args.noise_window,
        "noise_band": noise_band_label(args) if args.method == "noise" else None,
        "pmcc_step": args.step,
        "family_min_pixels": args.family_min_pixels,
        "family_baz_tol": args.family_baz_tol,
        "family_vel_tol": args.family_vel_tol,
        "readiness": [{"severity": f.severity, "title": f.title}
                      for f in findings],
        "cli_args": vars(args),
    }
    results_io.write_metadata(args.out_dir, run_tag, metadata)

    total_det = 0
    noise_parts = []            # compact levels for the whole-run statistics
    noise_unit_set = set()
    status_counts = collections.Counter()
    summary_parts = []          # compact detections for the run figure
    summary_families = []
    total_families = 0
    figure_dir = os.path.join(args.out_dir, "figures")
    if args.save_figures:
        os.makedirs(figure_dir, exist_ok=True)
        result_plots.apply_plot_style()      # same typography as the GUI
    n_chunks = int(np.ceil((t_end - t_start) / args.chunk))
    print(f"[INFO] Processing {n_chunks} chunk(s) of {args.chunk:.0f}s "
          f"with method='{args.method}' → {args.out_dir}/")

    cstart = t_start
    chunk_i = 0
    while cstart < t_end:
        cend = min(cstart + args.chunk, t_end)
        last = cend >= t_end
        # Analyse one window past the chunk end (see rows_in_chunk).
        lookahead = (args.noise_window if args.method == "noise"
                     else args.win_length)
        aend = cend if last else min(cend + lookahead, t_end)
        chunk_i += 1
        chunk_id = f"{run_tag}_c{chunk_i:05d}"
        # Every chunk gets an index row, so gaps in a long catalogue can be
        # told apart: ok | no_detections | no_data | gap_skipped | failed.
        row = {"run_id": chunk_id, "array": args.array_name,
               "method": args.method, "chunk_start": str(cstart),
               "chunk_end": str(cend), "freq_min": args.fmin,
               "freq_max": args.fmax, "n_detections": 0}
        df = None
        try:
            # Read only this chunk (+ pad) into memory, process, then free it.
            st, gap_frac = read_window(index, cstart, aend, pad=pad)
            if st is None:
                row["status"] = "no_data"
                if args.verbose:
                    print(f"[INFO] chunk {chunk_i}/{n_chunks}: no data")
            else:
                row["gap_fraction"] = round(gap_frac, 4)
            if st is not None and gap_frac > args.max_gap:
                row["status"] = "gap_skipped"
                print(f"[WARN] chunk {chunk_i}/{n_chunks} ({cstart}): "
                      f"{gap_frac:.0%} of samples missing > --max-gap "
                      f"{args.max_gap:.0%}; skipped")
            elif st is not None and args.method == "noise":
                chunk_inv = fk.resolve_epochs(inv, st)
                df, units, calibration = noise_chunk(
                    st, chunk_inv, args, cstart, cend, aend, last)
                row.update(n_channels=len(st), units=units,
                           calibration=calibration)
                row["status"] = ("ok" if df["level_db"].notna().any()
                                 else "no_detections")
                if args.save_figures and row["status"] == "ok":
                    fig_path = os.path.join(
                        figure_dir, f"{chunk_id}.{args.figure_format}")
                    save_chunk_figure(
                        fig_path, df, st, args, cstart, cend, units,
                        title=(f"{args.array_name}  ·  noise levels  ·  "
                               f"{noise_band_label(args)}  ·  "
                               f"{cstart.strftime('%Y-%m-%d %H:%M')} – "
                               f"{cend.strftime('%H:%M')} UTC"))
                    row["figure_path"] = os.path.relpath(fig_path, args.out_dir)
            elif st is not None:
                # One channel epoch per trace (latest calibration wins);
                # resolved per chunk as a long run can cross a calibration.
                chunk_inv = fk.resolve_epochs(inv, st)
                if args.no_response:
                    units, calibration = ["counts"] * len(st), "none"
                else:
                    st, units, calibration = preprocess(
                        st, chunk_inv, output=args.resp_output,
                        verbose=args.verbose)
                row.update(n_channels=len(st), units=summarise_units(units),
                           calibration=calibration)
                fk.attach_coordinates(st, chunk_inv, verbose=args.verbose)
                diagnostics = {} if args.save_figures else None
                fk_result = None
                if args.method == "fk":
                    windows, fk_result = analyse_fk_chunk(
                        st, chunk_inv, args, cstart, aend)
                    fk_result["beam_units"] = row["units"]
                else:
                    windows = analyse_chunk(st, args, cstart, aend,
                                            diagnostics)
                all_windows = rows_in_chunk(windows, cstart, cend,
                                            args.win_length, last)
                families = None
                if args.method == "pmcc" and all_windows is not None:
                    # Families are formed within a chunk; an arrival that
                    # straddles a chunk boundary can be split in two.
                    all_windows, families = pmcc_mod.pmcc_families(
                        all_windows, time_tol=args.win_length,
                        baz_tol=args.family_baz_tol,
                        vel_tol=args.family_vel_tol,
                        min_pixels=args.family_min_pixels,
                        pixel_overlap=args.win_length / args.step)
                    row["n_families"] = len(families)
                    row["n_mixed_families"] = int(families["mixed"].sum())
                df = all_windows
                if df is not None and not df.empty:
                    if args.method == "fk":
                        df = df[df["semblance"] >= args.semb_thresh]
                    elif args.method != "pmcc" and args.min_coherence:
                        df = df[df["coherence"] >= args.min_coherence]
                    df = df.reset_index(drop=True)
                row["status"] = ("ok" if df is not None and not df.empty
                                 else "no_detections")
                if args.save_figures:
                    fig_path = os.path.join(
                        figure_dir, f"{chunk_id}.{args.figure_format}")
                    save_chunk_figure(
                        fig_path,
                        fk_result if fk_result is not None else all_windows,
                        st, args, cstart, cend,
                        row["units"],
                        windows=(diagnostics or {}).get("windows"),
                        title=(f"{args.array_name}  ·  {args.method.upper()}"
                               f"  ·  chunk {chunk_i}/{n_chunks}  ·  "
                               f"{cstart.strftime('%Y-%m-%d %H:%M:%S')} – "
                               f"{cend.strftime('%H:%M:%S')} UTC"),
                        families=families)
                    row["figure_path"] = os.path.relpath(fig_path, args.out_dir)
            del st                       # release chunk memory before the next
        except Exception as e:
            print(f"[WARN] chunk {chunk_i}/{n_chunks} ({cstart}): {e}")
            row.update(status="failed", error=str(e)[:300])
            df = None
            if args.verbose:
                traceback.print_exc()

        if row["status"] == "no_detections" and args.verbose:
            print(f"[INFO] chunk {chunk_i}/{n_chunks}: no detections")
        if row["status"] == "ok" and args.method == "noise":
            lv_path, fmt = results_io.write_detections(
                args.out_dir, chunk_id, df, prefer_parquet=not args.csv,
                subdir="noise")
            stats = nl.noise_statistics(df, warn_db=args.noise_warn_db,
                                        fault_db=args.noise_fault_db,
                                        units=row["units"])
            flagged = stats[stats["flag"] != "ok"]
            row.update(
                n_windows=int(df["level_db"].notna().sum()),
                array_L50_db=round(float(stats["L50"].median()), 2),
                array_Leq_db=round(float(stats["Leq"].median()), 2),
                noise_flags=";".join(f"{r.station}:{r.flag}"
                                     for r in flagged.itertuples()),
                detections_path=os.path.relpath(lv_path, args.out_dir),
                detections_fmt=fmt)
            noise_parts.append(df[["time", "trace_id", "station", "rms",
                                   "level_db"]])
            noise_unit_set.add(row["units"])
            print(f"[INFO] chunk {chunk_i}/{n_chunks} ({cstart} → {cend}): "
                  f"array L50 {row['array_L50_db']:.1f} dB, "
                  f"{row['n_windows']} levels"
                  + (f", flagged {row['noise_flags']}" if len(flagged) else ""))
        elif row["status"] == "ok":
            det_path, fmt = results_io.write_detections(
                args.out_dir, chunk_id, df, prefer_parquet=not args.csv)
            row.update(
                n_detections=len(df),
                median_baz=round(pmcc_mod.circular_median_deg(df["baz"]), 2),
                median_app_vel=round(float(np.median(df["app_vel"])), 1),
                detections_path=os.path.relpath(det_path, args.out_dir),
                detections_fmt=fmt)
            if families is not None and not families.empty:
                # Family numbers restart per chunk; run_id makes them unique.
                families.insert(0, "run_id", chunk_id)
                fam_path, _ = results_io.write_detections(
                    args.out_dir, chunk_id, families,
                    prefer_parquet=not args.csv, subdir="families")
                row["families_path"] = os.path.relpath(fam_path, args.out_dir)
                total_families += len(families)
                if args.save_figures:
                    summary_families.append(families)

            # Optional bootstrap uncertainty for FK methods (on the windows
            # kept by --min-coherence).
            if args.bootstrap and args.method != "pmcc" and len(df) >= 2:
                try:
                    # Overlapping windows are correlated: resample blocks of
                    # window/hop consecutive windows (moving-block bootstrap).
                    boot = fk.bootstrap_beam(
                        df["baz"].to_numpy(), df["app_vel"].to_numpy(),
                        df["slowness"].to_numpy(), n_boot=args.n_boot,
                        block=max(1, round(1.0 / (1.0 - args.overlap))),
                        slowness_step=args.ds)
                    row.update({
                        "baz_mean": round(boot["baz_mean"], 3),
                        "baz_std": round(boot["baz_std"], 3),
                        "baz_ci95_lo": round(boot["baz_ci"][0], 3),
                        "baz_ci95_hi": round(boot["baz_ci"][1], 3),
                        "vel_mean": round(boot["vel_mean"], 2),
                        "vel_std": round(boot["vel_std"], 2),
                        "vel_ci95_lo": round(boot["vel_ci"][0], 2),
                        "vel_ci95_hi": round(boot["vel_ci"][1], 2),
                        "boot_block": boot["block"],
                        "boot_grid_limited": boot["resolution_limited"],
                    })
                except Exception as e:
                    if args.verbose:
                        print(f"[WARN] bootstrap failed: {e}")
            if args.save_figures:
                colour = {"pmcc": "f_center", "fk": "semblance"}.get(
                    args.method, "coherence")
                summary_parts.append(df[["time", "baz", "app_vel", colour]])
            total_det += len(df)
            print(f"[INFO] chunk {chunk_i}/{n_chunks} ({cstart} → {cend}): "
                  f"{len(df)} detections → {os.path.basename(det_path)}")

        results_io.append_index(args.out_dir, row)
        status_counts[row["status"]] += 1
        cstart = cend

    noise_stats_path = None
    if args.method == "noise" and noise_parts:
        run_levels = pd.concat(noise_parts, ignore_index=True)
        # dB re 20 µPa only if every chunk was calibrated to Pa.
        noise_units = "Pa" if noise_unit_set == {"Pa"} else "counts"
        run_stats = nl.noise_statistics(
            run_levels, warn_db=args.noise_warn_db,
            fault_db=args.noise_fault_db, units=noise_units)
        noise_stats_path = os.path.join(args.out_dir,
                                        f"{run_tag}_noise_stats.csv")
        run_stats.to_csv(noise_stats_path, index=False)
        nl.hourly_levels(run_levels).to_csv(
            os.path.join(args.out_dir, f"{run_tag}_noise_hourly.csv"),
            index_label="hour_utc")

    summary_path = None
    if args.save_figures and args.method == "noise" and noise_parts:
        from matplotlib.figure import Figure

        summary_path = os.path.join(
            args.out_dir, f"{run_tag}_summary.{args.figure_format}")
        fig = Figure(figsize=(12, 8), constrained_layout=True)
        result_plots.noise_figure(
            fig, run_levels, run_stats, units=noise_units,
            band_label=noise_band_label(args), t0=t_start, t1=t_end,
            warn_db=args.noise_warn_db, fault_db=args.noise_fault_db,
            diurnal=(t_end - t_start) >= 86400,
            title=(f"{args.array_name}  ·  noise levels  ·  "
                   f"{noise_band_label(args)}  ·  "
                   f"{args.noise_window:g} s windows"))
        fig.savefig(summary_path, dpi=args.figure_dpi)
    elif args.save_figures:
        summary_path = os.path.join(
            args.out_dir, f"{run_tag}_summary.{args.figure_format}")
        try:
            save_summary_figure(summary_path, summary_parts, args,
                                t_start, t_end, families=summary_families)
        except Exception as e:
            print(f"[WARN] summary figure failed: {e}")
            summary_path = None

    what = (f"{sum(len(p) for p in noise_parts)} window levels"
            if args.method == "noise" else f"{total_det} total detections")
    print(f"\n[DONE] {what} across {chunk_i} chunk(s): "
          + ", ".join(f"{n} {k}" for k, n in sorted(status_counts.items())))
    if args.method == "pmcc":
        print(f"       Families: {total_families} → {args.out_dir}/families/")
    if noise_stats_path:
        print(f"       Noise   : {args.out_dir}/noise/  ·  statistics "
              f"{noise_stats_path}")
        for r in run_stats.itertuples():
            print(f"         {r.station:>8}: L90 {r.L90:5.1f}  L50 {r.L50:5.1f}"
                  f"  L10 {r.L10:5.1f}  Leq {r.Leq:5.1f}  offset "
                  f"{r.offset_db:+.2f} dB  {r.flag.upper()}")
    if args.method != "noise":
        print(f"       Archive : {args.out_dir}/detections/")
    print(f"       Index   : {os.path.join(args.out_dir, 'index.csv')}")
    print(f"       Metadata: {os.path.join(args.out_dir, run_tag + '.json')}")
    if args.save_figures:
        print(f"       Figures : {figure_dir}/")
        if summary_path:
            print(f"       Summary : {summary_path}")
    if status_counts["failed"]:
        sys.exit(f"[ERROR] {status_counts['failed']} chunk(s) failed; "
                 "see warnings above.")


if __name__ == "__main__":
    main()
