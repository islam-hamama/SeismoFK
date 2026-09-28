"""
results_io.py — long-term result storage for SeismoFK v1.2.1.

Copyright (c) 2024-2026 Islam Hamama
Contact: islam.hamama@nriag.sci.eg
Licensed under the MIT License.

Storage layout (see also README_v1.2.1.md):

    <out_dir>/
        detections/<run_id>.parquet   # one detection table per processing chunk
        families/<run_id>.parquet     # PMCC families (one row per arrival)
        <run_id>.json                 # run-level metadata + provenance (sidecar)
        index.csv                     # master catalogue, one row per chunk

Why this layout for long-term analysis
---------------------------------------
Detection rows (time, baz, app_vel, semblance, ...) are tabular, append-heavy
and typed — ideal for a *columnar* format.  Parquet stores them ~5-10× smaller
than CSV, preserves dtypes (no string round-tripping of timestamps/floats), and
is fast to scan/filter with pandas or DuckDB even across years of runs:

    import pandas as pd
    df = pd.read_parquet("results/detections")        # whole archive
    df = df[df.semblance > 0.5]                        # cheap column filter

JSON is kept for the small, *nested* run metadata (parameters, station list,
software version) where its flexibility helps and its verbosity does not hurt.
If pyarrow/fastparquet is unavailable we transparently fall back to CSV.
"""

import csv
import json
import os
import tempfile
from datetime import datetime, timezone


def _detections_dir(out_dir, subdir="detections"):
    d = os.path.join(out_dir, subdir)
    os.makedirs(d, exist_ok=True)
    return d


def write_detections(out_dir, run_id, df, prefer_parquet=True,
                     subdir="detections"):
    """
    Persist a detection DataFrame under ``<out_dir>/<subdir>/``.  Returns
    (path, fmt) where fmt is 'parquet' or 'csv'.  Falls back to CSV if no
    Parquet engine is installed.  PMCC families use ``subdir="families"``.
    """
    ddir = _detections_dir(out_dir, subdir)
    if prefer_parquet:
        path = os.path.join(ddir, f"{run_id}.parquet")
        try:
            df.to_parquet(path, index=False)
            return path, "parquet"
        except Exception as e:  # pyarrow/fastparquet missing or write error
            print(f"[WARN] Parquet write failed ({e}); falling back to CSV.")
    path = os.path.join(ddir, f"{run_id}.csv")
    df.to_csv(path, index=False)
    return path, "csv"


def write_metadata(out_dir, run_id, metadata):
    """Write the run-level JSON sidecar.  Returns the path."""
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{run_id}.json")
    meta = dict(metadata)
    meta.setdefault("written_utc", datetime.now(timezone.utc).isoformat())
    with open(path, "w") as fh:
        json.dump(meta, fh, indent=2, default=str)
    return path


def append_index(out_dir, row):
    """
    Append one row to the master ``index.csv`` (creating it with a header on
    first use). Extend the header atomically when a new column appears.
    """
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "index.csv")
    exists = os.path.exists(path)

    # Keep a stable, readable column order; extend with any extra keys.
    base_cols = ["run_id", "array", "method", "chunk_start", "chunk_end",
                 "freq_min", "freq_max", "n_detections", "median_baz",
                 "median_app_vel", "detections_path", "detections_fmt"]
    fieldnames = base_cols + [k for k in row if k not in base_cols]
    if exists:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            header = reader.fieldnames or []
            added = [k for k in fieldnames if k not in header]
            if added:
                fieldnames = header + added
                old_rows = list(reader)
            else:
                fieldnames = header
                old_rows = None
        if old_rows is not None:
            with tempfile.NamedTemporaryFile(
                    mode="w", newline="", dir=out_dir, delete=False) as fh:
                tmp_path = fh.name
                writer = csv.DictWriter(fh, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(old_rows)
            try:
                os.replace(tmp_path, path)
            finally:
                if os.path.exists(tmp_path):
                    os.unlink(tmp_path)
    with open(path, "a", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        if not exists:
            writer.writeheader()
        writer.writerow({k: row.get(k, "") for k in fieldnames})
    return path
