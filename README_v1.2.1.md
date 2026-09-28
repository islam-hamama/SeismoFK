# SeismoFK v1.2.1 — Enhanced Array Analysis & Long-Term Batch CLI

Release notes and archive-format details for v1.2.1: the updated
`Infra_Analysis.py` GUI, analysis modules, and batch CLI. Run commands from the
repository root after installing `requirements.txt`. The previous release is
available as the `v1.2.0` tag.

## What's new in v1.2.1

| Feature | Module | Notes |
|---|---|---|
| **Capon (MVDR) & MUSIC** high-resolution beamforming | `fk_analysis.fk_beamform` | Sharper slowness peaks than conventional FK; resolves close/multiple arrivals. Bartlett also available for comparison. |
| **Windowed FK scanner** | `fk_analysis.fk_scan` | One detection row per sliding window — workhorse for batch runs. |
| **Bootstrap uncertainty** | `fk_analysis.bootstrap_beam` | `baz ± σ`, `vel ± σ`, confidence intervals, and a slowness-space confidence ellipse. |
| **Array Response Function** | `fk_analysis.array_response` | Theoretical ARF for the current geometry (aliasing & resolution). |
| **PMCC detector + families** | `pmcc.pmcc`, `pmcc.pmcc_families` | Multi-band, triplet-consistency detection pixels and least-squares slowness, aggregated into families (one per coherent arrival) by linking pixels adjacent in time and frequency with matching back-azimuth (±10°) and trace velocity (±15%). Pixels and families carry 95% confidence intervals (arrival-time error model with a Cramér–Rao floor, validated by a synthetic coverage test), and families that mix two directions are flagged. |
| **Long-term batch CLI** | `seismofk_cli.py` | Processes long time ranges in chunks; archives detections to Parquet (or CSV) with JSON metadata and `index.csv`. |
| **Parametrised classic FK** | `fk_analysis.fk_array` | Slowness grid & prewhitening are now arguments (defaults reproduce v1.2.0 exactly). |

All routines are verified against a synthetic plane wave in `_selftest.py`
(`python _selftest.py`): every method recovers the true back-azimuth (120°) and
trace velocity (340 m/s).
Run `python -m unittest test_release.py` for archive and CLI safety checks.

## Install

```bash
pip install -r requirements.txt      # core dependencies
pip install ".[parquet]"             # optional: pyarrow for Parquet archives
```

## Command-line usage (long-term analysis)

```bash
# One month of data, Capon high-resolution FK, hourly chunks, with bootstrap CI
python seismofk_cli.py \
    --mseed "data/2025-05/*.mseed" --inventory XML/ \
    --start 2025-05-01T00:00:00 --end 2025-06-01T00:00:00 \
    --method capon --fmin 0.5 --fmax 4.0 \
    --win-length 30 --overlap 0.5 --chunk 3600 \
    --array-name I31 --out-dir results/ --bootstrap

# PMCC-style detection pixels with 6 log-spaced bands
python seismofk_cli.py --mseed data/ --inventory XML/ \
    --start 2025-05-01 --end 2025-05-02 \
    --method pmcc --fmin 0.1 --fmax 8.0 --pmcc-bands 6 \
    --win-length 30 --array-name I31 --out-dir results/
```

`--method` accepts `bartlett`, `capon`, `music`, or `pmcc`.
Run `python seismofk_cli.py --help` for the full option list.

## Output layout & the JSON-vs-CSV question

For long-term work you accumulate **millions of typed, append-heavy detection
rows** (time, baz, app_vel, slowness, semblance/coherence, …). That is exactly
what a **columnar** format is built for, so the recommended layout is
**Parquet for the detection tables + a small JSON sidecar for run metadata**:

```
results/
├── detections/<run_id>.parquet   # detection table per chunk (typed, ~5-10× smaller than CSV)
├── families/<run_id>.parquet  # PMCC families (one row per arrival)
├── <run_id>.json                 # run metadata: array, band, params, channels, version
└── index.csv                     # master catalogue, one row per chunk (human-browsable)
```

**Why this beats the alternatives**

- **Parquet vs CSV** — Parquet is typically 5–10× smaller, preserves dtypes
  (timestamps/floats are not re-parsed from strings), and supports fast
  *column* filtering, so querying years of detections stays cheap:

  ```python
  import pandas as pd
  df = pd.read_parquet("results/detections")   # whole archive at once
  df = df[df.app_vel.between(280, 450)]         # cheap columnar filter
  ```

- **Parquet vs pure JSON** — JSON has no columnar layout; a single file holding
  millions of detection objects is bloated and slow to scan. JSON *is* the right
  tool for the small, nested **run metadata** (parameters, station list,
  software version), which is why it's kept as a per-run sidecar.

- **`index.csv`** stays CSV on purpose — it's tiny (one row per chunk) and you'll
  often open it by eye or in a spreadsheet.

If `pyarrow` is not installed, the CLI automatically falls back to CSV for the
detection tables (or force it with `--csv`).
Each CLI invocation gets a unique run ID, so repeating the same arguments keeps
earlier outputs. A run with failed chunks exits nonzero; completed chunks remain
in the archive.

## Library usage

```python
import fk_analysis as fk
from obspy import read, read_inventory

st = read("array.mseed"); inv = read_inventory("inv.xml")
fk.attach_coordinates(st, inv)

# High-resolution slowness map
sx, sy, P = fk.fk_beamform(st, fmin=0.5, fmax=4.0, method="capon")
peak = fk.estimate_slowness_peak(sx, sy, P)
print(peak["baz"], peak["app_vel"])

# Array response function
sx, sy, arf = fk.array_response(st, 0.5, 4.0)
```
```python
import pmcc
bands = pmcc.log_bands(0.1, 8.0, 6)
det = pmcc.pmcc(st, bands, window_sec=30)   # DataFrame of detection pixels
```
