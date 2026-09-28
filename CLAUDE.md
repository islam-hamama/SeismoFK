# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

SeismoFK v1.2.1 is an infrasound/seismic array-analysis tool: a PyQt5 desktop GUI plus a batch CLI, sharing one set of numerical modules. `AGENTS.md` holds the repository guidelines (style, commits, PRs); this file doesn't repeat them.

## Commands

Run everything from the repository root. The modules are flat top-level files that import each other by name; `pyproject.toml` lists them as `py-modules`, so add any new module there too. The previous release is the `v1.2.0` git tag.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt          # or: pip install -e ".[parquet]"
python _selftest.py                      # synthetic plane-wave check of FK/ARF/bootstrap/PMCC
python -m unittest test_release.py       # release checks (archive integrity, CLI validation, readiness)
python -m unittest test_release.ReleaseChecks.test_cli_archives_synthetic_pmcc_run   # single test
python Infra_Analysis.py                 # GUI (entry point: seismofk)
python seismofk_cli.py --help            # batch CLI (entry point: seismofk-cli)
```

There is no linter, formatter or build step. `_selftest.py` asserts that every method recovers back-azimuth 120° (±5°) and 340 m/s (±40 m/s). Run it after any change to `fk_analysis.py` or `pmcc.py`.

## Architecture

**Analysis core (no Qt imports):**
- `fk_analysis.py` is the base module. `station_xy_km` projects station coordinates to local flat-Earth x/y in km. The obspy-based `fk_array` is the classic FK that the GUI uses. `fk_beamform` (Bartlett/Capon/MUSIC over a CSDM from `cross_spectral_matrix`) and its sliding-window wrapper `fk_scan` are what the CLI uses. `bootstrap_beam` and `array_response` are also here.
- `pmcc.py` is a PMCC-style detector. It measures pairwise delays per log-spaced band (`log_bands`), keeps triplet-consistent pixels and solves for slowness by least squares. It imports `station_xy_km` from `fk_analysis`. `pmcc_families(pixels, time_tol=window)` then groups pixels into families (one per arrival) by single linkage: |Δt| ≤ window, bands at most one apart, |Δbaz| ≤ 10°, |Δv| ≤ 15%, and at least 5 pixels. Uncertainty model: each sensor's arrival time has an error σₐ, so Cov(s) = σₐ²(XcᵀXc)⁻¹ with Xc the centred positions; the naive σ²(GᵀG)⁻¹ is N/2 too small. σₐ² is the largest of the fit residual RSS/(N(N−3)), half the delay Cramér–Rao bound computed from the overlap-corrected correlation, and a timing floor. Families combine pixels by inverse-variance weighting, inflated by `pixel_overlap` (window/hop) and the Birge ratio. `mixed` is a split-half direction-change test, not the χ² p-value, which rejects most real arrivals. `test_pmcc_uncertainty_covers_truth` guards the calibration, so re-run it after touching `_pair_delays`, `_delay_crlb`, `_slowness_uncertainty` or `_combine_slowness`. Families only form reliably at hop = window/4, which is the GUI and CLI default; window/2 leaves too few pixels per arrival. The CLI forms families per chunk, so an arrival that straddles a chunk boundary can be split.
- `data_readiness.assess_stream` runs pre-flight checks (collinear geometry, band vs Nyquist, missing metadata) and returns `Finding` objects with a severity. The GUI calls it before processing.

**Units:** slowness in s/km, positions in km, `app_vel` in m/s in outputs, back-azimuth in degrees. Azimuth statistics must wrap correctly at north, and the tests check this.

**Coordinates and epochs:** IMS StationXML (for example `XML_IM/I48TN.xml`) contains a long nominal epoch that overlaps the dated calibration history. For I48 in 2026 that means 10000 counts/Pa where the newer epoch says 7194–7874 counts/Pa. ObsPy's `get_response`, `get_coordinates` and `remove_response` silently use the *first* match. So every pipeline (GUI `ProcessThread` and `_resolve_inventory`, CLI per chunk, readiness) first calls `fk_analysis.resolve_epochs(inv, stream)`, which keeps one epoch per trace, the one with the latest start. Do the same in any new code path before reading coordinates or responses. IMS responses also often have no stages, only a sensitivity. `require_response_stages` turns ObsPy's IndexError into a clear message before the scalar-sensitivity fallback.

**Two front ends, one pipeline.** `Infra_Analysis.py` (`ProcessThread`, a QThread) and `seismofk_cli.py` (`read_window` → `preprocess` → `analyse_chunk`) both do the same steps: load MiniSEED, remove the response with `pre_filt=(0.1, 0.5, 9.0, 10.0)` (falling back to scalar sensitivity), then run the analysis. The CLI comment says it "mirrors the GUI ProcessThread pipeline", so keep the two in sync when you change preprocessing. Both get output units from `fk_analysis.calibrated_units`: ObsPy applies a Pa response as-is even with `output='VEL'`, so infrasound stays in Pa. Never hard-code the units. Both also run `data_readiness.assess_stream`: the GUI through **Check data readiness**, and the CLI as a preflight on the first chunk (`--check-only`, `--skip-checks`).

**GUI layering:** `Infra_Analysis.py` holds the main window, results, spectrogram, waveform viewer and event-database dialogs. `gui_methods.py` holds the matplotlib-in-Qt dialogs for ARF, slowness map and PMCC. `result_ui.py` is the single theme for every window, including the StationXML editor. It defines colours once in `PALETTE`, and the one style sheet comes from it (`apply_main_style` for main windows, `apply_result_style` for dialogs). `apply_app_theme(app)` in each `main()` sets the platform UI font: SF Pro on macOS, Segoe UI on Windows, or Inter if installed. Widgets choose a role with `setObjectName`: `primaryAction`, `methodAction`, `dangerAction`, `heroAction`, `resultCard`, `fieldLabel`, `mutedText`. Don't add inline colours or emoji button labels. Keep font weights at 500–600; 700–800 render as Heavy/Black in system fonts. Figures get the matching typography from `result_plots.apply_plot_style()`, which runs when `Infra_Analysis` is imported and in the CLI when saving figures. Embed plots with `matplotlib.figure.Figure`, not `plt.subplots`, so figures aren't leaked through pyplot. `plot_spectrogram_window.compute_spectrogram` is kept free of the GUI so tests can import it. `result_plots.py` is Qt-free too. It draws the PMCC and FK-scan figures (waveform, correlation or peak/mean power, azimuth, velocity over the full interval) into any `Figure`, and both `gui_methods.PMCCWindow` and `seismofk_cli --save-figures` use it. `fk_results_figure` is the six-panel conventional-FK figure. The GUI's `_show_results` and `seismofk_cli --method fk` both draw it from the `fk_array` result dict. Change a figure there, not in two places. PMCC's `diagnostics['windows']` includes rejected windows, which provide the noise floor on the correlation panel. The module deliberately doesn't call `matplotlib.use` at import time.

**Noise levels:** `noise_levels.py` is Qt-free and shared by `gui_methods.NoiseLevelsDialog` and `seismofk_cli --method noise`. Windows lie on an absolute UTC grid (multiples of the window), so chunking never changes the results. Leq averages power, and L90/L50/L10 are exceedance percentiles. Sensor offsets are measured against the median of the *other* sensors; including the sensor itself hides part of its own fault. Filtering uses padded `sosfiltfilt`, not a taper, to keep the first and last windows exact, and zero-filled gaps are excluded through `coverage`. Calibration uses the StationXML sensitivity, with no response pre-filter that would cut the band.

**Persistence:**
- `results_io.py` writes CLI archives to `<out_dir>/detections/<run_id>.parquet` (CSV if pyarrow is missing or `--csv` is given), with a `<run_id>.json` metadata sidecar and an `index.csv` catalogue. Writes are atomic, and the index header can grow without dropping rows. Run IDs come from `make_run_tag` and are unique per invocation.
- `db_manager.py` is the SQLite event catalogue used by the GUI. It migrates the schema with `ALTER TABLE ADD COLUMN` inside `init_db`. `app_paths.database_path()` prefers a legacy `fk_events.db` in the checkout; otherwise it uses `~/.seismofk/` (override with `SEISMOFK_DATA_DIR`).

**CLI validation:** invalid parameters, such as `--chunk 0`, `--overlap 1`, `fmin >= fmax` or a path-traversing `--array-name`, must exit with code 2 before any file is read. Blocking readiness findings also exit with code 2, before anything is written. `read_window` returns `(stream, gap_fraction)`, and chunks above `--max-gap` are skipped because zero-filled gaps are coherent across sensors. Each chunk is analysed one window past its end, and `rows_in_chunk` assigns each window to the chunk where it starts, so boundary windows are not lost. Every chunk writes an `index.csv` row with a `status` (`ok`, `no_detections`, `no_data`, `gap_skipped`, `failed`) plus `units`, `calibration`, `gap_fraction` and `n_channels`. FK `coherence` is peak/mean slowness power (≥ 1); `--min-coherence` gates it. `test_release.py` enforces all of this.

## Data and generated files

- StationXML inventories go in `XML/` or `XML_IM/`. Read `XML/README.md` before adding metadata. `xml_creator.py`, `convert_ims.py` and `export_stations.py` are the metadata utilities.
- A GUI run writes `FK_<name>.jpg` and `Output_<name>.csv` to the current working directory. These, `fk_events.db`, personal `XML/` inventories and exported `.npz` files are ignored by `.gitignore`. `XML_IM/`, the IMS StationXML inventories, is tracked and published (decided for v1.2.1).
