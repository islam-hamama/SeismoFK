# Repository Guidelines

## Project Structure & Module Organization

This release is a flat Python project. `Infra_Analysis.py` starts the PyQt5 GUI; `seismofk_cli.py` runs long-term batch analysis. Array algorithms live in `fk_analysis.py` and `pmcc.py`, with output handling in `results_io.py` and database access in `db_manager.py`. Supporting GUI and utility scripts include `gui_methods.py`, `xml_creator.py`, `convert_ims.py`, and `plot_spectrogram_window.py`. StationXML inventories belong in `XML/` or `XML_IM/`; consult `XML/README.md` before adding metadata. `_selftest.py` is the current synthetic analysis check. User instructions are in `README.md` and `USAGE.md`.

## Build, Test, and Development Commands

There is no build step. From this directory, create a virtual environment and install the declared dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python _selftest.py
python Infra_Analysis.py
python seismofk_cli.py --help
```

`_selftest.py` exercises FK, array response, bootstrap, and PMCC on a synthetic plane wave. The GUI command opens the desktop application; the CLI help command lists batch options. A real batch run needs MiniSEED waveforms and a StationXML inventory; see `README.md` for an example.

## Coding Style & Naming Conventions

Follow the existing Python style: four-space indentation, `snake_case` functions and modules, `PascalCase` classes, and descriptive constants. Keep numerical analysis in reusable functions rather than GUI callbacks, and document units for slowness, velocity, frequency, and time. No formatter or linter configuration is committed, so match nearby code when editing.

## Testing Guidelines

Run `python _selftest.py` after changing array algorithms and inspect its reported back-azimuth and velocity against the synthetic target of 120° and 340 m/s. The repository has no configured pytest suite or coverage threshold. For GUI, XML, or archive changes, manually exercise the affected workflow with representative local data and describe that check in the pull request.

## Commit & Pull Request Guidelines

Recent commit subjects are short, imperative descriptions such as `Add screenshots: ...` and `Update README.md`. Use the same style and keep each commit focused. Pull requests should explain the behavior changed, list commands or manual checks performed, link relevant issues when available, and include screenshots for visible GUI changes. Do not commit private station metadata or generated analysis outputs without confirming redistribution rights.
