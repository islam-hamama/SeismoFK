"""Writable application paths shared by the GUI and database layer.

Copyright (c) 2024-2026 Islam Hamama
Contact: islam.hamama@nriag.sci.eg
Licensed under the MIT License.
"""

import os
from pathlib import Path


def data_dir(create=False):
    """Return the user's SeismoFK directory (override with SEISMOFK_DATA_DIR)."""
    root = Path(os.environ.get("SEISMOFK_DATA_DIR", "~/.seismofk")).expanduser()
    if create:
        root.mkdir(parents=True, exist_ok=True)
    return root


def inventory_dir(create=False):
    path = data_dir(create=create) / "XML"
    if create:
        path.mkdir(parents=True, exist_ok=True)
    return path


def database_path():
    """Keep an existing checkout database accessible; use user data for new ones."""
    legacy = Path(__file__).resolve().parent / "fk_events.db"
    if legacy.is_file():
        return legacy
    return data_dir() / "fk_events.db"
