"""
db_manager.py — SQLite persistence layer for SeismoFK event results.

Copyright (c) 2024-2025 Islam Hamama
Contact: islam.hamama@nriag.sci.eg

Licensed under the MIT License — see LICENSE for details.
"""

import os
import sqlite3
from datetime import datetime, timezone

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'fk_events.db')

# ── Classification vocabulary ─────────────────────────────────────────────────
CLASSIFICATIONS = [
    "Unknown",
    "Explosion / Blast",
    "Mining",
    "Volcanic",
    "Earthquake",
    "Meteor / Bolide",
    "Aircraft / Sonic Boom",
    "Ocean / Microbaroms",
    "Industrial",
    "Noise / Artifact",
]

# ── Schema ────────────────────────────────────────────────────────────────────
_SCHEMA = """
CREATE TABLE IF NOT EXISTS events (
    id              INTEGER PRIMARY KEY AUTOINCREMENT,
    saved_at        TEXT    NOT NULL,
    event_name      TEXT,
    classification  TEXT,
    note            TEXT,
    -- array / analysis
    station         TEXT,
    start_time      TEXT,
    duration        REAL,
    min_freq        REAL,
    max_freq        REAL,
    win_length      REAL,
    overlap         REAL,
    semb_thresh     REAL,
    -- FK results
    med_baz         REAL,
    med_vel         REAL,
    expected_baz    REAL,
    n_detections    INTEGER,
    n_windows       INTEGER,
    -- source location
    event_lat       REAL,
    event_lon       REAL,
    -- optional event physics
    origin_time     TEXT,       -- known event origin time (ISO string)
    celerity        REAL,       -- infrasound propagation speed (m/s)
    -- file references
    figure_path     TEXT,
    csv_path        TEXT,
    figure_blob     BLOB        -- PNG bytes of the analysis figure
);
"""


def _connect():
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    return conn


def init_db():
    with _connect() as conn:
        conn.executescript(_SCHEMA)
        # Migrate existing databases that pre-date newer columns
        for col_def in (
            "figure_blob  BLOB",
            "origin_time  TEXT",
            "celerity     REAL",
        ):
            try:
                conn.execute(f"ALTER TABLE events ADD COLUMN {col_def}")
            except sqlite3.OperationalError:
                pass  # column already exists


# ── CRUD ──────────────────────────────────────────────────────────────────────

def save_event(data: dict) -> int:
    """Insert a new event row and return the new id."""
    init_db()
    data = dict(data)                           # don't mutate caller's dict
    data.setdefault('saved_at', datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S'))
    cols   = ', '.join(data.keys())
    placeh = ', '.join('?' * len(data))
    with _connect() as conn:
        cur = conn.execute(f"INSERT INTO events ({cols}) VALUES ({placeh})",
                           list(data.values()))
        return cur.lastrowid


def fetch_all() -> list:
    init_db()
    with _connect() as conn:
        return conn.execute(
            "SELECT * FROM events ORDER BY saved_at DESC"
        ).fetchall()


def fetch_by_classification(cls: str) -> list:
    init_db()
    with _connect() as conn:
        return conn.execute(
            "SELECT * FROM events WHERE classification = ? ORDER BY saved_at DESC",
            (cls,)
        ).fetchall()


def update_event(event_id: int, classification: str, note: str):
    with _connect() as conn:
        conn.execute(
            "UPDATE events SET classification=?, note=? WHERE id=?",
            (classification, note, event_id)
        )


def delete_event(event_id: int):
    with _connect() as conn:
        conn.execute("DELETE FROM events WHERE id=?", (event_id,))


def export_csv(path: str):
    """Dump entire events table to a CSV file."""
    import csv
    rows = fetch_all()
    if not rows:
        return 0
    with open(path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows([dict(r) for r in rows])
    return len(rows)
