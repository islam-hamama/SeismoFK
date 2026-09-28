"""Preflight checks for array waveforms and StationXML metadata.

Copyright (c) 2024-2026 Islam Hamama
Contact: islam.hamama@nriag.sci.eg
Licensed under the MIT License.
"""

from dataclasses import dataclass

import numpy as np
from obspy.core.util import AttribDict

import fk_analysis as fk


@dataclass(frozen=True)
class Finding:
    severity: str
    title: str
    detail: str


def assess_stream(stream, inventory, *, fmin, fmax, start=None, duration=None):
    """Return actionable findings without changing waveform or inventory data.

    Errors block FK/PMCC use. Warnings describe calibration, coverage, or
    signal-quality conditions that merit inspection before interpretation.
    """
    findings = []

    def add(severity, title, detail):
        findings.append(Finding(severity, title, detail))

    if not stream:
        add("error", "No waveform loaded", "Select a MiniSEED file first.")
        return findings

    if len(stream) < 3:
        add("error", "Too few sensors", "Array methods need at least three traces.")
    ids = [tr.id for tr in stream]
    if len(set(ids)) != len(ids):
        add("warning", "Repeated channel IDs",
            "Merge or select one trace per channel before array processing.")

    rates = np.array([float(tr.stats.sampling_rate) for tr in stream])
    if np.any(~np.isfinite(rates)) or np.any(rates <= 0):
        add("error", "Invalid sample rate", "Every trace needs a positive sample rate.")
    else:
        if len(set(np.round(rates, 6))) > 1:
            add("warning", "Mixed sample rates",
                f"Rates span {rates.min():g}–{rates.max():g} Hz; resampling is required.")
        if fmin <= 0 or fmax <= fmin:
            add("error", "Invalid frequency band", "Set 0 < minimum < maximum frequency.")
        elif fmax >= rates.min() / 2:
            add("error", "Band exceeds Nyquist",
                f"Maximum frequency must be below {rates.min() / 2:g} Hz for the slowest trace.")

    common_start = max(tr.stats.starttime for tr in stream)
    common_end = min(tr.stats.endtime for tr in stream)
    if common_end <= common_start:
        add("error", "No common recording interval",
            "Choose traces whose time ranges overlap.")
    elif start is not None and duration is not None:
        requested_end = start + float(duration)
        if requested_end <= common_start or start >= common_end:
            add("error", "Selected range misses the data",
                "Pick a start time inside the common recording interval.")
        elif start < common_start or requested_end > common_end:
            add("warning", "Selected range extends beyond data",
                "Part of the requested duration has no samples on every trace.")

    for tr in stream:
        data = np.asarray(tr.data)
        if tr.stats.npts < 32:
            add("warning", f"Short trace · {tr.id}", "Fewer than 32 samples are available.")
        elif not np.all(np.isfinite(data)):
            add("warning", f"Nonfinite samples · {tr.id}",
                "Repair or exclude NaN/Inf samples before analysis.")
        elif np.ptp(data) == 0:
            add("warning", f"Flat trace · {tr.id}",
                "A constant channel cannot contribute timing information.")

    if inventory is None:
        add("error", "No StationXML selected",
            "Select matching station metadata for FK and PMCC. Spectrograms can use raw counts.")
        return findings

    from plot_spectrogram_window import get_sensitivity_counts_per_pa

    overlapping = [(tid, c, k) for tid, (_, _, c, k)
                   in fk.epoch_choices(inventory, stream).items() if k > 1]
    if overlapping:
        chosen = "\n".join(
            f"{tid}: {k} epochs, using {c.start_date.date}"
            + (f" ({c.response.instrument_sensitivity.value:g} "
               f"{c.response.instrument_sensitivity.output_units or ''}/"
               f"{c.response.instrument_sensitivity.input_units or ''})"
               if c.response and c.response.instrument_sensitivity else "")
            for tid, c, k in overlapping)
        add("warning", "Overlapping StationXML epochs",
            "Several channel epochs cover the data; SeismoFK uses the most "
            f"recent calibration. Check this is intended.\n{chosen}")
    inventory = fk.resolve_epochs(inventory, stream)

    located = stream.copy()
    positions = []
    for tr in located:
        try:
            try:
                inventory.get_coordinates(tr.id, datetime=tr.stats.starttime)
            except Exception:
                add("warning", f"Coordinate fallback · {tr.id}",
                    "Exact channel and epoch metadata were unavailable; check the selected inventory.")
            coords = fk.get_coordinates_safe(
                inventory, tr.id, datetime=tr.stats.starttime, verbose=False)
            tr.stats.coordinates = AttribDict(coords)
            positions.append(tr)
        except Exception:
            add("error", f"Missing coordinates · {tr.id}",
                "Add matching channel coordinates to StationXML.")

        try:
            get_sensitivity_counts_per_pa(inventory, tr)
        except Exception:
            add("warning", f"Pressure calibration unavailable · {tr.id}",
                "Spectrograms will use raw counts for this array.")

    if len(positions) >= 3:
        xy, _, _ = fk.station_xy_km(positions)
        if np.linalg.matrix_rank(xy - xy.mean(axis=0), tol=1e-3) < 2:
            add("error", "Collinear array geometry",
                "PMCC direction and velocity need non-collinear sensor coordinates.")

    if not findings:
        add("info", "Ready", "Waveforms, timing, geometry, and metadata passed the checks.")
    return findings
