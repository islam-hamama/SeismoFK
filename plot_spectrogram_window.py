"""
plot_spectrogram_window.py — Spectrogram plotting utility for SeismoFK.

Copyright (c) 2024-2025 Islam Hamama
Contact: islam.hamama@nriag.sci.eg

Licensed under the MIT License — see LICENSE for details.
"""
import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path("tmp_mplconfig").resolve()))

# NOTE: matplotlib.use("Agg") is intentionally NOT called at import time.
# This module is imported by the PyQt5 GUI (Infra_Analysis.py), which selects
# the "Qt5Agg" backend. Forcing "Agg" here would break interactive GUI canvases.
# The non-interactive "Agg" backend is selected only on the CLI path — see main().

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
from obspy import UTCDateTime, read, read_inventory
from scipy.signal import spectrogram


# Default configuration values. Paths are deliberately neutral placeholders —
# the script must not ship frozen to one private dataset. Supply your own
# waveform and inventory files via the command-line flags (see --help).
DEFAULT_OUTPUT_DIR = "spectrogram_exports"
DEFAULT_WINDOW_START = "2026-04-01T23:40:00"
DEFAULT_WINDOW_END = "2026-04-02T00:20:00"
DEFAULT_DPI = 300
DEFAULT_REFERENCE_PRESSURE_PA = 20e-6
DEFAULT_FREQ_MAX_HZ = 5.0
DEFAULT_NPERSEG = 512
DEFAULT_NOVERLAP = 460
DEFAULT_FILTER_FREQMIN = 0.5
DEFAULT_FILTER_FREQMAX = 6.0


def parse_args(argv=None):
    """Parse command-line arguments for the spectrogram plotting utility."""
    parser = argparse.ArgumentParser(
        description=(
            "Plot infrasound/seismic spectrograms and filtered waveforms "
            "for a time window of an array dataset."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-i", "--input-file", required=True,
        help="Waveform file to read (any format ObsPy can read, e.g. MiniSEED).",
    )
    parser.add_argument(
        "-x", "--inventory-file", required=True,
        help="StationXML inventory file providing instrument response.",
    )
    parser.add_argument(
        "-o", "--output-dir", default=DEFAULT_OUTPUT_DIR,
        help="Directory where the spectrogram and waveform PNGs are written.",
    )
    parser.add_argument(
        "--window-start", default=DEFAULT_WINDOW_START,
        help="UTC start time of the analysis window (ISO 8601).",
    )
    parser.add_argument(
        "--window-end", default=DEFAULT_WINDOW_END,
        help="UTC end time of the analysis window (ISO 8601).",
    )
    parser.add_argument(
        "--dpi", type=int, default=DEFAULT_DPI,
        help="Resolution (dots per inch) of the exported figures.",
    )
    parser.add_argument(
        "--reference-pressure", type=float, default=DEFAULT_REFERENCE_PRESSURE_PA,
        help="Reference pressure in Pa for dB-scaled PSD (20 uPa is standard).",
    )
    parser.add_argument(
        "--freq-max", type=float, default=DEFAULT_FREQ_MAX_HZ,
        help="Maximum frequency (Hz) displayed on the spectrogram.",
    )
    parser.add_argument(
        "--nperseg", type=int, default=DEFAULT_NPERSEG,
        help="Spectrogram segment length (samples) passed to scipy.signal.spectrogram.",
    )
    parser.add_argument(
        "--noverlap", type=int, default=DEFAULT_NOVERLAP,
        help="Spectrogram segment overlap (samples) passed to scipy.signal.spectrogram.",
    )
    parser.add_argument(
        "--filter-freqmin", type=float, default=DEFAULT_FILTER_FREQMIN,
        help="Lower corner frequency (Hz) of the bandpass filter.",
    )
    parser.add_argument(
        "--filter-freqmax", type=float, default=DEFAULT_FILTER_FREQMAX,
        help="Upper corner frequency (Hz) of the bandpass filter.",
    )
    return parser.parse_args(argv)


def utc_to_matplotlib(times_sec, window_start):
    return [
        (window_start + float(offset)).datetime.replace(tzinfo=None)
        for offset in times_sec
    ]


def get_sensitivity_counts_per_pa(inventory, trace):
    response = inventory.get_response(trace.id, trace.stats.starttime)
    sensitivity = response.instrument_sensitivity
    if sensitivity is None:
        raise ValueError("No instrument sensitivity found for {0}".format(trace.id))
    if sensitivity.input_units.upper() != "PA" or sensitivity.output_units.upper() != "COUNTS":
        raise ValueError(
            "Unexpected sensitivity units for {0}: {1} -> {2}".format(
                trace.id, sensitivity.input_units, sensitivity.output_units
            )
        )
    return float(sensitivity.value)


def compute_spectrogram(
    trace,
    inventory,
    *,
    filter_freqmin,
    filter_freqmax,
    freq_max,
    nperseg,
    noverlap,
    reference_pressure=DEFAULT_REFERENCE_PRESSURE_PA,
):
    """Compute a dB-scaled PSD spectrogram for a single ObsPy trace.

    This is a pure, importable function — no plotting, no file I/O. Both the
    CLI ``main()`` and the SeismoFK GUI (``SpectrogramWindow``) call it.

    Pipeline (identical to the original standalone script):
      1. detrend("linear") -> detrend("demean")
      2. counts -> Pa using the instrument sensitivity from ``inventory``
      3. zero-phase 4-corner bandpass [filter_freqmin, filter_freqmax]
      4. scipy.signal.spectrogram (Hann window, density/PSD)
      5. crop to ``freq <= freq_max`` and convert to dB re ``reference_pressure``

    Parameters
    ----------
    trace : obspy.Trace
        A single trace, ALREADY sliced to the desired time window by the
        caller. (The CLI slices before calling; the GUI passes a windowed
        stream.) Must contain at least one sample.
    inventory : obspy.Inventory
        Inventory providing the instrument response for ``trace.id``.
    filter_freqmin, filter_freqmax : float
        Bandpass corner frequencies (Hz). ``filter_freqmin`` must be > 0 and
        strictly less than ``filter_freqmax``.
    freq_max : float
        Upper frequency bound (Hz) retained in the returned arrays.
    nperseg, noverlap : int
        Segment length / overlap (samples) for ``scipy.signal.spectrogram``.
        ``noverlap`` must be < ``nperseg``.
    reference_pressure : float, optional
        Reference pressure (Pa) for the dB conversion. Defaults to 20 µPa.

    Returns
    -------
    (freqs, time_nums, power_db) : tuple of numpy.ndarray
        ``freqs``      — frequency bin centres (Hz), cropped to ``freq_max``.
        ``time_nums``  — segment time centres as matplotlib date numbers,
                         anchored to ``trace.stats.starttime``.
        ``power_db``   — 2-D PSD in dB, shape ``(len(freqs), len(time_nums))``.

    Raises
    ------
    ValueError
        On an empty trace, invalid filter band, or invalid spectrogram
        segmentation parameters.
    """
    if trace.stats.npts == 0:
        raise ValueError(
            "Trace {0} has no samples in the requested window.".format(trace.id)
        )
    if not (filter_freqmin > 0):
        raise ValueError("filter_freqmin must be > 0 Hz.")
    if not (filter_freqmin < filter_freqmax):
        raise ValueError("filter_freqmin must be < filter_freqmax.")
    if not (freq_max > 0):
        raise ValueError("freq_max must be > 0 Hz.")
    if nperseg <= 0:
        raise ValueError("nperseg must be a positive integer.")
    if noverlap < 0 or noverlap >= nperseg:
        raise ValueError("noverlap must satisfy 0 <= noverlap < nperseg.")

    work = trace.copy()
    work.detrend("linear")
    work.detrend("demean")

    data_counts = work.data.astype(np.float64)
    sensitivity_counts_per_pa = get_sensitivity_counts_per_pa(inventory, work)
    data_pa = data_counts / sensitivity_counts_per_pa

    fs = work.stats.sampling_rate
    nyquist = fs / 2.0
    if filter_freqmax >= nyquist:
        raise ValueError(
            "filter_freqmax ({0} Hz) must be below the Nyquist "
            "frequency ({1} Hz) for fs={2} Hz.".format(
                filter_freqmax, nyquist, fs)
        )

    filtered = work.copy()
    filtered.data = data_pa.copy()
    filtered.filter(
        "bandpass",
        freqmin=filter_freqmin,
        freqmax=filter_freqmax,
        corners=4,
        zerophase=True,
    )
    effective_nperseg = min(nperseg, filtered.stats.npts)
    effective_noverlap = min(noverlap, max(effective_nperseg - 1, 0))

    freqs, times, power = spectrogram(
        filtered.data.astype(np.float64),
        fs=fs,
        window="hann",
        nperseg=effective_nperseg,
        noverlap=effective_noverlap,
        detrend=False,
        scaling="density",
        mode="psd",
    )

    band_mask = freqs <= freq_max
    if not band_mask.any():
        raise ValueError(
            "freq_max ({0} Hz) is below the lowest spectrogram "
            "frequency bin ({1:.3f} Hz).".format(freq_max, float(freqs[0]))
        )
    freqs = freqs[band_mask]
    power = power[band_mask, :]

    power_db = 10 * np.log10(
        power / (reference_pressure ** 2) + np.finfo(float).eps
    )
    time_nums = mdates.date2num(
        utc_to_matplotlib(times, work.stats.starttime)
    )
    return freqs, time_nums, power_db


def style_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, color="white", linewidth=0.6, alpha=0.25)


def main(argv=None):
    # The CLI path is headless — select the non-interactive Agg backend here,
    # not at module import time, so GUI use (Qt5Agg) is not disturbed.
    import matplotlib
    matplotlib.use("Agg")

    args = parse_args(argv)

    input_file = Path(args.input_file)
    inventory_file = Path(args.inventory_file)
    output_dir = Path(args.output_dir)
    window_start = UTCDateTime(args.window_start)
    window_end = UTCDateTime(args.window_end)
    dpi = args.dpi
    reference_pressure_pa = args.reference_pressure
    freq_max_hz = args.freq_max
    nperseg = args.nperseg
    noverlap = args.noverlap
    filter_freqmin = args.filter_freqmin
    filter_freqmax = args.filter_freqmax

    stream = read(str(input_file))
    inventory = read_inventory(str(inventory_file))
    output_dir.mkdir(parents=True, exist_ok=True)
    Path("tmp_mplconfig").mkdir(exist_ok=True)

    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "#f6f7f9",
            "axes.edgecolor": "#4a4a4a",
            "axes.labelcolor": "#202124",
            "xtick.color": "#202124",
            "ytick.color": "#202124",
            "font.size": 10,
            "axes.titlesize": 13,
            "axes.labelsize": 10,
        }
    )

    panel_data = []
    global_min = []
    global_max = []

    waveform_panels = []

    for trace in stream:
        sliced = trace.copy().slice(window_start, window_end)
        if sliced.stats.npts == 0:
            continue

        # ── Spectrogram PSD via the shared pure function ──────────────────
        freqs, time_nums, power_db = compute_spectrogram(
            sliced,
            inventory,
            filter_freqmin=filter_freqmin,
            filter_freqmax=filter_freqmax,
            freq_max=freq_max_hz,
            nperseg=nperseg,
            noverlap=noverlap,
            reference_pressure=reference_pressure_pa,
        )
        panel_data.append(
            {
                "trace": trace,
                "freqs": freqs,
                "time_nums": time_nums,
                "power_db": power_db,
            }
        )
        global_min.append(np.percentile(power_db, 15))
        global_max.append(np.percentile(power_db, 98))

        # ── Filtered-waveform panel (in Pa, same band) ────────────────────
        detrended = sliced.copy()
        detrended.detrend("linear")
        detrended.detrend("demean")
        sensitivity_counts_per_pa = get_sensitivity_counts_per_pa(
            inventory, detrended
        )
        filtered = detrended.copy()
        filtered.data = (
            detrended.data.astype(np.float64) / sensitivity_counts_per_pa
        )
        filtered.filter(
            "bandpass",
            freqmin=filter_freqmin,
            freqmax=filter_freqmax,
            corners=4,
            zerophase=True,
        )
        waveform_panels.append(
            {
                "trace": trace,
                "times": mdates.date2num(
                    utc_to_matplotlib(filtered.times(), filtered.stats.starttime)
                ),
                "data_pa": filtered.data.astype(np.float64),
            }
        )

    fig = plt.figure(figsize=(15, 10.5), facecolor="white")
    gs = GridSpec(
        len(panel_data),
        2,
        width_ratios=[40, 1.8],
        hspace=0.12,
        wspace=0.08,
        left=0.07,
        right=0.92,
        top=0.9,
        bottom=0.09,
    )

    axes = []
    cax = fig.add_subplot(gs[:, 1])
    vmin = min(global_min)
    vmax = max(global_max)
    last_mesh = None

    for idx, panel in enumerate(panel_data):
        ax = fig.add_subplot(gs[idx, 0], sharex=axes[0] if axes else None)
        axes.append(ax)
        trace = panel["trace"]

        last_mesh = ax.pcolormesh(
            panel["time_nums"],
            panel["freqs"],
            panel["power_db"],
            shading="auto",
            cmap="magma",
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_ylabel("{0}\nHz".format(trace.stats.station), rotation=0, labelpad=28)
        ax.set_ylim(0, freq_max_hz)
        ax.set_yticks(list(range(0, int(freq_max_hz) + 1)))
        style_axes(ax)
        ax.text(
            0.015,
            0.92,
            "{0}".format(trace.id),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            color="white",
            bbox=dict(facecolor="black", alpha=0.28, edgecolor="none", pad=3.0),
        )

    # Derive a human-readable label and time strings from the data / args
    # so the figures are not frozen to one named dataset.
    array_label = panel_data[0]["trace"].stats.station
    window_label = "{0} to {1} UTC".format(
        window_start.strftime("%Y-%m-%d %H:%M"),
        window_end.strftime("%Y-%m-%d %H:%M"),
    )
    file_stamp = "{0}_{1}".format(
        window_start.strftime("%Y%m%dT%H%M"),
        window_end.strftime("%Y%m%dT%H%M"),
    )
    ref_uref_label = "{0:g} µPa".format(reference_pressure_pa * 1e6)

    fig.suptitle(
        "{0} Spectrograms\nPSD referenced to {1}, {2}".format(
            array_label, ref_uref_label, window_label
        ),
        fontsize=15,
        fontweight="semibold",
        y=0.965,
    )
    axes[-1].set_xlabel("Time (UTC)")

    locator = mdates.AutoDateLocator()
    formatter = mdates.DateFormatter("%H:%M")
    axes[-1].xaxis.set_major_locator(locator)
    axes[-1].xaxis.set_major_formatter(formatter)
    for ax in axes[:-1]:
        plt.setp(ax.get_xticklabels(), visible=False)

    if last_mesh is not None:
        cbar = fig.colorbar(last_mesh, cax=cax)
        cbar.set_label("PSD (dB re 20 µPa^2/Hz)")
        cbar.outline.set_linewidth(0.8)

    out_name = "combined_spectrogram_{0}.png".format(file_stamp)
    fig.savefig(str(output_dir / out_name), dpi=dpi)
    plt.close(fig)

    wave_fig, wave_axes = plt.subplots(
        len(waveform_panels),
        1,
        figsize=(15, 9.5),
        sharex=True,
        facecolor="white",
    )
    if not isinstance(wave_axes, np.ndarray):
        wave_axes = np.array([wave_axes])

    for ax, panel in zip(wave_axes, waveform_panels):
        trace = panel["trace"]
        ax.plot(panel["times"], panel["data_pa"], color="#0b3954", linewidth=0.8)
        ax.set_ylabel("{0}\nPa".format(trace.stats.station), rotation=0, labelpad=28)
        style_axes(ax)
        ax.text(
            0.015,
            0.9,
            "{0}".format(trace.id),
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            color="#202124",
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none", pad=2.5),
        )

    wave_fig.suptitle(
        "{0} Filtered Waveforms\n{1:g}-{2:g} Hz bandpass, {3}".format(
            array_label, filter_freqmin, filter_freqmax, window_label
        ),
        fontsize=15,
        fontweight="semibold",
        y=0.965,
    )
    wave_axes[-1].set_xlabel("Time (UTC)")
    wave_axes[-1].xaxis.set_major_locator(locator)
    wave_axes[-1].xaxis.set_major_formatter(formatter)
    for ax in wave_axes[:-1]:
        plt.setp(ax.get_xticklabels(), visible=False)

    wave_out = "filtered_waveforms_{0:g}_{1:g}Hz_{2}.png".format(
        filter_freqmin, filter_freqmax, file_stamp
    )
    wave_fig.savefig(str(output_dir / wave_out), dpi=dpi)
    plt.close(wave_fig)


if __name__ == "__main__":
    main()
