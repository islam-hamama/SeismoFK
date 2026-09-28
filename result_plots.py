"""
result_plots.py — Qt-free result figures shared by the GUI and the batch CLI.

Copyright (c) 2024-2026 Islam Hamama
Contact: islam.hamama@nriag.sci.eg
Licensed under the MIT License.

Every function draws into a caller-supplied ``matplotlib.figure.Figure`` so the
GUI can reuse its embedded canvas and the CLI can save straight to disk with no
display.  Each figure spans the full analysed interval and starts with the
band-passed waveform, so detections can be judged against the surrounding
noise rather than seen in isolation.

Units: time axes are UTC; back-azimuth in degrees; trace velocity in m/s;
frequency in Hz.
"""

import matplotlib.dates as mdates
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm
from obspy import UTCDateTime

import fk_analysis as fk
from matplotlib.ticker import (FixedFormatter, FixedLocator, FuncFormatter,
                               NullFormatter, NullLocator)

import pmcc as pmcc_mod

_GRID = dict(color="#e3e9ee", linewidth=0.6)

# Typography and ink shared with the desktop theme (result_ui.PALETTE).
PLOT_STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Inter", "Helvetica Neue", "Helvetica", "Segoe UI",
                        "Arial", "DejaVu Sans"],
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.titleweight": "semibold",
    "axes.titlecolor": "#14263a",
    "axes.labelsize": 9,
    "axes.labelcolor": "#34495c",
    "axes.edgecolor": "#a9b7c2",
    "axes.linewidth": 0.8,
    "xtick.color": "#4b6072",
    "ytick.color": "#4b6072",
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "grid.color": "#e3e9ee",
    "legend.fontsize": 8,
    "legend.edgecolor": "#dde4ea",
    "legend.framealpha": 0.92,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "figure.titleweight": "semibold",
}


def apply_plot_style():
    """Apply the SeismoFK figure style (GUI and CLI figures match)."""
    import logging

    import matplotlib

    # Missing optional families (e.g. Inter) are expected; do not warn.
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    matplotlib.rcParams.update(PLOT_STYLE)
_NOISE = "#b7c3cc"
_LINE = "#1d506f"


def _utc(t):
    """UTCDateTime / datetime → matplotlib date number."""
    return mdates.date2num(getattr(t, "datetime", t))


def reference_waveform(st, fmin, fmax, t0, t1):
    """
    Band-passed waveform of the first trace over [t0, t1].

    Returns (trace_id, time_nums, data) or None if no data overlap.  A single
    sensor is used on purpose: it shows the signal-to-noise ratio a reader can
    check against the raw record, without the gain of an array beam.
    """
    if not st:
        return None
    tr = st[0].copy()
    tr.trim(t0, t1)
    if tr.stats.npts < 8:
        return None
    tr.detrend("demean")
    tr.taper(0.02)
    nyquist = tr.stats.sampling_rate / 2.0
    tr.filter("bandpass", freqmin=fmin, freqmax=min(fmax, 0.95 * nyquist),
              corners=4, zerophase=True)
    times = _utc(tr.stats.starttime) + tr.times() / 86400.0
    return tr.id, times, tr.data


def _waveform_panel(ax, st, fmin, fmax, t0, t1, units):
    ref = reference_waveform(st, fmin, fmax, t0, t1)
    if ref is None:
        ax.text(0.5, 0.5, "No waveform in this interval", ha="center",
                va="center", transform=ax.transAxes)
        return
    trace_id, times, data = ref
    ax.plot(times, data, color="#34495e", lw=0.5)
    ax.set_ylabel(units, fontsize=8)
    ax.text(0.005, 0.05, f"{trace_id}  ·  {fmin:g}–{fmax:g} Hz",
            transform=ax.transAxes, fontsize=8, va="bottom",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none",
                      alpha=0.85))


def _time_axis(ax, t0, t1):
    ax.set_xlim(_utc(t0), _utc(t1))
    span = float(t1 - t0)
    fmt = "%H:%M" if span > 6 * 3600 else "%H:%M:%S"
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mdates.DateFormatter(fmt))
    ax.set_xlabel(f"UTC time  ({getattr(t0, 'datetime', t0):%Y-%m-%d})")
    for label in ax.get_xticklabels():
        label.set_rotation(20)
        label.set_ha("right")


def _frequency_colorbar(fig, artist, axes, df):
    colorbar = fig.colorbar(artist, ax=axes, pad=0.015, fraction=0.03)
    colorbar.set_label("Center frequency (Hz)")
    centres = sorted(df["f_center"].unique())
    if len(centres) <= 10:
        colorbar.ax.yaxis.set_major_locator(FixedLocator(centres))
        colorbar.ax.yaxis.set_major_formatter(
            FixedFormatter([f"{f:.2g}" for f in centres]))
        colorbar.ax.yaxis.set_minor_locator(NullLocator())


def pmcc_figure(fig, df, windows, st, *, fmin, fmax, t0, t1, corr_min,
                min_vel, max_vel, units="counts", title=None, families=None,
                window_sec=0.0):
    """
    Four-panel PMCC display over the full interval [t0, t1]:

    1. band-passed reference waveform (signal vs noise context);
    2. mean correlation of *every* window and band — rejected windows in grey
       form the noise floor, accepted pixels are coloured by frequency;
    3. back-azimuth and 4. trace velocity of accepted pixels.

    *windows* is ``diagnostics['windows']`` from :func:`pmcc.pmcc` (may be
    None).  With *families* (from :func:`pmcc.pmcc_families`, and *df*
    carrying its ``family`` column) each family's time span is shaded and
    labelled, family pixels are coloured by frequency and isolated pixels
    are drawn faint; direction and velocity summaries then use family pixels.
    Returns a summary dict (n, n_families, direction, concentration,
    median_vel).
    """
    fig.clear()
    axes = fig.subplots(4, 1, sharex=True,
                        gridspec_kw=dict(height_ratios=[1.1, 1, 1.25, 1]))
    ax_w, ax_c, ax_b, ax_v = axes
    _waveform_panel(ax_w, st, fmin, fmax, t0, t1, units)
    ax_w.set_title(title or "PMCC detector", fontsize=11, fontweight="bold",
                   loc="left")

    if windows is not None and len(windows):
        rejected = windows[~windows["accepted"]]
        ax_c.scatter(mdates.date2num(rejected["time"]), rejected["mean_corr"],
                     s=6, color=_NOISE, lw=0, label="rejected (noise)",
                     zorder=1)
    ax_c.axhline(corr_min, color="crimson", ls="--", lw=1,
                 label=f"threshold {corr_min:.2f}")
    ax_c.set_ylim(0, 1.02)
    ax_c.set_ylabel("Mean corr.", fontsize=8)

    summary = dict(n=0 if df is None else len(df), n_families=0,
                   direction=np.nan, concentration=0.0, median_vel=np.nan)
    use_families = (families is not None and df is not None
                    and "family" in df)
    if use_families and not df.empty:
        isolated = df[df["family"] == 0]
        if len(isolated):
            ti = mdates.date2num(isolated["time"])
            faint = dict(s=12, facecolors="none", edgecolors="#9fb0bd",
                         linewidths=0.6, zorder=2)
            ax_c.scatter(ti, isolated["mean_corr"], label="isolated pixel",
                         **faint)
            ax_b.scatter(ti, isolated["baz"], **faint)
            ax_v.scatter(ti, isolated["app_vel"], **faint)
        _family_spans(axes, families, window_sec, t0, t1)
        summary["n_families"] = len(families)
        df = df[df["family"] > 0]
    if df is None or df.empty:
        for ax in (ax_b, ax_v):
            ax.text(0.5, 0.5, "No PMCC detections in this interval",
                    ha="center", va="center", transform=ax.transAxes)
    else:
        # Larger low-frequency pixels drawn first keep every band visible
        # where several bands agree on the same time, direction and velocity.
        plot_df = df.sort_values("f_center")
        t = mdates.date2num(plot_df["time"])
        fraction = np.clip(np.log(plot_df["f_center"] / fmin)
                           / np.log(fmax / fmin), 0.0, 1.0)
        style = dict(c=plot_df["f_center"], cmap="viridis",
                     norm=LogNorm(vmin=fmin, vmax=fmax),
                     s=10 + 36 * (1 - fraction), edgecolors="#21354c",
                     linewidths=0.25, zorder=3)
        ax_c.scatter(t, plot_df["mean_corr"],
                     label="family pixel" if use_families else "accepted",
                     **style)
        sc = ax_b.scatter(t, plot_df["baz"], **style)
        ax_v.scatter(t, plot_df["app_vel"], **style)

        radians = np.deg2rad(df["baz"].to_numpy())
        weights = df["mean_corr"].to_numpy()
        summary.update(
            direction=pmcc_mod.circular_mean_deg(df["baz"], weights=weights),
            concentration=float(np.hypot(np.sum(weights * np.sin(radians)),
                                         np.sum(weights * np.cos(radians)))
                                / np.sum(weights)),
            median_vel=float(df["app_vel"].median()))
        if summary["concentration"] >= 0.5:
            ax_b.axhline(summary["direction"], color=_LINE, ls="--", lw=1,
                         alpha=0.75, label=f"mean {summary['direction']:.1f}°")
            ax_b.legend(loc="upper right", fontsize=7, framealpha=0.9)
        ax_v.axhline(summary["median_vel"], color=_LINE, ls="--", lw=1,
                     alpha=0.75, label=f"median {summary['median_vel']:.0f} m/s")
        ax_v.legend(loc="upper right", fontsize=7, framealpha=0.9)
        _frequency_colorbar(fig, sc, list(axes), df)
    legend = ax_c.legend(loc="upper right", fontsize=7, framealpha=0.9,
                         ncol=3)
    for handle in getattr(legend, "legend_handles",
                          getattr(legend, "legendHandles", [])):
        if hasattr(handle, "set_sizes"):
            handle.set_sizes([18])       # uniform key, whatever the data size

    ax_b.set_ylim(-5, 365)
    ax_b.set_yticks([0, 90, 180, 270, 360])
    ax_b.set_ylabel("Baz (°)", fontsize=8)
    ax_v.set_ylim(min_vel, max_vel)
    ax_v.set_ylabel("Vel. (m/s)", fontsize=8)
    for ax in axes:
        ax.set_axisbelow(True)
        ax.grid(**_GRID)
        ax.tick_params(labelsize=8)
    fig.align_ylabels(list(axes))
    _time_axis(ax_v, t0, t1)
    return summary


def _interval(lo, hi):
    """'[lo, hi]' with just enough decimals that the two ends differ, so a
    narrow interval is never displayed as zero width."""
    width = abs(hi - lo)
    decimals = 1 if width >= 1 else 2 if width >= 0.1 else 3
    if round(lo, decimals) == round(hi, decimals):
        return f"[{lo:.{decimals}f}, {hi:.{decimals}f}] (< {10 ** -decimals:g} wide)"
    return f"[{lo:.{decimals}f}, {hi:.{decimals}f}]"


def _family_spans(axes, families, window_sec, t0, t1):
    """Shade each PMCC family's time span (pixel centres ± half a window)
    on every panel, label it on the waveform panel, and mark its estimate
    with a 95 % confidence error bar on the back-azimuth and velocity
    panels.  Labels of families close in time are stacked on separate rows
    so they never overlap; a family flagged ``mixed`` gets a red label."""
    if families is None or families.empty:
        return
    half = pd.to_timedelta(window_sec / 2.0, unit="s")
    lo, hi = _utc(t0), _utc(t1)
    label_width = 0.16                  # axes fraction reserved per label
    row_end = []                        # right edge of the last label per row
    for fam in families.itertuples():
        x0 = mdates.date2num(fam.t_start - half)
        x1 = mdates.date2num(fam.t_end + half)
        for ax in axes:
            ax.axvspan(x0, x1, color="#0b7a75", alpha=0.09, lw=0, zorder=0)
        frac = (x0 - lo) / (hi - lo) if hi > lo else 0.0
        row = next((r for r, end in enumerate(row_end) if frac >= end),
                   len(row_end))
        if row == len(row_end):
            row_end.append(0.0)
        row_end[row] = frac + label_width
        ci_baz = getattr(fam, "baz_ci95", np.nan)
        ci_vel = getattr(fam, "vel_ci95", np.nan)
        mixed = bool(getattr(fam, "mixed", False))
        label = (f"F{fam.family}  {fam.baz:.0f}°" +
                 (f"±{ci_baz:.1f}" if np.isfinite(ci_baz) else "") +
                 f"  {fam.app_vel:.0f} m/s" + ("  ⚠ mixed" if mixed else ""))
        ink, edge = ("#b3261e", "#e6bdb9") if mixed else ("#07524f", "#b9dcd9")
        axes[0].annotate(
            label, xy=(x0, 1.0), xycoords=("data", "axes fraction"),
            xytext=(2, -3 - 15 * row), textcoords="offset points", va="top",
            fontsize=7.5, fontweight="semibold", color=ink,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=edge,
                      lw=0.6, alpha=0.95), zorder=6)
        # Family estimate with its 95 % interval on the baz / velocity panels.
        if len(axes) >= 4 and np.isfinite(ci_baz):
            mid = (x0 + x1) / 2.0
            bar = dict(fmt="D", ms=5, color=ink, mfc="white", mew=1.2,
                       elinewidth=1.4, capsize=4, zorder=7)
            axes[2].errorbar(mid, fam.baz, yerr=ci_baz, **bar)
            axes[3].errorbar(mid, fam.app_vel, yerr=ci_vel, **bar)


def fk_scan_figure(fig, df, st, *, fmin, fmax, t0, t1, method,
                   units="counts", coherence_min=None, min_vel=200.0,
                   max_vel=600.0, title=None):
    """
    Four-panel display of an FK scan (:func:`fk_analysis.fk_scan`) over
    [t0, t1]: waveform, peak-to-mean slowness power (a sharpness measure,
    high where one coherent arrival dominates), back-azimuth and trace
    velocity — the latter two coloured by that sharpness.  The velocity
    axis is fixed to [min_vel, max_vel]: a noise window whose slowness peak
    sits at the grid origin has a near-infinite velocity and would otherwise
    flatten the panel.
    """
    fig.clear()
    axes = fig.subplots(4, 1, sharex=True,
                        gridspec_kw=dict(height_ratios=[1.1, 1, 1.25, 1]))
    ax_w, ax_c, ax_b, ax_v = axes
    _waveform_panel(ax_w, st, fmin, fmax, t0, t1, units)
    ax_w.set_title(title or f"{method.upper()} FK scan", fontsize=11,
                   fontweight="bold", loc="left")
    if df is None or df.empty:
        for ax in (ax_c, ax_b, ax_v):
            ax.text(0.5, 0.5, "No FK windows in this interval", ha="center",
                    va="center", transform=ax.transAxes)
    else:
        t = mdates.date2num(df["time"])
        coherence = df["coherence"].to_numpy(dtype=float)
        ax_c.plot(t, coherence, color=_LINE, lw=0.8)
        if coherence_min:
            ax_c.axhline(coherence_min, color="crimson", ls="--", lw=1,
                         label=f"threshold {coherence_min:g}")
            ax_c.legend(loc="upper right", fontsize=7)
        order = np.argsort(coherence)            # sharpest peaks on top
        style = dict(c=coherence[order], cmap="plasma", s=16,
                     edgecolors="#21354c", linewidths=0.3, zorder=3)
        sc = ax_b.scatter(t[order], df["baz"].to_numpy()[order], **style)
        ax_v.scatter(t[order], df["app_vel"].to_numpy()[order], **style)
        colorbar = fig.colorbar(sc, ax=list(axes), pad=0.015, fraction=0.03)
        colorbar.set_label("Peak / mean power")
    ax_c.set_ylabel("Peak/mean", fontsize=8)
    ax_b.set_ylim(-5, 365)
    ax_b.set_yticks([0, 90, 180, 270, 360])
    ax_b.set_ylabel("Baz (°)", fontsize=8)
    ax_v.set_ylim(min_vel, max_vel)
    ax_v.set_ylabel("Vel. (m/s)", fontsize=8)
    for ax in axes:
        ax.set_axisbelow(True)
        ax.grid(**_GRID)
        ax.tick_params(labelsize=8)
    fig.align_ylabels(list(axes))
    _time_axis(ax_v, t0, t1)


def run_summary_figure(fig, df, *, method, t0, t1, fmin, fmax,
                       min_vel=None, max_vel=None, title=None, families=None):
    """
    Whole-run overview for the CLI: back-azimuth and trace velocity of every
    archived detection between t0 and t1, coloured by frequency (PMCC) or by
    peak-to-mean power (FK methods).
    """
    fig.clear()
    ax_b, ax_v = fig.subplots(2, 1, sharex=True)
    ax_b.set_title(title or f"{method.upper()} detections", fontsize=11,
                   fontweight="bold", loc="left")
    if df is None or df.empty:
        for ax in (ax_b, ax_v):
            ax.text(0.5, 0.5, "No detections", ha="center", va="center",
                    transform=ax.transAxes)
    else:
        t = mdates.date2num(df["time"])
        if method == "pmcc":
            style = dict(c=df["f_center"], cmap="viridis",
                         norm=LogNorm(vmin=fmin, vmax=fmax), s=10)
            label = "Center frequency (Hz)"
        elif method == "fk":
            style = dict(c=df["semblance"], cmap="plasma", vmin=0, vmax=1,
                         s=8)
            label = "Semblance"
        else:
            style = dict(c=df["coherence"], cmap="plasma", s=8)
            label = "Peak / mean power"
        style.update(edgecolors="none", zorder=3)
        sc = ax_b.scatter(t, df["baz"], **style)
        ax_v.scatter(t, df["app_vel"], **style)
        colorbar = fig.colorbar(sc, ax=[ax_b, ax_v], pad=0.015,
                                fraction=0.03, label=label)
        if method == "pmcc":             # plain Hz ticks on the log scale
            colorbar.ax.yaxis.set_major_formatter(
                FuncFormatter(lambda v, _: f"{v:g}"))
            colorbar.ax.yaxis.set_minor_formatter(NullFormatter())
    if families is not None and not families.empty:
        # One outlined marker per PMCC family at its start time.
        tf = mdates.date2num(families["t_start"])
        ring = dict(s=70, facecolors="none", edgecolors="#b3261e",
                    linewidths=1.2, zorder=5)
        ax_b.scatter(tf, families["baz"], label=f"{len(families)} families",
                     **ring)
        ax_v.scatter(tf, families["app_vel"], **ring)
        ax_b.legend(loc="upper right", fontsize=8)
    ax_b.set_ylim(-5, 365)
    ax_b.set_yticks([0, 90, 180, 270, 360])
    ax_b.set_ylabel("Back-azimuth (°)")
    ax_v.set_ylabel("Trace velocity (m/s)")
    if min_vel is not None and max_vel is not None:
        ax_v.set_ylim(min_vel, max_vel)
    for ax in (ax_b, ax_v):
        ax.set_axisbelow(True)
        ax.grid(**_GRID)
    _time_axis(ax_v, t0, t1)


def fk_results_figure(fig, r, params, *, event_name, fmin, fmax, semb_thresh):
    """
    The six-panel conventional-FK result figure shared by the GUI results
    window and ``seismofk_cli --method fk --save-figures``.

    Left column: Fisher statistic, back-azimuth, trace velocity and the beam
    over time (noise windows grey, detections coloured by semblance).  Right
    column: polar back-azimuth/velocity detection map and array geometry.

    *r* is the dict from :func:`fk_analysis.fk_array` (``beam_units`` is used
    for the beam axis).  *params* may hold ``bootstrap`` (bool),
    ``origin_time`` + ``celerity`` (expected-arrival line) and
    ``event_lat``/``event_lon``.  The expected direction is drawn only when
    ``r['expected_bazi']`` is finite.  Returns ``{'n_det', 'boot'}``.
    """
    n_sta       = int(r.get('n_stations', 0))
    exp_baz     = r.get('expected_bazi', np.nan)
    has_exp     = exp_baz is not None and np.isfinite(exp_baz)


    gs = gridspec.GridSpec(
        4, 2,
        width_ratios=[2.5, 1.2],
        hspace=0.45, wspace=0.34,
        left=0.08, right=0.97,
        top=0.80,  bottom=0.14,
    )

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    ax3 = fig.add_subplot(gs[2, 0], sharex=ax1)
    ax4 = fig.add_subplot(gs[3, 0], sharex=ax1)
    ax5 = fig.add_subplot(gs[0:2, 1], projection='polar')  # BAZ detection map
    ax6 = fig.add_subplot(gs[2:4, 1])   # Array geometry

    semb     = r['semblance']
    det_mask = semb >= semb_thresh
    noi_mask = ~det_mask
    n_det    = int(np.sum(det_mask))

    # ── v1.2.1 — optional bootstrap uncertainty on baz / app-vel ───────
    boot = None
    if params.get('bootstrap') and n_det >= 2:
        try:
            # fk_array's ``overlap`` is the step fraction: 0.1 → windows
            # overlap 90 %, so resample blocks of 1/step consecutive windows.
            step = float(params.get('overlap') or 1.0)
            boot = fk.bootstrap_beam(
                np.asarray(r['bazi']), np.asarray(r['app_vel']),
                np.asarray(r['slowness']), mask=det_mask, n_boot=1000,
                block=max(1, round(1.0 / max(step, 1e-3))),
                slowness_step=float(params.get('slowness_step', 0.16)))
        except Exception as e:
            print(f"[WARN] bootstrap failed: {e}")

    # ── Visual encoding ────────────────────────────────────────────────
    # Noise:      tiny, light-grey, very transparent → recedes to background
    # Detection:  larger (scaled by semblance), plasma colormap, black edge
    det_sizes = 6 + 20 * ((semb[det_mask] - semb_thresh) /
                          max(1 - semb_thresh, 1e-6))

    kw_noise = dict(color='#cccccc', alpha=0.25, s=4,
                    linewidths=0, zorder=1)
    kw_det   = dict(c=semb[det_mask], cmap='plasma',
                    vmin=semb_thresh, vmax=1.0,
                    s=det_sizes, alpha=0.95,
                    edgecolors='black', linewidths=0.5, zorder=3)

    # ── Panel 1 — Fisher ──────────────────────────────────────────────
    ax1.scatter(r['time'][noi_mask], r['fisher'][noi_mask], **kw_noise)
    sc = ax1.scatter(r['time'][det_mask], r['fisher'][det_mask], **kw_det)
    # Threshold line
    fisher_thresh = (n_sta - 1) * semb_thresh / (1 - semb_thresh + 1e-9)
    ax1.axhline(fisher_thresh, color='crimson', ls='--', lw=1.2,
                label=f"Threshold (semb={semb_thresh:.2f})")
    ax1.set_ylabel('Fisher', fontsize=8)
    ax1.legend(loc='upper right', fontsize=7)
    ax1.set_facecolor('#fafafa')
    ax1.grid(True, alpha=0.3)

    # ── Panel 2 — Back-azimuth ─────────────────────────────────────────
    ax2.scatter(r['time'][noi_mask], r['bazi'][noi_mask], **kw_noise)
    ax2.scatter(r['time'][det_mask], r['bazi'][det_mask], **kw_det)
    if has_exp:
        ax2.axhline(exp_baz, color='royalblue', ls='--', lw=1.6,
                    label=f"Expected  {int(exp_baz)}°", zorder=4)
        # ±20° acceptance band
        baz_lo = (exp_baz - 20) % 360
        baz_hi = (exp_baz + 20) % 360
        if baz_lo < baz_hi:
            ax2.axhspan(baz_lo, baz_hi, color='royalblue', alpha=0.08, zorder=0)
    ax2.set_ylim(0, 360); ax2.set_yticks([0, 90, 180, 270, 360])
    ax2.set_ylabel('Azimuth (°)', fontsize=8)
    if has_exp:
        ax2.legend(loc='upper right', fontsize=7)
    ax2.set_facecolor('#fafafa')
    ax2.grid(True, alpha=0.3)

    # ── Panel 3 — Apparent velocity ────────────────────────────────────
    ax3.scatter(r['time'][noi_mask], r['app_vel'][noi_mask], **kw_noise)
    ax3.scatter(r['time'][det_mask], r['app_vel'][det_mask], **kw_det)
    ax3.axhspan(300, 380, color='limegreen', alpha=0.08, zorder=0,
                label='300–380 m/s')
    ax3.set_ylabel('Velocity (m/s)', fontsize=8)
    ax3.set_ylim(200, 450)
    ax3.legend(loc='upper right', fontsize=7)
    ax3.set_facecolor('#fafafa')
    ax3.grid(True, alpha=0.3)

    # ── Panel 4 — Beam waveform ────────────────────────────────────────
    bwave   = r['beam_waveform']
    btim    = r['beam_times']
    max_val = float(np.max(np.abs(bwave))) if len(bwave) > 0 else 1.0

    ax4.plot(btim, bwave, color='#333333', lw=0.8, zorder=2)
    ax4.fill_between(btim, 0, bwave, where=bwave >= 0,
                     color='steelblue', alpha=0.5, zorder=1)
    ax4.fill_between(btim, 0, bwave, where=bwave <  0,
                     color='tomato',    alpha=0.5, zorder=1)
    beam_units = r.get('beam_units', 'unverified units')
    ax4.set_ylabel(f'Beam ({beam_units})', fontsize=8)
    ax4.set_facecolor('#fafafa')
    ax4.grid(True, alpha=0.3)
    ax4.text(0.02, 0.94, f"Max |beam| = {max_val:.4g} {beam_units}",
             transform=ax4.transAxes, fontsize=8,
             bbox=dict(boxstyle='round,pad=0.3',
                       facecolor='white', edgecolor='grey', alpha=0.85))

    # ── Celerity arrival line (optional) ──────────────────────────────
    origin_time_str = params.get('origin_time')
    celerity_ms     = params.get('celerity')
    if origin_time_str and celerity_ms:
        try:
            from obspy.geodetics import gps2dist_azimuth
            arr_lat  = r.get('array_lat', 0)
            arr_lon  = r.get('array_lon', 0)
            src_lat  = params.get('event_lat', 0)
            src_lon  = params.get('event_lon', 0)

            dist_m    = gps2dist_azimuth(src_lat, src_lon,
                                         arr_lat,  arr_lon)[0]
            travel_s  = dist_m / celerity_ms
            arrival   = UTCDateTime(origin_time_str) + travel_s
            arr_mpl   = mdates.date2num(arrival.datetime)

            ax4.axvline(arr_mpl, color='darkorange', lw=2.0,
                        ls='--', zorder=5,
                        label=(f'Expected arrival\n'
                               f'c = {celerity_ms:.0f} m/s\n'
                               f'Δ = {dist_m/1000:.1f} km   '
                               f'Δt = {travel_s:.0f} s'))
            ax4.legend(loc='upper right', fontsize=7,
                       framealpha=0.9, edgecolor='darkorange')
        except Exception as _cel_err:
            print(f"[WARN] Could not plot celerity line: {_cel_err}")

    # ── Shared time axis formatting ────────────────────────────────────
    fmt = mdates.DateFormatter('%H:%M:%S')
    for ax in (ax1, ax2, ax3):
        ax.tick_params(labelbottom=False)
    for ax in (ax1, ax2, ax3, ax4):
        ax.xaxis.set_major_locator(mdates.AutoDateLocator())
        ax.xaxis.set_major_formatter(fmt)

    fig.autofmt_xdate(rotation=25, ha='right')

    # ── Panel 5 — Back-Azimuth Detection Map (polar) ──────────────────
    ax5.set_theta_zero_location('N')   # North at top
    ax5.set_theta_direction(-1)         # Clockwise like a compass

    baz_det  = r['bazi'][det_mask]
    vel_det  = r['app_vel'][det_mask]
    semb_det = semb[det_mask]

    # Rose histogram bars — semblance-weighted count per 10° bin
    n_bins     = 36
    bin_edges  = np.linspace(0, 360, n_bins + 1)
    if len(baz_det) > 0:
        counts, _ = np.histogram(baz_det, bins=bin_edges, weights=semb_det)
        counts_n  = counts / max(counts.max(), 1e-9)
    else:
        counts_n  = np.zeros(n_bins)

    bin_centers = np.deg2rad(0.5 * (bin_edges[:-1] + bin_edges[1:]))
    bin_width   = np.deg2rad(360 / n_bins)
    bars = ax5.bar(bin_centers, counts_n,
                   width=bin_width * 0.9, bottom=0,
                   color='steelblue', alpha=0.22,
                   edgecolor='steelblue', linewidth=0.4, zorder=1)

    # Detection scatter: theta=baz, r=app_vel, colour=semblance
    if len(baz_det) > 0:
        sc5 = ax5.scatter(
            np.deg2rad(baz_det), vel_det,
            c=semb_det, cmap='plasma',
            vmin=semb_thresh, vmax=1.0,
            s=6 + 20 * (semb_det - semb_thresh) /
              max(1 - semb_thresh, 1e-6),
            alpha=0.80, edgecolors='black',
            linewidths=0.5, zorder=3)
    else:
        sc5 = ax5.scatter([], [], c=[], cmap='plasma',
                          vmin=semb_thresh, vmax=1.0)

    # Expected BAZ line
    if has_exp:
        exp_rad = np.deg2rad(exp_baz)
        ax5.plot([exp_rad, exp_rad], [200, 450],
                 color='royalblue', lw=2.0, ls='--',
                 label=f'Expected {int(exp_baz)}°', zorder=4)

    # ── v1.2.1 bootstrap uncertainty overlay (green cross-hair) ─────────
    if boot is not None:
        bm, vm = boot['baz_mean'], boot['vel_mean']
        b_lo, b_hi = boot['baz_ci']
        v_lo, v_hi = boot['vel_ci']
        vm_c = min(max(vm, 200), 450)
        # Angular arc across the 95 % interval (unwrapped around the mean).
        lo_off = (b_lo - bm + 180.0) % 360.0 - 180.0
        hi_off = (b_hi - bm + 180.0) % 360.0 - 180.0
        arc = np.deg2rad(bm + np.linspace(lo_off, hi_off, 40))
        ax5.plot(arc, np.full_like(arc, vm_c), color='limegreen',
                 lw=3.0, solid_capstyle='round', zorder=6)
        # Radial line across the velocity interval.
        ax5.plot([np.deg2rad(bm)] * 2,
                 [max(v_lo, 200), min(v_hi, 450)],
                 color='limegreen', lw=3.0, zorder=6)
        ax5.plot(np.deg2rad(bm), vm_c, marker='o', color='limegreen',
                 mec='black', mew=0.8, ms=9, zorder=7,
                 label=f"Bootstrap {bm:.1f}° {_interval(b_lo, b_hi)}")

    # Radial axis = apparent velocity
    ax5.set_ylim(200, 450)
    ax5.set_yticks([250, 300, 350, 400, 450])
    ax5.set_yticklabels(['250', '300', '350', '400', '450\nm/s'],
                        fontsize=6, color='#444')
    # Put the radial tick labels opposite the detection cluster (or the
    # expected direction when nothing is detected) so they never overlap.
    # Diagonals keep the stacked labels from running together along a
    # horizontal or vertical radius.
    if len(baz_det) > 0:
        det_rad = np.deg2rad(baz_det)
        label_dir = np.rad2deg(np.arctan2(np.sin(det_rad).mean(),
                                          np.cos(det_rad).mean()))
    else:
        label_dir = exp_baz if has_exp else 0.0
    ax5.set_rlabel_position(max(
        (45, 135, 225, 315),
        key=lambda d: abs((d - label_dir + 180.0) % 360.0 - 180.0)))

    # Cardinal labels
    ax5.set_thetagrids(range(0, 360, 45),
                       ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'],
                       fontsize=8)
    ax5.set_title('Back-Azimuth Detection Map',
                  fontsize=10, fontweight='bold', pad=12)
    # Key for the expected-azimuth line and bootstrap marker, kept in the
    # gap under the polar axes so it cannot cover detections.
    if ax5.get_legend_handles_labels()[0]:
        ax5.legend(loc='upper center', bbox_to_anchor=(0.5, -0.07),
                   ncol=2, fontsize=7, frameon=False)
    ax5.grid(True, alpha=0.3)

    # Panel 6 — Array Geometry ─────────────────────────────────────────
    x_km    = r.get('station_x_km', np.array([]))
    y_km    = r.get('station_y_km', np.array([]))
    sta_ids = r.get('station_ids',  [])

    ax6.scatter(x_km, y_km, s=60, color='steelblue',
                edgecolors='k', linewidths=0.7, zorder=3)
    # Push each label away from its nearest neighbour so tightly spaced
    # sensors (e.g. an inner triangle) do not print over each other.
    pts = np.column_stack([x_km, y_km]) if len(x_km) else np.empty((0, 2))
    for i, sid in enumerate(sta_ids):
        others = np.delete(pts, i, axis=0)
        away = pts[i] - (others[np.argmin(np.hypot(*(others - pts[i]).T))]
                         if len(others) else 0.0)
        if not np.any(away):
            away = pts[i] if np.any(pts[i]) else np.array([1.0, 1.0])
        ux, uy = away / np.hypot(*away)
        ax6.annotate(sid, pts[i], textcoords='offset points',
                     xytext=(8 * ux, 8 * uy),
                     ha='left' if ux > 0.3 else 'right' if ux < -0.3 else 'center',
                     va='bottom' if uy > 0.3 else 'top' if uy < -0.3 else 'center',
                     fontsize=7, color='#2c3e50')
    ax6.plot(0, 0, 'k+', ms=12, mew=2, zorder=4, label='Centroid')
    ax6.set_xlabel('East offset  (km)', fontsize=9)
    ax6.set_ylabel('North offset  (km)', fontsize=9)
    ax6.set_title('Array Geometry', fontsize=10, fontweight='bold')
    ax6.set_aspect('equal', adjustable='datalim')
    # Pad the data limits (which 'datalim' aspect works from, unlike
    # set_xlim/set_ylim) to leave room for the outward labels.
    if len(pts):
        reach = max(np.abs(pts).max() * 1.35, 0.05)
        ax6.update_datalim([(-reach, -reach), (reach, reach)])
        ax6.autoscale_view()
    ax6.legend(fontsize=7, loc='upper right')
    ax6.grid(True, alpha=0.25, ls='--')

    # The polar tick labels and lower map title need their own space at
    # desktop window heights. Keep the data limits while shortening the
    # occupied boxes on the right column.
    polar_box = ax5.get_position()
    ax5.set_position([polar_box.x0, polar_box.y0 + 0.035,
                      polar_box.width, polar_box.height - 0.035])
    geometry_box = ax6.get_position()
    ax6.set_position([geometry_box.x0, geometry_box.y0,
                      geometry_box.width, geometry_box.height - 0.035])

    # Suptitle ─────────────────────────────────────────────────────────
    suptitle = (
        f"Event: {event_name}    "
        f"Filter: {fmin}–{fmax} Hz    "
        f"Sensors: {n_sta}    "
        + (f"Expected BAZ: {int(exp_baz)}°    " if has_exp else "")
        + f"Detections: {n_det}/{len(semb)}")
    if boot is not None:
        suptitle += (f"\nBootstrap {int(boot['ci'])}% CI "
                     f"({boot['n_used']} windows, blocks of {boot['block']}):"
                     f"   baz {boot['baz_mean']:.1f}° "
                     f"{_interval(*boot['baz_ci'])}"
                     f"    app. vel {boot['vel_mean']:.0f} m/s "
                     f"{_interval(*boot['vel_ci'])}"
                     + ("  (limited by slowness grid)"
                        if boot.get('resolution_limited') else ""))
    fig.suptitle(suptitle, fontsize=11, y=0.975, fontweight='bold')

    # Semblance colorbar — span only the LEFT (time-series) column and sit
    # below the suptitle, so it no longer overruns the polar-plot title.
    cbar_ax = fig.add_axes([0.08, 0.885, 0.555, 0.012])
    cb = fig.colorbar(sc, cax=cbar_ax, orientation='horizontal', extend='min')
    cb.set_label(f'Semblance  (detections ≥ {semb_thresh:.2f}  |  grey = noise)',
                 fontsize=9, fontweight='bold', labelpad=-1)
    cb.ax.tick_params(labelsize=8)

    return {'n_det': n_det, 'boot': boot}


_SENSOR_COLOURS = ["#3b6e8f", "#c77c2e", "#4f8f5b", "#8a5a9e", "#b8474f",
                   "#2f8f8a", "#7a6a3a", "#5a6f9e", "#9e5a7a", "#6b8f2f"]


def noise_figure(fig, levels, stats, *, units, band_label, t0, t1,
                 warn_db=3.0, fault_db=6.0, title=None, diurnal=False):
    """
    Noise-level report for :mod:`noise_levels` results.

    Top: level of every sensor per window and the array median over [t0, t1]
    (long runs are summarised as medians per time bin so millions of windows
    stay readable).  Bottom left: per-sensor L90–L50–L10 box with Leq
    (diamond), or with ``diurnal=True`` the median level per hour of day.
    Bottom right: each sensor's median offset from the array median with
    the ±warn / ±fault bands.
    """
    import noise_levels as nl

    _, ref_label = nl.reference(units)
    fig.clear()
    grid = fig.add_gridspec(2, 2, height_ratios=[1.35, 1])
    ax_t = fig.add_subplot(grid[0, :])
    ax_d = fig.add_subplot(grid[1, 0])
    ax_o = fig.add_subplot(grid[1, 1])
    ax_t.set_title(title or f"Noise levels  ·  {band_label}", loc="left",
                   fontsize=11, fontweight="semibold")
    valid = levels.dropna(subset=["level_db"])
    if valid.empty:
        for ax in (ax_t, ax_d, ax_o):
            ax.text(0.5, 0.5, "No complete windows", ha="center",
                    va="center", transform=ax.transAxes)
        return
    traces = sorted(valid["trace_id"].unique())
    colour = {tid: _SENSOR_COLOURS[i % len(_SENSOR_COLOURS)]
              for i, tid in enumerate(traces)}
    label = {tid: tid.split(".")[1] or tid for tid in traces}

    # Time series (binned for long runs).
    grid_lv = valid.pivot_table(index="time", columns="trace_id",
                                values="level_db")
    if len(grid_lv) > 3000:
        span = grid_lv.index[-1] - grid_lv.index[0]
        grid_lv = grid_lv.resample(max(span / 1500, pd.Timedelta("1min"))) \
            .median().dropna(how="all")
    for tid in traces:
        if tid in grid_lv:
            ax_t.plot(grid_lv.index, grid_lv[tid], lw=0.7, alpha=0.75,
                      color=colour[tid], label=label[tid])
    ax_t.plot(grid_lv.index, grid_lv.median(axis=1), lw=1.8,
              color="#14263a", label="array median")
    ax_t.set_ylabel(f"Level ({ref_label})")
    ax_t.legend(loc="upper right", fontsize=7, ncol=min(len(traces) + 1, 8))
    _time_axis(ax_t, t0, t1)

    stats = stats.set_index("trace_id").reindex(traces)
    if diurnal:
        hourly = nl.hourly_levels(valid)
        for tid in traces:
            if tid in hourly:
                ax_d.plot(hourly.index, hourly[tid], marker="o", ms=3,
                          lw=1, color=colour[tid])
        ax_d.plot(hourly.index, hourly.median(axis=1), lw=2,
                  color="#14263a")
        ax_d.set_xlim(-0.5, 23.5)
        ax_d.set_xticks(range(0, 24, 3))
        ax_d.set_xlabel("Hour of day (UTC)")
        ax_d.set_title("Median level by hour of day", loc="left",
                       fontsize=10)
    else:
        data = [valid.loc[valid["trace_id"] == tid, "level_db"].to_numpy()
                for tid in traces]
        box = ax_d.boxplot(data, whis=(10, 90), showfliers=False,
                           patch_artist=True, widths=0.55,
                           medianprops=dict(color="#14263a", lw=1.5))
        for patch, tid in zip(box["boxes"], traces):
            patch.set(facecolor=colour[tid], alpha=0.35, edgecolor=colour[tid])
        ax_d.scatter(range(1, len(traces) + 1), stats["Leq"], marker="D",
                     s=26, color="#b3261e", zorder=4, label="Leq")
        ax_d.set_xticks(range(1, len(traces) + 1))
        ax_d.set_xticklabels([label[t] for t in traces], fontsize=8)
        ax_d.legend(loc="upper right", fontsize=7)
        ax_d.set_title("L90 · L50 · L10 (whiskers, box = quartiles)",
                       loc="left", fontsize=10)
    ax_d.set_ylabel(ref_label)

    offsets = stats["offset_db"].to_numpy(float)
    bar_colour = ["#b3261e" if f == "fault" else "#c77c2e" if f == "warning"
                  else "#4f8f5b" for f in stats["flag"]]
    ypos = np.arange(len(traces))
    lim = max(fault_db * 1.25, np.nanmax(np.abs(offsets)) * 1.15)
    ax_o.axvspan(-lim, -fault_db, color="#b3261e", alpha=0.08)
    ax_o.axvspan(fault_db, lim, color="#b3261e", alpha=0.08)
    ax_o.axvspan(-fault_db, -warn_db, color="#c77c2e", alpha=0.08)
    ax_o.axvspan(warn_db, fault_db, color="#c77c2e", alpha=0.08)
    ax_o.barh(ypos, offsets, xerr=stats["offset_mad_db"], color=bar_colour,
              height=0.55, error_kw=dict(ecolor="#4b6072", lw=0.8, capsize=2))
    ax_o.axvline(0, color="#14263a", lw=0.8)
    ax_o.set_xlim(-lim, lim)
    ax_o.set_yticks(ypos)
    ax_o.set_yticklabels([label[t] for t in traces], fontsize=8)
    ax_o.invert_yaxis()
    ax_o.set_xlabel("Offset from the other sensors' median (dB)")
    ax_o.set_title(f"Sensor offsets  ·  warn ±{warn_db:g} dB, "
                   f"fault ±{fault_db:g} dB", loc="left", fontsize=10)
    for ax in (ax_t, ax_d, ax_o):
        ax.set_axisbelow(True)
        ax.grid(**_GRID)
