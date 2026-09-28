"""
gui_methods.py — v1.2.1 interactive method dialogs for the SeismoFK GUI.

Copyright (c) 2024-2026 Islam Hamama
Contact: islam.hamama@nriag.sci.eg
Licensed under the MIT License.

These Qt dialogs expose the new v1.2.1 array-analysis methods from the main
window without disturbing the conventional-FK run pipeline:

    • ArrayResponseDialog  — theoretical array response function (ARF).
    • SlownessMapDialog     — Capon / MUSIC / Bartlett slowness map for one window.
    • PMCCWindow            — Progressive Multi-Channel Correlation detector.

Each dialog receives a stream that already has `tr.stats.coordinates` attached
(see FKAnalysisGUI._array_stream_with_coords) and does its own detrend/filter,
so instrument-response removal is not required for direction estimation.
"""

import numpy as np
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QApplication, QCheckBox, QDialog, QVBoxLayout, QHBoxLayout, QLabel,
    QPushButton, QComboBox, QDoubleSpinBox, QSpinBox, QMessageBox, QWidget,
    QGridLayout, QFrame, QFileDialog,
)
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from mpl_toolkits.axes_grid1 import make_axes_locatable

import fk_analysis as fk
import pmcc as pmcc_mod
import result_plots
from result_ui import apply_result_style, result_header, style_plot_toolbar


class _MplDialog(QDialog):
    """Base dialog with a Matplotlib canvas, toolbar and a controls row."""

    def __init__(self, title, parent=None, figsize=(7.5, 6.5), subtitle="Array analysis result"):
        super().__init__(parent)
        apply_result_style(self)
        self.setWindowTitle(title)
        self.setMinimumSize(760, 600)
        self.resize(1050, 790)
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(16, 16, 16, 12)
        self._layout.setSpacing(10)
        self._layout.addWidget(result_header(title, subtitle))

        self.controls_frame = QFrame()
        self.controls_frame.setObjectName("resultCard")
        self.controls = QHBoxLayout(self.controls_frame)
        self.controls.setContentsMargins(14, 10, 14, 10)
        self.controls.setSpacing(10)
        self.controls_frame.hide()
        self._layout.addWidget(self.controls_frame)

        self.fig = Figure(figsize=figsize, facecolor="white",
                          constrained_layout=True)
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        style_plot_toolbar(self.toolbar)
        plot_card = QFrame()
        plot_card.setObjectName("resultCard")
        plot_layout = QVBoxLayout(plot_card)
        plot_layout.setContentsMargins(10, 7, 10, 10)
        plot_layout.setSpacing(3)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas, stretch=1)
        self._layout.addWidget(plot_card, stretch=1)

        info_card = QFrame()
        info_card.setObjectName("resultCard")
        info_layout = QVBoxLayout(info_card)
        info_layout.setContentsMargins(15, 10, 15, 10)
        self.info = QLabel("")
        self.info.setObjectName("resultInfo")
        self.info.setWordWrap(True)
        info_layout.addWidget(self.info)
        self._layout.addWidget(info_card)

        actions = QHBoxLayout()
        self.actions_row = actions
        actions.addStretch()
        export = QPushButton("Export figure")
        export.clicked.connect(self._export_figure)
        actions.addWidget(export)
        close = QPushButton("Close")
        close.clicked.connect(self.close)
        actions.addWidget(close)
        self._layout.addLayout(actions)

    def _export_figure(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export analysis figure", "analysis_figure.png",
            "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)")
        if not filename:
            return
        if not filename.lower().endswith((".png", ".pdf", ".svg")):
            filename += ".png"
        try:
            self.fig.savefig(filename, dpi=300, bbox_inches="tight")
            QMessageBox.information(self, "Exported", f"Figure saved:\n{filename}")
        except Exception as error:
            QMessageBox.critical(self, "Export error", str(error))

    def _map_colorbar(self, ax, artist, label):
        divider = make_axes_locatable(ax)
        color_axis = divider.append_axes("right", size="4%", pad=0.12)
        self.fig.colorbar(artist, cax=color_axis, label=label)


# ─────────────────────────────────────────────────────────────────────────────
#  Array Response Function
# ─────────────────────────────────────────────────────────────────────────────

class ArrayResponseDialog(_MplDialog):
    """Show the theoretical ARF for the array geometry over [fmin, fmax]."""

    def __init__(self, stream, fmin, fmax, parent=None):
        super().__init__("Array response", parent,
                         subtitle="Array geometry and spatial resolution")
        self.stream = stream
        self.fmin, self.fmax = fmin, fmax
        try:
            self._draw()
        except Exception as e:
            QMessageBox.critical(self, "ARF error", str(e))

    def _draw(self):
        sx, sy, arf = fk.array_response(self.stream, self.fmin, self.fmax,
                                        smax=4.0, ds=0.05)
        ax = self.fig.add_subplot(111)
        m = ax.pcolormesh(sx, sy, arf, shading="auto", cmap="viridis")
        ax.contour(sx, sy, arf, levels=[0.5], colors="white",
                   linewidths=0.8, alpha=0.7)
        ax.plot(0, 0, "r+", ms=12, mew=2)
        ax.set_xlabel("Slowness East  (s/km)")
        ax.set_ylabel("Slowness North  (s/km)")
        ax.set_aspect("equal")
        self._map_colorbar(ax, m, "Normalised response")
        ax.set_title(f"ARF   {self.fmin:.2f}–{self.fmax:.2f} Hz   "
                     f"({len(self.stream)} sensors)")
        self.info.setText(
            "White contour = half-power (resolution) limit.  Secondary peaks "
            "away from the centre indicate spatial aliasing.")
        self.canvas.draw()


# ─────────────────────────────────────────────────────────────────────────────
#  High-resolution slowness map (Capon / MUSIC / Bartlett)
# ─────────────────────────────────────────────────────────────────────────────

class SlownessMapDialog(_MplDialog):
    """Adaptive-beamforming slowness map for a single analysis window."""

    _METHODS = [("Capon (MVDR)", "capon"),
                ("MUSIC", "music"),
                ("Bartlett (conventional)", "bartlett")]

    def __init__(self, stream, fmin, fmax, win_start, win_length,
                 expected_baz=None, parent=None):
        super().__init__("Slowness map", parent,
                         subtitle="Compare Bartlett, Capon, and MUSIC estimates")
        self.stream = stream
        self.fmin, self.fmax = fmin, fmax
        self.win_start, self.win_length = win_start, win_length
        self.expected_baz = expected_baz

        self.controls.addWidget(QLabel("Method:"))
        self.method_combo = QComboBox()
        for label, _ in self._METHODS:
            self.method_combo.addItem(label)
        self.controls.addWidget(self.method_combo)

        self.controls.addWidget(QLabel("Sources (MUSIC):"))
        self.nsrc_spin = QSpinBox()
        self.nsrc_spin.setRange(1, max(1, len(stream) - 1))
        self.nsrc_spin.setValue(1)
        self.controls.addWidget(self.nsrc_spin)

        recompute = QPushButton("↻ Recompute")
        recompute.clicked.connect(self._draw)
        self.controls.addWidget(recompute)
        self.controls.addStretch()
        self.controls_frame.show()

        self._draw()

    def _draw(self):
        method = self._METHODS[self.method_combo.currentIndex()][1]
        try:
            win = self.stream.slice(self.win_start,
                                    self.win_start + self.win_length)
            if len(win) < 2 or any(tr.stats.npts < 16 for tr in win):
                raise ValueError("Picked window has too few samples; "
                                 "pick a start time inside the data and/or "
                                 "increase the window length.")
            win = win.copy()
            win.detrend("demean")
            win.filter("bandpass", freqmin=self.fmin, freqmax=self.fmax,
                       corners=4, zerophase=True)
            sx, sy, P = fk.fk_beamform(
                win, self.fmin, self.fmax, method=method,
                smax=4.0, ds=0.05, n_sources=self.nsrc_spin.value())
            pk = fk.estimate_slowness_peak(sx, sy, P)
        except Exception as e:
            QMessageBox.critical(self, "Slowness-map error", str(e))
            return

        self.fig.clear()
        ax = self.fig.add_subplot(111)
        m = ax.pcolormesh(sx, sy, P, shading="auto", cmap="turbo")
        ax.plot(pk["sx"], pk["sy"], "w*", ms=15, mec="k", mew=0.8,
                label="Peak")

        # slowness circles for reference apparent velocities
        for v in (300, 340, 400):
            ang = np.linspace(0, 2 * np.pi, 200)
            s = 1000.0 / v
            ax.plot(s * np.cos(ang), s * np.sin(ang), color="white",
                    lw=0.5, alpha=0.5)
            ax.text(0, s, f"{v} m/s", color="white", fontsize=6,
                    ha="center", va="bottom", alpha=0.7)

        if self.expected_baz is not None:
            br = np.deg2rad(self.expected_baz)
            ax.plot([0, -np.sin(br) * 4], [0, -np.cos(br) * 4],
                    color="cyan", ls="--", lw=1.4,
                    label=f"Expected {int(self.expected_baz)}°")

        ax.set_xlabel("Slowness East  (s/km)")
        ax.set_ylabel("Slowness North  (s/km)")
        ax.set_aspect("equal")
        self._map_colorbar(ax, m, f"{method.capitalize()} power")
        ax.set_title(f"{method.upper()} slowness map   "
                     f"{self.fmin:.2f}–{self.fmax:.2f} Hz")
        ax.legend(loc="upper right", fontsize=7)
        self.info.setText(
            f"Peak →  back-azimuth = {pk['baz']:.1f}°,   "
            f"apparent velocity = {pk['app_vel']:.0f} m/s,   "
            f"sharpness (peak/mean) = {pk['coherence']:.1f}")
        self.canvas.draw()


# ─────────────────────────────────────────────────────────────────────────────
#  PMCC detector window
# ─────────────────────────────────────────────────────────────────────────────

class PMCCWindow(_MplDialog):
    """Run the PMCC detector and show waveform, correlation (including the
    rejected noise windows), direction and velocity over the full span."""

    def __init__(self, stream, fmin, fmax, win_start, duration,
                 window_sec, units="counts", parent=None):
        super().__init__("PMCC detector", parent, figsize=(10, 8.5),
                         subtitle="Signal against noise: waveform, correlation, "
                                  "direction and trace velocity")
        self.setMinimumSize(980, 780)
        self.resize(1180, 900)
        self.stream = stream
        self.fmin, self.fmax = fmin, fmax
        self.units = units
        self.win_start = win_start
        self.record_start = max(tr.stats.starttime for tr in stream)
        self.record_end = min(tr.stats.endtime for tr in stream)

        panel = QWidget()
        grid = QGridLayout(panel)
        grid.setHorizontalSpacing(16)
        grid.setVerticalSpacing(7)

        self.span_spin = QDoubleSpinBox()
        self.span_spin.setRange(10, 7 * 86400)
        self.span_spin.setDecimals(0)
        self.span_spin.setSingleStep(300)
        self.span_spin.setSuffix(" s")
        # At least 20 windows so an arrival is seen against its noise.
        self.span_spin.setValue(max(float(duration), 20 * float(window_sec)))
        self.span_spin.setToolTip(
            "Length analysed from the picked start time. Long spans show "
            "arrivals against the background noise.")
        self.whole_chk = QCheckBox("Whole record")
        self.whole_chk.setToolTip(
            f"Analyse {self.record_start.strftime('%H:%M:%S')}–"
            f"{self.record_end.strftime('%H:%M:%S')} UTC "
            f"({self.record_end - self.record_start:.0f} s).")
        self.whole_chk.toggled.connect(self.span_spin.setDisabled)
        self.bands_spin = QSpinBox()
        self.bands_spin.setRange(2, 20)
        self.bands_spin.setValue(6)
        self.win_spin = QDoubleSpinBox()
        self.win_spin.setRange(2, 600)
        self.win_spin.setValue(float(window_sec))
        self.win_spin.setSuffix(" s")
        self.corr_spin = QDoubleSpinBox()
        self.corr_spin.setRange(0.0, 0.99)
        self.corr_spin.setSingleStep(0.05)
        self.corr_spin.setValue(0.5)
        self.corr_spin.setToolTip("Minimum mean positive pairwise correlation (0–1).")
        self.closure_spin = QDoubleSpinBox()
        self.closure_spin.setRange(0.001, 0.5)
        self.closure_spin.setDecimals(3)
        self.closure_spin.setSingleStep(0.01)
        self.closure_spin.setValue(0.05)
        self.closure_spin.setSuffix(" s")
        self.closure_spin.setToolTip("Maximum RMS triplet delay mismatch; smaller is stricter.")
        self.min_vel_spin = QSpinBox()
        self.min_vel_spin.setRange(50, 5000)
        self.min_vel_spin.setValue(200)
        self.min_vel_spin.setToolTip("Reject apparent velocities below this value.")
        self.max_vel_spin = QSpinBox()
        self.max_vel_spin.setRange(50, 5000)
        self.max_vel_spin.setValue(600)
        self.max_vel_spin.setToolTip("Reject apparent velocities above this value.")
        self.family_px_spin = QSpinBox()
        self.family_px_spin.setRange(2, 200)
        self.family_px_spin.setValue(5)
        self.family_px_spin.setToolTip(
            "Minimum number of linked pixels for a family (one detection).")
        self.family_baz_spin = QDoubleSpinBox()
        self.family_baz_spin.setRange(1.0, 45.0)
        self.family_baz_spin.setValue(10.0)
        self.family_baz_spin.setSuffix("°")
        self.family_baz_spin.setToolTip(
            "Pixels adjacent in time and frequency join the same family when "
            "their back-azimuths differ by at most this much (velocity: 15%).")
        for col, (label, widget) in enumerate((
                ("Span", self.span_spin), ("", self.whole_chk),
                ("Bands", self.bands_spin), ("Window", self.win_spin),
                ("Family min pixels", self.family_px_spin))):
            grid.addWidget(QLabel(label), 0, col)
            grid.addWidget(widget, 1, col)
        for col, (label, widget) in enumerate((
                ("Min correlation", self.corr_spin),
                ("Max closure", self.closure_spin),
                ("Min velocity (m/s)", self.min_vel_spin),
                ("Max velocity (m/s)", self.max_vel_spin),
                ("Family Δ azimuth", self.family_baz_spin))):
            grid.addWidget(QLabel(label), 2, col)
            grid.addWidget(widget, 3, col)
        self.controls.addWidget(panel, stretch=1)

        recompute = QPushButton("Run PMCC")
        recompute.setObjectName("primaryAction")
        recompute.clicked.connect(self._draw)
        self.controls.addWidget(recompute)
        self.controls_frame.show()

        self.families = None
        self.export_families_btn = QPushButton("Export families…")
        self.export_families_btn.setEnabled(False)
        self.export_families_btn.clicked.connect(self._export_families)
        self.actions_row.insertWidget(1, self.export_families_btn)

        if len(stream) < 3:
            QMessageBox.warning(self, "PMCC", "PMCC needs at least 3 sensors.")
        else:
            self._draw()

    def _interval(self):
        """UTC [start, end] to analyse, clipped to the common record."""
        if self.whole_chk.isChecked():
            return self.record_start, self.record_end
        start = max(self.win_start, self.record_start)
        end = min(start + self.span_spin.value(), self.record_end)
        if end - start < self.win_spin.value():
            raise ValueError("The span from the picked start time is shorter "
                             "than one PMCC window. Pick an earlier start or "
                             "use the whole record.")
        return start, end

    def _draw(self):
        try:
            if self.min_vel_spin.value() >= self.max_vel_spin.value():
                raise ValueError("Minimum velocity must be below maximum velocity.")
            stime, etime = self._interval()
            bands = pmcc_mod.log_bands(self.fmin, self.fmax,
                                       self.bands_spin.value())
            stats = {}
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                # A quarter-window hop gives each arrival enough pixels to
                # form a family (window/2 often leaves only 3-4 per arrival).
                df = pmcc_mod.pmcc(
                    self.stream, bands, window_sec=self.win_spin.value(),
                    step_sec=self.win_spin.value() / 4.0,
                    corr_min=self.corr_spin.value(),
                    consistency_max=self.closure_spin.value(),
                    min_app_vel=self.min_vel_spin.value(),
                    max_app_vel=self.max_vel_spin.value(),
                    stime=stime, etime=etime, diagnostics=stats)
            finally:
                QApplication.restoreOverrideCursor()
        except Exception as e:
            QMessageBox.critical(self, "PMCC error", str(e))
            return

        df, self.families = pmcc_mod.pmcc_families(
            df, time_tol=self.win_spin.value(),
            baz_tol=self.family_baz_spin.value(),
            min_pixels=self.family_px_spin.value(),
            pixel_overlap=4.0)          # window / hop (hop = window / 4)
        self.export_families_btn.setEnabled(not self.families.empty)
        summary = result_plots.pmcc_figure(
            self.fig, df, stats.get("windows"), self.stream,
            fmin=self.fmin, fmax=self.fmax, t0=stime, t1=etime,
            corr_min=self.corr_spin.value(),
            min_vel=self.min_vel_spin.value(),
            max_vel=self.max_vel_spin.value(), units=self.units,
            families=self.families, window_sec=self.win_spin.value(),
            title=f"PMCC  ·  {stime.strftime('%Y-%m-%d %H:%M:%S')} – "
                  f"{etime.strftime('%H:%M:%S')} UTC  ({etime - stime:.0f} s)")
        self.canvas.draw()

        rejected = (f"{stats['low_correlation']} low correlation  ·  "
                    f"{stats['inconsistent']} inconsistent  ·  "
                    f"{stats['outside_velocity']} outside velocity range")
        if df.empty:
            self.info.setText(
                f"0 of {stats['candidates']} windows accepted  ·  {rejected}. "
                "Adjust the gate responsible for most rejections.")
            return
        in_family = int((df["family"] > 0).sum())
        if self.families.empty:
            family_text = (f"No families: no group of ≥ "
                           f"{self.family_px_spin.value()} linked pixels. "
                           "Lower Family min pixels or widen the gates.")
        else:
            family_text = "\n".join(
                f"F{f.family}  {f.t_start:%H:%M:%S}–{f.t_end:%H:%M:%S}   "
                f"{f.baz:.1f}° ± {f.baz_ci95:.1f}°   "
                f"{f.app_vel:.0f} ± {f.vel_ci95:.0f} m/s   {f.n_pixels} px"
                + (f"   ⚠ direction changes {f.mixed_dbaz:.1f}° within the "
                   "family: possibly two sources" if f.mixed else "")
                for f in self.families.itertuples())
        self.info.setText(
            f"{len(self.families)} famil{'y' if len(self.families) == 1 else 'ies'}"
            f"  ·  {in_family} of {len(df)} accepted pixels in families  ·  "
            f"{stats['candidates']} windows tested  ·  ± values are 95% "
            f"confidence intervals\n{family_text}\n"
            f"Rejected: {rejected}. Grey dots are rejected windows (noise); "
            "hollow circles are isolated pixels.")

    def _export_families(self):
        if self.families is None or self.families.empty:
            return
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export PMCC families", "pmcc_families.csv", "CSV (*.csv)")
        if not filename:
            return
        if not filename.lower().endswith(".csv"):
            filename += ".csv"
        try:
            self.families.to_csv(filename, index=False)
            QMessageBox.information(self, "Exported",
                                    f"{len(self.families)} families saved:\n{filename}")
        except Exception as error:
            QMessageBox.critical(self, "Export error", str(error))


# ─────────────────────────────────────────────────────────────────────────────
#  Noise levels (time-domain RMS, dB re 20 µPa)
# ─────────────────────────────────────────────────────────────────────────────

class NoiseLevelsDialog(_MplDialog):
    """RMS noise level of every sensor in consecutive windows over the whole
    record, with L90 / L50 / L10 / Leq statistics and sensor-offset checks
    (see :mod:`noise_levels`)."""

    _STAT_HEADERS = ["Sensor", "Windows", "L90", "L50", "L10", "Leq",
                     "Min", "Max", "Offset", "Offset MAD", "Status"]

    def __init__(self, stream, fmin, fmax, units="counts", parent=None):
        import noise_levels as nl
        from PyQt5.QtWidgets import QHeaderView, QTableWidget

        _, ref_label = nl.reference(units)
        super().__init__("Noise levels", parent, figsize=(11, 8),
                         subtitle=f"Time-domain RMS per window  ·  {ref_label}")
        self.setMinimumSize(1000, 820)
        self.resize(1200, 940)
        self.stream, self.units = stream, units
        self.levels = self.stats = None
        self.t0 = min(tr.stats.starttime for tr in stream)
        self.t1 = max(tr.stats.endtime for tr in stream)

        panel = QWidget()
        grid = QGridLayout(panel)
        grid.setHorizontalSpacing(16)
        grid.setVerticalSpacing(7)
        self.window_spin = QDoubleSpinBox()
        self.window_spin.setRange(1, 3600)
        self.window_spin.setDecimals(0)
        self.window_spin.setValue(1)
        self.window_unit = QComboBox()
        self.window_unit.addItems(["minutes", "seconds"])
        self.window_unit.setToolTip("RMS window: e.g. 1 minute, 5 minutes, 30 seconds.")
        self.fmin_spin = QDoubleSpinBox()
        self.fmin_spin.setRange(0.001, 50)
        self.fmin_spin.setDecimals(3)
        self.fmin_spin.setValue(fmin)
        self.fmax_spin = QDoubleSpinBox()
        self.fmax_spin.setRange(0.01, 100)
        self.fmax_spin.setValue(fmax)
        self.broadband_chk = QCheckBox("Broadband (no filter)")
        self.broadband_chk.toggled.connect(
            lambda on: (self.fmin_spin.setDisabled(on),
                        self.fmax_spin.setDisabled(on)))
        self.warn_spin = QDoubleSpinBox()
        self.warn_spin.setRange(0.5, 20)
        self.warn_spin.setValue(3.0)
        self.warn_spin.setSuffix(" dB")
        self.fault_spin = QDoubleSpinBox()
        self.fault_spin.setRange(1, 40)
        self.fault_spin.setValue(6.0)
        self.fault_spin.setSuffix(" dB")
        self.diurnal_chk = QCheckBox("Hour-of-day view")
        self.diurnal_chk.setToolTip(
            "Replace the distribution panel with median level per hour (UTC); "
            "most useful for records longer than a day.")
        self.diurnal_chk.setChecked(self.t1 - self.t0 >= 86400)
        for col, (label, widget) in enumerate((
                ("Window", self.window_spin), ("", self.window_unit),
                ("Band min (Hz)", self.fmin_spin),
                ("Band max (Hz)", self.fmax_spin),
                ("", self.broadband_chk))):
            grid.addWidget(QLabel(label), 0, col)
            grid.addWidget(widget, 1, col)
        for col, (label, widget) in enumerate((
                ("Sensor warning ±", self.warn_spin),
                ("Sensor fault ±", self.fault_spin),
                ("", self.diurnal_chk))):
            grid.addWidget(QLabel(label), 2, col)
            grid.addWidget(widget, 3, col)
        self.controls.addWidget(panel, stretch=1)
        run = QPushButton("Compute levels")
        run.setObjectName("primaryAction")
        run.clicked.connect(self._draw)
        self.controls.addWidget(run)
        self.controls_frame.show()

        self.table = QTableWidget(0, len(self._STAT_HEADERS))
        self.table.setHorizontalHeaderLabels(self._STAT_HEADERS)
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setMinimumHeight(110)
        # Plot and statistics share a draggable split; the plot gets most of
        # the height by default.
        from PyQt5.QtWidgets import QSplitter
        plot_card = self._layout.itemAt(2).widget()
        self._layout.removeWidget(plot_card)
        splitter = QSplitter(Qt.Vertical)
        splitter.setChildrenCollapsible(False)
        splitter.addWidget(plot_card)
        splitter.addWidget(self.table)
        splitter.setStretchFactor(0, 4)
        splitter.setStretchFactor(1, 1)
        splitter.setSizes([640, 170])
        self._layout.insertWidget(2, splitter, stretch=1)

        self.export_levels_btn = QPushButton("Export levels…")
        self.export_levels_btn.setEnabled(False)
        self.export_levels_btn.clicked.connect(self._export_levels)
        self.actions_row.insertWidget(1, self.export_levels_btn)
        self._draw()

    def _window_seconds(self):
        factor = 60.0 if self.window_unit.currentText() == "minutes" else 1.0
        return self.window_spin.value() * factor

    def _draw(self):
        import noise_levels as nl
        from PyQt5.QtGui import QColor
        from PyQt5.QtWidgets import QTableWidgetItem

        window = self._window_seconds()
        broadband = self.broadband_chk.isChecked()
        fmin = None if broadband else self.fmin_spin.value()
        fmax = None if broadband else self.fmax_spin.value()
        try:
            if fmin is not None and fmin >= fmax:
                raise ValueError("Band min must be below band max.")
            if window > self.t1 - self.t0:
                raise ValueError("The window is longer than the record.")
            if self.fault_spin.value() <= self.warn_spin.value():
                raise ValueError("The fault threshold must exceed the warning.")
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                self.levels = nl.noise_levels(self.stream, window, fmin, fmax,
                                              units=self.units)
                self.stats = nl.noise_statistics(
                    self.levels, warn_db=self.warn_spin.value(),
                    fault_db=self.fault_spin.value(), units=self.units)
            finally:
                QApplication.restoreOverrideCursor()
        except Exception as e:
            QMessageBox.critical(self, "Noise levels", str(e))
            return

        band = "broadband" if broadband else f"{fmin:g}–{fmax:g} Hz"
        span = (f"{window / 60:g} min" if window >= 60 else f"{window:g} s")
        result_plots.noise_figure(
            self.fig, self.levels, self.stats, units=self.units,
            band_label=band, t0=self.t0, t1=self.t1,
            warn_db=self.warn_spin.value(), fault_db=self.fault_spin.value(),
            diurnal=self.diurnal_chk.isChecked(),
            title=f"Noise levels  ·  {band}  ·  {span} windows")
        self.canvas.draw()

        self.table.setRowCount(len(self.stats))
        tint = {"fault": QColor("#fbeceb"), "warning": QColor("#fdf3e6")}
        for r, row in enumerate(self.stats.itertuples()):
            values = [row.station, f"{row.n_windows}", f"{row.L90:.1f}",
                      f"{row.L50:.1f}", f"{row.L10:.1f}", f"{row.Leq:.1f}",
                      f"{row.min_db:.1f}", f"{row.max_db:.1f}",
                      f"{row.offset_db:+.2f}", f"{row.offset_mad_db:.2f}",
                      row.flag.upper()]
            for c, value in enumerate(values):
                item = QTableWidgetItem(value)
                item.setTextAlignment(Qt.AlignCenter)
                if row.flag in tint:
                    item.setBackground(tint[row.flag])
                self.table.setItem(r, c, item)
        self.export_levels_btn.setEnabled(not self.levels.empty)

        _, ref_label = nl.reference(self.units)
        n_bad = int((self.stats["flag"] != "ok").sum())
        array_l50 = float(self.stats["L50"].median()) if len(self.stats) else float("nan")
        note = ("" if self.units == "Pa" else
                "  ·  ⚠ Not calibrated to Pa: levels are relative (dB re 1 count).")
        self.info.setText(
            f"{len(self.levels)} windows  ·  array L50 {array_l50:.1f} {ref_label}"
            f"  ·  {n_bad} sensor(s) flagged  ·  levels in {ref_label}; "
            f"L90/L50/L10 are exceeded 90/50/10% of the time; Leq averages "
            f"power.{note}")

    def _export_levels(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export noise levels", "noise_levels.csv", "CSV (*.csv)")
        if not filename:
            return
        base = filename[:-4] if filename.lower().endswith(".csv") else filename
        try:
            self.levels.to_csv(f"{base}.csv", index=False)
            self.stats.to_csv(f"{base}_stats.csv", index=False)
            QMessageBox.information(
                self, "Exported",
                f"Levels: {base}.csv\nStatistics: {base}_stats.csv")
        except Exception as error:
            QMessageBox.critical(self, "Export error", str(error))
