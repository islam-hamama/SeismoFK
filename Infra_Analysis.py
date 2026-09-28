"""
Infra_Analysis.py — SeismoFK Infrasound FK Array Analysis GUI

Copyright (c) 2024-2026 Islam Hamama
Contact: islam.hamama@nriag.sci.eg

Licensed under the MIT License — see LICENSE for details.
"""
import sys
import os
import shutil
import json
import html

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.dates as mdates
from matplotlib.figure import Figure
import numpy as np

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QLineEdit, QFileDialog, QMessageBox,
    QDoubleSpinBox, QDateTimeEdit, QProgressBar, QDialog, QFrame,
    QComboBox, QStyledItemDelegate, QGroupBox, QSizePolicy, QScrollArea,
    QTextEdit, QTableWidget, QTableWidgetItem, QHeaderView, QCheckBox,
)
from PyQt5.QtCore  import Qt, QDateTime, QThread, pyqtSignal
from PyQt5.QtGui   import QFontDatabase, QPixmap, QStandardItem, QStandardItemModel

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)

from obspy import read, read_inventory, UTCDateTime
import fk_analysis as fk
from fk_analysis import fk_array
import db_manager
from app_paths import inventory_dir
import result_plots

result_plots.apply_plot_style()
from matplotlib.ticker import FuncFormatter
from result_ui import (PALETTE, _track, apply_app_theme, apply_main_style,
                       apply_result_style, metric_card, result_header,
                       style_plot_toolbar)


# ═══════════════════════════════════════════════════════════════════════════
#  Background worker thread
# ═══════════════════════════════════════════════════════════════════════════
def _envelope(data, delta, offset, max_points=6000):
    """
    Display envelope of a trace: (times_s, values), demeaned.  Long traces
    are reduced to per-bin min/max pairs, which keeps every spike visible
    while drawing a fixed number of points.
    """
    data = np.asarray(data, dtype=float)
    data = data - data.mean() if data.size else data
    n = data.size
    if n <= max_points:
        return offset + np.arange(n) * delta, data
    bins = max_points // 2
    step = n // bins
    blocks = data[:bins * step].reshape(bins, step)
    centres = offset + (np.arange(bins) * step + step / 2.0) * delta
    return (np.repeat(centres, 2),
            np.column_stack([blocks.min(axis=1), blocks.max(axis=1)]).ravel())


class ProcessThread(QThread):
    finished = pyqtSignal(dict)
    progress = pyqtSignal(int)
    error    = pyqtSignal(str)
    status   = pyqtSignal(str)

    def __init__(self, params):
        super().__init__()
        self.params = params

    def run(self):
        try:
            # ── Step 1: Use pre-loaded stream or read from path ───────────
            self.status.emit("Loading data ...")
            self.progress.emit(5)
            if self.params.get('stream') is not None:
                st = self.params['stream'].copy()
            else:
                st = read(self.params['mseed_file'])
            self.progress.emit(20)

            # ── Step 2: Load inventory ────────────────────────────────────
            inv_path = self.params['inventory_file']
            if os.path.isdir(inv_path):
                # Merge every .xml file in the directory into one inventory
                from obspy.core.inventory import Inventory
                inv = Inventory()
                xml_files = sorted(
                    os.path.join(inv_path, f)
                    for f in os.listdir(inv_path) if f.endswith('.xml'))
                for xf in xml_files:
                    try:
                        inv += read_inventory(xf)
                    except Exception:
                        pass
                self.status.emit(f"Loaded {len(xml_files)} station files ...")
            else:
                inv = read_inventory(inv_path)
            # One channel epoch per trace (latest calibration wins) so response
            # removal, sensitivity and coordinates all use the same metadata.
            inv = fk.resolve_epochs(inv, st)
            self.progress.emit(35)

            # ── Step 3: Merge and clean masked arrays ─────────────────────
            st.merge(fill_value=0)
            for tr in st:
                if isinstance(tr.data, np.ma.MaskedArray):
                    tr.data = tr.data.filled(0)
            self.progress.emit(45)

            # ── Step 4: Full instrument response removal ──────────────────
            # Pressure sensors come out in Pa (ObsPy applies a PA response
            # as-is even with output='VEL'); units are read from StationXML.
            self.status.emit("Removing instrument response ...")
            try:
                fk.require_response_stages(inv, st)
                calibrated = st.copy()
                calibrated.remove_response(
                    inventory=inv,
                    output='VEL',
                    pre_filt=(0.1, 0.5, 9.0, 10.0),
                    water_level=60,
                )
                st = calibrated
                trace_units = [fk.calibrated_units(inv, tr, 'VEL') for tr in st]
            except Exception as e:
                # Fallback: scalar sensitivity division if response removal fails
                self.status.emit(f"Response removal failed ({e}), using scalar sensitivity ...")
                trace_units = []
                for tr in st:
                    unit = "counts"
                    try:
                        resp = inv.get_response(tr.id, datetime=tr.stats.starttime)
                        sens = float(resp.instrument_sensitivity.value)
                        if np.isfinite(sens) and sens > 0:
                            tr.data = tr.data.astype(float) / sens
                            unit = fk.calibrated_units(inv, tr, 'DEF')
                    except Exception:
                        pass   # keep raw counts if no response found
                    trace_units.append(unit)
            uncalibrated = [tr.id for tr, unit in zip(st, trace_units)
                            if unit == "counts"]
            if uncalibrated:
                self.status.emit("No response for " + ", ".join(uncalibrated)
                                 + " — beam amplitude is in raw counts.")
            self.progress.emit(60)

            # ── Step 5: FK analysis ───────────────────────────────────────
            self.status.emit("Running FK analysis ...")
            self.progress.emit(65)

            result = fk_array(
                st,
                inv,
                self.params['min_freq'],
                self.params['max_freq'],
                self.params['win_length'],
                self.params['overlap'],
                self.params['start_time'],
                self.params['duration'],
                self.params['event_name'],
                self.params['event_lat'],
                self.params['event_lon'],
            )
            self.progress.emit(100)
            self.status.emit("Done.")
            # Pass semb_thresh through so _show_results can use it
            result['semb_thresh'] = self.params.get('semb_thresh', 0.3)
            result['beam_units'] = (trace_units[0] if len(set(trace_units)) == 1
                                    else "mixed units")
            self.finished.emit(result)

        except Exception as e:
            import traceback
            self.status.emit("Error.")
            self.error.emit(f"{e}\n\n{traceback.format_exc()}")


# ═══════════════════════════════════════════════════════════════════════════
#  Shared footer
# ═══════════════════════════════════════════════════════════════════════════
class FooterWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(4, 2, 4, 2)

        logo_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logo.png")
        if os.path.exists(logo_path):
            logo_label = QLabel()
            pix = QPixmap(logo_path).scaled(28, 28, Qt.KeepAspectRatio,
                                             Qt.SmoothTransformation)
            logo_label.setPixmap(pix)
            layout.addWidget(logo_label)
        layout.addStretch()

        copy_label = QLabel("© 2024-2026 Islam Hamama  |  islam.hamama@nriag.sci.eg")
        copy_label.setObjectName("footerText")
        copy_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(copy_label)


# ═══════════════════════════════════════════════════════════════════════════
#  Save-to-database dialog
# ═══════════════════════════════════════════════════════════════════════════
class SaveEventDialog(QDialog):
    """
    Shown after FK analysis to let the user classify the event, add a note,
    and persist everything to the SQLite database.
    """


    def __init__(self, result: dict, params: dict, fig=None, parent=None):
        super().__init__(parent)
        self.result = result
        self.params = params
        self.fig    = fig
        self.setWindowTitle("Save Event to Database")
        self.setMinimumWidth(620)
        self.resize(680, 650)
        apply_result_style(self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 12)
        layout.setSpacing(10)
        layout.addWidget(result_header(
            "Save event", "Classify this analysis and add a note to the archive",
            kicker="SEISMOFK / EVENT ARCHIVE", badge="REVIEW"))

        # ── Summary box ───────────────────────────────────────────────────
        summ_grp    = QGroupBox("Analysis Summary")
        summ_layout = QVBoxLayout()
        semb  = result.get('semblance', np.array([]))
        thresh = params.get('semb_thresh', 0.3)
        n_det  = int(np.sum(semb >= thresh)) if len(semb) else 0
        lines = [
            f"Event:        {params.get('event_name', '—')}",
            f"Array:        {result.get('station_name', '—')}",
            f"Start time:   {params.get('start_time', '—')}",
            f"Duration:     {params.get('duration', '—')} s",
            f"Filter:       {params.get('min_freq', '—')} – {params.get('max_freq', '—')} Hz",
            f"Med BAZ:      {result.get('med_baz', 0):.1f}°   "
            f"(expected {result.get('expected_bazi', 0):.1f}°)",
            f"Med Vel:      {result.get('med_vel', 0):.0f} m/s",
            f"Detections:   {n_det} / {len(semb)} windows",
        ]
        if params.get('origin_time'):
            lines.append(f"Origin Time:  {params['origin_time']}")
        if params.get('celerity'):
            lines.append(f"Celerity:     {params['celerity']:.1f} m/s")
        lbl = QLabel('\n'.join(lines))
        lbl.setObjectName("monoSummary")
        lbl.setFont(QFontDatabase.systemFont(QFontDatabase.FixedFont))
        summ_layout.addWidget(lbl)
        summ_grp.setLayout(summ_layout)
        layout.addWidget(summ_grp)

        # ── Classification ─────────────────────────────────────────────────
        cls_grp    = QGroupBox("Event Classification")
        cls_layout = QHBoxLayout()
        cls_layout.addWidget(QLabel("Type:"))
        self.cls_combo = QComboBox()
        self.cls_combo.addItems(db_manager.CLASSIFICATIONS)
        cls_layout.addWidget(self.cls_combo, stretch=1)
        cls_grp.setLayout(cls_layout)
        layout.addWidget(cls_grp)

        # ── Note ──────────────────────────────────────────────────────────
        note_grp    = QGroupBox("Note  (optional)")
        note_layout = QVBoxLayout()
        self.note_edit = QTextEdit()
        self.note_edit.setPlaceholderText(
            "Add observations, context, or any remarks …")
        self.note_edit.setFixedHeight(90)
        note_layout.addWidget(self.note_edit)
        note_grp.setLayout(note_layout)
        layout.addWidget(note_grp)

        # ── Buttons ───────────────────────────────────────────────────────
        btn_row  = QHBoxLayout()
        save_btn = QPushButton("Save to database")
        save_btn.setObjectName("primaryAction")
        save_btn.clicked.connect(self._save)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        btn_row.addStretch()
        btn_row.addWidget(cancel_btn)
        btn_row.addWidget(save_btn)
        layout.addLayout(btn_row)

    def _save(self):
        r  = self.result
        p  = self.params
        semb   = r.get('semblance', np.array([]))
        thresh = p.get('semb_thresh', 0.3)
        n_det  = int(np.sum(semb >= thresh)) if len(semb) else 0

        safe_name = "".join(c if c.isalnum() or c in '-_' else '_'
                            for c in str(p.get('event_name', 'event')))

        # Render figure to PNG bytes
        fig_blob = None
        if self.fig is not None:
            import io
            buf = io.BytesIO()
            try:
                self.fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
                fig_blob = buf.getvalue()
            except Exception:
                pass

        data = dict(
            event_name     = p.get('event_name', ''),
            classification = self.cls_combo.currentText(),
            note           = self.note_edit.toPlainText().strip(),
            station        = str(r.get('station_name', '')),
            start_time     = str(p.get('start_time', '')),
            duration       = float(p.get('duration', 0)),
            min_freq       = float(p.get('min_freq', 0)),
            max_freq       = float(p.get('max_freq', 0)),
            win_length     = float(p.get('win_length', 0)),
            overlap        = float(p.get('overlap', 0)),
            semb_thresh    = float(thresh),
            med_baz        = float(r.get('med_baz', 0)),
            med_vel        = float(r.get('med_vel', 0)),
            expected_baz   = float(r.get('expected_bazi', 0)),
            n_detections   = n_det,
            n_windows      = int(len(semb)),
            event_lat      = float(p.get('event_lat', 0)),
            event_lon      = float(p.get('event_lon', 0)),
            origin_time    = p.get('origin_time') or None,
            celerity       = float(p['celerity']) if p.get('celerity') else None,
            figure_path    = str(r.get('figure_path', '')),
            csv_path       = f"Output_{safe_name}.csv",
            figure_blob    = fig_blob,
        )
        try:
            row_id = db_manager.save_event(data)
            QMessageBox.information(
                self, "Saved",
                f"Event saved to database  (ID {row_id}).")
            self.accept()
        except Exception as e:
            QMessageBox.critical(self, "Save Error", str(e))


# ═══════════════════════════════════════════════════════════════════════════
#  Database browser dialog
# ═══════════════════════════════════════════════════════════════════════════
class DatabaseBrowserDialog(QDialog):
    """Browse, filter, edit and export the saved FK events database."""

    _HDR = ['ID', 'Saved At', 'Event', 'Array', 'Classification',
            'Med BAZ', 'Med Vel', 'Detections', 'Freq (Hz)', 'Note']
    _COL_KEYS = ['id', 'saved_at', 'event_name', 'station', 'classification',
                 'med_baz', 'med_vel', 'n_detections', '_freq', 'note']

    def __init__(self, parent=None):
        super().__init__(parent)
        apply_result_style(self)
        self.setWindowTitle("SeismoFK — Event Database")
        self.setMinimumSize(900, 520)
        self.resize(1150, 660)
        self._all_rows = []
        self._build_ui()
        self._load()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 12)
        layout.setSpacing(10)
        layout.addWidget(result_header(
            "Event database", "Review, classify and export archived analyses",
            kicker="SEISMOFK / EVENT ARCHIVE", badge="DATABASE"))

        # ── Filter bar ────────────────────────────────────────────────────
        filter_row = QHBoxLayout()
        filter_row.addWidget(QLabel("Filter by classification:"))
        self.filter_combo = QComboBox()
        self.filter_combo.addItem("All")
        self.filter_combo.addItems(db_manager.CLASSIFICATIONS)
        self.filter_combo.currentTextChanged.connect(self._apply_filter)
        filter_row.addWidget(self.filter_combo)
        filter_row.addStretch()

        refresh_btn = QPushButton("Refresh")
        refresh_btn.clicked.connect(self._load)
        filter_row.addWidget(refresh_btn)
        layout.addLayout(filter_row)

        # ── Table ─────────────────────────────────────────────────────────
        self.table = QTableWidget()
        self.table.setColumnCount(len(self._HDR))
        self.table.setHorizontalHeaderLabels(self._HDR)
        self.table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.table.setSelectionBehavior(QTableWidget.SelectRows)
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.verticalHeader().setVisible(False)
        layout.addWidget(self.table, stretch=1)

        # ── Row count label ───────────────────────────────────────────────
        self.count_lbl = QLabel("")
        self.count_lbl.setObjectName("mutedText")
        layout.addWidget(self.count_lbl)

        # ── Action buttons ────────────────────────────────────────────────
        btn_row = QHBoxLayout()

        view_fig_btn = QPushButton("View figure")
        view_fig_btn.setObjectName("primaryAction")
        view_fig_btn.clicked.connect(self._view_figure)
        edit_btn = QPushButton("Edit classification…")
        edit_btn.clicked.connect(self._edit_row)
        export_btn = QPushButton("Export CSV…")
        export_btn.clicked.connect(self._export)
        delete_btn = QPushButton("Delete…")
        delete_btn.setObjectName("dangerAction")
        delete_btn.clicked.connect(self._delete_row)
        close_btn  = QPushButton("Close")
        close_btn.clicked.connect(self.close)

        for b in (view_fig_btn, edit_btn, export_btn):
            btn_row.addWidget(b)
        btn_row.addStretch()
        btn_row.addWidget(delete_btn)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)

    # ── Data loading ──────────────────────────────────────────────────────
    def _load(self):
        self._all_rows = [dict(r) for r in db_manager.fetch_all()]
        self._apply_filter(self.filter_combo.currentText())

    def _apply_filter(self, cls_text):
        if cls_text == "All":
            rows = self._all_rows
        else:
            rows = [r for r in self._all_rows if r.get('classification') == cls_text]
        self._populate(rows)

    def _populate(self, rows):
        self.table.setRowCount(len(rows))
        for ri, row in enumerate(rows):
            freq_str = (f"{row.get('min_freq','?')}–{row.get('max_freq','?')}")
            row['_freq'] = freq_str
            for ci, key in enumerate(self._COL_KEYS):
                val = row.get(key, '')
                if isinstance(val, float):
                    val = f"{val:.1f}"
                item = QTableWidgetItem(str(val) if val is not None else '')
                item.setData(Qt.UserRole, row.get('id'))
                self.table.setItem(ri, ci, item)
        self.count_lbl.setText(f"{len(rows)} event(s) shown")

    # ── Actions ───────────────────────────────────────────────────────────
    def _selected_id(self):
        row = self.table.currentRow()
        if row < 0:
            QMessageBox.warning(self, "No selection", "Select a row first.")
            return None
        item = self.table.item(row, 0)
        return int(item.text()) if item else None

    def _view_figure(self):
        event_id = self._selected_id()
        if event_id is None:
            return
        row_data = next((r for r in self._all_rows if r.get('id') == event_id), {})
        blob = row_data.get('figure_blob')
        if not blob:
            QMessageBox.information(self, "No Figure",
                                    "No figure image was stored for this event.")
            return

        from PyQt5.QtGui import QPixmap
        from PyQt5.QtCore import QByteArray
        pixmap = QPixmap()
        pixmap.loadFromData(QByteArray(blob))

        dlg = QDialog(self)
        apply_result_style(dlg)
        title = (f"Figure — {row_data.get('event_name','?')}  "
                 f"[{row_data.get('saved_at','?')}]")
        dlg.setWindowTitle(title)
        dlg.setMinimumSize(900, 650)
        vlay = QVBoxLayout(dlg)
        vlay.setContentsMargins(16, 16, 16, 12)
        vlay.setSpacing(10)
        vlay.addWidget(result_header(
            row_data.get('event_name') or f"Event #{event_id}",
            f"Saved {row_data.get('saved_at', '?')}  ·  "
            f"{row_data.get('classification') or 'Unclassified'}",
            kicker="SEISMOFK / EVENT ARCHIVE", badge=f"#{event_id}"))

        scroll_lbl = QLabel()
        scroll_lbl.setPixmap(
            pixmap.scaled(1200, 900, Qt.KeepAspectRatio, Qt.SmoothTransformation))
        scroll_lbl.setAlignment(Qt.AlignCenter)

        from PyQt5.QtWidgets import QScrollArea
        scroll = QScrollArea()
        scroll.setWidget(scroll_lbl)
        scroll.setWidgetResizable(True)
        vlay.addWidget(scroll)

        # Save button
        btn_row = QHBoxLayout()
        save_btn = QPushButton("Save PNG…")
        save_btn.clicked.connect(lambda: self._save_blob_png(blob, row_data))
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dlg.close)
        btn_row.addStretch()
        btn_row.addWidget(save_btn)
        btn_row.addWidget(close_btn)
        vlay.addLayout(btn_row)
        dlg.exec_()

    def _save_blob_png(self, blob, row_data):
        default = (f"FK_{row_data.get('event_name','event')}_"
                   f"{row_data.get('id','')}.png")
        fname, _ = QFileDialog.getSaveFileName(
            self, "Save Figure", default, "PNG Files (*.png)")
        if not fname:
            return
        try:
            with open(fname, 'wb') as f:
                f.write(blob)
            QMessageBox.information(self, "Saved", f"Figure saved:\n{fname}")
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))

    def _edit_row(self):
        event_id = self._selected_id()
        if event_id is None:
            return
        # Find the row dict
        row_data = next((r for r in self._all_rows if r.get('id') == event_id), {})

        dlg = QDialog(self)
        apply_result_style(dlg)
        dlg.setWindowTitle(f"Edit Event #{event_id}")
        dlg.setMinimumWidth(480)
        vlay = QVBoxLayout(dlg)
        vlay.setContentsMargins(16, 16, 16, 12)
        vlay.setSpacing(8)
        vlay.addWidget(result_header(
            f"Edit event #{event_id}", row_data.get('event_name') or "",
            kicker="SEISMOFK / EVENT ARCHIVE", badge="EDIT"))

        vlay.addWidget(QLabel("Classification:"))
        cls_combo = QComboBox()
        cls_combo.addItems(db_manager.CLASSIFICATIONS)
        cur_cls = row_data.get('classification', '')
        if cur_cls in db_manager.CLASSIFICATIONS:
            cls_combo.setCurrentText(cur_cls)
        vlay.addWidget(cls_combo)

        vlay.addWidget(QLabel("Note:"))
        note_edit = QTextEdit()
        note_edit.setPlainText(row_data.get('note', '') or '')
        note_edit.setFixedHeight(100)
        vlay.addWidget(note_edit)

        btn_row = QHBoxLayout()
        ok_btn  = QPushButton("Update")
        ok_btn.setObjectName("primaryAction")
        ok_btn.clicked.connect(dlg.accept)
        cn_btn = QPushButton("Cancel")
        cn_btn.clicked.connect(dlg.reject)
        btn_row.addStretch(); btn_row.addWidget(cn_btn); btn_row.addWidget(ok_btn)
        vlay.addLayout(btn_row)

        if dlg.exec_() == QDialog.Accepted:
            try:
                db_manager.update_event(
                    event_id,
                    cls_combo.currentText(),
                    note_edit.toPlainText().strip()
                )
                self._load()
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    def _delete_row(self):
        event_id = self._selected_id()
        if event_id is None:
            return
        if QMessageBox.question(
                self, "Confirm Delete",
                f"Permanently delete event #{event_id}?",
                QMessageBox.Yes | QMessageBox.No) == QMessageBox.Yes:
            try:
                db_manager.delete_event(event_id)
                self._load()
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    def _export(self):
        fname, _ = QFileDialog.getSaveFileName(
            self, "Export Database to CSV", "fk_events_export.csv",
            "CSV Files (*.csv)")
        if not fname:
            return
        try:
            n = db_manager.export_csv(fname)
            QMessageBox.information(self, "Exported",
                                    f"Exported {n} event(s) to:\n{fname}")
        except Exception as e:
            QMessageBox.critical(self, "Export Error", str(e))


# ═══════════════════════════════════════════════════════════════════════════
#  Results window
# ═══════════════════════════════════════════════════════════════════════════
class ResultsWindow(QDialog):
    def __init__(self, fig, result=None, params=None, parent=None):
        super().__init__(parent)
        apply_result_style(self)
        self.setWindowTitle("SeismoFK — Analysis Results")
        self.setMinimumSize(860, 620)
        self.resize(1220, 840)
        self.fig    = fig
        self.result = result or {}
        self.params = params or {}

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 10)
        layout.setSpacing(10)
        event_name = self.params.get("event_name") or "Untitled event"
        sensors = self.result.get("n_stations")
        start = self.params.get("start_time") or ""
        layout.addWidget(result_header(
            f"FK results · {event_name}",
            f"{sensors} sensors  ·  {start}" if sensors else str(start),
            kicker="SEISMOFK / ARRAY ANALYSIS", badge="RESULTS"))

        def metric(value, spec, suffix=""):
            try:
                number = float(value)
                return f"{number:{spec}}{suffix}" if np.isfinite(number) else "—"
            except (TypeError, ValueError):
                return "—"

        semblance = np.asarray(self.result.get("semblance", []))
        threshold = float(self.params.get("semb_thresh", 0.3))
        detections = int(np.sum(semblance >= threshold))
        summary = QHBoxLayout()
        summary.setSpacing(10)
        for label, value in (
                ("Back-azimuth", metric(self.result.get("med_baz"), ".1f", "°")),
                ("Trace velocity", metric(self.result.get("med_vel"), ".0f", " m/s")),
                ("Detected windows", f"{detections} / {len(semblance)}"),
                ("Expected direction", metric(
                    self.result.get("expected_bazi"), ".1f", "°"))):
            summary.addWidget(metric_card(label, value), stretch=1)
        layout.addLayout(summary)

        plot_card = QFrame()
        plot_card.setObjectName("resultCard")
        plot_layout = QVBoxLayout(plot_card)
        plot_layout.setContentsMargins(10, 7, 10, 10)
        plot_layout.setSpacing(3)
        self.canvas  = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        style_plot_toolbar(self.toolbar)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas, stretch=1)
        layout.addWidget(plot_card, stretch=1)

        btn_row = QHBoxLayout()
        save_fig_btn = QPushButton("Export figure")
        save_fig_btn.clicked.connect(self.save_figure)

        save_db_btn = QPushButton("Save event")
        save_db_btn.setObjectName("primaryAction")
        save_db_btn.clicked.connect(self._save_to_db)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)

        btn_row.addStretch()
        btn_row.addWidget(save_fig_btn)
        btn_row.addWidget(save_db_btn)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)
        layout.addWidget(FooterWidget(self))

    def save_figure(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save Figure", "",
            "PNG Files (*.png);;PDF Files (*.pdf);;SVG Files (*.svg);;All Files (*)")
        if not filename:
            return
        if not any(filename.endswith(e) for e in ('.png', '.pdf', '.svg')):
            filename += '.png'
        try:
            self.fig.savefig(filename, dpi=300, bbox_inches='tight')
            QMessageBox.information(self, "Saved", f"Figure saved:\n{filename}")
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))

    def _save_to_db(self):
        dlg = SaveEventDialog(self.result, self.params, fig=self.fig, parent=self)
        dlg.exec_()


# ═══════════════════════════════════════════════════════════════════════════
#  Spectrogram window — multi-panel PSD spectrograms with in-window recompute
# ═══════════════════════════════════════════════════════════════════════════
class SpectrogramWindow(QDialog):
    """Multi-panel spectrogram viewer (one panel per trace, magma cmap).

    Modelled on ``ResultsWindow``: embeds a matplotlib canvas + navigation
    toolbar + a "Save Figure" button. Adds an in-window control row
    (bandpass min/max, freq-max, nperseg/noverlap) with a "Recompute" button
    that re-runs ``compute_spectrogram()`` and redraws the canvas in place.

    Parameters
    ----------
    stream : obspy.Stream
        The data to display — ALREADY sliced to the desired time window by
        the caller.
    inventory : obspy.Inventory
        Resolved inventory providing instrument responses for every trace.
    params : dict
        Initial parameters. Recognised keys (all optional, with defaults):
        ``filter_freqmin``, ``filter_freqmax``, ``freq_max``, ``nperseg``,
        ``noverlap``, ``smooth_bins``, ``reference_pressure``.
    """

    # Fixed defaults mirror plot_spectrogram_window.py's standalone values.
    _DEFAULTS = {
        "filter_freqmin": 0.5,
        "filter_freqmax": 6.0,
        "freq_max": 5.0,
        "nperseg": 512,
        "noverlap": 460,
        "smooth_bins": 1.0,
        "reference_pressure": 20e-6,
    }

    def __init__(self, stream, inventory, params=None, parent=None):
        super().__init__(parent)
        apply_result_style(self)
        self.setWindowTitle("SeismoFK — Spectrograms")
        self.setMinimumSize(900, 680)
        self.resize(1250, 900)

        self.stream    = stream
        self.inventory = inventory
        cfg = dict(self._DEFAULTS)
        if params:
            cfg.update({k: v for k, v in params.items() if v is not None})
        self.reference_pressure = float(cfg["reference_pressure"])
        from plot_spectrogram_window import get_sensitivity_counts_per_pa
        try:
            if inventory is None:
                raise ValueError("No StationXML inventory selected")
            for trace in stream:
                get_sensitivity_counts_per_pa(inventory, trace)
            self.units_mode = "pressure"
            self.calibration_note = "Calibrated pressure"
        except Exception as error:
            self.units_mode = "counts"
            self.calibration_note = f"Raw counts · pressure calibration unavailable: {error}"

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 10)
        layout.setSpacing(10)
        if stream:
            first_sample = min(tr.stats.starttime for tr in stream)
            last_sample = max(tr.stats.endtime for tr in stream)
            data_range = (
                f"{len(stream)} channels  ·  "
                f"{first_sample.strftime('%Y-%m-%d %H:%M:%S')} to "
                f"{last_sample.strftime('%H:%M:%S')} UTC  "
                f"({last_sample - first_sample:.0f} s)")
        else:
            data_range = "No waveform data loaded"
        layout.addWidget(result_header(
            "Spectrograms", data_range,
            kicker="SEISMOFK / TIME–FREQUENCY",
            badge="Pa" if self.units_mode == "pressure" else "COUNTS"))

        # ── In-window control row (Option C) ──────────────────────────────
        ctrl_group  = QFrame()
        ctrl_group.setObjectName("resultCard")
        ctrl_layout = QHBoxLayout(ctrl_group)
        ctrl_layout.setContentsMargins(14, 10, 14, 10)
        ctrl_layout.setSpacing(8)
        self.freqmin_input  = QLineEdit(f"{cfg['filter_freqmin']:g}")
        self.freqmax_input  = QLineEdit(f"{cfg['filter_freqmax']:g}")
        self.freq_max_input = QLineEdit(f"{cfg['freq_max']:g}")
        self.nperseg_input  = QLineEdit(str(int(cfg["nperseg"])))
        self.noverlap_input = QLineEdit(str(int(cfg["noverlap"])))
        self.smoothing_input = QDoubleSpinBox()
        self.smoothing_input.setRange(0.0, 3.0)
        self.smoothing_input.setDecimals(1)
        self.smoothing_input.setSingleStep(0.2)
        self.smoothing_input.setValue(float(cfg["smooth_bins"]))
        self.smoothing_input.setToolTip(
            "Smooth adjacent PSD bins before display. Set 0 for the original bins.")
        for w in (self.freqmin_input, self.freqmax_input, self.freq_max_input,
                  self.nperseg_input, self.noverlap_input, self.smoothing_input):
            w.setMaximumWidth(80)
        recompute_btn = QPushButton("Recompute")
        recompute_btn.setObjectName("primaryAction")
        recompute_btn.clicked.connect(self._recompute)
        for w in (QLabel("Bandpass min (Hz):"), self.freqmin_input,
                  QLabel("max (Hz):"),          self.freqmax_input,
                  QLabel("Freq. max (Hz):"),    self.freq_max_input,
                  QLabel("nperseg:"),           self.nperseg_input,
                  QLabel("noverlap:"),          self.noverlap_input,
                  QLabel("Smooth:"),            self.smoothing_input,
                  recompute_btn):
            ctrl_layout.addWidget(w)
        ctrl_layout.addStretch()
        layout.addWidget(ctrl_group)

        self.feedback = QLabel(self.calibration_note)
        self.feedback.setObjectName("mutedText")
        self.feedback.setWordWrap(True)
        layout.addWidget(self.feedback)

        # ── Canvas + toolbar ──────────────────────────────────────────────
        # Use Figure() directly (not plt.figure()) so the figure is NOT
        # registered in pyplot's global figure manager — otherwise every
        # "Plot Spectrogram" invocation would orphan and leak a figure for
        # the process lifetime. This is the documented matplotlib-in-Qt
        # pattern; the canvas/toolbar own the figure's lifecycle.
        self.fig    = Figure(figsize=(13, 9), facecolor="white")
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
        layout.addWidget(plot_card, stretch=1)

        # ── Action buttons ────────────────────────────────────────────────
        btn_row = QHBoxLayout()
        save_fig_btn = QPushButton("Export figure")
        save_fig_btn.clicked.connect(self.save_figure)
        save_data_btn = QPushButton("Save plotted data")
        save_data_btn.setToolTip(
            "Export UTC times, frequencies, raw PSD, displayed PSD, and processing settings.")
        save_data_btn.clicked.connect(self.save_data)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        btn_row.addStretch()
        btn_row.addWidget(save_data_btn)
        btn_row.addWidget(save_fig_btn)
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)
        layout.addWidget(FooterWidget(self))

        # Initial render using the seeded parameters.
        self._render(cfg["filter_freqmin"], cfg["filter_freqmax"],
                     cfg["freq_max"], int(cfg["nperseg"]), int(cfg["noverlap"]),
                     float(cfg["smooth_bins"]))

    # ── Parameter parsing ─────────────────────────────────────────────────
    def _read_params(self):
        """Validate and return the control-row params, or raise ValueError."""
        try:
            freqmin  = float(self.freqmin_input.text())
            freqmax  = float(self.freqmax_input.text())
            freq_max = float(self.freq_max_input.text())
            nperseg  = int(self.nperseg_input.text())
            noverlap = int(self.noverlap_input.text())
        except ValueError:
            raise ValueError("All parameter fields must be numeric.")
        if freqmin <= 0:
            raise ValueError("Bandpass min must be > 0 Hz.")
        if freqmin >= freqmax:
            raise ValueError("Bandpass min must be < bandpass max.")
        if freq_max <= 0:
            raise ValueError("Freq. max must be > 0 Hz.")
        if nperseg <= 0:
            raise ValueError("nperseg must be a positive integer.")
        if noverlap < 0 or noverlap >= nperseg:
            raise ValueError("noverlap must satisfy 0 <= noverlap < nperseg.")
        return (freqmin, freqmax, freq_max, nperseg, noverlap,
                self.smoothing_input.value())

    def _recompute(self):
        try:
            params = self._read_params()
        except ValueError as e:
            QMessageBox.warning(self, "Invalid parameters", str(e))
            return
        self._render(*params)

    # ── Rendering ─────────────────────────────────────────────────────────
    def _render(self, filter_freqmin, filter_freqmax, freq_max,
                nperseg, noverlap, smooth_bins):
        """Compute spectrograms for every trace and (re)draw the canvas."""
        from plot_spectrogram_window import compute_spectrogram

        self.fig.clear()

        panels = []
        skipped = []
        for trace in self.stream:
            try:
                freqs, time_nums, power_db = compute_spectrogram(
                    trace,
                    self.inventory,
                    filter_freqmin=filter_freqmin,
                    filter_freqmax=filter_freqmax,
                    freq_max=freq_max,
                    nperseg=nperseg,
                    noverlap=noverlap,
                    reference_pressure=self.reference_pressure,
                    units=self.units_mode,
                    smooth_bins=smooth_bins,
                )
            except Exception as e:               # noqa: BLE001 — per-trace isolation
                skipped.append(f"{trace.id}: {e}")
                continue
            panels.append((trace, freqs, time_nums, power_db))

        if not panels:
            self._panels = []
            self.fig.text(0.5, 0.52, "NO SPECTROGRAM DATA",
                          ha="center", va="center", fontsize=16,
                          color="#34556a", fontweight="bold")
            self.fig.text(0.5, 0.45,
                          "Check the selected window, filter band, and sample rate.",
                          ha="center", va="center", fontsize=10,
                          color="#78909e")
            self.feedback.setText("No traces could be plotted. " + "; ".join(skipped))
            self.canvas.draw()
            return

        self._panels = panels
        self._last_plot_params = dict(
            filter_freqmin=filter_freqmin, filter_freqmax=filter_freqmax,
            freq_max=freq_max, nperseg=nperseg, noverlap=noverlap,
            smooth_bins=smooth_bins, units=self.units_mode,
            reference_pressure_pa=self.reference_pressure,
            psd_reference=(f"({self.reference_pressure * 1e6:g} µPa)²/Hz"
                           if self.units_mode == "pressure" else "1 count²/Hz"),
            smoothing_method="Gaussian on linear PSD before dB conversion")

        vmin = min(np.percentile(p[3], 15) for p in panels)
        vmax = max(np.percentile(p[3], 98) for p in panels)

        gs = gridspec.GridSpec(
            len(panels), 2, width_ratios=[40, 1.8],
            hspace=0.14, wspace=0.08,
            left=0.08, right=0.91, top=0.88, bottom=0.08, figure=self.fig)
        cax = self.fig.add_subplot(gs[:, 1])

        axes = []
        last_mesh = None
        for idx, (trace, freqs, time_nums, power_db) in enumerate(panels):
            ax = self.fig.add_subplot(
                gs[idx, 0], sharex=axes[0] if axes else None)
            axes.append(ax)
            last_mesh = ax.pcolormesh(
                time_nums, freqs, power_db,
                shading="gouraud" if min(len(freqs), len(time_nums)) > 1
                else "auto", cmap="magma", vmin=vmin, vmax=vmax)
            ax.set_ylabel("Hz", rotation=0, labelpad=18)
            ax.set_ylim(0, freq_max)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.grid(True, color="white", linewidth=0.6, alpha=0.25)
            ax.text(
                0.015, 0.92, trace.id, transform=ax.transAxes,
                ha="left", va="top", fontsize=9, color="white",
                bbox=dict(facecolor="black", alpha=0.28,
                          edgecolor="none", pad=3.0))

        self.fig.suptitle(
            f"Time–frequency power  ·  {len(panels)} channels  ·  "
            f"{filter_freqmin:g}–{filter_freqmax:g} Hz bandpass  ·  "
            + (f"Gaussian σ={smooth_bins:g} bins" if smooth_bins
               else "unsmoothed"),
            fontsize=13, fontweight="semibold", y=0.965)

        axes[-1].set_xlabel("Time (UTC)")
        locator   = mdates.AutoDateLocator()
        formatter = mdates.DateFormatter("%H:%M:%S")
        axes[-1].xaxis.set_major_locator(locator)
        axes[-1].xaxis.set_major_formatter(formatter)
        for ax in axes[:-1]:
            plt.setp(ax.get_xticklabels(), visible=False)

        if last_mesh is not None:
            cbar = self.fig.colorbar(last_mesh, cax=cax)
            # compute_spectrogram divides by reference², i.e. dB re (p_ref)²/Hz.
            cbar.set_label(f"PSD (dB re ({self.reference_pressure * 1e6:g} µPa)²/Hz)"
                           if self.units_mode == "pressure"
                           else "PSD (dB re 1 count²/Hz)")
            cbar.outline.set_linewidth(0.8)

        self.canvas.draw()

        self.feedback.setText(
            f"{self.calibration_note} · {len(panels)} trace(s) plotted · "
            f"smoothing {smooth_bins:g} bins"
            + (f" · {len(skipped)} skipped (hover for details)" if skipped else ""))
        self.feedback.setToolTip("\n".join(skipped))

    def save_data(self):
        if not getattr(self, "_panels", None):
            QMessageBox.warning(self, "No plotted data", "Recompute a valid spectrogram first.")
            return
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save spectrogram data", "spectrogram_data.npz",
            "NumPy archive (*.npz)")
        if not filename:
            return
        if not filename.lower().endswith(".npz"):
            filename += ".npz"
        from plot_spectrogram_window import compute_spectrogram

        try:
            arrays = {}
            metadata = dict(self._last_plot_params)
            metadata["traces"] = []
            for index, (trace, freqs, time_nums, displayed_db) in enumerate(self._panels):
                _, _, raw_db = compute_spectrogram(
                    trace, self.inventory,
                    filter_freqmin=metadata["filter_freqmin"],
                    filter_freqmax=metadata["filter_freqmax"],
                    freq_max=metadata["freq_max"],
                    nperseg=metadata["nperseg"], noverlap=metadata["noverlap"],
                    reference_pressure=self.reference_pressure,
                    units=self.units_mode, smooth_bins=0)
                key = f"trace_{index}"
                metadata["traces"].append({"key": key, "id": trace.id})
                arrays[f"{key}_frequency_hz"] = freqs
                arrays[f"{key}_time_utc"] = np.array(
                    [mdates.num2date(value).isoformat() for value in time_nums])
                arrays[f"{key}_raw_psd_db"] = raw_db
                arrays[f"{key}_display_psd_db"] = displayed_db
            arrays["metadata_json"] = np.array(json.dumps(metadata, sort_keys=True))
            np.savez_compressed(filename, **arrays)
            QMessageBox.information(self, "Saved", f"Spectrogram data saved:\n{filename}")
        except Exception as error:
            QMessageBox.critical(self, "Save error", str(error))

    def save_figure(self):
        filename, _ = QFileDialog.getSaveFileName(
            self, "Save Figure", "",
            "PNG Files (*.png);;PDF Files (*.pdf);;SVG Files (*.svg);;All Files (*)")
        if not filename:
            return
        if not any(filename.endswith(e) for e in ('.png', '.pdf', '.svg')):
            filename += '.png'
        try:
            self.fig.savefig(filename, dpi=300, bbox_inches='tight')
            QMessageBox.information(self, "Saved", f"Figure saved:\n{filename}")
        except Exception as e:
            QMessageBox.critical(self, "Error", str(e))


# ═══════════════════════════════════════════════════════════════════════════
#  Waveform viewer / time-picker dialog
# ═══════════════════════════════════════════════════════════════════════════
class WaveformViewer(QDialog):
    time_selected = pyqtSignal(str)

    def __init__(self, stream, parent=None):
        super().__init__(parent)
        apply_result_style(self)
        self.setWindowTitle("SeismoFK — Waveform Viewer")
        self.setMinimumSize(900, 560)
        self.resize(1250, 720)
        self.stream          = stream
        self.original_stream = stream.copy()
        self.selected_time   = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 10)
        layout.setSpacing(10)
        layout.addWidget(result_header(
            "Pick analysis start",
            "Click a waveform to set the analysis window start time",
            kicker="SEISMOFK / WAVEFORMS", badge=f"{len(stream)} TRACES"))

        # ── Filter controls ───────────────────────────────────────────────
        filter_group  = QFrame()
        filter_group.setObjectName("resultCard")
        filter_layout = QHBoxLayout(filter_group)
        filter_layout.setContentsMargins(14, 10, 14, 10)
        self.low_freq_input  = QLineEdit("0.5")
        self.high_freq_input = QLineEdit("5.0")
        apply_btn = QPushButton("Apply Filter")
        reset_btn = QPushButton("Reset")
        apply_btn.clicked.connect(self.apply_filter)
        reset_btn.clicked.connect(self.reset_filter)
        for w in (QLabel("Band-pass low (Hz):"), self.low_freq_input,
                  QLabel("high (Hz):"), self.high_freq_input,
                  apply_btn, reset_btn):
            filter_layout.addWidget(w)
        for w in (self.low_freq_input, self.high_freq_input):
            w.setMaximumWidth(80)
        filter_layout.addStretch()
        layout.addWidget(filter_group)

        # Figure() rather than plt.subplots() so each viewer does not leave a
        # figure registered in pyplot for the life of the process.
        self.fig_wv = Figure(figsize=(13, 5), facecolor='white')
        self.ax_wv = self.fig_wv.add_subplot(111)
        self.canvas_wv  = FigureCanvas(self.fig_wv)
        self.toolbar_wv = NavigationToolbar(self.canvas_wv, self)
        style_plot_toolbar(self.toolbar_wv)
        plot_card = QFrame()
        plot_card.setObjectName("resultCard")
        plot_layout = QVBoxLayout(plot_card)
        plot_layout.setContentsMargins(10, 7, 10, 10)
        plot_layout.addWidget(self.toolbar_wv)
        plot_layout.addWidget(self.canvas_wv, stretch=1)
        layout.addWidget(plot_card, stretch=1)

        self._draw_waveforms(self.stream)
        self.canvas_wv.mpl_connect('button_press_event', self._on_click)

        btn_row = QHBoxLayout()
        self.time_label = QLabel("Selected time: —")
        self.time_label.setObjectName("resultInfo")
        self.confirm_btn = QPushButton("Use this start time")
        self.confirm_btn.setObjectName("primaryAction")
        self.confirm_btn.setEnabled(False)
        self.confirm_btn.clicked.connect(self.confirm_selection)
        cancel_btn  = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)
        btn_row.addWidget(self.time_label, stretch=1)
        btn_row.addWidget(cancel_btn)
        btn_row.addWidget(self.confirm_btn)
        layout.addLayout(btn_row)
        layout.addWidget(FooterWidget(self))

    def _draw_waveforms(self, stream, title="Waveforms"):
        """Same style as the main waveform explorer: one ink colour, traces
        stacked top to bottom with station labels, UTC time axis; x values
        stay seconds from the first trace start (used by the click pick)."""
        ax = self.ax_wv
        ax.clear()
        t0 = stream[0].stats.starttime
        n = len(stream)
        for i, tr in enumerate(stream):
            x, y = _envelope(tr.data, tr.stats.delta,
                             float(tr.stats.starttime - t0))
            ax.plot(x, 0.45 * y / (np.max(np.abs(y)) or 1.0) + (n - 1 - i),
                    color=PALETTE['ink_2'], lw=0.55)
        ax.set_yticks(range(n))
        ax.set_yticklabels([f"{tr.stats.station} {tr.stats.channel}"
                            for tr in reversed(stream)], fontsize=8)
        ax.set_ylim(-0.7, n - 0.3)
        for side in ('top', 'right', 'left'):
            ax.spines[side].set_visible(False)
        ax.tick_params(axis='y', length=0)
        ax.grid(axis='x', color='#eef2f5', lw=0.6)
        ax.xaxis.set_major_formatter(
            FuncFormatter(lambda x, _: (t0 + x).strftime('%H:%M:%S')))
        ax.set_xlabel(f"UTC  ·  {t0.strftime('%Y-%m-%d')}", fontsize=8)
        ax.set_title(title, loc='left', fontsize=10)
        if self.selected_time is not None:
            ax.axvline(self.selected_time - t0, color=PALETTE['accent'],
                       lw=1.4)
        self.fig_wv.tight_layout(pad=0.6)
        self.canvas_wv.draw()

    def _on_click(self, event):
        if event.inaxes != self.ax_wv or event.button != 1:
            return
        self.selected_time = self.stream[0].stats.starttime + event.xdata
        self.time_label.setText(f"Selected time: {self.selected_time}")
        self.confirm_btn.setEnabled(True)
        self._draw_waveforms(self.stream,
                             title=f"Waveforms  [pick: {self.selected_time}]")

    def apply_filter(self):
        try:
            lo = float(self.low_freq_input.text())
            hi = float(self.high_freq_input.text())
            if lo >= hi:
                raise ValueError("Low must be < High")
            self.stream = self.original_stream.copy()
            self.stream.filter('bandpass', freqmin=lo, freqmax=hi, corners=4)
            self._draw_waveforms(self.stream, f"Filtered {lo}–{hi} Hz")
        except ValueError as e:
            QMessageBox.critical(self, "Filter Error", str(e))

    def reset_filter(self):
        self.stream = self.original_stream.copy()
        self._draw_waveforms(self.stream, "Waveforms (original)")

    def confirm_selection(self):
        if self.selected_time is None:
            QMessageBox.warning(self, "No pick",
                                "Click the waveform to pick a time first.")
            return
        self.time_selected.emit(str(self.selected_time))
        self.accept()


# ═══════════════════════════════════════════════════════════════════════════
#  Main window
# ═══════════════════════════════════════════════════════════════════════════
class FKAnalysisGUI(QMainWindow):


    def __init__(self):
        super().__init__()
        apply_main_style(self)
        self.setWindowTitle("SeismoFK  —  Infrasound FK Array Analysis")
        self.setGeometry(80, 60, 1280, 920)

        self.stream            = None
        self.inventory_file    = None
        self.process_thread    = None
        self.selected_time     = None
        self._xml_creator_win  = None

        self._build_ui()
        self.refresh_inventory_list()

    # ── UI construction ───────────────────────────────────────────────────
    def _build_ui(self):
        root   = QWidget()
        layout = QVBoxLayout(root)
        layout.setSpacing(12)
        layout.setContentsMargins(16, 14, 16, 10)
        self.setCentralWidget(root)

        hero = QFrame()
        hero.setObjectName("hero")
        hero_row = QHBoxLayout(hero)
        hero_row.setContentsMargins(21, 13, 21, 13)
        hero_text = QVBoxLayout()
        title = QLabel("SeismoFK")
        title.setObjectName("heroTitle")
        subtitle = QLabel("INFRASOUND ARRAY ANALYSIS")
        subtitle.setObjectName("heroSubtitle")
        _track(subtitle, 1.6)
        hero_text.setSpacing(1)
        hero_text.addWidget(title)
        hero_text.addWidget(subtitle)
        hero_row.addLayout(hero_text)
        hero_row.addStretch()

        # Workspace tools live in the header band, away from analysis actions.
        xml_creator_btn = QPushButton("StationXML editor")
        xml_creator_btn.clicked.connect(self._open_xml_creator)
        db_btn = QPushButton("Event database")
        db_btn.clicked.connect(self._open_db_browser)
        self.spectro_btn = QPushButton("Spectrogram")
        self.spectro_btn.clicked.connect(self._open_spectrogram)
        for button in (self.spectro_btn, db_btn, xml_creator_btn):
            button.setObjectName("heroAction")
            hero_row.addWidget(button)
        hero_row.addSpacing(6)
        badge = QLabel("v1.2.1")
        badge.setObjectName("versionBadge")
        hero_row.addWidget(badge)
        layout.addWidget(hero)

        content = QHBoxLayout()
        content.setSpacing(14)
        layout.addLayout(content, stretch=1)

        sidebar_scroll = QScrollArea()
        sidebar_scroll.setWidgetResizable(True)
        sidebar_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        sidebar_scroll.setFixedWidth(390)
        sidebar = QWidget()
        sidebar.setObjectName("sidebar")
        side_layout = QVBoxLayout(sidebar)
        side_layout.setContentsMargins(0, 0, 5, 0)
        side_layout.setSpacing(12)
        sidebar_scroll.setWidget(sidebar)
        # Let the window background show behind the section titles.
        sidebar_scroll.viewport().setAutoFillBackground(False)
        sidebar.setAutoFillBackground(False)
        content.addWidget(sidebar_scroll)

        workspace = QVBoxLayout()
        workspace.setSpacing(11)
        content.addLayout(workspace, stretch=1)

        # ── 1. File selection ─────────────────────────────────────────────
        file_grp    = QGroupBox("DATA SOURCE")
        file_layout = QVBoxLayout()
        file_layout.setSpacing(8)

        mseed_row = QHBoxLayout()
        self.mseed_edit = QLineEdit()
        self.mseed_edit.setReadOnly(True)
        self.mseed_edit.setPlaceholderText("Select a MiniSEED file …")
        browse_mseed = QPushButton("Browse")
        browse_mseed.clicked.connect(self.load_mseed)
        mseed_caption = QLabel("Waveform  ·  MiniSEED")
        mseed_caption.setObjectName("fieldLabel")
        file_layout.addWidget(mseed_caption)
        mseed_row.addWidget(self.mseed_edit, stretch=1)
        mseed_row.addWidget(browse_mseed)
        file_layout.addLayout(mseed_row)

        inv_row = QHBoxLayout()
        self.inv_combo = QComboBox()
        self.inv_combo.setMinimumWidth(0)
        self.inv_combo.setPlaceholderText("Select matching StationXML inventory")
        self.inv_combo.setItemDelegate(QStyledItemDelegate())
        refresh_btn = QPushButton("Refresh")
        refresh_btn.clicked.connect(self.refresh_inventory_list)
        add_xml_btn = QPushButton("+ Add XML")
        add_xml_btn.clicked.connect(self.add_inventory_file)
        inv_caption = QLabel("Station metadata  ·  StationXML")
        inv_caption.setObjectName("fieldLabel")
        file_layout.addWidget(inv_caption)
        inv_row.addWidget(self.inv_combo, stretch=1)
        file_layout.addLayout(inv_row)
        inv_actions = QHBoxLayout()
        inv_actions.addWidget(refresh_btn)
        inv_actions.addWidget(add_xml_btn)
        file_layout.addLayout(inv_actions)
        self.check_data_btn = QPushButton("Check data readiness")
        self.check_data_btn.setToolTip(
            "Check timing, sample rates, array geometry, and StationXML before analysis.")
        self.check_data_btn.clicked.connect(self._show_readiness)
        file_layout.addWidget(self.check_data_btn)
        file_grp.setLayout(file_layout)
        side_layout.addWidget(file_grp)

        # ── 2. Waveform preview ───────────────────────────────────────────
        wave_grp    = QGroupBox("WAVEFORM EXPLORER")
        wave_layout = QVBoxLayout()
        wave_hint = QLabel("Click a trace to set the analysis start time. Drag to zoom or use the toolbar.")
        wave_hint.setObjectName("sectionHint")
        wave_layout.addWidget(wave_hint)

        self.figure = Figure(figsize=(11, 3), facecolor='white')
        self.figure.text(0.5, 0.54, "No waveform loaded",
                         ha='center', va='center', color=PALETTE['ink'],
                         fontsize=14, fontweight='semibold')
        self.figure.text(0.5, 0.45,
                         "Browse for a MiniSEED file to inspect traces and pick an analysis window.",
                         ha='center', va='center', color=PALETTE['muted'], fontsize=10)
        self._preview_axes = []
        self._window_artists = []
        self.canvas  = FigureCanvas(self.figure)
        # Connected once; _update_preview only redraws.
        self.canvas.mpl_connect('button_press_event', self._on_preview_click)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.toolbar = NavigationToolbar(self.canvas, self)
        style_plot_toolbar(self.toolbar)
        wave_layout.addWidget(self.toolbar)
        wave_layout.addWidget(self.canvas)

        pick_row = QHBoxLayout()
        self.time_label = QLabel("Pick: —")
        self.time_label.setObjectName("pickLabel")
        self.clear_pick_btn = QPushButton("Clear pick")
        self.clear_pick_btn.clicked.connect(self.clear_time_pick)
        self.clear_pick_btn.setEnabled(False)
        self.view_window_btn = QPushButton("Open larger")
        self.view_window_btn.clicked.connect(self.view_waveform)
        self.view_window_btn.setEnabled(False)
        pick_row.addWidget(self.time_label, stretch=1)
        pick_row.addWidget(self.clear_pick_btn)
        pick_row.addWidget(self.view_window_btn)
        wave_layout.addLayout(pick_row)
        wave_grp.setLayout(wave_layout)
        workspace.addWidget(wave_grp, stretch=1)

        # ── 3. Parameters ─────────────────────────────────────────────────
        param_grp    = QGroupBox("FK CONFIGURATION")
        param_layout = QVBoxLayout()
        param_layout.setSpacing(10)

        self.min_freq = QDoubleSpinBox()
        self.min_freq.setRange(0.01, 100); self.min_freq.setValue(0.5)
        self.max_freq = QDoubleSpinBox()
        self.max_freq.setRange(0.01, 100); self.max_freq.setValue(6.0)
        self.win_len = QDoubleSpinBox()
        self.win_len.setRange(1, 3600); self.win_len.setValue(20.0)
        self.win_len.setSuffix(" s")
        self.overlap = QDoubleSpinBox()
        self.overlap.setRange(0.01, 0.99); self.overlap.setValue(0.1)
        self.overlap.setSingleStep(0.05)
        self.semb_thresh = QDoubleSpinBox()
        self.semb_thresh.setRange(0.0, 1.0); self.semb_thresh.setValue(0.3)
        self.semb_thresh.setSingleStep(0.05); self.semb_thresh.setDecimals(2)
        def fields_row(fields):
            row = QHBoxLayout()
            row.setSpacing(8)
            for label, widget in fields:
                column = QVBoxLayout()
                column.setSpacing(4)
                caption = QLabel(label)
                caption.setObjectName("fieldLabel")
                column.addWidget(caption)
                widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
                column.addWidget(widget)
                row.addLayout(column, stretch=1)
            return row

        param_layout.addLayout(fields_row([
            ("Minimum frequency · Hz", self.min_freq),
            ("Maximum frequency · Hz", self.max_freq),
        ]))
        param_layout.addLayout(fields_row([
            ("Window length", self.win_len),
            ("Window step", self.overlap),
        ]))
        param_layout.addLayout(fields_row([
            ("Semblance threshold", self.semb_thresh),
        ]))

        self.start_time_input = QDateTimeEdit(QDateTime.currentDateTimeUtc())
        self.start_time_input.setTimeSpec(Qt.UTC)
        self.start_time_input.setDisplayFormat("yyyy-MM-dd HH:mm:ss")
        self.duration = QDoubleSpinBox()
        self.duration.setRange(1, 86400); self.duration.setValue(900)
        self.duration.setSuffix(" s")
        self.event_name = QLineEdit("Event")
        self.event_lat  = QDoubleSpinBox()
        self.event_lat.setRange(-90, 90); self.event_lat.setDecimals(4)
        self.event_lon  = QDoubleSpinBox()
        self.event_lon.setRange(-180, 180); self.event_lon.setDecimals(4)
        param_layout.addLayout(fields_row([("Start time · UTC", self.start_time_input)]))
        param_layout.addLayout(fields_row([
            ("Duration", self.duration), ("Event name", self.event_name),
        ]))
        param_layout.addLayout(fields_row([
            ("Event latitude", self.event_lat),
            ("Event longitude", self.event_lon),
        ]))
        # Optional origin time and propagation speed.

        self.origin_time_chk = QCheckBox("Origin Time:")
        self.origin_time_chk.setToolTip(
            "Known event origin time — used together with Celerity\n"
            "to draw the expected infrasound arrival on the beam.")
        self.origin_time_input = QDateTimeEdit(QDateTime.currentDateTimeUtc())
        self.origin_time_input.setTimeSpec(Qt.UTC)
        self.origin_time_input.setDisplayFormat("yyyy-MM-dd HH:mm:ss")
        self.origin_time_input.setEnabled(False)
        self.origin_time_chk.toggled.connect(self.origin_time_input.setEnabled)

        self.celerity_chk = QCheckBox("Celerity (m/s):")
        self.celerity_chk.setToolTip(
            "Infrasound propagation speed.  Typical range: 220–340 m/s.\n"
            "Requires Origin Time to be set.")
        self.celerity_spin = QDoubleSpinBox()
        self.celerity_spin.setRange(100.0, 500.0)
        self.celerity_spin.setValue(340.0)
        self.celerity_spin.setDecimals(1)
        self.celerity_spin.setSingleStep(5.0)
        self.celerity_spin.setEnabled(False)
        self.celerity_chk.toggled.connect(self.celerity_spin.setEnabled)

        param_layout.addWidget(self.origin_time_chk)
        param_layout.addWidget(self.origin_time_input)
        param_layout.addWidget(self.celerity_chk)
        param_layout.addWidget(self.celerity_spin)
        arrival_note = QLabel("Optional: show expected arrival on the beam trace.")
        arrival_note.setObjectName("sectionHint")
        param_layout.addWidget(arrival_note)

        # Bootstrap uncertainty in the result summary.
        self.bootstrap_chk = QCheckBox("Estimate uncertainty (bootstrap)")
        self.bootstrap_chk.setToolTip(
            "Bootstrap the detected back-azimuth / apparent velocity to report\n"
            "baz ± σ and vel ± σ, and overlay the uncertainty on the polar map.")
        param_layout.addWidget(self.bootstrap_chk)

        param_grp.setLayout(param_layout)
        side_layout.addWidget(param_grp)
        side_layout.addStretch()

        # ── 4. Status + Progress + Run ────────────────────────────────────
        run_frame = QFrame()
        run_frame.setObjectName("resultCard")
        run_layout = QHBoxLayout(run_frame)
        run_layout.setContentsMargins(15, 12, 15, 12)
        run_info = QVBoxLayout()
        run_title = QLabel("FK BEAMFORMING")
        run_title.setObjectName("sectionTitle")
        run_info.addWidget(run_title)
        self.status_label = QLabel("Load a waveform and select StationXML")
        self.status_label.setObjectName("mutedText")
        run_info.addWidget(self.status_label)
        run_layout.addLayout(run_info, stretch=1)

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        workspace.addWidget(self.progress_bar)

        self.process_btn = QPushButton("Run FK analysis  →")
        self.process_btn.setObjectName("primaryAction")
        self.process_btn.clicked.connect(self.process_data)
        run_layout.addWidget(self.process_btn)
        workspace.addWidget(run_frame)

        # ── Advanced array methods (v1.2.1) ───────────────────────────────
        methods_frame = QGroupBox("ARRAY METHODS")
        methods_layout = QVBoxLayout(methods_frame)
        methods_hint = QLabel("Explore array geometry and time–frequency detections")
        methods_hint.setObjectName("sectionHint")
        methods_layout.addWidget(methods_hint)
        methods_row = QHBoxLayout()

        self.arf_btn = QPushButton("Array response")
        arf_btn = self.arf_btn
        arf_btn.setToolTip("Theoretical array response function (ARF) for the "
                           "current geometry — shows aliasing and resolution.")
        arf_btn.setObjectName("methodAction")
        arf_btn.clicked.connect(self._open_array_response)
        methods_row.addWidget(arf_btn)

        self.slowness_btn = QPushButton("Slowness map · Capon / MUSIC")
        slowness_btn = self.slowness_btn
        slowness_btn.setToolTip("High-resolution adaptive beamforming slowness "
                                "map for the picked window.")
        slowness_btn.setObjectName("methodAction")
        slowness_btn.clicked.connect(self._open_slowness_map)
        methods_row.addWidget(slowness_btn)

        self.pmcc_btn = QPushButton("PMCC detector")
        pmcc_btn = self.pmcc_btn
        pmcc_btn.setToolTip("PMCC-style time–frequency detection pixels with "
                            "back-azimuth and trace velocity estimates.")
        pmcc_btn.setObjectName("methodAction")
        pmcc_btn.clicked.connect(self._open_pmcc)
        methods_row.addWidget(pmcc_btn)

        self.noise_btn = QPushButton("Noise levels")
        self.noise_btn.setToolTip(
            "RMS noise level per sensor in 1-minute (adjustable) windows over "
            "the whole record, in dB re 20 µPa, with L90/L50/L10/Leq and "
            "sensor checks.")
        self.noise_btn.setObjectName("methodAction")
        self.noise_btn.clicked.connect(self._open_noise_levels)
        methods_row.addWidget(self.noise_btn)
        methods_layout.addLayout(methods_row)
        workspace.addWidget(methods_frame)

        self.inv_combo.currentIndexChanged.connect(self._update_analysis_actions)
        self.start_time_input.dateTimeChanged.connect(self._draw_analysis_window)
        self.duration.valueChanged.connect(self._draw_analysis_window)
        self._update_analysis_actions()

        layout.addWidget(FooterWidget(self))

    # ── XML Creator launcher ──────────────────────────────────────────────
    def _open_xml_creator(self):
        from xml_creator import XMLCreatorGUI
        if self._xml_creator_win is None or not self._xml_creator_win.isVisible():
            self._xml_creator_win = XMLCreatorGUI()
            # Refresh inventory list when creator is closed
            self._xml_creator_win.destroyed.connect(self.refresh_inventory_list)
        self._xml_creator_win.show()
        self._xml_creator_win.raise_()

    def _open_db_browser(self):
        dlg = DatabaseBrowserDialog(parent=self)
        dlg.exec_()

    # ── Spectrogram launcher ──────────────────────────────────────────────
    def _resolve_inventory(self, resolve_epochs=True):
        """Return the selected inventory; with a stream loaded and
        *resolve_epochs* set, reduced to one channel epoch per trace (see
        ``fk_analysis.resolve_epochs``)."""
        inv = self._read_selected_inventory()
        if resolve_epochs and self.stream:
            inv = fk.resolve_epochs(inv, self.stream)
        return inv

    def _read_selected_inventory(self):
        """Read the inventory selected in inv_combo into an Inventory.

        The combo's UserRole data is either a single .xml file path or a
        directory ("★ All IMS Stations") whose .xml files are merged.
        Returns an obspy.Inventory, or raises ValueError with a friendly
        message describing what is missing/invalid.
        """
        if self.inv_combo.currentIndex() < 0:
            raise ValueError("No inventory is selected. Pick one from the "
                             "inventory list (or add an XML file).")
        item = self.inv_combo.model().item(self.inv_combo.currentIndex())
        inv_path = item.data(Qt.UserRole) if item is not None else None
        if not inv_path:
            raise ValueError("The selected inventory entry has no file path.")

        if os.path.isdir(inv_path):
            from obspy.core.inventory import Inventory
            inv = Inventory()
            xml_files = sorted(
                os.path.join(inv_path, f)
                for f in os.listdir(inv_path) if f.endswith('.xml'))
            if not xml_files:
                raise ValueError(f"No .xml files found in:\n{inv_path}")
            loaded = 0
            for xf in xml_files:
                try:
                    inv += read_inventory(xf)
                    loaded += 1
                except Exception:
                    pass
            if loaded == 0:
                raise ValueError(f"None of the .xml files in\n{inv_path}\n"
                                 "could be read as a valid inventory.")
            return inv

        if not os.path.isfile(inv_path):
            raise ValueError(f"Inventory file not found:\n{inv_path}")
        try:
            return read_inventory(inv_path)
        except Exception as e:
            raise ValueError(f"Could not read inventory:\n{inv_path}\n\n{e}")

    def _open_spectrogram(self):
        # ── 1. Require loaded data ────────────────────────────────────────
        if not self.stream:
            QMessageBox.warning(
                self, "No data",
                "Load a MiniSEED file before plotting a spectrogram.")
            return

        # Pressure calibration is optional; the viewer labels raw counts.
        inventory = None
        if self.inv_combo.currentIndex() >= 0:
            try:
                inventory = self._resolve_inventory()
            except ValueError as e:
                QMessageBox.warning(self, "Inventory unavailable", str(e))

        # ── 3. Time range — the FK duration, not the short FK window ──────
        # self.selected_time is a relative-seconds offset into the preview;
        # absolute start = stream[0].starttime + selected_time.
        stream = self.stream
        if self.selected_time is not None:
            try:
                t0 = self.stream[0].stats.starttime + float(self.selected_time)
                t1 = t0 + float(self.duration.value())
                windowed = self.stream.slice(t0, t1)
                if windowed and any(tr.stats.npts > 0 for tr in windowed):
                    stream = windowed
                else:
                    QMessageBox.warning(
                        self, "Empty window",
                        "The picked window contains no samples — "
                        "using the full loaded stream instead.")
            except Exception as e:
                QMessageBox.warning(
                    self, "Window error",
                    f"Could not slice the picked window ({e}).\n"
                    "Using the full loaded stream instead.")

        # ── 4. Seed parameters ────────────────────────────────────────────
        sample_counts = [tr.stats.npts for tr in stream if tr.stats.npts >= 32]
        if not sample_counts:
            QMessageBox.warning(
                self, "No spectrogram data",
                "The selected time range has too few samples. Pick an earlier "
                "start time or increase Duration in FK Configuration.")
            return
        nyquist = min(tr.stats.sampling_rate / 2 for tr in stream
                      if tr.stats.npts >= 32)
        band_high = min(self.max_freq.value(), nyquist * 0.9)
        band_low = min(self.min_freq.value(), band_high * 0.5)
        segment = min(512, max(32, int(np.median(sample_counts)) // 4))
        params = {
            "filter_freqmin": band_low,
            "filter_freqmax": band_high,
            "freq_max": min(self.max_freq.value(), nyquist),
            "nperseg": segment,
            "noverlap": int(segment * 0.75),
        }

        # ── 5. Open the window ────────────────────────────────────────────
        try:
            win = SpectrogramWindow(stream, inventory, params, parent=self)
        except Exception as e:
            QMessageBox.critical(
                self, "Spectrogram error",
                f"Could not build the spectrogram window:\n\n{e}")
            return
        win.exec_()

    def _show_readiness(self):
        from data_readiness import Finding, assess_stream

        inventory = None
        inventory_error = None
        if self.inv_combo.currentIndex() >= 0:
            try:
                # Unresolved, so the check can report overlapping epochs.
                inventory = self._resolve_inventory(resolve_epochs=False)
            except ValueError as error:
                inventory_error = str(error)
        findings = assess_stream(
            self.stream, inventory, fmin=self.min_freq.value(),
            fmax=self.max_freq.value(), start=self._window_start(),
            duration=self.duration.value())
        if inventory_error:
            findings.append(Finding("error", "StationXML could not be read",
                                    inventory_error))

        errors = sum(item.severity == "error" for item in findings)
        warnings = sum(item.severity == "warning" for item in findings)
        dialog = QDialog(self)
        apply_result_style(dialog)
        dialog.setWindowTitle("Data readiness · SeismoFK")
        dialog.resize(680, 520)
        column = QVBoxLayout(dialog)
        column.setContentsMargins(16, 16, 16, 12)
        column.setSpacing(10)
        verdict = ("Not ready — fix blocking issues first" if errors else
                   "Ready, with warnings" if warnings else "Ready for analysis")
        column.addWidget(result_header(
            verdict, f"{errors} blocking issue(s)  ·  {warnings} warning(s)",
            kicker="SEISMOFK / DATA READINESS",
            badge="BLOCKED" if errors else "READY"))
        colours = {"error": "#a33a2c", "warning": "#9a6a00"}
        rank = {"error": 0, "warning": 1}
        detail = QTextEdit()
        detail.setReadOnly(True)
        detail.setHtml("".join(
            f"<p><b style='color:{colours.get(item.severity, '#24455c')}'>"
            f"{html.escape(item.severity.upper())}</b> &nbsp;"
            f"<b>{html.escape(item.title)}</b><br>"
            f"<span style='color:#4a6474'>"
            f"{html.escape(item.detail).replace(chr(10), '<br>')}</span></p>"
            for item in sorted(findings, key=lambda f: rank.get(f.severity, 2))))
        column.addWidget(detail, stretch=1)
        close = QPushButton("Close")
        close.setObjectName("primaryAction")
        close.clicked.connect(dialog.accept)
        column.addWidget(close, alignment=Qt.AlignRight)
        dialog.exec_()

    # ── v1.2.1 advanced array methods ─────────────────────────────────────
    def _array_stream_with_coords(self):
        """Return a copy of the loaded stream with coordinates attached and a
        common sampling rate — shared by the ARF / slowness-map / PMCC dialogs.
        Raises ValueError with a user-friendly message on any problem."""
        if not self.stream:
            raise ValueError("Load a MiniSEED file before running this method.")
        inv = self._resolve_inventory()
        st = self.stream.copy()
        rates = {tr.stats.sampling_rate for tr in st}
        if len(rates) > 1:
            target = min(rates)
            for tr in st:
                if tr.stats.sampling_rate != target:
                    tr.resample(target)
        fk.attach_coordinates(st, inv, verbose=False)
        if len(st) < 2:
            raise ValueError("Array analysis needs at least 2 channels with "
                             "coordinates in the inventory.")
        # Physical units for waveform panels; direction and velocity do not
        # depend on this scaling.
        st.units = fk.scalar_calibrate(st, inv)
        return st

    def _window_start(self):
        """Absolute UTC start for single-window methods: the picked time if set,
        otherwise the Start Time field."""
        if self.selected_time is not None and self.stream:
            return self.stream[0].stats.starttime + float(self.selected_time)
        return UTCDateTime(self.start_time_input.dateTime()
                           .toString("yyyy-MM-dd HH:mm:ss"))

    def _open_array_response(self):
        from gui_methods import ArrayResponseDialog
        try:
            st = self._array_stream_with_coords()
        except ValueError as e:
            QMessageBox.warning(self, "Array Response", str(e))
            return
        ArrayResponseDialog(st, self.min_freq.value(), self.max_freq.value(),
                            parent=self).exec_()

    def _open_slowness_map(self):
        from gui_methods import SlownessMapDialog
        try:
            st = self._array_stream_with_coords()
        except ValueError as e:
            QMessageBox.warning(self, "Slowness Map", str(e))
            return
        exp_baz = None
        try:
            from obspy.geodetics import gps2dist_azimuth
            exp_baz = gps2dist_azimuth(
                self.event_lat.value(), self.event_lon.value(),
                st[0].stats.coordinates.latitude,
                st[0].stats.coordinates.longitude)[2]
        except Exception:
            pass
        SlownessMapDialog(st, self.min_freq.value(), self.max_freq.value(),
                          self._window_start(), self.win_len.value(),
                          expected_baz=exp_baz, parent=self).exec_()

    def _open_noise_levels(self):
        """Noise levels over the whole loaded record; calibrated to Pa with
        the StationXML sensitivities when an inventory is selected."""
        from gui_methods import NoiseLevelsDialog
        if not self.stream:
            return
        st = self.stream.copy()
        units = "counts"
        if self.inv_combo.currentIndex() >= 0:
            try:
                units = fk.scalar_calibrate(st, self._resolve_inventory())
            except ValueError as e:
                QMessageBox.warning(self, "Noise levels", str(e))
        NoiseLevelsDialog(st, self.min_freq.value(), self.max_freq.value(),
                          units=units, parent=self).exec_()

    def _open_pmcc(self):
        from gui_methods import PMCCWindow
        try:
            st = self._array_stream_with_coords()
        except ValueError as e:
            QMessageBox.warning(self, "PMCC", str(e))
            return
        PMCCWindow(st, self.min_freq.value(), self.max_freq.value(),
                   self._window_start(), self.duration.value(),
                   window_sec=self.win_len.value(),
                   units=getattr(st, "units", "counts"), parent=self).exec_()

    # ── inventory helpers ─────────────────────────────────────────────────
    def _set_start_time_utc(self, when):
        timestamp = QDateTime(when.year, when.month, when.day, when.hour,
                              when.minute, when.second,
                              int(when.microsecond / 1000))
        timestamp.setTimeSpec(Qt.UTC)
        self.start_time_input.setDateTime(timestamp)

    def _update_analysis_actions(self):
        ready = bool(self.stream) and self.inv_combo.currentIndex() >= 0
        for button in (self.process_btn, self.arf_btn, self.slowness_btn,
                       self.pmcc_btn):
            button.setEnabled(ready)
            button.setToolTip("" if ready else
                              "Load a MiniSEED waveform and select StationXML first.")
        self.spectro_btn.setEnabled(bool(self.stream))
        self.noise_btn.setEnabled(bool(self.stream))
        self.check_data_btn.setEnabled(bool(self.stream))
        self.spectro_btn.setToolTip(
            "Plot loaded waveforms; without pressure metadata, values use raw counts."
            if self.stream else "Load a MiniSEED waveform first.")
        if not ready:
            self.status_label.setText("Load a waveform and select StationXML")
        elif self.process_thread is None or not self.process_thread.isRunning():
            self.status_label.setText("Ready to analyze")

    def refresh_inventory_list(self):
        self.inv_combo.clear()
        base_dir = os.path.dirname(os.path.abspath(__file__))
        xml_dir    = os.path.join(base_dir, 'XML')
        xml_im_dir = os.path.join(base_dir, 'XML_IM')
        user_xml_dir = inventory_dir()

        model = QStandardItemModel()
        self.inv_combo.setModel(model)

        # ── All IMS stations shortcut ──────────────────────────────────────
        if os.path.isdir(xml_im_dir) and any(
                f.endswith('.xml') for f in os.listdir(xml_im_dir)):
            item = QStandardItem("★ All IMS Stations (XML_IM/)")
            item.setData(xml_im_dir, Qt.UserRole)   # directory path → load all
            model.appendRow(item)

        # ── Individual XML_IM/ files ───────────────────────────────────────
        if os.path.isdir(xml_im_dir):
            for f in sorted(f for f in os.listdir(xml_im_dir) if f.endswith('.xml')):
                item = QStandardItem(f"  {f.replace('.xml', '')}  [IMS]")
                item.setData(os.path.join(xml_im_dir, f), Qt.UserRole)
                model.appendRow(item)

        # ── Legacy XML/ files ──────────────────────────────────────────────
        if os.path.isdir(xml_dir):
            for f in sorted(f for f in os.listdir(xml_dir) if f.endswith('.xml')):
                display = f.replace('.txt.xml', '').replace('.xml', '')
                item = QStandardItem(display)
                item.setData(os.path.join(xml_dir, f), Qt.UserRole)
                model.appendRow(item)
        if user_xml_dir.is_dir() and str(user_xml_dir.resolve()) != os.path.realpath(xml_dir):
            for path in sorted(user_xml_dir.glob('*.xml')):
                item = QStandardItem(f"{path.stem}  [User]")
                item.setData(str(path), Qt.UserRole)
                model.appendRow(item)
        self.inv_combo.setCurrentIndex(-1)

    def add_inventory_file(self):
        fname, _ = QFileDialog.getOpenFileName(
            self, "Select XML Inventory", "", "XML Files (*.xml);;All Files (*)")
        if not fname:
            return
        xml_dir = inventory_dir(create=True)
        dest = str(xml_dir / os.path.basename(fname))
        if os.path.exists(dest):
            if QMessageBox.question(
                    self, "Overwrite?",
                    f"{os.path.basename(fname)} already exists. Overwrite?",
                    QMessageBox.Yes | QMessageBox.No) == QMessageBox.No:
                return
        shutil.copy2(fname, dest)
        self.refresh_inventory_list()
        QMessageBox.information(self, "Added",
                                f"Inventory added:\n{os.path.basename(fname)}")

    # ── MiniSEED loading ──────────────────────────────────────────────────
    def load_mseed(self):
        fnames, _ = QFileDialog.getOpenFileNames(
            self, "Select MiniSEED file(s)", "",
            "MiniSEED (*.mseed *.msd *.ms);;All Files (*)")
        if not fnames:
            return
        try:
            from obspy import Stream
            combined = Stream()
            for f in sorted(fnames):
                combined += read(f)
            combined.merge(fill_value=0)
            for tr in combined:
                if isinstance(tr.data, np.ma.MaskedArray):
                    tr.data = tr.data.filled(0)
            self.stream = combined
            self.selected_time = None
            self.time_label.setText("Pick: —")
            self.clear_pick_btn.setEnabled(False)
            first = self.stream[0].stats.starttime
            self._set_start_time_utc(first)
            if len(fnames) == 1:
                self.mseed_edit.setText(fnames[0])
            else:
                self.mseed_edit.setText(
                    f"{len(fnames)} files merged  [{', '.join(os.path.basename(f) for f in fnames)}]")
            self._update_preview()
            self._update_analysis_actions()
        except Exception as e:
            QMessageBox.critical(self, "Load Error", str(e))

    # ── Waveform preview ──────────────────────────────────────────────────
    def _update_preview(self):
        """Draw each trace as a min/max envelope (instant even for long
        records) on a shared UTC time axis, labelled inside its panel."""
        if not self.stream:
            return
        t0 = self.stream[0].stats.starttime
        self.figure.clear()
        axes = self.figure.subplots(len(self.stream), 1, sharex=True,
                                    squeeze=False)[:, 0]
        self._preview_axes = list(axes)
        self._window_artists = []
        for ax, tr in zip(axes, self.stream):
            x, y = _envelope(tr.data, tr.stats.delta,
                             float(tr.stats.starttime - t0))
            ax.plot(x, y / (np.max(np.abs(y)) or 1.0),
                    color=PALETTE['ink_2'], lw=0.55)
            ax.set_ylim(-1.12, 1.12)
            ax.set_yticks([])
            for side in ('top', 'right', 'left'):
                ax.spines[side].set_visible(False)
            ax.spines['bottom'].set_color('#e6ebef')
            ax.tick_params(axis='x', length=0)
            ax.grid(axis='x', color='#eef2f5', lw=0.6)
            ax.text(0.004, 0.9, f"{tr.stats.station}  {tr.stats.channel}",
                    transform=ax.transAxes, va='top', fontsize=8,
                    fontweight='semibold', color=PALETTE['ink_2'],
                    bbox=dict(boxstyle='round,pad=0.2', fc='white',
                              ec='none', alpha=0.85))
        last = axes[-1]
        last.spines['bottom'].set_color('#a9b7c2')
        last.tick_params(axis='x', length=3, labelsize=8)
        last.xaxis.set_major_formatter(
            FuncFormatter(lambda x, _: (t0 + x).strftime('%H:%M:%S')))
        last.set_xlabel(
            f"UTC  ·  {t0.strftime('%Y-%m-%d')}  ·  {len(self.stream)} traces"
            f"  ·  {self.stream[0].stats.sampling_rate:g} Hz", fontsize=8)
        self.figure.tight_layout(pad=0.4, h_pad=0.05)
        self._draw_analysis_window()
        self.view_window_btn.setEnabled(True)

    def _draw_analysis_window(self, *_):
        """Shade [start, start + Duration] on every trace; follows the start
        time and Duration fields as they change."""
        if not (self.stream and self._preview_axes):
            return
        for artist in self._window_artists:
            artist.remove()
        self._window_artists = []
        t0 = self.stream[0].stats.starttime
        start = float(self._window_start() - t0)
        xlim = self._preview_axes[0].get_xlim()      # keep the user's zoom
        for ax in self._preview_axes:
            self._window_artists += [
                ax.axvspan(start, start + self.duration.value(),
                           color=PALETTE['accent'], alpha=0.08, lw=0, zorder=0),
                ax.axvline(start, color=PALETTE['accent'], lw=1.1, zorder=3)]
        self._preview_axes[0].set_xlim(xlim)
        self.canvas.draw_idle()

    def _on_preview_click(self, event):
        if not (event.inaxes and event.button == 1):
            return
        self.selected_time = event.xdata
        # Always compute absolute time from first trace start
        pt = self.stream[0].stats.starttime + self.selected_time
        self.time_label.setText(
            f"Pick: {pt.strftime('%Y-%m-%d %H:%M:%S.%f')[:-4]}")
        self._set_start_time_utc(pt)
        self.clear_pick_btn.setEnabled(True)
        self._draw_analysis_window()

    def clear_time_pick(self):
        self.selected_time = None
        self.time_label.setText("Pick: —")
        self.clear_pick_btn.setEnabled(False)
        if self.stream:
            first = self.stream[0].stats.starttime
            self._set_start_time_utc(first)
        self._draw_analysis_window()

    def view_waveform(self):
        if not self.stream:
            return
        viewer = WaveformViewer(self.stream, self)
        viewer.time_selected.connect(self._apply_viewer_pick)
        viewer.exec_()

    def _apply_viewer_pick(self, time_str):
        try:
            t = UTCDateTime(time_str)
            self._set_start_time_utc(t)
        except Exception as e:
            QMessageBox.critical(self, "Time Error", str(e))

    # ── Processing ────────────────────────────────────────────────────────
    def process_data(self):
        if not self.mseed_edit.text():
            QMessageBox.warning(self, "Missing input", "Please load a MiniSEED file.")
            return
        if self.inv_combo.currentIndex() < 0:
            QMessageBox.warning(self, "Missing input", "Please select an inventory file.")
            return

        item = self.inv_combo.model().item(self.inv_combo.currentIndex())
        self.inventory_file = item.data(Qt.UserRole)
        if not self.inventory_file:
            QMessageBox.warning(self, "Error", "Invalid inventory selection.")
            return

        self.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        self.status_label.setText("Starting ...")

        params = dict(
            mseed_file     = self.mseed_edit.text(),
            stream         = self.stream,           # pass pre-merged stream
            inventory_file = self.inventory_file,
            min_freq       = self.min_freq.value(),
            max_freq       = self.max_freq.value(),
            win_length     = self.win_len.value(),
            overlap        = self.overlap.value(),
            start_time     = self.start_time_input.dateTime().toString(
                                 "yyyy-MM-dd HH:mm:ss"),
            duration       = self.duration.value(),
            event_name     = self.event_name.text(),
            event_lat      = self.event_lat.value(),
            event_lon      = self.event_lon.value(),
            semb_thresh    = self.semb_thresh.value(),
            origin_time    = (self.origin_time_input.dateTime()
                              .toString("yyyy-MM-dd HH:mm:ss")
                              if self.origin_time_chk.isChecked() else None),
            celerity       = (self.celerity_spin.value()
                              if self.celerity_chk.isChecked() else None),
            bootstrap      = self.bootstrap_chk.isChecked(),
        )

        self._last_params = params          # kept for DB save dialog

        if self.process_thread:
            self.process_thread.quit()
            self.process_thread.wait()
        self.process_thread = ProcessThread(params)
        self.process_thread.progress.connect(self._on_progress)
        self.process_thread.finished.connect(self._on_finished)
        self.process_thread.error.connect(self._on_error)
        self.process_thread.status.connect(self.status_label.setText)
        self.process_thread.start()

    def _on_progress(self, val):
        self.progress_bar.setValue(val)

    def _on_finished(self, result):
        self.progress_bar.setVisible(False)
        self.setEnabled(True)
        self._show_results(result, self._last_params)

    def _on_error(self, msg):
        self.setEnabled(True)
        self.progress_bar.setVisible(False)
        QMessageBox.critical(self, "Processing Error", msg)

    # ── Results figure — 6-panel 2-column layout ──────────────────────────
    def _show_results(self, r, params=None):
        """
        Left column  (4 rows): Fisher | BAZ | App. Velocity | Beam waveform
        Right column (2 spans): FK slowness map | Array geometry
        """
        if params is None:
            params = getattr(self, '_last_params', {})
        semb_thresh = r.get('semb_thresh', self.semb_thresh.value())
        event_name  = self.event_name.text()
        fmin        = self.min_freq.value()
        fmax        = self.max_freq.value()

        # Figure() rather than plt.figure() so results windows are not kept
        # alive by pyplot; the same figure is drawn by the batch CLI.
        fig = Figure(figsize=(15, 11))
        result_plots.fk_results_figure(
            fig, r, params, event_name=event_name, fmin=fmin, fmax=fmax,
            semb_thresh=semb_thresh)

        # Auto-save JPG ────────────────────────────────────────────────────
        safe_name = "".join(c if c.isalnum() or c in '-_' else '_'
                            for c in event_name)
        jpg_file  = f"FK_{safe_name}.jpg"
        try:
            # matplotlib ≥3.3 moved JPEG quality into pil_kwargs (the old
            # `quality=` kwarg was removed and raised on newer versions).
            fig.savefig(jpg_file, dpi=300, bbox_inches='tight',
                        format='jpeg', pil_kwargs={'quality': 92})
            print(f"[INFO] Figure saved → {jpg_file}")
        except Exception as e:
            print(f"[WARN] Could not auto-save figure: {e}")

        # Attach figure path to result so SaveEventDialog can reference it
        r['figure_path'] = os.path.abspath(jpg_file)

        win = ResultsWindow(fig, result=r, params=params, parent=self)
        win.show()


# ═══════════════════════════════════════════════════════════════════════════
def main():
    app = QApplication(sys.argv)
    apply_app_theme(app)
    win = FKAnalysisGUI()
    win.show()
    return app.exec_()


if __name__ == '__main__':
    sys.exit(main())
