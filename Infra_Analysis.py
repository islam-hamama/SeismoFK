"""
Infra_Analysis.py — SeismoFK Infrasound FK Array Analysis GUI

Copyright (c) 2024-2025 Islam Hamama
Contact: islam.hamama@nriag.sci.eg

Licensed under the MIT License — see LICENSE for details.
"""
import sys
import os
import shutil

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
    QComboBox, QStyledItemDelegate, QGroupBox, QSizePolicy,
    QTextEdit, QTableWidget, QTableWidgetItem, QHeaderView, QCheckBox,
)
from PyQt5.QtCore  import Qt, QDateTime, QThread, pyqtSignal
from PyQt5.QtGui   import QPixmap, QStandardItem, QStandardItemModel

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)

from obspy import read, read_inventory, UTCDateTime
from fk_analysis import fk_array
import db_manager


# ═══════════════════════════════════════════════════════════════════════════
#  Background worker thread
# ═══════════════════════════════════════════════════════════════════════════
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
            self.progress.emit(35)

            # ── Step 3: Merge and clean masked arrays ─────────────────────
            st.merge(fill_value=0)
            for tr in st:
                if isinstance(tr.data, np.ma.MaskedArray):
                    tr.data = tr.data.filled(0)
            self.progress.emit(45)

            # ── Step 4: Full instrument response removal → output in Pa ──────
            self.status.emit("Removing instrument response ...")
            try:
                st.remove_response(
                    inventory=inv,
                    output='VEL',
                    pre_filt=(0.1, 0.5, 9.0, 10.0),
                    water_level=60,
                )
            except Exception as e:
                # Fallback: scalar sensitivity division if response removal fails
                self.status.emit(f"Response removal failed ({e}), using scalar sensitivity ...")
                for tr in st:
                    try:
                        resp = inv.get_response(tr.id, datetime=tr.stats.starttime)
                        sens = resp.instrument_sensitivity.value
                        if sens and sens != 0:
                            tr.data = tr.data.astype(float) / sens
                    except Exception:
                        pass   # keep raw counts if no response found
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

        logo_label = QLabel()
        logo_label.setStyleSheet(
            "border-radius:8px; background:rgba(0,0,0,0.05); padding:4px;")
        logo_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logo.png")
        if os.path.exists(logo_path):
            pix = QPixmap(logo_path).scaled(28, 28, Qt.KeepAspectRatio,
                                             Qt.SmoothTransformation)
            logo_label.setPixmap(pix)
        layout.addWidget(logo_label)
        layout.addStretch()

        copy_label = QLabel("© 2024-2025 Islam Hamama  |  islam.hamama@nriag.sci.eg")
        copy_label.setStyleSheet(
            "font-weight:bold; color:#666; padding:4px;"
            "background:rgba(0,0,0,0.05); border-radius:4px;")
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

    _STYLE = """
        QGroupBox { font-weight:bold; border:1px solid #ccc;
                    border-radius:6px; margin-top:8px; padding-top:6px; }
        QGroupBox::title { subcontrol-origin:margin; left:10px; color:#2c3e50; }
        QTextEdit { border:1px solid #bdc3c7; border-radius:4px;
                    background:white; padding:4px; }
        QComboBox { padding:5px; border:1px solid #bdc3c7;
                    border-radius:4px; background:white; }
    """

    def __init__(self, result: dict, params: dict, fig=None, parent=None):
        super().__init__(parent)
        self.result = result
        self.params = params
        self.fig    = fig
        self.setWindowTitle("Save Event to Database")
        self.setMinimumWidth(520)
        self.setStyleSheet(self._STYLE)

        layout = QVBoxLayout(self)
        layout.setSpacing(10)

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
        lbl.setStyleSheet(
            "font-family:monospace; font-size:9pt; color:#2c3e50;"
            "background:#f0f4f8; border-radius:4px; padding:8px;")
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
        save_btn = QPushButton("💾  Save")
        save_btn.setStyleSheet(
            "QPushButton{background:#27ae60;color:white;font-weight:bold;"
            "padding:8px 20px;border-radius:5px;}"
            "QPushButton:hover{background:#1e8449;}")
        save_btn.clicked.connect(self._save)
        cancel_btn = QPushButton("Cancel")
        cancel_btn.setStyleSheet(
            "QPushButton{background:#95a5a6;color:white;padding:8px 20px;"
            "border-radius:5px;} QPushButton:hover{background:#7f8c8d;}")
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
        self.setWindowTitle("SeismoFK — Event Database")
        self.setGeometry(120, 80, 1150, 620)
        self._all_rows = []
        self._build_ui()
        self._load()

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # ── Filter bar ────────────────────────────────────────────────────
        filter_row = QHBoxLayout()
        filter_row.addWidget(QLabel("Filter by classification:"))
        self.filter_combo = QComboBox()
        self.filter_combo.addItem("All")
        self.filter_combo.addItems(db_manager.CLASSIFICATIONS)
        self.filter_combo.currentTextChanged.connect(self._apply_filter)
        filter_row.addWidget(self.filter_combo)
        filter_row.addStretch()

        refresh_btn = QPushButton("⟳ Refresh")
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
        self.table.setStyleSheet(
            "QTableWidget { font-size:9pt; }"
            "QHeaderView::section { background:#2c3e50; color:white;"
            "  font-weight:bold; padding:4px; }")
        layout.addWidget(self.table)

        # ── Row count label ───────────────────────────────────────────────
        self.count_lbl = QLabel("")
        self.count_lbl.setStyleSheet("color:#555; font-style:italic;")
        layout.addWidget(self.count_lbl)

        # ── Action buttons ────────────────────────────────────────────────
        btn_row = QHBoxLayout()

        edit_btn = QPushButton("✏  Edit Classification / Note")
        edit_btn.clicked.connect(self._edit_row)
        view_fig_btn = QPushButton("🖼  View Figure")
        view_fig_btn.setStyleSheet(
            "QPushButton{background:#2980b9;color:white;border-radius:5px;"
            "padding:6px 14px;} QPushButton:hover{background:#1f618d;}")
        view_fig_btn.clicked.connect(self._view_figure)
        delete_btn = QPushButton("🗑  Delete Selected")
        delete_btn.setStyleSheet(
            "QPushButton{background:#e74c3c;color:white;border-radius:5px;"
            "padding:6px 14px;} QPushButton:hover{background:#c0392b;}")
        delete_btn.clicked.connect(self._delete_row)
        export_btn = QPushButton("📤  Export to CSV")
        export_btn.clicked.connect(self._export)
        close_btn  = QPushButton("✖  Close")
        close_btn.clicked.connect(self.close)

        for b in (edit_btn, view_fig_btn, delete_btn, export_btn):
            btn_row.addWidget(b)
        btn_row.addStretch()
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
        title = (f"Figure — {row_data.get('event_name','?')}  "
                 f"[{row_data.get('saved_at','?')}]")
        dlg.setWindowTitle(title)
        dlg.setMinimumSize(900, 650)
        vlay = QVBoxLayout(dlg)

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
        save_btn = QPushButton("💾  Save as PNG")
        save_btn.clicked.connect(lambda: self._save_blob_png(blob, row_data))
        close_btn = QPushButton("✖  Close")
        close_btn.clicked.connect(dlg.close)
        btn_row.addWidget(save_btn)
        btn_row.addStretch()
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
        dlg.setWindowTitle(f"Edit Event #{event_id}")
        dlg.setMinimumWidth(420)
        vlay = QVBoxLayout(dlg)

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
        ok_btn  = QPushButton("✔  Update")
        ok_btn.setStyleSheet(
            "QPushButton{background:#27ae60;color:white;border-radius:5px;"
            "padding:6px 16px;} QPushButton:hover{background:#1e8449;}")
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
        self.setWindowTitle("SeismoFK — Analysis Results")
        self.setGeometry(150, 100, 1300, 1050)
        self.fig    = fig
        self.result = result or {}
        self.params = params or {}

        layout = QVBoxLayout(self)
        self.canvas  = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        btn_row = QHBoxLayout()
        save_fig_btn = QPushButton("💾  Save Figure (300 DPI)")
        save_fig_btn.clicked.connect(self.save_figure)

        save_db_btn = QPushButton("🗄  Save to Database")
        save_db_btn.setStyleSheet(
            "QPushButton{background:#8e44ad;color:white;font-weight:bold;"
            "padding:7px 16px;border-radius:5px;}"
            "QPushButton:hover{background:#6c3483;}")
        save_db_btn.clicked.connect(self._save_to_db)

        close_btn = QPushButton("✖  Close")
        close_btn.clicked.connect(self.close)

        btn_row.addWidget(save_fig_btn)
        btn_row.addWidget(save_db_btn)
        btn_row.addStretch()
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
        ``noverlap``, ``reference_pressure``.
    """

    # Fixed defaults mirror plot_spectrogram_window.py's standalone values.
    _DEFAULTS = {
        "filter_freqmin": 0.5,
        "filter_freqmax": 6.0,
        "freq_max": 5.0,
        "nperseg": 512,
        "noverlap": 460,
        "reference_pressure": 20e-6,
    }

    def __init__(self, stream, inventory, params=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle("SeismoFK — Spectrograms")
        self.setGeometry(140, 90, 1300, 1000)

        self.stream    = stream
        self.inventory = inventory
        cfg = dict(self._DEFAULTS)
        if params:
            cfg.update({k: v for k, v in params.items() if v is not None})
        self.reference_pressure = float(cfg["reference_pressure"])

        layout = QVBoxLayout(self)

        # ── In-window control row (Option C) ──────────────────────────────
        ctrl_group  = QGroupBox("Spectrogram Parameters")
        ctrl_layout = QHBoxLayout()
        self.freqmin_input  = QLineEdit(f"{cfg['filter_freqmin']:g}")
        self.freqmax_input  = QLineEdit(f"{cfg['filter_freqmax']:g}")
        self.freq_max_input = QLineEdit(f"{cfg['freq_max']:g}")
        self.nperseg_input  = QLineEdit(str(int(cfg["nperseg"])))
        self.noverlap_input = QLineEdit(str(int(cfg["noverlap"])))
        for w in (self.freqmin_input, self.freqmax_input, self.freq_max_input,
                  self.nperseg_input, self.noverlap_input):
            w.setMaximumWidth(80)
        recompute_btn = QPushButton("↻  Recompute")
        recompute_btn.setStyleSheet(
            "QPushButton{background:#2980b9;color:white;font-weight:bold;"
            "padding:6px 14px;border-radius:5px;}"
            "QPushButton:hover{background:#1f618d;}")
        recompute_btn.clicked.connect(self._recompute)
        for w in (QLabel("Bandpass min (Hz):"), self.freqmin_input,
                  QLabel("max (Hz):"),          self.freqmax_input,
                  QLabel("Freq. max (Hz):"),    self.freq_max_input,
                  QLabel("nperseg:"),           self.nperseg_input,
                  QLabel("noverlap:"),          self.noverlap_input,
                  recompute_btn):
            ctrl_layout.addWidget(w)
        ctrl_layout.addStretch()
        ctrl_group.setLayout(ctrl_layout)
        layout.addWidget(ctrl_group)

        # ── Canvas + toolbar ──────────────────────────────────────────────
        # Use Figure() directly (not plt.figure()) so the figure is NOT
        # registered in pyplot's global figure manager — otherwise every
        # "Plot Spectrogram" invocation would orphan and leak a figure for
        # the process lifetime. This is the documented matplotlib-in-Qt
        # pattern; the canvas/toolbar own the figure's lifecycle.
        self.fig    = Figure(figsize=(13, 9), facecolor="white")
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas, stretch=1)

        # ── Action buttons ────────────────────────────────────────────────
        btn_row = QHBoxLayout()
        save_fig_btn = QPushButton("💾  Save Figure (300 DPI)")
        save_fig_btn.clicked.connect(self.save_figure)
        close_btn = QPushButton("✖  Close")
        close_btn.clicked.connect(self.close)
        btn_row.addWidget(save_fig_btn)
        btn_row.addStretch()
        btn_row.addWidget(close_btn)
        layout.addLayout(btn_row)
        layout.addWidget(FooterWidget(self))

        # Initial render using the seeded parameters.
        self._render(cfg["filter_freqmin"], cfg["filter_freqmax"],
                     cfg["freq_max"], int(cfg["nperseg"]), int(cfg["noverlap"]))

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
        return freqmin, freqmax, freq_max, nperseg, noverlap

    def _recompute(self):
        try:
            params = self._read_params()
        except ValueError as e:
            QMessageBox.warning(self, "Invalid parameters", str(e))
            return
        self._render(*params)

    # ── Rendering ─────────────────────────────────────────────────────────
    def _render(self, filter_freqmin, filter_freqmax, freq_max,
                nperseg, noverlap):
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
                )
            except Exception as e:               # noqa: BLE001 — per-trace isolation
                skipped.append(f"{trace.id}: {e}")
                continue
            panels.append((trace, freqs, time_nums, power_db))

        if not panels:
            self.canvas.draw()
            QMessageBox.critical(
                self, "Spectrogram failed",
                "No trace could be processed.\n\n" + "\n".join(skipped))
            return

        vmin = min(np.percentile(p[3], 15) for p in panels)
        vmax = max(np.percentile(p[3], 98) for p in panels)

        gs = gridspec.GridSpec(
            len(panels), 2, width_ratios=[40, 1.8],
            hspace=0.14, wspace=0.08,
            left=0.08, right=0.91, top=0.92, bottom=0.08, figure=self.fig)
        cax = self.fig.add_subplot(gs[:, 1])

        axes = []
        last_mesh = None
        for idx, (trace, freqs, time_nums, power_db) in enumerate(panels):
            ax = self.fig.add_subplot(
                gs[idx, 0], sharex=axes[0] if axes else None)
            axes.append(ax)
            last_mesh = ax.pcolormesh(
                time_nums, freqs, power_db,
                shading="auto", cmap="magma", vmin=vmin, vmax=vmax)
            ax.set_ylabel(f"{trace.stats.station}\nHz", rotation=0, labelpad=28)
            ax.set_ylim(0, freq_max)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.grid(True, color="white", linewidth=0.6, alpha=0.25)
            ax.text(
                0.015, 0.92, trace.id, transform=ax.transAxes,
                ha="left", va="top", fontsize=9, color="white",
                bbox=dict(facecolor="black", alpha=0.28,
                          edgecolor="none", pad=3.0))

        ref_label = f"{self.reference_pressure * 1e6:g} µPa"
        self.fig.suptitle(
            f"{panels[0][0].stats.station} Spectrograms\n"
            f"PSD referenced to {ref_label}  ·  "
            f"bandpass {filter_freqmin:g}-{filter_freqmax:g} Hz",
            fontsize=14, fontweight="semibold", y=0.975)

        axes[-1].set_xlabel("Time (UTC)")
        locator   = mdates.AutoDateLocator()
        formatter = mdates.DateFormatter("%H:%M:%S")
        axes[-1].xaxis.set_major_locator(locator)
        axes[-1].xaxis.set_major_formatter(formatter)
        for ax in axes[:-1]:
            plt.setp(ax.get_xticklabels(), visible=False)

        if last_mesh is not None:
            cbar = self.fig.colorbar(last_mesh, cax=cax)
            cbar.set_label("PSD (dB re 20 µPa²/Hz)")
            cbar.outline.set_linewidth(0.8)

        self.canvas.draw()

        if skipped:
            QMessageBox.warning(
                self, "Some traces skipped",
                "These traces could not be processed and were omitted:\n\n"
                + "\n".join(skipped))

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
        self.setWindowTitle("SeismoFK — Waveform Viewer")
        self.setGeometry(80, 80, 1280, 650)
        self.stream          = stream
        self.original_stream = stream.copy()
        self.selected_time   = None

        layout = QVBoxLayout(self)

        # ── Filter controls ───────────────────────────────────────────────
        filter_group  = QGroupBox("Band-pass Filter")
        filter_layout = QHBoxLayout()
        self.low_freq_input  = QLineEdit("0.5")
        self.high_freq_input = QLineEdit("5.0")
        apply_btn = QPushButton("Apply Filter")
        reset_btn = QPushButton("Reset")
        apply_btn.clicked.connect(self.apply_filter)
        reset_btn.clicked.connect(self.reset_filter)
        for w in (QLabel("Low (Hz):"), self.low_freq_input,
                  QLabel("High (Hz):"), self.high_freq_input,
                  apply_btn, reset_btn):
            filter_layout.addWidget(w)
        filter_group.setLayout(filter_layout)
        layout.addWidget(filter_group)

        instr = QLabel("⬤  Left-click on waveform to pick analysis start time")
        instr.setStyleSheet("color:#2980b9; font-weight:bold; padding:4px;")
        layout.addWidget(instr)

        self.fig_wv, self.ax_wv = plt.subplots(figsize=(13, 5))
        self.fig_wv.patch.set_facecolor('white')
        self.canvas_wv  = FigureCanvas(self.fig_wv)
        self.toolbar_wv = NavigationToolbar(self.canvas_wv, self)
        layout.addWidget(self.toolbar_wv)
        layout.addWidget(self.canvas_wv)

        self._draw_waveforms(self.stream)
        self.canvas_wv.mpl_connect('button_press_event', self._on_click)

        btn_row = QHBoxLayout()
        self.time_label = QLabel("Selected time: —")
        self.time_label.setStyleSheet("font-weight:bold; color:#c0392b;")
        confirm_btn = QPushButton("✔  Confirm")
        confirm_btn.clicked.connect(self.confirm_selection)
        cancel_btn  = QPushButton("✖  Cancel")
        cancel_btn.clicked.connect(self.reject)
        btn_row.addWidget(self.time_label, stretch=1)
        btn_row.addWidget(confirm_btn)
        btn_row.addWidget(cancel_btn)
        layout.addLayout(btn_row)
        layout.addWidget(FooterWidget(self))

    def _draw_waveforms(self, stream, title="Waveforms"):
        self.ax_wv.clear()
        for i, tr in enumerate(stream):
            t    = tr.times()
            norm = np.max(np.abs(tr.data)) or 1
            self.ax_wv.plot(t, tr.data / norm + i, lw=0.7, label=tr.id)
        self.ax_wv.set_xlabel("Time (s)")
        self.ax_wv.set_ylabel("Norm. Amplitude + offset")
        self.ax_wv.set_title(title)
        self.ax_wv.legend(fontsize=7, loc='upper right')
        self.ax_wv.grid(True, alpha=0.3)
        if self.selected_time is not None:
            rel = self.selected_time - stream[0].stats.starttime
            self.ax_wv.axvline(rel, color='red', lw=1.4, ls='--', alpha=0.85)
        self.canvas_wv.draw()

    def _on_click(self, event):
        if event.inaxes != self.ax_wv or event.button != 1:
            return
        self.selected_time = self.stream[0].stats.starttime + event.xdata
        self.time_label.setText(f"Selected time: {self.selected_time}")
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

    _STYLE = """
        QMainWindow, QDialog { background:#f4f6f8; }
        QWidget              { font-family:'Segoe UI',Arial,sans-serif;
                               font-size:10pt; }
        QGroupBox            { font-weight:bold; border:1px solid #ccc;
                               border-radius:6px; margin-top:8px; padding-top:6px; }
        QGroupBox::title     { subcontrol-origin:margin; left:10px; color:#2c3e50; }
        QPushButton          { background:#2980b9; color:white; border:none;
                               padding:7px 16px; border-radius:5px;
                               font-weight:bold; min-width:80px; }
        QPushButton:hover    { background:#1f618d; }
        QPushButton:disabled { background:#bdc3c7; color:#7f8c8d; }
        QLineEdit, QComboBox, QDoubleSpinBox, QDateTimeEdit {
                               padding:5px; border:1px solid #bdc3c7;
                               border-radius:4px; background:white; }
        QLabel               { color:#2c3e50; }
        QProgressBar         { border:2px solid #bdc3c7; border-radius:8px;
                               text-align:center; height:18px; }
    """

    def __init__(self):
        super().__init__()
        self.setStyleSheet(self._STYLE)
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
        layout.setSpacing(8)
        self.setCentralWidget(root)

        # ── 1. File selection ─────────────────────────────────────────────
        file_grp    = QGroupBox("File Selection")
        file_layout = QVBoxLayout()

        mseed_row = QHBoxLayout()
        self.mseed_edit = QLineEdit()
        self.mseed_edit.setPlaceholderText("Select a MiniSEED file …")
        browse_mseed = QPushButton("Browse")
        browse_mseed.clicked.connect(self.load_mseed)
        mseed_row.addWidget(QLabel("MiniSEED:"))
        mseed_row.addWidget(self.mseed_edit, stretch=1)
        mseed_row.addWidget(browse_mseed)
        file_layout.addLayout(mseed_row)

        inv_row = QHBoxLayout()
        self.inv_combo = QComboBox()
        self.inv_combo.setMinimumWidth(280)
        self.inv_combo.setItemDelegate(QStyledItemDelegate())
        refresh_btn = QPushButton("⟳ Refresh")
        refresh_btn.clicked.connect(self.refresh_inventory_list)
        add_xml_btn = QPushButton("+ Add XML")
        add_xml_btn.clicked.connect(self.add_inventory_file)
        inv_row.addWidget(QLabel("Inventory:"))
        inv_row.addWidget(self.inv_combo, stretch=1)
        inv_row.addWidget(refresh_btn)
        inv_row.addWidget(add_xml_btn)
        file_layout.addLayout(inv_row)
        file_grp.setLayout(file_layout)
        layout.addWidget(file_grp)

        # ── 2. Waveform preview ───────────────────────────────────────────
        wave_grp    = QGroupBox("Waveform Preview  (click to pick analysis start)")
        wave_layout = QVBoxLayout()

        self.figure = plt.figure(figsize=(11, 3))
        self.figure.patch.set_facecolor('white')
        self.canvas  = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.toolbar = NavigationToolbar(self.canvas, self)
        wave_layout.addWidget(self.toolbar)
        wave_layout.addWidget(self.canvas)

        pick_row = QHBoxLayout()
        self.time_label = QLabel("Pick: —")
        self.time_label.setStyleSheet("font-weight:bold; color:#c0392b;")
        self.clear_pick_btn = QPushButton("✖ Clear Pick")
        self.clear_pick_btn.clicked.connect(self.clear_time_pick)
        self.clear_pick_btn.setEnabled(False)
        self.view_window_btn = QPushButton("⤢ Open in Window")
        self.view_window_btn.clicked.connect(self.view_waveform)
        self.view_window_btn.setEnabled(False)
        pick_row.addWidget(self.time_label, stretch=1)
        pick_row.addWidget(self.clear_pick_btn)
        pick_row.addWidget(self.view_window_btn)
        wave_layout.addLayout(pick_row)
        wave_grp.setLayout(wave_layout)
        layout.addWidget(wave_grp, stretch=1)

        # ── 3. Parameters ─────────────────────────────────────────────────
        param_grp    = QGroupBox("Analysis Parameters")
        param_layout = QVBoxLayout()

        row1 = QHBoxLayout()
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
        for label, widget in [
            ("Min Freq (Hz):", self.min_freq),
            ("Max Freq (Hz):", self.max_freq),
            ("Window (s):",    self.win_len),
            ("Overlap:",       self.overlap),
            ("Semb. Threshold:", self.semb_thresh),
        ]:
            row1.addWidget(QLabel(label)); row1.addWidget(widget)
        param_layout.addLayout(row1)

        row2 = QHBoxLayout()
        self.start_time_input = QDateTimeEdit(QDateTime.currentDateTime())
        self.start_time_input.setDisplayFormat("yyyy-MM-dd HH:mm:ss")
        self.duration = QDoubleSpinBox()
        self.duration.setRange(1, 86400); self.duration.setValue(900)
        self.duration.setSuffix(" s")
        self.event_name = QLineEdit("Event")
        self.event_lat  = QDoubleSpinBox()
        self.event_lat.setRange(-90, 90); self.event_lat.setDecimals(4)
        self.event_lon  = QDoubleSpinBox()
        self.event_lon.setRange(-180, 180); self.event_lon.setDecimals(4)
        for label, widget in [
            ("Start Time:", self.start_time_input),
            ("Duration:",   self.duration),
            ("Event Name:", self.event_name),
            ("Lat:",        self.event_lat),
            ("Lon:",        self.event_lon),
        ]:
            row2.addWidget(QLabel(label)); row2.addWidget(widget)
        param_layout.addLayout(row2)
        # ── row3 — optional origin time & celerity ────────────────────────
        row3 = QHBoxLayout()

        self.origin_time_chk = QCheckBox("Origin Time:")
        self.origin_time_chk.setToolTip(
            "Known event origin time — used together with Celerity\n"
            "to draw the expected infrasound arrival on the beam.")
        self.origin_time_input = QDateTimeEdit(QDateTime.currentDateTime())
        self.origin_time_input.setDisplayFormat("yyyy-MM-dd HH:mm:ss")
        self.origin_time_input.setEnabled(False)
        self.origin_time_chk.toggled.connect(self.origin_time_input.setEnabled)

        self.celerity_chk = QCheckBox("Celerity (m/s):")
        self.celerity_chk.setToolTip(
            "Infrasound propagation speed.  Typical range: 300–360 m/s.\n"
            "Requires Origin Time to be set.")
        self.celerity_spin = QDoubleSpinBox()
        self.celerity_spin.setRange(100.0, 500.0)
        self.celerity_spin.setValue(340.0)
        self.celerity_spin.setDecimals(1)
        self.celerity_spin.setSingleStep(5.0)
        self.celerity_spin.setEnabled(False)
        self.celerity_chk.toggled.connect(self.celerity_spin.setEnabled)

        row3.addWidget(self.origin_time_chk)
        row3.addWidget(self.origin_time_input)
        row3.addSpacing(24)
        row3.addWidget(self.celerity_chk)
        row3.addWidget(self.celerity_spin)
        row3.addStretch()

        note = QLabel("  ← optional: draw expected arrival on beam waveform")
        note.setStyleSheet("color:#888; font-style:italic; font-size:8pt;")
        row3.addWidget(note)

        param_layout.addLayout(row3)

        param_grp.setLayout(param_layout)
        layout.addWidget(param_grp)

        # ── 4. Status + Progress + Run ────────────────────────────────────
        self.status_label = QLabel("Ready.")
        self.status_label.setStyleSheet(
            "font-style:italic; color:#555; padding:2px 6px;"
            "background:rgba(0,0,0,0.04); border-radius:4px;")
        layout.addWidget(self.status_label)

        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)

        self.process_btn = QPushButton("▶  Run FK Analysis")
        self.process_btn.setStyleSheet(
            "QPushButton{background:#27ae60;font-size:11pt;padding:10px;}"
            "QPushButton:hover{background:#1e8449;}")
        self.process_btn.clicked.connect(self.process_data)
        layout.addWidget(self.process_btn)

        sep = QFrame(); sep.setFrameShape(QFrame.HLine); sep.setFrameShadow(QFrame.Sunken)
        layout.addWidget(sep)

        tools_row = QHBoxLayout()
        xml_creator_btn = QPushButton("🛠  XML Creator / Editor")
        xml_creator_btn.setStyleSheet(
            "QPushButton{background:#8e44ad;padding:6px 14px;}"
            "QPushButton:hover{background:#6c3483;}")
        xml_creator_btn.clicked.connect(self._open_xml_creator)
        tools_row.addWidget(xml_creator_btn)

        db_btn = QPushButton("🗄  Event Database")
        db_btn.setStyleSheet(
            "QPushButton{background:#16a085;padding:6px 14px;}"
            "QPushButton:hover{background:#1abc9c;}")
        db_btn.clicked.connect(self._open_db_browser)
        tools_row.addWidget(db_btn)

        spectro_btn = QPushButton("📊  Plot Spectrogram")
        spectro_btn.setStyleSheet(
            "QPushButton{background:#d35400;padding:6px 14px;}"
            "QPushButton:hover{background:#a04000;}")
        spectro_btn.clicked.connect(self._open_spectrogram)
        tools_row.addWidget(spectro_btn)

        tools_row.addStretch()
        layout.addLayout(tools_row)

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
    def _resolve_inventory(self):
        """Resolve the inventory selected in inv_combo into an Inventory.

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

        # ── 2. Resolve the inventory robustly ─────────────────────────────
        try:
            inventory = self._resolve_inventory()
        except ValueError as e:
            QMessageBox.warning(self, "Inventory required", str(e))
            return

        # ── 3. Time window — picked window, else whole stream ─────────────
        # self.selected_time is a relative-seconds offset into the preview;
        # absolute start = stream[0].starttime + selected_time.
        stream = self.stream
        if self.selected_time is not None:
            try:
                t0 = self.stream[0].stats.starttime + float(self.selected_time)
                t1 = t0 + float(self.win_len.value())
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
        params = {
            "filter_freqmin": self.min_freq.value(),
            "filter_freqmax": self.max_freq.value(),
            "freq_max":   5.0,
            "nperseg":    512,
            "noverlap":   460,
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

    # ── inventory helpers ─────────────────────────────────────────────────
    def refresh_inventory_list(self):
        self.inv_combo.clear()
        base_dir = os.path.dirname(os.path.abspath(__file__))
        xml_dir    = os.path.join(base_dir, 'XML')
        xml_im_dir = os.path.join(base_dir, 'XML_IM')

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
        elif not os.path.isdir(xml_im_dir):
            QMessageBox.warning(self, "Warning", "No XML or XML_IM directory found.")

    def add_inventory_file(self):
        fname, _ = QFileDialog.getOpenFileName(
            self, "Select XML Inventory", "", "XML Files (*.xml);;All Files (*)")
        if not fname:
            return
        xml_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'XML')
        os.makedirs(xml_dir, exist_ok=True)
        dest = os.path.join(xml_dir, os.path.basename(fname))
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
            if len(fnames) == 1:
                self.mseed_edit.setText(fnames[0])
            else:
                self.mseed_edit.setText(
                    f"{len(fnames)} files merged  [{', '.join(os.path.basename(f) for f in fnames)}]")
            self._update_preview()
        except Exception as e:
            QMessageBox.critical(self, "Load Error", str(e))

    # ── Waveform preview ──────────────────────────────────────────────────
    def _update_preview(self):
        if not self.stream:
            return

        n     = len(self.stream)
        colors = plt.cm.tab10.colors

        self.figure.clear()
        axes = self.figure.subplots(n, 1, sharex=True)
        if n == 1:
            axes = [axes]

        for i, tr in enumerate(self.stream):
            ax    = axes[i]
            times = np.arange(len(tr.data)) * tr.stats.delta
            norm  = np.max(np.abs(tr.data)) or 1
            data  = tr.data / norm
            col   = colors[i % len(colors)]

            ax.plot(times, data, color=col, lw=0.7)
            ax.fill_between(times, data, 0, where=data >= 0,
                            color=col, alpha=0.20)
            ax.fill_between(times, data, 0, where=data <  0,
                            color=col, alpha=0.20)
            ax.set_ylim(-1.15, 1.15)
            ax.set_ylabel("Norm.", fontsize=7)
            ax.grid(True, alpha=0.25, ls='--')
            ax.tick_params(labelsize=7)

            label = (f"{tr.stats.network}.{tr.stats.station}."
                     f"{tr.stats.channel}  "
                     f"Fs={tr.stats.sampling_rate:.0f} Hz  "
                     f"[{tr.stats.starttime.strftime('%Y-%m-%d %H:%M:%S')}]")
            ax.set_title(label, fontsize=8, loc='left', pad=2)

            if self.selected_time is not None:
                ax.axvline(self.selected_time, color='red',
                           lw=1.4, ls='--', alpha=0.85)
                if i == 0:
                    pt = tr.stats.starttime + self.selected_time
                    ax.annotate(
                        f"Start: {pt.strftime('%H:%M:%S')}",
                        xy=(self.selected_time, 0.9), xytext=(6, 0),
                        textcoords='offset points', color='red', fontsize=8,
                        bbox=dict(boxstyle='round,pad=0.3',
                                  fc='white', ec='red', alpha=0.85))

        axes[-1].set_xlabel("Time (s)", fontsize=8)
        self.figure.tight_layout(pad=0.4, h_pad=0.3)
        self.canvas.mpl_connect('button_press_event', self._on_preview_click)
        self.view_window_btn.setEnabled(True)
        self.canvas.draw()

    def _on_preview_click(self, event):
        if not (event.inaxes and event.button == 1):
            return
        self.selected_time = event.xdata
        # Always compute absolute time from first trace start
        pt = self.stream[0].stats.starttime + self.selected_time
        self.time_label.setText(
            f"Pick: {pt.strftime('%Y-%m-%d %H:%M:%S.%f')[:-4]}")
        self.start_time_input.setDateTime(QDateTime(
            pt.year, pt.month, pt.day, pt.hour, pt.minute, pt.second,
            int(pt.microsecond / 1000)))
        self.clear_pick_btn.setEnabled(True)
        self._update_preview()

    def clear_time_pick(self):
        self.selected_time = None
        self.time_label.setText("Pick: —")
        self.clear_pick_btn.setEnabled(False)
        self.start_time_input.setDateTime(QDateTime.currentDateTime())
        self._update_preview()

    def view_waveform(self):
        if not self.stream:
            return
        viewer = WaveformViewer(self.stream, self)
        viewer.time_selected.connect(self._apply_viewer_pick)
        viewer.exec_()

    def _apply_viewer_pick(self, time_str):
        try:
            t = UTCDateTime(time_str)
            self.start_time_input.setDateTime(QDateTime(
                t.year, t.month, t.day, t.hour, t.minute, t.second))
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
        clr = "#e74c3c" if val < 34 else "#f39c12" if val < 67 else "#27ae60"
        self.progress_bar.setStyleSheet(f"""
            QProgressBar {{ border:2px solid {clr}; border-radius:8px;
                            text-align:center; }}
            QProgressBar::chunk {{ background:{clr}; border-radius:6px; }}""")

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
        n_sta       = int(r.get('n_stations', 0))
        exp_baz     = r['expected_bazi']

        fig = plt.figure(figsize=(15, 11))

        gs = gridspec.GridSpec(
            4, 2,
            width_ratios=[2.5, 1.2],
            hspace=0.45, wspace=0.32,
            left=0.08, right=0.97,
            top=0.84,  bottom=0.07,
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
        ax1.set_ylabel('Fisher', fontsize=9)
        ax1.legend(loc='upper right', fontsize=7)
        ax1.set_facecolor('#fafafa')
        ax1.grid(True, alpha=0.3)

        # ── Panel 2 — Back-azimuth ─────────────────────────────────────────
        ax2.scatter(r['time'][noi_mask], r['bazi'][noi_mask], **kw_noise)
        ax2.scatter(r['time'][det_mask], r['bazi'][det_mask], **kw_det)
        ax2.axhline(exp_baz, color='royalblue', ls='--', lw=1.6,
                    label=f"Expected  {int(exp_baz)}°", zorder=4)
        # ±20° acceptance band
        baz_lo = (exp_baz - 20) % 360
        baz_hi = (exp_baz + 20) % 360
        if baz_lo < baz_hi:
            ax2.axhspan(baz_lo, baz_hi, color='royalblue', alpha=0.08, zorder=0)
        ax2.set_ylim(0, 360); ax2.set_yticks([0, 90, 180, 270, 360])
        ax2.set_ylabel('Back-Az (°)', fontsize=9)
        ax2.legend(loc='upper right', fontsize=7)
        ax2.set_facecolor('#fafafa')
        ax2.grid(True, alpha=0.3)

        # ── Panel 3 — Apparent velocity ────────────────────────────────────
        ax3.scatter(r['time'][noi_mask], r['app_vel'][noi_mask], **kw_noise)
        ax3.scatter(r['time'][det_mask], r['app_vel'][det_mask], **kw_det)
        ax3.axhspan(300, 380, color='limegreen', alpha=0.08, zorder=0,
                    label='300–380 m/s')
        ax3.set_ylabel('App. Vel. (m/s)', fontsize=9)
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
        ax4.set_ylabel('Pressure (Pa)', fontsize=9)
        ax4.set_facecolor('#fafafa')
        ax4.grid(True, alpha=0.3)
        ax4.text(0.02, 0.94, f"Max |p| = {max_val:.4f} Pa",
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

        # Detection count annotation on Fisher panel
        ax1.text(0.99, 0.94,
                 f"{n_det}/{len(semb)} windows detected",
                 transform=ax1.transAxes, fontsize=8, ha='right',
                 color='crimson', fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.3',
                           facecolor='white', edgecolor='crimson', alpha=0.85))

        # ── Shared time axis formatting ────────────────────────────────────
        fmt = mdates.DateFormatter('%H:%M:%S')
        for ax in (ax1, ax2, ax3, ax4):
            ax.xaxis.set_major_locator(mdates.AutoDateLocator())
            ax.xaxis.set_major_formatter(fmt)
        for ax in (ax1, ax2, ax3):
            plt.setp(ax.get_xticklabels(), visible=False)
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
        exp_rad = np.deg2rad(exp_baz)
        ax5.plot([exp_rad, exp_rad], [200, 450],
                 color='royalblue', lw=2.0, ls='--',
                 label=f'Expected {int(exp_baz)}°', zorder=4)

        # Radial axis = apparent velocity
        ax5.set_ylim(200, 450)
        ax5.set_yticks([250, 300, 350, 400, 450])
        ax5.set_yticklabels(['250', '300', '350', '400', '450\nm/s'],
                            fontsize=6, color='#444')
        ax5.set_rlabel_position(45)

        # Cardinal labels
        ax5.set_thetagrids(range(0, 360, 45),
                           ['N', 'NE', 'E', 'SE', 'S', 'SW', 'W', 'NW'],
                           fontsize=8)
        ax5.set_title('Back-Azimuth Detection Map',
                      fontsize=10, fontweight='bold', pad=14)
        ax5.legend(loc='upper right', fontsize=7,
                   bbox_to_anchor=(1.35, 1.12))
        ax5.grid(True, alpha=0.3)

        # Panel 6 — Array Geometry ─────────────────────────────────────────
        x_km    = r.get('station_x_km', np.array([]))
        y_km    = r.get('station_y_km', np.array([]))
        sta_ids = r.get('station_ids',  [])

        ax6.scatter(x_km, y_km, s=60, color='steelblue',
                    edgecolors='k', linewidths=0.7, zorder=3)
        for xi, yi, sid in zip(x_km, y_km, sta_ids):
            ax6.annotate(sid, (xi, yi),
                         textcoords='offset points', xytext=(5, 5),
                         fontsize=7, color='#2c3e50')
        ax6.plot(0, 0, 'k+', ms=12, mew=2, zorder=4, label='Centroid')
        ax6.set_xlabel('East offset  (km)', fontsize=9)
        ax6.set_ylabel('North offset  (km)', fontsize=9)
        ax6.set_title('Array Geometry', fontsize=10, fontweight='bold')
        ax6.set_aspect('equal', adjustable='datalim')
        ax6.legend(fontsize=7, loc='upper right')
        ax6.grid(True, alpha=0.25, ls='--')

        # Suptitle ─────────────────────────────────────────────────────────
        fig.suptitle(
            f"Event: {event_name}    "
            f"Filter: {fmin}–{fmax} Hz    "
            f"Sensors: {n_sta}    "
            f"Expected BAZ: {int(exp_baz)}°    "
            f"Detections: {n_det}/{len(semb)}",
            fontsize=11, y=0.98, fontweight='bold')

        # Semblance colorbar (detections only, plasma scale) ──────────────
        cbar_ax = fig.add_axes([0.08, 0.895, 0.89, 0.012])
        cb = fig.colorbar(sc, cax=cbar_ax, orientation='horizontal', extend='min')
        cb.set_label(f'Semblance  (detections ≥ {semb_thresh:.2f}  |  grey = noise)',
                     fontsize=9, fontweight='bold', labelpad=-1)
        cb.ax.tick_params(labelsize=8)

        # Auto-save JPG ────────────────────────────────────────────────────
        safe_name = "".join(c if c.isalnum() or c in '-_' else '_'
                            for c in event_name)
        jpg_file  = f"FK_{safe_name}.jpg"
        try:
            fig.savefig(jpg_file, dpi=300, bbox_inches='tight',
                        format='jpeg', quality=92)
            print(f"[INFO] Figure saved → {jpg_file}")
        except Exception as e:
            print(f"[WARN] Could not auto-save figure: {e}")

        # Attach figure path to result so SaveEventDialog can reference it
        r['figure_path'] = os.path.abspath(jpg_file)

        win = ResultsWindow(fig, result=r, params=params, parent=self)
        win.show()


# ═══════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    app = QApplication(sys.argv)
    win = FKAnalysisGUI()
    win.show()
    sys.exit(app.exec_())
