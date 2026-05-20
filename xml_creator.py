"""
xml_creator.py — StationXML Creator / Editor Tool for SeismoFK

Copyright (c) 2024-2025 Islam Hamama
Contact: islam.hamama@nriag.sci.eg

Licensed under the MIT License — see LICENSE for details.

Creates FDSN StationXML files compatible with obspy read_inventory().
Supports multi-station arrays, per-channel sensitivity, sensor presets,
loading/editing existing XMLs, CSV import, and one-click channel copy.
"""

import sys
import os
import csv
import xml.dom.minidom

from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QLineEdit, QFileDialog, QMessageBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QGroupBox,
    QComboBox, QDialog, QTextEdit, QSplitter, QAbstractItemView,
    QDialogButtonBox, QDoubleSpinBox, QDateTimeEdit, QFrame,
    QSizePolicy,
)
from PyQt5.QtCore import Qt, QDateTime
from PyQt5.QtGui  import QFont, QColor

# ─────────────────────────────────────────────────────────────────────────────
#  Sensor presets  {name: (channel_code, location, fs, sensitivity, freq, input_units, output_units)}
# ─────────────────────────────────────────────────────────────────────────────
PRESETS = {
    "Custom": ("BDF", "", 20.0, 1.0, 1.0, "PA", "COUNTS"),
    "Generic Infrasound — 10000 ct/Pa @ 20 Hz": ("BDF", "", 20.0, 10000.0, 1.0, "PA", "COUNTS"),
    "Generic Infrasound — 10000 ct/Pa @ 100 Hz": ("BDF", "", 100.0, 10000.0, 1.0, "PA", "COUNTS"),
    "HLW Array — 2.945e-5 ct/Pa @ 100 Hz": ("BDF", "", 100.0, 2.945e-5, 1.0, "PA", "COUNTS"),
    "HLW Array — 2.945e-5 ct/Pa @ 50 Hz": ("BDF", "", 50.0, 2.945e-5, 1.0, "PA", "COUNTS"),
    "IMS/CTBTO — 8000 ct/Pa @ 20 Hz": ("BDF", "", 20.0, 8000.0, 0.25, "PA", "COUNTS"),
    "IMS/CTBTO — 10000 ct/Pa @ 20 Hz": ("BDF", "", 20.0, 10000.0, 1.0, "PA", "COUNTS"),
    "Seismic — velocity (HHZ)": ("HHZ", "", 100.0, 1500.0, 1.0, "M/S", "COUNTS"),
    "Seismic — broadband (BHZ)": ("BHZ", "", 40.0, 1500.0, 1.0, "M/S", "COUNTS"),
}

# Column indices for station table
S_CODE, S_LAT, S_LON, S_ELEV, S_SITE = range(5)
STATION_HEADERS = ["Station Code", "Latitude", "Longitude", "Elevation (m)", "Site Name"]

# Column indices for channel table
C_CODE, C_LOC, C_FS, C_SENS, C_FREQ, C_INUNIT, C_OUTUNIT = range(7)
CHANNEL_HEADERS = ["Ch. Code", "Location", "Sample Rate (Hz)",
                   "Sensitivity", "Ref. Freq (Hz)", "Input Units", "Output Units"]


# ─────────────────────────────────────────────────────────────────────────────
#  XML builder
# ─────────────────────────────────────────────────────────────────────────────

def _e(val):
    """Escape XML characters in a string value."""
    return str(val).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def build_xml(net_code, net_start, net_end, source, stations):
    """
    Build a StationXML string from the provided data.

    stations : list of dicts with keys:
        code, lat, lon, elev, site, start, end,
        channels: list of dicts with keys:
            code, location, fs, sensitivity, freq, input_units, output_units,
            start, end, azimuth, dip, depth
    """
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<FDSNStationXML xmlns="http://www.fdsn.org/xml/station/1" schemaVersion="1">',
        f'    <Source>{_e(source)}</Source>',
        f'    <Created>{QDateTime.currentDateTimeUtc().toString("yyyy-MM-ddTHH:mm:ss")}</Created>',
        f'    <Network code="{_e(net_code)}" '
        f'startDate="{_e(net_start)}" endDate="{_e(net_end)}">',
        f'        <TotalNumberStations>{len(stations)}</TotalNumberStations>',
    ]

    for sta in stations:
        lines += [
            f'        <Station code="{_e(sta["code"])}" '
            f'startDate="{_e(sta["start"])}" endDate="{_e(sta["end"])}">',
            f'            <Latitude>{_e(sta["lat"])}</Latitude>',
            f'            <Longitude>{_e(sta["lon"])}</Longitude>',
            f'            <Elevation>{_e(sta["elev"])}</Elevation>',
            f'            <Site>',
            f'                <Name>{_e(sta["site"])}</Name>',
            f'            </Site>',
        ]
        for ch in sta.get('channels', []):
            lines += [
                f'            <Channel code="{_e(ch["code"])}" '
                f'locationCode="{_e(ch["location"])}" '
                f'startDate="{_e(ch["start"])}" endDate="{_e(ch["end"])}">',
                f'                <Latitude>{_e(sta["lat"])}</Latitude>',
                f'                <Longitude>{_e(sta["lon"])}</Longitude>',
                f'                <Elevation>{_e(sta["elev"])}</Elevation>',
                f'                <Depth>{_e(ch.get("depth", 0))}</Depth>',
                f'                <Azimuth>{_e(ch.get("azimuth", 0))}</Azimuth>',
                f'                <Dip>{_e(ch.get("dip", 0))}</Dip>',
                f'                <SampleRate>{_e(ch["fs"])}</SampleRate>',
                f'                <Response>',
                f'                    <InstrumentSensitivity>',
                f'                        <Value>{_e(ch["sensitivity"])}</Value>',
                f'                        <Frequency>{_e(ch["freq"])}</Frequency>',
                f'                        <InputUnits>',
                f'                            <Name>{_e(ch["input_units"])}</Name>',
                f'                        </InputUnits>',
                f'                        <OutputUnits>',
                f'                            <Name>{_e(ch["output_units"])}</Name>',
                f'                        </OutputUnits>',
                f'                    </InstrumentSensitivity>',
                f'                </Response>',
                f'            </Channel>',
            ]
        lines.append('        </Station>')

    lines += ['    </Network>', '</FDSNStationXML>']
    raw = '\n'.join(lines)
    # Pretty-print
    try:
        pretty = xml.dom.minidom.parseString(raw).toprettyxml(indent='    ')
        # Remove the extra XML declaration added by toprettyxml
        pretty = '\n'.join(pretty.split('\n')[1:])
        return '<?xml version="1.0" encoding="UTF-8"?>\n' + pretty
    except Exception:
        return raw


# ─────────────────────────────────────────────────────────────────────────────
#  XML Preview dialog
# ─────────────────────────────────────────────────────────────────────────────

class PreviewDialog(QDialog):
    def __init__(self, xml_text, parent=None):
        super().__init__(parent)
        self.setWindowTitle("XML Preview")
        self.setGeometry(200, 100, 900, 700)
        layout = QVBoxLayout(self)
        self.editor = QTextEdit()
        self.editor.setReadOnly(True)
        self.editor.setFont(QFont("Courier New", 9))
        self.editor.setPlainText(xml_text)
        layout.addWidget(self.editor)
        btns = QDialogButtonBox(QDialogButtonBox.Close)
        btns.rejected.connect(self.reject)
        layout.addWidget(btns)


# ─────────────────────────────────────────────────────────────────────────────
#  Main window
# ─────────────────────────────────────────────────────────────────────────────

class XMLCreatorGUI(QMainWindow):

    _STYLE = """
        QMainWindow  { background:#f4f6f8; }
        QWidget      { font-family:'Segoe UI',Arial,sans-serif; font-size:10pt; }
        QGroupBox    { font-weight:bold; border:1px solid #ccc; border-radius:6px;
                       margin-top:8px; padding-top:6px; }
        QGroupBox::title { subcontrol-origin:margin; left:10px; color:#2c3e50; }
        QPushButton  { background:#2980b9; color:white; border:none;
                       padding:6px 14px; border-radius:4px;
                       font-weight:bold; min-width:70px; }
        QPushButton:hover    { background:#1f618d; }
        QPushButton:disabled { background:#bdc3c7; color:#7f8c8d; }
        QTableWidget { border:1px solid #ccc; gridline-color:#e0e0e0; }
        QHeaderView::section { background:#2c3e50; color:white;
                               padding:4px; font-weight:bold; }
        QLineEdit, QComboBox, QDoubleSpinBox, QDateTimeEdit {
                       padding:4px; border:1px solid #bdc3c7;
                       border-radius:4px; background:white; }
        QLabel       { color:#2c3e50; }
    """

    def __init__(self):
        super().__init__()
        self.setStyleSheet(self._STYLE)
        self.setWindowTitle("SeismoFK — StationXML Creator")
        self.setGeometry(100, 60, 1400, 820)
        self._build_ui()

    # ── UI ────────────────────────────────────────────────────────────────
    def _build_ui(self):
        root = QWidget()
        main = QVBoxLayout(root)
        main.setSpacing(6)
        self.setCentralWidget(root)

        # ── Network header ────────────────────────────────────────────────
        net_grp    = QGroupBox("Network")
        net_layout = QHBoxLayout()

        self.net_code  = QLineEdit("EN");  self.net_code.setMaximumWidth(80)
        self.net_src   = QLineEdit("LOCAL"); self.net_src.setMaximumWidth(120)
        self.net_start = QDateTimeEdit(QDateTime(2000, 1, 1, 0, 0, 0))
        self.net_start.setDisplayFormat("yyyy-MM-dd HH:mm:ss")
        self.net_end   = QDateTimeEdit(QDateTime(2099, 12, 31, 23, 59, 59))
        self.net_end.setDisplayFormat("yyyy-MM-dd HH:mm:ss")

        for lbl, w in [("Network Code:", self.net_code),
                       ("Source:",       self.net_src),
                       ("Start:",        self.net_start),
                       ("End:",          self.net_end)]:
            net_layout.addWidget(QLabel(lbl))
            net_layout.addWidget(w)
        net_layout.addStretch()
        net_grp.setLayout(net_layout)
        main.addWidget(net_grp)

        # ── Splitter: stations (left) | channels (right) ──────────────────
        splitter = QSplitter(Qt.Horizontal)

        # ── Station panel ─────────────────────────────────────────────────
        sta_widget = QWidget()
        sta_layout = QVBoxLayout(sta_widget)
        sta_layout.setContentsMargins(0, 0, 0, 0)

        sta_grp    = QGroupBox("Stations")
        sta_inner  = QVBoxLayout()

        # station date range (shared default)
        date_row = QHBoxLayout()
        self.sta_start = QDateTimeEdit(QDateTime(2000, 1, 1, 0, 0, 0))
        self.sta_start.setDisplayFormat("yyyy-MM-dd HH:mm:ss")
        self.sta_end   = QDateTimeEdit(QDateTime(2099, 12, 31, 23, 59, 59))
        self.sta_end.setDisplayFormat("yyyy-MM-dd HH:mm:ss")
        date_row.addWidget(QLabel("Default station start:"))
        date_row.addWidget(self.sta_start)
        date_row.addWidget(QLabel("end:"))
        date_row.addWidget(self.sta_end)
        date_row.addStretch()
        sta_inner.addLayout(date_row)

        self.sta_table = QTableWidget(0, len(STATION_HEADERS))
        self.sta_table.setHorizontalHeaderLabels(STATION_HEADERS)
        self.sta_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.sta_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.sta_table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.sta_table.itemSelectionChanged.connect(self._on_station_selected)
        sta_inner.addWidget(self.sta_table)

        sta_btn_row = QHBoxLayout()
        add_sta  = QPushButton("+ Add Station")
        add_sta.clicked.connect(self.add_station)
        dup_sta  = QPushButton("⧉ Duplicate")
        dup_sta.clicked.connect(self.duplicate_station)
        rem_sta  = QPushButton("− Remove")
        rem_sta.clicked.connect(self.remove_station)
        imp_csv  = QPushButton("📂 Import CSV")
        imp_csv.clicked.connect(self.import_csv)
        for b in (add_sta, dup_sta, rem_sta, imp_csv):
            sta_btn_row.addWidget(b)
        sta_btn_row.addStretch()
        sta_inner.addLayout(sta_btn_row)

        # CSV format hint
        hint = QLabel("CSV columns: code, lat, lon, elev, site  (header row optional)")
        hint.setStyleSheet("color:#888; font-size:8pt; font-style:italic;")
        sta_inner.addWidget(hint)

        sta_grp.setLayout(sta_inner)
        sta_layout.addWidget(sta_grp)
        splitter.addWidget(sta_widget)

        # ── Channel panel ─────────────────────────────────────────────────
        ch_widget = QWidget()
        ch_layout = QVBoxLayout(ch_widget)
        ch_layout.setContentsMargins(0, 0, 0, 0)

        ch_grp   = QGroupBox("Channels  (select a station first)")
        ch_inner = QVBoxLayout()

        # Preset row
        preset_row = QHBoxLayout()
        self.preset_combo = QComboBox()
        self.preset_combo.addItems(PRESETS.keys())
        apply_preset = QPushButton("Apply Preset to New Channel")
        apply_preset.clicked.connect(self.apply_preset)
        copy_all_btn = QPushButton("Copy Channels → All Stations")
        copy_all_btn.clicked.connect(self.copy_channels_to_all)
        preset_row.addWidget(QLabel("Sensor Preset:"))
        preset_row.addWidget(self.preset_combo, stretch=1)
        preset_row.addWidget(apply_preset)
        preset_row.addWidget(copy_all_btn)
        ch_inner.addLayout(preset_row)

        # channel date range
        ch_date_row = QHBoxLayout()
        self.ch_start = QDateTimeEdit(QDateTime(2000, 1, 1, 0, 0, 0))
        self.ch_start.setDisplayFormat("yyyy-MM-dd HH:mm:ss")
        self.ch_end   = QDateTimeEdit(QDateTime(2099, 12, 31, 23, 59, 59))
        self.ch_end.setDisplayFormat("yyyy-MM-dd HH:mm:ss")
        ch_date_row.addWidget(QLabel("Default channel start:"))
        ch_date_row.addWidget(self.ch_start)
        ch_date_row.addWidget(QLabel("end:"))
        ch_date_row.addWidget(self.ch_end)
        ch_date_row.addStretch()
        ch_inner.addLayout(ch_date_row)

        self.ch_table = QTableWidget(0, len(CHANNEL_HEADERS))
        self.ch_table.setHorizontalHeaderLabels(CHANNEL_HEADERS)
        self.ch_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.ch_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        ch_inner.addWidget(self.ch_table)

        ch_btn_row = QHBoxLayout()
        add_ch  = QPushButton("+ Add Channel")
        add_ch.clicked.connect(self.add_channel)
        dup_ch  = QPushButton("⧉ Duplicate")
        dup_ch.clicked.connect(self.duplicate_channel)
        rem_ch  = QPushButton("− Remove")
        rem_ch.clicked.connect(self.remove_channel)
        for b in (add_ch, dup_ch, rem_ch):
            ch_btn_row.addWidget(b)
        ch_btn_row.addStretch()
        ch_inner.addLayout(ch_btn_row)

        ch_grp.setLayout(ch_inner)
        ch_layout.addWidget(ch_grp)
        splitter.addWidget(ch_widget)

        splitter.setSizes([500, 900])
        main.addWidget(splitter, stretch=1)

        # ── Bottom action bar ─────────────────────────────────────────────
        sep = QFrame(); sep.setFrameShape(QFrame.HLine); sep.setFrameShadow(QFrame.Sunken)
        main.addWidget(sep)

        action_row = QHBoxLayout()
        load_btn    = QPushButton("📂  Load XML")
        load_btn.clicked.connect(self.load_xml)
        preview_btn = QPushButton("🔍  Preview XML")
        preview_btn.clicked.connect(self.preview_xml)
        validate_btn = QPushButton("✔  Validate")
        validate_btn.clicked.connect(self.validate_xml)
        save_btn    = QPushButton("💾  Save XML")
        save_btn.setStyleSheet(
            "QPushButton{background:#27ae60;font-size:10pt;padding:8px 20px;}"
            "QPushButton:hover{background:#1e8449;}")
        save_btn.clicked.connect(self.save_xml)

        self.save_path_label = QLabel("Not saved yet.")
        self.save_path_label.setStyleSheet("color:#888; font-style:italic;")

        for w in (load_btn, preview_btn, validate_btn, save_btn):
            action_row.addWidget(w)
        action_row.addStretch()
        action_row.addWidget(self.save_path_label)
        main.addLayout(action_row)

        # Internal: map station row → channel list
        # Channels are stored per station row in self._channels dict
        self._channels = {}   # {sta_row: [ {code, location, fs, ...}, ... ]}
        self._current_sta_row = None

    # ── Station helpers ───────────────────────────────────────────────────

    def _new_station_row(self, code="STA", lat="0.0", lon="0.0",
                         elev="0.0", site=""):
        row = self.sta_table.rowCount()
        self.sta_table.insertRow(row)
        for col, val in enumerate([code, lat, lon, elev, site]):
            item = QTableWidgetItem(str(val))
            item.setTextAlignment(Qt.AlignCenter)
            self.sta_table.setItem(row, col, item)
        self._channels[row] = []
        return row

    def add_station(self):
        row = self._new_station_row()
        self.sta_table.selectRow(row)

    def duplicate_station(self):
        row = self._selected_sta_row()
        if row is None:
            return
        vals = [self.sta_table.item(row, c).text()
                for c in range(self.sta_table.columnCount())]
        new_row = self._new_station_row(*vals)
        import copy
        self._channels[new_row] = copy.deepcopy(self._channels.get(row, []))
        self.sta_table.selectRow(new_row)

    def remove_station(self):
        row = self._selected_sta_row()
        if row is None:
            return
        if QMessageBox.question(
                self, "Remove Station",
                f"Remove station at row {row + 1}?",
                QMessageBox.Yes | QMessageBox.No) == QMessageBox.No:
            return
        self.sta_table.removeRow(row)
        # Re-key channels dict
        new_ch = {}
        for k, v in self._channels.items():
            if k < row:
                new_ch[k] = v
            elif k > row:
                new_ch[k - 1] = v
        self._channels = new_ch
        self._current_sta_row = None
        self.ch_table.setRowCount(0)

    def _selected_sta_row(self):
        rows = self.sta_table.selectedItems()
        if not rows:
            return None
        return self.sta_table.currentRow()

    def _on_station_selected(self):
        row = self._selected_sta_row()
        if row is None or row == self._current_sta_row:
            return
        # Save current channel table back to storage
        self._save_current_channels()
        self._current_sta_row = row
        self._load_channels(row)

    def _save_current_channels(self):
        if self._current_sta_row is None:
            return
        chs = []
        for r in range(self.ch_table.rowCount()):
            ch = {}
            keys = ['code', 'location', 'fs', 'sensitivity',
                    'freq', 'input_units', 'output_units']
            for c, k in enumerate(keys):
                item = self.ch_table.item(r, c)
                ch[k] = item.text() if item else ''
            chs.append(ch)
        self._channels[self._current_sta_row] = chs

    def _load_channels(self, sta_row):
        self.ch_table.setRowCount(0)
        for ch in self._channels.get(sta_row, []):
            self._add_channel_row(ch)

    # ── Channel helpers ───────────────────────────────────────────────────

    def _add_channel_row(self, ch):
        row = self.ch_table.rowCount()
        self.ch_table.insertRow(row)
        keys = ['code', 'location', 'fs', 'sensitivity',
                'freq', 'input_units', 'output_units']
        defaults = ['BDF', '', '20', '10000', '1.0', 'PA', 'COUNTS']
        for col, (k, default) in enumerate(zip(keys, defaults)):
            val = ch.get(k, default)
            item = QTableWidgetItem(str(val))
            item.setTextAlignment(Qt.AlignCenter)
            self.ch_table.setItem(row, col, item)

    def add_channel(self):
        if self._current_sta_row is None:
            QMessageBox.warning(self, "No station", "Select a station first.")
            return
        preset_name = self.preset_combo.currentText()
        code, loc, fs, sens, freq, inunit, outunit = PRESETS.get(
            preset_name, PRESETS["Custom"])
        self._add_channel_row({
            'code': code, 'location': loc,
            'fs': str(fs), 'sensitivity': str(sens),
            'freq': str(freq),
            'input_units': inunit, 'output_units': outunit,
        })

    def duplicate_channel(self):
        row = self.ch_table.currentRow()
        if row < 0:
            return
        ch = {}
        keys = ['code', 'location', 'fs', 'sensitivity',
                'freq', 'input_units', 'output_units']
        for c, k in enumerate(keys):
            item = self.ch_table.item(row, c)
            ch[k] = item.text() if item else ''
        self._add_channel_row(ch)

    def remove_channel(self):
        row = self.ch_table.currentRow()
        if row < 0:
            return
        self.ch_table.removeRow(row)

    def apply_preset(self):
        """Apply selected preset to every selected channel row (or add new)."""
        if self._current_sta_row is None:
            QMessageBox.warning(self, "No station", "Select a station first.")
            return
        preset_name = self.preset_combo.currentText()
        code, loc, fs, sens, freq, inunit, outunit = PRESETS[preset_name]
        selected = self.ch_table.selectedItems()
        if selected:
            for row in set(i.row() for i in selected):
                vals = [code, loc, str(fs), str(sens), str(freq), inunit, outunit]
                for col, v in enumerate(vals):
                    self.ch_table.setItem(row, col, QTableWidgetItem(v))
        else:
            self.add_channel()

    def copy_channels_to_all(self):
        """Copy channels of the selected station to every other station."""
        if self._current_sta_row is None:
            QMessageBox.warning(self, "No station", "Select a source station first.")
            return
        self._save_current_channels()
        src = self._channels.get(self._current_sta_row, [])
        if not src:
            QMessageBox.warning(self, "No channels",
                                "The selected station has no channels to copy.")
            return
        import copy
        for r in range(self.sta_table.rowCount()):
            if r != self._current_sta_row:
                self._channels[r] = copy.deepcopy(src)
        QMessageBox.information(
            self, "Done",
            f"Channels copied to all {self.sta_table.rowCount() - 1} other stations.")

    # ── CSV import ────────────────────────────────────────────────────────

    def import_csv(self):
        fname, _ = QFileDialog.getOpenFileName(
            self, "Import Stations from CSV", "",
            "CSV Files (*.csv);;All Files (*)")
        if not fname:
            return
        try:
            with open(fname, newline='', encoding='utf-8') as f:
                reader = csv.reader(f)
                rows = [r for r in reader if r]

            # Auto-skip header if first cell is non-numeric
            start = 0
            try:
                float(rows[0][1])
            except (ValueError, IndexError):
                start = 1

            added = 0
            for row in rows[start:]:
                if len(row) < 4:
                    continue
                code = row[0].strip()
                lat  = row[1].strip()
                lon  = row[2].strip()
                elev = row[3].strip()
                site = row[4].strip() if len(row) > 4 else code
                self._new_station_row(code, lat, lon, elev, site)
                added += 1

            QMessageBox.information(self, "Imported",
                                    f"{added} stations imported from CSV.")
        except Exception as e:
            QMessageBox.critical(self, "Import Error", str(e))

    # ── Collect all data ──────────────────────────────────────────────────

    def _collect(self):
        """Return (net_meta, stations_list) from current UI state."""
        self._save_current_channels()

        net_meta = {
            'code':  self.net_code.text().strip() or "XX",
            'start': self.net_start.dateTime().toString("yyyy-MM-ddTHH:mm:ss"),
            'end':   self.net_end.dateTime().toString("yyyy-MM-ddTHH:mm:ss"),
            'source': self.net_src.text().strip() or "LOCAL",
        }

        sta_default_start = self.sta_start.dateTime().toString("yyyy-MM-ddTHH:mm:ss")
        sta_default_end   = self.sta_end.dateTime().toString("yyyy-MM-ddTHH:mm:ss")
        ch_default_start  = self.ch_start.dateTime().toString("yyyy-MM-ddTHH:mm:ss")
        ch_default_end    = self.ch_end.dateTime().toString("yyyy-MM-ddTHH:mm:ss")

        stations = []
        for r in range(self.sta_table.rowCount()):
            def cell(c):
                item = self.sta_table.item(r, c)
                return item.text().strip() if item else ''
            sta = {
                'code':  cell(S_CODE) or f"STA{r+1}",
                'lat':   cell(S_LAT)  or '0.0',
                'lon':   cell(S_LON)  or '0.0',
                'elev':  cell(S_ELEV) or '0.0',
                'site':  cell(S_SITE) or cell(S_CODE),
                'start': sta_default_start,
                'end':   sta_default_end,
                'channels': [],
            }
            for ch_dict in self._channels.get(r, []):
                sta['channels'].append({
                    'code':         ch_dict.get('code', 'BDF'),
                    'location':     ch_dict.get('location', ''),
                    'fs':           ch_dict.get('fs', '20'),
                    'sensitivity':  ch_dict.get('sensitivity', '10000'),
                    'freq':         ch_dict.get('freq', '1.0'),
                    'input_units':  ch_dict.get('input_units', 'PA'),
                    'output_units': ch_dict.get('output_units', 'COUNTS'),
                    'depth':        0,
                    'azimuth':      0,
                    'dip':          0,
                    'start':        ch_default_start,
                    'end':          ch_default_end,
                })
            stations.append(sta)
        return net_meta, stations

    # ── Actions ───────────────────────────────────────────────────────────

    def preview_xml(self):
        net, stations = self._collect()
        if not stations:
            QMessageBox.warning(self, "Empty", "Add at least one station first.")
            return
        xml_text = build_xml(net['code'], net['start'], net['end'],
                             net['source'], stations)
        dlg = PreviewDialog(xml_text, self)
        dlg.exec_()

    def validate_xml(self):
        net, stations = self._collect()
        errors = []
        if not net['code']:
            errors.append("Network code is empty.")
        for i, sta in enumerate(stations):
            if not sta['code']:
                errors.append(f"Station {i+1}: code is empty.")
            try:
                float(sta['lat']); float(sta['lon']); float(sta['elev'])
            except ValueError:
                errors.append(f"Station {i+1} ({sta['code']}): invalid lat/lon/elev.")
            if not sta['channels']:
                errors.append(f"Station {sta['code']}: has no channels.")
            for j, ch in enumerate(sta['channels']):
                try:
                    float(ch['fs'])
                except ValueError:
                    errors.append(f"Station {sta['code']} Ch{j+1}: invalid sample rate.")
                try:
                    float(ch['sensitivity'])
                except ValueError:
                    errors.append(f"Station {sta['code']} Ch{j+1}: invalid sensitivity.")
        if errors:
            QMessageBox.warning(self, "Validation Errors",
                                "\n".join(f"• {e}" for e in errors))
        else:
            # Also try parsing with obspy
            try:
                import tempfile
                from obspy import read_inventory
                xml_text = build_xml(net['code'], net['start'], net['end'],
                                     net['source'], stations)
                with tempfile.NamedTemporaryFile(suffix='.xml', delete=False,
                                                 mode='w', encoding='utf-8') as tf:
                    tf.write(xml_text)
                    tmp_path = tf.name
                read_inventory(tmp_path)
                os.unlink(tmp_path)
                QMessageBox.information(self, "Valid",
                                        "XML is valid and loads correctly with obspy.")
            except Exception as e:
                QMessageBox.warning(self, "obspy Validation Failed", str(e))

    def save_xml(self):
        net, stations = self._collect()
        if not stations:
            QMessageBox.warning(self, "Empty", "Add at least one station first.")
            return

        default_name = f"{net['code']}_array.xml"
        xml_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'XML')
        os.makedirs(xml_dir, exist_ok=True)

        fname, _ = QFileDialog.getSaveFileName(
            self, "Save StationXML", os.path.join(xml_dir, default_name),
            "XML Files (*.xml);;All Files (*)")
        if not fname:
            return
        if not fname.endswith('.xml'):
            fname += '.xml'

        try:
            xml_text = build_xml(net['code'], net['start'], net['end'],
                                 net['source'], stations)
            with open(fname, 'w', encoding='utf-8') as f:
                f.write(xml_text)
            self.save_path_label.setText(f"Saved: {os.path.basename(fname)}")
            QMessageBox.information(
                self, "Saved",
                f"StationXML saved to:\n{fname}\n\n"
                f"Stations: {len(stations)}\n"
                f"Total channels: {sum(len(s['channels']) for s in stations)}")
        except Exception as e:
            QMessageBox.critical(self, "Save Error", str(e))

    def load_xml(self):
        """Load an existing StationXML and populate the editor."""
        fname, _ = QFileDialog.getOpenFileName(
            self, "Load StationXML", "",
            "XML Files (*.xml);;All Files (*)")
        if not fname:
            return
        try:
            from obspy import read_inventory
            inv = read_inventory(fname)
        except Exception as e:
            QMessageBox.critical(self, "Load Error", str(e))
            return

        # Clear current state
        self.sta_table.setRowCount(0)
        self._channels = {}
        self._current_sta_row = None
        self.ch_table.setRowCount(0)

        for network in inv:
            self.net_code.setText(network.code)
            if network.start_date:
                d = network.start_date
                self.net_start.setDateTime(QDateTime(
                    d.year, d.month, d.day, d.hour, d.minute, d.second))
            if network.end_date:
                d = network.end_date
                self.net_end.setDateTime(QDateTime(
                    d.year, d.month, d.day, d.hour, d.minute, d.second))

            for station in network:
                row = self._new_station_row(
                    code=station.code,
                    lat=str(station.latitude),
                    lon=str(station.longitude),
                    elev=str(station.elevation),
                    site=station.site.name if station.site else station.code,
                )
                chs = []
                for ch in station:
                    sens_val  = ''
                    sens_freq = ''
                    in_unit   = ''
                    out_unit  = ''
                    if ch.response and ch.response.instrument_sensitivity:
                        s = ch.response.instrument_sensitivity
                        sens_val  = str(s.value)
                        sens_freq = str(s.frequency)
                        in_unit   = s.input_units or ''
                        out_unit  = s.output_units or ''
                    chs.append({
                        'code':         ch.code,
                        'location':     ch.location_code or '',
                        'fs':           str(ch.sample_rate),
                        'sensitivity':  sens_val,
                        'freq':         sens_freq,
                        'input_units':  in_unit,
                        'output_units': out_unit,
                    })
                self._channels[row] = chs

        if self.sta_table.rowCount() > 0:
            self.sta_table.selectRow(0)

        QMessageBox.information(self, "Loaded",
                                f"Loaded {self.sta_table.rowCount()} stations "
                                f"from {os.path.basename(fname)}.")


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    app = QApplication(sys.argv)
    win = XMLCreatorGUI()
    win.show()
    sys.exit(app.exec_())
