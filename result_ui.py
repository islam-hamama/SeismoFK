"""
result_ui.py — the SeismoFK desktop theme.

Copyright (c) 2024-2026 Islam Hamama
Contact: islam.hamama@nriag.sci.eg
Licensed under the MIT License.

One palette, one font and one style sheet for the main window, the result
dialogs and the StationXML editor, so every window looks like the same
application.  Colours are defined once in ``PALETTE``; widgets opt into a role
with ``setObjectName`` (e.g. ``primaryAction``, ``dangerAction``,
``resultCard``) rather than inline colours.

Typography: the platform UI font (SF Pro on macOS, Segoe UI on Windows,
the desktop default on Linux), or Inter when installed.  Weights stay at
medium/semibold — heavier weights render as "Black" in several system fonts.
"""

from PyQt5.QtWidgets import QFrame, QHBoxLayout, QLabel, QVBoxLayout

PALETTE = dict(
    canvas="#f2f5f8",        # window background
    surface="#ffffff",       # cards
    surface_alt="#f7f9fb",   # inputs, table stripes
    line="#dde4ea",          # card borders
    line_strong="#c6d2db",   # input borders
    ink="#14263a",           # primary text
    ink_2="#34495c",         # secondary text
    muted="#6a7d8c",         # hints, captions
    navy="#0f2a3d",          # header bands
    navy_2="#1b3d55",
    accent="#0b7a75",        # primary actions, focus
    accent_hover="#09625e",
    accent_soft="#e4f2f1",
    accent_line="#b9dcd9",
    danger="#b3261e",
    danger_soft="#fbeceb",
    header_text="#f3f8fb",
    header_muted="#9fc3cc",
)

_BASE = """
QMainWindow, QDialog#analysisDialog {{ background:{canvas}; }}
QWidget {{ color:{ink}; }}
QToolTip {{ background:{navy}; color:{header_text}; border:0;
            padding:6px 8px; border-radius:4px; }}

/* Cards and sections: the group title sits above its card. */
QGroupBox {{ background:{surface}; border:1px solid {line}; border-radius:10px;
             margin-top:22px; padding:14px 14px 12px; }}
QGroupBox::title {{ subcontrol-origin:margin; subcontrol-position:top left;
                    left:4px; top:0; padding:0 0 5px 0; color:{muted};
                    font-size:9pt; font-weight:600; }}
QFrame#resultCard {{ background:{surface}; border:1px solid {line};
                     border-radius:10px; }}
QFrame#metricCard {{ background:{surface}; border:1px solid {line};
                     border-radius:8px; }}

/* Header band */
QFrame#resultHeader, QFrame#hero {{ background:{navy}; border-radius:10px; }}
QLabel#resultKicker, QLabel#heroSubtitle {{ color:{header_muted};
                                            font-size:9pt; font-weight:600; }}
QLabel#resultTitle {{ color:{header_text}; font-size:18pt; font-weight:600; }}
QLabel#heroTitle {{ color:{header_text}; font-size:20pt; font-weight:700; }}
QLabel#resultSubtitle {{ color:#c4d8e1; font-size:10pt; }}
QLabel#resultBadge, QLabel#versionBadge {{ color:#cdeeea; background:{navy_2};
    padding:6px 11px; border-radius:6px; font-size:9pt; font-weight:600; }}

/* Text roles */
QLabel#metricLabel {{ color:{muted}; font-size:9pt; font-weight:600; }}
QLabel#metricValue {{ color:{ink}; font-size:17pt; font-weight:600; }}
QLabel#resultInfo {{ color:{ink_2}; font-size:10pt; font-weight:500; }}
QLabel#mutedText, QLabel#sectionHint {{ color:{muted}; font-size:9pt; }}
QLabel#fieldLabel {{ color:{ink_2}; font-size:9pt; font-weight:500; }}
QLabel#sectionTitle {{ color:{ink}; font-weight:600; font-size:10pt; }}
QLabel#pickLabel {{ color:{accent}; font-weight:600; }}
QLabel#footerText {{ color:#8b9aa6; font-size:8pt; }}
QLabel#monoSummary {{ color:{ink_2}; background:{surface_alt};
                      border:1px solid {line}; border-radius:6px; padding:10px; }}

/* Buttons */
QPushButton {{ background:{surface}; color:{ink_2}; border:1px solid {line_strong};
               padding:7px 13px; border-radius:6px; font-weight:500; }}
QPushButton:hover {{ background:{surface_alt}; border-color:#9fb2bf; color:{ink}; }}
QPushButton:pressed {{ background:#e9eef2; }}
QPushButton:focus {{ border:1px solid {accent}; }}
QPushButton:disabled {{ color:#a3b0ba; background:#f4f6f8; border-color:#e3e8ec; }}
QPushButton#primaryAction {{ background:{accent}; color:white;
                             border:1px solid {accent}; font-weight:600; }}
QPushButton#primaryAction:hover {{ background:{accent_hover};
                                   border-color:{accent_hover}; }}
QPushButton#primaryAction:disabled {{ background:#d5e3e6; border-color:#d5e3e6;
                                      color:#7d939c; }}
QPushButton#methodAction {{ background:{accent_soft}; color:#07524f;
                            border:1px solid {accent_line}; font-weight:600;
                            padding:9px 13px; }}
QPushButton#methodAction:hover {{ background:#d3ebe9; border-color:{accent}; }}
QPushButton#methodAction:disabled {{ background:#f1f4f6; color:#9aabb5;
                                     border-color:#e1e7eb; }}
QPushButton#dangerAction {{ background:{surface}; color:{danger};
                            border:1px solid #e6bdb9; }}
QPushButton#dangerAction:hover {{ background:{danger_soft}; border-color:{danger}; }}
QPushButton#heroAction {{ background:transparent; color:#d7e7ee;
                          border:1px solid #37596f; font-weight:500; }}
QPushButton#heroAction:hover {{ background:{navy_2}; border-color:#5d8196;
                                color:white; }}
QPushButton#heroAction:disabled {{ color:#5f7b8c; border-color:#284a60; }}

/* Inputs */
QLineEdit, QComboBox, QDoubleSpinBox, QSpinBox, QDateTimeEdit, QTextEdit,
QPlainTextEdit {{ background:{surface_alt}; border:1px solid {line_strong};
                  border-radius:6px; padding:6px 7px;
                  selection-background-color:{accent}; selection-color:white; }}
QLineEdit:focus, QComboBox:focus, QDoubleSpinBox:focus, QSpinBox:focus,
QDateTimeEdit:focus, QTextEdit:focus, QPlainTextEdit:focus {{
                  border:1px solid {accent}; background:{surface}; }}
QLineEdit:disabled, QComboBox:disabled, QDoubleSpinBox:disabled,
QSpinBox:disabled, QDateTimeEdit:disabled {{ color:#9aa9b4; background:#f1f4f6; }}
QComboBox::drop-down {{ border:0; width:24px; }}
QComboBox QAbstractItemView {{ background:{surface}; border:1px solid {line};
                               selection-background-color:{accent_soft};
                               selection-color:{ink}; outline:0; }}
QCheckBox {{ spacing:8px; color:{ink_2}; }}

/* Tables */
QTableWidget, QTableView {{ background:{surface}; alternate-background-color:{surface_alt};
               border:1px solid {line}; border-radius:8px; gridline-color:#e8edf1;
               selection-background-color:{accent_soft}; selection-color:{ink}; }}
QHeaderView::section {{ background:{surface_alt}; color:{ink_2}; font-weight:600;
                        border:0; border-bottom:1px solid {line}; padding:7px 8px; }}

/* Progress, scrolling, toolbars */
QProgressBar {{ border:0; border-radius:3px; background:#dde6ec; height:6px;
                text-align:center; color:transparent; }}
QProgressBar::chunk {{ background:{accent}; border-radius:3px; }}
QScrollArea {{ border:0; background:transparent; }}
QWidget#sidebar {{ background:transparent; }}
QScrollArea#framedScroll {{ background:{surface}; border:1px solid {line};
                            border-radius:8px; }}
QScrollBar:vertical {{ background:transparent; width:10px; margin:2px; }}
QScrollBar:horizontal {{ background:transparent; height:10px; margin:2px; }}
QScrollBar::handle:vertical, QScrollBar::handle:horizontal {{
                  background:#c7d2da; border-radius:3px; min-height:28px;
                  min-width:28px; }}
QScrollBar::handle:hover {{ background:#9fb0bd; }}
QScrollBar::add-line, QScrollBar::sub-line {{ width:0; height:0; }}
QScrollBar::add-page, QScrollBar::sub-page {{ background:transparent; }}
QToolBar {{ background:transparent; border:0; spacing:4px; }}
QToolButton {{ border:0; padding:3px; min-width:24px; min-height:24px;
               border-radius:5px; }}
QToolButton:hover {{ background:#e8eef2; }}
QToolButton:checked {{ background:{accent_soft}; }}
QSplitter::handle {{ background:{line}; }}
"""

#: Style sheet for dialogs and any window that sets its own base font size.
RESULT_DIALOG_STYLE = _BASE.format(**PALETTE)

#: Style sheet for the main application windows (adds the base text size).
MAIN_WINDOW_STYLE = RESULT_DIALOG_STYLE + "\nQWidget { font-size:11pt; }\n"


def ui_font_family():
    """Inter when installed, otherwise the platform's UI font."""
    from PyQt5.QtGui import QFontDatabase

    if "Inter" in QFontDatabase().families():
        return "Inter"
    return QFontDatabase.systemFont(QFontDatabase.GeneralFont).family()


def apply_app_theme(app):
    """Set the application font and the shared style sheet on *app*."""
    from PyQt5.QtGui import QFont

    font = QFont(ui_font_family())
    font.setHintingPreference(QFont.PreferNoHinting)
    font.setStyleStrategy(QFont.PreferAntialias)
    app.setFont(font)
    app.setStyleSheet(RESULT_DIALOG_STYLE)


def apply_result_style(dialog):
    dialog.setObjectName("analysisDialog")
    dialog.setStyleSheet(RESULT_DIALOG_STYLE)


def apply_main_style(window):
    window.setStyleSheet(MAIN_WINDOW_STYLE)


def style_plot_toolbar(toolbar):
    from PyQt5.QtCore import QSize
    from PyQt5.QtGui import QColor, QIcon, QPainter

    toolbar.setIconSize(QSize(18, 18))
    toolbar.setMinimumHeight(34)
    # Matplotlib may choose white SVG icons from the native dark palette even
    # when the toolbar sits on a white result card. Recolor their existing
    # alpha masks to the secondary text color so every action stays visible.
    for action in toolbar.actions():
        if action.icon().isNull():
            continue
        icon = action.icon().pixmap(24, 24)
        painter = QPainter(icon)
        painter.setCompositionMode(QPainter.CompositionMode_SourceIn)
        painter.fillRect(icon.rect(), QColor(PALETTE["ink_2"]))
        painter.end()
        action.setIcon(QIcon(icon))


def result_header(title, subtitle, *, kicker="SEISMOFK / ANALYSIS", badge="v1.2.1"):
    frame = QFrame()
    frame.setObjectName("resultHeader")
    row = QHBoxLayout(frame)
    row.setContentsMargins(20, 14, 20, 14)
    row.setSpacing(12)
    heading = QVBoxLayout()
    heading.setSpacing(2)
    eyebrow = QLabel(kicker)
    eyebrow.setObjectName("resultKicker")
    _track(eyebrow, 1.2)
    name = QLabel(title)
    name.setObjectName("resultTitle")
    name.setWordWrap(True)
    detail = QLabel(subtitle)
    detail.setObjectName("resultSubtitle")
    detail.setWordWrap(True)
    heading.addWidget(eyebrow)
    heading.addWidget(name)
    heading.addWidget(detail)
    row.addLayout(heading, stretch=1)
    version = QLabel(badge)
    version.setObjectName("resultBadge")
    row.addWidget(version)
    return frame


def metric_card(label, value):
    frame = QFrame()
    frame.setObjectName("metricCard")
    column = QVBoxLayout(frame)
    column.setContentsMargins(15, 10, 15, 10)
    column.setSpacing(2)
    caption = QLabel(label.upper())
    caption.setObjectName("metricLabel")
    _track(caption, 0.8)
    number = QLabel(value)
    number.setObjectName("metricValue")
    column.addWidget(caption)
    column.addWidget(number)
    return frame


def _track(label, spacing):
    """Letter-spacing (px) for uppercase eyebrow text; QSS cannot set it."""
    from PyQt5.QtGui import QFont

    font = label.font()
    font.setLetterSpacing(QFont.AbsoluteSpacing, spacing)
    label.setFont(font)
