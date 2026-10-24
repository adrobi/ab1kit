"""Batch panel: a QC table over a folder of .ab1 files."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QCheckBox,
    QHBoxLayout,
    QHeaderView,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from chrom_core.batch import FileReport
from i18n import tr

STATUS_TEXT = {"ok": "✓", "warn": "⚠", "bad": "✗", "error": "✗"}
STATUS_COLOR = {"ok": "#1e9e1e", "warn": "#e0a000", "bad": "#e02020", "error": "#e02020"}


class BatchPanel(QWidget):
    """Table of scan results; double-click opens a file in the viewer."""

    open_requested = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.reports: list[FileReport] = []

        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(
            [
                tr("batch.col.file"),
                tr("batch.col.length"),
                tr("batch.col.meanq"),
                tr("batch.col.q20"),
                tr("batch.col.trim"),
                tr("batch.col.status"),
            ]
        )
        self.table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.Stretch
        )
        self.table.verticalHeader().setVisible(False)
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.table.itemDoubleClicked.connect(self._double_clicked)
        layout.addWidget(self.table)

        buttons = QHBoxLayout()
        self.hide_bad_box = QCheckBox(tr("batch.hide_bad"))
        self.hide_bad_box.toggled.connect(self._apply_filter)
        buttons.addWidget(self.hide_bad_box)
        buttons.addStretch(1)
        self.csv_button = QPushButton(tr("batch.csv"))
        self.export_button = QPushButton(tr("batch.export"))
        buttons.addWidget(self.csv_button)
        buttons.addWidget(self.export_button)
        layout.addLayout(buttons)

    # ----------------------------------------------------------------- data
    def show_reports(self, reports: list[FileReport]) -> None:
        self.reports = list(reports)
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(reports))
        for row, report in enumerate(reports):
            self._set_item(row, 0, report.path.name, data=str(report.path))
            self._set_item(row, 1, report.length if not report.error else "")
            self._set_item(
                row, 2,
                "" if report.error else "—" if report.mean_q is None else report.mean_q,
            )
            self._set_item(row, 3, report.q20 if not report.error else "")
            self._set_item(row, 4, report.trimmed_length if not report.error else "")
            status = QTableWidgetItem(STATUS_TEXT[report.status])
            status.setForeground(QColor(STATUS_COLOR[report.status]))
            if report.error:
                status.setToolTip(report.error)
            self.table.setItem(row, 5, status)
        self.table.setSortingEnabled(True)
        self._apply_filter()

    def _set_item(self, row: int, col: int, value, data=None) -> None:
        item = QTableWidgetItem()
        item.setData(Qt.ItemDataRole.DisplayRole, value)
        if data is not None:
            item.setData(Qt.ItemDataRole.UserRole, data)
        self.table.setItem(row, col, item)

    def _apply_filter(self) -> None:
        hide = self.hide_bad_box.isChecked()
        for row in range(self.table.rowCount()):
            status = self.table.item(row, 5)
            bad = status is not None and status.text() == "✗"
            self.table.setRowHidden(row, hide and bad)

    def _double_clicked(self, item: QTableWidgetItem) -> None:
        path_item = self.table.item(item.row(), 0)
        if path_item is not None:
            path = path_item.data(Qt.ItemDataRole.UserRole)
            if path:
                self.open_requested.emit(str(path))
