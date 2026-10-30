"""Batch dock: a QC table over a folder of .ab1 files, with scan and export."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QCheckBox,
    QDockWidget,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from chrom_core import batch, export
from chrom_core.batch import FileReport
from chrom_core.config import AppConfig
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


class BatchDock(QDockWidget):
    """Dock around :class:`BatchPanel`; owns folder scanning and batch export.

    ``start_task`` is the shell's background-task runner; ``message`` is
    emitted for status-bar text (text, timeout in ms).
    """

    open_requested = Signal(str)
    message = Signal(str, int)

    def __init__(self, config: AppConfig, start_task: Callable, parent: QWidget | None = None) -> None:
        super().__init__(tr("batch.title"), parent)
        self.config = config
        self._start_task = start_task
        self.reports: list[batch.FileReport] = []
        self.batch_panel = BatchPanel()
        self.setWidget(self.batch_panel)
        self.batch_panel.open_requested.connect(self.open_requested)
        self.batch_panel.export_button.clicked.connect(self.export_trimmed)
        self.batch_panel.csv_button.clicked.connect(self.save_csv)

    def scan_folder(self, folder: Path) -> None:
        window, min_q = self.config.trim_window, self.config.trim_min_q

        def done(reports: object) -> None:
            self.reports = reports  # type: ignore[assignment]
            self.batch_panel.show_reports(self.reports)
            self.show()
            self.message.emit(tr("batch.scan_done", n=len(self.reports)), 5000)

        # scan_folder reports no progress and is not interruptible: on cancel
        # its late result is dropped and the table stays as it was.
        self._start_task(
            lambda progress, is_cancelled: batch.scan_folder(folder, window, min_q),
            done,
            tr("menu.file.open_folder"),
            key="scan",
        )

    def export_trimmed(self) -> None:
        if not self.reports:
            return
        out_dir = QFileDialog.getExistingDirectory(self, tr("batch.export_dir"))
        if not out_dir:
            return
        results = batch.batch_export(
            self.reports,
            Path(out_dir),
            window=self.config.trim_window,
            min_q=self.config.trim_min_q,
        )
        exported = sum(1 for _, dest, error in results if dest and not error)
        errors = sum(1 for _, _, error in results if error)
        self.message.emit(tr("batch.exported", n=exported, err=errors), 8000)

    def save_csv(self) -> None:
        if not self.reports:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, tr("batch.csv"), "batch_summary.csv", "CSV (*.csv)"
        )
        if path:
            export.write_text(path, batch.summary_csv(self.reports))
            self.message.emit(tr("export.done", path=path), 5000)
