"""Expected-SNP verification dock: verdict table with CSV/XLSX buttons."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDockWidget,
    QHBoxLayout,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from i18n import tr
from ui.tables import configure_table, paint_status, set_full_text_tooltip


VERIFY_VERDICT_STATUS = {
    "confirmed": "ok",
    "heterozygous": "mixed",
    "not_found": "errors",
    "ref_mismatch": "errors",
    "low_quality": "incomplete",
    "no_data": "no_data",
    "unsupported": "no_data",
}


class VerifyDock(QDockWidget):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(tr("verify.title"), parent)
        self.verify_table = QTableWidget(0, 8)
        self.verify_table.setHorizontalHeaderLabels(
            [
                tr("verify.col.name"),
                tr("verify.col.pos"),
                tr("verify.col.ref"),
                tr("verify.col.alt"),
                tr("verify.col.verdict"),
                tr("verify.col.observed"),
                tr("verify.col.q"),
                tr("verify.col.note"),
            ]
        )
        configure_table(self.verify_table)
        self.verify_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.verify_table.setSortingEnabled(True)
        buttons = QWidget()
        row = QHBoxLayout(buttons)
        self.csv_button = QPushButton(tr("verify.save_csv"))
        self.xlsx_button = QPushButton(tr("verify.save_xlsx"))
        row.addWidget(self.csv_button)
        row.addWidget(self.xlsx_button)
        row.addStretch(1)
        container = QWidget()
        layout = QVBoxLayout(container)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self.verify_table)
        layout.addWidget(buttons)
        self.setWidget(container)

    def fill(self, verdicts: list) -> None:
        labels = {
            "confirmed": tr("verify.confirmed"),
            "heterozygous": tr("verify.heterozygous"),
            "not_found": tr("verify.not_found"),
            "no_data": tr("verify.no_data"),
            "ref_mismatch": tr("verify.ref_mismatch"),
            "low_quality": tr("verify.low_quality"),
            "unsupported": tr("verify.unsupported"),
        }
        self.verify_table.setSortingEnabled(False)
        self.verify_table.setRowCount(len(verdicts))
        for row, v in enumerate(verdicts):
            values = [
                v.variant.name,
                v.variant.position,
                v.variant.ref,
                v.variant.alt,
                labels.get(v.verdict, v.verdict),
                v.observed,
                v.quality if v.quality is not None else "",
                tr(f"verify.note.{v.note_code}", **v.note_args) if v.note_code else v.note,
            ]
            for col, value in enumerate(values):
                item = QTableWidgetItem()
                item.setData(Qt.ItemDataRole.DisplayRole, value)
                set_full_text_tooltip(item)
                if col == 4:
                    paint_status(item, VERIFY_VERDICT_STATUS.get(v.verdict, "no_data"))
                self.verify_table.setItem(row, col, item)
        self.verify_table.setSortingEnabled(True)
