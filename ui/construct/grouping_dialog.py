"""Grouping dialog: which read belongs to which clone (and primer).

A table "file / clone / primer" with editable clone and primer cells; reads
without a clone are highlighted. A custom regex (named groups ``clone`` and
``primer``, matched against the file stem) regroups all reads at once as a
live preview; an invalid pattern gets a red frame and disables OK.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from chrom_core.clones import DEFAULT_PRIMER_TOKENS, Grouping, ReadRef, group_reads
from i18n import tr
from ui import status
from ui.tables import configure_table
from ui.theme import theme as palette

COL_FILE, COL_CLONE, COL_PRIMER = range(3)
_PATH_ROLE = Qt.ItemDataRole.UserRole + 2
_INVALID_STYLE = "QLineEdit { border: 2px solid #d32f2f; }"


class GroupingDialog(QDialog):
    def __init__(self, grouping: Grouping, pattern: str, parent: QWidget | None = None,
                 primer_tokens: Iterable[str] = DEFAULT_PRIMER_TOKENS) -> None:
        super().__init__(parent)
        self.setWindowTitle(tr("grouping.title"))
        self.resize(640, 480)
        self._tokens = tuple(primer_tokens)
        self._paths = [r.path for r in grouping.reads] + list(grouping.unassigned)
        self._invalid = False

        self.pattern_edit = QLineEdit(pattern)
        self.pattern_edit.setPlaceholderText(tr("grouping.pattern_placeholder"))
        self.pattern_edit.setToolTip(tr("grouping.pattern_tip"))
        self._error = QLabel()
        self._error.setStyleSheet("color: #d32f2f;")
        self._summary = QLabel()

        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(
            [tr("grouping.col.file"), tr("grouping.col.clone"), tr("grouping.col.primer")])
        configure_table(self.table, stretch_column=COL_FILE)
        self.table.setAlternatingRowColors(False)   # the highlight must stay readable
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.itemChanged.connect(self._item_changed)

        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)

        row = QHBoxLayout()
        row.addWidget(QLabel(tr("grouping.pattern")))
        row.addWidget(self.pattern_edit, 1)
        layout = QVBoxLayout(self)
        layout.addLayout(row)
        layout.addWidget(self._error)
        layout.addWidget(self.table, 1)
        layout.addWidget(self._summary)
        layout.addWidget(self.buttons)

        self._fill(grouping)
        self._error.hide()
        self.pattern_edit.textChanged.connect(self._pattern_changed)

    # ------------------------------------------------------------- public API
    def result_grouping(self) -> Grouping:
        grouping = Grouping()
        for row in range(self.table.rowCount()):
            path = Path(self.table.item(row, COL_FILE).data(_PATH_ROLE))
            clone = self.table.item(row, COL_CLONE).text().strip()
            primer = self.table.item(row, COL_PRIMER).text().strip()
            if clone:
                grouping.reads.append(ReadRef(path, clone, primer))
            else:
                grouping.unassigned.append(path)
        return grouping

    def result_pattern(self) -> str:
        return self.pattern_edit.text()

    def pattern_invalid(self) -> bool:
        return self._invalid

    def unassigned_rows(self) -> list[int]:
        return [row for row in range(self.table.rowCount())
                if not self.table.item(row, COL_CLONE).text().strip()]

    # --------------------------------------------------------------- private
    def _pattern_changed(self, text: str) -> None:
        try:
            grouping = group_reads(self._paths, text or None, self._tokens)
        except ValueError as exc:
            self._set_invalid(str(exc))
            return
        self._set_invalid("")
        self._fill(grouping)

    def _set_invalid(self, message: str) -> None:
        self._invalid = bool(message)
        self.pattern_edit.setStyleSheet(_INVALID_STYLE if message else "")
        self._error.setText(message)
        self._error.setVisible(bool(message))
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(not message)

    def _fill(self, grouping: Grouping) -> None:
        refs = {r.path: r for r in grouping.reads}
        self.table.blockSignals(True)
        self.table.setRowCount(len(self._paths))
        for row, path in enumerate(self._paths):
            ref = refs.get(path)
            name = QTableWidgetItem(path.name)
            name.setData(_PATH_ROLE, str(path))
            name.setToolTip(str(path))
            name.setFlags(name.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self.table.setItem(row, COL_FILE, name)
            self.table.setItem(row, COL_CLONE, QTableWidgetItem(ref.clone if ref else ""))
            self.table.setItem(row, COL_PRIMER, QTableWidgetItem(ref.primer if ref else ""))
        self.table.blockSignals(False)
        self._highlight()

    def _item_changed(self, _item: QTableWidgetItem) -> None:
        self._highlight()

    def _highlight(self) -> None:
        """Rows without a clone get the theme highlight; the summary counts them."""
        colors = palette(status.current_theme())
        mark = QBrush(QColor(colors["highlight"]))
        unassigned = set(self.unassigned_rows())
        self.table.blockSignals(True)
        for row in range(self.table.rowCount()):
            for col in range(3):
                item = self.table.item(row, col)
                item.setBackground(mark if row in unassigned else QBrush())
                if row in unassigned and col == COL_CLONE:
                    item.setToolTip(tr("grouping.no_clone"))
                elif col == COL_CLONE:
                    item.setToolTip("")
        self.table.blockSignals(False)
        clones = {self.table.item(r, COL_CLONE).text().strip()
                  for r in range(self.table.rowCount())} - {""}
        self._summary.setText(tr("grouping.summary", clones=len(clones),
                                 reads=self.table.rowCount(), unassigned=len(unassigned)))
