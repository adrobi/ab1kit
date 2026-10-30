"""Shared table behaviour: content-sized columns, stretching last column, tooltips."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QHeaderView, QTableWidget, QTableWidgetItem

from ui.status import current_theme, status_color, status_text_color


STATUS_ROLE = Qt.ItemDataRole.UserRole + 1
TOOLTIP_MIN_CHARS = 40


def configure_table(table: QTableWidget, stretch_column: int | None = None) -> None:
    """Columns fit their content; ``stretch_column`` (default: the last) takes the rest."""
    header = table.horizontalHeader()
    header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
    if stretch_column is None:
        header.setStretchLastSection(True)
    else:
        header.setStretchLastSection(False)
        header.setSectionResizeMode(stretch_column, QHeaderView.ResizeMode.Stretch)
    header.setMinimumSectionSize(40)
    table.verticalHeader().setVisible(False)
    table.setAlternatingRowColors(True)
    table.setWordWrap(False)


def set_full_text_tooltip(item: QTableWidgetItem) -> None:
    """Show the whole text on hover for long cells (a narrow column may elide them)."""
    text = item.text()
    if len(text) > TOOLTIP_MIN_CHARS:
        item.setToolTip(text)


def paint_status(item: QTableWidgetItem, verdict: str) -> None:
    """Colour a cell with the shared status palette (current theme)."""
    item.setData(STATUS_ROLE, verdict)
    theme = current_theme()
    item.setBackground(status_color(verdict, theme))
    item.setForeground(status_text_color(verdict, theme))


def repaint_statuses(table: QTableWidget) -> None:
    """Re-apply the status palette (after a theme switch) to every status cell."""
    for row in range(table.rowCount()):
        for col in range(table.columnCount()):
            item = table.item(row, col)
            verdict = item.data(STATUS_ROLE) if item is not None else None
            if verdict:
                paint_status(item, verdict)
