"""Shared table behaviour: content-sized columns, stretching last column, tooltips."""

from __future__ import annotations

from PySide6.QtWidgets import QHeaderView, QTableWidget, QTableWidgetItem

from ui.status import current_theme, status_color, status_text_color


def configure_table(table: QTableWidget) -> None:
    """Columns fit their content and the last one takes the remaining width."""
    header = table.horizontalHeader()
    header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
    header.setStretchLastSection(True)
    header.setMinimumSectionSize(40)
    table.verticalHeader().setVisible(False)
    table.setAlternatingRowColors(True)
    table.setWordWrap(False)


def set_full_text_tooltip(item: QTableWidgetItem) -> None:
    """Show the whole text on hover (it may be elided by a narrow column)."""
    text = item.text()
    if text:
        item.setToolTip(text)


def paint_status(item: QTableWidgetItem, verdict: str) -> None:
    """Colour a cell with the shared status palette (current theme)."""
    theme = current_theme()
    item.setBackground(status_color(verdict, theme))
    item.setForeground(status_text_color(verdict, theme))
