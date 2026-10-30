"""Difference detail: every read of a clone centred on one reference position,
plus the clone's differences table, so the user can judge an error by eye."""

from __future__ import annotations

import struct
from pathlib import Path

from PySide6.QtCore import QEvent, Qt, Signal
from PySide6.QtGui import QColor, QFont
from PySide6.QtWidgets import (
    QAbstractItemView, QGraphicsOpacityEffect, QLabel, QScrollArea, QSplitter, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget,
)

from chrom_core.align import ReadAlignment
from chrom_core.clones import CloneReport, load_read
from chrom_core.reference import Construct
from i18n import tr
from ui import status
from ui.chrom_view import MAX_ZOOM, MIN_ZOOM, ChromatogramWidget
from ui.read_view.analysis_dock import reference_overlay
from ui.tables import STATUS_ROLE, configure_table
from ui.theme import theme as palette

CONTEXT = 10             # reference bases shown on each side of the position
VIEW_HEIGHT = 190        # compact height of one chromatogram
ROW_STATUS = {"confirmed": "errors", "disputed": "errors", "low_quality": "no_data"}
STATUS_KEY = {"confirmed": "matrix.legend.confirmed", "disputed": "matrix.legend.disputed",
              "low_quality": "matrix.legend.low_quality"}
READ_ERRORS = (OSError, ValueError, KeyError, TypeError, struct.error)
DISPUTED_TINT = 0.35     # share of the errors red mixed into the background
UNCOVERED_OPACITY = 0.35
EFFECT_CODES = ("synonymous", "frameshift", "inframe_ins", "inframe_del")


def displayed_index(alignment: ReadAlignment, ref_pos: int, read_length: int) -> int | None:
    """Index in the *displayed* read of the base aligned to ``ref_pos``.

    ``ref_to_read`` indexes the aligned orientation (the reverse complement on
    '-'), so the displayed index is ``n-1-idx`` there. A position deleted in the
    read falls back to the nearest mapped neighbour; ``None`` when not covered.
    """
    if not alignment.covers(ref_pos):
        return None
    mapping, length = alignment.ref_to_read, alignment.ref_length
    idx = None
    for step in range(len(mapping) + (length or 0) + 1):
        for pos in ((ref_pos + step, ref_pos - step) if step else (ref_pos,)):
            if length:
                pos %= length
            if pos in mapping:
                idx = mapping[pos]
                break
        if idx is not None:
            break
    if idx is None:
        return None
    return idx if alignment.strand == "+" else read_length - 1 - idx


class DetailPane(QWidget):
    open_read = Signal(Path)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._report: CloneReport | None = None
        self._construct: Construct | None = None
        self._views: list[ChromatogramWidget] = []
        self._results: list = []
        self._rows: list[QWidget] = []     # one container per read (view or placeholder)
        self._notes: list[QLabel] = []     # 'position not covered' label per view
        self._paths: dict[ChromatogramWidget, Path] = {}
        self._ref_pos: int | None = None
        self._zoom: float | None = None
        self._zoom_chosen = False     # the user zoomed: keep it across resizes
        self._syncing = False
        self._theme_name = status.current_theme()

        self._ref_label = QLabel()
        mono = QFont("Consolas")
        mono.setStyleHint(QFont.StyleHint.Monospace)
        self._ref_label.setFont(mono)
        self._ref_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        self._stack = QWidget()
        self._stack_layout = QVBoxLayout(self._stack)
        self._stack_layout.setContentsMargins(0, 0, 0, 0)
        self._stack_layout.setSpacing(2)
        self._stack_layout.addStretch(1)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self._stack)

        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.addWidget(self._ref_label)
        left_layout.addWidget(scroll, 1)

        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels([
            tr("detail.col.pos"), tr("detail.col.change"), tr("detail.col.effect"),
            tr("detail.col.support"), tr("detail.col.status"),
        ])
        configure_table(self.table, stretch_column=2)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.itemSelectionChanged.connect(self._row_selected)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(left)
        splitter.addWidget(self.table)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(splitter)

    # ---------------------------------------------------------------- API
    def views(self) -> list[ChromatogramWidget]:
        return list(self._views)

    def reference_text(self) -> str:
        return self._ref_label.text()

    def set_theme(self, name: str) -> None:
        self._theme_name = name
        for view in self._views:
            view.set_theme(palette(name))
        self._paint_table()

    def show_clone(self, report: CloneReport, construct: Construct) -> None:
        self._report, self._construct = report, construct
        self._ref_pos, self._zoom, self._zoom_chosen = None, None, False
        self._ref_label.clear()
        for row in self._rows:
            row.hide()
            row.deleteLater()
        self._views, self._results, self._paths, self._rows, self._notes = [], [], {}, [], []
        for result in report.reads:
            name = result.ref.path.name
            if result.alignment is None:
                error = tr("detail.err.no_alignment") if result.error == "no_alignment" else result.error
                self._add_row(QLabel(f"{result.ref.primer} · {name}: {error}"))
                continue
            try:
                chrom = load_read(result.ref.path)
            except READ_ERRORS as exc:
                self._add_row(QLabel(f"{result.ref.primer} · {name}: {exc}"))
                continue
            view = ChromatogramWidget()
            view.setMinimumSize(0, VIEW_HEIGHT)
            view.setMaximumHeight(VIEW_HEIGHT)
            view.set_theme(palette(self._theme_name))
            view.set_chromatogram(chrom)
            letters, mapping = reference_overlay(result.alignment, construct.seq, len(chrom.seq))
            view.set_reference(letters, mapping)
            title = f"{result.ref.primer} · {name} ({result.alignment.strand})"
            if result.ambiguous:
                title += " — " + tr("detail.ambiguous")
            view.setToolTip(title)
            view.installEventFilter(self)
            view.view_changed.connect(lambda v=view: self._view_changed(v))
            note = QLabel(tr("detail.uncovered"))
            note.hide()
            container = QWidget()
            box = QVBoxLayout(container)
            box.setContentsMargins(0, 0, 0, 0)
            box.setSpacing(0)
            box.addWidget(QLabel(title))
            box.addWidget(note)
            box.addWidget(view)
            self._paths[view] = result.ref.path
            self._views.append(view)
            self._results.append(result)
            self._notes.append(note)
            self._add_row(container)
        self._reset_default_zoom()
        self._fill_table(report)

    def _add_row(self, widget: QWidget) -> None:
        self._rows.append(widget)
        self._stack_layout.insertWidget(self._stack_layout.count() - 1, widget)

    def show_position(self, ref_pos: int) -> None:
        if self._report is None or self._construct is None:
            return
        self._ref_pos = ref_pos
        self._select_row(ref_pos)
        self._ref_label.setText(self._reference_context(ref_pos))
        for view, result, note in zip(self._views, self._results, self._notes):
            exact = ref_pos in result.alignment.ref_to_read
            idx = displayed_index(result.alignment, ref_pos, len(view.chrom.seq))
            view.set_highlights([(idx, idx + 1)] if exact and idx is not None else [])
            note.setVisible(idx is None)
            effect = None
            if idx is None:
                effect = QGraphicsOpacityEffect(view)
                effect.setOpacity(UNCOVERED_OPACITY)
            view.setGraphicsEffect(effect)
        self._apply_view()

    # ----------------------------------------------------------- internals
    def _reference_context(self, ref_pos: int) -> str:
        seq, n = self._construct.seq, len(self._construct.seq)

        def around(a: int, b: int) -> str:
            if self._construct.circular:
                return "".join(seq[i % n] for i in range(a, b))
            return seq[max(0, a):max(0, b)]

        return (f"{ref_pos + 1}:  {around(ref_pos - CONTEXT, ref_pos)}[{seq[ref_pos % n]}]"
                f"{around(ref_pos + 1, ref_pos + CONTEXT + 1)}")

    def _reset_default_zoom(self) -> None:
        """Default scale (about 21 bases visible), unless the user chose one."""
        if not self._views or self._zoom_chosen:
            return
        spacing = self._views[0]._avg_spacing
        width = max(self._views[0].width(), 200)
        self._zoom = min(MAX_ZOOM, max(MIN_ZOOM, width / ((2 * CONTEXT + 1) * spacing)))
        self._apply_view()

    def _apply_view(self) -> None:
        """Give every view the shared zoom and centre it on the position."""
        if self._zoom is None:
            return
        self._syncing = True
        try:
            for view, result in zip(self._views, self._results):
                view.zoom = self._zoom
                if self._ref_pos is not None:
                    idx = displayed_index(result.alignment, self._ref_pos, len(view.chrom.seq))
                    if idx is not None:
                        view.center_on_point(view.chrom.peak_indices[idx])
                view.update()
        finally:
            self._syncing = False

    def _view_changed(self, source: ChromatogramWidget) -> None:
        """A zoom change in one view is applied to all (shared scale)."""
        if self._syncing or source.zoom == self._zoom:
            return
        self._zoom, self._zoom_chosen = source.zoom, True
        self._apply_view()

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._reset_default_zoom()
        if self._zoom_chosen:
            self._apply_view()

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        if event.type() == QEvent.Type.MouseButtonDblClick and obj in self._paths:
            self.open_read.emit(self._paths[obj])
            return True
        return super().eventFilter(obj, event)

    def _fill_table(self, report: CloneReport) -> None:
        self.table.blockSignals(True)
        self.table.setRowCount(len(report.differences))
        for row, d in enumerate(report.differences):
            effect = tr(f"effect.{d.effect}") if d.effect in EFFECT_CODES else d.effect
            cells = [str(d.ref_pos + 1), f"{d.ref_base or '-'}→{d.read_base or '-'}", effect,
                     f"{d.support} / {d.against}", tr(STATUS_KEY[d.status])]
            for col, text in enumerate(cells):
                item = QTableWidgetItem(text)
                if col == 0:
                    item.setData(Qt.ItemDataRole.UserRole, d.ref_pos)
                if col == 4:
                    item.setData(STATUS_ROLE, d.status)
                if not d.in_region:
                    item.setToolTip(tr("matrix.legend.outside"))
                    font = item.font()
                    font.setItalic(True)
                    item.setFont(font)
                self.table.setItem(row, col, item)
        self.table.blockSignals(False)
        self._paint_table()

    def _paint_table(self) -> None:
        """Status cells: confirmed solid red, disputed a weak red tint, low quality
        grey (shared palette); rows outside the region dimmed. Current theme."""
        theme, colours = self._theme_name, palette(self._theme_name)
        base = QColor(colours["bg"])
        red = status.status_color("errors", theme)

        def mix(a: int, b: int) -> int:
            return round(a * DISPUTED_TINT + b * (1 - DISPUTED_TINT))

        for row, d in enumerate(self._report.differences if self._report else ()):
            for col in range(self.table.columnCount()):
                item = self.table.item(row, col)
                if item is None:
                    continue
                if col == 4 and d.status == "disputed":
                    item.setBackground(QColor(mix(red.red(), base.red()), mix(red.green(), base.green()),
                                              mix(red.blue(), base.blue())))
                    item.setForeground(QColor(colours["text"]))
                elif col == 4:
                    verdict = ROW_STATUS[d.status]
                    item.setBackground(status.status_color(verdict, theme))
                    item.setForeground(status.status_text_color(verdict, theme))
                elif not d.in_region:
                    item.setForeground(QColor(colours["dim"]))

    def _row_selected(self) -> None:
        rows = self.table.selectionModel().selectedRows()
        if rows:
            self.show_position(self.table.item(rows[0].row(), 0).data(Qt.ItemDataRole.UserRole))

    def _select_row(self, ref_pos: int) -> None:
        self.table.blockSignals(True)
        self.table.clearSelection()
        for row in range(self.table.rowCount()):
            if self.table.item(row, 0).data(Qt.ItemDataRole.UserRole) == ref_pos:
                self.table.selectRow(row)
                break
        self.table.blockSignals(False)
