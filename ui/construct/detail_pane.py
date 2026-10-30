"""Difference detail: every read of a clone centred on one reference position,
plus the clone's differences table, so the user can judge an error by eye."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QEvent, Qt, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QAbstractItemView, QLabel, QScrollArea, QSplitter, QTableWidget, QTableWidgetItem,
    QVBoxLayout, QWidget,
)

from chrom_core.align import ReadAlignment
from chrom_core.clones import CloneReport, load_read
from chrom_core.reference import Construct
from i18n import tr
from ui import status
from ui.chrom_view import ChromatogramWidget
from ui.read_view.analysis_dock import reference_overlay
from ui.tables import configure_table, paint_status, repaint_statuses
from ui.theme import theme as palette

CONTEXT = 10             # reference bases shown on each side of the position
VIEW_HEIGHT = 190        # compact height of one chromatogram
ROW_STATUS = {"confirmed": "errors", "disputed": "incomplete", "low_quality": "no_data"}
STATUS_KEY = {"confirmed": "matrix.legend.confirmed", "disputed": "matrix.legend.disputed",
              "low_quality": "matrix.legend.low_quality"}
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
        self._paths: dict[ChromatogramWidget, Path] = {}
        self._ref_pos: int | None = None
        self._zoom: float | None = None
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
        repaint_statuses(self.table)

    def show_clone(self, report: CloneReport, construct: Construct) -> None:
        self._report, self._construct = report, construct
        self._ref_pos, self._zoom = None, None
        self._ref_label.clear()
        for view in self._views:
            view.hide()
            view.deleteLater()
        self._views, self._results, self._paths = [], [], {}
        for result in report.reads:
            if result.alignment is None:
                continue
            try:
                chrom = load_read(result.ref.path)
            except Exception:  # unreadable file: nothing to draw for this read
                continue
            view = ChromatogramWidget()
            view.setMinimumSize(0, VIEW_HEIGHT)
            view.setMaximumHeight(VIEW_HEIGHT)
            view.set_theme(palette(self._theme_name))
            view.set_chromatogram(chrom)
            letters, mapping = reference_overlay(result.alignment, construct.seq, len(chrom.seq))
            view.set_reference(letters, mapping)
            view.setToolTip(f"{result.ref.primer} · {result.ref.path.name} ({result.alignment.strand})")
            view.installEventFilter(self)
            view.view_changed.connect(lambda v=view: self._view_changed(v))
            self._paths[view] = result.ref.path
            self._views.append(view)
            self._results.append(result)
            self._stack_layout.insertWidget(self._stack_layout.count() - 1, view)
        self._fill_table(report)

    def show_position(self, ref_pos: int) -> None:
        if self._report is None or self._construct is None:
            return
        self._ref_pos = ref_pos
        self._zoom = None
        self._select_row(ref_pos)
        self._ref_label.setText(self._reference_context(ref_pos))
        for view, result in zip(self._views, self._results):
            exact = ref_pos in result.alignment.ref_to_read
            idx = displayed_index(result.alignment, ref_pos, len(view.chrom.seq))
            view.set_highlights([(idx, idx + 1)] if exact and idx is not None else [])
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

    def _apply_view(self) -> None:
        if self._ref_pos is None or not self._views:
            return
        if self._zoom is None:
            spacing = self._views[0]._avg_spacing
            width = max(self._views[0].width(), 200)
            self._zoom = width / ((2 * CONTEXT + 1) * spacing)
        self._syncing = True
        try:
            for view, result in zip(self._views, self._results):
                view.zoom = self._zoom
                idx = displayed_index(result.alignment, self._ref_pos, len(view.chrom.seq))
                if idx is not None:
                    view.center_on_point(view.chrom.peak_indices[idx])
                view.update()
        finally:
            self._syncing = False

    def _view_changed(self, source: ChromatogramWidget) -> None:
        """A zoom change in one view is applied to all (shared scale)."""
        if self._syncing or self._ref_pos is None or source.zoom == self._zoom:
            return
        self._zoom = source.zoom
        self._apply_view()

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        if self._ref_pos is not None:
            self._zoom = None
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
                    paint_status(item, ROW_STATUS[d.status])
                self.table.setItem(row, col, item)
        self.table.blockSignals(False)

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
