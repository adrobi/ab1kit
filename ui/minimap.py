"""Whole-read overview strip: quality map, optional quality profile, viewport."""

from __future__ import annotations

from PySide6.QtCore import QPointF, Qt, Signal
from PySide6.QtGui import QColor, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import QWidget

from chrom_core import Chromatogram
from ui.theme import LIGHT

STRIP_H = 18
PROFILE_H = 40
PAD = 2


class MiniMap(QWidget):
    """One-line map of the whole read; click or drag to navigate.

    Shows a per-base quality strip (colored by phred bucket), optionally a
    quality profile line above it, plus a rectangle marking the main view's
    viewport.
    """

    request_center = Signal(float)  # trace point the main view should center on

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.chrom: Chromatogram | None = None
        self.theme = LIGHT
        self.show_profile = False
        self._offset = 0.0
        self._visible = 0.0
        self._total = 1.0
        self.highlights: list[tuple[int, int]] = []
        self.setMouseTracking(True)
        self._apply_height()

    # ------------------------------------------------------------- state API
    def set_chromatogram(self, chrom: Chromatogram) -> None:
        self.chrom = chrom
        self.highlights = []
        self._total = max(1, max(len(t) for t in chrom.traces.values()))
        self.update()

    def set_theme(self, theme: dict) -> None:
        self.theme = theme
        self.update()

    def set_profile(self, on: bool) -> None:
        self.show_profile = on
        self._apply_height()
        self.update()

    def set_view(self, offset: float, visible: float, total: float) -> None:
        self._offset = offset
        self._visible = visible
        self._total = max(1.0, total)
        self.update()

    def set_highlights(self, ranges: list[tuple[int, int]]) -> None:
        self.highlights = ranges
        self.update()

    def _apply_height(self) -> None:
        self.setFixedHeight(STRIP_H + PAD + (PROFILE_H if self.show_profile else 0))

    # --------------------------------------------------------- coordinate map
    def _h(self) -> int:
        return self.height()

    def _strip_y(self) -> int:
        return self._h() - STRIP_H - PAD

    def _x(self, point: float) -> float:
        return point / self._total * self.width()

    # ----------------------------------------------------------- interaction
    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self._navigate(event.position().x())

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        if event.buttons() & Qt.MouseButton.LeftButton:
            self._navigate(event.position().x())

    def _navigate(self, x: float) -> None:
        if self.chrom is not None:
            self.request_center.emit(x / max(self.width(), 1) * self._total)

    # -------------------------------------------------------------- painting
    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.fillRect(0, 0, self.width(), self._h(), QColor(self.theme["panel"]))
        if self.chrom is not None:
            if self.show_profile:
                self._paint_profile(painter)
            self._paint_strip(painter)
            self._paint_highlights(painter)
        self._paint_viewport(painter)
        painter.end()

    def _paint_strip(self, painter: QPainter) -> None:
        chrom = self.chrom
        assert chrom is not None
        painter.setPen(Qt.PenStyle.NoPen)
        y = self._strip_y()
        peaks = list(chrom.peak_indices) + [self._total]
        for i, quality in enumerate(chrom.quality):
            bucket = "low_q" if quality < 10 else "mid_q" if quality < 20 else "high_q"
            painter.setBrush(QColor(self.theme[bucket]))
            x0 = self._x(peaks[i])
            x1 = self._x(peaks[i + 1])
            if x1 - x0 < 0.5:
                x1 = x0 + 0.5
            painter.drawRect(int(x0), y, int(x1 - x0) + 1, STRIP_H - 2)

    def _paint_profile(self, painter: QPainter) -> None:
        chrom = self.chrom
        assert chrom is not None
        max_q = max(60, max(chrom.quality) if chrom.quality else 60)
        poly = QPolygonF()
        for i, quality in enumerate(chrom.quality):
            x = self._x(chrom.peak_indices[i])
            y = PROFILE_H - PAD - quality / max_q * (PROFILE_H - 2 * PAD)
            poly.append(QPointF(x, y))
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setPen(QPen(QColor(self.theme["dim"]), 1))
        painter.drawPolyline(poly)

    def _paint_highlights(self, painter: QPainter) -> None:
        if not self.highlights or not self.chrom:
            return
        color = QColor(self.theme["low_q"])
        color.setAlpha(140)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(color)
        peaks = self.chrom.peak_indices
        for start, end in self.highlights:
            x0 = self._x(peaks[max(0, start)])
            x1 = self._x(peaks[min(len(peaks) - 1, max(0, end - 1))])
            painter.drawRect(int(x0), 0, max(1, int(x1 - x0)), self._h())

    def _paint_viewport(self, painter: QPainter) -> None:
        painter.setPen(QPen(QColor(self.theme["viewport"]), 1))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        x = self._x(self._offset)
        w = max(8.0, self._visible / self._total * self.width())
        painter.drawRect(int(x), 0, int(w), self._h() - 1)
