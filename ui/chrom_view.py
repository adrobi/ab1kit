"""The chromatogram view widget: traces, bases, quality, navigation."""

from __future__ import annotations

import bisect

from PySide6.QtCore import QPointF, QRect, Qt, Signal
from PySide6.QtGui import QColor, QFont, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import QRubberBand, QToolTip, QWidget

from chrom_core import Chromatogram
from i18n import tr
from ui.theme import LIGHT

MIN_ZOOM = 0.02
MAX_ZOOM = 64.0

RULE_H = 20          # base-number ruler at the top
LETTER_H = 24        # called-base letters below the traces
Q_STRIP_H = 56       # per-base quality bars at the bottom


class ChromatogramWidget(QWidget):
    """Paints trace polylines, called bases and per-base quality bars.

    The horizontal axis is trace points; zoom is pixels per trace point.
    """

    view_changed = Signal()      # offset/zoom moved (scrollbar & minimap sync)
    cursor_info = Signal(str)    # description of the hovered base ("" on leave)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.chrom: Chromatogram | None = None
        self.zoom = 8.0        # pixels per trace point
        self.offset = 0.0      # leftmost visible trace point
        self._avg_spacing = 12.0  # trace points between called bases
        self.theme = LIGHT
        self.visible_channels: set[str] = set("GATC")
        self.smoothing = False
        self.highlights: list[tuple[int, int]] = []  # base ranges, [start, end)
        self._hover_x: float | None = None
        self._drag_x: float | None = None
        self._select_origin: float | None = None
        self._rubber = QRubberBand(QRubberBand.Shape.Rectangle, self)
        self.setMinimumSize(400, 300)
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

    # ------------------------------------------------------------- state API
    def set_chromatogram(self, chrom: Chromatogram) -> None:
        self.chrom = chrom
        self.offset = 0.0
        self.zoom = 8.0
        self.highlights = []
        self._avg_spacing = self._total_points() / max(1, len(chrom.seq))
        self.update()
        self.view_changed.emit()

    def set_theme(self, theme: dict) -> None:
        self.theme = theme
        self.update()

    def set_visible_channel(self, base: str, visible: bool) -> None:
        (self.visible_channels.add if visible else self.visible_channels.discard)(base)
        self.update()

    def set_smoothing(self, on: bool) -> None:
        self.smoothing = on
        self.update()

    def set_highlights(self, ranges: list[tuple[int, int]]) -> None:
        self.highlights = ranges
        self.update()

    # ------------------------------------------------------------ view helpers
    def _total_points(self) -> int:
        if not self.chrom:
            return 0
        return max(len(t) for t in self.chrom.traces.values())

    def _visible_points(self) -> float:
        return self.width() / self.zoom

    def _clamp_offset(self, offset: float) -> float:
        return max(0.0, min(offset, max(0.0, self._total_points() - self._visible_points())))

    def _x(self, point: float) -> float:
        return (point - self.offset) * self.zoom

    def _point_at(self, x: float) -> float:
        return self.offset + x / self.zoom

    def set_offset(self, offset: float) -> None:
        self.offset = self._clamp_offset(offset)
        self.update()
        self.view_changed.emit()

    def center_on_point(self, point: float) -> None:
        self.set_offset(point - self._visible_points() / 2)

    def goto_base(self, number: int) -> None:
        """Scroll so that 1-based base ``number`` is in the middle."""
        if not self.chrom or not (1 <= number <= len(self.chrom.seq)):
            return
        self.center_on_point(self.chrom.peak_indices[number - 1])

    def _zoom(self, factor: float, anchor_ratio: float = 0.5) -> None:
        anchor = self.offset + anchor_ratio * self._visible_points()
        self.zoom = min(MAX_ZOOM, max(MIN_ZOOM, self.zoom * factor))
        self.offset = self._clamp_offset(anchor - anchor_ratio * self._visible_points())
        self.update()
        self.view_changed.emit()

    def zoom_in(self) -> None:
        self._zoom(1.25)

    def zoom_out(self) -> None:
        self._zoom(1 / 1.25)

    def zoom_fit(self) -> None:
        if not self.chrom:
            return
        self.zoom = min(MAX_ZOOM, max(MIN_ZOOM, self.width() / self._total_points()))
        self.offset = 0.0
        self.update()
        self.view_changed.emit()

    def zoom_reset(self) -> None:
        self.zoom = 8.0
        self.offset = 0.0
        self.update()
        self.view_changed.emit()

    def zoom_to_points(self, p0: float, p1: float) -> None:
        span = max(p1 - p0, self._total_points() / 2000.0)
        self.zoom = min(MAX_ZOOM, max(MIN_ZOOM, self.width() / span))
        self.offset = self._clamp_offset(p0 - self._visible_points() * 0.02)
        self.update()
        self.view_changed.emit()

    # ----------------------------------------------------------- interaction
    def wheelEvent(self, event) -> None:  # noqa: N802 (Qt naming)
        if self.chrom is None:
            return
        anchor_ratio = event.position().x() / max(self.width(), 1)
        self._zoom(1.25 if event.angleDelta().y() > 0 else 1 / 1.25, anchor_ratio)

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if self.chrom is None:
            return
        if event.button() == Qt.MouseButton.LeftButton and event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            self._select_origin = event.position().x()
            self._rubber.setGeometry(QRect(int(self._select_origin), 0, 0, self.height()))
            self._rubber.show()
        elif event.button() == Qt.MouseButton.LeftButton:
            self._drag_x = event.position().x()

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        x = event.position().x()
        if self._select_origin is not None:
            left = int(min(self._select_origin, x))
            self._rubber.setGeometry(QRect(left, 0, int(abs(x - self._select_origin)), self.height()))
            return
        if self._drag_x is not None:
            dx = x - self._drag_x
            self._drag_x = x
            self.set_offset(self.offset - dx / self.zoom)
            return
        if self.chrom is not None:
            self._hover_x = x
            self._report_cursor(event.globalPosition().toPoint())
            self.update()

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        if self._select_origin is not None and self.chrom is not None:
            self._rubber.hide()
            x0, x1 = sorted((self._select_origin, event.position().x()))
            self._select_origin = None
            if x1 - x0 > 12:
                self.zoom_to_points(self._point_at(x0), self._point_at(x1))
        self._drag_x = None

    def leaveEvent(self, event) -> None:  # noqa: N802
        self._hover_x = None
        self.cursor_info.emit("")
        self.update()

    def keyPressEvent(self, event) -> None:  # noqa: N802
        step = self._avg_spacing
        page = self._visible_points() * 0.9
        key = event.key()
        if key in (Qt.Key.Key_Left, Qt.Key.Key_Right):
            direction = -1 if key == Qt.Key.Key_Left else 1
            self.set_offset(self.offset + direction * step)
        elif key == Qt.Key.Key_PageUp:
            self.set_offset(self.offset - page)
        elif key == Qt.Key.Key_PageDown:
            self.set_offset(self.offset + page)
        elif key == Qt.Key.Key_Home:
            self.set_offset(0.0)
        elif key == Qt.Key.Key_End:
            self.set_offset(self._total_points())
        elif key in (Qt.Key.Key_Plus, Qt.Key.Key_Equal):
            self.zoom_in()
        elif key == Qt.Key.Key_Minus:
            self.zoom_out()
        else:
            super().keyPressEvent(event)

    def _report_cursor(self, global_pos) -> None:
        point = self._point_at(self._hover_x or 0.0)
        i = self._nearest_base(point)
        if i is None:
            self.cursor_info.emit("")
            return
        base = self.chrom.seq[i]
        quality = self.chrom.quality[i]
        text = tr("cursor.info", n=i + 1, base=base, q=quality)
        self.cursor_info.emit(text)
        QToolTip.showText(global_pos, text, self)

    def _nearest_base(self, point: float) -> int | None:
        if not self.chrom:
            return None
        peaks = self.chrom.peak_indices
        i = bisect.bisect_left(peaks, point)
        best: int | None = None
        best_dist: float | None = None
        for cand in (i - 1, i):
            if 0 <= cand < len(peaks):
                dist = abs(peaks[cand] - point)
                if best_dist is None or dist < best_dist:
                    best, best_dist = cand, dist
        return best

    # -------------------------------------------------------------- painting
    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        painter.fillRect(0, 0, w, h, QColor(self.theme["bg"]))
        if self.chrom is None:
            painter.setPen(QColor(self.theme["dim"]))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, tr("view.empty"))
            return

        trace_top = RULE_H + 2
        trace_bottom = h - Q_STRIP_H - LETTER_H - 4

        self._paint_highlights(painter, w, h)
        self._paint_grid(painter, w, trace_top, trace_bottom)
        self._paint_ruler(painter, w)
        self._paint_traces(painter, w, trace_top, trace_bottom)
        self._paint_bases(painter, w, h)
        self._paint_cursor(painter, h)
        painter.end()

    def _paint_highlights(self, painter: QPainter, w: int, h: int) -> None:
        if not self.highlights or not self.chrom:
            return
        color = QColor(self.theme["highlight"])
        color.setAlpha(90)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(color)
        peaks = self.chrom.peak_indices
        for start, end in self.highlights:
            x0 = self._x(peaks[max(0, start)])
            x1 = self._x(peaks[min(len(peaks) - 1, max(0, end - 1))]) + self._avg_spacing
            if x1 < 0 or x0 > w:
                continue
            painter.drawRect(int(x0), 0, int(x1 - x0), h)

    def _paint_grid(self, painter: QPainter, w: int, top: int, bottom: int) -> None:
        if not self.chrom or self.zoom * self._avg_spacing < 7:
            return
        painter.setPen(QPen(QColor(self.theme["grid"]), 1))
        for peak in self.chrom.peak_indices:
            x = self._x(peak)
            if x < 0 or x > w:
                continue
            painter.drawLine(QPointF(x, top), QPointF(x, bottom))

    def _paint_ruler(self, painter: QPainter, w: int) -> None:
        chrom = self.chrom
        assert chrom is not None
        painter.setPen(QColor(self.theme["dim"]))
        font = QFont(self.font())
        font.setPixelSize(11)
        painter.setFont(font)
        step = next(
            (s for s in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
             if s * self._avg_spacing * self.zoom >= 100),
            1000,
        )
        for i in range(0, len(chrom.seq), step):
            x = self._x(chrom.peak_indices[i])
            if x < -30 or x > w + 30:
                continue
            painter.drawLine(QPointF(x, RULE_H - 4), QPointF(x, RULE_H))
            painter.drawText(QPointF(x + 2, 13), str(i + 1))

    @staticmethod
    def _smooth(values: list[int], k: int = 3) -> list[float]:
        out = []
        for i in range(len(values)):
            lo = max(0, i - k // 2)
            hi = min(len(values), i + k // 2 + 1)
            out.append(sum(values[lo:hi]) / (hi - lo))
        return out

    def _paint_traces(self, painter: QPainter, w: int, top: int, bottom: int) -> None:
        chrom = self.chrom
        assert chrom is not None
        start = max(0, int(self.offset) - 2)
        end = min(self._total_points(), int(self.offset + self._visible_points()) + 3)
        # normalize to the visible window, not the global max: early peaks are
        # often much taller than the rest of the read (FinchTV behaves the same)
        max_val = max(max(t[start:end]) for t in chrom.traces.values()) or 1
        usable = bottom - top - 4  # keep the tallest peak off the ruler
        for base in "GATC":
            if base not in self.visible_channels:
                continue
            values = chrom.traces[base][start:end]
            if self.smoothing:
                values = self._smooth(list(values))
            poly = QPolygonF()
            for j, value in enumerate(values):
                y = bottom - 2 - value / max_val * usable
                poly.append(QPointF(self._x(start + j), y))
            painter.setPen(QPen(QColor(self.theme["bases"][base]), 2 if self.zoom >= 4 else 1))
            painter.drawPolyline(poly)

    def _paint_bases(self, painter: QPainter, w: int, h: int) -> None:
        chrom = self.chrom
        assert chrom is not None
        letter_y = h - Q_STRIP_H - 6
        q_bottom = h - 2
        max_q = max(60, max(chrom.quality) if chrom.quality else 60)
        bar_w = max(2, min(int(self.zoom * 6), 16))
        font = QFont(self.font())
        font.setFamilies(["Consolas", "Courier New", "monospace"])
        font.setPixelSize(15)
        painter.setFont(font)
        prev_peak = 0.0
        for i, (base, q) in enumerate(zip(chrom.seq, chrom.quality)):
            peak = chrom.peak_indices[i]
            x = self._x(peak)
            if x > w + 30:
                break
            if x < -30 or (i and x - self._x(prev_peak) < 4):
                prev_peak = peak
                continue
            prev_peak = peak
            bucket = "low_q" if q < 10 else "mid_q" if q < 20 else "high_q"
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(QColor(self.theme[bucket]))
            bar_h = int((Q_STRIP_H - 6) * q / max_q)
            painter.drawRect(int(x), q_bottom - bar_h, bar_w, bar_h)
            painter.setPen(QColor(self.theme["bases"].get(base.upper(), QColor(self.theme["text"]))))
            painter.drawText(QPointF(x - 4, letter_y), base)

    def _paint_cursor(self, painter: QPainter, h: int) -> None:
        if self._hover_x is None or self.chrom is None:
            return
        painter.setPen(QPen(QColor(self.theme["cursor"]), 1, Qt.PenStyle.DashLine))
        painter.drawLine(QPointF(self._hover_x, RULE_H), QPointF(self._hover_x, h - 2))
