"""Clone matrix: every clone as a row on one shared coordinate axis under the
reference track (layout "B" of the construct screen).

Colours
-------
- Coverage (per region position, reads with Q >= min_q): both strands -> theme
  ``cov_both`` (saturated), one strand -> ``cov_one`` (light), none -> the empty
  ``track`` background. Positions outside the region carry no coverage.
- Difference marks (in-region differences only) use the saturated variant of
  the status palette (``STATUS_HEX`` dark tone) in both themes, so a thin mark
  stays visible on light coverage: ``confirmed`` -> errors red, ``disputed`` ->
  red/pale-red horizontal stripes, ``low_quality`` -> no-data grey; mixed
  (double-peak) positions -> mixed purple.
- Verdict chips use ``status_color`` of the current theme.

The axis is linear 0..L also for circular constructs; anything crossing the
origin (``start > end``) is drawn as two segments.
"""

from __future__ import annotations

import math

from PySide6.QtCore import QPoint, QPointF, QRect, QRectF, QSize, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontMetrics, QMouseEvent, QPainter, QPen, QWheelEvent
from PySide6.QtWidgets import QWidget

from chrom_core.clones import CloneReport, Region
from chrom_core.reference import Construct
from i18n import tr
from ui import status
from ui.construct import VERDICTS
from ui.theme import theme as palette

MARGIN = 6
RULER_H = 16
REF_ROW_H = 28
ROW_H = 24
TRACK_H = 10        # clone coverage bar
REF_TRACK_H = 14    # reference / feature bar
MARK_H = 16
LEGEND_H = 20
CHIP_W = 84
CHIP_H = 16
GAP = 6
NAME_MIN, NAME_MAX = 44, 160
ZOOM_STEP = 1.25    # per wheel notch
MIN_SPAN = 20.0     # bp visible at most zoom
DRAG_SLOP = 3       # px before a press becomes a pan
HIT_SLOP = 4        # px around a mark that still hits it

_STATUS_ORDER = ("low_quality", "disputed", "confirmed")   # drawn bottom -> top


def segments(start: int, end: int, length: int) -> list[tuple[int, int]]:
    """Linear pieces of a half-open interval; ``start > end`` wraps the origin."""
    if start <= end:
        return [(start, end)]
    return [(start, length), (0, end)]


def mark_colors(theme: str) -> dict[str, QColor]:
    """Mark colour per difference status, plus ``mixed`` and the pale stripe of
    ``disputed`` (``disputed_stripe``)."""
    del theme  # same saturated tones in both themes (see module docstring)
    return {
        "confirmed": status.status_color("errors", "dark"),
        "disputed": status.status_color("errors", "dark"),
        "disputed_stripe": status.status_color("errors", "light"),
        "low_quality": status.status_color("no_data", "dark"),
        "mixed": status.status_color("mixed", "dark"),
    }


def _coverage_runs(report: CloneReport) -> list[tuple[int, int, int]]:
    """``(start, end, level)`` runs of consecutive positions; level 2 = both
    strands, 1 = one strand (positions without coverage are left out)."""
    runs: list[list[int]] = []
    for pos in sorted(report.coverage):
        plus, minus = report.coverage[pos]
        level = 2 if plus and minus else 1 if plus + minus else 0
        if level == 0:
            continue
        if runs and runs[-1][1] == pos and runs[-1][2] == level:
            runs[-1][1] = pos + 1
        else:
            runs.append([pos, pos + 1, level])
    return [tuple(r) for r in runs]


class _Row:
    def __init__(self, report: CloneReport) -> None:
        self.report = report
        self.runs = _coverage_runs(report)
        self.differences = sorted(   # stronger statuses painted last, on top
            report.region_differences,
            key=lambda d: _STATUS_ORDER.index(d.status) if d.status in _STATUS_ORDER else 0)
        self.mixed = report.mixed_positions


class CloneMatrix(QWidget):
    """Reference track plus one row per clone; custom painted.

    Wheel zooms the shared axis around the cursor, dragging pans it, a Shift+drag
    on the reference track picks a region. A click on a mark emits
    ``difference_clicked(clone, ref_pos)``, elsewhere on a row ``clone_clicked``.
    """

    difference_clicked = Signal(str, int)
    clone_clicked = Signal(str)
    region_changed = Signal(object)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._theme = status.current_theme()
        self._construct: Construct | None = None
        self._region: Region | None = None
        self._rows: list[_Row] = []
        self._filter = set(VERDICTS)
        self._view_start = 0.0
        self._view_span = 1.0
        self._press: QPoint | None = None
        self._press_view = 0.0
        self._mode = ""            # "", "click", "pan", "region"
        self._band: tuple[float, float] | None = None
        self.setToolTip(tr("matrix.hint"))
        self.setMinimumWidth(320)

    # ------------------------------------------------------------- public API
    def set_data(self, construct: Construct, region: Region | None,
                 reports: list[CloneReport]) -> None:
        self._construct = construct
        self._region = region
        self._rows = [_Row(r) for r in reports]
        self._view_start, self._view_span = 0.0, float(max(1, len(construct.seq)))
        self.updateGeometry()
        self.update()

    def set_filter(self, verdicts: set[str]) -> None:
        self._filter = set(verdicts)
        self.updateGeometry()
        self.update()

    def set_theme(self, name: str) -> None:
        self._theme = "dark" if name == "dark" else "light"
        self.update()

    def visible_clones(self) -> list[str]:
        return [row.report.clone for row in self._visible()]

    def view_range(self) -> tuple[float, float]:
        """``(first position, number of positions)`` shown on the axis."""
        return (self._view_start, self._view_span)

    def x_for(self, pos: float) -> float:
        left, width = self._track_left(), self._track_width()
        return left + (pos - self._view_start) / self._view_span * width

    def pos_at(self, x: float) -> float:
        return self._view_start + (x - self._track_left()) / self._track_width() * self._view_span

    def reference_rect(self) -> QRect:
        return QRect(round(self._track_left()), MARGIN + RULER_H,
                     round(self._track_width()), REF_ROW_H)

    def row_rect(self, clone: str) -> QRect:
        """Full-width rectangle of a shown clone's row; null when hidden or unknown."""
        for i, row in enumerate(self._visible()):
            if row.report.clone == clone:
                return QRect(0, self._rows_top() + i * ROW_H, self.width(), ROW_H)
        return QRect()

    def chip_rect(self, clone: str) -> QRect:
        row = self.row_rect(clone)
        if row.isNull():
            return QRect()
        return QRect(self.width() - MARGIN - CHIP_W, row.center().y() - CHIP_H // 2 + 1,
                     CHIP_W, CHIP_H)

    def mark_point(self, clone: str, ref_pos: int) -> QPoint:
        """Centre of the mark drawn for ``ref_pos`` in ``clone``'s row."""
        data = self._row(clone)
        diffs = data.differences if data is not None else []
        kind = next((d.kind for d in diffs if d.ref_pos == ref_pos), "substitution")
        return QPoint(round(self._mark_x(ref_pos, kind)), self.row_rect(clone).center().y())

    def hit_test(self, pos: QPoint) -> tuple[str, int] | None:
        """``(clone, ref_pos)`` of the difference or mixed mark under ``pos``."""
        index = self._row_index_at(pos.y())
        if index is None or not self._track_left() - HIT_SLOP <= pos.x() <= self._track_right() + HIT_SLOP:
            return None
        row = self._visible()[index]
        slop = max(HIT_SLOP, self._px_per_base() / 2 + 1)
        best: tuple[float, int, int] | None = None   # (distance, priority, ref_pos)
        candidates = [(d.ref_pos, d.kind, 0) for d in row.differences]
        candidates += [(p, "substitution", 1) for p in row.mixed]
        for ref_pos, kind, priority in candidates:
            distance = abs(self._mark_x(ref_pos, kind) - pos.x())
            if distance <= slop and (best is None or (distance, priority) < best[:2]):
                best = (distance, priority, ref_pos)
        return (row.report.clone, best[2]) if best else None

    def sizeHint(self) -> QSize:  # noqa: N802
        return QSize(720, self._content_height())

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        return QSize(320, self._content_height())

    # --------------------------------------------------------------- geometry
    def _visible(self) -> list[_Row]:
        return [r for r in self._rows if r.report.verdict in self._filter]

    def _row(self, clone: str) -> _Row | None:
        return next((r for r in self._rows if r.report.clone == clone), None)

    def _length(self) -> int:
        return len(self._construct.seq) if self._construct else 1

    def _name_width(self) -> int:
        fm = QFontMetrics(self._bold_font())
        widest = max((fm.horizontalAdvance(r.report.clone) for r in self._rows), default=0)
        return min(NAME_MAX, max(NAME_MIN, widest + GAP))

    def _track_left(self) -> float:
        return float(MARGIN + self._name_width())

    def _track_right(self) -> float:
        return float(self.width() - MARGIN - CHIP_W - GAP)

    def _track_width(self) -> float:
        return max(1.0, self._track_right() - self._track_left())

    def _px_per_base(self) -> float:
        return self._track_width() / self._view_span

    def _rows_top(self) -> int:
        return MARGIN + RULER_H + REF_ROW_H

    def _content_height(self) -> int:
        rows = max(1, len(self._visible()))
        return self._rows_top() + rows * ROW_H + GAP + LEGEND_H + MARGIN

    def _row_index_at(self, y: int) -> int | None:
        offset = y - self._rows_top()
        if offset < 0:
            return None
        index = offset // ROW_H
        return index if index < len(self._visible()) else None

    def _mark_x(self, ref_pos: int, kind: str) -> float:
        # an insertion sits before ``ref_pos``: on the boundary, not the base centre
        return self.x_for(ref_pos if kind == "insertion" else ref_pos + 0.5)

    def _bold_font(self) -> QFont:
        font = QFont(self.font())
        font.setBold(True)
        return font

    def _small_font(self) -> QFont:
        font = QFont(self.font())
        font.setPointSizeF(max(6.0, font.pointSizeF() - 1.5))
        return font

    def _set_view(self, start: float, span: float) -> None:
        length = float(self._length())
        span = min(length, max(min(MIN_SPAN, length), span))
        start = min(max(0.0, start), length - span)
        if (start, span) != (self._view_start, self._view_span):
            self._view_start, self._view_span = start, span
            self.update()

    # ------------------------------------------------------------------ input
    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802
        if self._construct is None:
            return
        steps = event.angleDelta().y() / 120.0
        if not steps:
            return
        x = event.position().x()
        if not self._track_left() <= x <= self._track_right():
            x = (self._track_left() + self._track_right()) / 2
        anchor = self.pos_at(x)
        fraction = (x - self._track_left()) / self._track_width()
        span = self._view_span / ZOOM_STEP ** steps
        length = float(self._length())
        span = min(length, max(min(MIN_SPAN, length), span))
        self._set_view(anchor - fraction * span, span)
        event.accept()

    def mousePressEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if event.button() != Qt.MouseButton.LeftButton or self._construct is None:
            return
        point = event.position().toPoint()
        self._press, self._press_view = point, self._view_start
        on_reference = self.reference_rect().adjusted(-HIT_SLOP, 0, HIT_SLOP, 0).contains(point)
        if event.modifiers() & Qt.KeyboardModifier.ShiftModifier and on_reference:
            self._mode = "region"
            x = self._clamp_x(point.x())
            self._band = (x, x)
            self.update()
        else:
            self._mode = "click"

    def mouseMoveEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if self._press is None:
            return
        x = event.position().x()
        if self._mode == "region":
            self._band = (self._band[0], self._clamp_x(x))
            self.update()
            return
        dx = x - self._press.x()
        if self._mode == "click" and abs(dx) > DRAG_SLOP:
            self._mode = "pan"
        if self._mode == "pan":
            self._set_view(self._press_view - dx / self._track_width() * self._view_span,
                           self._view_span)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if self._press is None or event.button() != Qt.MouseButton.LeftButton:
            return
        mode, press = self._mode, self._press
        self._press, self._mode = None, ""
        if mode == "region":
            a, b = sorted(self._band or (0.0, 0.0))
            self._band = None
            self.update()
            start, end = round(self.pos_at(a)), round(self.pos_at(b))
            start, end = max(0, start), min(self._length(), end)
            if end > start:
                self.region_changed.emit((start, end))
        elif mode == "click":
            hit = self.hit_test(press)
            if hit is not None:
                self.difference_clicked.emit(*hit)
                return
            index = self._row_index_at(press.y())
            if index is not None:
                self.clone_clicked.emit(self._visible()[index].report.clone)

    def _clamp_x(self, x: float) -> float:
        lo, hi = self.x_for(0), self.x_for(self._length())
        return min(max(x, max(lo, self._track_left())), min(hi, self._track_right()))

    # --------------------------------------------------------------- painting
    def paintEvent(self, event) -> None:  # noqa: N802, ARG002
        colors = palette(self._theme)
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor(colors["bg"]))
        if self._construct is None:
            self._paint_empty(painter, colors)
            painter.end()
            return
        self._paint_ruler(painter, colors)
        self._paint_reference(painter, colors)
        rows = self._visible()
        for i, row in enumerate(rows):
            self._paint_row(painter, colors, row, self._rows_top() + i * ROW_H)
        if not rows:
            self._paint_empty(painter, colors)
        self._paint_region_lines(painter, colors, len(rows))
        self._paint_legend(painter, colors, max(1, len(rows)))
        painter.end()

    def _track_clip(self) -> QRectF:
        return QRectF(self._track_left(), 0, self._track_width(), self.height())

    def _spans(self, start: int, end: int) -> list[QRectF]:
        """Pixel x-extents (as zero-height rects) of an interval, split at the origin."""
        return [QRectF(self.x_for(a), 0, self.x_for(b) - self.x_for(a), 0)
                for a, b in segments(start, end, self._length())]

    def _paint_empty(self, painter: QPainter, colors: dict) -> None:
        painter.setPen(QColor(colors["dim"]))
        rect = QRect(0, self._rows_top(), self.width(), ROW_H)
        painter.drawText(rect, Qt.AlignmentFlag.AlignCenter, tr("matrix.empty"))

    def _paint_ruler(self, painter: QPainter, colors: dict) -> None:
        painter.setFont(self._small_font())
        fm = painter.fontMetrics()
        label_w = fm.horizontalAdvance(f"{self._length():,}") + 16
        raw = self._view_span * label_w / self._track_width()
        magnitude = 10 ** math.floor(math.log10(max(raw, 1)))
        step = next(m * magnitude for m in (1, 2, 5, 10) if m * magnitude >= raw)
        y = MARGIN + RULER_H
        painter.setPen(QPen(QColor(colors["dim"]), 1))
        painter.save()
        painter.setClipRect(self._track_clip().adjusted(-label_w, 0, label_w, 0))
        first = max(step, math.ceil(self._view_start / step) * step)
        v = first
        while v <= self._view_start + self._view_span:
            x = self.x_for(v)
            painter.drawLine(QPointF(x, y - 4), QPointF(x, y))
            painter.drawText(QRectF(x - label_w / 2, MARGIN - 2, label_w, RULER_H - 3),
                             Qt.AlignmentFlag.AlignCenter, f"{int(v):,}".replace(",", " "))
            v += step
        painter.restore()
        painter.drawLine(QPointF(self._track_left(), y), QPointF(self._track_right(), y))

    def _paint_reference(self, painter: QPainter, colors: dict) -> None:
        top = MARGIN + RULER_H
        cy = top + REF_ROW_H / 2
        track = QRectF(self._track_left(), cy - REF_TRACK_H / 2, self._track_width(), REF_TRACK_H)
        painter.setFont(self._small_font())
        painter.setPen(QColor(colors["dim"]))
        painter.drawText(QRectF(self.width() - MARGIN - CHIP_W, top, CHIP_W, REF_ROW_H),
                         Qt.AlignmentFlag.AlignCenter, tr("matrix.reference"))
        painter.save()
        painter.setClipRect(self._track_clip())
        self._fill_axis(painter, track, 0, self._length(), QColor(colors["track"]))
        seen = set()
        for feature in self._construct.features:
            key = (feature.name, feature.start, feature.end)
            if feature.kind == "source" or key in seen:
                continue
            seen.add(key)
            colour = QColor(colors["feature_cds" if feature.kind in ("CDS", "gene") else "feature"])
            for span in self._spans(feature.start, feature.end):
                rect = QRectF(span.left(), track.top(), span.width(), track.height())
                painter.fillRect(rect, colour)
                if rect.width() > 24:
                    painter.setPen(QColor("#ffffff"))
                    text = painter.fontMetrics().elidedText(
                        feature.name, Qt.TextElideMode.ElideRight, int(rect.width()) - 6)
                    visible = rect.intersected(self._track_clip())
                    painter.drawText(visible.adjusted(6, 0, -3, 0),
                                     Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, text)
        start, end = self._region if self._region is not None else (0, self._length())
        painter.setPen(QPen(QColor(colors["region"]), 2))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        for span in self._spans(start, end):
            painter.drawRoundedRect(QRectF(span.left(), track.top() - 3, span.width(),
                                           track.height() + 6), 3, 3)
        if self._band is not None:
            a, b = sorted(self._band)
            band = QColor(colors["region"])
            band.setAlpha(80)
            painter.fillRect(QRectF(a, track.top() - 3, b - a, track.height() + 6), band)
        painter.restore()

    def _fill_axis(self, painter: QPainter, track: QRectF, start: int, end: int,
                   colour: QColor) -> None:
        for span in self._spans(start, end):
            painter.fillRect(QRectF(span.left(), track.top(), span.width(), track.height()), colour)

    def _paint_row(self, painter: QPainter, colors: dict, row: _Row, top: int) -> None:
        report = row.report
        cy = top + ROW_H / 2
        painter.setFont(self._bold_font())
        painter.setPen(QColor(colors["text"]))
        name = painter.fontMetrics().elidedText(report.clone, Qt.TextElideMode.ElideRight,
                                                self._name_width() - GAP)
        painter.drawText(QRectF(MARGIN, top, self._name_width() - GAP, ROW_H),
                         Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, name)

        track = QRectF(self._track_left(), cy - TRACK_H / 2, self._track_width(), TRACK_H)
        painter.save()
        painter.setClipRect(self._track_clip())
        self._fill_axis(painter, track, 0, self._length(), QColor(colors["track"]))
        shades = {2: QColor(colors["cov_both"]), 1: QColor(colors["cov_one"])}
        for start, end, level in row.runs:
            self._fill_axis(painter, track, start, end, shades[level])
        marks = mark_colors(self._theme)
        width = max(3.0, self._px_per_base())
        for pos in row.mixed:
            self._paint_mark(painter, self._mark_x(pos, "substitution"), cy, width, marks["mixed"])
        for d in row.differences:
            x = self._mark_x(d.ref_pos, d.kind)
            if d.status == "disputed":
                self._paint_mark(painter, x, cy, width, marks["disputed_stripe"])
                for k in range(0, MARK_H, 4):
                    painter.fillRect(QRectF(x - width / 2, cy - MARK_H / 2 + k, width, 2),
                                     marks["disputed"])
            else:
                self._paint_mark(painter, x, cy, width, marks.get(d.status, marks["low_quality"]))
        painter.restore()
        self._paint_chip(painter, report)

    @staticmethod
    def _paint_mark(painter: QPainter, x: float, cy: float, width: float, colour: QColor) -> None:
        painter.fillRect(QRectF(x - width / 2, cy - MARK_H / 2, width, MARK_H), colour)

    def _chip_text(self, report: CloneReport) -> str:
        if report.verdict == "errors":
            n = sum(1 for d in report.region_differences if d.status != "low_quality")
            return tr("matrix.chip.errors", n=n)
        if report.verdict == "incomplete":
            return tr("matrix.chip.incomplete", pct=math.floor(report.region_coverage_pct))
        return tr(f"matrix.chip.{report.verdict}") if report.verdict in VERDICTS else report.verdict

    def _paint_chip(self, painter: QPainter, report: CloneReport) -> None:
        rect = QRectF(self.chip_rect(report.clone))
        painter.save()
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(status.status_color(report.verdict, self._theme))
        painter.drawRoundedRect(rect, CHIP_H / 2, CHIP_H / 2)
        painter.setPen(status.status_text_color(report.verdict, self._theme))
        painter.setFont(self._small_font())
        text = painter.fontMetrics().elidedText(self._chip_text(report),
                                                Qt.TextElideMode.ElideRight, CHIP_W - 8)
        painter.drawText(rect, Qt.AlignmentFlag.AlignCenter, text)
        painter.restore()

    def _paint_region_lines(self, painter: QPainter, colors: dict, n_rows: int) -> None:
        """Thin dashed region boundaries through the clone rows."""
        if self._region is None or n_rows == 0:
            return
        top, bottom = self._rows_top(), self._rows_top() + n_rows * ROW_H
        pen = QPen(QColor(colors["region"]), 1, Qt.PenStyle.DashLine)
        painter.save()
        painter.setClipRect(self._track_clip())
        painter.setPen(pen)
        for x in {self.x_for(self._region[0]), self.x_for(self._region[1])}:
            painter.drawLine(QPointF(x, top), QPointF(x, bottom))
        painter.restore()

    def _paint_legend(self, painter: QPainter, colors: dict, n_rows: int) -> None:
        top = self._rows_top() + n_rows * ROW_H + GAP
        marks = mark_colors(self._theme)
        items = [
            ("block", QColor(colors["cov_both"]), "matrix.legend.cov_both"),
            ("block", QColor(colors["cov_one"]), "matrix.legend.cov_one"),
            ("block", QColor(colors["track"]), "matrix.legend.cov_none"),
            ("mark", marks["confirmed"], "matrix.legend.confirmed"),
            ("stripe", marks["disputed"], "matrix.legend.disputed"),
            ("mark", marks["low_quality"], "matrix.legend.low_quality"),
            ("mark", marks["mixed"], "matrix.legend.mixed"),
        ]
        painter.setFont(self._small_font())
        fm = painter.fontMetrics()
        x, cy = float(self._track_left()), top + LEGEND_H / 2
        for kind, colour, key in items:
            if kind == "block":
                painter.fillRect(QRectF(x, cy - 4, 14, 8), colour)
                x += 18
            else:
                if kind == "stripe":
                    painter.fillRect(QRectF(x, cy - 6, 4, 12), marks["disputed_stripe"])
                    for k in range(0, 12, 4):
                        painter.fillRect(QRectF(x, cy - 6 + k, 4, 2), colour)
                else:
                    painter.fillRect(QRectF(x, cy - 6, 4, 12), colour)
                x += 8
            text = tr(key)
            painter.setPen(QColor(colors["dim"]))
            painter.drawText(QPointF(x, cy + fm.ascent() / 2 - 1), text)
            x += fm.horizontalAdvance(text) + 12
