"""Main ab1kit window: open a .ab1 file and render its chromatogram."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QPointF, Qt
from PySide6.QtGui import QAction, QColor, QFont, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import QApplication, QFileDialog, QMainWindow, QWidget

from chrom_core import Chromatogram, read_ab1

BASE_COLORS = {
    "A": QColor("#1e9e1e"),  # green
    "C": QColor("#1478ff"),  # blue
    "G": QColor("#101010"),  # black
    "T": QColor("#e02020"),  # red
}
LOW_Q_COLOR = QColor("#e02020")
MID_Q_COLOR = QColor("#e0a000")
HIGH_Q_COLOR = QColor("#1e9e1e")

RULE_H = 20          # base-number ruler at the top
LETTER_H = 24        # called-base letters below the traces
Q_STRIP_H = 56       # per-base quality bars at the bottom


class ChromatogramWidget(QWidget):
    """Paints trace polylines, called bases and per-base quality bars.

    The horizontal axis is trace points; zoom is pixels per trace point.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.chrom: Chromatogram | None = None
        self.zoom = 8.0        # pixels per trace point
        self.offset = 0.0      # leftmost visible trace point
        self._avg_spacing = 12.0  # trace points between called bases
        self.setMinimumSize(400, 300)
        self.setMouseTracking(True)
        self._drag_x: float | None = None

    def set_chromatogram(self, chrom: Chromatogram) -> None:
        self.chrom = chrom
        self.offset = 0.0
        self._avg_spacing = self._total_points() / max(1, len(chrom.seq))
        self.update()

    # ------------------------------------------------------------------ view
    def _total_points(self) -> int:
        if not self.chrom:
            return 0
        return max(len(t) for t in self.chrom.traces.values())

    def _visible_points(self) -> float:
        return self.width() / self.zoom

    def _x(self, point: float) -> float:
        return (point - self.offset) * self.zoom

    # ---------------------------------------------------------- interaction
    def wheelEvent(self, event) -> None:  # noqa: N802 (Qt naming)
        if self.chrom is None:
            return
        anchor_ratio = event.position().x() / max(self.width(), 1)
        anchor = self.offset + anchor_ratio * self._visible_points()
        factor = 1.25 if event.angleDelta().y() > 0 else 1 / 1.25
        self.zoom = min(64.0, max(1.0, self.zoom * factor))
        self.offset = max(0.0, anchor - anchor_ratio * self._visible_points())
        self.update()

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self._drag_x = event.position().x()

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        if self._drag_x is None:
            return
        dx = event.position().x() - self._drag_x
        self._drag_x = event.position().x()
        self.offset = max(0.0, self.offset - dx / self.zoom)
        self.update()

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        self._drag_x = None

    # -------------------------------------------------------------- painting
    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        w, h = self.width(), self.height()
        painter.fillRect(0, 0, w, h, QColor("#ffffff"))
        if self.chrom is None:
            painter.setPen(QColor("#808080"))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Open a .ab1 file to begin")
            return

        chrom = self.chrom
        trace_top = RULE_H + 2
        trace_bottom = h - Q_STRIP_H - LETTER_H - 4
        trace_h = max(10, trace_bottom - trace_top)

        self._paint_ruler(painter, w)
        self._paint_traces(painter, w, trace_top, trace_bottom)
        self._paint_bases(painter, w, h)
        painter.end()

    def _paint_ruler(self, painter: QPainter, w: int) -> None:
        chrom = self.chrom
        assert chrom is not None
        painter.setPen(QColor("#555555"))
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
            values = chrom.traces[base]
            poly = QPolygonF()
            for j in range(start, end):
                y = bottom - 2 - values[j] / max_val * usable
                poly.append(QPointF(self._x(j), y))
            painter.setPen(QPen(BASE_COLORS[base], 1))
            painter.drawPolyline(poly)

    def _paint_bases(self, painter: QPainter, w: int, h: int) -> None:
        chrom = self.chrom
        assert chrom is not None
        letter_y = h - Q_STRIP_H - 6
        q_bottom = h - 2
        max_q = max(60, max(chrom.quality) if chrom.quality else 60)
        bar_w = max(2, min(int(self.zoom * 6), 16))
        font = QFont(self.font())
        font.setPixelSize(15)
        painter.setFont(font)
        prev_peak = 0.0
        for i, (base, q) in enumerate(zip(chrom.seq, chrom.quality)):
            peak = chrom.peak_indices[i]
            x = self._x(peak)
            next_peak = chrom.peak_indices[i + 1] if i + 1 < len(chrom.peak_indices) else peak + 10
            if x > w + 30:
                break
            if x < -30 or i and x - self._x(prev_peak) < 4:
                prev_peak = peak
                continue
            prev_peak = peak
            color = LOW_Q_COLOR if q < 10 else MID_Q_COLOR if q < 20 else HIGH_Q_COLOR
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(color)
            bar_h = int((Q_STRIP_H - 6) * q / max_q)
            painter.drawRect(int(x), q_bottom - bar_h, bar_w, bar_h)
            painter.setPen(BASE_COLORS.get(base.upper(), QColor("#404040")))
            painter.drawText(QPointF(x - 4, letter_y), base)


class MainWindow(QMainWindow):
    """Top-level ab1kit window: menu + status bar around the chromatogram."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("ab1kit")
        self.resize(1100, 700)
        self.chrom_view = ChromatogramWidget()
        self.setCentralWidget(self.chrom_view)

        file_menu = self.menuBar().addMenu("&File")
        open_action = QAction("&Open .ab1…", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self._open_file)
        file_menu.addAction(open_action)
        quit_action = QAction("E&xit", self)
        quit_action.setShortcut("Ctrl+Q")
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)

    def _open_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open chromatogram",
            "",
            "Sanger chromatograms (*.ab1);;All files (*)",
        )
        if not path:
            return
        chrom = read_ab1(Path(path))
        self.chrom_view.set_chromatogram(chrom)
        self.setWindowTitle(f"ab1kit — {Path(path).name}")
        n = len(chrom.seq)
        mean_q = sum(chrom.quality) / max(n, 1)
        self.statusBar().showMessage(f"{chrom.sample_id}: {n} bases, mean Q {mean_q:.1f}")


def run(argv: list[str] | None = None) -> int:
    """Application entry point used by the ``ab1kit`` console script."""
    app = QApplication(argv if argv is not None else [])
    window = MainWindow()
    window.show()
    return app.exec()
