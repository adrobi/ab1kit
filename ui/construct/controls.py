"""Construct-screen controls: drop zones of the empty state and the row of
region / parameter inputs above the matrix."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QDragEnterEvent, QDropEvent, QMouseEvent
from PySide6.QtWidgets import (
    QCheckBox,
    QDoubleSpinBox,
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from chrom_core.clones import Region, VerifyParams
from i18n import tr
from ui.theme import theme as palette


def dropped_paths(event: QDropEvent | QDragEnterEvent) -> list[Path]:
    return [Path(u.toLocalFile()) for u in event.mimeData().urls() if u.toLocalFile()]


class DropZone(QFrame):
    """Big dashed area: accepts dropped files/folders and works as a button."""

    clicked = Signal()
    dropped = Signal(list)       # list[Path]

    def __init__(self, text: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setAcceptDrops(True)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setMinimumSize(280, 180)
        self._prompt = text
        self._done = ""
        self._hover = False
        self._theme = "light"
        self._label = QLabel(text)
        self._label.setWordWrap(True)
        self._label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._detail = QLabel()
        self._detail.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._detail.setWordWrap(True)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(24, 24, 24, 24)
        layout.addStretch(1)
        layout.addWidget(self._label)
        layout.addWidget(self._detail)
        self._extra = QVBoxLayout()
        layout.addLayout(self._extra)
        layout.addStretch(1)
        self._restyle()

    def text(self) -> str:
        return self._label.text()

    def add_widget(self, widget: QWidget) -> None:
        """Extra control under the texts (e.g. a secondary "choose files" button)."""
        self._extra.addWidget(widget, 0, Qt.AlignmentFlag.AlignHCenter)

    def detail_text(self) -> str:
        return self._detail.text()

    def set_done(self, text: str) -> None:
        """Mark the zone as filled (``text`` names what was loaded); "" resets it."""
        self._done = text
        self._detail.setText(f"✓ {text}" if text else tr("workspace.click_to_choose"))
        self._restyle()

    def set_theme(self, name: str) -> None:
        self._theme = name
        self._restyle()

    def _restyle(self) -> None:
        colors = palette(self._theme)
        accent = colors["region"]
        border = accent if (self._hover or self._done) else colors["dim"]
        self.setStyleSheet(
            f"DropZone {{ border: 2px dashed {border}; border-radius: 12px;"
            f" background-color: {colors['panel']}; }}"
        )
        self._label.setStyleSheet(f"font-size: 15pt; color: {colors['text']};")
        self._detail.setStyleSheet(f"color: {accent if self._done else colors['dim']};")
        if not self._detail.text():
            self._detail.setText(tr("workspace.click_to_choose"))

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton and self.rect().contains(
                event.position().toPoint()):
            self.clicked.emit()

    def dragEnterEvent(self, event: QDragEnterEvent) -> None:  # noqa: N802
        if dropped_paths(event):
            self._hover = True
            self._restyle()
            event.acceptProposedAction()

    def dragLeaveEvent(self, event) -> None:  # noqa: N802
        self._hover = False
        self._restyle()

    def dropEvent(self, event: QDropEvent) -> None:  # noqa: N802
        self._hover = False
        self._restyle()
        paths = dropped_paths(event)
        if paths:
            event.acceptProposedAction()
            self.dropped.emit(paths)


class EmptyState(QWidget):
    """Shown until both a construct and reads are loaded."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.construct_zone = DropZone(tr("workspace.drop_construct"))
        self.reads_zone = DropZone(tr("workspace.drop_reads"))
        self.reads_files = QPushButton(tr("workspace.choose_files"))
        self.reads_files.setFlat(True)
        self.reads_files.setCursor(Qt.CursorShape.PointingHandCursor)
        title = QLabel(tr("workspace.empty_title"))
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title.setStyleSheet("font-size: 18pt; font-weight: 600;")
        hint = QLabel(tr("workspace.empty_hint"))
        hint.setAlignment(Qt.AlignmentFlag.AlignCenter)
        hint.setWordWrap(True)
        self._hint = hint

        zones = QHBoxLayout()
        zones.setSpacing(24)
        zones.addStretch(1)
        zones.addWidget(self.construct_zone, 3)
        self.reads_zone.add_widget(self.reads_files)
        zones.addWidget(self.reads_zone, 3)
        zones.addStretch(1)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(24, 24, 24, 24)
        layout.addStretch(1)
        layout.addWidget(title)
        layout.addWidget(hint)
        layout.addSpacing(16)
        layout.addLayout(zones, 2)
        layout.addStretch(1)

    def zones(self) -> list[DropZone]:
        return [self.construct_zone, self.reads_zone]

    def set_theme(self, name: str) -> None:
        for zone in self.zones():
            zone.set_theme(name)
        self._hint.setStyleSheet(f"color: {palette(name)['dim']};")


class ControlsBar(QWidget):
    """Region (1-based, inclusive), circular flag, thresholds and data buttons.

    Spin boxes report only finished input (no keyboard tracking). Programmatic
    updates (``set_region`` etc.) emit nothing.
    """

    region_edited = Signal(object)     # Region (0-based half-open)
    whole_requested = Signal()
    params_changed = Signal()
    circular_toggled = Signal(bool)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.start_spin = self._spin(1, 1)
        self.end_spin = self._spin(1, 1)
        self.whole_button = QPushButton(tr("workspace.whole"))
        self.whole_button.setToolTip(tr("workspace.whole_tip"))
        self.circular = QCheckBox(tr("workspace.circular"))
        self.circular.setToolTip(tr("workspace.circular_tip"))
        defaults = VerifyParams()
        self.min_q = self._spin(0, 60, defaults.min_q)
        self.min_q.setToolTip(tr("workspace.min_q_tip"))
        self.het_ratio = QDoubleSpinBox()
        self.het_ratio.setRange(0.05, 0.95)
        self.het_ratio.setSingleStep(0.05)
        self.het_ratio.setDecimals(2)
        self.het_ratio.setValue(defaults.het_ratio)
        self.het_ratio.setKeyboardTracking(False)
        self.het_ratio.setToolTip(tr("workspace.het_ratio_tip"))
        self.mixed_min = self._spin(1, 99, defaults.mixed_min_positions)
        self.mixed_min.setToolTip(tr("workspace.mixed_min_tip"))
        self.grouping_button = QPushButton(tr("workspace.grouping"))
        self.construct_button = QPushButton(tr("workspace.construct_button"))
        self.reads_button = QPushButton(tr("workspace.add_reads"))

        layout = QHBoxLayout(self)
        layout.setContentsMargins(8, 2, 8, 2)
        layout.setSpacing(6)
        for widget in (QLabel(tr("workspace.region_from")), self.start_spin,
                       QLabel(tr("workspace.region_to")), self.end_spin, self.whole_button):
            layout.addWidget(widget)
        layout.addSpacing(10)
        layout.addWidget(self.circular)
        layout.addSpacing(10)
        for text, widget in ((tr("workspace.min_q"), self.min_q),
                             (tr("workspace.het_ratio"), self.het_ratio),
                             (tr("workspace.mixed_min"), self.mixed_min)):
            layout.addWidget(QLabel(text))
            layout.addWidget(widget)
        layout.addStretch(1)
        # data buttons live in their own strip (placed by the workspace next to
        # the summary) so the parameter row stays narrow
        self.data_buttons = QWidget()
        buttons = QHBoxLayout(self.data_buttons)
        buttons.setContentsMargins(0, 0, 8, 0)
        for button in (self.grouping_button, self.reads_button, self.construct_button):
            buttons.addWidget(button)

        self.start_spin.valueChanged.connect(self._region_from_spins)
        self.end_spin.valueChanged.connect(self._region_from_spins)
        self.whole_button.clicked.connect(self.whole_requested)
        self.circular.toggled.connect(self.circular_toggled)
        for spin in (self.min_q, self.het_ratio, self.mixed_min):
            spin.valueChanged.connect(lambda _v: self.params_changed.emit())

    @staticmethod
    def _spin(lo: int, hi: int, value: int | None = None) -> QSpinBox:
        spin = QSpinBox()
        spin.setRange(lo, hi)
        spin.setKeyboardTracking(False)
        if value is not None:
            spin.setValue(value)
        return spin

    # ------------------------------------------------------------- public API
    def set_region(self, region: Region | None, length: int) -> None:
        """Show ``region`` (``None``: 1..length) without emitting ``region_edited``."""
        start, end = (0, length) if region is None else region
        for spin, value in ((self.start_spin, start + 1), (self.end_spin, end or length)):
            spin.blockSignals(True)
            spin.setRange(1, max(1, length))
            spin.setValue(value)
            spin.blockSignals(False)

    def set_circular(self, flag: bool) -> None:
        self.circular.blockSignals(True)
        self.circular.setChecked(flag)
        self.circular.blockSignals(False)

    def params(self) -> VerifyParams:
        return VerifyParams(self.min_q.value(), round(self.het_ratio.value(), 2),
                            self.mixed_min.value())

    def set_params(self, params: VerifyParams) -> None:
        for spin, value in ((self.min_q, params.min_q), (self.het_ratio, params.het_ratio),
                            (self.mixed_min, params.mixed_min_positions)):
            spin.blockSignals(True)
            spin.setValue(value)
            spin.blockSignals(False)

    def _region_from_spins(self) -> None:
        self.region_edited.emit((self.start_spin.value() - 1, self.end_spin.value()))
