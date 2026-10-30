"""One-line summary of a construct check: name, length, topology, region and
per-verdict clone counters that double as row filters for the matrix."""

from __future__ import annotations

from collections import Counter

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QHBoxLayout, QLabel, QPushButton, QWidget

from chrom_core.clones import CloneReport, Region
from chrom_core.reference import Construct
from i18n import tr
from ui import status
from ui.construct import VERDICTS
from ui.theme import theme as palette


def region_text(region: Region | None) -> str:
    """``region 7–647`` (1-based, inclusive); ``None`` is the whole construct."""
    if region is None:
        return tr("construct.region_all")
    start, end = region
    return tr("construct.region", start=start + 1, end=end)


def _thousands(n: int) -> str:
    return f"{n:,}".replace(",", " ")


class SummaryBar(QWidget):
    """Construct name and facts, then a counter chip per verdict.

    Clicking a counter toggles that verdict in the shown set and emits
    ``filter_changed`` with the verdicts still shown.
    """

    filter_changed = Signal(set)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._theme = status.current_theme()
        self._shown = set(VERDICTS)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(8, 4, 8, 4)
        layout.setSpacing(8)
        self._title = QLabel()
        self._title.setStyleSheet("font-weight: 600;")
        self._info = QLabel()
        layout.addWidget(self._title)
        layout.addWidget(self._info)
        self._counters: dict[str, QPushButton] = {}
        for verdict in VERDICTS:
            button = QPushButton()
            button.setCheckable(True)
            button.setChecked(True)
            button.setFlat(True)
            button.setToolTip(tr("construct.filter_tip"))
            button.toggled.connect(lambda on, v=verdict: self._toggle(v, on))
            button.hide()
            layout.addWidget(button)
            self._counters[verdict] = button
        layout.addStretch(1)
        self._restyle()

    # ------------------------------------------------------------- public API
    def set_data(self, construct: Construct, region: Region | None,
                 reports: list[CloneReport]) -> None:
        self._title.setText(construct.name)
        topology = tr("construct.circular" if construct.circular else "construct.linear")
        self._info.setText(" · ".join((
            f"{tr('construct.bp', n=_thousands(len(construct.seq)))}, {topology}",
            region_text(region),
        )))
        counts = Counter(r.verdict for r in reports)
        for verdict, button in self._counters.items():
            n = counts.get(verdict, 0)
            button.setText(tr(f"construct.count.{verdict}", n=n))
            button.setVisible(n > 0)

    def counter(self, verdict: str) -> QPushButton:
        return self._counters[verdict]

    def shown_verdicts(self) -> set[str]:
        return set(self._shown)

    def title_text(self) -> str:
        return self._title.text()

    def info_text(self) -> str:
        return self._info.text()

    def set_theme(self, name: str) -> None:
        self._theme = "dark" if name == "dark" else "light"
        self._restyle()

    # --------------------------------------------------------------- private
    def _toggle(self, verdict: str, on: bool) -> None:
        if on:
            self._shown.add(verdict)
        else:
            self._shown.discard(verdict)
        self._restyle()
        self.filter_changed.emit(set(self._shown))

    def _restyle(self) -> None:
        dim = palette(self._theme)["dim"]
        for verdict, button in self._counters.items():
            bg = status.status_color(verdict, self._theme).name()
            fg = status.status_text_color(verdict, self._theme).name()
            # shown: filled chip; hidden: outline only, struck through
            button.setStyleSheet(
                "QPushButton { border-radius: 9px; padding: 1px 10px; font-weight: 600;"
                f" border: 1px solid {bg}; background-color: transparent; color: {dim};"
                " text-decoration: line-through; }"
                f"QPushButton:checked {{ background-color: {bg}; color: {fg};"
                " text-decoration: none; }"
            )
