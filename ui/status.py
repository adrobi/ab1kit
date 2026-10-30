"""Shared visual language: status palette, quality colours, chips and icons.

Status colours come from ``chrom_core.clone_report.STATUS_HEX`` (single source of
truth, also used by the XLSX report). Quality colours come from the theme dicts;
``quality_bucket`` is the one place that defines the Q thresholds.
"""

from __future__ import annotations

from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QAction, QColor, QIcon, QPainter, QPixmap
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtWidgets import QLabel

from chrom_core.clone_report import STATUS_HEX
from ui.theme import resource_dir, theme as palette

Q_LOW = 20    # Q below this is red
Q_HIGH = 30   # Q at or above this is green; in between is yellow

_ICON_SIZES = (16, 20, 24, 32, 48)
_current_theme = "light"
_icon_cache: dict[tuple[str, str], QIcon] = {}
_bound: list[tuple[QAction, str]] = []


def quality_bucket(q: int) -> str:
    """Theme key for a phred score: ``low_q`` (<20), ``mid_q`` (20-29), ``high_q`` (>=30)."""
    return "low_q" if q < Q_LOW else "mid_q" if q < Q_HIGH else "high_q"


def quality_color(q: int, theme: str) -> QColor:
    return QColor(palette(theme)[quality_bucket(q)])


def status_color(verdict: str, theme: str) -> QColor:
    """Fill colour of a verdict (``ok/errors/incomplete/mixed/no_data``); unknown -> no_data."""
    light, dark = STATUS_HEX.get(verdict, STATUS_HEX["no_data"])
    return QColor(dark if theme == "dark" else light)


def status_text_color(verdict: str, theme: str) -> QColor:
    """Readable text colour on top of ``status_color``."""
    bg = status_color(verdict, theme)
    luminance = 0.299 * bg.red() + 0.587 * bg.green() + 0.114 * bg.blue()
    return QColor("#1a1a1a" if luminance > 150 else "#ffffff")


def current_theme() -> str:
    return _current_theme


def status_chip(verdict: str, text: str) -> QLabel:
    """Rounded pill in the current theme's status colours."""
    chip = QLabel(text)
    bg = status_color(verdict, _current_theme).name()
    fg = status_text_color(verdict, _current_theme).name()
    chip.setStyleSheet(
        f"QLabel {{ background-color: {bg}; color: {fg}; border-radius: 9px;"
        " padding: 1px 10px; font-weight: 600; }"
    )
    chip.setAlignment(Qt.AlignmentFlag.AlignCenter)
    return chip


def icon(name: str) -> QIcon:
    """Lucide icon ``ui/icons/<name>.svg`` tinted with the current theme's text colour."""
    key = (name, _current_theme)
    cached = _icon_cache.get(key)
    if cached is not None:
        return cached
    svg = (resource_dir() / "icons" / f"{name}.svg").read_text(encoding="utf-8")
    svg = svg.replace("currentColor", palette(_current_theme)["text"])
    renderer = QSvgRenderer(svg.encode("utf-8"))
    result = QIcon()
    for size in _ICON_SIZES:
        pixmap = QPixmap(size, size)
        pixmap.fill(Qt.GlobalColor.transparent)
        painter = QPainter(pixmap)
        renderer.render(painter, QRectF(0, 0, size, size))
        painter.end()
        result.addPixmap(pixmap)
    _icon_cache[key] = result
    return result


def bind_icon(action: QAction, name: str) -> None:
    """Give ``action`` the icon and keep it tinted when the theme changes."""
    action.setIcon(icon(name))
    _bound.append((action, name))


def set_theme(name: str) -> None:
    """Switch the icon tint and re-icon every bound action."""
    global _current_theme
    _current_theme = "dark" if name == "dark" else "light"
    alive = []
    for action, icon_name in _bound:
        try:
            action.setIcon(icon(icon_name))
        except RuntimeError:  # underlying QAction already deleted
            continue
        alive.append((action, icon_name))
    _bound[:] = alive
