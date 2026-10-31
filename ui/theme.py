"""Light/dark palettes shared by all ab1kit widgets."""

from __future__ import annotations

import sys
import tempfile
from functools import lru_cache
from pathlib import Path

from PySide6.QtGui import QIcon

LIGHT = {
    "bg": "#ffffff",
    "panel": "#f3f3f3",
    "text": "#1a1a1a",
    "dim": "#808080",
    "grid": "#e8e8e8",
    "cursor": "#b04040",
    "highlight": "#ffe98a",
    "viewport": "#606060",
    "bases": {"A": "#1e9e1e", "C": "#1478ff", "G": "#333333", "T": "#e02020"},
    "low_q": "#e02020",
    "mid_q": "#e0a000",
    "high_q": "#1e9e1e",
    # construct matrix: coverage on both strands / one strand, empty track,
    # region frame, features (CDS/gene and others)
    "cov_both": "#4caf50",
    "cov_one": "#b5e0b7",
    "track": "#e8e8e8",
    "region": "#2563eb",
    "feature_cds": "#2e8b57",
    "feature": "#7a8290",
}

DARK = {
    "bg": "#1e1e1e",
    "panel": "#252526",
    "text": "#d4d4d4",
    "dim": "#8a8a8a",
    "grid": "#3a3a3a",
    "cursor": "#d47070",
    "highlight": "#5a5220",
    "viewport": "#9a9a9a",
    "bases": {"A": "#37d437", "C": "#4a9bff", "G": "#e0e0e0", "T": "#ff5252"},
    "low_q": "#ff5252",
    "mid_q": "#e0a000",
    "high_q": "#37d437",
    "cov_both": "#43a047",
    "cov_one": "#3d6e41",
    "track": "#333333",
    "region": "#5b8def",
    "feature_cds": "#2e7d4f",
    "feature": "#5c6370",
}


def theme(name: str) -> dict:
    return DARK if name == "dark" else LIGHT


def resource_dir() -> Path:
    """Folder holding style*.qss and icons/: next to this module, or under
    ``sys._MEIPASS/ui`` in a PyInstaller bundle (see ab1kit.spec datas)."""
    frozen = getattr(sys, "_MEIPASS", None)
    return Path(frozen) / "ui" if frozen else Path(__file__).resolve().parent


def app_icon() -> QIcon:
    """The application icon (``ui/icons/app/ab1kit.png``), frozen-safe."""
    return QIcon(str(resource_dir() / "icons" / "app" / "ab1kit.png"))


@lru_cache(maxsize=2)
def _read_qss(dark: bool) -> str:
    name = "style_dark.qss" if dark else "style.qss"
    return (resource_dir() / name).read_text(encoding="utf-8")


def tinted_icon(name: str, color: str, stroke: str = "2") -> str:
    """Path of icon ``name`` (``ui/icons/<name>.svg``) drawn in ``color``, for QSS
    ``image: url(...)``. Written once per icon and colour to the temp folder;
    ``""`` when that folder cannot be written (see ``image_rule``)."""
    folder = Path(tempfile.gettempdir()) / "ab1kit_icons"
    target = folder / f"{name}_{color.lstrip('#')}_{stroke}.svg"
    try:
        folder.mkdir(exist_ok=True)
        if not target.exists():
            svg = (resource_dir() / "icons" / f"{name}.svg").read_text(encoding="utf-8")
            target.write_text(svg.replace("currentColor", color)
                              .replace('stroke-width="2"', f'stroke-width="{stroke}"'),
                              encoding="utf-8")
    except OSError:
        return ""
    return target.as_posix()


def image_rule(path: str) -> str:
    """QSS ``image: url(path);``, or nothing for a missing (``""``) image."""
    return f"image: url({path});" if path else ""


def tick_image(background: str) -> str:
    """Path of a tick (Lucide ``check``) SVG readable on ``background``."""
    r, g, b = (int(background.lstrip("#")[i:i + 2], 16) for i in (0, 2, 4))
    tick = "#1a1a1a" if 0.299 * r + 0.587 * g + 0.114 * b > 150 else "#ffffff"
    return tinted_icon("check", tick, "3")


def app_qss(dark: bool) -> str:
    """Application stylesheet (ui/style.qss or ui/style_dark.qss). Both themes set
    explicit colours so the chrome never depends on the OS light/dark mode."""
    arrow = theme("dark" if dark else "light")["text"]
    return (_read_qss(dark)
            .replace("image: url(@TICK@);", image_rule(tick_image("#2f6fed")))
            .replace("image: url(@UP@);", image_rule(tinted_icon("chevron-up", arrow, "3")))
            .replace("image: url(@DOWN@);",
                     image_rule(tinted_icon("chevron-down", arrow, "3"))))
