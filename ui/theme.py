"""Light/dark palettes shared by all ab1kit widgets."""

from __future__ import annotations

import sys
from functools import lru_cache
from pathlib import Path

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
}


def theme(name: str) -> dict:
    return DARK if name == "dark" else LIGHT


def resource_dir() -> Path:
    """Folder holding style*.qss and icons/: next to this module, or under
    ``sys._MEIPASS/ui`` in a PyInstaller bundle (see ab1kit.spec datas)."""
    frozen = getattr(sys, "_MEIPASS", None)
    return Path(frozen) / "ui" if frozen else Path(__file__).resolve().parent


@lru_cache(maxsize=2)
def _read_qss(dark: bool) -> str:
    name = "style_dark.qss" if dark else "style.qss"
    return (resource_dir() / name).read_text(encoding="utf-8")


def app_qss(dark: bool) -> str:
    """Application stylesheet (ui/style.qss or ui/style_dark.qss). Both themes set
    explicit colours so the chrome never depends on the OS light/dark mode."""
    return _read_qss(dark)
