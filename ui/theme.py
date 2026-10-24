"""Light/dark palettes shared by all ab1kit widgets."""

from __future__ import annotations

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


def app_qss(dark: bool) -> str:
    """Application stylesheet. Both themes set explicit colors so the chrome
    never depends on the OS light/dark mode (native dark Windows makes our
    fixed light-theme text colors unreadable)."""
    if dark:
        return (
            "QWidget { background-color: #252526; color: #d4d4d4; }"
            "QMenuBar, QMenu, QToolBar { background-color: #2d2d30; }"
            "QLineEdit, QSpinBox { background-color: #1e1e1e; }"
            "QDockWidget::title { background-color: #2d2d30; }"
        )
    return (
        "QWidget { background-color: #f3f3f3; color: #1a1a1a; }"
        "QMenuBar, QMenu, QToolBar { background-color: #e9e9e9; }"
        "QLineEdit, QSpinBox { background-color: #ffffff; }"
        "QDockWidget::title { background-color: #e9e9e9; }"
    )
