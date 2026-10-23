"""GUI smoke tests: run Qt offscreen, no real display needed."""

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from pathlib import Path

import pytest

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"
AB1 = sorted(EXAMPLES.glob("*.ab1"))


@pytest.fixture(scope="module")
def qapp():
    from PySide6.QtWidgets import QApplication

    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def no_config_writes(monkeypatch):
    monkeypatch.setattr("ui.main_window.save_config", lambda *a, **k: None)


def test_entry_point_importable():
    from chrom_cli import main
    from ui.main_window import MainWindow, run

    assert callable(main)
    assert callable(run)
    assert MainWindow is not None


@pytest.mark.skipif(not AB1, reason="no example .ab1 files")
def test_main_window_opens_file_and_renders(qapp):
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    window = MainWindow(AppConfig(recent_files=[]))
    assert window.open_path(AB1[0])
    assert window.chrom_view.chrom is not None
    assert window.minimap.chrom is not None

    pixmap = window.grab()
    assert not pixmap.isNull()
    assert pixmap.width() > 100 and pixmap.height() > 100

    window.set_theme("dark")
    window.minimap.set_profile(True)
    window.chrom_view.set_highlights([(10, 40)])
    window.chrom_view.goto_base(50)
    assert not window.grab().isNull()

    window.close()


@pytest.mark.skipif(not AB1, reason="no example .ab1 files")
def test_open_broken_file_shows_error(qapp, monkeypatch):
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    shown = {}
    monkeypatch.setattr(
        "ui.main_window.QMessageBox.critical",
        lambda *a, **k: shown.setdefault("called", True),
    )
    window = MainWindow(AppConfig(recent_files=[]))
    assert window.open_path("no-such-file.ab1") is False
    assert shown.get("called")
    window.close()
