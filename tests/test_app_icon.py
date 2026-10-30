"""App icon: committed assets, regeneration script, window icon."""

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import struct
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
APP_ICONS = ROOT / "ui" / "icons" / "app"


def test_png_and_ico_exist():
    assert (APP_ICONS / "ab1kit.png").is_file()
    assert (APP_ICONS / "ab1kit.ico").is_file()


def test_ico_has_required_sizes():
    data = (APP_ICONS / "ab1kit.ico").read_bytes()
    reserved, kind, count = struct.unpack_from("<HHH", data, 0)
    assert (reserved, kind) == (0, 1)
    sizes = set()
    for i in range(count):
        w, h = struct.unpack_from("<BB", data, 6 + 16 * i)
        sizes.add(w or 256)
        assert (w or 256) == (h or 256)
    assert {16, 24, 32, 48, 64, 128, 256} <= sizes


def test_svg_has_no_text_elements():
    svg = (APP_ICONS / "ab1kit.svg").read_text(encoding="utf-8")
    assert "<text" not in svg


def test_make_icon_regenerates_identical_png(tmp_path):
    env = {k: v for k, v in os.environ.items() if k != "QT_QPA_FONTDIR"}
    env["QT_QPA_PLATFORM"] = "offscreen"
    subprocess.run([sys.executable, str(ROOT / "scripts" / "make_icon.py"),
                    "--out", str(tmp_path)], check=True, env=env, cwd=ROOT)
    assert (tmp_path / "ab1kit.png").read_bytes() == (APP_ICONS / "ab1kit.png").read_bytes()
    assert (tmp_path / "ab1kit.ico").read_bytes() == (APP_ICONS / "ab1kit.ico").read_bytes()


def test_main_window_sets_application_icon(monkeypatch):
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    monkeypatch.setattr("ui.main_window.save_config", lambda *a, **k: None)
    from ui.main_window import MainWindow

    app.setWindowIcon(__import__("PySide6.QtGui", fromlist=["QIcon"]).QIcon())
    window = MainWindow()
    icon = app.windowIcon()
    assert not icon.isNull()
    pix = icon.pixmap(32, 32)
    assert (pix.width(), pix.height()) == (32, 32)
    window.close()
