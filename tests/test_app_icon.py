"""App icon: committed assets, regeneration script, window icon."""

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import struct
import subprocess
import sys
from pathlib import Path

import pytest
from PySide6.QtGui import QIcon, QImage

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


def _ico_entries(data: bytes) -> list[tuple[int, bytes]]:
    count = struct.unpack_from("<H", data, 4)[0]
    out = []
    for i in range(count):
        w, _h, _c, _r, _p, _b, size, offset = struct.unpack_from("<BBBBHHII", data, 6 + 16 * i)
        out.append((w or 256, data[offset:offset + size]))
    return out


def _decode(blob: bytes) -> QImage:
    image = QImage.fromData(blob, "PNG")
    assert not image.isNull()
    return image.convertToFormat(QImage.Format.Format_ARGB32)


def _assert_pixels_close(a: QImage, b: QImage) -> None:
    """Max per-channel diff <= 2, except at most 0.1% outlier pixels."""
    assert (a.width(), a.height()) == (b.width(), b.height())
    outliers = 0
    for y in range(a.height()):
        for x in range(a.width()):
            pa, pb = a.pixel(x, y), b.pixel(x, y)
            if any(abs(((pa >> s) & 255) - ((pb >> s) & 255)) > 2 for s in (0, 8, 16, 24)):
                outliers += 1
    assert outliers <= max(1, a.width() * a.height() // 1000)


def test_make_icon_regenerates_same_pixels(tmp_path):
    env = {k: v for k, v in os.environ.items() if k != "QT_QPA_FONTDIR"}
    env["QT_QPA_PLATFORM"] = "offscreen"
    subprocess.run([sys.executable, str(ROOT / "scripts" / "make_icon.py"),
                    "--out", str(tmp_path)], check=True, env=env, cwd=ROOT)
    _assert_pixels_close(_decode((tmp_path / "ab1kit.png").read_bytes()),
                         _decode((APP_ICONS / "ab1kit.png").read_bytes()))
    new = _ico_entries((tmp_path / "ab1kit.ico").read_bytes())
    old = _ico_entries((APP_ICONS / "ab1kit.ico").read_bytes())
    assert [size for size, _ in new] == [size for size, _ in old]
    for (_, blob_new), (_, blob_old) in zip(new, old):
        _assert_pixels_close(_decode(blob_new), _decode(blob_old))


def test_main_window_sets_application_icon(monkeypatch):
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    monkeypatch.setattr("ui.main_window.save_config", lambda *a, **k: None)
    from ui.main_window import MainWindow

    app.setWindowIcon(QIcon())
    window = MainWindow()
    icon = app.windowIcon()
    assert not icon.isNull()
    pix = icon.pixmap(32, 32)
    assert (pix.width(), pix.height()) == (32, 32)
    assert pix.toImage().convertToFormat(QImage.Format.Format_ARGB32).pixelColor(16, 16).alpha() > 0
    window.close()
