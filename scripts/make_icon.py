"""Regenerate ui/icons/app/ab1kit.png and ab1kit.ico from ab1kit.svg.

Usage: python scripts/make_icon.py [--out DIR]

The SVG is the source of truth (letters are vector paths, so no fonts are
needed). The ICO holds PNG-compressed entries, assembled here without extra
dependencies.

Note: the SVG letter outlines were extracted once from Consolas Bold 22px and are
baked in as paths; this script does not touch them.
"""

from __future__ import annotations

import argparse
import os
import struct
import sys
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QBuffer, QByteArray, QIODevice
from PySide6.QtGui import QGuiApplication, QImage, QPainter
from PySide6.QtSvg import QSvgRenderer

APP_ICONS = Path(__file__).resolve().parent.parent / "ui" / "icons" / "app"
ICO_SIZES = (16, 24, 32, 48, 64, 128, 256)


def render(renderer: QSvgRenderer, size: int) -> QImage:
    image = QImage(size, size, QImage.Format.Format_ARGB32)
    image.fill(0)
    painter = QPainter(image)
    painter.setRenderHints(QPainter.RenderHint.Antialiasing
                           | QPainter.RenderHint.SmoothPixmapTransform)
    renderer.render(painter)
    painter.end()
    return image


def png_bytes(image: QImage) -> bytes:
    data = QByteArray()
    buf = QBuffer(data)
    buf.open(QIODevice.OpenModeFlag.WriteOnly)
    image.save(buf, "PNG")
    buf.close()
    return bytes(data)


def build_ico(entries: list[tuple[int, bytes]]) -> bytes:
    header = struct.pack("<HHH", 0, 1, len(entries))
    offset = 6 + 16 * len(entries)
    directory, blobs = b"", b""
    for size, blob in entries:
        dim = 0 if size >= 256 else size
        directory += struct.pack("<BBBBHHII", dim, dim, 0, 0, 1, 32, len(blob), offset)
        blobs += blob
        offset += len(blob)
    return header + directory + blobs


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=APP_ICONS, help="output folder")
    args = parser.parse_args(argv)

    QGuiApplication.instance() or QGuiApplication([])
    renderer = QSvgRenderer(str(APP_ICONS / "ab1kit.svg"))
    if not renderer.isValid():
        print("invalid SVG", file=sys.stderr)
        return 1
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "ab1kit.png").write_bytes(png_bytes(render(renderer, 256)))
    entries = [(s, png_bytes(render(renderer, s))) for s in ICO_SIZES]
    (args.out / "ab1kit.ico").write_bytes(build_ico(entries))
    return 0


if __name__ == "__main__":
    sys.exit(main())
