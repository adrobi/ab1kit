"""Render the chromatogram widget offscreen to a PNG (dev smoke check).

Usage: python scripts/render_preview.py [sample.ab1] [out.png]
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication  # noqa: E402

from chrom_core import read_ab1  # noqa: E402
from ui.main_window import ChromatogramWidget  # noqa: E402


def main() -> None:
    examples = sorted((Path(__file__).resolve().parent.parent / "examples").glob("*.ab1"))
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else examples[0]
    out = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(tempfile.gettempdir()) / "ab1kit_preview.png"

    app = QApplication([])
    view = ChromatogramWidget()
    view.resize(1200, 700)
    view.set_chromatogram(read_ab1(src))
    half = view._total_points() / 2
    view.offset = max(0.0, min(half, half - view._visible_points() / 2))
    view.grab().save(str(out))
    print(f"{src.name} -> {out}")


if __name__ == "__main__":
    main()
