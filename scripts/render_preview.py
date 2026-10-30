"""Render the chromatogram widget offscreen to a PNG (dev smoke check).

Usage: python scripts/render_preview.py [sample.ab1] [out.png] [--theme light|dark] [--window]

--window renders the whole main window (toolbar, docks, tables) instead of the bare widget.
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication  # noqa: E402

from chrom_core import batch, read_ab1  # noqa: E402
from ui.chrom_view import ChromatogramWidget  # noqa: E402


def render_window(src: Path, out: Path, theme: str) -> None:
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    app = QApplication.instance()
    app.setStyle("Fusion")
    window = MainWindow(AppConfig(recent_files=[], theme=theme))
    window.resize(1300, 820)
    window.show()
    window.open_path(src)
    seq = "".join(window.read_model.seq)
    ref = list(seq[100:400])
    ref[50] = "A" if ref[50] != "A" else "C"
    ref[120] = "G" if ref[120] != "G" else "T"
    window._set_reference_seq("".join(ref))
    rows = ["name;pos;ref;alt", f"snp1;51;{seq[150]};T", f"snp2;121;{seq[220]};G"]
    window._verify_variants("\n".join(rows) + "\n")
    window.batch_dock.show()
    window.batch_panel.show_reports(batch.scan_folder(src.parent))
    for _ in range(20):
        app.processEvents()
    window.grab().save(str(out))
    print(f"{src.name} [{theme}] -> {out}")
    window.close()


def main() -> None:
    examples = sorted((Path(__file__).resolve().parent.parent / "examples").glob("*.ab1"))
    argv = sys.argv[1:]
    theme = "light"
    if "--theme" in argv:
        i = argv.index("--theme")
        theme = argv[i + 1]
        del argv[i:i + 2]
    window_mode = "--window" in argv
    argv = [a for a in argv if a != "--window"]
    src = Path(argv[0]) if argv else examples[0]
    out = Path(argv[1]) if len(argv) > 1 else Path(tempfile.gettempdir()) / "ab1kit_preview.png"

    app = QApplication([])
    if window_mode:
        render_window(src, out, theme)
        return
    view = ChromatogramWidget()
    view.resize(1200, 700)
    view.set_chromatogram(read_ab1(src))
    half = view._total_points() / 2
    view.offset = max(0.0, min(half, half - view._visible_points() / 2))
    view.grab().save(str(out))
    print(f"{src.name} -> {out}")


if __name__ == "__main__":
    main()
