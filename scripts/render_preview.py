"""Render the chromatogram widget offscreen to a PNG (dev smoke check).

Usage: python scripts/render_preview.py [sample.ab1] [out.png] [--theme light|dark] [--window]
       python scripts/render_preview.py --construct [out.png] [--theme light|dark] [--lang ru|en] [--empty]

--window renders the whole main window (toolbar, docks, tables) instead of the bare widget.
--construct renders the "construct check" tab: the COI reference of tests/data with one
substitution introduced (so one clone shows an error) against the 8 real reads, with the
detail pane open on that difference; --empty renders the empty state instead.
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


def render_window(src: Path, out: Path, theme: str, uncheck: str = "") -> None:
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
    for base in uncheck:
        window.read_view._channel_boxes[base].setChecked(False)
    window.batch_dock.show()
    window.batch_panel.show_reports(batch.scan_folder(src.parent))
    for _ in range(20):
        app.processEvents()
    window.grab().save(str(out))
    print(f"{src.name} [{theme}] -> {out}")
    window.close()


ROOT = Path(__file__).resolve().parent.parent
MUTATION_AT = 250  # 0-based reference position of the introduced substitution


def _mutated_reference(folder: Path) -> Path:
    from Bio import SeqIO
    from Bio.Seq import Seq

    record = SeqIO.read(ROOT / "tests" / "data" / "refs" / "coi_HM417896.gb", "genbank")
    seq = str(record.seq).upper()
    base = "G" if seq[MUTATION_AT] != "G" else "T"
    record.seq = Seq(seq[:MUTATION_AT] + base + seq[MUTATION_AT + 1:])
    target = folder / "COI_construct.gb"
    SeqIO.write(record, target, "genbank")
    return target


def render_construct(out: Path, theme: str, lang: str, empty: bool) -> None:
    import time

    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow
    import ui.main_window as main_window
    import ui.read_view.read_view as read_view

    main_window.save_config = read_view.save_config = lambda *a, **k: None  # keep the user's config
    app = QApplication.instance()
    app.setStyle("Fusion")
    window = MainWindow(AppConfig(recent_files=[], theme=theme, language=lang))
    window.resize(1300, 820)
    window.show()
    ws = window.workspace
    if not empty:
        folder = Path(tempfile.mkdtemp(prefix="ab1kit_construct_"))
        ws.handle_drop([_mutated_reference(folder), ROOT / "tests" / "data" / "real"])
        deadline = time.monotonic() + 60
        while ws.is_verifying() and time.monotonic() < deadline:
            app.processEvents()
            time.sleep(0.01)
        ws.set_region((6, 647))
        while ws.is_verifying() and time.monotonic() < deadline:
            app.processEvents()
            time.sleep(0.01)
        clone = next(r.clone for r in ws.reports() if r.verdict == "errors")
        ws.matrix.difference_clicked.emit(clone, MUTATION_AT)
    for _ in range(20):
        app.processEvents()
    window.grab().save(str(out))
    print(f"construct [{theme}, {lang}{', empty' if empty else ''}] -> {out}")
    window.close()


def main() -> None:
    examples = sorted((Path(__file__).resolve().parent.parent / "examples").glob("*.ab1"))
    argv = sys.argv[1:]
    theme = "light"
    if "--theme" in argv:
        i = argv.index("--theme")
        theme = argv[i + 1]
        del argv[i:i + 2]
    uncheck = ""
    if "--uncheck" in argv:
        i = argv.index("--uncheck")
        uncheck = argv[i + 1]
        del argv[i:i + 2]
    lang = "ru"
    if "--lang" in argv:
        i = argv.index("--lang")
        lang = argv[i + 1]
        del argv[i:i + 2]
    if "--construct" in argv:
        rest = [a for a in argv if a not in ("--construct", "--empty")]
        out = Path(rest[0]) if rest else Path(tempfile.gettempdir()) / "ab1kit_construct.png"
        app = QApplication([])  # noqa: F841 - must outlive the window
        render_construct(out, theme, lang, "--empty" in argv)
        return
    window_mode = "--window" in argv
    argv = [a for a in argv if a != "--window"]
    src = Path(argv[0]) if argv else examples[0]
    out = Path(argv[1]) if len(argv) > 1 else Path(tempfile.gettempdir()) / "ab1kit_preview.png"

    app = QApplication([])
    if window_mode:
        render_window(src, out, theme, uncheck)
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
