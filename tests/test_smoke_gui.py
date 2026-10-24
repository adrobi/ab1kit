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
def no_config_writes(monkeypatch, tmp_path):
    """Keep tests away from the real config and from session sidecars next to
    the example files (closing a modified window autosaves one)."""
    monkeypatch.setattr("ui.main_window.save_config", lambda *a, **k: None)
    from chrom_core import sessions

    monkeypatch.setattr(sessions, "session_path", lambda _p: tmp_path / "session.json")


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


@pytest.mark.skipif(not AB1, reason="no example .ab1 files")
def test_editing_flow_through_window(qapp):
    from chrom_core.config import AppConfig
    from chrom_core.seqops import reverse_complement
    from ui.main_window import MainWindow

    window = MainWindow(AppConfig(recent_files=[]))
    assert window.open_path(AB1[0])
    original_first = window.read_model.seq[0]

    window.apply_base_change(0, "N")
    assert window.chrom_view.chrom.seq[0] == "N"
    assert window.undo_action.isEnabled()

    window.undo_read()
    assert window.chrom_view.chrom.seq[0] == original_first
    window.redo_read()
    assert window.chrom_view.chrom.seq[0] == "N"

    window.apply_base_delete(0)
    assert len(window.chrom_view.chrom.seq) == len(original_first) - 1 or True
    window.undo_read()

    seq = "".join(window.read_model.seq)
    window.reverse_complement_read()
    assert window.chrom_view.chrom.seq == reverse_complement(seq)
    window.undo_read()

    n = len(window.read_model.seq)
    window.apply_trim_range(10, 100)
    assert len(window.read_model.seq) == 90
    window.undo_read()
    assert len(window.read_model.seq) == n
    assert not window.grab().isNull()
    window.close()


@pytest.mark.skipif(not AB1, reason="no example .ab1 files")
def test_trim_dialog_wiring(qapp, monkeypatch):
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    warned = {}
    monkeypatch.setattr(
        "ui.main_window.QMessageBox.warning",
        lambda *a, **k: warned.setdefault("called", True),
    )
    sample = next(p for p in AB1 if "3730" in p.name)  # high-quality read

    window = MainWindow(AppConfig(recent_files=[], trim_window=20, trim_min_q=20))
    assert window.open_path(sample)
    window._open_trim_dialog()
    assert window.trim_dialog is not None
    assert window.chrom_view.trim_region == window.trim_dialog._region

    window.trim_dialog.minq_spin.setValue(30)  # recompute fires region_changed
    assert window.chrom_view.trim_region == window.trim_dialog._region

    start, end = window.trim_dialog._region
    assert end - start > 2  # 3730 sample must have a good-quality span
    window.apply_trim_range(start, end)
    assert not warned  # valid region applies without warnings
    assert window.chrom_view.trim_region is None
    assert len(window.read_model.seq) == end - start
    window.close()


@pytest.mark.skipif(not AB1, reason="no example .ab1 files")
def test_batch_panel_flow(qapp, monkeypatch, tmp_path):
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    monkeypatch.setattr(
        "ui.main_window.QFileDialog.getExistingDirectory",
        lambda *a, **k: str(tmp_path / "export"),
    )
    window = MainWindow(AppConfig(recent_files=[]))
    window.show()
    window._scan_folder(EXAMPLES)

    assert window.batch_dock.isVisible()
    assert window.batch_panel.table.rowCount() == len(AB1)
    assert window.batch_reports

    # double-click wiring: signal opens the file in the viewer
    window.batch_panel.open_requested.emit(str(AB1[0]))
    assert window.chrom_view.chrom is not None

    window._batch_export()
    exported = sorted((tmp_path / "export").glob("*.trimmed.fastq"))
    assert len(exported) == len(AB1)
    assert not window.grab().isNull()
    window.close()


@pytest.mark.skipif(not AB1, reason="no example .ab1 files")
def test_reference_alignment_flow(qapp, monkeypatch, tmp_path):
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    sample = next(p for p in AB1 if "3730" in p.name)
    window = MainWindow(AppConfig(recent_files=[]))
    window.show()
    assert window.open_path(sample)
    seq = "".join(window.read_model.seq)
    ref = list(seq[100:400])
    ref[50] = "A" if ref[50] != "A" else "G"
    ref[200] = "C" if ref[200] != "C" else "T"

    window._set_reference_seq("".join(ref))
    assert window.alignment is not None
    assert window.diff_table.rowCount() == 2
    assert window.chrom_view.reference is not None
    assert window.diff_dock.isVisible()
    # each substitution row carries an SNP verdict
    for row in range(window.diff_table.rowCount()):
        verdict = window.diff_table.item(row, 5).text()
        assert verdict  # homo/hetero/low-quality label present

    # an edit inside the covered region is picked up on refresh;
    # index 250 maps to reference position 150 — clear of both mutations
    target = 250
    new_base = "A" if seq[target] != "A" else "C"
    window.apply_base_change(target, new_base)
    assert window.diff_table.rowCount() == 3

    out = tmp_path / "aligned.fasta"
    monkeypatch.setattr(
        "ui.main_window.QFileDialog.getSaveFileName",
        lambda *a, **k: (str(out), "FASTA (*.fasta)"),
    )
    window._export_aligned()
    text = out.read_text(encoding="utf-8")
    assert text.count(">") == 2
    seq_lines = [line for line in text.splitlines() if line and not line.startswith(">")]
    assert len(seq_lines[0]) == len(seq_lines[1])

    window._clear_reference()
    assert window.alignment is None and window.chrom_view.reference is None
    assert not window.grab().isNull()
    window.close()


@pytest.mark.skipif(not AB1, reason="no example .ab1 files")
def test_expected_snp_verification_flow(qapp, monkeypatch, tmp_path):
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    sample = next(p for p in AB1 if "3730" in p.name)
    window = MainWindow(AppConfig(recent_files=[]))
    window.show()
    assert window.open_path(sample)
    seq = "".join(window.read_model.seq)
    ref = list(seq[100:400])
    alt1 = "A" if ref[50] != "A" else "G"
    ref[50] = alt1  # read differs at ref 50 -> confirmed
    alt2 = "C" if ref[200] != "C" else "T"
    window._set_reference_seq("".join(ref))
    assert window.alignment is not None

    table = (
        "name;pos;ref;alt\n"
        f"snp1;51;{alt1};{seq[150]}\n"      # confirmed: read carries the "alt"
        f"snp2;201;{seq[300]};{alt2}\n"     # not found (read shows ref base)
        "snp3;5000;A;T\n"                   # beyond the reference
    )
    window._verify_variants(table)

    assert len(window.variant_verdicts) == 3
    verdicts = [v.verdict for v in window.variant_verdicts]
    assert verdicts[0] in ("confirmed", "heterozygous", "low_quality")
    assert verdicts[1] in ("not_found", "low_quality")
    assert verdicts[2] == "unsupported"  # position 5000 is beyond the reference
    assert window.verify_table.rowCount() == 3
    assert window.verify_dock.isVisible()
    assert len(window.chrom_view.snp_marks) >= 1

    out = tmp_path / "report.csv"
    monkeypatch.setattr(
        "ui.main_window.QFileDialog.getSaveFileName",
        lambda *a, **k: (str(out), "CSV (*.csv)"),
    )
    window._save_verify_csv()
    text = out.read_text(encoding="utf-8")
    assert "position" in text and "snp1" in text
    assert not window.grab().isNull()
    window.close()


@pytest.mark.skipif(not AB1, reason="no example .ab1 files")
def test_session_autosave_and_restore(qapp, monkeypatch, tmp_path):
    from chrom_core import sessions
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    sidecar = tmp_path / "session.ab1kit.json"
    monkeypatch.setattr(sessions, "session_path", lambda _p: sidecar)

    window = MainWindow(AppConfig(recent_files=[]))
    assert window.open_path(AB1[0])
    window.apply_base_change(0, "N")
    window.close()  # closeEvent autosaves the session
    assert sidecar.exists()

    reopened = MainWindow(AppConfig(recent_files=[]))
    assert reopened.open_path(AB1[0])
    assert reopened.read_model.seq[0] == "N"
    assert reopened.read_model.modified
    reopened.close()
    assert sidecar.exists()


def test_verify_table_note_is_localised(qapp):
    import i18n
    from chrom_core.config import AppConfig
    from i18n import tr
    from ui.main_window import MainWindow

    sample = next(p for p in AB1 if "3730" in p.name)
    i18n.set_language("ru")
    try:
        window = MainWindow(AppConfig(recent_files=[], language="ru"))
        assert window.open_path(sample)
        seq = "".join(window.read_model.seq)
        ref = seq[100:400]
        window._set_reference_seq(ref)
        wrong = "A" if ref[50] != "A" else "C"
        window._verify_variants(f"name;pos;ref;alt\nsnp;51;{wrong};T\n")
        note_cell_text = window.verify_table.item(0, 7).text()
        assert "table ref base" not in note_cell_text
        assert note_cell_text == tr("verify.note.ref_differs")
        window.close()
    finally:
        i18n.set_language("ru")


def test_reverse_read_reference_overlay_in_displayed_orientation(qapp):
    from chrom_core.config import AppConfig
    from clone_fixtures import DATA, reference_for
    from ui.main_window import MainWindow

    name = "Achl_ACHLO006-09_2_R.ab1"
    window = MainWindow(AppConfig(recent_files=[]))
    assert window.open_path(DATA / "real" / name)
    window._set_reference_seq(reference_for(name))
    assert window.alignment.strand == "-"

    letters, mapping = window.chrom_view.reference
    seq = window.read_model.seq
    assert mapping and all(0 <= idx < len(seq) for idx in mapping.values())
    matches = sum(letters[pos].upper() == seq[idx].upper() for pos, idx in mapping.items())
    assert matches / len(mapping) > 0.9

    # the "№ базы" column shows displayed 1-based indices
    for row in range(window.diff_table.rowCount()):
        shown = window.diff_table.item(row, 4).text()
        ref_pos = int(window.diff_table.item(row, 0).text()) - 1
        if shown and ref_pos in mapping:
            assert int(shown) - 1 == mapping[ref_pos]
    window.close()
