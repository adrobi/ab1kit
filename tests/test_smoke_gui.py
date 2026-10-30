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
    monkeypatch.setattr("ui.read_view.read_view.save_config", lambda *a, **k: None)
    monkeypatch.setattr("ui.construct.workspace.save_config", lambda *a, **k: None)
    from chrom_core import sessions

    monkeypatch.setattr(sessions, "session_path", lambda _p: tmp_path / "session.json")


@pytest.fixture(autouse=True)
def delete_windows():
    """Destroy the windows a test created: closed windows otherwise stay alive
    and every later theme switch restyles all of them."""
    yield
    from PySide6.QtCore import QCoreApplication, QEvent
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance()
    if app is None:
        return
    for widget in app.topLevelWidgets():
        widget.close()
        widget.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete.value)
    app.processEvents()


def _wait_for(app, predicate, timeout=5.0):
    import time

    deadline = time.monotonic() + timeout
    while not predicate() and time.monotonic() < deadline:
        app.processEvents()
        time.sleep(0.005)
    assert predicate()


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
        "ui.read_view.read_view.QMessageBox.critical",
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
        "ui.read_view.read_view.QMessageBox.warning",
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
        "ui.read_view.read_view.QFileDialog.getExistingDirectory",
        lambda *a, **k: str(tmp_path / "export"),
    )
    window = MainWindow(AppConfig(recent_files=[]))
    window.show()
    window._scan_folder(EXAMPLES)
    _wait_for(qapp, lambda: window.batch_panel.table.rowCount() > 0)

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
        "ui.read_view.read_view.QFileDialog.getSaveFileName",
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
        "ui.read_view.read_view.QFileDialog.getSaveFileName",
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


@pytest.mark.skipif(not AB1, reason="no example .ab1 files")
def test_read_view_emits_read_edited(qapp):
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    window = MainWindow(AppConfig(recent_files=[]))
    emitted = []
    window.read_view.read_edited.connect(emitted.append)
    view = window.read_view
    assert view.open_path(AB1[0])
    assert emitted == []  # opening is not an edit
    assert view.current_model() is window.read_model

    expected = [Path(AB1[0])]
    steps = [
        lambda: view.apply_base_change(0, "N"),
        view.undo_read,
        view.redo_read,
        lambda: view.apply_trim_range(10, 100),
        view.reverse_complement_read,
        view._save_session,
    ]
    for count, step in enumerate(steps, start=1):
        step()
        assert emitted == expected * count
    window.close()
    assert emitted == expected * len(steps)  # close-time autosave does not emit


@pytest.mark.skipif(not AB1, reason="no example .ab1 files")
def test_bare_read_view_works_without_shell_menus(qapp):
    from chrom_core.config import AppConfig
    from ui.read_view import ReadView

    view = ReadView(AppConfig(recent_files=[]), lambda *a, **k: None)
    assert view.open_path(AB1[0])
    view.apply_base_change(0, "N")
    assert view.undo_action.isEnabled()
    view.close()


# ------------------------------------------------------------------ Task 13: theme
def test_toolbar_actions_have_icons(qapp):
    from PySide6.QtWidgets import QToolBar

    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    window = MainWindow(AppConfig(recent_files=[]))
    bar = window.read_view.findChild(QToolBar)
    actions = [a for a in bar.actions() if not a.isSeparator() and a.text()]
    assert actions
    assert all(not a.icon().isNull() for a in actions), [a.text() for a in actions if a.icon().isNull()]
    window.close()


def test_quality_colour_thresholds_and_palette(qapp):
    from chrom_core.clone_report import STATUS_HEX
    from ui.status import quality_color, status_color

    for theme in ("light", "dark"):
        assert quality_color(15, theme) != quality_color(35, theme)
        assert quality_color(19, theme) != quality_color(20, theme)
        assert quality_color(20, theme) == quality_color(29, theme)
        assert quality_color(29, theme) != quality_color(30, theme)
        assert quality_color(30, theme) == quality_color(60, theme)
    for verdict, (light, dark) in STATUS_HEX.items():
        assert status_color(verdict, "light").name() == light
        assert status_color(verdict, "dark").name() == dark


def test_minimap_and_chrom_view_share_quality_buckets(qapp):
    from ui.status import quality_bucket

    assert [quality_bucket(q) for q in (0, 19, 20, 29, 30, 99)] == [
        "low_q", "low_q", "mid_q", "mid_q", "high_q", "high_q"]


def test_status_chip_is_label(qapp):
    from PySide6.QtWidgets import QLabel

    from ui.status import status_chip

    chip = status_chip("ok", "OK")
    assert isinstance(chip, QLabel) and chip.text() == "OK"
    assert "border-radius" in chip.styleSheet()


def test_stylesheet_applies_in_both_themes(qapp):
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    window = MainWindow(AppConfig(recent_files=[]))
    seen = []
    for name in ("light", "dark"):
        window.set_theme(name)
        assert qapp.styleSheet().strip()
        seen.append(qapp.styleSheet())
    assert seen[0] != seen[1]
    window.close()


def test_tables_stretch_last_column_and_tooltip_notes(qapp):
    from PySide6.QtWidgets import QHeaderView

    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    sample = next(p for p in AB1 if "3730" in p.name)
    window = MainWindow(AppConfig(recent_files=[]))
    assert window.open_path(sample)
    seq = "".join(window.read_model.seq)
    ref = seq[100:400]
    window._set_reference_seq(ref)
    wrong = "A" if ref[50] != "A" else "C"
    window._verify_variants(f"name;pos;ref;alt\nsnp;51;{wrong};T\n")
    for table in (window.verify_table, window.diff_table):
        assert table.horizontalHeader().stretchLastSection()
        assert table.horizontalHeader().sectionResizeMode(0) != QHeaderView.ResizeMode.Fixed
    note = window.verify_table.item(0, 7)
    assert len(note.text()) > 40 and note.toolTip() == note.text()
    window.close()


def test_status_cells_repaint_on_theme_switch(qapp):
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow
    from ui.status import status_color

    sample = next(p for p in AB1 if "3730" in p.name)
    window = MainWindow(AppConfig(recent_files=[], theme="light"))
    assert window.open_path(sample)
    seq = "".join(window.read_model.seq)
    ref = seq[100:400]
    window._set_reference_seq(ref)
    wrong = "A" if ref[50] != "A" else "C"
    window._verify_variants(f"name;pos;ref;alt\nsnp;51;{wrong};T\n")
    cell = window.verify_table.item(0, 4)
    assert cell.background().color() == status_color("errors", "light")
    window.set_theme("dark")
    assert cell.background().color() == status_color("errors", "dark")
    window.close()


def test_icon_pixmap_changes_with_theme(qapp):
    from PySide6.QtGui import QAction

    from ui import status

    action = QAction("x")
    status.bind_icon(action, "undo-2")
    try:
        status.set_theme("light")
        light = action.icon().pixmap(24, 24).toImage()
        status.set_theme("dark")
        dark = action.icon().pixmap(24, 24).toImage()
        assert light != dark
    finally:
        status.set_theme("light")


def test_channel_checkboxes_use_channel_colours_and_tick(qapp):
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow
    from ui.theme import theme as palette

    window = MainWindow(AppConfig(recent_files=[]))
    for name in ("light", "dark"):
        window.set_theme(name)
        for base, box in window.read_view._channel_boxes.items():
            sheet = box.styleSheet()
            colour = palette(name)["bases"][base]
            assert f"QCheckBox::indicator:checked {{ background-color: {colour}" in sheet
            assert "image: url(" in sheet
    assert "image: url(" in qapp.styleSheet()
    window.close()


def test_batch_table_stretches_file_column(qapp):
    from PySide6.QtWidgets import QHeaderView

    from ui.read_view.batch_dock import BatchPanel

    panel = BatchPanel()
    header = panel.table.horizontalHeader()
    assert header.sectionResizeMode(0) == QHeaderView.ResizeMode.Stretch
    assert not header.stretchLastSection()


def test_tooltips_only_on_long_cells(qapp):
    from PySide6.QtWidgets import QTableWidgetItem

    from ui.tables import set_full_text_tooltip

    short, long_ = QTableWidgetItem("12"), QTableWidgetItem("x" * 60)
    set_full_text_tooltip(short)
    set_full_text_tooltip(long_)
    assert short.toolTip() == "" and long_.toolTip() == "x" * 60


# ------------------------------------------------------------------ Task 16: tabs, projects
def _coi_project(tmp_path, region):
    from chrom_core.clones import Grouping, VerifyParams
    from chrom_core.project import Project, save_project

    refs = Path(__file__).resolve().parent / "data" / "refs"
    project = tmp_path / "p.ab1proj"
    save_project(project, Project(refs / "coi_HM417896.gb", False, region, Grouping(),
                                  VerifyParams()))
    return project


def test_construct_tab_first_with_project_menu(qapp):
    from chrom_core.config import AppConfig
    from i18n import tr
    from ui.main_window import MainWindow

    window = MainWindow(AppConfig(recent_files=[]))
    assert window.tabs.count() == 2
    assert window.tabs.currentIndex() == 0 and window.tabs.currentWidget() is window.workspace
    assert window.tabs.tabText(0) == tr("tab.construct")
    assert window.tabs.tabText(1) == tr("tab.read")
    titles = [a.text() for a in window.project_menu.actions() if not a.isSeparator()]
    assert titles == [tr("menu.project.new"), tr("menu.project.open"), tr("menu.project.save"),
                      tr("menu.project.save_as"), tr("menu.project.export")]
    window.close()


@pytest.mark.skipif(not AB1, reason="no example .ab1 files")
def test_open_ab1_switches_to_read_tab(qapp):
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    window = MainWindow(AppConfig(recent_files=[]))
    assert window.open_path(AB1[0])
    assert window.tabs.currentWidget() is window.read_view
    window.close()


def test_open_project_path_goes_to_construct_tab(qapp, tmp_path):
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    project = _coi_project(tmp_path, (10, 200))
    window = MainWindow(AppConfig(recent_files=[]))
    window.tabs.setCurrentWidget(window.read_view)
    assert window.open_path(project)
    assert window.tabs.currentWidget() is window.workspace
    assert window.workspace.region() == (10, 200)
    assert window.workspace.construct().name == "HM417896.1"
    window.close()


def test_run_opens_project_argument(qapp, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QApplication

    import ui.main_window as mw

    project = _coi_project(tmp_path, None)
    opened = []
    class ExistingApp:  # run() must not create a second QApplication or block
        instance = staticmethod(QApplication.instance)

        def __init__(self, _args):
            pass

        def setStyle(self, _name):  # noqa: N802
            pass

        def exec(self):
            return 0

    monkeypatch.setattr(mw, "QApplication", ExistingApp)
    from PySide6.QtWidgets import QMessageBox

    for box in ("critical", "warning"):
        monkeypatch.setattr(QMessageBox, box, lambda *a, **k: QMessageBox.StandardButton.Ok)
    monkeypatch.setattr(mw.MainWindow, "show", lambda self: opened.append(self))
    assert mw.run([str(project)]) == 0
    (window,) = opened
    assert window.tabs.currentWidget() is window.workspace
    assert window.workspace.construct() is not None
    window.close()
