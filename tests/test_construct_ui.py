"""Construct screen widgets: summary bar and clone matrix (offscreen Qt).

Reports come from ``verify_all`` on the real pair ``Achl_ACHLO006-09`` F+R
against its COI reference, with errors / trimming / double peaks introduced
through editing sessions (``write_read``), as in ``test_clone_verify``.
"""

from __future__ import annotations

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PySide6.QtCore import QEvent, QPoint, QPointF, Qt
from PySide6.QtGui import QColor, QImage, QMouseEvent, QWheelEvent
from PySide6.QtWidgets import QApplication

from chrom_core.align import align_to_construct
from chrom_core.clones import CloneDifference, CloneReport, Grouping, ReadRef, verify_all
from chrom_core.editing import EditableRead
from chrom_core.reference import load_construct
from clone_fixtures import COI_GB, DATA, real_read, with_secondary_peak, write_read

F_NAME = "Achl_ACHLO006-09_1_F.ab1"
R_NAME = "Achl_ACHLO006-09_2_R.ab1"
F = DATA / "real" / F_NAME
R = DATA / "real" / R_NAME
REGION = (6, 647)
P = 250          # substitution site inside the region
MIXED_AT = (200, 300, 400)


@pytest.fixture(scope="module")
def qapp():
    return QApplication.instance() or QApplication([])


@pytest.fixture(scope="module")
def construct():
    return load_construct(COI_GB)


def _f_index(construct, ref_pos: int) -> int:
    al = align_to_construct(real_read(F_NAME).seq, construct)
    assert al.strand == "+"
    return al.ref_to_read[ref_pos]


@pytest.fixture(scope="module")
def reports(construct, tmp_path_factory):
    """ok (c1), errors confirmed (c2, F only), errors disputed (c3), incomplete (c4), mixed (c5)."""
    folder = tmp_path_factory.mktemp("clones")
    i = _f_index(construct, P)
    # a C peak 4x the called base, then the call edited to C: the trace agrees
    # with the edit (no artificial double peak at P)
    assert real_read(F_NAME).seq[i] != "C"
    edited = EditableRead(with_secondary_peak(real_read(F_NAME), i, "C", 4.0))
    edited.set_base(i, "C")
    sub_dir = folder / "sub"
    sub_dir.mkdir()
    f_sub = write_read(edited.current(), sub_dir, F_NAME)

    trimmed = EditableRead(real_read(F_NAME))
    trimmed.apply_trim(100, 400)
    trim_dir = folder / "trim"
    trim_dir.mkdir()
    f_trim = write_read(trimmed.current(), trim_dir, F_NAME)

    chrom = real_read(F_NAME)
    for pos in MIXED_AT:
        j = _f_index(construct, pos)
        chrom = with_secondary_peak(chrom, j, next(b for b in "ACGT" if b != chrom.seq[j]), 0.5)
    mix_dir = folder / "mix"
    mix_dir.mkdir()
    f_mix = write_read(chrom, mix_dir, F_NAME)

    grouping = Grouping([
        ReadRef(F, "c1", "F"), ReadRef(R, "c1", "R"),
        ReadRef(f_sub, "c2", "F"),
        ReadRef(f_sub, "c3", "F"), ReadRef(R, "c3", "R"),
        ReadRef(f_trim, "c4", "F"),
        ReadRef(f_mix, "c5", "F"), ReadRef(R, "c5", "R"),
    ])
    out = verify_all(construct, REGION, grouping)
    assert [r.verdict for r in out] == ["ok", "errors", "errors", "incomplete", "mixed"]
    assert [d.status for d in out[1].region_differences] == ["confirmed"]
    assert [d.status for d in out[2].region_differences] == ["disputed"]
    return out


def _matrix(qapp, construct, reports, width=900, theme="light"):
    from ui.construct.matrix import CloneMatrix

    m = CloneMatrix()
    m.set_theme(theme)
    m.set_data(construct, REGION, reports)
    m.resize(width, m.sizeHint().height())
    return m


def _send_mouse(widget, kind, pos: QPoint, button, buttons, modifiers=Qt.KeyboardModifier.NoModifier):
    local = QPointF(pos)
    event = QMouseEvent(kind, local, widget.mapToGlobal(local), button, buttons, modifiers)
    QApplication.sendEvent(widget, event)


# ------------------------------------------------------------------ summary bar


def test_summary_shows_construct_and_counts(qapp, construct, reports):
    from i18n import tr
    from ui.construct.summary_bar import SummaryBar

    bar = SummaryBar()
    bar.set_data(construct, REGION, reports)
    assert construct.name in bar.title_text()
    info = bar.info_text()
    assert "658" in info and "7–647" in info
    assert tr("construct.linear") in info
    counts = {v: bar.counter(v) for v in ("ok", "errors", "incomplete", "mixed")}
    assert "1" in counts["ok"].text() and "2" in counts["errors"].text()
    assert "1" in counts["incomplete"].text() and "1" in counts["mixed"].text()
    assert bar.counter("no_data").isHidden()   # zero counters are not shown


def test_summary_counter_toggles_filter(qapp, construct, reports):
    from ui.construct.summary_bar import SummaryBar

    bar = SummaryBar()
    bar.set_data(construct, REGION, reports)
    seen: list[set] = []
    bar.filter_changed.connect(seen.append)
    bar.counter("ok").click()
    assert seen[-1] == {"errors", "incomplete", "mixed", "no_data"}
    bar.counter("errors").click()
    assert seen[-1] == {"incomplete", "mixed", "no_data"}
    bar.counter("ok").click()
    assert seen[-1] == {"ok", "incomplete", "mixed", "no_data"}


def test_summary_whole_construct_and_wrapping_region(qapp, construct, reports):
    from dataclasses import replace

    from i18n import tr
    from ui.construct.summary_bar import SummaryBar

    bar = SummaryBar()
    bar.set_data(construct, None, reports)
    assert tr("construct.region_all") in bar.info_text()
    circ = replace(construct, circular=True)
    bar.set_data(circ, (600, 50), reports)
    assert "601–50" in bar.info_text() and tr("construct.circular") in bar.info_text()


# ------------------------------------------------------------------ matrix


def test_filter_hides_rows(qapp, construct, reports):
    m = _matrix(qapp, construct, reports)
    assert m.visible_clones() == ["c1", "c2", "c3", "c4", "c5"]
    m.set_filter({"errors"})
    assert m.visible_clones() == ["c2", "c3"]
    assert m.row_rect("c2").top() < m.row_rect("c3").top()
    assert m.row_rect("c1").isNull()
    m.set_filter({"ok", "mixed"})
    assert m.visible_clones() == ["c1", "c5"]


def test_hit_test_on_difference_mark(qapp, construct, reports):
    m = _matrix(qapp, construct, reports)
    point = m.mark_point("c2", P)
    assert m.hit_test(point) == ("c2", P)
    assert m.hit_test(m.mark_point("c3", P)) == ("c3", P)
    # the same x on the clean clone, and far from any mark: nothing
    assert m.hit_test(QPoint(point.x(), m.row_rect("c1").center().y())) is None
    away = QPoint(round(m.x_for(500)), m.row_rect("c2").center().y())
    assert m.hit_test(away) is None
    # mixed positions are marks too
    assert m.hit_test(m.mark_point("c5", 300)) == ("c5", 300)


def test_click_on_mark_and_on_row(qapp, construct, reports):
    m = _matrix(qapp, construct, reports)
    diffs, clones = [], []
    m.difference_clicked.connect(lambda c, p: diffs.append((c, p)))
    m.clone_clicked.connect(clones.append)
    left, none = Qt.MouseButton.LeftButton, Qt.MouseButton.NoButton
    point = m.mark_point("c2", P)
    _send_mouse(m, QEvent.Type.MouseButtonPress, point, left, left)
    _send_mouse(m, QEvent.Type.MouseButtonRelease, point, left, none)
    assert diffs == [("c2", P)]
    row = QPoint(round(m.x_for(500)), m.row_rect("c4").center().y())
    _send_mouse(m, QEvent.Type.MouseButtonPress, row, left, left)
    _send_mouse(m, QEvent.Type.MouseButtonRelease, row, left, none)
    assert clones == ["c4"]


def test_shift_drag_on_reference_sets_region(qapp, construct, reports):
    m = _matrix(qapp, construct, reports)
    regions = []
    m.region_changed.connect(regions.append)
    y = m.reference_rect().center().y()
    a, b = QPoint(round(m.x_for(300)), y), QPoint(round(m.x_for(100)), y)
    left, none = Qt.MouseButton.LeftButton, Qt.MouseButton.NoButton
    shift = Qt.KeyboardModifier.ShiftModifier
    _send_mouse(m, QEvent.Type.MouseButtonPress, a, left, left, shift)
    _send_mouse(m, QEvent.Type.MouseMove, QPoint((a.x() + b.x()) // 2, y), none, left, shift)
    _send_mouse(m, QEvent.Type.MouseMove, b, none, left, shift)
    _send_mouse(m, QEvent.Type.MouseButtonRelease, b, left, none, shift)
    assert len(regions) == 1
    start, end = regions[0]
    assert abs(start - 100) <= 2 and abs(end - 300) <= 2
    assert m.view_range() == (0.0, float(len(construct.seq)))   # no pan happened


def test_wheel_zooms_and_drag_pans(qapp, construct, reports):
    m = _matrix(qapp, construct, reports)
    L = float(len(construct.seq))
    assert m.view_range() == (0.0, L)
    centre = QPointF(m.x_for(329), m.row_rect("c1").center().y())
    wheel = QWheelEvent(centre, m.mapToGlobal(centre), QPoint(), QPoint(0, 240),
                        Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier,
                        Qt.ScrollPhase.NoScrollPhase, False)
    QApplication.sendEvent(m, wheel)
    start, span = m.view_range()
    assert span < L and start > 0
    assert abs(m.x_for(329) - centre.x()) < 1.5   # zoom keeps the point under the cursor
    left, none = Qt.MouseButton.LeftButton, Qt.MouseButton.NoButton
    p = centre.toPoint()
    _send_mouse(m, QEvent.Type.MouseButtonPress, p, left, left)
    _send_mouse(m, QEvent.Type.MouseMove, p + QPoint(60, 0), none, left)
    _send_mouse(m, QEvent.Type.MouseButtonRelease, p + QPoint(60, 0), left, none)
    assert m.view_range()[0] < start and m.view_range()[1] == span
    # zooming out again never shows more than the construct
    out = QWheelEvent(centre, m.mapToGlobal(centre), QPoint(), QPoint(0, -2400),
                      Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier,
                      Qt.ScrollPhase.NoScrollPhase, False)
    QApplication.sendEvent(m, out)
    assert m.view_range() == (0.0, L)


@pytest.mark.parametrize("theme", ["light", "dark"])
def test_renders_in_both_themes(qapp, construct, reports, theme):
    from ui.construct.matrix import mark_colors
    from ui.status import status_color

    m = _matrix(qapp, construct, reports, theme=theme)
    image = QImage(m.size(), QImage.Format.Format_ARGB32)
    image.fill(QColor("#123456"))
    m.render(image)
    assert image.pixelColor(2, 2) != QColor("#123456")   # background painted
    point = m.mark_point("c2", P)
    assert image.pixelColor(point) == mark_colors(theme)["confirmed"]
    chip = m.chip_rect("c1")   # sample the pill's padding, left of the text
    assert image.pixelColor(chip.left() + 8, chip.center().y()) == status_color("ok", theme)


def test_wrapping_region_draws_two_segments():
    from ui.construct.matrix import segments

    assert segments(10, 40, 100) == [(10, 40)]
    assert segments(90, 5, 100) == [(90, 100), (0, 5)]
    assert segments(0, 100, 100) == [(0, 100)]


def _synthetic(clone: str, verdict: str, length: int) -> CloneReport:
    coverage = {p: ((1, 1) if p % 50 else (1, 0)) for p in range(length)}
    diffs = (CloneDifference(100, "substitution", "A", "G", 2, 0, 0, "confirmed", True),)
    return CloneReport(clone, (), coverage, diffs if verdict == "errors" else (), (), verdict)


def test_twelve_rows_fit_in_600_px(qapp, construct):
    from ui.construct.matrix import CloneMatrix

    verdicts = ["ok", "errors", "incomplete", "mixed", "no_data", "ok"] * 2
    reps = [_synthetic(f"clone_{i + 1}", v, len(construct.seq)) for i, v in enumerate(verdicts)]
    m = CloneMatrix()
    m.set_data(construct, None, reps)
    assert m.sizeHint().height() <= 600
    m.resize(900, 600)
    assert m.row_rect("clone_12").bottom() < 600
    assert m.hit_test(m.mark_point("clone_2", 100)) == ("clone_2", 100)
