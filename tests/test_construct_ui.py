"""Construct screen widgets: summary bar and clone matrix (offscreen Qt).

Reports come from ``verify_all`` on the real pair ``Achl_ACHLO006-09`` F+R
against its COI reference, with errors / trimming / double peaks introduced
through editing sessions (``write_read``), as in ``test_clone_verify``.
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from PySide6.QtCore import QEvent, QPoint, QPointF, Qt
from PySide6.QtGui import QAction, QColor, QImage, QMouseEvent, QWheelEvent
from PySide6.QtWidgets import QApplication

from chrom_core.align import align_to_construct
from chrom_core.clones import (
    CloneDifference, CloneReport, Grouping, ReadRef, VerifyParams, verify_all,
)
from chrom_core.editing import EditableRead
from chrom_core.reference import load_construct
from clone_fixtures import COI_GB, DATA, perf_factor, real_read, with_secondary_peak, write_read

F_NAME = "Achl_ACHLO006-09_1_F.ab1"
R_NAME = "Achl_ACHLO006-09_2_R.ab1"
F = DATA / "real" / F_NAME
R = DATA / "real" / R_NAME
REGION = (6, 647)
P = 250          # substitution site inside the region
OUTSIDE = 652    # covered by F only, outside the region
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
    # with the edit (no artificial double peak); same at OUTSIDE
    chrom = real_read(F_NAME)
    o = _f_index(construct, OUTSIDE)
    edits = {i: "C" if chrom.seq[i] != "C" else "G", o: "C" if chrom.seq[o] != "C" else "G"}
    for index, base in edits.items():
        chrom = with_secondary_peak(chrom, index, base, 4.0)
    edited = EditableRead(chrom)
    for index, base in edits.items():
        edited.set_base(index, base)
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
    assert [d.ref_pos for d in out[1].differences if not d.in_region] == [OUTSIDE]
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
    image = _render(m, QColor("#123456"))
    assert image.pixelColor(2, 2) != QColor("#123456")   # background painted
    point = m.mark_point("c2", P)
    assert image.pixelColor(point) == mark_colors()["confirmed"]
    chip = m.chip_rect("c1")   # sample the pill's padding, left of the text
    assert image.pixelColor(chip.left() + 8, chip.center().y()) == status_color("ok", theme)


def _render(m, fill=QColor("#123456")) -> QImage:
    image = QImage(m.size(), QImage.Format.Format_ARGB32)
    image.fill(fill)
    m.render(image)
    return image


def test_disputed_mark_is_hollow(qapp, construct, reports):
    from ui.construct.matrix import mark_colors

    m = _matrix(qapp, construct, reports)
    image = _render(m)
    red = mark_colors()["disputed"]
    centre = m.mark_point("c3", P)
    assert image.pixelColor(centre) != red                       # hollow inside
    edges = [image.pixelColor(centre + QPoint(dx, 0)) for dx in range(-4, 5)]
    assert edges.count(red) == 2                                 # left and right outline
    assert sum(1 for dx in range(-4, 5) if image.pixelColor(centre + QPoint(dx, 0)) != red
               and abs(dx) <= 1) == 3                            # >= 5 px wide outline
    assert image.pixelColor(m.mark_point("c2", P)) == mark_colors()["confirmed"]   # solid


def test_out_of_region_difference_drawn_and_clickable(qapp, construct, reports):
    from ui.theme import theme as palette

    m = _matrix(qapp, construct, reports)
    image = _render(m)
    point = m.mark_point("c2", OUTSIDE)
    colour = image.pixelColor(point)
    assert colour.red() > colour.green() + 40 and colour.red() > colour.blue() + 40
    assert m.hit_test(point) == ("c2", OUTSIDE)
    clicked = []
    m.difference_clicked.connect(lambda c, p: clicked.append((c, p)))
    left, none = Qt.MouseButton.LeftButton, Qt.MouseButton.NoButton
    _send_mouse(m, QEvent.Type.MouseButtonPress, point, left, left)
    _send_mouse(m, QEvent.Type.MouseButtonRelease, point, left, none)
    assert clicked == [("c2", OUTSIDE)]
    # outside the region the track is hatched ("not checked"); a gap in the
    # region's coverage stays the plain track colour ("no coverage")
    track = QColor(palette("light")["track"])
    y = m.row_rect("c1").center().y()
    hatched = [image.pixelColor(x, y + dy) for dx in range(10)
               for x in [round(m.x_for(647)) + 2 + dx] for dy in range(-4, 5)
               if x < round(m.x_for(len(construct.seq)))]
    assert hatched and any(c != track for c in hatched)
    y4 = m.row_rect("c4").center().y()
    assert all(image.pixelColor(x, y4) == track
               for x in range(round(m.x_for(450)), round(m.x_for(470))))


def test_hit_test_while_zoomed(qapp, construct, reports):
    m = _matrix(qapp, construct, reports)
    anchor = QPointF(m.x_for(P + 0.5), m.row_rect("c2").center().y())
    zoom = QWheelEvent(anchor, m.mapToGlobal(anchor), QPoint(), QPoint(0, 120 * 12),
                       Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier,
                       Qt.ScrollPhase.NoScrollPhase, False)
    QApplication.sendEvent(m, zoom)
    assert m.view_range()[1] < 100
    point = m.mark_point("c2", P)
    assert m.hit_test(point) == ("c2", P)
    neighbour = QPoint(round(m.x_for(P + 3.5)), point.y())     # three bases away
    assert m.hit_test(neighbour) is None
    assert m.hit_test(QPoint(point.x(), m.row_rect("c1").center().y())) is None


def test_ruler_labels_are_one_based(qapp, construct, reports):
    m = _matrix(qapp, construct, reports)
    ticks = m.ruler_ticks()
    assert ticks and all(v % ticks[0][0] == 0 for v, _ in ticks)
    for v, x in ticks:
        assert int(m.pos_at(x)) + 1 == v   # the tick sits on base number v (1-based)


def test_mark_at_last_position(qapp, construct):
    from ui.construct.matrix import CloneMatrix, mark_colors

    last = len(construct.seq) - 1
    rep = _synthetic("end", "errors", len(construct.seq))
    rep = CloneReport(rep.clone, (), rep.coverage, (
        CloneDifference(last, "substitution", "A", "G", 2, 0, 0, "confirmed", True),), (), "errors")
    m = CloneMatrix()
    m.set_data(construct, None, [rep])
    m.resize(900, m.sizeHint().height())
    point = m.mark_point("end", last)
    assert point.x() < m.chip_rect("end").left()
    assert m.hit_test(point) == ("end", last)
    assert _render(m).pixelColor(point) == mark_colors()["confirmed"]


# measured ~12 ms (best of 5) on the dev machine for 12 rows x 200 runs on 10 kb;
# the bound is about 2x that (it was ~440 ms before name widths were cached)
PAINT_BUDGET_S = 0.025


@pytest.mark.perf
def test_paint_cost_10kb_12_rows(qapp):
    import time

    from chrom_core.reference import Construct
    from ui.construct.matrix import CloneMatrix

    length = 10_000
    big = Construct("big", "ACGT" * (length // 4), False, ())
    coverage = {p: ((1, 1) if (p // 50) % 2 else (1, 0)) for p in range(length)}   # 200 runs
    reps = [CloneReport(f"clone_{i}", (), coverage, (), (), "ok") for i in range(12)]
    m = CloneMatrix()
    m.set_data(big, None, reps)
    m.resize(1200, m.sizeHint().height())
    image = QImage(m.size(), QImage.Format.Format_ARGB32)
    m.render(image)   # warm-up
    best = min(_timed(m, image, time) for _ in range(3))
    assert best < PAINT_BUDGET_S * perf_factor(), f"paint took {best * 1000:.0f} ms"


def _timed(m, image, time) -> float:
    start = time.perf_counter()
    m.render(image)
    return time.perf_counter() - start


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


# ------------------------------------------------------------------ detail pane


@pytest.fixture(scope="module")
def pane_report(construct, tmp_path_factory):
    """ACHLO006-09: F (+) and R (-) with a substitution at P in F's session."""
    folder = tmp_path_factory.mktemp("pane")
    i = _f_index(construct, P)
    chrom = real_read(F_NAME)
    base = "C" if chrom.seq[i] != "C" else "G"
    chrom = with_secondary_peak(chrom, i, base, 4.0)
    edited = EditableRead(chrom)
    edited.set_base(i, base)
    f_sub = write_read(edited.current(), folder, F_NAME)
    out = verify_all(construct, REGION, Grouping([ReadRef(f_sub, "c1", "F"), ReadRef(R, "c1", "R")]))
    report = out[0]
    assert [r.alignment.strand for r in report.reads] == ["+", "-"]
    assert [d.ref_pos for d in report.region_differences] == [P]
    return report


def _pane(qapp, construct, report):
    from ui.construct.detail_pane import DetailPane

    pane = DetailPane()
    pane.resize(1000, 500)
    pane.show_clone(report, construct)
    return pane


def test_displayed_index_converts_once_for_both_strands():
    from chrom_core.align import ReadAlignment
    from ui.construct.detail_pane import displayed_index

    def al(strand, mapping):
        return ReadAlignment(strand=strand, score=0, ref_start=0, ref_end=10, read_start=0, read_end=10,
                             ref_to_read=mapping, differences=[], identity=1.0, alternatives=1,
                             oriented_read="", blocks=[], ref_length=100)

    plus, minus = al("+", {5: 3, 6: 4}), al("-", {5: 3, 6: 4})
    assert displayed_index(plus, 5, 10) == 3
    assert displayed_index(minus, 5, 10) == 10 - 1 - 3
    assert displayed_index(minus, 50, 10) is None        # not covered
    # a deleted position falls back to the nearest mapped neighbour
    gap = al("+", {5: 3, 7: 4})
    assert displayed_index(gap, 6, 10) == 4


def test_displayed_index_out_of_range_is_not_covered():
    from chrom_core.align import ReadAlignment
    from ui.construct.detail_pane import displayed_index

    def al(strand):
        return ReadAlignment(strand=strand, score=0, ref_start=0, ref_end=10, read_start=0,
                             read_end=10, ref_to_read={5: 8}, differences=[], identity=1.0,
                             alternatives=1, oriented_read="", blocks=[], ref_length=100)

    # the chromatogram was shortened by an edit after the alignment was stored
    assert displayed_index(al("+"), 5, 6) is None      # index 8 >= 6
    assert displayed_index(al("-"), 5, 6) is None      # 6 - 1 - 8 < 0


def test_pane_survives_alignment_past_a_shortened_read(qapp, construct, pane_report):
    from dataclasses import replace

    first = pane_report.reads[0]
    far = {p: i + 100_000 for p, i in first.alignment.ref_to_read.items()}
    stale = replace(pane_report, reads=(
        replace(first, alignment=replace(first.alignment, ref_to_read=far)),
        *pane_report.reads[1:]))
    pane = _pane(qapp, construct, stale)
    pane.show_position(P)
    assert not pane._notes[0].isHidden() and pane._notes[1].isHidden()


def test_pane_one_chromatogram_per_aligned_read(qapp, construct, pane_report):
    pane = _pane(qapp, construct, pane_report)
    assert len(pane.views()) == 2
    assert pane.table.rowCount() == len(pane_report.differences)


def test_pane_centres_every_read_on_the_position(qapp, construct, pane_report):
    pane = _pane(qapp, construct, pane_report)
    pane.show_position(P)
    for result, view in zip(pane_report.reads, pane.views()):
        n = len(view.chrom.seq)
        idx = result.alignment.ref_to_read[P]
        want = idx if result.alignment.strand == "+" else n - 1 - idx
        assert view.highlights == [(want, want + 1)]
        centre = view.offset + view.width() / view.zoom / 2
        assert abs(centre - view.chrom.peak_indices[want]) <= 1.0
        letters, mapping = view.reference
        assert mapping[P] == want
        assert len({v.zoom for v in pane.views()}) == 1   # shared scale


def test_pane_reference_strip_shows_neighbourhood(qapp, construct, pane_report):
    pane = _pane(qapp, construct, pane_report)
    pane.show_position(P)
    text = pane.reference_text()
    assert construct.seq[P - 10:P] in text and construct.seq[P + 1:P + 11] in text
    assert f"[{construct.seq[P]}]" in text


def test_pane_table_row_selection_shows_position(qapp, construct, pane_report):
    pane = _pane(qapp, construct, pane_report)
    seen = []
    orig = pane.show_position
    pane.show_position = lambda p: (seen.append(p), orig(p))
    pane.table.selectRow(0)
    assert seen[-1] == P == pane_report.differences[0].ref_pos
    assert pane.table.item(0, 0).text() == str(P + 1)
    assert "→" in pane.table.item(0, 1).text()


def test_pane_double_click_opens_read(qapp, construct, pane_report):
    pane = _pane(qapp, construct, pane_report)
    opened = []
    pane.open_read.connect(opened.append)
    view = pane.views()[1]
    _send_mouse(view, QEvent.Type.MouseButtonDblClick, QPoint(20, 20),
                Qt.MouseButton.LeftButton, Qt.MouseButton.LeftButton)
    assert opened == [pane_report.reads[1].ref.path]


def test_pane_effect_labels_translated(qapp, construct, pane_report):
    from dataclasses import replace
    from i18n import tr

    d = pane_report.differences[0]
    report = replace(pane_report, differences=(replace(d, effect="frameshift"), replace(d, effect="D45N")))
    pane = _pane(qapp, construct, report)
    assert pane.table.item(0, 2).text() == tr("effect.frameshift")
    assert pane.table.item(1, 2).text() == "D45N"


def test_pane_status_cells_keep_one_colour_per_meaning(qapp, construct, pane_report):
    from dataclasses import replace
    from ui import status

    d = pane_report.differences[0]
    diffs = tuple(replace(d, status=s) for s in ("confirmed", "disputed", "low_quality"))
    pane = _pane(qapp, construct, replace(pane_report, differences=diffs))
    bg = [pane.table.item(r, 4).background().color() for r in range(3)]
    red = status.status_color("errors", "light")
    assert bg[0] == red
    assert bg[1] != red and bg[1].red() > bg[1].green() and bg[1].red() > bg[1].blue()   # weaker red tint
    assert bg[1] != status.status_color("incomplete", "light")
    assert bg[2] == status.status_color("no_data", "light")
    pane.set_theme("dark")
    assert pane.table.item(0, 4).background().color() == status.status_color("errors", "dark")


def test_pane_dims_rows_outside_region(qapp, construct, pane_report):
    from dataclasses import replace
    from i18n import tr

    d = pane_report.differences[0]
    pane = _pane(qapp, construct, replace(pane_report, differences=(d, replace(d, in_region=False))))
    inside, outside = pane.table.item(0, 0), pane.table.item(1, 0)
    assert not inside.font().italic() and inside.toolTip() == ""
    assert outside.font().italic() and outside.toolTip() == tr("matrix.legend.outside")
    assert outside.foreground().color() != inside.foreground().color()


def test_pane_zoom_is_shared_and_kept(qapp, construct, pane_report):
    pane = _pane(qapp, construct, pane_report)
    views = pane.views()
    assert len({v.zoom for v in views}) == 1
    views[0].zoom_in()                           # before any position is shown
    assert views[1].zoom == views[0].zoom
    pane.show_position(P)
    chosen = views[0].zoom
    assert views[1].zoom == chosen               # show_position keeps the user's zoom
    views[0].zoom_in()
    assert views[1].zoom == views[0].zoom > chosen
    for result, view in zip(pane_report.reads, views):
        n = len(view.chrom.seq)
        idx = result.alignment.ref_to_read[P]
        want = idx if result.alignment.strand == "+" else n - 1 - idx
        centre = view.offset + view.width() / view.zoom / 2
        assert abs(centre - view.chrom.peak_indices[want]) <= 1.0
    pane.resize(1100, 500)
    assert views[0].zoom == views[1].zoom and views[0].zoom != chosen


def test_pane_failed_and_ambiguous_reads_are_shown(qapp, construct, pane_report):
    from dataclasses import replace
    from PySide6.QtWidgets import QLabel
    from chrom_core.clones import ReadResult
    from i18n import tr

    f, r = pane_report.reads
    bad = ReadResult(ReadRef(Path("bad_read.ab1"), "c1", "F"), None, error="no_alignment")
    report = replace(pane_report, reads=(f, bad, replace(r, ambiguous=True)))
    pane = _pane(qapp, construct, report)
    assert len(pane.views()) == 2
    texts = [lab.text() for lab in pane.findChildren(QLabel)]
    assert any("bad_read.ab1" in s and tr("detail.err.no_alignment") in s for s in texts)
    assert tr("detail.ambiguous") in pane.views()[1].toolTip()
    assert tr("detail.ambiguous") not in pane.views()[0].toolTip()
    gone = ReadResult(ReadRef(Path("missing.ab1"), "c1", "F"), r.alignment)
    pane = _pane(qapp, construct, replace(pane_report, reads=(f, gone)))
    assert len(pane.views()) == 1
    assert any("missing.ab1" in lab.text() for lab in pane.findChildren(QLabel))


def test_pane_dims_read_that_does_not_cover_position(qapp, construct, pane_report):
    pane = _pane(qapp, construct, pane_report)
    covers = pane_report.reads[1].alignment
    outside = next(p for p in range(len(construct.seq)) if not covers.covers(p)
                   and pane_report.reads[0].alignment.covers(p))
    pane.show_position(outside)
    assert pane.views()[1].graphicsEffect() is not None
    assert pane.views()[1].highlights == []
    assert pane.views()[0].graphicsEffect() is None
    pane.show_position(P)
    assert pane.views()[1].graphicsEffect() is None


def test_pane_renders_png(qapp, construct, pane_report):
    target = os.environ.get("AB1KIT_DETAIL_PNG")
    if not target:
        pytest.skip("set AB1KIT_DETAIL_PNG to render the pane")
    pane = _pane(qapp, construct, pane_report)
    pane.show_position(P)
    assert pane.grab().save(target)


def test_pane_survives_corrupt_ab1(qapp, construct, pane_report, tmp_path):
    from dataclasses import replace
    from PySide6.QtWidgets import QLabel
    from chrom_core.abif import read_ab1
    from chrom_core.clones import ReadResult

    raw = bytearray((Path(__file__).resolve().parents[1] / "examples" / "310.ab1").read_bytes())
    directory = int.from_bytes(raw[26:30], "big")
    raw[directory + 12] ^= 0xFF          # deterministic: BioPython's ABIF parser asserts
    corrupt = tmp_path / "corrupt.ab1"
    corrupt.write_bytes(raw)
    with pytest.raises(AssertionError):
        read_ab1(corrupt)

    f, r = pane_report.reads
    broken = ReadResult(ReadRef(corrupt, "c1", "F"), f.alignment)
    pane = _pane(qapp, construct, replace(pane_report, reads=(broken, r)))
    assert len(pane.views()) == 1
    assert any("corrupt.ab1" in lab.text() for lab in pane.findChildren(QLabel))
    pane.show_position(P)
    assert pane.views()[0].highlights


# ------------------------------------------------------------------ workspace (Task 16)

REAL_DIR = DATA / "real"


@pytest.fixture
def window(qapp, monkeypatch):
    """Main window with config writes disabled and message boxes recorded."""
    from PySide6.QtWidgets import QMessageBox

    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    for target in ("ui.main_window.save_config", "ui.read_view.read_view.save_config",
                   "ui.construct.workspace.save_config"):
        monkeypatch.setattr(target, lambda *a, **k: None)
    messages: list[str] = []

    def record(_parent, title, text, *a, **k):
        messages.append(f"{title}: {text}")
        return QMessageBox.StandardButton.Ok

    monkeypatch.setattr(QMessageBox, "warning", record)
    monkeypatch.setattr(QMessageBox, "critical", record)
    w = MainWindow(AppConfig(recent_files=[]))
    w.messages = messages
    yield w
    w._cancel_tasks()
    w.close()
    _delete(qapp, w)


def _delete(qapp, widget):
    """Really destroy a test window (closed windows otherwise live on and make
    every later theme switch restyle them)."""
    from PySide6.QtCore import QCoreApplication, QEvent

    widget.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete.value)
    qapp.processEvents()


def _wait_idle(qapp, ws, timeout=30.0):
    import time

    deadline = time.monotonic() + timeout
    while ws.is_verifying() and time.monotonic() < deadline:
        qapp.processEvents()
        time.sleep(0.005)
    qapp.processEvents()
    assert not ws.is_verifying()


def _drop(widget, paths):
    from PySide6.QtCore import QMimeData, QUrl
    from PySide6.QtGui import QDropEvent

    mime = QMimeData()
    mime.setUrls([QUrl.fromLocalFile(str(p)) for p in paths])
    event = QDropEvent(QPointF(10, 10), Qt.DropAction.CopyAction, mime,
                       Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier)
    widget.dropEvent(event)


def _count_starts(ws):
    calls = []
    original = ws._start_task

    def counting(*a, **k):
        calls.append(k.get("key"))
        return original(*a, **k)

    ws._start_task = counting
    return calls


def _loaded(qapp, window, construct_path=COI_GB, reads=(REAL_DIR,)):
    ws = window.workspace
    _drop(ws, [construct_path, *reads])
    _wait_idle(qapp, ws)
    return ws


def test_workspace_is_first_tab_and_shows_empty_state(window):
    ws = window.workspace
    assert window.tabs.indexOf(ws) == 0
    assert window.tabs.currentWidget() is ws
    assert window.tabs.indexOf(window.read_view) == 1
    window.show()
    assert ws.empty_state.isVisible()
    assert not ws.matrix.isVisible()
    texts = " ".join(z.text() for z in ws.empty_state.zones())
    assert "FASTA" in texts and ".ab1" in texts


def test_drop_construct_and_reads_folder(qapp, window, construct):
    ws = _loaded(qapp, window)
    window.show()
    assert not ws.empty_state.isVisible() and ws.matrix.isVisible()
    expected = verify_all(construct, None, ws.grouping())
    assert ws.matrix.visible_clones() == [r.clone for r in expected]
    assert len(ws.matrix.visible_clones()) == 4
    assert [r.verdict for r in ws.reports()] == [r.verdict for r in expected]
    # same-lineage clones against the unmodified reference: no errors
    same = [r for r in ws.reports() if r.clone.startswith(("Achl_ACHLO006", "Achl_ACHLO007"))]
    assert len(same) == 2 and all(r.verdict in ("ok", "incomplete") for r in same)
    assert ws.summary.title_text() == construct.name
    assert window.messages == []


def test_matrix_click_shows_details_and_open_read_switches_tab(qapp, window, tmp_path, construct):
    seq = construct.seq
    mutant = "G" if seq[P] != "G" else "T"
    fasta = tmp_path / "mut.fa"
    fasta.write_text(f">mut\n{seq[:P]}{mutant}{seq[P + 1:]}\n")
    f_copy = tmp_path / F_NAME
    f_copy.write_bytes(F.read_bytes())
    ws = _loaded(qapp, window, fasta, [f_copy])
    ws.set_region((100, 600))
    _wait_idle(qapp, ws)
    (report,) = ws.reports()
    assert report.verdict == "errors" and [d.ref_pos for d in report.region_differences] == [P]

    ws.matrix.difference_clicked.emit(report.clone, P)
    assert len(ws.detail.views()) == 1

    ws.detail.open_read.emit(f_copy)
    assert window.tabs.currentWidget() is window.read_view
    assert window.read_model is not None

    # edit the read at the difference: the check re-runs and the verdict changes
    window.apply_base_change(_f_index(construct, P), mutant)
    assert ws.is_verifying()
    _wait_idle(qapp, ws)
    assert [r.verdict for r in ws.reports()] == ["ok"]

    window.undo_read()
    _wait_idle(qapp, ws)
    assert [r.verdict for r in ws.reports()] == ["errors"]


def test_read_edit_outside_project_does_not_reverify(qapp, window):
    ws = _loaded(qapp, window)
    starts = _count_starts(ws)
    examples = Path(__file__).resolve().parents[1] / "examples"
    window.read_view.read_edited.emit(examples / "310.ab1")
    assert starts == []


def test_project_roundtrip_restores_region_grouping_params(qapp, window, tmp_path):
    ws = _loaded(qapp, window)
    clone = ws.grouping().reads[0]
    grouping = ws.grouping()
    grouping.assign(clone.path, "renamed", clone.primer)
    ws.set_grouping(grouping)
    ws.set_region((100, 600))
    ws.controls.min_q.setValue(25)
    _wait_idle(qapp, ws)
    project = tmp_path / "check.ab1proj"
    assert ws.save_project(project)

    ws.clear()
    assert ws.reports() == [] and ws.region() is None
    assert ws.load_project(project)
    _wait_idle(qapp, ws)
    assert ws.region() == (100, 600)
    assert (ws.controls.start_spin.value(), ws.controls.end_spin.value()) == (101, 600)
    assert ws.grouping().to_dict() == grouping.to_dict()
    assert ws.params().min_q == 25 and ws.controls.min_q.value() == 25
    assert "renamed" in ws.matrix.visible_clones()
    assert window.messages == []


def test_broken_project_shows_message_and_empty_workspace(qapp, window, tmp_path):
    ws = _loaded(qapp, window)
    bad = tmp_path / "bad.ab1proj"
    bad.write_text("{not json", encoding="utf-8")
    assert ws.load_project(bad) is False
    assert len(window.messages) == 1 and "bad.ab1proj" in window.messages[0]
    window.show()
    assert ws.empty_state.isVisible()
    assert ws.reports() == [] and ws.grouping().reads == []


def test_invalid_region_shows_message_and_does_not_run(qapp, window):
    ws = _loaded(qapp, window)       # COI .gb is linear
    starts = _count_starts(ws)
    assert ws.set_region((500, 100)) is False
    assert starts == [] and ws.region() is None
    assert len(window.messages) == 1
    assert (ws.controls.start_spin.value(), ws.controls.end_spin.value()) == (1, 658)
    # the same through the spin boxes (1-based, inclusive): 500..658 is fine ...
    ws.controls.start_spin.setValue(500)
    assert starts == ["verify"] and ws.region() == (499, 658)
    # ... 500..100 wraps the origin of a linear construct: refused, spins restored
    ws.controls.end_spin.setValue(100)
    assert starts == ["verify"] and ws.region() == (499, 658) and len(window.messages) == 2
    assert (ws.controls.start_spin.value(), ws.controls.end_spin.value()) == (500, 658)
    _wait_idle(qapp, ws)


def test_region_spins_and_matrix_drag_stay_in_sync(qapp, window):
    ws = _loaded(qapp, window)
    starts = _count_starts(ws)
    ws.matrix.region_changed.emit((50, 400))
    assert ws.region() == (50, 400) and starts == ["verify"]
    assert (ws.controls.start_spin.value(), ws.controls.end_spin.value()) == (51, 400)
    ws.controls.whole_button.click()
    assert ws.region() is None and starts == ["verify", "verify"]
    _wait_idle(qapp, ws)


def test_cancelled_reverification_keeps_previous_verdicts(qapp, window):
    ws = _loaded(qapp, window)
    before = ws.reports()
    shown = ws.matrix.visible_clones()
    info = ws.summary.info_text()
    ws.set_region((100, 600))
    assert ws.is_verifying()
    window._cancel_tasks()
    _wait_idle(qapp, ws)
    assert ws.reports() is before
    assert ws.matrix.visible_clones() == shown and ws.summary.info_text() == info


def test_failed_reverification_keeps_previous_verdicts(qapp, window, monkeypatch):
    ws = _loaded(qapp, window)
    before = ws.reports()

    def boom(*a, **k):
        raise RuntimeError("boom")

    monkeypatch.setattr("ui.construct.workspace.verify_all", boom)
    ws.controls.min_q.setValue(30)
    _wait_idle(qapp, ws)
    assert ws.reports() is before
    assert any("boom" in m for m in window.messages)


def test_params_and_circular_trigger_reverification(qapp, window):
    ws = _loaded(qapp, window)
    starts = _count_starts(ws)
    ws.controls.het_ratio.setValue(0.5)
    ws.controls.mixed_min.setValue(4)
    ws.controls.circular.setChecked(True)
    assert starts == ["verify"] * 3
    _wait_idle(qapp, ws)
    assert ws.params().het_ratio == 0.5 and ws.params().mixed_min_positions == 4
    assert ws.construct().circular is True


def test_export_report_csv_and_xlsx(qapp, window, tmp_path):
    ws = _loaded(qapp, window)
    assert ws.export_report(tmp_path / "report.csv")
    assert (tmp_path / "report.csv").read_text(encoding="utf-8").startswith("clone,verdict")
    assert (tmp_path / "report_differences.csv").read_text(encoding="utf-8").startswith("clone,")
    assert ws.export_report(tmp_path / "report.xlsx")
    assert (tmp_path / "report.xlsx").stat().st_size > 0


def test_grouping_dialog_pattern_preview_and_invalid_regex(qapp):
    from PySide6.QtWidgets import QDialogButtonBox

    from chrom_core.clones import group_reads
    from ui.construct.grouping_dialog import GroupingDialog

    grouping = group_reads(sorted(REAL_DIR.glob("*.ab1")))
    dialog = GroupingDialog(grouping, "", None)
    ok = dialog.buttons.button(QDialogButtonBox.StandardButton.Ok)
    assert dialog.table.rowCount() == 8 and ok.isEnabled()

    dialog.pattern_edit.setText("([")
    assert not ok.isEnabled() and dialog.pattern_invalid()
    dialog.pattern_edit.setText("(.*)")              # no named groups
    assert not ok.isEnabled() and dialog.pattern_invalid()

    pattern = r"Achl_ACHLO00(?P<clone>\d-09)_\d_(?P<primer>[FR])"
    dialog.pattern_edit.setText(pattern)
    assert ok.isEnabled() and not dialog.pattern_invalid()
    result = dialog.result_grouping()
    assert sorted(result.clones()) == ["6-09", "7-09"]
    assert len(result.unassigned) == 4
    empty_rows = [r for r in range(8) if dialog.table.item(r, 1).text() == ""]
    assert dialog.unassigned_rows() == empty_rows and len(empty_rows) == 4
    assert dialog.result_pattern() == pattern

    # manual edit: clear a clone cell -> that read becomes unassigned
    first = next(r for r in range(8) if dialog.table.item(r, 1).text())
    dialog.table.item(first, 1).setText("")
    assert first in dialog.unassigned_rows()
    assert len(dialog.result_grouping().unassigned) == 5


def test_new_ui_keys_exist_in_both_languages():
    import re as _re

    import i18n

    root = Path(__file__).resolve().parents[1]
    files = [root / "ui" / "main_window.py", *sorted((root / "ui" / "construct").glob("*.py"))]
    keys = set()
    for path in files:
        text = path.read_text(encoding="utf-8")
        keys |= set(_re.findall(r"""\btr\(\s*["']([a-z_.]+)["']""", text))
    assert {"tab.construct", "menu.project.open", "workspace.drop_construct"} <= keys
    assert not (keys - set(i18n.RU)), keys - set(i18n.RU)
    assert not (keys - set(i18n.EN)), keys - set(i18n.EN)


# ------------------------------------------------------------------ Task 16 review fixes


def test_hidden_viewer_actions_do_nothing_on_construct_tab(qapp, window, tmp_path, construct):
    f_copy = tmp_path / F_NAME
    f_copy.write_bytes(F.read_bytes())
    ws = _loaded(qapp, window, COI_GB, [f_copy])
    ws.detail.open_read.emit(f_copy)
    assert window.tabs.currentWidget() is window.read_view
    window.apply_base_change(_f_index(construct, P), "N")
    _wait_idle(qapp, ws)
    edited = "".join(window.read_model.seq)
    assert window.undo_action.isEnabled()

    window.tabs.setCurrentWidget(ws)
    starts = _count_starts(ws)
    assert not window.undo_action.isEnabled()
    window.undo_action.trigger()
    for action in window.read_view.findChildren(QAction):
        if action.parent() is window.read_view:
            assert not action.isEnabled(), action.text()
    assert "".join(window.read_model.seq) == edited and starts == []

    window.tabs.setCurrentWidget(window.read_view)
    assert window.undo_action.isEnabled() and not window.redo_action.isEnabled()


def test_clear_resets_matrix_and_summary(qapp, window, tmp_path):
    ws = _loaded(qapp, window)
    assert ws.matrix.visible_clones()
    ws.clear()
    assert ws.matrix.visible_clones() == [] and ws.summary.title_text() == ""
    assert not ws.detail.isVisible()
    _drop(ws, [COI_GB, REAL_DIR])
    assert ws.is_verifying()
    window._cancel_tasks()
    _wait_idle(qapp, ws)
    assert ws.matrix.visible_clones() == [] and ws.reports() == []
    assert ws.summary.title_text() == ""


def test_project_with_missing_construct_is_not_reported_open(qapp, window, tmp_path):
    from chrom_core.project import Project, save_project

    project = tmp_path / "lost.ab1proj"
    save_project(project, Project(tmp_path / "gone.gb", False, None,
                                  Grouping([ReadRef(F, "c1", "F")]), VerifyParams()))
    said = []
    window.workspace.status_message.connect(lambda text, _ms: said.append(text))
    assert window.workspace.load_project(project) is False
    assert window.workspace.project_path() is None
    assert len(window.messages) == 1 and "gone.gb" in window.messages[0]
    assert not any("lost.ab1proj" in s for s in said)


def test_drop_without_changes_does_not_reverify(qapp, window, tmp_path):
    ws = _loaded(qapp, window)
    starts = _count_starts(ws)
    _drop(ws, [REAL_DIR])                      # nothing new
    bad = tmp_path / "broken.gb"
    bad.write_text("LOCUS nonsense\n", encoding="utf-8")
    _drop(ws, [bad])                           # construct fails to load
    assert starts == []
    assert len(window.messages) == 1


def test_wrapping_region_on_circular_construct_survives_project(qapp, window, tmp_path, construct):
    fasta = tmp_path / "coi.fa"                # FASTA constructs are circular
    fasta.write_text(f">coi\n{construct.seq}\n")
    ws = _loaded(qapp, window, fasta, [REAL_DIR])
    assert ws.construct().circular
    ws.controls.start_spin.setValue(600)
    ws.controls.end_spin.setValue(50)
    assert ws.region() == (599, 50)
    _wait_idle(qapp, ws)
    project = tmp_path / "wrap.ab1proj"
    assert ws.save_project(project)
    ws.clear()
    assert ws.load_project(project)
    _wait_idle(qapp, ws)
    assert ws.region() == (599, 50) and ws.construct().circular
    assert (ws.controls.start_spin.value(), ws.controls.end_spin.value()) == (600, 50)
    assert window.messages == []


def test_grouping_dialog_flags_invalid_initial_pattern(qapp):
    from PySide6.QtWidgets import QDialogButtonBox

    from chrom_core.clones import group_reads
    from ui.construct.grouping_dialog import GroupingDialog

    grouping = group_reads(sorted(REAL_DIR.glob("*.ab1")))
    dialog = GroupingDialog(grouping, "(?P<clone>[", None)
    assert dialog.pattern_invalid()
    assert not dialog.buttons.button(QDialogButtonBox.StandardButton.Ok).isEnabled()
    assert dialog.table.rowCount() == 8 and dialog.result_grouping().to_dict() == grouping.to_dict()


def test_double_spin_follows_ui_language(qapp):
    from i18n import set_language
    from ui.construct.controls import ControlsBar

    try:
        set_language("en")
        assert ControlsBar().het_ratio.text() == "0.40"
        set_language("ru")
        assert ControlsBar().het_ratio.text() == "0,40"
    finally:
        set_language("ru")
