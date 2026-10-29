"""Per-clone verification: coverage, merged differences and verdicts.

The base clone is the real pair ``Achl_ACHLO006-09`` F+R against its own COI
reference (HM417896, 658 bp). Errors are introduced into the *reference*
(``mutate``), so the unmodified reads report them, or into a read through an
editing session (``write_read``), which is how a user's edits reach the check.
"""

from __future__ import annotations

import random
import time
from concurrent.futures import CancelledError
from dataclasses import replace

import pytest

from chrom_core import clones
from chrom_core.clones import (
    MIN_IDENTITY,
    CloneReport,
    VerifyParams,
    load_read,
    region_positions,
    verify_all,
    verify_clone,
)
from chrom_core.align import align_to_construct
from chrom_core.clones import Grouping, ReadRef
from chrom_core.editing import EditableRead
from chrom_core.peaks import secondary_peak
from chrom_core.reference import Construct, load_construct
from clone_fixtures import (
    DATA,
    PUC19_GB,
    coi_reference,
    mutate,
    real_read,
    with_secondary_peak,
    write_read,
)

F_NAME = "Achl_ACHLO006-09_1_F.ab1"
R_NAME = "Achl_ACHLO006-09_2_R.ab1"
F_PATH = DATA / "real" / F_NAME
R_PATH = DATA / "real" / R_NAME
# The part of COI covered by both trimmed reads: F aligns 6..658, R 0..647.
REGION = (6, 647)
P = 250            # inside the region; context AAC[A]TGA, no homopolymer
OUTSIDE = 652      # covered by F only, outside the region


@pytest.fixture(scope="module")
def coi() -> str:
    return coi_reference()


def _reads(f=F_PATH, r=R_PATH) -> list[ReadRef]:
    out = []
    if f is not None:
        out.append(ReadRef(f, "c1", "F"))
    if r is not None:
        out.append(ReadRef(r, "c1", "R"))
    return out


def _construct(seq: str, circular: bool = False) -> Construct:
    return Construct("c", seq, circular)


def _f_index(coi: str, ref_pos: int) -> int:
    """Displayed index in the F read of reference position ``ref_pos``."""
    al = align_to_construct(real_read(F_NAME).seq, _construct(coi))
    assert al.strand == "+"
    return al.ref_to_read[ref_pos]


def test_ok_clone(coi):
    report = verify_clone(_construct(coi), REGION, _reads())
    assert report.verdict == "ok" and report.region_differences == ()
    assert report.region_coverage_pct == 100.0
    assert 0 < report.both_strands_pct <= 100.0
    assert set(report.coverage) == set(range(*REGION))
    assert all(r.error == "" and not r.ambiguous for r in report.reads)


def test_substitution_found_with_support(coi):
    ref = mutate(coi, P, "sub", "G")
    r = verify_clone(_construct(ref), REGION, _reads())
    d, = r.region_differences
    assert (d.ref_pos, d.kind, d.status, d.support, d.against) == (P, "substitution", "confirmed", 2, 0)
    assert (d.ref_base, d.read_base, d.in_region) == ("G", "A", True)
    assert r.verdict == "errors"


def test_insertion_and_deletion_found(coi):
    # base missing from the reference -> the reads show an insertion
    r = verify_clone(_construct(mutate(coi, P, "del")), REGION, _reads())
    d, = r.region_differences
    assert (d.ref_pos, d.kind, d.read_base, d.status, d.support) == (P, "insertion", "A", "confirmed", 2)
    assert r.verdict == "errors"
    # extra base in the reference -> the reads show a deletion
    r = verify_clone(_construct(mutate(coi, P, "ins", "G")), REGION, _reads())
    d, = r.region_differences
    assert (d.ref_pos, d.kind, d.ref_base, d.status, d.support) == (P, "deletion", "G", "confirmed", 2)
    assert r.verdict == "errors"


def test_disputed_when_reads_disagree(coi, tmp_path):
    # A C peak 4x the called A under P, then the call edited to C: the trace
    # agrees with the edit (no artificial double peak at P).
    i = _f_index(coi, P)
    model = EditableRead(with_secondary_peak(real_read(F_NAME), i, "C", 4.0))
    model.set_base(i, "C")
    assert secondary_peak(model.current(), i).ratio < VerifyParams().het_ratio
    f = write_read(model.current(), tmp_path, F_NAME)
    r = verify_clone(_construct(coi), REGION, _reads(f=f))
    d, = r.region_differences
    assert (d.ref_pos, d.kind, d.read_base, d.support, d.against, d.status) == (
        P, "substitution", "C", 1, 1, "disputed")
    assert r.verdict == "errors"


def test_outside_region_does_not_affect_verdict(coi):
    r = verify_clone(_construct(mutate(coi, OUTSIDE, "sub", "A")), REGION, _reads())
    assert r.verdict == "ok" and r.region_differences == ()
    d, = r.differences
    assert (d.ref_pos, d.in_region, d.status) == (OUTSIDE, False, "confirmed")


def test_incomplete_when_read_trimmed(coi, tmp_path):
    model = EditableRead(real_read(F_NAME))
    model.apply_trim(100, 400)
    f = write_read(model.current(), tmp_path, F_NAME)
    r = verify_clone(_construct(coi), REGION, _reads(f=f, r=None))
    assert r.verdict == "incomplete"
    assert r.region_coverage_pct < 100
    assert r.both_strands_pct == 0.0


def _mixed_read(coi, tmp_path, ratio: float, positions) -> object:
    chrom = real_read(F_NAME)
    for pos in positions:
        i = _f_index(coi, pos)
        other = next(b for b in "ACGT" if b != chrom.seq[i])
        chrom = with_secondary_peak(chrom, i, other, ratio)
        assert secondary_peak(chrom, i).ratio == pytest.approx(ratio, abs=0.05)
    return write_read(chrom, tmp_path, F_NAME)


MIXED_AT = (200, 300, 400)


def test_mixed_50_percent_flagged(coi, tmp_path):
    # The real R read itself has two second peaks in the region (595, 596):
    # below the 3-position threshold, so the unmodified clone stays "ok".
    baseline = verify_clone(_construct(coi), REGION, _reads()).mixed_positions
    assert len(baseline) < VerifyParams().mixed_min_positions
    f = _mixed_read(coi, tmp_path, 0.5, MIXED_AT)
    r = verify_clone(_construct(coi), REGION, _reads(f=f))
    assert r.verdict == "mixed"
    added = tuple(p for p in r.mixed_positions if p not in baseline)
    assert added == MIXED_AT and len(added) == 3


def test_mixed_50_percent_single_read_counts_exactly(coi, tmp_path):
    f = _mixed_read(coi, tmp_path, 0.5, MIXED_AT)
    r = verify_clone(_construct(coi), REGION, _reads(f=f, r=None))
    assert r.verdict == "mixed" and r.mixed_positions == MIXED_AT


def test_mixed_30_percent_not_flagged(coi, tmp_path):
    f = _mixed_read(coi, tmp_path, 0.3, MIXED_AT)
    r = verify_clone(_construct(coi), REGION, _reads(f=f))
    assert r.verdict == "ok"
    assert not set(MIXED_AT) & set(r.mixed_positions)


def test_no_data_for_foreign_read():
    puc19 = load_construct(PUC19_GB)
    raw = align_to_construct(real_read(F_NAME).seq, puc19)
    assert raw.identity < MIN_IDENTITY or raw.covered < clones.MIN_ALIGNED_BP
    r = verify_clone(puc19, None, [ReadRef(F_PATH, "c1", "F")])
    assert r.verdict == "no_data"
    res, = r.reads
    assert (res.error, res.alignment) == ("no_alignment", None)


def test_reverse_complemented_session_same_result(coi, tmp_path):
    construct = _construct(mutate(coi, P, "sub", "G"))
    plain = verify_clone(construct, REGION, _reads())
    paths = []
    for path, name in ((F_PATH, F_NAME), (R_PATH, R_NAME)):
        model = EditableRead(real_read(name))
        model.reverse_complement()
        paths.append(write_read(model.current(), tmp_path, name))
    flipped = verify_clone(construct, REGION, _reads(*paths))
    assert flipped.verdict == plain.verdict == "errors"
    assert flipped.differences == plain.differences
    assert flipped.coverage == plain.coverage   # strands are physical, not displayed
    assert flipped.mixed_positions == plain.mixed_positions
    assert [r.alignment.strand for r in flipped.reads] == ["-", "+"]


def test_region_across_origin(coi, tmp_path):
    # COI + a random spacer, so the circle is longer than a read; rotated so
    # the origin falls inside COI and the region wraps it. The reads are cut
    # to their COI part: their real flanks would align to the random spacer
    # as noise.
    paths = []
    for name in (F_NAME, R_NAME):
        model = EditableRead(real_read(name))
        model.apply_trim(40, 640)
        paths.append(write_read(model.current(), tmp_path, name))
    rng = random.Random(2)
    plasmid = mutate(coi, P, "sub", "G") + "".join(rng.choice("ACGT") for _ in range(400))
    L, k = len(plasmid), 300
    rotated = plasmid[k:] + plasmid[:k]
    region = ((REGION[0] - k) % L, (REGION[1] - k) % L)
    assert region[0] > region[1]
    linear = verify_clone(_construct(plasmid), REGION, _reads(*paths))
    circ = verify_clone(_construct(rotated, circular=True), region, _reads(*paths))
    assert all(x.alignment.wraps and not x.ambiguous for x in circ.reads)
    assert circ.verdict == linear.verdict == "errors"
    assert [replace(d, ref_pos=(d.ref_pos + k) % L) for d in circ.differences] == list(
        linear.differences)
    assert circ.region_differences[0].ref_pos == (P - k) % L
    assert circ.region_coverage_pct == linear.region_coverage_pct > 90.0
    assert circ.both_strands_pct == linear.both_strands_pct > 50.0
    assert {(p + k) % L: v for p, v in circ.coverage.items()} == linear.coverage


def test_region_positions():
    assert region_positions((2, 5), 10) == [2, 3, 4]
    assert region_positions((8, 2), 10) == [8, 9, 0, 1]
    assert region_positions(None, 4) == [0, 1, 2, 3]


def test_read_without_quality_is_incomplete_not_errors(coi, tmp_path):
    model = EditableRead(real_read(F_NAME))
    model.apply_trim(20, 680)   # drop the noisy ends: no quality, no trimming
    chrom = model.current()
    chrom = replace(chrom, quality=(0,) * len(chrom.seq))
    f = write_read(chrom, tmp_path, F_NAME)
    r = verify_clone(_construct(mutate(coi, P, "sub", "G")), REGION, _reads(f=f, r=None))
    assert r.verdict == "incomplete"
    assert r.region_coverage_pct == 0.0
    d, = r.region_differences
    assert (d.ref_pos, d.status, d.support, d.low_quality) == (P, "low_quality", 0, 1)


def test_reference_N_is_not_a_difference(coi):
    ref = mutate(coi, P, "sub", "N")
    ref = mutate(ref, P + 10, "sub", "R")   # IUPAC too (coi[260] is compatible or not, never counted)
    r = verify_clone(_construct(ref), REGION, _reads())
    assert r.verdict == "ok" and r.differences == ()


def test_unreadable_file_marked_not_fatal(coi, tmp_path):
    broken = tmp_path / "c1_X.ab1"
    broken.write_bytes(b"definitely not an ABIF file")
    reads = _reads() + [ReadRef(broken, "c1", "X")]
    r = verify_clone(_construct(coi), REGION, reads)
    assert r.verdict == "ok"
    bad = r.reads[-1]
    assert bad.error and bad.error != "no_alignment" and bad.alignment is None
    assert all(x.alignment is not None for x in r.reads[:2])


def test_ambiguous_read_excluded(coi, tmp_path):
    # reads trimmed well inside COI: no overhang can favour one of the copies
    paths = []
    for name in (F_NAME, R_NAME):
        model = EditableRead(real_read(name))
        model.apply_trim(100, 500)
        paths.append(write_read(model.current(), tmp_path, name))
    r = verify_clone(_construct(coi + coi), None, _reads(*paths))
    assert all(x.ambiguous and x.alignment.placements == 2 for x in r.reads)
    assert r.verdict == "no_data"


def test_genuine_repeat_single_read_is_no_data(coi, tmp_path):
    model = EditableRead(real_read(F_NAME))
    model.apply_trim(100, 500)
    f = write_read(model.current(), tmp_path, F_NAME)
    rng = random.Random(3)
    spacer = "".join(rng.choice("ACGT") for _ in range(300))
    r = verify_clone(_construct(coi + spacer + coi), None, _reads(f=f, r=None))
    res, = r.reads
    assert res.ambiguous and res.alignment.placements == 2
    assert r.verdict == "no_data"
    assert all(v == (0, 0) for v in r.coverage.values())


def test_homopolymer_deletion_is_an_error_not_no_data(coi):
    # co-optimal gap placements inside a run are one placement, not ambiguity (R15)
    r = verify_clone(_construct(mutate(coi, 200, "ins", "T")), REGION, _reads())  # TTTT -> TTTTT
    assert all(x.alignment.alternatives > 1 and not x.ambiguous for x in r.reads)
    assert r.verdict == "errors"
    assert [(d.kind, d.status) for d in r.region_differences] == [("deletion", "confirmed")]


@pytest.mark.parametrize("region, circular", [
    ((10, 10), True),        # empty
    ((0, 659), False),       # past the end
    ((-1, 10), False),
    ((658, 5), True),        # start outside the construct
    ((500, 100), False),     # wraps on a linear construct
])
def test_invalid_region_rejected(coi, monkeypatch, region, circular):
    def boom(*a, **k):
        raise AssertionError("aligned before validating the region")
    monkeypatch.setattr(clones, "align_to_construct", boom)
    monkeypatch.setattr(clones, "_load", boom)
    with pytest.raises(ValueError):
        verify_clone(_construct(coi, circular), region, _reads())
    with pytest.raises(ValueError):
        verify_all(_construct(coi, circular), region, Grouping(_reads()), VerifyParams())


def test_verify_all_progress_and_cancel(coi):
    grouping = Grouping([
        ReadRef(F_PATH, "c2", "F"), ReadRef(R_PATH, "c2", "R"),
        ReadRef(F_PATH, "c1", "F"),
    ])
    calls = []
    reports = verify_all(_construct(coi), REGION, grouping, VerifyParams(),
                         progress=lambda d, t: calls.append((d, t)))
    assert [r.clone for r in reports] == ["c1", "c2"]
    assert all(isinstance(r, CloneReport) for r in reports)
    assert [r.verdict for r in reports] == ["incomplete", "ok"]
    assert calls == [(1, 2), (2, 2)]

    calls.clear()
    with pytest.raises(CancelledError):
        verify_all(_construct(coi), REGION, grouping, VerifyParams(),
                   progress=lambda d, t: calls.append((d, t)),
                   cancelled=lambda: bool(calls))
    assert calls == [(1, 2)]


def test_annotate_sets_effect(coi):
    seen = []

    def annotate(construct, diff):
        seen.append((construct.name, diff.ref_pos))
        return "effect!"

    r = verify_clone(_construct(mutate(coi, P, "sub", "G")), REGION, _reads(), annotate=annotate)
    assert [d.effect for d in r.differences] == ["effect!"]
    assert seen == [("c", P)]


def test_load_read_applies_session(tmp_path):
    model = EditableRead(real_read(F_NAME))
    model.set_base(10, "A")
    model.apply_trim(5, 300)
    path = write_read(model.current(), tmp_path, F_NAME)
    chrom = load_read(path)
    assert chrom.seq == "".join(model.seq)
    assert chrom.quality == tuple(model.quality)


# Measured 2026-10-30 on the developer machine, after placements (R15) were
# added: 15.1-15.5 s (3 runs; the exact DP alignment against 10 kb dominates).
# LIMIT = 2 x 15.5, rounded up.
LIMIT = 32.0


def test_performance_12_clones_x_4_reads(coi):
    rng = random.Random(1)
    flank = "".join(rng.choice("ACGT") for _ in range(10_000 - len(coi)))
    construct = _construct(coi + flank)   # ~10 kb
    names = ["Achl_ACHLO006-09_1_F.ab1", "Achl_ACHLO006-09_2_R.ab1",
             "Achl_ACHLO007-09_1_F.ab1", "Achl_ACHLO007-09_2_R.ab1"]
    grouping = Grouping([
        ReadRef(DATA / "real" / n, f"clone{c}", n[-5]) for c in range(12) for n in names
    ])
    t = time.perf_counter()
    reports = verify_all(construct, None, grouping, VerifyParams())
    elapsed = time.perf_counter() - t
    print(f"verify_all 12x4 on {len(construct.seq)} bp: {elapsed:.2f} s")
    assert len(reports) == 12
    assert all(r.verdict != "no_data" for r in reports)
    assert not any(x.ambiguous or x.error for r in reports for x in r.reads)
    assert elapsed < LIMIT
