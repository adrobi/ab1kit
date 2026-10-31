import random

import pytest

from clone_fixtures import coi_reference, real_read
from chrom_core.align import (
    CIRCULAR_MARGIN,
    MAX_REFERENCE_LENGTH,
    align_read,
    align_to_construct,
    gapped_sequences,
)
from chrom_core.reference import Construct

READ_NAME = "Achl_ACHLO006-09_1_F.ab1"


@pytest.fixture(scope="module")
def read_seq() -> str:
    return real_read(READ_NAME).seq


@pytest.fixture(scope="module")
def ref() -> str:
    return coi_reference()


@pytest.fixture(scope="module")
def insert_read(read_seq, ref) -> str:
    """The part of the real read that lies inside the COI reference.

    The raw read (703 bp) overhangs BOTH ends of the 658 bp reference. Biopython's
    end gaps are free on one sequence per end only, so on an unrolled circle the
    44-base junk tail is cheaper aligned onto the start of the circle than gapped
    — the read then legitimately "covers" the whole ring. A clone read sits
    inside its construct, so the rotation tests use the inside part.
    """
    full = align_read(read_seq, ref)
    return read_seq[full.read_start:full.read_end]


@pytest.fixture(scope="module")
def aln0(insert_read, ref):
    return align_read(insert_read, ref)


@pytest.fixture(scope="module")
def mid(aln0) -> int:
    m = (aln0.ref_start + aln0.ref_end) // 2
    # the new origin must fall well inside the read's covered span
    assert aln0.ref_start + 100 < m < aln0.ref_end - 100
    return m


@pytest.fixture(scope="module")
def rotated(ref, mid) -> str:
    return ref[mid:] + ref[:mid]


def _random_seq(n: int, seed: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_margin_constant() -> None:
    assert CIRCULAR_MARGIN == 50


def test_linear_matches_align_read(read_seq, ref) -> None:
    aln0 = align_read(read_seq, ref)
    aln = align_to_construct(read_seq, Construct("coi", ref, False))
    assert aln.differences == aln0.differences
    assert aln.identity == aln0.identity
    assert (aln.ref_start, aln.ref_end) == (aln0.ref_start, aln0.ref_end)
    assert aln.ref_to_read == aln0.ref_to_read
    assert aln.ref_length == len(ref)
    assert aln.covered == aln0.covered
    assert aln.covered_positions() == list(range(aln0.ref_start, aln0.ref_end))
    assert aln0.ref_start > 0
    assert aln.covers(aln0.ref_start) and not aln.covers(aln0.ref_start - 1)


def test_read_across_origin_of_rotated_reference(insert_read, ref, aln0, mid, rotated) -> None:
    aln = align_to_construct(insert_read, Construct("rot", rotated, True))
    L = len(ref)
    assert aln.ref_length == L
    assert aln.ref_end < aln.ref_start  # wrapped
    assert len(aln.covered_positions()) == aln0.covered
    assert aln.covered == aln0.covered
    assert [(d.ref_pos - (L - mid)) % L for d in aln.differences] == [
        d.ref_pos for d in aln0.differences
    ]
    assert aln.covers(0) and aln.covers(L - 1)
    assert aln.identity == pytest.approx(aln0.identity)
    # every key is a real construct position; the mapping is the rotation of aln0's
    assert all(0 <= p < L for p in aln.ref_to_read)
    assert {(p + mid) % L: q for p, q in aln.ref_to_read.items()} == aln0.ref_to_read
    # gapped reconstruction over the unrolled reference gives the same pair
    assert gapped_sequences(aln, rotated * 2) == gapped_sequences(aln0, ref)


def test_same_rotation_linear_breaks_read(insert_read, aln0, rotated) -> None:
    aln = align_to_construct(insert_read, Construct("rot", rotated, False))
    assert aln.covered < aln0.covered  # without the ring the read is cut at the origin
    # only one side of the origin can be aligned: roughly half of the read's span
    assert aln.covered <= aln0.covered * 0.6
    assert aln.ref_end >= aln.ref_start


def test_real_read_longer_than_circular_construct(read_seq, ref) -> None:
    assert len(read_seq) > len(ref)
    aln = align_to_construct(read_seq, Construct("coi", ref, True))
    positions = aln.covered_positions()
    assert len(positions) == len(set(positions)) <= len(ref)
    assert all(0 <= p < len(ref) for p in positions)
    assert set(aln.ref_to_read) <= set(positions)
    assert all(0 <= d.ref_pos < len(ref) for d in aln.differences)
    assert aln.covered == len(positions)


def test_tandem_read_longer_than_circle_reports_each_position_once() -> None:
    circle = _random_seq(300, seed=7)
    read = (circle * 3)[50:50 + 400]  # wraps the circle more than once
    aln = align_to_construct(read, Construct("small", circle, True))
    L = len(circle)
    positions = aln.covered_positions()
    assert len(positions) == L and sorted(positions) == list(range(L))
    assert aln.covered == L
    assert (aln.ref_start, aln.ref_end) == (0, L)
    assert all(aln.covers(p) for p in range(L))
    assert sorted(aln.ref_to_read) == list(range(L))
    # the first pass over the circle wins: read index 0 sits at construct position 50
    assert aln.ref_to_read[50] == 0
    assert aln.ref_to_read[49] == L - 1
    assert aln.differences == []
    assert aln.identity == 1.0
    ref_g, read_g = gapped_sequences(aln, circle * 2)
    assert len(ref_g) == L and ref_g == read_g == read[:L]


def test_wrapped_coordinates_minus_strand() -> None:
    from chrom_core.seqops import reverse_complement

    circle = _random_seq(500, seed=11)
    insert = circle[420:] + circle[:90]  # 170 bp across the origin
    aln = align_to_construct(reverse_complement(insert), Construct("c", circle, True))
    assert aln.strand == "-"
    assert (aln.ref_start, aln.ref_end) == (420, 90)
    assert aln.covered == 170
    assert aln.covered_positions() == list(range(420, 500)) + list(range(90))
    assert aln.ref_to_read[420] == 0 and aln.ref_to_read[0] == 80 and aln.ref_to_read[89] == 169
    assert not aln.covers(90) and not aln.covers(419)


def test_wrapped_differences_are_mod_length() -> None:
    circle = _random_seq(400, seed=3)
    insert = list(circle[350:] + circle[:60])  # covers 350..399, 0..59
    insert[55] = "A" if insert[55] != "A" else "C"  # construct position 5
    aln = align_to_construct("".join(insert), Construct("c", circle, True))
    assert [(d.ref_pos, d.kind) for d in aln.differences] == [(5, "substitution")]


def test_limit_applies_to_original_length() -> None:
    big = _random_seq(MAX_REFERENCE_LENGTH, seed=5)
    read = big[1000:1030]
    aln = align_to_construct(read, Construct("big", big, True))
    assert aln.ref_length == MAX_REFERENCE_LENGTH
    with pytest.raises(ValueError):
        align_to_construct(read, Construct("big", big + "A", True))
    with pytest.raises(ValueError):
        align_to_construct(read, Construct("big", big + "A", False))


def test_read_in_unrolled_prefix_is_not_ambiguous() -> None:
    circle = _random_seq(1000, seed=21)
    read = circle[10:210]  # also present in the unrolled copy of the start
    aln = align_to_construct(read, Construct("c", circle, True))
    assert (aln.ref_start, aln.ref_end) == (10, 210)
    assert aln.alternatives == 1
    assert aln.placements == 1
    linear = align_to_construct(read, Construct("c", circle, False))
    assert (linear.alternatives, linear.placements) == (1, 1)


def test_genuine_repeat_in_circle_stays_ambiguous() -> None:
    circle = list(_random_seq(1000, seed=23))
    repeat = _random_seq(100, seed=24)
    circle[20:120] = repeat   # inside the unrolled prefix as well
    circle[600:700] = repeat
    circle = "".join(circle)
    aln = align_to_construct(repeat, Construct("c", circle, True))
    assert aln.alternatives == 2
    assert aln.placements == 2


def test_homopolymer_indel_across_origin_is_one_placement() -> None:
    circle = _random_seq(1000, seed=31)
    circle = circle[:995] + "GTTTC"           # a TTT run just before the origin
    read = circle[900:996] + "TTTT" + circle[999:] + circle[:100]   # one extra T
    aln = align_to_construct(read, Construct("c", circle, True))
    assert aln.wraps
    assert aln.alternatives > 1
    assert aln.placements == 1
