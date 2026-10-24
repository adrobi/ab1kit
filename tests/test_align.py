from pathlib import Path

import pytest

from chrom_core import read_ab1
from chrom_core.align import Difference, align_read, gapped_sequences, parse_sequence
from chrom_core.seqops import reverse_complement

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


# ------------------------------------------------------------- synthetic ref
@pytest.fixture
def reference() -> str:
    # deterministic pseudo-random DNA
    import random

    rng = random.Random(42)
    return "".join(rng.choice("ACGT") for _ in range(400))


def test_perfect_read_inside_reference(reference: str) -> None:
    read = reference[10:60]
    a = align_read(read, reference)
    assert a.strand == "+"
    assert (a.ref_start, a.ref_end) == (10, 60)
    assert a.differences == []
    assert a.identity == 1.0
    assert a.alternatives == 1
    assert a.ref_to_read[10] == 0 and a.ref_to_read[59] == 49


def test_substitution_in_middle(reference: str) -> None:
    read = list(reference[10:60])
    read[25] = "N" if read[25] != "N" else "A"
    a = align_read("".join(read), reference)
    subs = [d for d in a.differences if d.kind == "substitution"]
    assert len(subs) == 1
    assert subs[0].ref_pos == 10 + 25
    assert subs[0].read_pos == 25


def test_substitution_at_read_end_not_dropped(reference: str) -> None:
    # the case that breaks naive local Smith-Waterman
    read = list(reference[10:60])
    read[-1] = "A" if read[-1] != "A" else "C"
    a = align_read("".join(read), reference)
    subs = [d for d in a.differences if d.kind == "substitution"]
    assert len(subs) == 1 and subs[0].ref_pos == 59


def test_substitution_at_read_start_not_dropped(reference: str) -> None:
    read = list(reference[10:60])
    read[0] = "A" if read[0] != "A" else "C"
    a = align_read("".join(read), reference)
    subs = [d for d in a.differences if d.kind == "substitution"]
    assert len(subs) == 1 and subs[0].ref_pos == 10


def test_single_base_deletion_in_read(reference: str) -> None:
    read = reference[10:30] + reference[31:60]  # base 30 deleted
    a = align_read(read, reference)
    dels = [d for d in a.differences if d.kind == "deletion"]
    assert len(dels) == 1 and dels[0].ref_pos == 30
    assert dels[0].read_base == "-" and dels[0].read_pos is None
    assert 30 not in a.ref_to_read
    assert a.ref_to_read[31] == 20


def test_single_base_insertion_in_read(reference: str) -> None:
    read = reference[10:30] + "G" + reference[30:60]
    a = align_read(read, reference)
    ins = [d for d in a.differences if d.kind == "insertion"]
    assert len(ins) == 1
    assert ins[0].read_base == "G" and ins[0].ref_base == "-"
    assert ins[0].ref_pos == 30  # anchored at the ref base after the insertion


def test_reverse_strand_is_detected(reference: str) -> None:
    read = reverse_complement(reference[100:160])
    a = align_read(read, reference)
    assert a.strand == "-"
    assert (a.ref_start, a.ref_end) == (100, 160)
    assert a.identity == 1.0


def test_free_end_gaps_on_both_sides(reference: str) -> None:
    flank = "TTTTTTTTTT" * 20
    read = reference[5:55]
    a = align_read(read, flank + read[:0] + reference + flank)  # noqa: E501
    assert (a.ref_start, a.ref_end) == (200 + 5, 200 + 55)


def test_no_similarity_gives_empty_span(reference: str) -> None:
    # 'A'*10 has no match on either strand of an all-C reference
    a = align_read("A" * 10, "C" * 100)
    assert a.covered == 0
    assert a.identity == 0.0


def test_ambiguous_repetitive_reference_flagged() -> None:
    a = align_read("ATATAT", "ATATATATATATATAT")
    assert a.alternatives > 1
    assert a.identity == 1.0


def test_empty_inputs_raise(reference: str) -> None:
    with pytest.raises(ValueError):
        align_read("", reference)
    with pytest.raises(ValueError):
        align_read("ACGT", "")


def test_huge_reference_is_refused() -> None:
    with pytest.raises(ValueError, match="too long"):
        align_read("ACGT", "A" * 100_001)


def test_identity_value(reference: str) -> None:
    read = list(reference[10:60])
    read[5] = "A" if read[5] != "A" else "C"
    a = align_read("".join(read), reference)
    assert a.identity == pytest.approx(49 / 50)


def test_gapped_sequences(reference: str) -> None:
    read = reference[10:30] + "G" + reference[30:60]  # one insertion
    a = align_read(read, reference)
    ref_aligned, read_aligned = gapped_sequences(a, reference)
    assert len(ref_aligned) == len(read_aligned) == 51
    assert read_aligned.count("-") == 0  # insertion is a gap in REFERENCE
    assert ref_aligned.count("-") == 1
    read_del = reference[10:30] + reference[31:60]  # one deletion
    a2 = align_read(read_del, reference)
    ref2, read2 = gapped_sequences(a2, reference)
    assert read2.count("-") == 1 and len(ref2) == len(read2)


# ---------------------------------------------------------------- real reads
@pytest.mark.skipif(
    not list(EXAMPLES.glob("*.ab1")), reason="no example .ab1 files"
)
def test_real_read_with_known_mutations() -> None:
    chrom = read_ab1(EXAMPLES / "3730.ab1")
    seq = chrom.seq[100:400]
    ref = list(seq)
    ref[50] = "A" if ref[50] != "A" else "G"
    ref[200] = "C" if ref[200] != "C" else "T"
    a = align_read(seq, "".join(ref))
    subs = [d for d in a.differences if d.kind == "substitution"]
    assert len(subs) == 2
    assert {s.ref_pos for s in subs} == {50, 200}
    assert a.identity > 0.99


# ------------------------------------------------------------------ parsing
def test_parse_plain_and_fasta() -> None:
    assert parse_sequence("acgt n") == "ACGTN"
    assert parse_sequence(">seq1 desc\nACGT\nTGCA\n>seq2\nTTTT") == "ACGTTGCA"


def test_parse_rejects_bad_chars() -> None:
    with pytest.raises(ValueError, match="non-IUPAC"):
        parse_sequence("ACGT123")
