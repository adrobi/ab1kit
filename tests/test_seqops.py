import pytest

from chrom_core.seqops import find_matches, iupac_compatible, reverse_complement


def test_plain_bases() -> None:
    assert reverse_complement("ATGC") == "GCAT"


def test_palindrome_is_invariant() -> None:
    assert reverse_complement("AAGCTT") == "AAGCTT"


def test_iupac_codes() -> None:
    assert reverse_complement("ACYGR") == "YCRGT"


def test_lowercase_and_n() -> None:
    assert reverse_complement("atgn") == "ncat"


def test_empty() -> None:
    assert reverse_complement("") == ""


# ---------------------------------------------------------------- find_matches
def test_exact_match() -> None:
    assert find_matches("GC", "ATGCGT") == [2]


def test_no_match() -> None:
    assert find_matches("AAA", "ATGCGT") == []


def test_overlapping_matches() -> None:
    assert find_matches("AA", "AAAA") == [0, 1, 2]


def test_iupac_pattern() -> None:
    # N matches anything
    assert find_matches("GN", "ATGC") == [2]
    # R = A/G: matches both AG (idx 1) and GG (idx 2)
    assert find_matches("RG", "AAGG") == [1, 2]


def test_reverse_match_maps_to_original_coordinates() -> None:
    # ATG is not present directly in seq; its rc (CAT) matches seq[6:9]
    assert find_matches("ATG", "TTAGCACAT", reverse=True) == [6]


def test_reverse_equals_direct_for_revcomp_palindrome() -> None:
    # AATT is its own reverse complement
    assert find_matches("AATT", "CAATTG") == [1]
    assert find_matches("AATT", "CAATTG", reverse=True) == [1]


def test_empty_pattern() -> None:
    assert find_matches("", "ATGC") == []


def test_bad_pattern_raises() -> None:
    with pytest.raises(ValueError):
        find_matches("AXT", "ATGC")


def test_iupac_compatible() -> None:
    assert iupac_compatible("R", "A")
    assert iupac_compatible("N", "T")
    assert not iupac_compatible("A", "G")
