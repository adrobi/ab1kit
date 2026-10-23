from chrom_core.seqops import reverse_complement


def test_plain_bases():
    assert reverse_complement("ATGC") == "GCAT"


def test_palindrome_is_invariant():
    assert reverse_complement("AAGCTT") == "AAGCTT"


def test_iupac_codes():
    assert reverse_complement("ACYGR") == "YCRGT"


def test_lowercase_and_n():
    assert reverse_complement("atgn") == "ncat"


def test_empty():
    assert reverse_complement("") == ""
