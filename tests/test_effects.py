"""Amino-acid effects of clone differences inside CDS features."""

from dataclasses import replace

from Bio.Seq import Seq

from chrom_core.clones import CloneDifference
from chrom_core.effects import difference_effect
from chrom_core.reference import Construct, Feature, load_construct
from chrom_core.seqops import reverse_complement
from tests.clone_fixtures import COI_GB

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def aa(codon: str) -> str:
    return str(Seq(codon).translate(table=11))


def sub(pos: int, ref: str, alt: str) -> CloneDifference:
    return CloneDifference(pos, "substitution", ref, alt, 2, 0, 0, "confirmed", True)


def coi():
    c = load_construct(COI_GB)
    cds = next(f for f in c.features if f.kind == "CDS")
    return c, cds


def codon_at(c, cds, n: int) -> tuple[int, str]:
    """Start position and codon of amino acid ``n`` (1-based) on the plus strand."""
    start = cds.start + cds.codon_start - 1 + (n - 1) * 3
    return start, c.seq[start:start + 3].upper()


def alt_changing(codon: str, within: int, synonymous: bool, avoid_stop=True) -> str:
    for b in "ACGT":
        if b == codon[within]:
            continue
        new = codon[:within] + b + codon[within + 1:]
        if (aa(new) == aa(codon)) == synonymous and not (avoid_stop and aa(new) == "*"):
            return b
    raise AssertionError("no such alternative")


def test_synonymous_third_codon_position():
    c, cds = coi()
    start, codon = codon_at(c, cds, 10)
    alt = alt_changing(codon, 2, synonymous=True)
    assert aa(codon[:2] + alt) == aa(codon)
    assert difference_effect(c, sub(start + 2, codon[2], alt)) == "synonymous"


def test_missense_label():
    c, cds = coi()
    n = 45
    start, codon = codon_at(c, cds, n)
    alt = alt_changing(codon, 0, synonymous=False)
    d = sub(start, codon[0], alt)
    aa_ref, aa_alt = aa(codon), aa(alt + codon[1:])
    assert difference_effect(c, d) == f"{aa_ref}{n}{aa_alt}"


def test_stop_gained():
    c, cds = coi()
    for n in range(2, 200):
        start, codon = codon_at(c, cds, n)
        for within in range(3):
            for b in "ACGT":
                new = codon[:within] + b + codon[within + 1:]
                if b != codon[within] and aa(new) == "*" and aa(codon) != "*":
                    got = difference_effect(c, sub(start + within, codon[within], b))
                    assert got == f"{aa(codon)}{n}*"
                    return
    raise AssertionError("no stop-gaining substitution found")


def test_frameshift_and_inframe():
    c, cds = coi()
    pos = cds.start + 30
    ins1 = CloneDifference(pos, "insertion", "-", "A", 2, 0, 0, "confirmed", True)
    ins3 = replace(ins1, read_base="ACG")
    del1 = CloneDifference(pos, "deletion", c.seq[pos], "-", 2, 0, 0, "confirmed", True)
    assert difference_effect(c, ins1) == "frameshift"
    assert difference_effect(c, ins3) == "inframe_ins"
    assert difference_effect(c, del1) == "frameshift"


def test_minus_strand_cds():
    """Mirrored COI on the minus strand (controller ruling R11).

    pUC19's GenBank has no CDS (only 'source'), so ``bla`` cannot be used: the
    COI construct is reverse-complemented and its CDS mirrored to strand -1.
    """
    c, cds = coi()
    L = len(c.seq)
    rc = Construct(
        "coi_rc", reverse_complement(c.seq.upper()), False,
        (replace(cds, start=L - cds.end, end=L - cds.start, strand=-1),),
    )
    for n in (1, 45, 100):
        start, codon = codon_at(c, cds, n)
        for within, synonymous in ((0, False), (2, True)):
            alt = alt_changing(codon, within, synonymous)
            plus = difference_effect(c, sub(start + within, codon[within], alt))
            mpos = L - 1 - (start + within)
            assert rc.seq[mpos] == COMP[codon[within]]
            minus = difference_effect(rc, sub(mpos, COMP[codon[within]], COMP[alt]))
            assert minus == plus


def test_wrapping_cds():
    c, cds = coi()
    # rotate so the CDS crosses the origin
    shift = 100
    seq = c.seq.upper()
    L = len(seq)
    rot = Construct(
        "rot", seq[shift:] + seq[:shift], True,
        (replace(cds, start=(cds.start - shift) % L, end=(cds.start + 600 - shift) % L),),
    )
    start, codon = codon_at(c, cds, 45)
    alt = alt_changing(codon, 0, synonymous=False)
    d = sub((start - shift) % L, codon[0], alt)
    assert difference_effect(rot, d) == f"{aa(codon)}45{aa(alt + codon[1:])}"


def test_outside_cds_empty():
    c, cds = coi()
    only = Construct("x", c.seq, False, (replace(cds, start=100, end=200),))
    assert difference_effect(only, sub(5, c.seq[5], "A" if c.seq[5] != "A" else "C")) == ""


def test_no_features_empty():
    c, _ = coi()
    plain = Construct("x", c.seq, False)
    assert difference_effect(plain, sub(10, "A", "C")) == ""
