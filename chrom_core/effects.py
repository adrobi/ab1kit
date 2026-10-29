"""Amino-acid consequences of clone differences inside CDS features."""

from __future__ import annotations

from typing import TYPE_CHECKING

from Bio.Seq import Seq


if TYPE_CHECKING:
    from chrom_core.clones import CloneDifference
    from chrom_core.reference import Construct, Feature

TRANSLATION_TABLE = 11
_BASES = frozenset("ACGT")
_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def _feature_positions(feature: Feature, length: int) -> list[int]:
    """Construct positions of ``feature`` in transcript order (wrap and strand aware)."""
    if feature.start <= feature.end:
        pos = list(range(feature.start, feature.end))
    else:
        pos = list(range(feature.start, length)) + list(range(0, feature.end))
    return pos[::-1] if feature.strand == -1 else pos


def _translate(codon: str) -> str:
    return str(Seq(codon).translate(table=TRANSLATION_TABLE))


def difference_effect(construct: Construct, diff: CloneDifference) -> str:
    """Effect of ``diff`` on the first CDS of ``construct`` that contains it.

    Substitution: ``"synonymous"``, ``"D45N"`` (one-letter codes, amino-acid
    number counted from the CDS start, frame from ``codon_start``) or
    ``"W12*"`` for a stop. Insertion/deletion: ``"frameshift"`` when the length
    is not a multiple of 3, else ``"inframe_ins"`` / ``"inframe_del"``. The
    length is that of this one ``CloneDifference``: a multi-base deletion
    arrives as several one-base differences, each of which is reported as
    ``"frameshift"``. ``""`` outside any CDS, with no features, or when the
    read base is not A/C/G/T or the codon is incomplete.
    """
    seq = construct.seq.upper()
    length = len(seq)
    if not length:
        return ""
    for feature in construct.features:
        if feature.kind != "CDS":
            continue
        pos = _feature_positions(feature, length)
        index = {p: i for i, p in enumerate(pos)}
        if diff.ref_pos not in index:
            continue
        if diff.kind == "insertion":
            if (diff.ref_pos - 1) % length not in index:
                continue
            return "frameshift" if len(diff.read_base) % 3 else "inframe_ins"
        if diff.kind == "deletion":
            return "frameshift" if len(diff.ref_base.replace("-", "")) % 3 else "inframe_del"
        if diff.kind != "substitution":
            return ""
        alt = diff.read_base.upper()
        if len(alt) != 1 or alt not in _BASES:
            return ""
        offset = index[diff.ref_pos] - (feature.codon_start - 1)
        if offset < 0:
            return ""
        number, within = divmod(offset, 3)
        first = number * 3 + feature.codon_start - 1
        codon_pos = pos[first:first + 3]
        if len(codon_pos) < 3:
            return ""
        if feature.strand == -1:   # positions are in transcript order: complement only
            alt = alt.translate(_COMPLEMENT)
            codon = "".join(seq[p].translate(_COMPLEMENT) for p in codon_pos)
        else:
            codon = "".join(seq[p] for p in codon_pos)
        mutated = codon[:within] + alt + codon[within + 1:]
        return _label(codon, mutated, number + 1)
    return ""


def _label(codon: str, mutated: str, number: int) -> str:
    ref_aa, alt_aa = _translate(codon), _translate(mutated)
    if ref_aa == alt_aa:
        return "synonymous"
    return f"{ref_aa}{number}{alt_aa}"
