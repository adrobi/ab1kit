"""SNP calling: classify alignment differences using trace peak shapes.

A substitution difference between read and reference becomes one of:
- ``homozygous``  — clean single peak of a non-reference base;
- ``heterozygous``— double peak: a second channel rises to at least
  ``het_ratio`` of the primary height at the same position (classic ~50/50
  Sanger het), or the basecaller itself emitted an IUPAC mixed base;
- ``low_quality`` — phred below ``min_q``; reported but flagged untrusted.
"""

from __future__ import annotations

from dataclasses import dataclass

from chrom_core.abif import Chromatogram
from chrom_core.align import ReadAlignment
from chrom_core.seqops import IUPAC

DEFAULT_MIN_Q = 20
DEFAULT_HET_RATIO = 0.40
PEAK_WINDOW = 5  # trace points around the called peak scanned for heights

_IUPAC_SETS = {
    frozenset(bases): code for code, bases in IUPAC.items() if len(bases) > 1
}


def alleles_to_iupac(alleles) -> str:
    """Smallest IUPAC code covering exactly these concrete bases."""
    concrete = {a.upper() for a in alleles if a.upper() in "ACGT"}
    if len(concrete) <= 1:
        return next(iter(concrete), "N")
    return _IUPAC_SETS.get(frozenset(concrete), "N")


@dataclass(frozen=True)
class SnpCall:
    ref_pos: int          # 0-based position in the reference
    read_pos: int         # 0-based index in the read as displayed
    ref_base: str
    read_base: str        # called base at the position (may itself be IUPAC)
    kind: str             # 'homozygous' | 'heterozygous' | 'low_quality'
    quality: int
    ratio: float | None   # secondary/primary peak height; None if unmeasurable
    iupac: str            # code summarizing the observed alleles
    note: str = ""


def peak_heights(
    chrom: Chromatogram, read_index: int, window: int = PEAK_WINDOW
) -> dict[str, float]:
    """Max trace intensity per channel around the called peak."""
    peak = chrom.peak_indices[read_index]
    heights = {}
    for base, trace in chrom.traces.items():
        lo = max(0, peak - window)
        hi = min(len(trace), peak + window + 1)
        heights[base] = float(max(trace[lo:hi]))
    return heights


def call_snps(
    alignment: ReadAlignment,
    chrom: Chromatogram,
    min_q: int = DEFAULT_MIN_Q,
    het_ratio: float = DEFAULT_HET_RATIO,
) -> list[SnpCall]:
    """Classify every substitution in ``alignment`` using ``chrom`` traces.

    ``chrom`` must be the read in its displayed orientation; the alignment's
    oriented-read coordinates are mapped back automatically for '-' strand.
    """
    calls: list[SnpCall] = []
    n = len(chrom.seq)
    for diff in alignment.differences:
        if diff.kind != "substitution":
            continue
        read_idx = diff.read_pos if alignment.strand == "+" else n - 1 - diff.read_pos
        if not 0 <= read_idx < n:
            continue
        quality = chrom.quality[read_idx]
        called = chrom.seq[read_idx].upper()
        heights = peak_heights(chrom, read_idx)

        called_concrete = set(IUPAC.get(called, "ACGT"))
        primary = max(heights.get(b, 0.0) for b in called_concrete) or max(
            heights.values()
        )
        secondary_base, secondary = max(
            ((b, h) for b, h in heights.items()), key=lambda kv: kv[1]
        )
        # the strongest channel should be one of the called alleles; the
        # secondary candidate is the strongest channel outside them
        outside = [(b, h) for b, h in heights.items() if b not in called_concrete]
        if outside:
            secondary_base, secondary = max(outside, key=lambda kv: kv[1])
        ratio = secondary / primary if primary > 0 else None

        alleles = set(called_concrete)
        if ratio is not None and ratio >= het_ratio:
            alleles.add(secondary_base)
        iupac = alleles_to_iupac(alleles)

        if quality < min_q:
            kind, note = "low_quality", f"Q{quality} < {min_q}"
        elif len(called_concrete) > 1:
            kind, note = "heterozygous", "mixed base call"
        elif ratio is not None and ratio >= het_ratio:
            kind, note = "heterozygous", f"second peak {ratio:.0%} of primary"
        else:
            kind, note = "homozygous", ""

        calls.append(
            SnpCall(
                ref_pos=diff.ref_pos,
                read_pos=read_idx,
                ref_base=diff.ref_base,
                read_base=chrom.seq[read_idx],
                kind=kind,
                quality=quality,
                ratio=ratio,
                iupac=iupac,
                note=note,
            )
        )
    calls.sort(key=lambda c: c.ref_pos)
    return calls
