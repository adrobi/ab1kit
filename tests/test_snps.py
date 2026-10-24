import random
from pathlib import Path

import pytest

from chrom_core import read_ab1
from chrom_core.abif import Chromatogram
from chrom_core.align import align_read
from chrom_core.seqops import IUPAC, reverse_complement
from chrom_core.snps import alleles_to_iupac, call_snps, peak_heights

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"

SPACING = 10  # trace points between called bases


def make_chrom(seq: str, quality: int = 40, heights: dict[tuple[int, str], float] | None = None) -> Chromatogram:
    """Synthetic chromatogram: clean primary peaks, optional extra heights.

    ``heights`` maps (base index, channel) -> trace value at that peak, e.g.
    {(10, "C"): 500} models a secondary C peak at read position 10.
    """
    heights = heights or {}
    n = len(seq)
    total = (n + 1) * SPACING
    traces = {base: [80] * total for base in "ACGT"}
    peaks = []
    for i, called in enumerate(seq.upper()):
        peak = (i + 1) * SPACING
        peaks.append(peak)
        # an IUPAC called base lights up all of its channels
        active = set(IUPAC[called]) if called in IUPAC else set("ACGT")
        for base in "ACGT":
            traces[base][peak] = 1000 if base in active else 60
    for (index, base), value in heights.items():
        traces[base][(index + 1) * SPACING] = value
    return Chromatogram(
        sample_id="synthetic",
        seq=seq,
        quality=(quality,) * n,
        traces={b: tuple(t) for b, t in traces.items()},
        peak_indices=tuple(peaks),
    )


def test_alleles_to_iupac() -> None:
    assert alleles_to_iupac({"A"}) == "A"
    assert alleles_to_iupac({"A", "G"}) == "R"
    assert alleles_to_iupac({"C", "T"}) == "Y"
    assert alleles_to_iupac({"A", "C", "G"}) == "V"
    assert alleles_to_iupac({"A", "C", "G", "T"}) == "N"


def test_peak_heights_reads_channels() -> None:
    chrom = make_chrom("ACGT")
    heights = peak_heights(chrom, 1)  # C at index 1
    assert heights["C"] == 1000
    assert heights["A"] < 200


def test_homozygous_snp() -> None:
    rng = random.Random(7)
    reference = "".join(rng.choice("ACGT") for _ in range(120))
    read = list(reference)
    read[50] = "A" if read[50] != "A" else "G"
    read = "".join(read)
    chrom = make_chrom(read)
    alignment = align_read(read, reference)
    calls = call_snps(alignment, chrom, min_q=20, het_ratio=0.4)
    assert len(calls) == 1
    call = calls[0]
    assert call.kind == "homozygous"
    assert call.ref_pos == 50 and call.read_pos == 50
    assert call.iupac == read[50]
    assert call.ratio is not None and call.ratio < 0.2


def test_heterozygous_double_peak() -> None:
    rng = random.Random(8)
    reference = "".join(rng.choice("ACGT") for _ in range(120))
    read = list(reference)
    original = read[50]
    read[50] = "A" if original != "A" else "G"  # called allele vs reference
    read = "".join(read)
    # second channel rises to 50% of the primary at the same position
    chrom = make_chrom(read, heights={(50, "C"): 500})
    alignment = align_read(read, reference)
    calls = call_snps(alignment, chrom, min_q=20, het_ratio=0.4)
    assert len(calls) == 1
    call = calls[0]
    assert call.kind == "heterozygous"
    assert call.iupac == alleles_to_iupac({read[50], "C"})
    assert call.ratio == pytest.approx(0.5, abs=0.05)


def test_small_secondary_peak_stays_homozygous() -> None:
    rng = random.Random(9)
    reference = "".join(rng.choice("ACGT") for _ in range(120))
    read = list(reference)
    read[50] = "A" if read[50] != "A" else "G"
    read = "".join(read)
    chrom = make_chrom(read, heights={(50, "C"): 150})  # 15% of primary
    alignment = align_read(read, reference)
    calls = call_snps(alignment, chrom)
    assert calls[0].kind == "homozygous"


def test_low_quality_is_flagged() -> None:
    rng = random.Random(10)
    reference = "".join(rng.choice("ACGT") for _ in range(120))
    read = list(reference)
    read[50] = "A" if read[50] != "A" else "G"
    read = "".join(read)
    chrom = make_chrom(read, quality=10)
    alignment = align_read(read, reference)
    calls = call_snps(alignment, chrom, min_q=20)
    assert calls[0].kind == "low_quality"
    assert "Q10" in calls[0].note


def test_basecaller_iupac_call_is_heterozygous() -> None:
    rng = random.Random(11)
    reference = "".join(rng.choice("ACGT") for _ in range(120))
    read = list(reference)
    read[50] = "R"  # basecaller already emitted a mixed base
    read = "".join(read)
    chrom = make_chrom(read)
    alignment = align_read(read, reference)
    calls = call_snps(alignment, chrom)
    assert calls[0].kind == "heterozygous"
    assert calls[0].iupac == "R"


def test_reverse_strand_coordinates_map() -> None:
    rng = random.Random(12)
    forward = "".join(rng.choice("ACGT") for _ in range(120))
    read = list(forward)
    read[50] = "A" if read[50] != "A" else "G"  # mutated read, forward display
    read = "".join(read)
    # the reference is given on the opposite strand
    reference = reverse_complement(forward)
    chrom = make_chrom(read)  # chromatogram in the displayed (forward) orientation
    alignment = align_read(read, reference)
    assert alignment.strand == "-"
    calls = call_snps(alignment, chrom)
    assert len(calls) == 1
    # oriented-read index maps back onto the displayed read correctly
    assert calls[0].read_pos == 50
    assert calls[0].read_base == read[50]
    assert calls[0].kind in ("homozygous", "heterozygous", "low_quality")


def test_no_differences_no_calls() -> None:
    chrom = make_chrom("ACGT" * 10)
    alignment = align_read("ACGT" * 10, "ACGT" * 10)
    assert call_snps(alignment, chrom) == []


@pytest.mark.skipif(not (EXAMPLES / "3730.ab1").exists(), reason="no example")
def test_real_read_clean_position_is_homozygous() -> None:
    chrom = read_ab1(EXAMPLES / "3730.ab1")
    reference = list(chrom.seq)
    # pick a clean high-quality position with a small secondary peak
    target = None
    for i in range(200, len(chrom.seq) - 200):
        if chrom.quality[i] >= 45 and i + 1 < len(chrom.seq):
            heights = peak_heights(chrom, i)
            called = chrom.seq[i]
            others = sorted(
                (h for b, h in heights.items() if b != called), reverse=True
            )
            if others and others[0] / heights[called] < 0.2:
                target = i
                break
    assert target is not None, "no clean position found in the example read"
    reference[target] = "A" if reference[target] != "A" else "G"
    alignment = align_read(chrom.seq, "".join(reference))
    calls = [c for c in call_snps(alignment, chrom) if c.ref_pos == target]
    assert len(calls) == 1
    assert calls[0].kind == "homozygous"
    assert calls[0].quality >= 45
