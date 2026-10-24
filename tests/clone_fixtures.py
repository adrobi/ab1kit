"""Shared fixtures for clone-verification tests: the real-data corpus in ``tests/data``."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from chrom_core import read_ab1
from chrom_core.abif import Chromatogram
from chrom_core.seqops import IUPAC

DATA = Path(__file__).resolve().parent / "data"

REAL_READS: list[Path] = sorted((DATA / "real").glob("*.ab1"))

COI_GB = DATA / "refs" / "coi_HM417896.gb"

# Sample prefix -> COI accession of the closest-lineage reference (A. chlorotica
# is a cryptic aggregate; see tests/data/README.md). HM417896 is the primary one.
SAMPLE_REFERENCE: dict[str, str] = {
    "Achl_ACHLO006-09": "HM417896",
    "Achl_ACHLO007-09": "HM417896",
    "Achl_ACHLO040-09": "KC569719",
    "Achl_ACHLO041-09": "HM417923",
}
PUC19_GB = DATA / "refs" / "pUC19_L09137.gb"
SNAPGENE_DNA = DATA / "refs" / "sample-d.dna"


def real_read(name: str) -> Chromatogram:
    return read_ab1(DATA / "real" / name)


def coi_reference() -> str:
    return str(SeqIO.read(COI_GB, "genbank").seq).upper()


def reference_for(read_name: str) -> str:
    """COI sequence (upper case) of the lineage-matched reference for a real read."""
    for prefix, accession in SAMPLE_REFERENCE.items():
        if read_name.startswith(prefix):
            return str(SeqIO.read(DATA / "refs" / f"coi_{accession}.gb", "genbank").seq).upper()
    raise KeyError(read_name)


SPACING = 10  # trace points between called bases in synthetic chromatograms


def make_chrom(
    seq: str,
    quality: int = 40,
    heights: dict[tuple[int, str], float] | None = None,
    spacing: int = SPACING,
    primary: dict[str, float] | None = None,
) -> Chromatogram:
    """Synthetic chromatogram: clean primary peaks, optional extra heights.

    ``heights`` maps (base index, channel) -> trace value at that peak, e.g.
    {(10, "C"): 500} models a secondary C peak at read position 10.
    ``primary`` overrides the height of a called peak per channel (default
    1000), e.g. {"G": 500} models a systematically weak G dye.
    """
    heights = heights or {}
    primary = primary or {}
    n = len(seq)
    total = (n + 1) * spacing
    traces = {base: [80] * total for base in "ACGT"}
    peaks = []
    for i, called in enumerate(seq.upper()):
        peak = (i + 1) * spacing
        peaks.append(peak)
        # an IUPAC called base lights up all of its channels
        active = set(IUPAC[called]) if called in IUPAC else set("ACGT")
        for base in "ACGT":
            traces[base][peak] = primary.get(base, 1000) if base in active else 60
    for (index, base), value in heights.items():
        traces[base][(index + 1) * spacing] = value
    return Chromatogram(
        sample_id="synthetic",
        seq=seq,
        quality=(quality,) * n,
        traces={b: tuple(t) for b, t in traces.items()},
        peak_indices=tuple(peaks),
    )
