"""Shared fixtures for clone-verification tests: the real-data corpus in ``tests/data``."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from chrom_core import read_ab1
from chrom_core.abif import Chromatogram

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
