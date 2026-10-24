"""Shared fixtures for clone-verification tests: the real-data corpus in ``tests/data``."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from chrom_core import read_ab1
from chrom_core.abif import Chromatogram

DATA = Path(__file__).resolve().parent / "data"

REAL_READS: list[Path] = sorted((DATA / "real").glob("*.ab1"))

COI_GB = DATA / "refs" / "coi_HM417896.gb"
PUC19_GB = DATA / "refs" / "pUC19_L09137.gb"
SNAPGENE_DNA = DATA / "refs" / "sample-d.dna"


def real_read(name: str) -> Chromatogram:
    return read_ab1(DATA / "real" / name)


def coi_reference() -> str:
    return str(SeqIO.read(COI_GB, "genbank").seq).upper()
