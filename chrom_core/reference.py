"""Reference constructs (plasmid / insert): loading from FASTA, GenBank, SnapGene."""
from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path

from Bio import SeqIO

from .align import parse_sequence

_FORMAT_BY_EXT = {
    ".fa": "fasta", ".fasta": "fasta", ".fna": "fasta", ".txt": "fasta",
    ".gb": "genbank", ".gbk": "genbank", ".genbank": "genbank",
    ".dna": "snapgene",
}
_NAME_QUALIFIERS = ("label", "gene", "product", "note")


@dataclass(frozen=True)
class Feature:
    """Annotated region; 0-based, half-open. start > end: wraps the origin."""

    name: str
    kind: str
    start: int
    end: int
    strand: int
    codon_start: int = 1


@dataclass(frozen=True)
class Construct:
    name: str
    seq: str
    circular: bool
    features: tuple[Feature, ...] = ()

    def with_circular(self, flag: bool) -> Construct:
        return replace(self, circular=flag)


def construct_from_text(text: str, name: str = "pasted") -> Construct:
    return Construct(name, parse_sequence(text), True)


def _detect_format(path: Path) -> str:
    fmt = _FORMAT_BY_EXT.get(path.suffix.lower())
    if fmt:
        return fmt
    with open(path, "rb") as fh:
        head = fh.read(512).lstrip()
    if head.startswith(b">"):
        return "fasta"
    if head.startswith(b"LOCUS"):
        return "genbank"
    raise ValueError(f"cannot determine reference format of {path.name}")


def _feature(f) -> Feature:
    loc = f.location
    parts = list(loc.parts)
    if loc.strand == -1:
        parts.reverse()  # Biopython lists complement(join()) parts in transcript order
    start, end = int(parts[0].start), int(parts[-1].end)
    q = f.qualifiers
    name = next((q[k][0] for k in _NAME_QUALIFIERS if q.get(k)), f.type)
    try:
        codon_start = int(q["codon_start"][0]) if f.type == "CDS" else 1
    except (KeyError, IndexError, ValueError):
        codon_start = 1
    return Feature(str(name), f.type, start, end, int(loc.strand or 0), codon_start)


def load_construct(path: str | Path) -> Construct:
    path = Path(path)
    fmt = _detect_format(path)
    try:
        record = next(SeqIO.parse(path, fmt), None)
    except Exception as exc:  # Biopython raises assorted parser errors
        raise ValueError(f"cannot parse {path.name} as {fmt}: {exc}") from exc
    if record is None:
        raise ValueError(f"no sequence records in {path.name}")
    seq = parse_sequence(str(record.seq))
    if not seq:
        raise ValueError(f"empty sequence in {path.name}")
    if fmt == "fasta":
        circular, features = True, ()
    else:
        circular = record.annotations.get("topology") == "circular"
        features = tuple(_feature(f) for f in record.features)
    return Construct(record.id or path.stem, seq, circular, features)
