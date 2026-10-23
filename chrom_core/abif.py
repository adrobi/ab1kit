"""Reading Sanger chromatograms (``.ab1``) via ``Bio.Sequencing.Abi``."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO

# DATA9..DATA12 carry the four dye channels. Which base a channel reports
# depends on the instrument's dye set; this default targets BigDye v3.x on
# ABI 31xx/35xx sequencers (the most common case) and can be overridden per
# file via ``dye_order``.
DEFAULT_DYE_ORDER = {"DATA9": "G", "DATA10": "A", "DATA11": "T", "DATA12": "C"}


@dataclass(frozen=True)
class Chromatogram:
    """A parsed Sanger read.

    Attributes:
        sample_id: sample/file identifier from the .ab1 metadata.
        seq: called bases (PBAS2), may include IUPAC ambiguity codes.
        quality: per-base phred quality (PCON2), aligned with ``seq``.
        traces: raw dye intensities keyed by base letter ('A'/'C'/'G'/'T').
        peak_indices: trace-space position of each called base (PLOC1).
    """

    sample_id: str
    seq: str
    quality: tuple[int, ...]
    traces: dict[str, tuple[int, ...]]
    peak_indices: tuple[int, ...]

    def __post_init__(self) -> None:
        n = len(self.seq)
        if len(self.quality) != n:
            raise ValueError(f"quality length {len(self.quality)} != seq length {n}")
        if len(self.peak_indices) != n:
            raise ValueError(f"peak count {len(self.peak_indices)} != seq length {n}")
        if not self.traces:
            raise ValueError("no dye channels found")


def read_ab1(path: str | Path, dye_order: dict[str, str] | None = None) -> Chromatogram:
    """Parse a ``.ab1`` file into a :class:`Chromatogram`."""
    with open(path, "rb") as fh:
        record = SeqIO.read(fh, "abi")
    raw = record.annotations.get("abif_raw", {})
    order = dye_order or DEFAULT_DYE_ORDER
    traces = {
        base: tuple(raw[channel])
        for channel, base in order.items()
        if channel in raw
    }
    if not traces:
        raise ValueError(f"{path}: no DATA9-12 dye channels in ABIF tags")
    return Chromatogram(
        sample_id=record.id,
        seq=str(record.seq),
        quality=tuple(record.letter_annotations.get("phred_quality", ())),
        traces=traces,
        peak_indices=tuple(raw.get("PLOC1", ())),
    )
