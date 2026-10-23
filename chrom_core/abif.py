"""Reading Sanger chromatograms (``.ab1``) via ``Bio.Sequencing.Abi``."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO

# DATA9..DATA12 carry the four dye channels. Which base a channel reports
# depends on the instrument's dye set; this default targets BigDye v3.x on
# ABI 31xx/35xx sequencers (the most common case) and can be overridden per
# file via ``dye_order``.
DEFAULT_DYE_ORDER = {"DATA9": "G", "DATA10": "A", "DATA11": "T", "DATA12": "C"}

# ABIF tags surfaced as human-readable metadata.
META_TAGS = {
    "SMPL1": "sample",
    "MODL1": "instrument",
    "RUND1": "run_date",
    "RUNT1": "run_time",
    "CMNT1": "comment",
    "PRJT1": "project",
    "LANE1": "capillary",
    "PDMF1": "mobility",
}

# Allowed called-base alphabet: IUPAC codes plus gap markers some
# basecallers emit.
_ALLOWED_ALPHABET = set("ACGTURYSWKMBDHVN-.")


@dataclass(frozen=True)
class Chromatogram:
    """A parsed Sanger read.

    Attributes:
        sample_id: sample/file identifier from the .ab1 metadata.
        seq: called bases (PBAS2), may include IUPAC ambiguity codes.
        quality: per-base phred quality (PCON2), aligned with ``seq``.
        traces: raw dye intensities keyed by base letter ('A'/'C'/'G'/'T').
        peak_indices: trace-space position of each called base (PLOC1).
        metadata: selected ABIF tags (sample, instrument, dates, ...).
    """

    sample_id: str
    seq: str
    quality: tuple[int, ...]
    traces: dict[str, tuple[int, ...]]
    peak_indices: tuple[int, ...]
    metadata: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        n = len(self.seq)
        if len(self.quality) != n:
            raise ValueError(f"quality length {len(self.quality)} != seq length {n}")
        if len(self.peak_indices) != n:
            raise ValueError(f"peak count {len(self.peak_indices)} != seq length {n}")
        if not self.traces:
            raise ValueError("no dye channels found")
        bad_bases = set(self.seq.upper()) - _ALLOWED_ALPHABET
        if bad_bases:
            raise ValueError(f"unexpected bases in called sequence: {sorted(bad_bases)}")
        if self.quality and (min(self.quality) < 0 or max(self.quality) > 99):
            raise ValueError("phred quality outside 0..99")
        trace_len = len(next(iter(self.traces.values())))
        if self.peak_indices and max(self.peak_indices) >= trace_len:
            raise ValueError("peak position beyond the end of the traces")


def _tag_text(value) -> str:
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    return str(value).replace("\x00", "").strip()


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
    metadata = {
        key: _tag_text(raw[tag])
        for tag, key in META_TAGS.items()
        if tag in raw and _tag_text(raw[tag])
    }
    return Chromatogram(
        sample_id=record.id,
        seq=str(record.seq),
        quality=tuple(record.letter_annotations.get("phred_quality", ())),
        traces=traces,
        peak_indices=tuple(raw.get("PLOC1", ())),
        metadata=metadata,
    )
