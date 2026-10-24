"""FASTA / FASTQ export helpers."""

from __future__ import annotations

from pathlib import Path


def to_fasta(identifier: str, seq: str, width: int = 60) -> str:
    body = "\n".join(seq[i : i + width] for i in range(0, len(seq), width))
    return f">{identifier}\n{body}\n"


def to_fastq(identifier: str, seq: str, quality) -> str:
    if len(seq) != len(quality):
        raise ValueError("sequence and quality lengths differ")
    quals = "".join(chr(33 + min(93, q)) for q in quality)
    return f"@{identifier}\n{seq}\n+\n{quals}\n"


def write_text(path: str | Path, text: str) -> None:
    """Write LF line endings regardless of platform (bio-tool convention)."""
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)
