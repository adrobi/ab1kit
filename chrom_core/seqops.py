"""Small sequence operations used across ab1kit."""

from __future__ import annotations

_COMPLEMENT = str.maketrans(
    "ACGTNRYSWKMDHVBacgtnryswkmdhvb",
    "TGCANYRSWMKHDVBtgcanyrswmkhdvb",
)


def reverse_complement(seq: str) -> str:
    """Reverse-complement a DNA string with IUPAC ambiguity codes."""
    return seq.translate(_COMPLEMENT)[::-1]
