"""Small sequence operations used across ab1kit."""

from __future__ import annotations

_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVBacgturyswkmbdhvb",
    "TGCAAYRSWMKVHDVBtgcaayrswmkvhdvb",
)

# IUPAC ambiguity: code -> the set of concrete bases it stands for.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
    "-": "-", ".": ".",
}


def reverse_complement(seq: str) -> str:
    """Reverse-complement a DNA string with IUPAC ambiguity codes."""
    return seq.translate(_COMPLEMENT)[::-1]


def iupac_compatible(a: str, b: str) -> bool:
    """True when IUPAC codes ``a`` and ``b`` share at least one concrete base."""
    return bool(set(IUPAC[a.upper()]) & set(IUPAC[b.upper()]))


def find_matches(pattern: str, seq: str, reverse: bool = False) -> list[int]:
    """Start indices (0-based, in original ``seq`` coordinates) of IUPAC-aware
    matches of ``pattern``.

    With ``reverse=True`` the pattern is matched against the reverse
    complement: a returned index ``i`` means ``seq[i:i+len(pattern)]``
    reverse-complements to the pattern.
    """
    if not pattern:
        return []
    probe = reverse_complement(pattern) if reverse else pattern
    bad = set(probe.upper()) - set(IUPAC)
    if bad:
        raise ValueError(f"pattern contains non-IUPAC characters: {sorted(bad)}")
    plen = len(probe)
    return [
        i
        for i in range(len(seq) - plen + 1)
        if all(iupac_compatible(seq[i + k], probe[k]) for k in range(plen))
    ]
