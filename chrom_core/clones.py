"""Clone verification: grouping reads into clones by file name."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

DEFAULT_PRIMER_TOKENS = (
    "F", "R", "FW", "FWD", "REV", "RV", "T7", "T7TERM", "SP6", "T3",
    "M13F", "M13R", "M13", "PUC-F", "PUC-R", "BGH", "CMV",
)

_SEP = r"[_\-. ]"
_TOKEN_RE = re.compile(r"[^_\-. ]+")
_WELL_RE = re.compile(r"[A-H]\d{1,2}")
_MAX_PRIMER_PARTS = 3  # "pUC-F" spans two tokens


def _norm(text: str) -> str:
    return re.sub(f"{_SEP}+", "-", text).upper()


@dataclass(frozen=True)
class ReadRef:
    path: Path
    clone: str
    primer: str


def _natural_key(text: str) -> list:
    return [
        (0, int(part), "") if part.isdigit() else (1, 0, part.casefold())
        for part in re.split(r"(\d+)", text)
        if part != ""
    ]


@dataclass
class Grouping:
    reads: list[ReadRef] = field(default_factory=list)
    unassigned: list[Path] = field(default_factory=list)

    def clones(self) -> dict[str, list[ReadRef]]:
        out: dict[str, list[ReadRef]] = {}
        for ref in self.reads:
            out.setdefault(ref.clone, []).append(ref)
        return {k: out[k] for k in sorted(out, key=_natural_key)}

    def assign(self, path: Path, clone: str, primer: str) -> None:
        """Put ``path`` into ``clone`` (moving it if already grouped or unassigned)."""
        path = Path(path)
        self.unassigned = [p for p in self.unassigned if p != path]
        ref = ReadRef(path, clone, primer)
        for i, old in enumerate(self.reads):
            if old.path == path:
                self.reads[i] = ref
                return
        self.reads.append(ref)

    def to_dict(self) -> dict:
        return {
            "reads": [
                {"path": str(r.path), "clone": r.clone, "primer": r.primer}
                for r in self.reads
            ],
            "unassigned": [str(p) for p in self.unassigned],
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Grouping":
        return cls(
            reads=[
                ReadRef(Path(r["path"]), r["clone"], r["primer"])
                for r in d.get("reads", [])
            ],
            unassigned=[Path(p) for p in d.get("unassigned", [])],
        )


def parse_name(
    stem: str, primer_tokens: Iterable[str] = DEFAULT_PRIMER_TOKENS
) -> tuple[str, str] | None:
    """Split a file stem into ``(clone, primer)``; ``None`` if no clone remains.

    The primer is the last dictionary token (else the last token). The clone is
    the text before it, minus trailing separators, a lone digit (primer order)
    and a leading plate well (``A05``) -- the well is stripped only when another
    clone token follows it, so ``B2_F`` is clone ``B2``. Anything after the
    primer (e.g. a repeat suffix ``_2``) is ignored.
    """
    tokens = list(_TOKEN_RE.finditer(stem))
    if not tokens:
        return None
    known = {_norm(t) for t in primer_tokens}
    best: tuple[int, int] | None = None  # (start, end) in ``stem``
    for i in range(len(tokens)):
        for j in range(i, min(i + _MAX_PRIMER_PARTS, len(tokens))):
            start, end = tokens[i].start(), tokens[j].end()
            if _norm(stem[start:end]) not in known:
                continue
            # the last match wins; on equal end the longer (earlier start) one
            if best is None or end > best[1] or (end == best[1] and start < best[0]):
                best = (start, end)
    if best is None:
        best = (tokens[-1].start(), tokens[-1].end())
    start, end = best
    primer = stem[start:end]
    clone = re.sub(f"{_SEP}+$", "", stem[:start])
    clone = re.sub(rf"(?:^|{_SEP})\d$", "", clone)          # primer order number
    clone = re.sub(f"{_SEP}+$", "", clone)
    # leading plate well, only if another clone token follows it
    clone = re.sub(rf"^{_WELL_RE.pattern}{_SEP}+(?=[^_\-. ])", "", clone)
    clone = re.sub(f"^{_SEP}+", "", clone)
    return (clone, primer) if clone else None


def group_reads(
    paths: Iterable[Path],
    pattern: str | None = None,
    primer_tokens: Iterable[str] = DEFAULT_PRIMER_TOKENS,
) -> Grouping:
    """Group reads into clones by file name (default rules or a custom regex)."""
    rx = None
    if pattern:
        try:
            rx = re.compile(pattern)
        except re.error as exc:
            raise ValueError(f"invalid grouping pattern: {exc}") from exc
        if not {"clone", "primer"} <= set(rx.groupindex):
            raise ValueError("grouping pattern needs named groups 'clone' and 'primer'")
    tokens = tuple(primer_tokens)
    grouping = Grouping()
    for p in map(Path, paths):
        if rx is not None:
            m = rx.fullmatch(p.stem)
            parsed = (m.group("clone") or "", m.group("primer") or "") if m else None
            if parsed is not None and not parsed[0]:
                parsed = None
        else:
            parsed = parse_name(p.stem, tokens)
        if parsed is None:
            grouping.unassigned.append(p)
        else:
            grouping.reads.append(ReadRef(p, parsed[0], parsed[1]))
    return grouping
