"""Clone verification: grouping reads into clones, per-clone checks and verdicts."""

from __future__ import annotations

import re
from concurrent.futures import CancelledError
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Iterable

from chrom_core import sessions
from chrom_core.abif import Chromatogram, read_ab1
from chrom_core.align import Difference, ReadAlignment, align_to_construct
from chrom_core.editing import EditableRead
from chrom_core.effects import difference_effect
from chrom_core.peaks import scan_secondary_peaks
from chrom_core.reference import Construct
from chrom_core.seqops import iupac_compatible, reverse_complement
from chrom_core.trim import trim_window

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


# ------------------------------------------------------------------ verification

MIN_ALIGNED_BP = 50     # a read counts as aligned only past both thresholds
MIN_IDENTITY = 0.80
TRIM_WINDOW = 20        # quality trimming before alignment (sliding window)

Region = tuple[int, int]  # 0-based half-open; start > end wraps the origin (circular only)


@dataclass(frozen=True)
class VerifyParams:
    min_q: int = 20
    het_ratio: float = 0.40
    mixed_min_positions: int = 3


@dataclass(frozen=True)
class ReadResult:
    ref: ReadRef
    alignment: ReadAlignment | None   # None when ``error`` is set
    error: str = ""                   # "no_alignment" or why the file could not be used
    ambiguous: bool = False           # placements > 1: excluded from the verdict


@dataclass(frozen=True)
class CloneDifference:
    ref_pos: int
    kind: str           # 'substitution' | 'deletion' | 'insertion' (before ref_pos)
    ref_base: str
    read_base: str
    support: int        # reads showing it with Q >= min_q
    against: int        # reads showing the reference there with Q >= min_q
    low_quality: int    # reads showing it with Q < min_q
    status: str         # 'confirmed' | 'disputed' | 'low_quality'
    in_region: bool
    effect: str = ""


@dataclass(frozen=True)
class CloneReport:
    clone: str
    reads: tuple[ReadResult, ...]
    coverage: dict[int, tuple[int, int]]   # region position -> (reads on +, reads on -)
    differences: tuple[CloneDifference, ...]
    mixed_positions: tuple[int, ...]
    verdict: str        # 'no_data' | 'mixed' | 'errors' | 'incomplete' | 'ok'

    def _pct(self, test) -> float:
        if not self.coverage:
            return 0.0
        hits = sum(1 for v in self.coverage.values() if test(v))
        return 100.0 * hits / len(self.coverage)

    @property
    def region_coverage_pct(self) -> float:
        return self._pct(lambda v: v[0] + v[1] > 0)

    @property
    def both_strands_pct(self) -> float:
        return self._pct(lambda v: v[0] > 0 and v[1] > 0)

    @property
    def region_differences(self) -> tuple[CloneDifference, ...]:
        return tuple(d for d in self.differences if d.in_region)


def check_region(region: Region | None, length: int, circular: bool) -> Region:
    """``region`` normalised (``None`` -> whole construct); ``ValueError`` if invalid."""
    if region is None:
        return (0, length)
    start, end = region
    if not (0 <= start < length and 0 <= end <= length):
        raise ValueError(f"region {start}..{end} lies outside the construct (0..{length})")
    if start == end:
        raise ValueError("empty region")
    if start > end and not circular:
        raise ValueError("a region across the origin needs a circular construct")
    return (start, end)


def region_positions(region: Region | None, length: int) -> list[int]:
    """Construct positions of ``region`` in order (``None``: the whole construct)."""
    start, end = check_region(region, length, circular=True)
    if start < end:
        return list(range(start, end))
    return list(range(start, length)) + list(range(end))


def _load(path: Path) -> tuple[Chromatogram, bool]:
    """Current chromatogram (editing session applied) and its reversed flag."""
    path = Path(path)
    model = EditableRead(read_ab1(path))
    data = sessions.load_session(sessions.session_path(path))
    if data and data.get("source") == path.name:
        try:
            model.restore_state(data)
        except (KeyError, TypeError, ValueError):
            pass
    return model.current(), model.reversed_complement


def load_read(path: Path) -> Chromatogram:
    """Parse a read and apply its editing session sidecar when present."""
    return _load(path)[0]


def _align_trimmed(chrom: Chromatogram, construct: Construct, min_q: int) -> ReadAlignment:
    """Align the quality-trimmed read; indices refer to the full displayed read.

    Trimming keeps low-quality overhangs from wrapping around a circle. A read
    without a usable window (e.g. no quality values at all) is aligned whole.
    """
    n = len(chrom.seq)
    start, end = trim_window(chrom.quality, window=TRIM_WINDOW, min_q=min_q)
    if end - start < 2:
        start, end = 0, n
    al = align_to_construct(chrom.seq[start:end], construct)
    if (start, end) == (0, n):
        return al
    off = start if al.strand == "+" else n - end
    return replace(
        al,
        read_start=al.read_start + off,
        read_end=al.read_end + off,
        ref_to_read={p: i + off for p, i in al.ref_to_read.items()},
        differences=[
            d if d.read_pos is None else replace(d, read_pos=d.read_pos + off)
            for d in al.differences
        ],
        oriented_read=chrom.seq.upper() if al.strand == "+"
        else reverse_complement(chrom.seq.upper()),
        blocks=[(t, (q0 + off, q1 + off)) for t, (q0, q1) in al.blocks],
    )


class _Read:
    """An aligned read used for the verdict: quality lookups at reference sites."""

    def __init__(self, chrom: Chromatogram, al: ReadAlignment, strand: str, circular: bool):
        self.chrom, self.al, self.circular = chrom, al, circular
        self.strand = strand  # physical strand: a reverse-complemented display is undone
        self.map = al.ref_to_read

    def q(self, idx: int) -> int:
        n = len(self.chrom.seq)
        d = idx if self.al.strand == "+" else n - 1 - idx
        return self.chrom.quality[d] if 0 <= d < n else 0

    def _mapped_from(self, pos: int, step: int) -> int | None:
        """Read index of the nearest mapped position beyond ``pos`` towards ``step``."""
        p = pos
        for _ in range(self.al.ref_length):
            p += step
            if self.circular:
                p %= self.al.ref_length
            if not self.al.covers(p):
                return None
            if p in self.map:
                return self.map[p]
        return None

    def gap_q(self, pos: int) -> int:
        """Q of the read bases flanking a deleted position or an insertion point."""
        left = self._mapped_from(pos, -1)
        right = self.map[pos] if pos in self.map else self._mapped_from(pos, +1)
        if left is None or right is None:
            return 0
        return min(self.q(left), self.q(right))

    def site_q(self, pos: int) -> int | None:
        """Q at a reference position (flanks if deleted); None if not covered."""
        if not self.al.covers(pos):
            return None
        return self.q(self.map[pos]) if pos in self.map else self.gap_q(pos)

    def difference_q(self, d: Difference) -> int:
        return self.q(d.read_pos) if d.kind == "substitution" else self.gap_q(d.ref_pos)

    def shows_reference(self, d: Difference, ref: str, min_q: int) -> bool:
        m, pos = self.map, d.ref_pos
        if d.kind == "insertion":
            prev = (pos - 1) % len(ref)
            return (prev in m and pos in m and m[pos] == m[prev] + 1
                    and min(self.q(m[prev]), self.q(m[pos])) >= min_q)
        return (pos in m and self.al.oriented_read[m[pos]] == ref[pos]
                and self.q(m[pos]) >= min_q)


def _informative(d: Difference) -> bool:
    """N/IUPAC in the reference is never a difference; an ambiguous read call
    compatible with the reference base (e.g. N) is a no-call, not a difference."""
    if d.kind != "insertion" and d.ref_base not in "ACGT":
        return False
    if d.kind == "substitution" and d.read_base not in "ACGT":
        return not iupac_compatible(d.read_base, d.ref_base)
    return True


def _status(support: int, against: int) -> str:
    if support == 0:
        return "low_quality"
    return "disputed" if against else "confirmed"


def _verdict(used, mixed, differences, coverage, params: VerifyParams) -> str:
    if not used:
        return "no_data"
    if len(mixed) >= params.mixed_min_positions:
        return "mixed"
    if any(d.in_region and d.status != "low_quality" for d in differences):
        return "errors"
    if any(plus + minus == 0 for plus, minus in coverage.values()):
        return "incomplete"
    return "ok"


def verify_clone(
    construct: Construct,
    region: Region | None,
    reads: list[ReadRef],
    params: VerifyParams = VerifyParams(),
    annotate=difference_effect,
) -> CloneReport:
    """Check one colony's reads against ``construct`` within ``region``.

    ``annotate(construct, difference) -> str`` fills ``CloneDifference.effect``
    (default: ``difference_effect``; ``None`` leaves it empty).
    """
    L = len(construct.seq)
    positions = region_positions(check_region(region, L, construct.circular), L)
    in_region = set(positions)
    ref = construct.seq.upper()

    results: list[ReadResult] = []
    used: list[_Read] = []
    for read_ref in reads:
        try:
            chrom, flipped = _load(read_ref.path)
            al = _align_trimmed(chrom, construct, params.min_q)
        except Exception as exc:  # unreadable / unalignable file: mark it, go on
            results.append(ReadResult(read_ref, None, error=str(exc) or type(exc).__name__))
            continue
        if al.covered < MIN_ALIGNED_BP or al.identity < MIN_IDENTITY:
            results.append(ReadResult(read_ref, None, error="no_alignment"))
            continue
        ambiguous = al.placements > 1   # distinct placements, not co-optimal paths
        results.append(ReadResult(read_ref, al, ambiguous=ambiguous))
        if not ambiguous:
            strand = al.strand if not flipped else {"+": "-", "-": "+"}[al.strand]
            used.append(_Read(chrom, al, strand, construct.circular))

    coverage: dict[int, tuple[int, int]] = {}
    for pos in positions:
        good = [r.strand for r in used
                if r.al.covers(pos) and r.site_q(pos) >= params.min_q]
        coverage[pos] = (good.count("+"), good.count("-"))

    # merge per (position, kind, read base): reads showing it, with their Q
    found: dict[tuple, tuple[Difference, list[_Read], list[int]]] = {}
    for r in used:
        for d in r.al.differences:
            if _informative(d):
                entry = found.setdefault((d.ref_pos, d.kind, d.read_base), (d, [], []))
                entry[1].append(r)
                entry[2].append(r.difference_q(d))

    differences = []
    for key in sorted(found):
        d, showing, quals = found[key]
        support = sum(1 for q in quals if q >= params.min_q)
        against = sum(1 for r in used
                      if r not in showing and r.shows_reference(d, ref, params.min_q))
        inside = d.ref_pos in in_region or (
            d.kind == "insertion" and (d.ref_pos - 1) % L in in_region)
        cd = CloneDifference(
            d.ref_pos, d.kind, d.ref_base, d.read_base, support, against,
            len(quals) - support, _status(support, against), inside,
        )
        if annotate is not None:
            cd = replace(cd, effect=annotate(construct, cd) or "")
        differences.append(cd)

    mixed = tuple(sorted({
        pos
        for r in used
        for pos, _ in scan_secondary_peaks(r.chrom, r.al, params.het_ratio, params.min_q)
        if pos in in_region
    }))

    return CloneReport(
        clone=reads[0].clone if reads else "",
        reads=tuple(results),
        coverage=coverage,
        differences=tuple(differences),
        mixed_positions=mixed,
        verdict=_verdict(used, mixed, differences, coverage, params),
    )


def verify_all(
    construct: Construct,
    region: Region | None,
    grouping: Grouping,
    params: VerifyParams = VerifyParams(),
    progress=None,
    cancelled=None,
    annotate=difference_effect,
) -> list[CloneReport]:
    """``verify_clone`` for every clone of ``grouping``, in natural name order.

    ``progress(done, total)`` is called after each clone; ``cancelled()`` is
    polled before each one and, when true, raises
    ``concurrent.futures.CancelledError``.
    """
    check_region(region, len(construct.seq), construct.circular)
    groups = grouping.clones()
    reports = []
    for done, reads in enumerate(groups.values(), start=1):
        if cancelled is not None and cancelled():
            raise CancelledError()
        reports.append(verify_clone(construct, region, reads, params, annotate))
        if progress is not None:
            progress(done, len(groups))
    return reports
