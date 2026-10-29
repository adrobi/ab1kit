"""Reference alignment for Sanger reads: semiglobal DP with affine gaps.

Built on ``Bio.Align.PairwiseAligner`` — exact dynamic programming, free end
gaps (the read is expected to sit inside the reference, ends may be ragged
after trimming), affine gap penalties tuned for single-base Sanger indels.
Both read orientations are tried; the best score wins.

Scoring notes: a naive *local* (classic Smith-Waterman) alignment silently
drops weak read ends — a substitution near the end would be reported as "no
data" instead of a mismatch. The free-end-gap global mode used here keeps
the whole read anchored. References beyond ``MAX_REFERENCE_LENGTH`` should
use a seeded aligner instead; we refuse rather than hang.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import TYPE_CHECKING

from Bio import Align

from chrom_core.seqops import IUPAC, reverse_complement

if TYPE_CHECKING:  # reference.py imports this module; avoid a runtime cycle
    from chrom_core.reference import Construct

MATCH_SCORE = 2.0
MISMATCH_SCORE = -3.0
GAP_OPEN = -10.0
GAP_EXTEND = -1.0
MAX_REFERENCE_LENGTH = 100_000
CIRCULAR_MARGIN = 50  # extra unrolled bases past len(read) for circular constructs


@dataclass(frozen=True)
class Difference:
    ref_pos: int            # 0-based; for insertions — the ref base AFTER the inserted bases
    kind: str               # 'substitution' | 'deletion' | 'insertion'
    ref_base: str           # '-' for insertions
    read_base: str          # '-' for deletions, inserted string for insertions
    read_pos: int | None    # 0-based index in the aligned read; None for deletions


@dataclass(frozen=True)
class ReadAlignment:
    strand: str                       # '+' or '-', relative to the read as passed in
    score: float
    ref_start: int                    # covered reference span, 0-based [start, end);
    ref_end: int                      # ref_end < ref_start: wraps the origin of a circle
    read_start: int
    read_end: int
    ref_to_read: dict[int, int]       # covered ref position -> read index
    differences: list[Difference]
    identity: float                   # matched / aligned columns
    alternatives: int                 # number of co-optimal alignments (1 = unique)
    oriented_read: str                # the read in its aligned orientation
    blocks: list                      # biopython aligned blocks, for gapped reconstruction
                                      # (unrolled coordinates for circular constructs)
    ref_length: int = 0               # construct length (0 = unknown, from align_read)

    @property
    def wraps(self) -> bool:
        return self.ref_end < self.ref_start

    @property
    def covered(self) -> int:
        if self.wraps:
            return self.ref_length - self.ref_start + self.ref_end
        return self.ref_end - self.ref_start

    def covers(self, pos: int) -> bool:
        if self.wraps:
            return self.ref_start <= pos < self.ref_length or 0 <= pos < self.ref_end
        return self.ref_start <= pos < self.ref_end

    def covered_positions(self) -> list[int]:
        if self.wraps:
            return list(range(self.ref_start, self.ref_length)) + list(range(self.ref_end))
        return list(range(self.ref_start, self.ref_end))


def _make_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = MATCH_SCORE
    aligner.mismatch_score = MISMATCH_SCORE
    # affine internal gaps; free end gaps on both sequences (the read sits
    # somewhere inside the reference, ends may be ragged after trimming).
    # Position-based attributes (internal/end) keep the aligner consistent —
    # mixing combined setters desyncs target/query scores in some versions.
    aligner.open_internal_gap_score = GAP_OPEN
    aligner.extend_internal_gap_score = GAP_EXTEND
    aligner.open_end_gap_score = 0.0
    aligner.extend_end_gap_score = 0.0
    return aligner


def parse_sequence(text: str) -> str:
    """Extract a single sequence from pasted text or FASTA (first record)."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if lines and lines[0].startswith(">"):
        lines = lines[1:]
        lines = lines[: next((i for i, l in enumerate(lines) if l.startswith(">")), len(lines))]
    seq = "".join("".join(lines).split()).upper()  # drop ALL whitespace incl. internal
    bad = set(seq) - set(IUPAC)
    if bad:
        raise ValueError(f"reference contains non-IUPAC characters: {sorted(bad)}")
    return seq


def _clean(read: str, reference: str) -> tuple[str, str]:
    read = "".join(read.split()).upper()
    reference = "".join(reference.split()).upper()
    if not read:
        raise ValueError("empty read")
    if not reference:
        raise ValueError("empty reference")
    return read, reference


def _check_length(length: int) -> None:
    if length > MAX_REFERENCE_LENGTH:
        raise ValueError(
            f"reference of {length} bp is too long for exact alignment "
            f"(limit {MAX_REFERENCE_LENGTH})"
        )


def align_read(read: str, reference: str) -> ReadAlignment:
    """Align a (possibly reverse-strand) read into the reference."""
    read, reference = _clean(read, reference)
    _check_length(len(reference))
    return _align(read, reference)


def _align(read: str, reference: str) -> ReadAlignment:
    """align_read without the size limit (callers check it on the real length)."""
    aligner = _make_aligner()
    best: tuple[float, object, str, str] | None = None
    for strand, oriented in (("+", read), ("-", reverse_complement(read))):
        alignments = aligner.align(reference, oriented)
        if best is None or alignments.score > best[0]:
            best = (alignments.score, alignments, strand, oriented)
    score, alignments, strand, oriented = best

    alignment = alignments[0]
    blocks = [(tuple(t), tuple(q)) for t, q in zip(alignment.aligned[0], alignment.aligned[1])]
    try:
        alternatives = len(alignments)
    except Exception:  # counting can blow up on pathological repeats
        alternatives = 0  # "unknown"

    ref_to_read, differences, identity = _walk_blocks(blocks, reference, oriented)
    differences.sort(key=lambda d: (d.ref_pos, d.kind))
    if blocks:
        ref_start, ref_end = int(blocks[0][0][0]), int(blocks[-1][0][1])
        read_start, read_end = int(blocks[0][1][0]), int(blocks[-1][1][1])
    else:  # no similarity at all: empty alignment (everything became free end gaps)
        ref_start = ref_end = read_start = read_end = 0

    return ReadAlignment(
        strand=strand,
        score=score,
        ref_start=ref_start,
        ref_end=ref_end,
        read_start=read_start,
        read_end=read_end,
        ref_to_read=ref_to_read,
        differences=differences,
        identity=identity,
        alternatives=alternatives,
        oriented_read=oriented,
        blocks=blocks,
    )


def _walk_blocks(
    blocks: list, reference: str, oriented: str
) -> tuple[dict[int, int], list[Difference], float]:
    """Per-position mapping, differences (unsorted) and identity of a block list."""
    ref_to_read: dict[int, int] = {}
    differences: list[Difference] = []
    matches = 0
    columns = 0
    prev_t: int | None = None
    prev_q = 0
    for (t0, t1), (q0, q1) in blocks:
        if prev_t is not None:
            ref_gap = t0 - prev_t
            read_gap = q0 - prev_q
            if ref_gap > 0 and read_gap == 0:
                for p in range(prev_t, t0):
                    differences.append(Difference(p, "deletion", reference[p], "-", None))
                    columns += 1
            elif read_gap > 0 and ref_gap == 0:
                differences.append(
                    Difference(int(t0), "insertion", "-", oriented[prev_q:q0], int(prev_q))
                )
                columns += read_gap
        for k in range(t1 - t0):
            ref_base = reference[t0 + k]
            read_base = oriented[q0 + k]
            ref_to_read[int(t0 + k)] = int(q0 + k)
            columns += 1
            if ref_base == read_base:
                matches += 1
            else:
                differences.append(
                    Difference(int(t0 + k), "substitution", ref_base, read_base, int(q0 + k))
                )
        prev_t, prev_q = t1, q1
    identity = matches / columns if columns else 0.0
    return ref_to_read, differences, identity


def align_to_construct(read: str, construct: Construct) -> ReadAlignment:
    """Align a read to a construct; circular constructs allow reads across the origin.

    Linear: ``align_read`` plus ``ref_length``. Circular: the read is aligned to
    the unrolled sequence ``seq + seq[:len(read) + CIRCULAR_MARGIN]`` and every
    ``ref_*`` coordinate, ``ref_to_read`` key and ``Difference.ref_pos`` is
    reported modulo ``L = len(seq)`` (``ref_end`` as ``((ref_end - 1) % L) + 1``),
    so ``ref_end < ref_start`` means the span wraps the origin. ``blocks`` stay
    in unrolled coordinates: pass ``construct.seq * 2`` to ``gapped_sequences``.

    A read longer than the circle can align over the same positions twice. The
    first pass (in alignment order, from the alignment start) is kept: blocks
    are clipped at ``start + L``, so at most L positions are reported, each
    once, and identity/differences are computed over the clipped alignment.
    Full coverage is normalised to ``(0, L)``. Differences are listed in
    alignment order (the part after the origin last). The size limit applies
    to ``L``, not to the unrolled sequence.
    """
    read, seq = _clean(read, construct.seq)
    L = len(seq)
    _check_length(L)
    if not construct.circular:
        return replace(_align(read, seq), ref_length=L)

    unrolled = seq + seq[: len(read) + CIRCULAR_MARGIN]
    raw = _align(read, unrolled)
    if not raw.blocks:
        return replace(raw, ref_length=L)

    shift = (raw.blocks[0][0][0] // L) * L  # make the alignment start in the first copy
    limit = raw.blocks[0][0][0] - shift + L
    clipped = []
    for (t0, t1), (q0, q1) in raw.blocks:
        t0, t1 = t0 - shift, t1 - shift
        if t0 >= limit:
            break
        if t1 > limit:
            q1 = q0 + (limit - t0)
            t1 = limit
        clipped.append(((int(t0), int(t1)), (int(q0), int(q1))))

    ref_to_read, differences, identity = _walk_blocks(clipped, unrolled, raw.oriented_read)
    differences.sort(key=lambda d: (d.ref_pos, d.kind))  # unrolled order = alignment order
    ref_start, ref_end = clipped[0][0][0], clipped[-1][0][1]
    if ref_end - ref_start >= L:
        ref_start, ref_end = 0, L
    else:
        ref_start, ref_end = ref_start % L, ((ref_end - 1) % L) + 1
    return replace(
        raw,
        ref_start=ref_start,
        ref_end=ref_end,
        read_start=clipped[0][1][0],
        read_end=clipped[-1][1][1],
        ref_to_read={p % L: q for p, q in ref_to_read.items()},
        differences=[replace(d, ref_pos=d.ref_pos % L) for d in differences],
        identity=identity,
        blocks=clipped,
        ref_length=L,
    )


def gapped_sequences(alignment: ReadAlignment, reference: str) -> tuple[str, str]:
    """Reconstruct the aligned pair over the covered span ('-' for gaps)."""
    reference = reference.upper()
    ref_aligned: list[str] = []
    read_aligned: list[str] = []
    if not alignment.blocks:
        return "", ""
    # start at the first block (== ref_start/read_start for align_read; for a
    # circular alignment the blocks are in unrolled coordinates)
    (prev_t, _), (prev_q, _) = alignment.blocks[0]
    for (t0, t1), (q0, q1) in alignment.blocks:
        ref_aligned.extend(reference[prev_t:t0])
        read_aligned.extend("-" * (t0 - prev_t))
        read_aligned.extend(alignment.oriented_read[prev_q:q0])
        ref_aligned.extend("-" * (q0 - prev_q))
        span = t1 - t0
        ref_aligned.extend(reference[t0:t1])
        read_aligned.extend(alignment.oriented_read[q0:q1])
        prev_t, prev_q = t1, q1
    return "".join(ref_aligned), "".join(read_aligned)
