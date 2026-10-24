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

from dataclasses import dataclass

from Bio import Align

from chrom_core.seqops import IUPAC, reverse_complement

MATCH_SCORE = 2.0
MISMATCH_SCORE = -3.0
GAP_OPEN = -10.0
GAP_EXTEND = -1.0
MAX_REFERENCE_LENGTH = 100_000


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
    ref_start: int                    # covered reference span, 0-based [start, end)
    ref_end: int
    read_start: int
    read_end: int
    ref_to_read: dict[int, int]       # covered ref position -> read index
    differences: list[Difference]
    identity: float                   # matched / aligned columns
    alternatives: int                 # number of co-optimal alignments (1 = unique)
    oriented_read: str                # the read in its aligned orientation
    blocks: list                      # biopython aligned blocks, for gapped reconstruction

    @property
    def covered(self) -> int:
        return self.ref_end - self.ref_start


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


def align_read(read: str, reference: str) -> ReadAlignment:
    """Align a (possibly reverse-strand) read into the reference."""
    read = "".join(read.split()).upper()
    reference = "".join(reference.split()).upper()
    if not read:
        raise ValueError("empty read")
    if not reference:
        raise ValueError("empty reference")
    if len(reference) > MAX_REFERENCE_LENGTH:
        raise ValueError(
            f"reference of {len(reference)} bp is too long for exact alignment "
            f"(limit {MAX_REFERENCE_LENGTH})"
        )

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

    differences.sort(key=lambda d: (d.ref_pos, d.kind))
    identity = matches / columns if columns else 0.0
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


def gapped_sequences(alignment: ReadAlignment, reference: str) -> tuple[str, str]:
    """Reconstruct the aligned pair over the covered span ('-' for gaps)."""
    reference = reference.upper()
    ref_aligned: list[str] = []
    read_aligned: list[str] = []
    prev_t = alignment.ref_start
    prev_q = alignment.read_start
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
