"""Normalised secondary-peak analysis: the single source of double-peak calls.

Raw channel heights are not comparable between dyes: one channel can be
systematically weaker than the others along the whole read, which inflates
(or hides) a "second peak" ratio measured on raw intensities. Here every
channel height at a position is divided by that channel's local baseline —
the median height of nearby peaks where that channel was the called base —
and the ratio is taken between normalised heights.

The measurement window follows the local peak spacing (half the distance to
the nearest neighbouring peak) so a neighbour's peak never leaks in.
"""

from __future__ import annotations

from dataclasses import dataclass
from statistics import median

from chrom_core.abif import Chromatogram
from chrom_core.align import ReadAlignment
from chrom_core.seqops import IUPAC

BASELINE_SPAN = 20   # nearest called bases used for a channel's baseline
MIN_WINDOW = 2       # smallest half-window, trace points
_LONE_PEAK_WINDOW = 5  # half-window for a read with a single peak


@dataclass(frozen=True)
class SecondaryPeak:
    primary_base: str     # strongest (normalised) channel among the called alleles
    base: str             # strongest channel outside the called alleles; '' if none
    ratio: float | None   # normalised secondary / primary; None if unmeasurable


def peak_window(chrom: Chromatogram, i: int) -> int:
    """Half the distance to the nearest neighbouring peak, at least 2."""
    peaks = chrom.peak_indices
    gaps = []
    if i > 0:
        gaps.append(peaks[i] - peaks[i - 1])
    if i + 1 < len(peaks):
        gaps.append(peaks[i + 1] - peaks[i])
    if not gaps:
        return _LONE_PEAK_WINDOW
    return max(MIN_WINDOW, min(gaps) // 2)


def _heights(chrom: Chromatogram, i: int) -> dict[str, float]:
    """Per-channel peak height within ``peak_window`` of called peak ``i``.

    A channel's height is its tallest *local maximum* strictly inside the
    window, or its value at the called peak position if that is higher. A
    plain window maximum is not used: at half the peak spacing a neighbour's
    peak is still ~half its height, so its rising/falling shoulder at the
    window edge would masquerade as a second peak (on real reads this alone
    produced dozens of false double peaks per read).
    """
    peak = chrom.peak_indices[i]
    window = peak_window(chrom, i)
    heights = {}
    for base, trace in chrom.traces.items():
        if not 0 <= peak < len(trace):
            heights[base] = 0.0
            continue
        best = trace[peak]
        lo = max(0, peak - window)
        hi = min(len(trace) - 1, peak + window)
        for k in range(lo + 1, hi):
            if trace[k - 1] <= trace[k] >= trace[k + 1] and trace[k] > best:
                best = trace[k]
        heights[base] = float(best)
    return heights


def _called(chrom: Chromatogram, i: int) -> str:
    return chrom.seq[i].upper()


def _neighbours(n: int, i: int, span: int) -> list[int]:
    """The ``span`` indices nearest to ``i`` (excluding ``i``) within [0, n)."""
    out: list[int] = []
    d = 1
    while len(out) < span and (i - d >= 0 or i + d < n):
        if i - d >= 0:
            out.append(i - d)
        if len(out) < span and i + d < n:
            out.append(i + d)
        d += 1
    return out


def _baseline(chrom: Chromatogram, i: int, span: int, heights_of) -> dict[str, float]:
    per_channel: dict[str, list[float]] = {base: [] for base in chrom.traces}
    primaries: list[float] = []
    for j in _neighbours(len(chrom.seq), i, span):
        heights = heights_of(j)
        called = _called(chrom, j)
        if called in per_channel:
            per_channel[called].append(heights[called])
            primaries.append(heights[called])
        else:  # mixed / N call: its tallest channel stands in for the peak height
            primaries.append(max(heights.values(), default=0.0))
    fallback = median(primaries) if primaries else 0.0
    return {
        base: median(values) if values else fallback
        for base, values in per_channel.items()
    }


def channel_baseline(
    chrom: Chromatogram, i: int, span: int = BASELINE_SPAN
) -> dict[str, float]:
    """Per-channel baseline height around read position ``i``.

    For each channel: the median height of that channel over the ``span``
    nearest called bases whose called base is that channel. A channel never
    called nearby falls back to the median peak height of all ``span`` peaks.
    """
    return _baseline(chrom, i, span, lambda j: _heights(chrom, j))


def _secondary(heights: dict[str, float], baseline: dict[str, float], called: str) -> SecondaryPeak:
    norm = {
        base: h / baseline[base] if baseline.get(base, 0) > 0 else h
        for base, h in heights.items()
    }
    called_alleles = set(IUPAC.get(called, "ACGT")) & set(norm)
    if not called_alleles:
        called_alleles = set(norm)
    primary_base = max(sorted(called_alleles), key=lambda b: norm[b])
    primary = norm[primary_base]
    outside = [b for b in sorted(norm) if b not in called_alleles]
    if not outside:
        return SecondaryPeak(primary_base, "", None)
    base = max(outside, key=lambda b: norm[b])
    ratio = norm[base] / primary if primary > 0 else None
    return SecondaryPeak(primary_base, base, ratio)


def secondary_peak(chrom: Chromatogram, i: int) -> SecondaryPeak:
    """Strongest second channel under the called peak ``i``, normalised per dye."""
    return _secondary(_heights(chrom, i), channel_baseline(chrom, i), _called(chrom, i))


def scan_secondary_peaks(
    chrom: Chromatogram,
    alignment: ReadAlignment,
    min_ratio: float,
    min_q: int,
) -> list[tuple[int, SecondaryPeak]]:
    """Every covered reference position whose read peak has a second peak.

    ``chrom`` is the read in its displayed orientation (as for ``call_snps``);
    positions with quality below ``min_q`` are skipped. Sorted by ``ref_pos``.
    """
    n = len(chrom.seq)
    cache: dict[int, dict[str, float]] = {}

    def heights_of(j: int) -> dict[str, float]:
        if j not in cache:
            cache[j] = _heights(chrom, j)
        return cache[j]

    hits: list[tuple[int, SecondaryPeak]] = []
    for ref_pos in sorted(alignment.ref_to_read):
        idx = alignment.ref_to_read[ref_pos]
        if alignment.strand == "-":
            idx = n - 1 - idx
        if not 0 <= idx < n or chrom.quality[idx] < min_q:
            continue
        baseline = _baseline(chrom, idx, BASELINE_SPAN, heights_of)
        peak = _secondary(heights_of(idx), baseline, _called(chrom, idx))
        if peak.ratio is not None and peak.ratio >= min_ratio:
            hits.append((ref_pos, peak))
    return hits
