"""Quality trimming of Sanger reads (sliding window, Chromas-style)."""

from __future__ import annotations

from collections.abc import Sequence


def trim_window(
    quality: Sequence[int], window: int = 20, min_q: int = 20
) -> tuple[int, int]:
    """Return the ``(start, end)`` span of good-quality bases.

    A window of ``window`` bases is slid along the read; a window is "good"
    when its mean quality is at least ``min_q``. The span starts at the first
    good window and ends at the last one. Returns ``(0, 0)`` when no window
    qualifies (including empty input).
    """
    n = len(quality)
    if n == 0 or window <= 0:
        return (0, 0)
    w = min(window, n)
    good = [
        sum(quality[i : i + w]) / w >= min_q
        for i in range(n - w + 1)
    ]
    if not any(good):
        return (0, 0)
    start = good.index(True)
    end = len(good) - 1 - good[::-1].index(True) + w
    return (start, end)
