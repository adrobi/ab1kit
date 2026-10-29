"""Editable called sequence with undo/redo on top of a parsed chromatogram."""

from __future__ import annotations

from dataclasses import dataclass

from chrom_core import Chromatogram
from chrom_core.seqops import reverse_complement

IUPAC_LETTERS = "ACGTNRYSWKMBDHV"
INSERT_LETTERS = "ACGTN"


@dataclass(frozen=True)
class _Snapshot:
    seq: str
    quality: tuple[int, ...]
    peaks: tuple[int, ...]
    edited: dict[int, str | None]
    reversed_complement: bool


class EditableRead:
    """Mutable editing state over an immutable parsed :class:`Chromatogram`.

    Every mutation pushes the previous state onto the undo stack (capped at
    ``MAX_HISTORY``); ``edited`` maps edited positions to the original base
    (``None`` for inserted bases).
    """

    MAX_HISTORY = 100

    def __init__(self, chrom: Chromatogram) -> None:
        self.original = chrom
        self.seq: list[str] = list(chrom.seq)
        self.quality: list[int] = list(chrom.quality)
        self.peaks: list[int] = list(chrom.peak_indices)
        self.edited: dict[int, str | None] = {}
        self.reversed_complement = False
        self.modified = False
        self._undo: list[_Snapshot] = []
        self._redo: list[_Snapshot] = []

    # ---------------------------------------------------------------- history
    def _snapshot(self) -> _Snapshot:
        return _Snapshot(
            "".join(self.seq),
            tuple(self.quality),
            tuple(self.peaks),
            dict(self.edited),
            self.reversed_complement,
        )

    def _push(self) -> None:
        self._undo.append(self._snapshot())
        if len(self._undo) > self.MAX_HISTORY:
            self._undo.pop(0)
        self._redo.clear()
        self.modified = True

    def undo(self) -> bool:
        if not self._undo:
            return False
        current = self._snapshot()
        self._restore(self._undo.pop())
        self._redo.append(current)
        return True

    def redo(self) -> bool:
        if not self._redo:
            return False
        current = self._snapshot()
        self._restore(self._redo.pop())
        self._undo.append(current)
        return True

    def _restore(self, snap: _Snapshot) -> None:
        self.seq = list(snap.seq)
        self.quality = list(snap.quality)
        self.peaks = list(snap.peaks)
        self.edited = dict(snap.edited)
        self.reversed_complement = snap.reversed_complement

    @property
    def can_undo(self) -> bool:
        return bool(self._undo)

    @property
    def can_redo(self) -> bool:
        return bool(self._redo)

    # --------------------------------------------------------------- mutations
    def set_base(self, i: int, base: str) -> None:
        self._push()
        if i not in self.edited:
            self.edited[i] = self.seq[i]
        self.seq[i] = base

    def insert_base(self, i: int, base: str) -> None:
        """Insert ``base`` before index ``i`` (``len(seq)`` appends)."""
        self._push()
        self.seq.insert(i, base)
        self.quality.insert(i, self.quality[i - 1] if 0 < i <= len(self.quality) - 1
                            else (self.quality[i] if self.quality else 20))
        self.peaks.insert(i, self._interp_peak(i))
        self.edited = {k + (k >= i): v for k, v in self.edited.items()}
        self.edited[i] = None

    def delete_base(self, i: int) -> None:
        self._push()
        del self.seq[i]
        del self.quality[i]
        del self.peaks[i]
        self.edited.pop(i, None)
        self.edited = {k - (k > i): v for k, v in self.edited.items()}

    def reverse_complement(self) -> None:
        """Flip the whole read (as Chromas does) as one undoable operation."""
        self._push()
        n = len(self.seq)
        self.seq = list(reverse_complement("".join(self.seq)))
        self.quality.reverse()
        # current() reverses the traces over their full length: map peaks the same way
        last = len(next(iter(self.original.traces.values()))) - 1
        self.peaks = [last - p for p in reversed(self.peaks)]
        self.edited = {n - 1 - k: v for k, v in self.edited.items()}
        self.reversed_complement = not self.reversed_complement

    def apply_trim(self, start: int, end: int) -> None:
        """Keep bases ``[start, end)`` as one undoable operation."""
        start, end = max(0, start), min(len(self.seq), end)
        if end - start < 2:
            raise ValueError("trim would leave fewer than 2 bases")
        self._push()
        self.seq = self.seq[start:end]
        self.quality = self.quality[start:end]
        self.peaks = self.peaks[start:end]
        self.edited = {k - start: v for k, v in self.edited.items()
                       if start <= k < end}

    # ------------------------------------------------------------------ misc
    def _interp_peak(self, i: int) -> int:
        trace_len = len(next(iter(self.original.traces.values())))
        if not self.peaks:
            return 0
        if i == 0:
            return max(0, self.peaks[0] - 10)
        if i >= len(self.peaks):
            return min(trace_len - 1, self.peaks[-1] + 10)
        mid = (self.peaks[i - 1] + self.peaks[i]) // 2
        mid = max(mid, self.peaks[i - 1] + 1)
        mid = min(mid, self.peaks[i] - 1)
        return max(0, min(trace_len - 1, mid))

    def current(self) -> Chromatogram:
        traces = self.original.traces
        if self.reversed_complement:
            traces = {
                base: tuple(reversed(traces[comp]))
                for base, comp in (("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"))
            }
        return Chromatogram(
            self.original.sample_id,
            "".join(self.seq),
            tuple(self.quality),
            traces,
            tuple(self.peaks),
            dict(self.original.metadata),
        )

    # --------------------------------------------------------------- sessions
    def to_dict(self) -> dict:
        return {
            "seq": "".join(self.seq),
            "quality": list(self.quality),
            "peaks": list(self.peaks),
            "edited": [[i, orig] for i, orig in sorted(self.edited.items())],
            "reversed_complement": self.reversed_complement,
        }

    def restore_state(self, data: dict) -> None:
        seq = data["seq"]
        quality = list(data["quality"])
        peaks = list(data["peaks"])
        if not (len(seq) == len(quality) == len(peaks)):
            raise ValueError("session arrays have inconsistent lengths")
        if seq and (min(quality) < 0 or max(quality) > 99):
            raise ValueError("session quality out of range")
        self.seq = list(seq)
        self.quality = quality
        self.peaks = peaks
        self.edited = {int(i): (orig if orig else None) for i, orig in data["edited"]}
        self.reversed_complement = bool(data["reversed_complement"])
        self._undo.clear()
        self._redo.clear()
        self.modified = True
