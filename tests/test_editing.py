from pathlib import Path

import pytest

from chrom_core import read_ab1
from chrom_core.editing import EditableRead
from chrom_core.seqops import reverse_complement

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"
AB1 = sorted(EXAMPLES.glob("*.ab1"))

pytestmark = pytest.mark.skipif(not AB1, reason="no example .ab1 files")


@pytest.fixture
def model() -> EditableRead:
    return EditableRead(read_ab1(AB1[0]))


def test_set_base_and_undo(model: EditableRead) -> None:
    original = model.seq[5]
    model.set_base(5, "N")
    assert model.seq[5] == "N"
    assert model.edited[5] == original
    assert model.modified and model.can_undo

    assert model.undo()
    assert model.seq[5] == original
    assert model.can_redo
    assert model.redo()
    assert model.seq[5] == "N"


def test_insert_delete_roundtrip(model: EditableRead) -> None:
    n = len(model.seq)
    model.insert_base(3, "A")
    assert len(model.seq) == n + 1
    assert model.seq[3] == "A"
    assert model.edited[3] is None  # inserted
    # edited markers after the insert shift right
    model.set_base(10, "N")
    assert model.edited[10] is not None

    model.delete_base(3)
    assert len(model.seq) == n
    assert model.undo() and model.undo() and model.undo()
    assert "".join(model.seq) == model.original.seq


def test_insert_peaks_stay_monotonic(model: EditableRead) -> None:
    for i in (0, 1, len(model.seq) // 2, len(model.seq) - 1, len(model.seq)):
        model.insert_base(i, "C")
    assert list(model.peaks) == sorted(model.peaks)
    # Chromatogram rebuild must pass validation (peaks within traces)
    chrom = model.current()
    assert len(chrom.peak_indices) == len(chrom.seq)


def test_reverse_complement(model: EditableRead) -> None:
    seq = "".join(model.seq)
    model.reverse_complement()
    assert "".join(model.seq) == reverse_complement(seq)
    assert model.quality == list(reversed(model.original.quality))
    assert model.reversed_complement
    chrom = model.current()
    # flipped traces: A channel is original T channel reversed
    assert chrom.traces["A"] == tuple(reversed(model.original.traces["T"]))
    # peaks follow the flipped traces (reversed over the full trace length)
    last = len(model.original.traces["A"]) - 1
    assert chrom.peak_indices == tuple(last - p for p in reversed(model.original.peak_indices))
    # rc twice restores sequence, quality and peaks
    model.reverse_complement()
    assert "".join(model.seq) == seq
    assert model.quality == list(model.original.quality)
    assert model.peaks == list(model.original.peak_indices)
    assert not model.reversed_complement


def test_rc_shifts_edited_markers(model: EditableRead) -> None:
    model.set_base(2, "N")
    n = len(model.seq)
    model.reverse_complement()
    assert 2 not in model.edited
    assert model.edited[n - 1 - 2] == model.original.seq[2]


def test_apply_trim_and_edited_remap(model: EditableRead) -> None:
    model.set_base(5, "N")
    model.set_base(50, "N")
    model.apply_trim(10, 100)
    assert len(model.seq) == 90
    assert 5 not in model.edited
    assert model.edited[50 - 10] == model.original.seq[50]
    assert model.undo()
    assert len(model.seq) == len(model.original.seq)


def test_trim_too_short_raises(model: EditableRead) -> None:
    with pytest.raises(ValueError):
        model.apply_trim(0, 1)


def test_history_cap(model: EditableRead) -> None:
    for i in range(EditableRead.MAX_HISTORY + 20):
        model.set_base(0, "ACGT"[i % 4])
    assert len(model._undo) <= EditableRead.MAX_HISTORY


def test_session_state_roundtrip(model: EditableRead) -> None:
    model.set_base(1, "N")
    model.insert_base(2, "G")
    state = model.to_dict()
    other = EditableRead(read_ab1(AB1[0]))
    other.restore_state(state)
    assert "".join(other.seq) == "".join(model.seq)
    assert other.quality == model.quality
    assert other.peaks == model.peaks
    assert other.edited == model.edited
    assert other.reversed_complement == model.reversed_complement


def test_restore_rejects_inconsistent(model: EditableRead) -> None:
    state = model.to_dict()
    state["quality"] = state["quality"][:-5]
    with pytest.raises(ValueError):
        model.restore_state(state)
