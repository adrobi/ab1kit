from pathlib import Path

import pytest

from chrom_core import read_ab1
from chrom_core.abif import Chromatogram

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"
AB1_FILES = sorted(EXAMPLES.glob("*.ab1"))

pytestmark = pytest.mark.skipif(not AB1_FILES, reason="no example .ab1 files in examples/")


@pytest.mark.parametrize("path", AB1_FILES, ids=lambda p: p.name)
def test_read_ab1_invariants(path: Path) -> None:
    chrom = read_ab1(path)

    assert len(chrom.seq) > 100
    assert len(chrom.quality) == len(chrom.seq)
    assert len(chrom.peak_indices) == len(chrom.seq)
    assert set(chrom.traces) == set("ACGT")

    lengths = {len(t) for t in chrom.traces.values()}
    assert len(lengths) == 1, "dye channels must be equal length"

    assert set(chrom.seq.upper()) <= set("ACGTURYSWKMBDHVN-.")
    assert min(chrom.quality) >= 0
    assert max(chrom.quality) <= 99
    assert list(chrom.peak_indices) == sorted(chrom.peak_indices)

    all_values = [v for t in chrom.traces.values() for v in t]
    assert min(all_values) >= 0
    assert max(all_values) > 1000, "expected real signal, not flat noise"


@pytest.mark.parametrize("path", AB1_FILES, ids=lambda p: p.name)
def test_metadata_is_strings(path: Path) -> None:
    chrom = read_ab1(path)
    assert isinstance(chrom.metadata, dict)
    for value in chrom.metadata.values():
        assert isinstance(value, str) and value.strip()
        assert "\x00" not in value


def _chrom(**overrides) -> Chromatogram:
    defaults = dict(
        sample_id="x",
        seq="ACGT",
        quality=(40, 40, 40, 40),
        traces={"A": (1, 2, 3, 4), "C": (1, 2, 3, 4),
                "G": (1, 2, 3, 4), "T": (1, 2, 3, 4)},
        peak_indices=(0, 1, 2, 3),
    )
    defaults.update(overrides)
    return Chromatogram(**defaults)


def test_quality_length_mismatch_rejected() -> None:
    with pytest.raises(ValueError, match="quality length"):
        _chrom(quality=(40, 40))


def test_peak_count_mismatch_rejected() -> None:
    with pytest.raises(ValueError, match="peak count"):
        _chrom(peak_indices=(0, 1))


def test_bad_alphabet_rejected() -> None:
    with pytest.raises(ValueError, match="unexpected bases"):
        _chrom(seq="ACGX")


def test_quality_out_of_range_rejected() -> None:
    with pytest.raises(ValueError, match="quality"):
        _chrom(quality=(40, 40, 40, 200))


def test_peak_beyond_traces_rejected() -> None:
    with pytest.raises(ValueError, match="beyond the end"):
        _chrom(peak_indices=(0, 1, 2, 99))


def test_empty_traces_rejected() -> None:
    with pytest.raises(ValueError, match="dye channels"):
        _chrom(traces={})
