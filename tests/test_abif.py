from pathlib import Path

import pytest

from chrom_core import read_ab1

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

    assert set(chrom.seq.upper()) <= set("ACGTNRYSWKMDHVB")
    assert min(chrom.quality) >= 0
    assert max(chrom.quality) <= 99
    assert list(chrom.peak_indices) == sorted(chrom.peak_indices)

    all_values = [v for t in chrom.traces.values() for v in t]
    assert min(all_values) >= 0
    assert max(all_values) > 1000, "expected real signal, not flat noise"
