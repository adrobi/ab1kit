"""Shared fixtures for clone-verification tests: the real-data corpus in ``tests/data``."""

from __future__ import annotations

import os
import shutil
import struct
from dataclasses import replace
from functools import cache
from pathlib import Path

from Bio import SeqIO

from chrom_core import read_ab1, sessions
from chrom_core.abif import DEFAULT_DYE_ORDER, Chromatogram
from chrom_core.peaks import channel_baseline, peak_window
from chrom_core.seqops import IUPAC

DATA = Path(__file__).resolve().parent / "data"

REAL_READS: list[Path] = sorted((DATA / "real").glob("*.ab1"))

COI_GB = DATA / "refs" / "coi_HM417896.gb"

# Sample prefix -> COI accession of the closest-lineage reference (A. chlorotica
# is a cryptic aggregate; see tests/data/README.md). HM417896 is the primary one.
SAMPLE_REFERENCE: dict[str, str] = {
    "Achl_ACHLO006-09": "HM417896",
    "Achl_ACHLO007-09": "HM417896",
    "Achl_ACHLO040-09": "KC569719",
    "Achl_ACHLO041-09": "HM417923",
}
PUC19_GB = DATA / "refs" / "pUC19_L09137.gb"
SNAPGENE_DNA = DATA / "refs" / "sample-d.dna"


def perf_factor() -> float:
    """Multiplier for timing bounds (``AB1KIT_PERF_FACTOR``, default 1) on slow machines."""
    return float(os.environ.get("AB1KIT_PERF_FACTOR", "1"))


def real_read(name: str) -> Chromatogram:
    return read_ab1(DATA / "real" / name)


def coi_reference() -> str:
    return str(SeqIO.read(COI_GB, "genbank").seq).upper()


def reference_for(read_name: str) -> str:
    """COI sequence (upper case) of the lineage-matched reference for a real read."""
    for prefix, accession in SAMPLE_REFERENCE.items():
        if read_name.startswith(prefix):
            return str(SeqIO.read(DATA / "refs" / f"coi_{accession}.gb", "genbank").seq).upper()
    raise KeyError(read_name)


SPACING = 10  # trace points between called bases in synthetic chromatograms


def make_chrom(
    seq: str,
    quality: int = 40,
    heights: dict[tuple[int, str], float] | None = None,
    spacing: int = SPACING,
    primary: dict[str, float] | None = None,
) -> Chromatogram:
    """Synthetic chromatogram: clean primary peaks, optional extra heights.

    ``heights`` maps (base index, channel) -> trace value at that peak, e.g.
    {(10, "C"): 500} models a secondary C peak at read position 10.
    ``primary`` overrides the height of a called peak per channel (default
    1000), e.g. {"G": 500} models a systematically weak G dye.
    """
    heights = heights or {}
    primary = primary or {}
    n = len(seq)
    total = (n + 1) * spacing
    traces = {base: [80] * total for base in "ACGT"}
    peaks = []
    for i, called in enumerate(seq.upper()):
        peak = (i + 1) * spacing
        peaks.append(peak)
        # an IUPAC called base lights up all of its channels
        active = set(IUPAC[called]) if called in IUPAC else set("ACGT")
        for base in "ACGT":
            traces[base][peak] = primary.get(base, 1000) if base in active else 60
    for (index, base), value in heights.items():
        traces[base][(index + 1) * spacing] = value
    return Chromatogram(
        sample_id="synthetic",
        seq=seq,
        quality=(quality,) * n,
        traces={b: tuple(t) for b, t in traces.items()},
        peak_indices=tuple(peaks),
    )


# ---------------------------------------------------------------- semi-synthetic


def mutate(seq: str, pos: int, kind: str, base: str = "") -> str:
    """``seq`` with a substitution / insertion (before ``pos``) / deletion at ``pos``."""
    if kind == "sub":
        return seq[:pos] + base + seq[pos + 1:]
    if kind == "ins":
        return seq[:pos] + base + seq[pos:]
    if kind == "del":
        return seq[:pos] + seq[pos + 1:]
    raise ValueError(kind)


def with_secondary_peak(chrom: Chromatogram, index: int, base: str, ratio: float) -> Chromatogram:
    """Add a second peak in channel ``base`` under called base ``index``.

    The called peak's shape is copied into ``base`` (so it is a real local
    maximum, as ``peaks`` measures it), scaled so the *dye-normalised*
    secondary/primary ratio is ``ratio``.
    """
    called = chrom.seq[index].upper()
    baseline = channel_baseline(chrom, index)
    scale = ratio * baseline[base] / baseline[called]
    centre = chrom.peak_indices[index]
    window = peak_window(chrom, index)
    source = chrom.traces[called]
    trace = list(chrom.traces[base])
    for k in range(max(0, centre - window + 1), min(len(trace), centre + window)):
        trace[k] = max(trace[k], round(source[k] * scale))
    traces = dict(chrom.traces)
    traces[base] = tuple(trace)
    return replace(chrom, traces=traces)


@cache
def _source_by_sample_id() -> dict[str, Path]:
    return {read_ab1(p).sample_id: p for p in REAL_READS}

_DYE_OF = {base: tag for tag, base in DEFAULT_DYE_ORDER.items()}
_RC_CHANNEL = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _flip(traces: dict[str, tuple[int, ...]]) -> dict[str, tuple[int, ...]]:
    """Display traces of a reverse-complemented read <-> file traces."""
    return {b: tuple(reversed(traces[_RC_CHANNEL[b]])) for b in traces}


def _mismatches(a: dict, b: dict) -> int:
    return sum(x != y for base in a for x, y in zip(a[base], b[base]))


def _patch_traces(path: Path, traces: dict[str, tuple[int, ...]]) -> None:
    """Overwrite the DATA9-12 dye channels of an ABIF file in place (same lengths)."""
    data = bytearray(path.read_bytes())
    _, _, _, _, count, _, offset = struct.unpack(">4siHHiii", data[6:30])
    wanted = {_DYE_OF[b]: values for b, values in traces.items()}
    for e in range(count):
        entry = offset + 28 * e
        name, number, _, _, n, _, where = struct.unpack(">4siHHiii", data[entry:entry + 24])
        tag = f"{name.decode('latin-1')}{number}"
        if tag in wanted:
            assert n == len(wanted[tag]), "trace length must not change"
            data[where:where + 2 * n] = struct.pack(f">{n}h", *wanted[tag])
    path.write_bytes(bytes(data))


def write_read(chrom: Chromatogram, folder: Path, name: str) -> Path:
    """Save a modified real read as ``folder/name`` + editing session.

    The source ``.ab1`` (found by sample id) is copied; the called bases,
    quality and peak positions go to the session sidecar, so ``load_read``
    returns ``chrom``. A reverse-complemented display is detected from the
    traces; traces changed by ``with_secondary_peak`` are written into the copy.
    """
    source = _source_by_sample_id()[chrom.sample_id]
    path = Path(folder) / (name if name.lower().endswith(".ab1") else name + ".ab1")
    shutil.copyfile(source, path)
    original = read_ab1(source).traces
    flipped = _mismatches(_flip(chrom.traces), original) < _mismatches(chrom.traces, original)
    file_traces = _flip(chrom.traces) if flipped else dict(chrom.traces)
    if file_traces != original:
        _patch_traces(path, file_traces)
    sessions.save_session(sessions.session_path(path), {
        "source": path.name,
        "seq": chrom.seq,
        "quality": list(chrom.quality),
        "peaks": list(chrom.peak_indices),
        "edited": [],
        "reversed_complement": flipped,
    })
    return path
