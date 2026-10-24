import shutil
from pathlib import Path

import pytest

from chrom_core import batch, read_ab1, sessions
from chrom_core.editing import EditableRead

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"
AB1 = sorted(EXAMPLES.glob("*.ab1"))

pytestmark = pytest.mark.skipif(not AB1, reason="no example .ab1 files")


@pytest.fixture
def folder(tmp_path: Path) -> Path:
    for sample in AB1:
        shutil.copy(sample, tmp_path / sample.name)
    (tmp_path / "broken.ab1").write_bytes(b"definitely not an ABIF file")
    return tmp_path


def test_scan_folder_qc(folder: Path) -> None:
    reports = batch.scan_folder(folder, window=20, min_q=20)
    by_name = {r.path.name: r for r in reports}

    assert set(by_name) == {"310.ab1", "3730.ab1", "broken.ab1"}

    good = by_name["3730.ab1"]
    assert good.length > 500
    assert good.mean_q > 20
    assert good.q20 > 100
    assert good.trimmed_length > 0
    assert good.status in ("ok", "warn")

    broken = by_name["broken.ab1"]
    assert broken.error
    assert broken.status == "error"


def test_scan_empty_and_missing(tmp_path: Path) -> None:
    empty = tmp_path / "empty"
    empty.mkdir()
    assert batch.scan_folder(empty) == []
    with pytest.raises(OSError):
        batch.scan_folder(tmp_path / "no-such-dir")


def test_batch_export_fastq(folder: Path, tmp_path: Path) -> None:
    reports = [r for r in batch.scan_folder(folder) if not r.error]
    out = tmp_path / "out"
    results = batch.batch_export(reports, out, window=20, min_q=20, fmt="fastq")

    assert len(results) == len(reports)
    for source, dest, error in results:
        assert not error
        assert dest is not None and dest.exists()
        lines = dest.read_text(encoding="utf-8").splitlines()
        assert lines[0].startswith("@")
        assert lines[2] == "+"
        assert len(lines[1]) == len(lines[3]) > 0
        assert dest.name.endswith(".trimmed.fastq")


def test_batch_export_fasta(folder: Path, tmp_path: Path) -> None:
    reports = [r for r in batch.scan_folder(folder) if not r.error][:1]
    results = batch.batch_export(reports, tmp_path / "fa", fmt="fasta")
    dest = results[0][1]
    assert dest and dest.exists()
    assert dest.read_text(encoding="utf-8").startswith(">")


def test_batch_export_respects_session(folder: Path, tmp_path: Path) -> None:
    sample = folder / "3730.ab1"
    edited_index = 400  # safely inside the trimmed region
    model = EditableRead(read_ab1(sample))
    model.set_base(edited_index, "N")
    sessions.save_session(
        sessions.session_path(sample), {"source": sample.name, **model.to_dict()}
    )
    reports = [r for r in batch.scan_folder(folder) if r.path == sample]
    results = batch.batch_export(reports, tmp_path / "out2")
    dest = results[0][1]
    assert dest is not None
    seq_line = dest.read_text(encoding="utf-8").splitlines()[1]
    offset = edited_index - reports[0].trim_start
    assert 0 <= offset < len(seq_line)
    assert seq_line[offset] == "N"


def test_summary_csv(folder: Path) -> None:
    reports = batch.scan_folder(folder)
    text = batch.summary_csv(reports)
    lines = text.strip().splitlines()
    assert lines[0] == "file;length;mean_q;q20;trimmed_length;status;error"
    assert len(lines) == 1 + len(reports)
    broken_row = next(line for line in lines if "broken" in line)
    assert broken_row.rstrip().endswith("error") or ";error" in broken_row


def test_csv_field_quoting() -> None:
    assert batch._csv_field("plain") == "plain"
    assert batch._csv_field('a"b;c') == '"a""b;c"'


def test_read_without_quality_reports_none_mean_q(folder: Path) -> None:
    by_name = {r.path.name: r for r in batch.scan_folder(folder)}
    assert by_name["310.ab1"].mean_q is None
    assert by_name["3730.ab1"].mean_q is not None
    row = next(
        line for line in batch.summary_csv(list(by_name.values())).splitlines()
        if line.startswith("310.ab1")
    )
    assert row.split(";")[2] == ""
