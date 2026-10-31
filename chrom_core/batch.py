"""Batch scanning and export of .ab1 folders."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from chrom_core import read_ab1, trim_window
from chrom_core.clones import load_read
from chrom_core.editing import EditableRead
from chrom_core.export import to_fasta, to_fastq, write_text

MIN_OK_LENGTH = 200    # trimmed length considered good
MIN_WARN_LENGTH = 100  # below this the read is "bad"
MIN_OK_MEAN_Q = 25.0   # mean quality required for "ok"


@dataclass(frozen=True)
class FileReport:
    path: Path
    length: int = 0
    mean_q: float | None = 0.0  # None: the read carries no quality values (all Q = 0)
    q20: int = 0
    trim_start: int = 0
    trim_end: int = 0
    error: str = ""

    @property
    def trimmed_length(self) -> int:
        return self.trim_end - self.trim_start

    @property
    def status(self) -> str:
        """One of 'ok', 'warn', 'bad', 'error'."""
        if self.error:
            return "error"
        if (
            self.trimmed_length >= MIN_OK_LENGTH
            and (self.mean_q or 0.0) >= MIN_OK_MEAN_Q
        ):
            return "ok"
        if self.trimmed_length >= MIN_WARN_LENGTH:
            return "warn"
        return "bad"


def scan_folder(folder: str | Path, window: int = 20, min_q: int = 20) -> list[FileReport]:
    """QC-scan every .ab1 in ``folder`` (non-recursive, sorted by name)."""
    folder = Path(folder)
    reports: list[FileReport] = []
    for path in sorted(p for p in folder.iterdir() if p.suffix.lower() == ".ab1"):
        try:
            chrom = read_ab1(path)
        except Exception as exc:
            reports.append(FileReport(path=path, error=str(exc)))
            continue
        start, end = trim_window(chrom.quality, window=window, min_q=min_q)
        reports.append(
            FileReport(
                path=path,
                length=len(chrom.seq),
                mean_q=(
                    round(sum(chrom.quality) / len(chrom.quality), 1)
                    if any(chrom.quality)
                    else None
                ),
                q20=sum(1 for q in chrom.quality if q >= 20),
                trim_start=start,
                trim_end=end,
            )
        )
    return reports


def batch_export(
    reports: list[FileReport],
    out_dir: str | Path,
    window: int = 20,
    min_q: int = 20,
    fmt: str = "fastq",
) -> list[tuple[Path, Path | None, str]]:
    """Trim and export every report; session edits are respected.

    Returns one ``(source, destination_or_None, error)`` tuple per report.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    results: list[tuple[Path, Path | None, str]] = []
    for report in reports:
        dest: Path | None = None
        error = ""
        try:
            model = EditableRead(load_read(report.path))
            start, end = trim_window(model.quality, window=window, min_q=min_q)
            if end - start >= 2:
                model.apply_trim(start, end)
            seq = "".join(model.seq)
            identifier = model.current().sample_id
            if fmt == "fasta":
                dest = out_dir / (report.path.stem + ".trimmed.fasta")
                write_text(dest, to_fasta(identifier, seq))
            else:
                dest = out_dir / (report.path.stem + ".trimmed.fastq")
                write_text(dest, to_fastq(identifier, seq, model.quality))
        except Exception as exc:
            error = str(exc)
        results.append((report.path, dest, error))
    return results


def csv_field(value, sep: str = ";") -> str:
    """One CSV field, quoted when it contains ``sep``, a quote or a line break."""
    text = str(value)
    if any(ch in text for ch in f'{sep}"\n\r'):
        text = '"' + text.replace('"', '""') + '"'
    return text


_csv_field = csv_field


def summary_csv(reports: list[FileReport]) -> str:
    """Semicolon-separated summary (opens cleanly in RU-locale Excel)."""
    lines = ["file;length;mean_q;q20;trimmed_length;status;error"]
    for report in reports:
        lines.append(
            ";".join(
                _csv_field(value)
                for value in (
                    report.path.name,
                    report.length,
                    "" if report.mean_q is None else report.mean_q,
                    report.q20,
                    report.trimmed_length,
                    report.status,
                    report.error,
                )
            )
        )
    return "\n".join(lines) + "\n"
