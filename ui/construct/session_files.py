"""File helpers of the construct workspace: read discovery, path identity, reports."""

from __future__ import annotations

import os
from collections.abc import Iterable
from pathlib import Path

from chrom_core.clone_report import clone_report_xlsx, clones_csv, differences_csv
from chrom_core.clones import CloneReport

PROJECT_SUFFIX = ".ab1proj"
CONSTRUCT_FILTER = "Constructs (*.fa *.fasta *.fna *.txt *.gb *.gbk *.genbank *.dna);;All files (*)"


def path_key(path: Path | str) -> str:
    """Comparable identity of a file path (absolute, case-normalised on Windows)."""
    return os.path.normcase(os.path.abspath(path))


def ab1_files(paths: Iterable[Path | str]) -> list[Path]:
    """``.ab1`` files among ``paths``; a folder contributes its own ``.ab1`` files."""
    out: list[Path] = []
    for p in map(Path, paths):
        if p.is_dir():
            out += sorted(x for x in p.iterdir() if x.suffix.lower() == ".ab1" and x.is_file())
        elif p.suffix.lower() == ".ab1":
            out.append(p)
    return out


def write_report(path: Path, reports: list[CloneReport], construct_name: str) -> None:
    """``.xlsx``: one workbook; otherwise ``path`` gets the per-clone CSV and
    ``<stem>_differences.csv`` next to it the per-difference CSV."""
    path = Path(path)
    if path.suffix.lower() == ".xlsx":
        clone_report_xlsx(path, reports, construct_name)
        return
    path.write_text(clones_csv(reports), encoding="utf-8")
    path.with_name(f"{path.stem}_differences.csv").write_text(
        differences_csv(reports), encoding="utf-8")
