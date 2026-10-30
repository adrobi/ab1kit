"""Reference-comparison dock: table of differences with SNP verdicts."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import QDockWidget, QTableWidget, QTableWidgetItem, QWidget

from chrom_core.align import ReadAlignment
from chrom_core.seqops import reverse_complement
from i18n import tr


def reference_overlay(
    alignment: ReadAlignment, reference_seq: str, read_length: int
) -> tuple[str, dict[int, int]]:
    """Reference letters and ref -> index mapping in the *displayed* read's
    orientation (the alignment indexes the reverse complement on '-')."""
    if alignment.strand == "+":
        return reference_seq, alignment.ref_to_read
    n = read_length
    return (
        reverse_complement(reference_seq)[::-1],  # complement, same ref indexing
        {pos: n - 1 - idx for pos, idx in alignment.ref_to_read.items()},
    )


def displayed_diff_index(alignment: ReadAlignment, diff, read_length: int) -> int | None:
    """0-based index of ``diff`` in the displayed read (None for deletions)."""
    if diff.read_pos is None:
        return None
    if alignment.strand == "+":
        return diff.read_pos
    n = read_length
    span = len(diff.read_base) if diff.kind == "insertion" else 1
    return n - diff.read_pos - span


class DiffDock(QDockWidget):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(tr("analysis.diffs.title"), parent)
        self.diff_table = QTableWidget(0, 8)
        self.diff_table.setHorizontalHeaderLabels(
            [
                tr("analysis.col.pos"),
                tr("analysis.col.kind"),
                tr("analysis.col.ref"),
                tr("analysis.col.read"),
                tr("analysis.col.readpos"),
                tr("snp.col.verdict"),
                tr("snp.col.q"),
                tr("snp.col.ratio"),
            ]
        )
        self.diff_table.verticalHeader().setVisible(False)
        self.diff_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self.diff_table.setSortingEnabled(True)
        self.setWidget(self.diff_table)

    def fill(self, alignment: ReadAlignment | None, snp_calls: list, read_length: int) -> None:
        kind_labels = {
            "substitution": tr("analysis.kind.substitution"),
            "deletion": tr("analysis.kind.deletion"),
            "insertion": tr("analysis.kind.insertion"),
        }
        verdict_labels = {
            "homozygous": tr("snp.homozygous"),
            "heterozygous": tr("snp.heterozygous"),
            "low_quality": tr("snp.low_quality"),
        }
        calls_by_ref = {c.ref_pos: c for c in snp_calls}
        diffs = alignment.differences if alignment else []
        self.diff_table.setSortingEnabled(False)
        self.diff_table.setRowCount(len(diffs))
        for row, diff in enumerate(diffs):
            call = calls_by_ref.get(diff.ref_pos)
            verdict = ""
            if call is not None:
                verdict = verdict_labels.get(call.kind, call.kind)
                if call.kind == "heterozygous":
                    verdict += f" ({call.iupac})"
            values = [
                diff.ref_pos + 1,
                kind_labels.get(diff.kind, diff.kind),
                diff.ref_base,
                diff.read_base,
                (call.read_pos + 1 if call is not None
                 else shown + 1 if (shown := displayed_diff_index(alignment, diff, read_length)) is not None
                 else ""),
                verdict,
                call.quality if call is not None else "",
                f"{call.ratio:.2f}" if call is not None and call.ratio is not None else "",
            ]
            for col, value in enumerate(values):
                item = QTableWidgetItem()
                item.setData(Qt.ItemDataRole.DisplayRole, value)
                if call is not None and col == 5 and call.kind == "heterozygous":
                    item.setForeground(QColor("#a05000"))
                self.diff_table.setItem(row, col, item)
        self.diff_table.setSortingEnabled(True)

    def clear(self) -> None:
        self.diff_table.setRowCount(0)
        self.hide()
