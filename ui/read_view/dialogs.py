"""Dialogs used by the single-read view."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from chrom_core import trim_window
from i18n import tr


class SearchDialog(QDialog):
    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle(tr("search.title"))
        layout = QVBoxLayout(self)
        self.edit = QLineEdit()
        self.edit.setPlaceholderText(tr("search.pattern"))
        self.reverse_box = QCheckBox(tr("search.reverse"))
        button = QPushButton(tr("search.button"))
        button.setDefault(True)
        button.clicked.connect(self.accept)
        layout.addWidget(self.edit)
        layout.addWidget(self.reverse_box)
        layout.addWidget(button)

    def pattern(self) -> str:
        return self.edit.text().strip().upper()

    def want_reverse(self) -> bool:
        return self.reverse_box.isChecked()


class TrimDialog(QDialog):
    """Non-modal quality-trim dialog; preview streams into the chromatogram."""

    region_changed = Signal(int, int)
    apply_requested = Signal(int, int)

    def __init__(self, quality_fn, window: int, min_q: int, total: int, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle(tr("trim.title"))
        self.quality_fn = quality_fn
        self.total = total
        self._region = (0, total)

        layout = QVBoxLayout(self)
        form = QFormLayout()
        self.window_spin = QSpinBox()
        self.window_spin.setRange(2, 100)
        self.window_spin.setValue(window)
        self.minq_spin = QSpinBox()
        self.minq_spin.setRange(0, 60)
        self.minq_spin.setValue(min_q)
        form.addRow(QLabel(tr("trim.window")), self.window_spin)
        form.addRow(QLabel(tr("trim.minq")), self.minq_spin)
        layout.addLayout(form)
        self.preview = QLabel()
        layout.addWidget(self.preview)
        buttons = QWidget()
        row = QVBoxLayout(buttons)
        apply_button = QPushButton(tr("trim.apply"))
        apply_button.setDefault(True)
        apply_button.clicked.connect(self._apply)
        close_button = QPushButton(tr("dialog.close"))
        close_button.clicked.connect(self.close)
        row.addWidget(apply_button)
        row.addWidget(close_button)
        layout.addWidget(buttons)

        self.window_spin.valueChanged.connect(self._recompute)
        self.minq_spin.valueChanged.connect(self._recompute)
        self._recompute()

    def window_value(self) -> int:
        return self.window_spin.value()

    def minq_value(self) -> int:
        return self.minq_spin.value()

    def _recompute(self) -> None:
        start, end = trim_window(self.quality_fn(), self.window_spin.value(), self.minq_spin.value())
        self._set_region(start, end)
        self.region_changed.emit(start, end)

    def _set_region(self, start: int, end: int) -> None:
        self._region = (start, end)
        self.preview.setText(tr("trim.preview", n=end - start, total=self.total))

    def set_region(self, start: int, end: int) -> None:
        """Sync preview when the user drags trim handles in the view."""
        self._set_region(start, end)

    def _apply(self) -> None:
        start, end = self._region
        self.apply_requested.emit(start, end)
        self.close()


def ask_reference_text(parent: QWidget) -> str | None:
    """Modal dialog for pasting/loading a reference; returns the raw text
    when accepted, ``None`` otherwise."""
    dialog = QDialog(parent)
    dialog.setWindowTitle(tr("analysis.ref.title"))
    dialog.resize(480, 240)
    layout = QVBoxLayout(dialog)
    edit = QPlainTextEdit()
    edit.setPlaceholderText(tr("analysis.ref.paste"))
    layout.addWidget(edit)

    def from_file() -> None:
        path, _ = QFileDialog.getOpenFileName(
            dialog, tr("analysis.ref.from_file"), "", "Sequence (*.fasta *.fa *.txt);;All files (*)"
        )
        if path:
            edit.setPlainText(Path(path).read_text(encoding="utf-8", errors="replace"))

    file_button = QPushButton(tr("analysis.ref.from_file"))
    file_button.clicked.connect(from_file)
    layout.addWidget(file_button)

    ok = QPushButton("OK")
    ok.clicked.connect(dialog.accept)
    cancel = QPushButton(tr("dialog.close"))
    cancel.clicked.connect(dialog.reject)
    row = QWidget()
    row_layout = QHBoxLayout(row)
    row_layout.addWidget(ok)
    row_layout.addWidget(cancel)
    layout.addWidget(row)

    if dialog.exec() != QDialog.DialogCode.Accepted:
        return None
    return edit.toPlainText()
