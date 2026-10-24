"""Main ab1kit window: toolbar, menus, overview map, chromatogram, docks."""

from __future__ import annotations

import sys
from importlib.metadata import PackageNotFoundError, version as pkg_version
from pathlib import Path

from PySide6.QtCore import QMarginsF, QSizeF, Qt, Signal
from PySide6.QtGui import QAction, QColor, QDragEnterEvent, QDropEvent, QPageSize, QPainter, QPdfWriter
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QDialog,
    QDockWidget,
    QFileDialog,
    QFormLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMenu,
    QMessageBox,
    QPushButton,
    QScrollBar,
    QSpinBox,
    QPlainTextEdit,
    QTableWidget,
    QTableWidgetItem,
    QToolBar,
    QVBoxLayout,
    QWidget,
)

from chrom_core import Chromatogram, EditableRead, batch, read_ab1, trim_window
from chrom_core import export, sessions
from chrom_core.align import ReadAlignment, align_read, gapped_sequences, parse_sequence
from chrom_core.config import AppConfig, load_config, save_config
from chrom_core.editing import INSERT_LETTERS, IUPAC_LETTERS
from chrom_core.snps import call_snps
from chrom_core.seqops import find_matches
from i18n import set_language, tr
from ui.batch_panel import BatchPanel
from ui.chrom_view import ChromatogramWidget
from ui.minimap import MiniMap
from ui.theme import app_qss, theme as palette

try:
    VERSION = pkg_version("ab1kit")
except PackageNotFoundError:  # running from source without install
    VERSION = "dev"


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


class MainWindow(QMainWindow):
    def __init__(self, config: AppConfig | None = None) -> None:
        super().__init__()
        self.config = config or load_config()
        set_language(self.config.language)

        self.setWindowTitle("ab1kit")
        self.resize(1200, 780)
        self._path: Path | None = None
        self.read_model: EditableRead | None = None
        self.trim_dialog: TrimDialog | None = None
        self.batch_reports: list[batch.FileReport] = []
        self.reference_seq: str = ""
        self.alignment: ReadAlignment | None = None
        self.snp_calls: list = []
        self._search_ranges: list[tuple[int, int]] = []
        self._search_pos = -1

        self.chrom_view = ChromatogramWidget()
        self.minimap = MiniMap()
        self.scrollbar = QScrollBar(Qt.Orientation.Horizontal)

        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self.minimap)
        layout.addWidget(self.chrom_view, 1)
        layout.addWidget(self.scrollbar)
        self.setCentralWidget(central)

        self._channel_boxes: dict[str, QCheckBox] = {}
        self._build_toolbar()
        self._build_info_dock()
        self._build_batch_dock()
        self._build_analysis_dock()
        self._build_menus()
        self._connect()
        self.setAcceptDrops(True)
        self.set_theme(self.config.theme)
        self.statusBar().showMessage(tr("status.ready"))
        self._fill_recent_menu()

    # ------------------------------------------------------------------ setup
    def _build_toolbar(self) -> None:
        bar = QToolBar()
        bar.setMovable(False)
        self.addToolBar(bar)

        def add(text: str, slot, shortcut: str | None = None) -> None:
            action = QAction(text, self)
            if shortcut:
                action.setShortcut(shortcut)
            action.triggered.connect(slot)
            bar.addAction(action)

        add(tr("tb.open"), self._open_dialog)
        bar.addSeparator()
        add(tr("tb.zoom_in"), self.chrom_view.zoom_in, "Ctrl++")
        add(tr("tb.zoom_out"), self.chrom_view.zoom_out, "Ctrl+-")
        add(tr("tb.fit"), self.chrom_view.zoom_fit)
        add(tr("tb.reset"), self.chrom_view.zoom_reset)
        bar.addSeparator()
        for base in "ACGT":
            box = QCheckBox(base)
            box.setChecked(True)
            box.toggled.connect(
                lambda on, b=base: self.chrom_view.set_visible_channel(b, on)
            )
            bar.addWidget(box)
            self._channel_boxes[base] = box

    def _build_menus(self) -> None:
        file_menu = self.menuBar().addMenu(tr("menu.file"))
        self._file_menu = file_menu
        open_action = QAction(tr("menu.file.open"), self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self._open_dialog)
        file_menu.addAction(open_action)
        open_folder = QAction(tr("menu.file.open_folder"), self)
        open_folder.setShortcut("Ctrl+Shift+O")
        open_folder.triggered.connect(self._open_folder_dialog)
        file_menu.addAction(open_folder)
        self.recent_menu = file_menu.addMenu(tr("menu.file.recent"))
        save_session = QAction(tr("menu.file.save_session"), self)
        save_session.setShortcut("Ctrl+S")
        save_session.triggered.connect(self._save_session)
        file_menu.addAction(save_session)
        export_menu = file_menu.addMenu(tr("menu.file.export"))
        for text, slot in (
            (tr("menu.file.export_fasta"), self._export_fasta),
            (tr("menu.file.export_fastq"), self._export_fastq),
            (tr("menu.file.export_png"), self._export_png),
            (tr("menu.file.export_pdf"), self._export_pdf),
        ):
            action = QAction(text, self)
            action.triggered.connect(slot)
            export_menu.addAction(action)
        file_menu.addSeparator()
        copy_action = QAction(tr("menu.file.copy_seq"), self)
        copy_action.setShortcut("Ctrl+Shift+C")
        copy_action.triggered.connect(self._copy_sequence)
        file_menu.addAction(copy_action)
        file_menu.addSeparator()
        exit_action = QAction(tr("menu.file.exit"), self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        edit_menu = self.menuBar().addMenu(tr("menu.edit"))
        self.undo_action = QAction(tr("menu.edit.undo"), self)
        self.undo_action.setShortcut("Ctrl+Z")
        self.undo_action.triggered.connect(self.undo_read)
        self.undo_action.setEnabled(False)
        edit_menu.addAction(self.undo_action)
        self.redo_action = QAction(tr("menu.edit.redo"), self)
        self.redo_action.setShortcut("Ctrl+Y")
        self.redo_action.triggered.connect(self.redo_read)
        self.redo_action.setEnabled(False)
        edit_menu.addAction(self.redo_action)
        edit_menu.addSeparator()
        find_action = QAction(tr("menu.edit.find"), self)
        find_action.setShortcut("Ctrl+F")
        find_action.triggered.connect(self._search)
        edit_menu.addAction(find_action)
        find_next = QAction(tr("menu.edit.find_next"), self)
        find_next.setShortcut("F3")
        find_next.triggered.connect(self._find_next)
        edit_menu.addAction(find_next)
        goto_action = QAction(tr("menu.edit.goto"), self)
        goto_action.setShortcut("Ctrl+G")
        goto_action.triggered.connect(self._goto)
        edit_menu.addAction(goto_action)
        edit_menu.addSeparator()
        revcomp_action = QAction(tr("menu.edit.revcomp"), self)
        revcomp_action.setShortcut("Ctrl+R")
        revcomp_action.triggered.connect(self.reverse_complement_read)
        edit_menu.addAction(revcomp_action)
        trim_action = QAction(tr("menu.edit.trim"), self)
        trim_action.setShortcut("Ctrl+T")
        trim_action.triggered.connect(self._open_trim_dialog)
        edit_menu.addAction(trim_action)

        view_menu = self.menuBar().addMenu(tr("menu.view"))
        self.dark_action = QAction(tr("menu.view.dark_theme"), self)
        self.dark_action.setCheckable(True)
        self.dark_action.setChecked(self.config.theme == "dark")
        self.dark_action.toggled.connect(
            lambda on: self.set_theme("dark" if on else "light")
        )
        view_menu.addAction(self.dark_action)
        smooth_action = QAction(tr("menu.view.smoothing"), self)
        smooth_action.setCheckable(True)
        smooth_action.toggled.connect(self.chrom_view.set_smoothing)
        view_menu.addAction(smooth_action)
        profile_action = QAction(tr("menu.view.quality_profile"), self)
        profile_action.setCheckable(True)
        profile_action.toggled.connect(self.minimap.set_profile)
        view_menu.addAction(profile_action)
        info_action = QAction(tr("menu.view.file_info"), self)
        info_action.setCheckable(True)
        info_action.toggled.connect(self.info_dock.setVisible)
        self.info_dock.visibilityChanged.connect(info_action.setChecked)
        view_menu.addAction(info_action)
        batch_action = QAction(tr("menu.view.batch_panel"), self)
        batch_action.setCheckable(True)
        batch_action.toggled.connect(self.batch_dock.setVisible)
        self.batch_dock.visibilityChanged.connect(batch_action.setChecked)
        view_menu.addAction(batch_action)

        analysis_menu = self.menuBar().addMenu(tr("menu.analysis"))
        reference_action = QAction(tr("menu.analysis.reference"), self)
        reference_action.setShortcut("Ctrl+L")
        reference_action.triggered.connect(self._reference_dialog)
        analysis_menu.addAction(reference_action)
        clear_ref = QAction(tr("menu.analysis.clear_reference"), self)
        clear_ref.triggered.connect(self._clear_reference)
        analysis_menu.addAction(clear_ref)
        analysis_menu.addSeparator()
        export_aligned = QAction(tr("menu.analysis.export_aligned"), self)
        export_aligned.triggered.connect(self._export_aligned)
        analysis_menu.addAction(export_aligned)

        help_menu = self.menuBar().addMenu(tr("menu.help"))
        about_action = QAction(tr("menu.help.about"), self)
        about_action.triggered.connect(
            lambda: QMessageBox.about(self, "ab1kit", tr("about.text", version=VERSION))
        )
        help_menu.addAction(about_action)

    def _build_info_dock(self) -> None:
        self.info_dock = QDockWidget(tr("info.title"), self)
        self.info_dock.setWidget(QWidget())
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.info_dock)
        self.info_dock.hide()

    def _build_batch_dock(self) -> None:
        self.batch_dock = QDockWidget(tr("batch.title"), self)
        self.batch_panel = BatchPanel()
        self.batch_dock.setWidget(self.batch_panel)
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.batch_dock)
        self.batch_dock.hide()
        self.batch_panel.open_requested.connect(self.open_path)
        self.batch_panel.export_button.clicked.connect(self._batch_export)
        self.batch_panel.csv_button.clicked.connect(self._batch_save_csv)

    def _build_analysis_dock(self) -> None:
        self.diff_dock = QDockWidget(tr("analysis.diffs.title"), self)
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
        self.diff_dock.setWidget(self.diff_table)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.diff_dock)
        self.diff_dock.hide()

    def _connect(self) -> None:
        self.chrom_view.view_changed.connect(self._sync_from_view)
        self.chrom_view.cursor_info.connect(self.statusBar().showMessage)
        self.chrom_view.base_clicked.connect(self._on_base_clicked)
        self.minimap.request_center.connect(self.chrom_view.center_on_point)
        self.scrollbar.valueChanged.connect(
            lambda value: self.chrom_view.set_offset(float(value))
        )

    # ---------------------------------------------------------------- opening
    def open_path(self, path: str | Path) -> bool:
        try:
            chrom = read_ab1(path)
        except Exception as exc:  # user-facing: any parse failure becomes a dialog
            QMessageBox.critical(self, tr("err.open_failed"), f"{path}:\n{exc}")
            return False
        self._path = Path(path)
        model = EditableRead(chrom)
        restored = self._try_restore_session(model)
        self.read_model = model
        self._search_ranges = []
        self._search_pos = -1
        self.chrom_view.set_chromatogram(model.current())
        self.chrom_view.set_edited(model.edited)
        self.minimap.set_chromatogram(model.current())
        self.minimap.set_highlights([])
        self._push_recent(str(self._path))
        self._refresh_info(model.current())
        self._sync_actions()
        if self.reference_seq:
            self._realign()
        n = len(model.seq)
        mean_q = sum(model.quality) / max(1, n)
        message = tr("status.file", n=n, q=mean_q)
        if restored:
            message += " · " + tr("session.loaded")
        self.statusBar().showMessage(message, 8000)
        self.setWindowTitle(f"ab1kit — {self._path.name}" + ("*" if model.modified else ""))
        return True

    def _try_restore_session(self, model: EditableRead) -> bool:
        sidecar = sessions.session_path(self._path)
        data = sessions.load_session(sidecar)
        if not data or data.get("source") != self._path.name:
            return False
        try:
            model.restore_state(data)
        except (KeyError, TypeError, ValueError):
            return False
        return True

    def _open_dialog(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            tr("menu.file.open"),
            str(self._path.parent) if self._path else "",
            "Sanger chromatograms (*.ab1);;All files (*)",
        )
        if path:
            self.open_path(path)

    # --------------------------------------------------------------- analysis
    def _reference_dialog(self) -> None:
        dialog = QDialog(self)
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
            return
        try:
            sequence = parse_sequence(edit.toPlainText())
        except ValueError as exc:
            QMessageBox.warning(self, tr("analysis.ref.title"), str(exc))
            return
        if not sequence:
            return
        self._set_reference_seq(sequence)

    def _set_reference_seq(self, sequence: str) -> None:
        self.reference_seq = sequence
        if self.read_model is not None:
            self._realign()
        else:
            self.statusBar().showMessage(tr("analysis.need_read"), 5000)

    def _clear_reference(self) -> None:
        self.reference_seq = ""
        self.alignment = None
        self.chrom_view.set_reference(None)
        self.diff_table.setRowCount(0)
        self.diff_dock.hide()

    def _realign(self) -> None:
        if not (self.read_model and self.reference_seq):
            return
        try:
            self.alignment = align_read("".join(self.read_model.seq), self.reference_seq)
        except ValueError as exc:
            QMessageBox.warning(self, tr("analysis.ref.title"), str(exc))
            self._clear_reference()
            return
        self.chrom_view.set_reference(self.reference_seq, self.alignment.ref_to_read)
        chrom = self.read_model.current()
        self.snp_calls = call_snps(
            self.alignment, chrom,
            min_q=self.config.snp_min_q, het_ratio=self.config.het_ratio,
        )
        self._fill_diff_table()
        self.diff_dock.setVisible(bool(self.alignment.differences))
        pct = f"{self.alignment.identity * 100:.1f}"
        if self.alignment.differences:
            message = tr(
                "analysis.status",
                n=len(self.alignment.differences),
                pct=pct,
                rstart=self.alignment.ref_start + 1,
                rend=self.alignment.ref_end,
            )
        else:
            message = tr("analysis.no_diffs")
        if self.snp_calls:
            counts = {k: sum(1 for c in self.snp_calls if c.kind == k)
                      for k in ("homozygous", "heterozygous", "low_quality")}
            message += " · " + tr(
                "snp.summary",
                homo=counts["homozygous"],
                het=counts["heterozygous"],
                low=counts["low_quality"],
            )
        if self.alignment.alternatives > 1:
            message += " · " + tr("analysis.ambiguous", n=self.alignment.alternatives)
        self.statusBar().showMessage(message, 8000)

    def _fill_diff_table(self) -> None:
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
        calls_by_ref = {c.ref_pos: c for c in getattr(self, "snp_calls", [])}
        diffs = self.alignment.differences if self.alignment else []
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
                 else diff.read_pos + 1 if diff.read_pos is not None else ""),
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

    def _export_aligned(self) -> None:
        if not (self.alignment and self.reference_seq and self.read_model):
            return
        path = self._export_target(
            tr("menu.analysis.export_aligned"), ".aligned.fasta", "FASTA (*.fasta)"
        )
        if not path:
            return
        ref_aligned, read_aligned = gapped_sequences(self.alignment, self.reference_seq)
        sample = self._export_id()
        strand = "reverse" if self.alignment.strand == "-" else "forward"
        content = (
            export.to_fasta(f"{sample}_reference[{self.alignment.ref_start + 1}-{self.alignment.ref_end}]", ref_aligned)
            + export.to_fasta(f"{sample}_read_{strand}", read_aligned)
        )
        export.write_text(path, content)
        self.statusBar().showMessage(tr("export.done", path=path), 5000)

    # ------------------------------------------------------------------ batch
    def _open_folder_dialog(self) -> None:
        start = str(self._path.parent) if self._path else ""
        folder = QFileDialog.getExistingDirectory(self, tr("menu.file.open_folder"), start)
        if folder:
            self._scan_folder(Path(folder))

    def _scan_folder(self, folder: Path) -> None:
        self.batch_reports = batch.scan_folder(
            folder, self.config.trim_window, self.config.trim_min_q
        )
        self.batch_panel.show_reports(self.batch_reports)
        self.batch_dock.show()
        self.statusBar().showMessage(tr("batch.scan_done", n=len(self.batch_reports)), 5000)

    def _batch_export(self) -> None:
        if not self.batch_reports:
            return
        out_dir = QFileDialog.getExistingDirectory(self, tr("batch.export_dir"))
        if not out_dir:
            return
        results = batch.batch_export(
            self.batch_reports,
            Path(out_dir),
            window=self.config.trim_window,
            min_q=self.config.trim_min_q,
        )
        exported = sum(1 for _, dest, error in results if dest and not error)
        errors = sum(1 for _, _, error in results if error)
        self.statusBar().showMessage(
            tr("batch.exported", n=exported, err=errors), 8000
        )

    def _batch_save_csv(self) -> None:
        if not self.batch_reports:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, tr("batch.csv"), "batch_summary.csv", "CSV (*.csv)"
        )
        if path:
            export.write_text(path, batch.summary_csv(self.batch_reports))
            self.statusBar().showMessage(tr("export.done", path=path), 5000)

    def dragEnterEvent(self, event: QDragEnterEvent) -> None:  # noqa: N802
        urls = event.mimeData().urls()
        if urls and urls[0].toLocalFile().lower().endswith(".ab1"):
            event.acceptProposedAction()

    def dropEvent(self, event: QDropEvent) -> None:  # noqa: N802
        path = event.mimeData().urls()[0].toLocalFile()
        if path:
            self.open_path(path)

    # ---------------------------------------------------------------- recent
    def _push_recent(self, path: str) -> None:
        recent = [p for p in self.config.recent_files if p != path]
        recent.insert(0, path)
        self.config.recent_files = recent[:10]
        save_config(self.config)
        self._fill_recent_menu()

    def _fill_recent_menu(self) -> None:
        self.recent_menu.clear()
        for path in self.config.recent_files:
            action = QAction(Path(path).name, self)
            action.setToolTip(path)
            action.triggered.connect(
                lambda _=False, p=path: Path(p).exists() and self.open_path(p)
            )
            self.recent_menu.addAction(action)

    # ------------------------------------------------------------------ view
    def _sync_from_view(self) -> None:
        total = self.chrom_view._total_points()
        visible = self.chrom_view._visible_points()
        self.scrollbar.blockSignals(True)
        self.scrollbar.setRange(0, max(0, int(total - visible)))
        self.scrollbar.setPageStep(max(1, int(visible)))
        self.scrollbar.setValue(int(self.chrom_view.offset))
        self.scrollbar.blockSignals(False)
        self.minimap.set_view(self.chrom_view.offset, visible, total)

    def set_theme(self, name: str) -> None:
        self.config.theme = name
        colors = palette(name)
        self.chrom_view.set_theme(colors)
        self.minimap.set_theme(colors)
        for base, box in self._channel_boxes.items():
            box.setStyleSheet(
                f"QCheckBox {{ color: {colors['bases'][base]}; font-weight: bold; padding: 0 4px; }}"
            )
        QApplication.instance().setStyleSheet(app_qss(name == "dark"))
        save_config(self.config)

    # --------------------------------------------------------------- editing
    def _sync_actions(self) -> None:
        model = self.read_model
        self.undo_action.setEnabled(bool(model and model.can_undo))
        self.redo_action.setEnabled(bool(model and model.can_redo))

    def refresh_read(self) -> None:
        if self.read_model is None:
            return
        chrom = self.read_model.current()
        self.chrom_view.set_chromatogram(chrom, keep_view=True)
        self.chrom_view.set_edited(self.read_model.edited)
        self.minimap.set_chromatogram(chrom)
        self.minimap.set_highlights([])
        self._search_ranges = []
        self._search_pos = -1
        self._refresh_info(chrom)
        self._sync_actions()
        self.setWindowTitle(
            f"ab1kit — {self._path.name}" + ("*" if self.read_model.modified else "")
        )
        if self.reference_seq:
            self._realign()

    def apply_base_change(self, index: int, base: str) -> None:
        if self.read_model:
            self.read_model.set_base(index, base)
            self.refresh_read()

    def apply_base_insert(self, index: int, base: str) -> None:
        if self.read_model:
            self.read_model.insert_base(index, base)
            self.refresh_read()

    def apply_base_delete(self, index: int) -> None:
        if self.read_model:
            self.read_model.delete_base(index)
            self.refresh_read()

    def apply_trim_range(self, start: int, end: int) -> None:
        if not self.read_model:
            return
        try:
            self.read_model.apply_trim(start, end)
        except ValueError as exc:
            QMessageBox.warning(self, tr("trim.title"), str(exc))
            return
        self.chrom_view.set_trim_region(None)
        self.refresh_read()
        self.statusBar().showMessage(tr("trim.status", n=len(self.read_model.seq)))

    def undo_read(self) -> None:
        if self.read_model and self.read_model.undo():
            self.refresh_read()

    def redo_read(self) -> None:
        if self.read_model and self.read_model.redo():
            self.refresh_read()

    def reverse_complement_read(self) -> None:
        if self.read_model:
            self.read_model.reverse_complement()
            self.refresh_read()

    def _on_base_clicked(self, index: int, global_pos) -> None:
        if self.read_model is None:
            return
        menu = QMenu(self)
        change = menu.addMenu(tr("edit.change"))
        for letter in IUPAC_LETTERS:
            change.addAction(letter, lambda l=letter: self.apply_base_change(index, l))
        insert_left = menu.addMenu(tr("edit.insert_left"))
        for letter in INSERT_LETTERS:
            insert_left.addAction(letter, lambda l=letter: self.apply_base_insert(index, l))
        insert_right = menu.addMenu(tr("edit.insert_right"))
        for letter in INSERT_LETTERS:
            insert_right.addAction(letter, lambda l=letter: self.apply_base_insert(index + 1, l))
        menu.addSeparator()
        menu.addAction(tr("edit.delete"), lambda: self.apply_base_delete(index))
        menu.exec(global_pos)

    # ------------------------------------------------------------------ trim
    def _open_trim_dialog(self) -> None:
        if self.read_model is None or self.trim_dialog is not None:
            return
        start, end = trim_window(
            self.read_model.quality, self.config.trim_window, self.config.trim_min_q
        )
        self.chrom_view.set_trim_region((start, end))
        dialog = TrimDialog(
            lambda: self.read_model.quality,
            self.config.trim_window,
            self.config.trim_min_q,
            len(self.read_model.seq),
            self,
        )
        self.trim_dialog = dialog
        dialog.region_changed.connect(
            lambda start, end: self.chrom_view.set_trim_region((start, end))
        )
        self.chrom_view.trim_range_changed.connect(dialog.set_region)
        dialog.apply_requested.connect(self._apply_trim_from_dialog)
        dialog.finished.connect(self._trim_dialog_closed)
        dialog.show()

    def _apply_trim_from_dialog(self, start: int, end: int) -> None:
        if self.trim_dialog is not None:
            self.config.trim_window = self.trim_dialog.window_value()
            self.config.trim_min_q = self.trim_dialog.minq_value()
            save_config(self.config)
        self.apply_trim_range(start, end)

    def _trim_dialog_closed(self) -> None:
        dialog = self.trim_dialog
        self.trim_dialog = None
        if dialog is not None:
            try:
                self.chrom_view.trim_range_changed.disconnect(dialog.set_region)
            except RuntimeError:
                pass
        self.chrom_view.set_trim_region(None)

    # ---------------------------------------------------------------- search
    def _search(self) -> None:
        if self.chrom_view.chrom is None:
            return
        dialog = SearchDialog(self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        pattern = dialog.pattern()
        if not pattern:
            return
        seq = self.chrom_view.chrom.seq
        try:
            fwd = find_matches(pattern, seq)
            rev = find_matches(pattern, seq, reverse=True) if dialog.want_reverse() else []
        except ValueError as exc:
            QMessageBox.warning(self, tr("search.title"), str(exc))
            return
        ranges = [(i, i + len(pattern)) for i in fwd + rev]
        self.chrom_view.set_highlights(ranges)
        self.minimap.set_highlights(ranges)
        self._search_ranges = ranges
        self._search_pos = -1
        if ranges:
            self.statusBar().showMessage(
                tr("search.found", n=len(ranges), f=len(fwd), r=len(rev))
            )
            self._find_next()
        else:
            self.statusBar().showMessage(tr("search.none"))

    def _find_next(self) -> None:
        if not self._search_ranges:
            return
        self._search_pos = (self._search_pos + 1) % len(self._search_ranges)
        self.chrom_view.goto_base(self._search_ranges[self._search_pos][0] + 1)

    # ----------------------------------------------------------- goto and info
    def _goto(self) -> None:
        chrom = self.chrom_view.chrom
        if chrom is None:
            return
        number, ok = QInputDialog.getInt(
            self, tr("goto.title"), tr("goto.label"), 1, 1, len(chrom.seq)
        )
        if ok:
            self.chrom_view.goto_base(number)

    def _refresh_info(self, chrom: Chromatogram) -> None:
        page = QWidget()
        form = QFormLayout(page)
        form.setContentsMargins(8, 8, 8, 8)

        def row(label: str, value) -> None:
            widget = QLabel(str(value))
            widget.setWordWrap(True)
            form.addRow(QLabel(label), widget)

        row(tr("info.file"), self._path)
        row(tr("info.length"), len(chrom.seq))
        row(tr("info.mean_q"), f"{sum(chrom.quality) / max(1, len(chrom.quality)):.1f}")
        row(tr("info.q20"), sum(1 for q in chrom.quality if q >= 20))
        meta_keys = (
            ("sample", "meta.sample"),
            ("instrument", "meta.instrument"),
            ("run_date", "meta.run_date"),
            ("run_time", "meta.run_time"),
            ("project", "meta.project"),
            ("capillary", "meta.capillary"),
            ("mobility", "meta.mobility"),
            ("comment", "meta.comment"),
        )
        for key, label_key in meta_keys:
            if key in chrom.metadata:
                row(tr(label_key), chrom.metadata[key])
        self.info_dock.setWidget(page)

    # ------------------------------------------------------- export and sessions
    def _export_target(self, caption: str, suffix: str, pattern: str) -> str | None:
        default = (self._path.stem if self._path else "read") + suffix
        path, _ = QFileDialog.getSaveFileName(self, caption, default, pattern)
        return path or None

    def _export_fasta(self) -> None:
        if not self.read_model:
            return
        path = self._export_target(tr("menu.file.export_fasta"), ".fasta", "FASTA (*.fasta)")
        if path:
            export.write_text(path, export.to_fasta(self._export_id(), "".join(self.read_model.seq)))
            self.statusBar().showMessage(tr("export.done", path=path), 5000)

    def _export_fastq(self) -> None:
        if not self.read_model:
            return
        path = self._export_target(tr("menu.file.export_fastq"), ".fastq", "FASTQ (*.fastq)")
        if path:
            content = export.to_fastq(
                self._export_id(), "".join(self.read_model.seq), self.read_model.quality
            )
            export.write_text(path, content)
            self.statusBar().showMessage(tr("export.done", path=path), 5000)

    def _export_id(self) -> str:
        return self.read_model.current().sample_id if self.read_model else "read"

    def _export_png(self) -> None:
        path = self._export_target(tr("menu.file.export_png"), ".png", "PNG (*.png)")
        if path and self.chrom_view.grab().save(path):
            self.statusBar().showMessage(tr("export.done", path=path), 5000)

    def _export_pdf(self) -> None:
        path = self._export_target(tr("menu.file.export_pdf"), ".pdf", "PDF (*.pdf)")
        if not path:
            return
        pixmap = self.chrom_view.grab()
        writer = QPdfWriter(path)
        writer.setPageSize(QPageSize(QSizeF(pixmap.width(), pixmap.height()), QPageSize.Unit.Point))
        writer.setPageMargins(QMarginsF(0, 0, 0, 0))
        painter = QPainter(writer)
        painter.drawPixmap(0, 0, writer.width(), writer.height(), pixmap)
        painter.end()
        self.statusBar().showMessage(tr("export.done", path=path), 5000)

    def _copy_sequence(self) -> None:
        if self.read_model:
            QApplication.clipboard().setText("".join(self.read_model.seq))
            self.statusBar().showMessage(
                tr("copy.done", n=len(self.read_model.seq)), 3000
            )

    def _save_session(self) -> None:
        if self.read_model is None or self._path is None:
            return
        sessions.save_session(
            sessions.session_path(self._path),
            {"source": self._path.name, **self.read_model.to_dict()},
        )
        self.statusBar().showMessage(tr("session.saved"), 3000)

    def closeEvent(self, event) -> None:  # noqa: N802
        if self.read_model is not None and self._path is not None:
            sidecar = sessions.session_path(self._path)
            try:
                if self.read_model.modified:
                    sessions.save_session(
                        sidecar, {"source": self._path.name, **self.read_model.to_dict()}
                    )
                elif sidecar.exists():
                    sidecar.unlink()
            except OSError:
                pass
        super().closeEvent(event)


def run(argv: list[str] | None = None) -> int:
    """Application entry point used by the ``ab1kit`` console script."""
    args = list(sys.argv[1:] if argv is None else argv)
    app = QApplication(args)
    app.setStyle("Fusion")  # chrome identical on every platform and OS theme
    window = MainWindow()
    window.show()
    for arg in args:
        if not arg.startswith("-") and Path(arg).exists():
            window.open_path(arg)
            break
    return app.exec()
