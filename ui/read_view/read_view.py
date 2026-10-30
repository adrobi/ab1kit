"""Single-read viewer: chromatogram, minimap, editing, trim, search, analysis docks."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import QMarginsF, QSizeF, Qt, QUrl, Signal
from PySide6.QtGui import QAction, QDesktopServices, QPageSize, QPainter, QPdfWriter
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QDialog,
    QDockWidget,
    QFileDialog,
    QFormLayout,
    QInputDialog,
    QLabel,
    QMainWindow,
    QMenu,
    QMessageBox,
    QScrollBar,
    QSizePolicy,
    QToolBar,
    QVBoxLayout,
    QWidget,
)

from chrom_core import Chromatogram, EditableRead, read_ab1, trim_window
from chrom_core import blast as blast_api
from chrom_core import export, sessions
from chrom_core.align import ReadAlignment, align_read, gapped_sequences, parse_sequence
from chrom_core.config import AppConfig, save_config
from chrom_core.editing import INSERT_LETTERS, IUPAC_LETTERS
from chrom_core.snps import call_snps
from chrom_core.verify import parse_expected, verify_expected, verdicts_csv, verdicts_xlsx
from chrom_core.seqops import find_matches
from i18n import tr
from ui import status
from ui.chrom_view import ChromatogramWidget
from ui.minimap import MiniMap
from ui.read_view import menus
from ui.read_view.analysis_dock import DiffDock, displayed_diff_index, reference_overlay
from ui.read_view.batch_dock import BatchDock
from ui.read_view.dialogs import SearchDialog, TrimDialog, ask_reference_text
from ui.read_view.verify_dock import VerifyDock
from ui.status import bind_icon, quality_bucket
from ui.tables import repaint_statuses
from ui.theme import theme as palette, tick_image


class ReadView(QMainWindow):
    """Everything that belongs to one open read. Embedded by the shell window.

    The shell supplies ``start_task`` (background-task runner, see
    ``MainWindow._start_task``) and listens to the signals below; the view
    never touches its parent.
    """

    read_edited = Signal(Path)  # after a base edit, undo/redo, trim, revcomp, session save
    file_opened = Signal(str)
    status_message = Signal(str, int)  # text, timeout in ms (0 = until replaced)
    title_changed = Signal(str)

    def __init__(self, config: AppConfig, start_task: Callable, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.config = config
        self._start_task = start_task
        self._path: Path | None = None
        self._active = True    # False while the shell shows another tab
        self.read_model: EditableRead | None = None
        self.trim_dialog: TrimDialog | None = None
        self.reference_seq: str = ""
        self.alignment: ReadAlignment | None = None
        self.snp_calls: list = []
        self.variant_verdicts: list = []
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

        self.undo_action = QAction(tr("menu.edit.undo"), self)
        self.undo_action.setShortcut("Ctrl+Z")
        self.undo_action.triggered.connect(self.undo_read)
        self.undo_action.setEnabled(False)
        self.redo_action = QAction(tr("menu.edit.redo"), self)
        self.redo_action.setShortcut("Ctrl+Y")
        self.redo_action.triggered.connect(self.redo_read)
        self.redo_action.setEnabled(False)

        self._channel_boxes: dict[str, QCheckBox] = {}
        self._build_toolbar()
        self._build_docks()
        self._connect()

    def _say(self, text: str, timeout: int = 0) -> None:
        self.status_message.emit(text, timeout)

    # ------------------------------------------------------------------ setup
    def _build_toolbar(self) -> None:
        bar = QToolBar()
        bar.setMovable(False)
        bar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.addToolBar(bar)

        def add(text: str, icon_name: str, slot, shortcut: str | None = None) -> None:
            action = QAction(text, self)
            bind_icon(action, icon_name)
            if shortcut:
                action.setShortcut(shortcut)
            action.triggered.connect(slot)
            bar.addAction(action)

        add(tr("tb.open"), "folder-open", self.open_dialog)
        bar.addSeparator()
        add(tr("tb.zoom_in"), "zoom-in", self.chrom_view.zoom_in, "Ctrl++")
        add(tr("tb.zoom_out"), "zoom-out", self.chrom_view.zoom_out, "Ctrl+-")
        add(tr("tb.fit"), "maximize", self.chrom_view.zoom_fit)
        add(tr("tb.reset"), "rotate-ccw", self.chrom_view.zoom_reset)
        bar.addSeparator()
        bind_icon(self.undo_action, "undo-2")
        bind_icon(self.redo_action, "redo-2")
        bar.addAction(self.undo_action)
        bar.addAction(self.redo_action)
        bar.addSeparator()
        for base in "ACGT":
            box = QCheckBox(base)
            box.setChecked(True)
            box.toggled.connect(
                lambda on, b=base: self.chrom_view.set_visible_channel(b, on)
            )
            bar.addWidget(box)
            self._channel_boxes[base] = box
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        bar.addWidget(spacer)
        self.quality_legend = QLabel()
        self.quality_legend.setToolTip(tr("tb.quality_legend.tip"))
        bar.addWidget(self.quality_legend)
        self._refresh_legend("light")

    def _build_docks(self) -> None:
        self.info_dock = QDockWidget(tr("info.title"), self)
        self.info_dock.setWidget(QWidget())
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.info_dock)
        self.info_dock.hide()

        self.batch_dock = BatchDock(self.config, self._start_task)
        self.batch_panel = self.batch_dock.batch_panel
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.batch_dock)
        self.batch_dock.hide()
        self.batch_dock.open_requested.connect(self.open_path)
        self.batch_dock.message.connect(self._say)

        self.diff_dock = DiffDock()
        self.diff_table = self.diff_dock.diff_table
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.diff_dock)
        self.diff_dock.hide()

        self.verify_dock = VerifyDock()
        self.verify_table = self.verify_dock.verify_table
        self.addDockWidget(Qt.DockWidgetArea.BottomDockWidgetArea, self.verify_dock)
        self.verify_dock.hide()
        self.verify_dock.csv_button.clicked.connect(self._save_verify_csv)
        self.verify_dock.xlsx_button.clicked.connect(self._save_verify_xlsx)

    @property
    def batch_reports(self) -> list:
        return self.batch_dock.reports

    def current_model(self) -> EditableRead | None:
        return self.read_model

    # ------------------------------------------------------------------ menus
    # The shell owns the menu bar and its order; it asks the view to append
    # the read-specific items (see ``menus.py``) at the right places.
    def add_file_actions(self, menu: QMenu) -> None:
        menus.add_file_actions(self, menu)

    def add_edit_menu(self, menu: QMenu) -> None:
        menus.add_edit_actions(self, menu)

    def add_view_actions(self, menu: QMenu) -> None:
        menus.add_view_actions(self, menu)

    def add_analysis_menu(self, menu: QMenu) -> None:
        menus.add_analysis_actions(self, menu)

    def _connect(self) -> None:
        self.chrom_view.view_changed.connect(self._sync_from_view)
        self.chrom_view.cursor_info.connect(self._say)
        self.chrom_view.base_clicked.connect(self._on_base_clicked)
        self.minimap.request_center.connect(self.chrom_view.center_on_point)
        self.scrollbar.valueChanged.connect(
            lambda value: self.chrom_view.set_offset(float(value))
        )

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
        self.file_opened.emit(str(self._path))
        self._refresh_info(model.current())
        self._sync_actions()
        if self.reference_seq:
            self._realign()
        n = len(model.seq)
        mean_q = sum(model.quality) / max(1, n)
        message = tr("status.file", n=n, q=mean_q)
        if restored:
            message += " · " + tr("session.loaded")
        self._say(message, 8000)
        self.title_changed.emit(f"ab1kit — {self._path.name}" + ("*" if model.modified else ""))
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

    def open_dialog(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            tr("menu.file.open"),
            str(self._path.parent) if self._path else "",
            "Sanger chromatograms (*.ab1);;All files (*)",
        )
        if path:
            self.open_path(path)

    # --------------------------------------------------------------- analysis

    # --------------------------------------------------------------- analysis
    def _reference_dialog(self) -> None:
        text = ask_reference_text(self)
        if text is None:
            return
        try:
            sequence = parse_sequence(text)
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
            self._say(tr("analysis.need_read"), 5000)

    def _reference_overlay(self) -> tuple[str, dict[int, int]]:
        return reference_overlay(self.alignment, self.reference_seq, len(self.read_model.seq))

    def _displayed_diff_index(self, diff) -> int | None:
        return displayed_diff_index(self.alignment, diff, len(self.read_model.seq))

    def _clear_reference(self) -> None:
        self.reference_seq = ""
        self.alignment = None
        self.chrom_view.set_reference(None)
        self.diff_dock.clear()

    def _realign(self) -> None:
        if not (self.read_model and self.reference_seq):
            return
        try:
            self.alignment = align_read("".join(self.read_model.seq), self.reference_seq)
        except ValueError as exc:
            QMessageBox.warning(self, tr("analysis.ref.title"), str(exc))
            self._clear_reference()
            return
        letters, ref_to_read = self._reference_overlay()
        self.chrom_view.set_reference(letters, ref_to_read)
        chrom = self.read_model.current()
        self.snp_calls = call_snps(
            self.alignment, chrom,
            min_q=self.config.snp_min_q, het_ratio=self.config.het_ratio,
        )
        self.diff_dock.fill(self.alignment, self.snp_calls, len(self.read_model.seq))
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
        self._say(message, 8000)

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
        self._say(tr("export.done", path=path), 5000)

    def _blast_submit(self) -> None:
        if not self.read_model:
            return
        seq = "".join(self.read_model.seq)
        try:
            params = blast_api.build_request(seq)
        except (OSError, RuntimeError, ValueError) as exc:
            QMessageBox.warning(self, tr("menu.analysis.blast"), str(exc))
            return

        def done(rid: object) -> None:
            QDesktopServices.openUrl(QUrl(blast_api.result_url(str(rid))))
            self._say(tr("blast.sent", rid=rid), 10000)

        # urllib cannot be interrupted: on cancel the late RID is dropped.
        self._start_task(
            lambda progress, is_cancelled: blast_api.submit(params, timeout=20),
            done,
            tr("menu.analysis.blast"),
            key="blast",
        )

    def _verify_dialog(self) -> None:
        if not (self.read_model and self.alignment and self.reference_seq):
            self._say(tr("verify.need_reference"), 6000)
            return
        path, _ = QFileDialog.getOpenFileName(
            self,
            tr("menu.analysis.verify"),
            str(self._path.parent) if self._path else "",
            "Variants (*.vcf *.csv *.tsv *.txt);;All files (*)",
        )
        if not path:
            return
        try:
            text = Path(path).read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            QMessageBox.warning(self, tr("verify.title"), str(exc))
            return
        self._verify_variants(text)

    def _verify_variants(self, text: str) -> None:
        try:
            variants, warnings = parse_expected(text)
        except ValueError as exc:
            QMessageBox.warning(self, tr("verify.title"), str(exc))
            return
        if not variants:
            QMessageBox.warning(self, tr("verify.title"), tr("search.none"))
            return
        self.variant_verdicts = verify_expected(
            self.alignment,
            self.read_model.current(),
            self.reference_seq,
            variants,
            min_q=self.config.snp_min_q,
            het_ratio=self.config.het_ratio,
        )
        self._fill_verify_table()
        self.verify_dock.show()
        self.chrom_view.set_snp_marks(
            [(v.read_pos, v.verdict) for v in self.variant_verdicts if v.read_pos is not None]
        )
        counts = {k: sum(1 for v in self.variant_verdicts if v.verdict == k)
                  for k in ("confirmed", "heterozygous", "not_found", "no_data")}
        message = tr(
            "verify.summary",
            n=len(self.variant_verdicts),
            confirmed=counts["confirmed"],
            het=counts["heterozygous"],
            notfound=counts["not_found"],
            nodata=counts["no_data"],
        )
        if warnings:
            message += f" · {len(warnings)} ⚠"
        self._say(message, 10000)

    def _fill_verify_table(self) -> None:
        self.verify_dock.fill(self.variant_verdicts)

    def _save_verify_csv(self) -> None:
        if not self.variant_verdicts:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, tr("verify.save_csv"), "snp_report.csv", "CSV (*.csv)"
        )
        if path:
            export.write_text(path, verdicts_csv(self.variant_verdicts))
            self._say(tr("export.done", path=path), 5000)

    def _save_verify_xlsx(self) -> None:
        if not self.variant_verdicts:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, tr("verify.save_xlsx"), "snp_report.xlsx", "Excel (*.xlsx)"
        )
        if not path:
            return
        try:
            verdicts_xlsx(path, self.variant_verdicts, sample_id=self._export_id())
        except RuntimeError as exc:
            QMessageBox.warning(self, tr("verify.save_xlsx"), str(exc))
            return
        self._say(tr("export.done", path=path), 5000)

    # ------------------------------------------------------------------ batch
    def open_folder_dialog(self) -> None:
        start = str(self._path.parent) if self._path else ""
        folder = QFileDialog.getExistingDirectory(self, tr("menu.file.open_folder"), start)
        if folder:
            self._scan_folder(Path(folder))

    def _scan_folder(self, folder: Path) -> None:
        self.batch_dock.scan_folder(folder)

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

    def _refresh_legend(self, name: str) -> None:
        """Small colour key for the quality bars (same thresholds everywhere)."""
        colors = palette(name)
        parts = [tr("tb.quality_legend")]
        for q, text in ((0, "&lt;20"), (20, "20–29"), (30, "≥30")):
            parts.append(f'<span style="color:{colors[quality_bucket(q)]}">■</span>&nbsp;{text}')
        self.quality_legend.setText("&nbsp;&nbsp;".join(parts))

    def apply_theme(self, name: str) -> None:
        status.set_theme(name)
        colors = palette(name)
        self._refresh_legend(name)
        self.chrom_view.set_theme(colors)
        self.minimap.set_theme(colors)
        for base, box in self._channel_boxes.items():
            c = colors["bases"][base]
            box.setStyleSheet(
                f"QCheckBox {{ color: {c}; font-weight: bold; padding: 0 4px; }}"
                f"QCheckBox::indicator {{ border: 1px solid {c}; }}"
                f"QCheckBox::indicator:checked {{ background-color: {c}; border-color: {c};"
                f" image: url({tick_image(c)}); }}"
            )
        for table in (self.diff_table, self.verify_table, self.batch_panel.table):
            repaint_statuses(table)

    # --------------------------------------------------------------- editing
    def set_active(self, active: bool) -> None:
        """Enable or disable all of the view's own actions (menus, toolbar,
        shortcuts), so a hidden view cannot be edited by a stray shortcut."""
        self._active = bool(active)
        for action in self.findChildren(QAction, options=Qt.FindChildOption.FindDirectChildrenOnly):
            action.setEnabled(self._active)
        self._sync_actions()

    def _sync_actions(self) -> None:
        model = self.read_model
        self.undo_action.setEnabled(self._active and bool(model and model.can_undo))
        self.redo_action.setEnabled(self._active and bool(model and model.can_redo))

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
        self.title_changed.emit(
            f"ab1kit — {self._path.name}" + ("*" if self.read_model.modified else "")
        )
        if self.reference_seq:
            self._realign()

    def _after_edit(self) -> None:
        self.refresh_read()
        self.read_edited.emit(self._path)

    def apply_base_change(self, index: int, base: str) -> None:
        if self.read_model:
            self.read_model.set_base(index, base)
            self._after_edit()

    def apply_base_insert(self, index: int, base: str) -> None:
        if self.read_model:
            self.read_model.insert_base(index, base)
            self._after_edit()

    def apply_base_delete(self, index: int) -> None:
        if self.read_model:
            self.read_model.delete_base(index)
            self._after_edit()

    def apply_trim_range(self, start: int, end: int) -> None:
        if not self.read_model:
            return
        try:
            self.read_model.apply_trim(start, end)
        except ValueError as exc:
            QMessageBox.warning(self, tr("trim.title"), str(exc))
            return
        self.chrom_view.set_trim_region(None)
        self._after_edit()
        self._say(tr("trim.status", n=len(self.read_model.seq)))

    def undo_read(self) -> None:
        if self.read_model and self.read_model.undo():
            self._after_edit()

    def redo_read(self) -> None:
        if self.read_model and self.read_model.redo():
            self._after_edit()

    def reverse_complement_read(self) -> None:
        if self.read_model:
            self.read_model.reverse_complement()
            self._after_edit()

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
            self._say(
                tr("search.found", n=len(ranges), f=len(fwd), r=len(rev))
            )
            self._find_next()
        else:
            self._say(tr("search.none"))

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
            self._say(tr("export.done", path=path), 5000)

    def _export_fastq(self) -> None:
        if not self.read_model:
            return
        path = self._export_target(tr("menu.file.export_fastq"), ".fastq", "FASTQ (*.fastq)")
        if path:
            content = export.to_fastq(
                self._export_id(), "".join(self.read_model.seq), self.read_model.quality
            )
            export.write_text(path, content)
            self._say(tr("export.done", path=path), 5000)

    def _export_id(self) -> str:
        return self.read_model.current().sample_id if self.read_model else "read"

    def _export_png(self) -> None:
        path = self._export_target(tr("menu.file.export_png"), ".png", "PNG (*.png)")
        if path and self.chrom_view.grab().save(path):
            self._say(tr("export.done", path=path), 5000)

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
        self._say(tr("export.done", path=path), 5000)

    def _copy_sequence(self) -> None:
        if self.read_model:
            QApplication.clipboard().setText("".join(self.read_model.seq))
            self._say(
                tr("copy.done", n=len(self.read_model.seq)), 3000
            )

    def _save_session(self) -> None:
        if self.read_model is None or self._path is None:
            return
        sessions.save_session(
            sessions.session_path(self._path),
            {"source": self._path.name, **self.read_model.to_dict()},
        )
        self._say(tr("session.saved"), 3000)
        self.read_edited.emit(self._path)

    def autosave_session(self) -> None:
        """Called by the shell on close: save the sidecar if modified, else delete a stale one."""
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
