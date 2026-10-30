"""The "construct check" workspace: construct + reads in, per-clone verdicts out.

Holds the session state (construct file, circular flag, region, grouping,
thresholds), re-runs ``verify_all`` in the background on every change (key
``"verify"``: a newer run cancels the older one; a cancelled or failed run
leaves the previous verdicts on screen) and reads/writes ``.ab1proj`` files.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QDragEnterEvent, QDropEvent
from PySide6.QtWidgets import (
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QMessageBox,
    QScrollArea,
    QSplitter,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from chrom_core import project as project_io
from chrom_core.align import MAX_REFERENCE_LENGTH
from chrom_core.clones import (
    CloneReport,
    Grouping,
    Region,
    VerifyParams,
    check_region,
    group_reads,
    verify_all,
)
from chrom_core.config import AppConfig, save_config
from chrom_core.reference import Construct, load_construct
from i18n import tr
from ui import status
from ui.construct.controls import ControlsBar, EmptyState, dropped_paths
from ui.construct.detail_pane import DetailPane
from ui.construct.grouping_dialog import GroupingDialog
from ui.construct.matrix import CloneMatrix
from ui.construct.session_files import (
    CONSTRUCT_FILTER,
    PROJECT_SUFFIX,
    ab1_files,
    path_key,
    write_report,
)
from ui.construct.summary_bar import SummaryBar



class ConstructWorkspace(QWidget):
    open_read = Signal(Path)
    status_message = Signal(str, int)

    def __init__(self, config: AppConfig, start_task: Callable,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.config = config
        self._start_task = start_task
        self._construct_path: Path | None = None
        self._base: Construct | None = None
        self._circular = False
        self._region: Region | None = None
        self._grouping = Grouping()
        self._pattern = config.group_pattern
        self._reports: list[CloneReport] = []
        self._shown_construct: Construct | None = None
        self._task = None
        self._generation = 0
        self._loading = False
        self._project_path: Path | None = None
        self._selected: tuple[str, int | None] | None = None   # (clone, ref_pos)

        self.empty_state = EmptyState()
        self.summary = SummaryBar()
        self.controls = ControlsBar()
        self.matrix = CloneMatrix()
        self.detail = DetailPane()
        self.detail.hide()
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setWidget(self.matrix)
        splitter = QSplitter(Qt.Orientation.Vertical)
        splitter.addWidget(scroll)
        splitter.addWidget(self.detail)
        splitter.setStretchFactor(0, 2)
        splitter.setStretchFactor(1, 3)
        main = QWidget()
        box = QVBoxLayout(main)
        box.setContentsMargins(0, 0, 0, 0)
        box.setSpacing(2)
        top = QHBoxLayout()
        top.addWidget(self.summary, 1)
        top.addWidget(self.controls.data_buttons)
        box.addLayout(top)
        box.addWidget(self.controls)
        box.addWidget(splitter, 1)
        self._main_page = main
        self._pages = QStackedWidget()
        self._pages.addWidget(self.empty_state)
        self._pages.addWidget(main)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._pages)
        self.setAcceptDrops(True)
        self._connect()
        self.set_theme(status.current_theme())
        self._refresh_zones()

    def _connect(self) -> None:
        empty = self.empty_state
        empty.construct_zone.clicked.connect(self.choose_construct)
        empty.reads_zone.clicked.connect(self.choose_reads_folder)
        empty.reads_files.clicked.connect(self.choose_read_files)
        for zone in empty.zones():
            zone.dropped.connect(self.handle_drop)
        c = self.controls
        c.region_edited.connect(self.set_region)
        c.whole_requested.connect(lambda: self.set_region(None))
        c.circular_toggled.connect(self.set_circular)
        c.params_changed.connect(lambda: self.set_params(c.params()))
        c.grouping_button.clicked.connect(self.edit_grouping)
        c.reads_button.clicked.connect(self.choose_reads_folder)
        c.construct_button.clicked.connect(self.choose_construct)
        self.summary.filter_changed.connect(self.matrix.set_filter)
        self.matrix.region_changed.connect(self.set_region)
        self.matrix.difference_clicked.connect(self._show_difference)
        self.matrix.clone_clicked.connect(lambda clone: self._show_difference(clone, None))
        self.detail.open_read.connect(self.open_read)

    # ------------------------------------------------------------- state API
    def construct(self) -> Construct | None:
        return self._base.with_circular(self._circular) if self._base else None

    def region(self) -> Region | None:
        return self._region

    def grouping(self) -> Grouping:
        return Grouping(list(self._grouping.reads), list(self._grouping.unassigned))

    def params(self) -> VerifyParams:
        return self.controls.params()

    def reports(self) -> list[CloneReport]:
        return self._reports

    def is_verifying(self) -> bool:
        return self._task is not None

    def project_path(self) -> Path | None:
        return self._project_path

    def has_data(self) -> bool:
        return self._base is not None or bool(self._grouping.reads or self._grouping.unassigned)

    def contains_read(self, path: Path | str) -> bool:
        key = path_key(path)
        return any(path_key(p) == key for p in self._all_paths())

    # -------------------------------------------------------------- changes
    def set_construct(self, path: Path | str) -> bool:
        path = Path(path)
        try:
            construct = load_construct(path)
        except (ValueError, OSError) as exc:
            self._warn(tr("workspace.construct_failed"), f"{path}:\n{exc}")
            return False
        if len(construct.seq) > MAX_REFERENCE_LENGTH:
            self._warn(tr("workspace.construct_failed"),
                       tr("workspace.construct_too_long", n=len(construct.seq),
                          limit=MAX_REFERENCE_LENGTH))
            return False
        self._construct_path, self._base, self._circular = path, construct, construct.circular
        if not self._region_ok(self._region, construct.circular):
            self._region = None
        self._sync_controls()
        self._reverify()
        return True

    def add_reads(self, paths: Iterable[Path | str]) -> int:
        known = {path_key(p) for p in self._all_paths()}
        new = [p for p in ab1_files(paths) if path_key(p) not in known]
        if not new:
            self.status_message.emit(tr("workspace.no_new_reads"), 5000)
            return 0
        try:
            added = group_reads(new, self._pattern or None, self.config.primer_tokens)
        except ValueError as exc:   # stale custom pattern in the config
            self._warn(tr("grouping.title"), str(exc))
            added = group_reads(new, None, self.config.primer_tokens)
        self._grouping.reads += added.reads
        self._grouping.unassigned += added.unassigned
        self._reverify()
        return len(new)

    def set_grouping(self, grouping: Grouping) -> None:
        self._grouping = Grouping(list(grouping.reads), list(grouping.unassigned))
        self._reverify()

    def set_region(self, region: Region | None) -> bool:
        if region is not None:
            region = (int(region[0]), int(region[1]))
        if not self._region_ok(region, self._circular, warn=True):
            self._sync_controls()
            return False
        if region == self._region:
            self._sync_controls()
            return True
        self._region = region
        self._sync_controls()
        self._reverify()
        return True

    def set_circular(self, flag: bool) -> bool:
        if not self._region_ok(self._region, flag, warn=True):
            self._sync_controls()
            return False
        self._circular = bool(flag)
        self._sync_controls()
        self._reverify()
        return True

    def set_params(self, params: VerifyParams) -> None:
        self.controls.set_params(params)
        self._reverify()

    def read_edited(self, path: Path | str) -> bool:
        """A read was edited elsewhere: re-run the check if it is part of this session."""
        if not self.contains_read(path):
            return False
        self._reverify()
        return True

    def clear(self) -> None:
        self._cancel_run()
        self._construct_path = self._base = None
        self._circular, self._region = False, None
        self._grouping, self._reports = Grouping(), []
        self._project_path = None
        self._selected = None
        self.detail.hide()
        self.controls.set_params(VerifyParams())
        self._sync_controls()
        self._refresh_zones()

    # --------------------------------------------------------------- project
    def load_project(self, path: Path | str) -> bool:
        path = Path(path)
        try:
            project = project_io.load_project(path)
        except ValueError as exc:
            self.clear()
            self._warn(tr("workspace.project_failed"), str(exc))
            return False
        self.clear()
        self._loading = True
        try:
            ok = self.set_construct(project.construct_path)
            self._grouping = project.grouping
            self.controls.set_params(project.params)
            if ok:
                self._circular = project.circular
                if self._region_ok(project.region, project.circular, warn=True):
                    self._region = project.region
            self._project_path = path
        finally:
            self._loading = False
        self._sync_controls()
        self._reverify()
        self.status_message.emit(tr("workspace.project_loaded", name=path.name), 5000)
        return ok

    def save_project(self, path: Path | str) -> bool:
        path = Path(path)
        if self._construct_path is None:
            self._warn(tr("workspace.project_save_failed"), tr("workspace.need_construct"))
            return False
        project = project_io.Project(self._construct_path, self._circular, self._region,
                                     self.grouping(), self.params())
        try:
            project_io.save_project(path, project)
        except OSError as exc:
            self._warn(tr("workspace.project_save_failed"), f"{path}:\n{exc}")
            return False
        self._project_path = path
        self.status_message.emit(tr("workspace.project_saved", name=path.name), 5000)
        return True

    def export_report(self, path: Path | str) -> bool:
        path = Path(path)
        if not self._reports or self._shown_construct is None:
            self._warn(tr("workspace.export_failed"), tr("workspace.no_results"))
            return False
        try:
            write_report(path, self._reports, self._shown_construct.name)
        except (OSError, RuntimeError) as exc:
            self._warn(tr("workspace.export_failed"), f"{path}:\n{exc}")
            return False
        self.status_message.emit(tr("workspace.exported", name=path.name), 5000)
        return True

    # --------------------------------------------------------------- dialogs
    def choose_construct(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, tr("workspace.construct_button"), self._start_dir(), CONSTRUCT_FILTER)
        if path:
            self.set_construct(path)

    def choose_reads_folder(self) -> None:
        folder = QFileDialog.getExistingDirectory(self, tr("workspace.add_reads"),
                                                  self._start_dir())
        if folder:
            self.add_reads([folder])

    def choose_read_files(self) -> None:
        paths, _ = QFileDialog.getOpenFileNames(
            self, tr("workspace.add_reads"), self._start_dir(),
            "Sanger chromatograms (*.ab1);;All files (*)")
        if paths:
            self.add_reads(paths)

    def edit_grouping(self) -> None:
        dialog = GroupingDialog(self._grouping, self._pattern, self, self.config.primer_tokens)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        if dialog.result_pattern() != self._pattern:
            self._pattern = self.config.group_pattern = dialog.result_pattern()
            save_config(self.config)
        self.set_grouping(dialog.result_grouping())

    def handle_drop(self, paths: list[Path]) -> None:
        """Projects open, ``.ab1`` files and folders become reads, anything else
        is taken as the construct."""
        paths = [Path(p) for p in paths]
        project = next((p for p in paths if p.suffix.lower() == PROJECT_SUFFIX), None)
        if project is not None:
            self.load_project(project)
            return
        constructs = [p for p in paths if p.is_file() and p.suffix.lower() != ".ab1"]
        reads = [p for p in paths if p.is_dir() or p.suffix.lower() == ".ab1"]
        self._loading = True
        try:
            if constructs:
                self.set_construct(constructs[0])
            if reads:
                self.add_reads(reads)
        finally:
            self._loading = False
        self._reverify()

    def dragEnterEvent(self, event: QDragEnterEvent) -> None:  # noqa: N802
        if dropped_paths(event):
            event.acceptProposedAction()

    def dropEvent(self, event: QDropEvent) -> None:  # noqa: N802
        paths = dropped_paths(event)
        if paths:
            event.acceptProposedAction()
            self.handle_drop(paths)

    def set_theme(self, name: str) -> None:
        self.summary.set_theme(name)
        self.matrix.set_theme(name)
        self.detail.set_theme(name)
        self.empty_state.set_theme(name)

    # ------------------------------------------------------------ verification
    def _reverify(self) -> None:
        if self._loading:
            return
        self._refresh_zones()
        construct = self.construct()
        if construct is None:
            return
        self._cancel_run()
        if not self._grouping.reads:
            self._show([], construct, self._region)
            return
        region, grouping, params = self._region, self.grouping(), self.params()
        generation = self._generation

        def work(progress, cancelled):
            return verify_all(construct, region, grouping, params, progress, cancelled)

        def done(reports):
            if generation == self._generation:
                self._task = None
                self._show(reports, construct, region)

        task = self._start_task(work, done, tr("workspace.verify_failed"), key="verify")
        self._task = task
        task.cancelled.connect(lambda: self._task_ended(task))
        task.failed.connect(lambda _m: self._task_ended(task))

    def _task_ended(self, task) -> None:
        if self._task is task:
            self._task = None

    def _cancel_run(self) -> None:
        self._generation += 1
        task, self._task = self._task, None
        if task is not None:
            task.cancel()

    def _show(self, reports: list[CloneReport], construct: Construct,
              region: Region | None) -> None:
        self._reports, self._shown_construct = reports, construct
        self.summary.set_data(construct, region, reports)
        self.matrix.set_data(construct, region, reports)
        self.matrix.set_filter(self.summary.shown_verdicts())
        if self._selected is not None:
            clone, pos = self._selected
            if any(r.clone == clone for r in reports):
                self._show_difference(clone, pos)
            else:
                self._selected = None
                self.detail.hide()
        if reports:
            self.status_message.emit(tr("workspace.verified", n=len(reports)), 5000)

    def _show_difference(self, clone: str, ref_pos: int | None) -> None:
        report = next((r for r in self._reports if r.clone == clone), None)
        if report is None or self._shown_construct is None:
            return
        self._selected = (clone, ref_pos)
        self.detail.show()
        self.detail.show_clone(report, self._shown_construct)
        if ref_pos is not None:
            self.detail.show_position(ref_pos)

    # --------------------------------------------------------------- helpers
    def _region_ok(self, region: Region | None, circular: bool, warn: bool = False) -> bool:
        if region is None or self._base is None:
            return True
        try:
            check_region(region, len(self._base.seq), circular)
        except ValueError as exc:
            if warn:
                start, end = region
                self._warn(tr("workspace.region_invalid"),
                           tr("workspace.region_invalid_text", start=start + 1, end=end,
                              reason=exc))
            return False
        return True

    def _sync_controls(self) -> None:
        length = len(self._base.seq) if self._base else 1
        self.controls.set_region(self._region, length)
        self.controls.set_circular(self._circular)

    def _refresh_zones(self) -> None:
        empty = self.empty_state
        empty.construct_zone.set_done(
            f"{self._base.name} · {self._construct_path.name}" if self._base else "")
        n = len(self._grouping.reads) + len(self._grouping.unassigned)
        empty.reads_zone.set_done(tr("workspace.reads_loaded", n=n) if n else "")
        ready = self._base is not None and n > 0
        self._pages.setCurrentWidget(self._main_page if ready else self.empty_state)

    def _all_paths(self) -> list[Path]:
        return [r.path for r in self._grouping.reads] + list(self._grouping.unassigned)

    def _start_dir(self) -> str:
        if self._project_path is not None:
            return str(self._project_path.parent)
        if self._construct_path is not None:
            return str(self._construct_path.parent)
        return ""

    def _warn(self, title: str, text: str) -> None:
        QMessageBox.warning(self, title, text)
