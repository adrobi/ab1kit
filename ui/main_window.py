"""ab1kit shell window: tabs "construct check" / "read viewer", File/Project/
View/Help menus, language, theme, recent files, background-task status bar.
The single-read viewer lives in ``ui.read_view``, the construct check in
``ui.construct.workspace``."""

from __future__ import annotations

import sys
from importlib.metadata import PackageNotFoundError, version as pkg_version
from pathlib import Path

from PySide6.QtGui import QAction, QDragEnterEvent, QDropEvent
from PySide6.QtWidgets import (
    QApplication,
    QFileDialog,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QTabWidget,
)

from chrom_core.config import AppConfig, load_config, save_config
from i18n import set_language, tr
from ui import status
from ui.construct.workspace import PROJECT_SUFFIX, ConstructWorkspace
from ui.read_view import ReadView
from ui.status import bind_icon
from ui.tasks import Task, run_task
from ui.theme import app_qss

try:
    VERSION = pkg_version("ab1kit")
except PackageNotFoundError:  # running from source without install
    VERSION = "dev"


def _delegate(name: str) -> property:
    """Read-only view onto ``ReadView.<name>`` (keeps the old MainWindow API)."""
    return property(lambda self: getattr(self.read_view, name))


class MainWindow(QMainWindow):
    # State and widgets that used to live on the window itself.
    read_model = _delegate("read_model")
    chrom_view = _delegate("chrom_view")
    minimap = _delegate("minimap")
    trim_dialog = _delegate("trim_dialog")
    reference_seq = _delegate("reference_seq")
    alignment = _delegate("alignment")
    snp_calls = _delegate("snp_calls")
    variant_verdicts = _delegate("variant_verdicts")
    batch_reports = _delegate("batch_reports")
    batch_panel = _delegate("batch_panel")
    batch_dock = _delegate("batch_dock")
    diff_dock = _delegate("diff_dock")
    diff_table = _delegate("diff_table")
    verify_dock = _delegate("verify_dock")
    verify_table = _delegate("verify_table")
    undo_action = _delegate("undo_action")
    redo_action = _delegate("redo_action")

    def __init__(self, config: AppConfig | None = None) -> None:
        super().__init__()
        self.config = config or load_config()
        set_language(self.config.language)

        self.setWindowTitle("ab1kit")
        self.resize(1200, 780)
        self.read_view = ReadView(self.config, self._start_task)
        self.workspace = ConstructWorkspace(self.config, self._start_task)
        self.tabs = QTabWidget()
        self.tabs.setDocumentMode(True)
        self.tabs.addTab(self.workspace, tr("tab.construct"))
        self.tabs.addTab(self.read_view, tr("tab.read"))
        self.setCentralWidget(self.tabs)
        self.read_view.file_opened.connect(self._push_recent)
        self.read_view.file_opened.connect(lambda _p: self.tabs.setCurrentWidget(self.read_view))
        self.read_view.status_message.connect(self.statusBar().showMessage)
        self.read_view.title_changed.connect(self.setWindowTitle)
        self.read_view.read_edited.connect(self._read_edited)
        self.workspace.status_message.connect(self.statusBar().showMessage)
        self.workspace.open_read.connect(self._open_read)

        self._build_menus()
        self.setAcceptDrops(True)
        self.set_theme(self.config.theme)
        self.statusBar().showMessage(tr("status.ready"))
        self._build_task_status()
        self._fill_recent_menu()

    # ------------------------------------------------------------------ tasks
    def _build_task_status(self) -> None:
        self._active_tasks: set[Task] = set()
        self._keyed_tasks: dict[str, Task] = {}
        self.task_progress = QProgressBar()
        self.task_progress.setMaximumWidth(180)
        self.task_progress.setFormat(tr("task.running"))
        self.task_cancel_button = QPushButton(tr("task.cancel"))
        self.task_cancel_button.clicked.connect(self._cancel_tasks)
        for widget in (self.task_progress, self.task_cancel_button):
            self.statusBar().addPermanentWidget(widget)
            widget.hide()

    def _cancel_tasks(self) -> None:
        for task in list(self._active_tasks):
            task.cancel()

    def _start_task(self, fn, on_finished, error_title: str, key: str = "") -> Task:
        """Run ``fn`` in the background; ``on_finished(result)`` runs in the GUI
        thread unless the task was cancelled or failed. A new task with the
        same ``key`` cancels the previous one (last start wins)."""
        if key and key in self._keyed_tasks:
            self._keyed_tasks.pop(key).cancel()
        task = run_task(fn, self)
        if key:
            self._keyed_tasks[key] = task
        self._active_tasks.add(task)
        self.task_progress.setRange(0, 0)  # busy until progress arrives
        self.task_progress.show()
        self.task_cancel_button.show()

        def release() -> None:
            self._active_tasks.discard(task)
            if key and self._keyed_tasks.get(key) is task:
                del self._keyed_tasks[key]
            if not self._active_tasks:
                self.task_progress.hide()
                self.task_cancel_button.hide()

        def on_progress(done: int, total: int) -> None:
            if total > 0:
                self.task_progress.setRange(0, total)
                self.task_progress.setValue(done)

        def on_failed(message: str) -> None:
            release()
            QMessageBox.warning(self, error_title, message)

        def on_cancelled() -> None:
            release()
            self.statusBar().showMessage(tr("task.cancelled"), 5000)

        def on_done(result: object) -> None:
            release()
            on_finished(result)

        task.progress.connect(on_progress)
        task.failed.connect(on_failed)
        task.cancelled.connect(on_cancelled)
        task.finished.connect(on_done)
        return task

    # ------------------------------------------------------------------ menus
    def _build_menus(self) -> None:
        file_menu = self.menuBar().addMenu(tr("menu.file"))
        self._file_menu = file_menu
        open_action = QAction(tr("menu.file.open"), self)
        bind_icon(open_action, "folder-open")
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(lambda: self._in_read_tab(self.read_view.open_dialog))
        file_menu.addAction(open_action)
        open_folder = QAction(tr("menu.file.open_folder"), self)
        open_folder.setShortcut("Ctrl+Shift+O")
        open_folder.triggered.connect(
            lambda: self._in_read_tab(self.read_view.open_folder_dialog))
        file_menu.addAction(open_folder)
        self.recent_menu = file_menu.addMenu(tr("menu.file.recent"))
        self.read_view.add_file_actions(file_menu)
        file_menu.addSeparator()
        exit_action = QAction(tr("menu.file.exit"), self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        self._build_project_menu(self.menuBar().addMenu(tr("menu.project")))
        self.read_view.add_edit_menu(self.menuBar().addMenu(tr("menu.edit")))

        view_menu = self.menuBar().addMenu(tr("menu.view"))
        self.dark_action = QAction(tr("menu.view.dark_theme"), self)
        bind_icon(self.dark_action, "moon")
        self.dark_action.setCheckable(True)
        self.dark_action.setChecked(self.config.theme == "dark")
        self.dark_action.toggled.connect(
            lambda on: self.set_theme("dark" if on else "light")
        )
        view_menu.addAction(self.dark_action)
        self.read_view.add_view_actions(view_menu)
        lang_menu = view_menu.addMenu(tr("menu.view.language"))
        self._lang_actions = {}
        for code, label in (("ru", tr("lang.ru")), ("en", tr("lang.en"))):
            action = QAction(label, self)
            action.setCheckable(True)
            action.setChecked(self.config.language == code)
            action.triggered.connect(lambda _=False, c=code: self._switch_language(c))
            lang_menu.addAction(action)
            self._lang_actions[code] = action

        self.read_view.add_analysis_menu(self.menuBar().addMenu(tr("menu.analysis")))

        help_menu = self.menuBar().addMenu(tr("menu.help"))
        about_action = QAction(tr("menu.help.about"), self)
        about_action.triggered.connect(
            lambda: QMessageBox.about(self, "ab1kit", tr("about.text", version=VERSION))
        )
        help_menu.addAction(about_action)

    def _build_project_menu(self, menu) -> None:
        self.project_menu = menu
        for key, shortcut, slot in (
            ("menu.project.new", "", self._project_new),
            ("menu.project.open", "Ctrl+Shift+P", self._project_open),
            ("menu.project.save", "Ctrl+Shift+S", self._project_save),
            ("menu.project.save_as", "", self._project_save_as),
            (None, "", None),
            ("menu.project.export", "Ctrl+Shift+E", self._project_export),
        ):
            if key is None:
                menu.addSeparator()
                continue
            action = QAction(tr(key), self)
            if shortcut:
                action.setShortcut(shortcut)
            action.triggered.connect(slot)
            menu.addAction(action)

    # --------------------------------------------------------------- projects
    def _show_workspace(self) -> None:
        self.tabs.setCurrentWidget(self.workspace)

    def _project_new(self) -> None:
        self.workspace.clear()
        self._show_workspace()

    def _project_open(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, tr("menu.project.open"), "", f"ab1kit project (*{PROJECT_SUFFIX})")
        if path:
            self.open_path(path)

    def _project_save(self) -> None:
        path = self.workspace.project_path()
        if path is None:
            self._project_save_as()
        elif self.workspace.save_project(path):
            self._push_recent(str(path))

    def _project_save_as(self) -> None:
        self._show_workspace()
        path, _ = QFileDialog.getSaveFileName(
            self, tr("menu.project.save_as"), "", f"ab1kit project (*{PROJECT_SUFFIX})")
        if not path:
            return
        if not path.lower().endswith(PROJECT_SUFFIX):
            path += PROJECT_SUFFIX
        if self.workspace.save_project(path):
            self._push_recent(path)

    def _project_export(self) -> None:
        self._show_workspace()
        path, chosen = QFileDialog.getSaveFileName(
            self, tr("menu.project.export"), "", "Excel (*.xlsx);;CSV (*.csv)")
        if not path:
            return
        if Path(path).suffix.lower() not in (".xlsx", ".csv"):
            path += ".csv" if "csv" in chosen.lower() else ".xlsx"
        self.workspace.export_report(path)

    def _in_read_tab(self, action) -> None:
        self.tabs.setCurrentWidget(self.read_view)
        action()

    def _open_read(self, path: Path) -> None:
        self.read_view.open_path(path)   # file_opened switches to the viewer tab

    def _read_edited(self, path: Path) -> None:
        """Edits reach the check through the read's session sidecar: write it
        now (or drop it after undoing everything), then re-run the check."""
        if self.workspace.contains_read(path):
            self.read_view.autosave_session()
            self.workspace.read_edited(path)

    def _switch_language(self, code: str) -> None:
        self.config.language = code
        save_config(self.config)
        for other, action in self._lang_actions.items():
            action.setChecked(other == code)
        QMessageBox.information(self, tr("menu.view.language"), tr("lang.restart"))

    def set_theme(self, name: str) -> None:
        self.config.theme = name
        status.set_theme(name)
        self.read_view.apply_theme(name)
        self.workspace.set_theme(name)
        QApplication.instance().setStyleSheet(app_qss(name == "dark"))
        save_config(self.config)

    # ----------------------------------------------------------- recent, drops
    def dragEnterEvent(self, event: QDragEnterEvent) -> None:  # noqa: N802
        urls = event.mimeData().urls()
        if urls and urls[0].toLocalFile().lower().endswith((".ab1", PROJECT_SUFFIX)):
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

    # ------------------------------------------------- delegates to the read view
    def open_path(self, path: str | Path) -> bool:
        """A ``.ab1proj`` opens in the construct check, anything else in the viewer."""
        if str(path).lower().endswith(PROJECT_SUFFIX):
            self._show_workspace()
            opened = self.workspace.load_project(path)
            if opened:
                self._push_recent(str(path))
            return opened
        return self.read_view.open_path(path)   # file_opened switches to the viewer tab

    def refresh_read(self) -> None:
        self.read_view.refresh_read()

    def apply_base_change(self, index: int, base: str) -> None:
        self.read_view.apply_base_change(index, base)

    def apply_base_insert(self, index: int, base: str) -> None:
        self.read_view.apply_base_insert(index, base)

    def apply_base_delete(self, index: int) -> None:
        self.read_view.apply_base_delete(index)

    def apply_trim_range(self, start: int, end: int) -> None:
        self.read_view.apply_trim_range(start, end)

    def undo_read(self) -> None:
        self.read_view.undo_read()

    def redo_read(self) -> None:
        self.read_view.redo_read()

    def reverse_complement_read(self) -> None:
        self.read_view.reverse_complement_read()

    def _set_reference_seq(self, sequence: str) -> None:
        self.read_view._set_reference_seq(sequence)

    def _reference_overlay(self) -> tuple[str, dict[int, int]]:
        return self.read_view._reference_overlay()

    def _clear_reference(self) -> None:
        self.read_view._clear_reference()

    def _export_aligned(self) -> None:
        self.read_view._export_aligned()

    def _blast_submit(self) -> None:
        self.read_view._blast_submit()

    def _verify_variants(self, text: str) -> None:
        self.read_view._verify_variants(text)

    def _fill_verify_table(self) -> None:
        self.read_view._fill_verify_table()

    def _save_verify_csv(self) -> None:
        self.read_view._save_verify_csv()

    def _scan_folder(self, folder: Path) -> None:
        self._in_read_tab(lambda: self.read_view._scan_folder(folder))

    def _batch_export(self) -> None:
        self.batch_dock.export_trimmed()

    def _open_trim_dialog(self) -> None:
        self.read_view._open_trim_dialog()

    def closeEvent(self, event) -> None:  # noqa: N802
        self._cancel_tasks()
        self.read_view.autosave_session()
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
