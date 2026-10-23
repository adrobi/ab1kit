"""Main ab1kit window: toolbar, menus, overview map, chromatogram, docks."""

from __future__ import annotations

import sys
from importlib.metadata import PackageNotFoundError, version as pkg_version
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QDragEnterEvent, QDropEvent
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
    QMessageBox,
    QPushButton,
    QScrollBar,
    QToolBar,
    QVBoxLayout,
    QWidget,
)

from chrom_core import Chromatogram, read_ab1
from chrom_core.config import AppConfig, load_config, save_config
from chrom_core.seqops import find_matches
from i18n import set_language, tr
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


class MainWindow(QMainWindow):
    def __init__(self, config: AppConfig | None = None) -> None:
        super().__init__()
        self.config = config or load_config()
        set_language(self.config.language)

        self.setWindowTitle("ab1kit")
        self.resize(1200, 780)
        self._path: Path | None = None
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
        open_action = QAction(tr("menu.file.open"), self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self._open_dialog)
        file_menu.addAction(open_action)
        self.recent_menu = file_menu.addMenu(tr("menu.file.recent"))
        file_menu.addSeparator()
        exit_action = QAction(tr("menu.file.exit"), self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        edit_menu = self.menuBar().addMenu(tr("menu.edit"))
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

    def _connect(self) -> None:
        self.chrom_view.view_changed.connect(self._sync_from_view)
        self.chrom_view.cursor_info.connect(self.statusBar().showMessage)
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
        self.chrom_view.set_chromatogram(chrom)
        self.minimap.set_chromatogram(chrom)
        self._search_ranges = []
        self._search_pos = -1
        self.setWindowTitle(f"ab1kit — {self._path.name}")
        self._push_recent(str(self._path))
        self._refresh_info(chrom)
        n = len(chrom.seq)
        mean_q = sum(chrom.quality) / max(1, len(chrom.quality))
        self.statusBar().showMessage(tr("status.file", n=n, q=mean_q), 5000)
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


def run(argv: list[str] | None = None) -> int:
    """Application entry point used by the ``ab1kit`` console script."""
    args = list(sys.argv[1:] if argv is None else argv)
    app = QApplication(args)
    window = MainWindow()
    window.show()
    for arg in args:
        if not arg.startswith("-") and Path(arg).exists():
            window.open_path(arg)
            break
    return app.exec()
