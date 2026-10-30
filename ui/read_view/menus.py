"""Read-specific menu items, appended by the shell into its own menus.

Shortcuts, order and enabled-state handling are the ones the single window
used to have; the actions are parented to the view so they live with it.
"""

from __future__ import annotations

from PySide6.QtGui import QAction
from PySide6.QtWidgets import QMenu

from i18n import tr


def add_file_actions(view, file_menu: QMenu) -> None:
    """Save session, export submenu, copy sequence (between 'Recent' and 'Exit')."""
    save_session = QAction(tr("menu.file.save_session"), view)
    save_session.setShortcut("Ctrl+S")
    save_session.triggered.connect(view._save_session)
    file_menu.addAction(save_session)
    export_menu = file_menu.addMenu(tr("menu.file.export"))
    for text, slot in (
        (tr("menu.file.export_fasta"), view._export_fasta),
        (tr("menu.file.export_fastq"), view._export_fastq),
        (tr("menu.file.export_png"), view._export_png),
        (tr("menu.file.export_pdf"), view._export_pdf),
    ):
        action = QAction(text, view)
        action.triggered.connect(slot)
        export_menu.addAction(action)
    file_menu.addSeparator()
    copy_action = QAction(tr("menu.file.copy_seq"), view)
    copy_action.setShortcut("Ctrl+Shift+C")
    copy_action.triggered.connect(view._copy_sequence)
    file_menu.addAction(copy_action)


def add_edit_actions(view, edit_menu: QMenu) -> None:
    edit_menu.addAction(view.undo_action)
    edit_menu.addAction(view.redo_action)
    edit_menu.addSeparator()
    find_action = QAction(tr("menu.edit.find"), view)
    find_action.setShortcut("Ctrl+F")
    find_action.triggered.connect(view._search)
    edit_menu.addAction(find_action)
    find_next = QAction(tr("menu.edit.find_next"), view)
    find_next.setShortcut("F3")
    find_next.triggered.connect(view._find_next)
    edit_menu.addAction(find_next)
    goto_action = QAction(tr("menu.edit.goto"), view)
    goto_action.setShortcut("Ctrl+G")
    goto_action.triggered.connect(view._goto)
    edit_menu.addAction(goto_action)
    edit_menu.addSeparator()
    revcomp_action = QAction(tr("menu.edit.revcomp"), view)
    revcomp_action.setShortcut("Ctrl+R")
    revcomp_action.triggered.connect(view.reverse_complement_read)
    edit_menu.addAction(revcomp_action)
    trim_action = QAction(tr("menu.edit.trim"), view)
    trim_action.setShortcut("Ctrl+T")
    trim_action.triggered.connect(view._open_trim_dialog)
    edit_menu.addAction(trim_action)


def add_view_actions(view, view_menu: QMenu) -> None:
    """Smoothing, quality profile, info and batch docks (after 'Dark theme')."""
    smooth_action = QAction(tr("menu.view.smoothing"), view)
    smooth_action.setCheckable(True)
    smooth_action.toggled.connect(view.chrom_view.set_smoothing)
    view_menu.addAction(smooth_action)
    profile_action = QAction(tr("menu.view.quality_profile"), view)
    profile_action.setCheckable(True)
    profile_action.toggled.connect(view.minimap.set_profile)
    view_menu.addAction(profile_action)
    info_action = QAction(tr("menu.view.file_info"), view)
    info_action.setCheckable(True)
    info_action.toggled.connect(view.info_dock.setVisible)
    view.info_dock.visibilityChanged.connect(info_action.setChecked)
    view_menu.addAction(info_action)
    batch_action = QAction(tr("menu.view.batch_panel"), view)
    batch_action.setCheckable(True)
    batch_action.toggled.connect(view.batch_dock.setVisible)
    view.batch_dock.visibilityChanged.connect(batch_action.setChecked)
    view_menu.addAction(batch_action)


def add_analysis_actions(view, analysis_menu: QMenu) -> None:
    reference_action = QAction(tr("menu.analysis.reference"), view)
    reference_action.setShortcut("Ctrl+L")
    reference_action.triggered.connect(view._reference_dialog)
    analysis_menu.addAction(reference_action)
    clear_ref = QAction(tr("menu.analysis.clear_reference"), view)
    clear_ref.triggered.connect(view._clear_reference)
    analysis_menu.addAction(clear_ref)
    analysis_menu.addSeparator()
    export_aligned = QAction(tr("menu.analysis.export_aligned"), view)
    export_aligned.triggered.connect(view._export_aligned)
    analysis_menu.addAction(export_aligned)
    analysis_menu.addSeparator()
    verify_action = QAction(tr("menu.analysis.verify"), view)
    verify_action.setShortcut("Ctrl+Shift+V")
    verify_action.triggered.connect(view._verify_dialog)
    analysis_menu.addAction(verify_action)
    analysis_menu.addSeparator()
    blast_action = QAction(tr("menu.analysis.blast"), view)
    blast_action.triggered.connect(view._blast_submit)
    analysis_menu.addAction(blast_action)
