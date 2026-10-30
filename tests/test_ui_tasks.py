"""Background task runner tests (offscreen, no pytest-qt)."""

import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import time
from concurrent.futures import CancelledError

from pathlib import Path

import pytest

AB1 = sorted((Path(__file__).resolve().parent.parent / "examples").glob("*.ab1"))


@pytest.fixture(scope="module")
def qapp():
    from PySide6.QtWidgets import QApplication

    return QApplication.instance() or QApplication([])


@pytest.fixture()
def parent(qapp):
    from PySide6.QtCore import QObject

    return QObject()


def _wait(qapp, predicate, timeout=5.0):
    deadline = time.monotonic() + timeout
    while not predicate() and time.monotonic() < deadline:
        qapp.processEvents()
        time.sleep(0.005)
    return predicate()


def _collect(task):
    events = []
    task.finished.connect(lambda v: events.append(("finished", v)))
    task.failed.connect(lambda m: events.append(("failed", m)))
    task.cancelled.connect(lambda: events.append(("cancelled", None)))
    task.progress.connect(lambda d, t: events.append(("progress", (d, t))))
    return events


def test_result_delivered_in_finished(qapp, parent):
    from ui.tasks import run_task

    task = run_task(lambda progress, is_cancelled: (progress(1, 2), 42)[1], parent)
    events = _collect(task)
    assert _wait(qapp, lambda: any(e[0] == "finished" for e in events))
    assert ("finished", 42) in events
    assert ("progress", (1, 2)) in events


def test_exception_goes_to_failed(qapp, parent):
    from ui.tasks import run_task

    def boom(progress, is_cancelled):
        raise ValueError("bad thing")

    task = run_task(boom, parent)
    events = _collect(task)
    assert _wait(qapp, lambda: bool(events))
    assert events == [("failed", "bad thing")]


def test_cancel_emits_cancelled_and_drops_result(qapp, parent):
    from ui.tasks import run_task

    def loop(progress, is_cancelled):
        for _ in range(500):
            if is_cancelled():
                raise CancelledError()
            time.sleep(0.01)
        return "late"

    task = run_task(loop, parent)
    events = _collect(task)
    task.cancel()
    assert task.is_cancelled()
    assert events == [("cancelled", None)]  # emitted immediately
    assert _wait(qapp, lambda: bool(events))
    assert events == [("cancelled", None)]


def test_cancel_drops_late_result_of_uninterruptible_fn(qapp, parent):
    from ui.tasks import run_task

    task = run_task(lambda p, c: (time.sleep(0.2), "late")[1], parent)
    events = _collect(task)
    task.cancel()
    assert _wait(qapp, lambda: bool(events))
    assert events == [("cancelled", None)]


def test_gui_thread_not_blocked(qapp, parent):
    from PySide6.QtCore import QTimer

    from ui.tasks import run_task

    ticks = []
    timer = QTimer()
    timer.timeout.connect(lambda: ticks.append(1))
    timer.start(20)
    task = run_task(lambda p, c: (time.sleep(0.3), "ok")[1], parent)
    events = _collect(task)
    assert _wait(qapp, lambda: bool(events))
    timer.stop()
    assert events == [("finished", "ok")]
    assert len(ticks) >= 5


def test_signals_delivered_in_gui_thread(qapp, parent):
    from PySide6.QtCore import QThread

    from ui.tasks import run_task

    threads = []
    task = run_task(lambda p, c: 1, parent)
    task.finished.connect(lambda v: threads.append(QThread.currentThread()))
    assert _wait(qapp, lambda: bool(threads))
    assert threads[0] is qapp.thread()


@pytest.mark.skipif(not AB1, reason="no example .ab1 files")
def test_main_window_blast_uses_task_and_cancel_drops_result(qapp, monkeypatch, tmp_path):
    from chrom_core import blast
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    monkeypatch.setattr("ui.main_window.save_config", lambda *a, **k: None)
    monkeypatch.setattr("ui.read_view.read_view.save_config", lambda *a, **k: None)
    opened = []
    monkeypatch.setattr(
        "ui.read_view.read_view.QDesktopServices.openUrl", lambda url: opened.append(url)
    )
    monkeypatch.setattr(blast, "submit", lambda params, timeout=20: (time.sleep(0.2), "RID1")[1])
    window = MainWindow(AppConfig(recent_files=[]))
    window.show()
    assert window.open_path(AB1[0])
    window._blast_submit()
    assert window.task_cancel_button.isVisible()
    assert window.task_progress.isVisible()
    assert _wait(qapp, lambda: bool(opened))
    assert _wait(qapp, lambda: not window.task_progress.isVisible())

    opened.clear()
    window._blast_submit()
    window.task_cancel_button.click()
    assert _wait(qapp, lambda: not window.task_progress.isVisible())
    time.sleep(0.3)
    qapp.processEvents()
    assert opened == []
    window.close()


def test_base_exception_still_terminal(qapp, parent):
    from ui.tasks import run_task

    def bail(progress, is_cancelled):
        raise KeyboardInterrupt()

    task = run_task(bail, parent)
    events = _collect(task)
    assert _wait(qapp, lambda: bool(events))
    assert events == [("failed", "KeyboardInterrupt")]


def test_parent_deleted_mid_run_is_clean(qapp, capfd):
    import shiboken6
    from PySide6.QtCore import QObject

    from ui import tasks
    from ui.tasks import run_task

    def slow(progress, is_cancelled):
        time.sleep(0.2)
        progress(1, 2)
        raise RuntimeError("after delete")

    parent = QObject()
    run_task(slow, parent)
    shiboken6.delete(parent)
    assert _wait(qapp, lambda: not tasks._running)
    time.sleep(0.05)
    assert "Traceback" not in capfd.readouterr().err
    assert not tasks._running


@pytest.mark.skipif(not AB1, reason="no example .ab1 files")
def test_scan_last_start_wins(qapp, monkeypatch, tmp_path):
    from chrom_core import batch
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    monkeypatch.setattr("ui.main_window.save_config", lambda *a, **k: None)
    monkeypatch.setattr("ui.read_view.read_view.save_config", lambda *a, **k: None)
    real = batch.scan_folder

    def fake(folder, *a, **k):
        if folder.name == "slow":
            time.sleep(0.4)
            return real(folder.parent, *a, **k)[:1]
        return real(folder.parent, *a, **k)[:2]

    monkeypatch.setattr(batch, "scan_folder", fake)
    window = MainWindow(AppConfig(recent_files=[]))
    window.show()
    ex = AB1[0].parent
    window._scan_folder(ex / "slow")
    window._scan_folder(ex / "fast")
    assert _wait(qapp, lambda: window.batch_panel.table.rowCount() == min(2, len(AB1)))
    time.sleep(0.6)
    qapp.processEvents()
    assert window.batch_panel.table.rowCount() == min(2, len(AB1))
    assert not window.task_progress.isVisible()
    window.close()


@pytest.mark.skipif(not AB1, reason="no example .ab1 files")
def test_second_blast_cancels_first(qapp, monkeypatch):
    from chrom_core import blast
    from chrom_core.config import AppConfig
    from ui.main_window import MainWindow

    monkeypatch.setattr("ui.main_window.save_config", lambda *a, **k: None)
    monkeypatch.setattr("ui.read_view.read_view.save_config", lambda *a, **k: None)
    opened = []
    monkeypatch.setattr(
        "ui.read_view.read_view.QDesktopServices.openUrl", lambda url: opened.append(url.toString())
    )
    rids = iter(["FIRST", "SECOND"])
    monkeypatch.setattr(blast, "submit", lambda params, timeout=20: (time.sleep(0.2), next(rids))[1])
    window = MainWindow(AppConfig(recent_files=[]))
    assert window.open_path(AB1[0])
    window._blast_submit()
    window._blast_submit()
    assert _wait(qapp, lambda: bool(opened))
    time.sleep(0.4)
    qapp.processEvents()
    assert len(opened) == 1
    window.close()
