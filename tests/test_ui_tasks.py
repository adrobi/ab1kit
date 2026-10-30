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
    from PySide6.QtCore import QObject

    from ui.tasks import run_task

    task = run_task(lambda progress, is_cancelled: (progress(1, 2), 42)[1], parent)
    events = _collect(task)
    assert _wait(qapp, lambda: any(e[0] == "finished" for e in events))
    assert ("finished", 42) in events
    assert ("progress", (1, 2)) in events


def test_exception_goes_to_failed(qapp, parent):
    from PySide6.QtCore import QObject

    from ui.tasks import run_task

    def boom(progress, is_cancelled):
        raise ValueError("bad thing")

    task = run_task(boom, parent)
    events = _collect(task)
    assert _wait(qapp, lambda: bool(events))
    assert events == [("failed", "bad thing")]


def test_cancel_emits_cancelled_and_drops_result(qapp, parent):
    from PySide6.QtCore import QObject

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
    assert _wait(qapp, lambda: bool(events))
    assert events == [("cancelled", None)]


def test_cancel_drops_late_result_of_uninterruptible_fn(qapp, parent):
    from PySide6.QtCore import QObject

    from ui.tasks import run_task

    task = run_task(lambda p, c: (time.sleep(0.2), "late")[1], parent)
    events = _collect(task)
    task.cancel()
    assert _wait(qapp, lambda: bool(events))
    assert events == [("cancelled", None)]


def test_gui_thread_not_blocked(qapp, parent):
    from PySide6.QtCore import QObject, QTimer

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
    from PySide6.QtCore import QObject, QThread

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
    opened = []
    monkeypatch.setattr(
        "ui.main_window.QDesktopServices.openUrl", lambda url: opened.append(url)
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
