"""Background tasks: run a function in the global thread pool with progress,
cancellation and GUI-thread signal delivery."""

from __future__ import annotations

from collections.abc import Callable
from concurrent.futures import CancelledError
import threading

from PySide6.QtCore import QObject, QRunnable, QThreadPool, Signal, Slot

ProgressFn = Callable[[int, int], None]
CancelFn = Callable[[], bool]
TaskFn = Callable[[ProgressFn, CancelFn], object]

# Strong references to running tasks (and, through them, their runnables) so
# they are not garbage-collected mid-run.
_running: set[Task] = set()


class Task(QObject):
    """Handle of a background task. Public signals are emitted in the GUI thread.

    Exactly one of ``finished`` / ``failed`` / ``cancelled`` is emitted. After
    ``cancel()`` the task always ends with ``cancelled``: a result that arrives
    late (an uninterruptible call such as urllib) is dropped.
    """

    progress = Signal(int, int)
    finished = Signal(object)
    failed = Signal(str)
    cancelled = Signal()

    # Worker -> GUI hand-off. Connected to a slot of this QObject (which lives
    # in the GUI thread), so emission from the worker is queued.
    _event = Signal(str, object)

    def __init__(self, fn: TaskFn, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._fn = fn
        self._cancel_event = threading.Event()
        self._done = False
        self._event.connect(self._dispatch)
        self._runnable = _Runnable(self)

    def cancel(self) -> None:
        """Request cancellation (call from the GUI thread). ``cancelled`` is
        emitted at once; whatever the worker produces later is ignored."""
        self._cancel_event.set()
        if not self._done:
            self._done = True
            _running.discard(self)
            self.cancelled.emit()

    def is_cancelled(self) -> bool:
        return self._cancel_event.is_set()

    def _start(self) -> None:
        _running.add(self)
        QThreadPool.globalInstance().start(self._runnable)

    @Slot(str, object)
    def _dispatch(self, kind: str, payload: object) -> None:
        if self._done:
            return
        if kind == "progress":
            if not self.is_cancelled():
                self.progress.emit(*payload)  # type: ignore[misc]
            return
        self._done = True
        _running.discard(self)
        if kind == "cancelled":
            self.cancelled.emit()
        elif kind == "failed":
            self.failed.emit(str(payload))
        else:
            self.finished.emit(payload)


class _Runnable(QRunnable):
    def __init__(self, task: Task) -> None:
        super().__init__()
        self.setAutoDelete(False)
        self._task = task

    def run(self) -> None:
        task = self._task
        def emit(kind: str, payload: object) -> None:
            try:
                task._event.emit(kind, payload)
            except RuntimeError:  # Task was destroyed with its parent
                _running.discard(task)

        try:
            result = task._fn(
                lambda done, total: emit("progress", (int(done), int(total))),
                task.is_cancelled,
            )
        except CancelledError:
            emit("cancelled", None)
        except BaseException as exc:  # noqa: BLE001 - reported through ``failed``
            emit("failed", str(exc) or type(exc).__name__)
        else:
            emit("finished", result)
        finally:
            _running.discard(task)


def run_task(fn: TaskFn, parent: QObject) -> Task:
    """Run ``fn(progress, is_cancelled)`` in the global thread pool."""
    task = Task(fn, parent)
    task._start()
    return task
