"""Timing tests are opt-in (``-m perf``) and their bounds scale with AB1KIT_PERF_FACTOR."""

import subprocess
import sys
from pathlib import Path

from clone_fixtures import perf_factor

ROOT = Path(__file__).resolve().parent.parent
PERF = {"tests/test_clone_verify.py::test_performance_12_clones_x_4_reads",
        "tests/test_construct_ui.py::test_paint_cost_10kb_12_rows"}


def _collected(*extra: str) -> set[str]:
    out = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-q", "-p", "no:cacheprovider",
         "tests/test_clone_verify.py", "tests/test_construct_ui.py", *extra],
        cwd=ROOT, capture_output=True, text=True, check=False).stdout
    return {line.strip().replace("\\", "/") for line in out.splitlines() if "::" in line}


def test_perf_tests_deselected_by_default_and_selected_with_marker():
    assert not PERF & _collected()
    assert _collected("-m", "perf") == PERF


def test_perf_factor(monkeypatch):
    monkeypatch.delenv("AB1KIT_PERF_FACTOR", raising=False)
    assert perf_factor() == 1.0
    monkeypatch.setenv("AB1KIT_PERF_FACTOR", "2.5")
    assert perf_factor() == 2.5
