from pathlib import Path

import pytest

from chrom_core import sessions


def test_session_path_sibling(tmp_path: Path) -> None:
    ab1 = tmp_path / "read.ab1"
    assert sessions.session_path(ab1) == tmp_path / "read.ab1kit.json"


def test_save_load_roundtrip(tmp_path: Path) -> None:
    path = tmp_path / "s.json"
    data = {
        "source": "read.ab1",
        "seq": "ACGT",
        "quality": [1, 2, 3, 4],
        "peaks": [0, 10, 20, 30],
        "edited": [[1, "T"]],
        "reversed_complement": False,
    }
    sessions.save_session(path, data)
    loaded = sessions.load_session(path)
    assert loaded is not None
    assert loaded["seq"] == "ACGT"
    assert loaded["app"] == "ab1kit" and loaded["version"] == 1


def test_load_missing_returns_none(tmp_path: Path) -> None:
    assert sessions.load_session(tmp_path / "nope.json") is None


def test_load_corrupt_returns_none(tmp_path: Path) -> None:
    path = tmp_path / "s.json"
    path.write_text("{oops", encoding="utf-8")
    assert sessions.load_session(path) is None


def test_load_foreign_json_returns_none(tmp_path: Path) -> None:
    path = tmp_path / "s.json"
    path.write_text('{"seq": "A"}', encoding="utf-8")
    assert sessions.load_session(path) is None


def test_load_missing_keys_returns_none(tmp_path: Path) -> None:
    path = tmp_path / "s.json"
    sessions.save_session(path, {
        "seq": "A", "quality": [1], "peaks": [0],
        "edited": [], "reversed_complement": False, "source": "x.ab1",
    })
    path.write_text(path.read_text(encoding="utf-8").replace('"peaks"', '"peakz"'), encoding="utf-8")
    assert sessions.load_session(path) is None
