"""Editing-session sidecar files (``<name>.ab1kit.json`` next to the .ab1)."""

from __future__ import annotations

import json
from pathlib import Path

REQUIRED_KEYS = ("seq", "quality", "peaks", "edited", "reversed_complement")


def session_path(ab1_path: str | Path) -> Path:
    ab1 = Path(ab1_path)
    return ab1.parent / (ab1.stem + ".ab1kit.json")


def save_session(path: str | Path, data: dict) -> None:
    payload = dict(data)
    payload["app"] = "ab1kit"
    payload["version"] = 1
    Path(path).write_text(
        json.dumps(payload, ensure_ascii=False), encoding="utf-8"
    )


def load_session(path: str | Path) -> dict | None:
    """Return the session dict, or None when missing/corrupt/not ours."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    if not isinstance(data, dict) or data.get("app") != "ab1kit":
        return None
    if any(key not in data for key in REQUIRED_KEYS):
        return None
    return data
