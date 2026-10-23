"""Application settings persisted as JSON in the user config directory."""

from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass, field
from pathlib import Path


@dataclass
class AppConfig:
    trim_window: int = 20
    trim_min_q: int = 20
    theme: str = "light"        # "light" | "dark"
    language: str = "ru"        # "ru" | "en"
    recent_files: list[str] = field(default_factory=list)


def config_dir() -> Path:
    appdata = os.environ.get("APPDATA")
    if appdata:
        return Path(appdata) / "ab1kit"
    xdg = os.environ.get("XDG_CONFIG_HOME")
    if xdg:
        return Path(xdg) / "ab1kit"
    return Path.home() / ".config" / "ab1kit"


def config_path() -> Path:
    return config_dir() / "config.json"


def load_config(path: Path | None = None) -> AppConfig:
    """Load settings; unknown keys are dropped, missing/corrupt file -> defaults."""
    path = path or config_path()
    data: dict = {}
    if path.exists():
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            data = {}
    if not isinstance(data, dict):
        data = {}
    cfg = AppConfig()
    for key, value in data.items():
        if hasattr(cfg, key):
            setattr(cfg, key, value)
    return cfg


def save_config(cfg: AppConfig, path: Path | None = None) -> None:
    path = path or config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(asdict(cfg), indent=2, ensure_ascii=False), encoding="utf-8"
    )
