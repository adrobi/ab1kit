"""``.ab1proj`` project files: a clone-verification session on disk."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

from chrom_core.clones import Grouping, Region, VerifyParams

PROJECT_VERSION = 1


@dataclass
class Project:
    construct_path: Path
    circular: bool
    region: Region | None
    grouping: Grouping
    params: VerifyParams


def _to_stored(path: Path, base: Path) -> str:
    """Relative POSIX path when ``path`` is inside ``base``, else absolute."""
    absolute = Path(os.path.abspath(path))
    try:
        return absolute.relative_to(base).as_posix()
    except ValueError:
        return str(absolute)


def _from_stored(text: str, base: Path) -> Path:
    p = Path(text)
    return p if p.is_absolute() else base / p


def save_project(path: Path, project: Project) -> None:
    path = Path(path)
    base = Path(os.path.abspath(path.parent))
    grouping = project.grouping.to_dict()
    for r in grouping["reads"]:
        r["path"] = _to_stored(Path(r["path"]), base)
    grouping["unassigned"] = [
        _to_stored(Path(p), base) for p in grouping["unassigned"]
    ]
    p = project.params
    payload = {
        "app": "ab1kit",
        "version": PROJECT_VERSION,
        "construct_path": _to_stored(project.construct_path, base),
        "circular": project.circular,
        "region": list(project.region) if project.region is not None else None,
        "grouping": grouping,
        "params": {
            "min_q": p.min_q,
            "het_ratio": p.het_ratio,
            "mixed_min_positions": p.mixed_min_positions,
        },
    }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def load_project(path: Path) -> Project:
    path = Path(path)
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as e:
        raise ValueError(f"Project file is not valid JSON: {path} ({e})") from e
    except (OSError, UnicodeDecodeError) as e:
        raise ValueError(f"Cannot read project file: {path} ({e})") from e
    if not isinstance(data, dict) or data.get("app") != "ab1kit":
        raise ValueError(f"Not an ab1kit project file: {path}")
    version = data.get("version")
    if not isinstance(version, int) or version > PROJECT_VERSION:
        raise ValueError(
            f"Unsupported project version {version!r} (this ab1kit supports "
            f"up to version {PROJECT_VERSION}): {path}"
        )
    base = Path(os.path.abspath(path.parent))
    try:
        grouping = Grouping.from_dict(
            {
                "reads": [
                    {**r, "path": str(_from_stored(r["path"], base))}
                    for r in data["grouping"].get("reads", [])
                ],
                "unassigned": [
                    str(_from_stored(p, base))
                    for p in data["grouping"].get("unassigned", [])
                ],
            }
        )
        region = data["region"]
        return Project(
            construct_path=_from_stored(data["construct_path"], base),
            circular=bool(data["circular"]),
            region=(int(region[0]), int(region[1])) if region is not None else None,
            grouping=grouping,
            params=VerifyParams(**data["params"]),
        )
    except (KeyError, TypeError, ValueError, IndexError) as e:
        raise ValueError(f"Project file is incomplete or damaged: {path} ({e!r})") from e
