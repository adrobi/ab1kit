import json
import shutil
from pathlib import Path

import pytest

from chrom_core.clones import Grouping, ReadRef, VerifyParams
from chrom_core.project import PROJECT_VERSION, Project, load_project, save_project


def _make(folder: Path, outside: Path | None = None) -> Project:
    reads = [
        ReadRef(folder / "reads" / "c1_F.ab1", "c1", "F"),
        ReadRef(folder / "reads" / "c1_R.ab1", "c1", "R"),
    ]
    unassigned = [folder / "x.ab1"]
    if outside is not None:
        reads.append(ReadRef(outside / "far.ab1", "c2", "F"))
    return Project(
        construct_path=folder / "plasmid.gb",
        circular=True,
        region=(900, 100),
        grouping=Grouping(reads=reads, unassigned=unassigned),
        params=VerifyParams(min_q=25, het_ratio=0.3, mixed_min_positions=4),
    )


def test_round_trip(tmp_path):
    proj = _make(tmp_path, outside=tmp_path.parent / "elsewhere")
    f = tmp_path / "s.ab1proj"
    save_project(f, proj)
    assert load_project(f) == proj


def test_whole_construct_region_none(tmp_path):
    proj = _make(tmp_path)
    proj.region = None
    f = tmp_path / "s.ab1proj"
    save_project(f, proj)
    assert load_project(f).region is None


def test_paths_relative_and_posix_inside_folder(tmp_path):
    outside = tmp_path.parent / "elsewhere"
    f = tmp_path / "s.ab1proj"
    save_project(f, _make(tmp_path, outside))
    data = json.loads(f.read_text(encoding="utf-8"))
    assert data["app"] == "ab1kit"
    assert data["version"] == PROJECT_VERSION
    assert data["construct_path"] == "plasmid.gb"
    paths = [r["path"] for r in data["grouping"]["reads"]]
    assert paths[0] == "reads/c1_F.ab1"
    assert Path(paths[2]).is_absolute()


def test_moved_folder_still_resolves(tmp_path):
    src = tmp_path / "a"
    src.mkdir()
    save_project(src / "s.ab1proj", _make(src))
    dst = tmp_path / "b" / "moved"
    dst.parent.mkdir()
    shutil.move(str(src), str(dst))
    loaded = load_project(dst / "s.ab1proj")
    assert loaded == _make(dst)


def test_cyrillic_folder_with_space(tmp_path):
    folder = tmp_path / "Данные" / "клоны 1"
    folder.mkdir(parents=True)
    f = folder / "проект.ab1proj"
    proj = _make(folder)
    save_project(f, proj)
    assert "Данные" not in f.read_text(encoding="utf-8")
    assert load_project(f) == proj
    moved = tmp_path / "Другие данные" / "клоны 1"
    moved.parent.mkdir()
    shutil.move(str(folder), str(moved))
    assert load_project(moved / "проект.ab1proj") == _make(moved)


def test_corrupt_json(tmp_path):
    f = tmp_path / "s.ab1proj"
    f.write_text("{not json", encoding="utf-8")
    with pytest.raises(ValueError, match="JSON"):
        load_project(f)


def test_foreign_file(tmp_path):
    f = tmp_path / "s.ab1proj"
    f.write_text(json.dumps({"app": "other", "version": 1}), encoding="utf-8")
    with pytest.raises(ValueError, match="ab1kit"):
        load_project(f)


def test_future_version(tmp_path):
    f = tmp_path / "s.ab1proj"
    f.write_text(
        json.dumps({"app": "ab1kit", "version": PROJECT_VERSION + 1}),
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="version"):
        load_project(f)


def test_missing_fields(tmp_path):
    f = tmp_path / "s.ab1proj"
    f.write_text(json.dumps({"app": "ab1kit", "version": 1}), encoding="utf-8")
    with pytest.raises(ValueError):
        load_project(f)
