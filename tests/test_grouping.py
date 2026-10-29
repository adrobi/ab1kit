from pathlib import Path

import pytest

from chrom_core.clones import (
    DEFAULT_PRIMER_TOKENS,
    Grouping,
    ReadRef,
    group_reads,
    parse_name,
)
from chrom_core.config import AppConfig, load_config
from tests.clone_fixtures import REAL_READS


@pytest.mark.parametrize("stem,expected", [
    ("Achl_ACHLO006-09_1_F", ("Achl_ACHLO006-09", "F")),
    ("Achl_ACHLO006-09_2_R", ("Achl_ACHLO006-09", "R")),
    ("cl3_T7", ("cl3", "T7")),
    ("cl3_T7_2", ("cl3", "T7")),          # repeat sequencing
    ("A05_K3-F", ("K3", "F")),
    ("clone12.rev", ("clone12", "rev")),
    ("клон 14_SP6", ("клон 14", "SP6")),
    ("T7", None),
])
def test_parse_name(stem, expected):
    assert parse_name(stem) == expected


def test_parse_name_fallback_to_last_token_and_custom_tokens():
    assert parse_name("cl3_xyz") == ("cl3", "xyz")
    assert parse_name("cl3_XYZ", ["xyz"]) == ("cl3", "XYZ")
    assert parse_name("single") is None


def test_parse_name_multi_part_primer_token():
    assert parse_name("cl5_pUC-F") == ("cl5", "pUC-F")


def test_default_tokens_contents():
    assert "T7TERM" in DEFAULT_PRIMER_TOKENS and len(DEFAULT_PRIMER_TOKENS) == 17


def test_group_real_reads_into_four_clones():
    g = group_reads(REAL_READS)
    assert len(g.clones()) == 4 and all(len(v) == 2 for v in g.clones().values())
    assert g.unassigned == []


def test_repeat_reads_both_kept():
    g = group_reads([Path("d/cl3_T7.ab1"), Path("d/cl3_T7_2.ab1")])
    assert [r.path.name for r in g.clones()["cl3"]] == ["cl3_T7.ab1", "cl3_T7_2.ab1"]


def test_custom_pattern_and_unassigned():
    g = group_reads(
        [Path("S1-fw.ab1"), Path("S12-rv.ab1"), Path("junk.ab1")],
        pattern=r"(?P<clone>S\d+)-(?P<primer>\w+)",
    )
    assert list(g.clones()) == ["S1", "S12"]
    assert g.clones()["S12"][0].primer == "rv"
    assert g.unassigned == [Path("junk.ab1")]


def test_bad_regex_raises():
    with pytest.raises(ValueError):
        group_reads([Path("a.ab1")], pattern="(unclosed")
    with pytest.raises(ValueError):
        group_reads([Path("a.ab1")], pattern=r"(\w+)")  # no named groups


def test_no_clone_goes_to_unassigned():
    g = group_reads([Path("T7.ab1")])
    assert g.reads == [] and g.unassigned == [Path("T7.ab1")]


def test_natural_sort():
    paths = [Path(f"cl{n}_T7.ab1") for n in (10, 2, 1)]
    assert list(group_reads(paths).clones()) == ["cl1", "cl2", "cl10"]


def test_grouping_roundtrip_and_manual_assign():
    g = group_reads([Path("cl1_T7.ab1"), Path("junk.ab1")])
    g.assign(Path("junk.ab1"), "cl1", "SP6")
    g.assign(Path("cl1_T7.ab1"), "cl2", "T7")  # reassign moves, never duplicates
    assert g.unassigned == []
    assert [(r.clone, r.primer) for r in g.reads] == [("cl2", "T7"), ("cl1", "SP6")]
    d = g.to_dict()
    assert all(isinstance(r["path"], str) for r in d["reads"])
    assert Grouping.from_dict(d) == g
    assert ReadRef(Path("a"), "c", "p") in Grouping.from_dict(
        {"reads": [{"path": "a", "clone": "c", "primer": "p"}], "unassigned": []}
    ).reads


def test_config_new_fields_and_old_files(tmp_path):
    cfg = AppConfig()
    assert cfg.primer_tokens == list(DEFAULT_PRIMER_TOKENS) and cfg.group_pattern == ""
    f = tmp_path / "c.json"
    f.write_text('{"theme": "dark"}', encoding="utf-8")
    assert load_config(f).primer_tokens == list(DEFAULT_PRIMER_TOKENS)
