from pathlib import Path

from chrom_core.config import AppConfig, load_config, save_config


def test_roundtrip(tmp_path: Path) -> None:
    path = tmp_path / "config.json"
    cfg = AppConfig(trim_window=30, theme="dark", recent_files=["a.ab1"])
    save_config(cfg, path)
    assert load_config(path) == cfg


def test_missing_file_gives_defaults(tmp_path: Path) -> None:
    assert load_config(tmp_path / "nope.json") == AppConfig()


def test_corrupt_file_gives_defaults(tmp_path: Path) -> None:
    path = tmp_path / "config.json"
    path.write_text("{broken json", encoding="utf-8")
    assert load_config(path) == AppConfig()


def test_unknown_keys_dropped(tmp_path: Path) -> None:
    path = tmp_path / "config.json"
    path.write_text('{"trim_window": 10, "mystery": 1}', encoding="utf-8")
    cfg = load_config(path)
    assert cfg.trim_window == 10
    assert not hasattr(cfg, "mystery")
