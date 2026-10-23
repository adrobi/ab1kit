import i18n
from i18n import set_language, tr


def setup_function() -> None:
    set_language("ru")


def test_russian_by_default() -> None:
    assert tr("menu.file.open") == "&Открыть .ab1…"


def test_english_after_switch() -> None:
    set_language("en")
    assert tr("menu.file.open") == "&Open .ab1…"


def test_fallback_to_russian_for_missing_key() -> None:
    set_language("en")
    assert tr("menu.file.open")  # key exists in both
    i18n.EN.pop("menu.file.exit")
    try:
        set_language("en")
        assert tr("menu.file.exit") == "В&ыход"
    finally:
        i18n.EN["menu.file.exit"] = "E&xit"


def test_unknown_key_returned_as_is() -> None:
    assert tr("no.such.key") == "no.such.key"


def test_formatting_kwargs() -> None:
    set_language("ru")
    assert tr("search.found", n=2, f=1, r=1) == "Найдено: 2 (прямой: 1, реверс: 1)"


def test_invalid_language_falls_back_to_ru() -> None:
    set_language("de")
    assert tr("menu.file.open") == "&Открыть .ab1…"
