"""Минималистичный i18n: русские строки здесь, английские в translations_en.py.

Использование: `tr("menu.file.open")`, с подстановкой: `tr("search.found", n=3)`.
"""

from __future__ import annotations

from translations_en import EN

RU = {
    "menu.file": "&Файл",
    "menu.file.open": "&Открыть .ab1…",
    "menu.file.recent": "Недавние",
    "menu.file.exit": "В&ыход",
    "menu.edit": "&Правка",
    "menu.edit.find": "&Поиск подпоследовательности…",
    "menu.edit.find_next": "Найти &далее",
    "menu.edit.goto": "&Перейти к базе…",
    "menu.view": "&Вид",
    "menu.view.dark_theme": "Тёмная тема",
    "menu.view.smoothing": "Сглаживание трассов",
    "menu.view.quality_profile": "Профиль качества",
    "menu.view.file_info": "Информация о файле",
    "menu.help": "&Справка",
    "menu.help.about": "О программе",
    "tb.open": "Открыть",
    "tb.zoom_in": "Приблизить",
    "tb.zoom_out": "Отдалить",
    "tb.fit": "Весь рид",
    "tb.reset": "Сброс вида",
    "goto.title": "Перейти к базе",
    "goto.label": "Номер базы",
    "search.title": "Поиск подпоследовательности",
    "search.pattern": "Последовательность",
    "search.reverse": "Искать и в реверс-комплементе",
    "search.button": "Найти",
    "search.found": "Найдено: {n} (прямой: {f}, реверс: {r})",
    "search.none": "Совпадений нет",
    "info.title": "Информация",
    "info.file": "Файл",
    "info.length": "Длина, баз",
    "info.mean_q": "Среднее качество",
    "info.q20": "Баз с Q≥20",
    "meta.sample": "Образец",
    "meta.instrument": "Прибор",
    "meta.run_date": "Дата прогона",
    "meta.run_time": "Время прогона",
    "meta.comment": "Комментарий",
    "meta.project": "Проект",
    "meta.capillary": "Капилляр",
    "meta.mobility": "Мобильность (dye set)",
    "about.text": "ab1kit — просмотр и анализ хроматограмм Sanger.\nВерсия {version}",
    "status.ready": "Готово",
    "status.file": "Загружено {n} баз, среднее качество {q:.1f}",
    "cursor.info": "База {n}: {base}, качество {q}",
    "view.empty": "Откройте файл .ab1",
    "err.open_failed": "Не удалось открыть файл",
}

_lang = "ru"


def set_language(lang: str) -> None:
    global _lang
    _lang = lang if lang in ("ru", "en") else "ru"


def language() -> str:
    return _lang


def tr(key: str, **kwargs) -> str:
    table = EN if _lang == "en" else RU
    text = table.get(key) or RU.get(key) or key
    return text.format(**kwargs) if kwargs else text
