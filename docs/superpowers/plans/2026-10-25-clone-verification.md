# Проверка клонов, точность пиков и новый интерфейс — план реализации

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** ab1kit проверяет 4–12 клонов одного конструкта по заданному региону и показывает вердикт по каждому клону на экране-матрице, с хроматограммами всех ридов в месте каждого расхождения.

**Architecture:** Вся логика — в UI-независимом `chrom_core/` (новые модули `peaks`, `reference`, `clones`, `effects`, `project`, `clone_report`; кольцевое выравнивание в `align`). Интерфейс на PySide6 Widgets: `main_window.py` разбивается на `ui/read_view/*`, новый экран — `ui/construct/*`, тяжёлые операции — через `ui/tasks.py` (`QThreadPool`). CLI — тонкая оболочка над ядром.

**Tech Stack:** Python 3.11+, PySide6 ≥ 6.6, BioPython ≥ 1.83 (`Bio.Align`, `Bio.SeqIO`), openpyxl, pytest.

**Spec:** `docs/superpowers/specs/2026-10-25-clone-verification-design.md`

## Global Constraints

- Python ≥ 3.11; новых runtime-зависимостей нет (только `biopython`, `PySide6`, `openpyxl`).
- `chrom_core/` не импортирует Qt.
- Все строки UI — через `tr("ключ")` из `i18n.py`; каждый новый ключ есть и в `RU` (`i18n.py`), и в `EN` (`translations_en.py`) — это проверяет `tests/test_i18n.py`. (В спецификации написано `_()` — в проекте функция называется `tr`.)
- Пороги по умолчанию: `min_q = 20`, `het_ratio = 0.40`, `mixed_min_positions = 3`; окно базовой линии `span = 20`; запас кольцевого выравнивания `+50` п.н.; лимит конструкта `MAX_REFERENCE_LENGTH = 100_000`.
- Тесты не ходят в сеть; данные из NCBI/GitHub скачиваются один раз скриптом и коммитятся.
- Данные с несовместимой с MIT лицензией (gear-genomics, GPL-3.0) в репозиторий не кладутся.
- Существующие тесты остаются зелёными после каждой задачи: `python -m pytest -q`.
- GUI-тесты — offscreen (`QT_QPA_PLATFORM=offscreen`), по образцу `tests/test_smoke_gui.py`.
- Палитра статусов: `ok` зелёный, `errors` красный, `incomplete` жёлтый, `mixed` фиолетовый, `no_data` серый — одни и те же цвета везде (матрица, сводка, отчёт XLSX).

## Review Focus

1. **Рид, развёрнутый реверс-комплементом в сессии правок** (`.ab1kit.json` с `reversed_complement: true`): покрытие и расхождения должны относиться к тем же позициям конструкта, что и у неразвёрнутого рида. Тест — Task 7.
2. **Регион через начало координат кольцевого конструкта** (`start > end`, например 4300–200 на 4420 п.н.): покрытие, расхождения и вердикт считаются по обёрнутому региону. Тест — Task 7.
3. **Рид без данных качества** (все Q = 0, как `examples/310.ab1`): расхождения не «подтверждаются», клон получает `incomplete`, а не `errors`; в батч-таблице «—» вместо 0. Тесты — Tasks 3 и 7.
4. **Конструкт с несколькими записями, строчными буквами или `N`**: берётся первая запись, верхний регистр, `N` в референсе не считается расхождением. Тесты — Tasks 4 и 7.
5. **Повторное секвенирование тем же праймером** (`cl3_T7.ab1` и `cl3_T7_2.ab1`) и **кириллица/пробелы в путях**: оба рида остаются в клоне `cl3`; проект и CLI работают с путями вида `D:\Данные\клоны 1\`. Тесты — Tasks 6, 10, 17.

---

## Фаза 1 — данные, анализ пиков, исправления

### Task 1: Корпус реальных тестовых данных

**Files:**
- Create: `scripts/fetch_test_data.py`
- Create: `tests/data/real/` (8 файлов `.ab1`), `tests/data/refs/`, `tests/data/README.md`
- Create: `tests/clone_fixtures.py`
- Test: `tests/test_data_corpus.py`

**Interfaces:**
- Produces (`tests/clone_fixtures.py`):
  - `DATA: Path` = `tests/data`
  - `REAL_READS: list[Path]` — 8 файлов, отсортированы
  - `real_read(name: str) -> Chromatogram`
  - `COI_GB: Path`, `PUC19_GB: Path`, `SNAPGENE_DNA: Path`
  - `coi_reference() -> str` — последовательность COI, верхний регистр

- [ ] **Step 1: Написать `scripts/fetch_test_data.py`.** Скачивает в `tests/data/`:
  - из `https://raw.githubusercontent.com/roblanf/sangeranalyseR/<commit>/inst/extdata/Allolobophora_chlorotica/ACHLO/` файлы `Achl_ACHLO006-09_{1_F,2_R}.ab1`, `Achl_ACHLO007-09_{1_F,2_R}.ab1`, `Achl_ACHLO040-09_{1_F,2_R}.ab1`, `Achl_ACHLO041-09_{1_F,2_R}.ab1` → `tests/data/real/`; `<commit>` — SHA HEAD на момент скачивания, зашит константой;
  - pUC19 `L09137` в формате GenBank через `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=L09137&rettype=gbwithparts&retmode=text` → `tests/data/refs/pUC19_L09137.gb`;
  - COI: accession передаётся аргументом `--coi ACC` → `tests/data/refs/coi_<ACC>.gb` (тот же efetch);
  - SnapGene-образец из `https://raw.githubusercontent.com/biopython/biopython/<commit>/Tests/SnapGene/sample-d.dna` → `tests/data/refs/sample-d.dna`.
  Только `urllib`; идемпотентен (существующие файлы не перекачивает без `--force`).
- [ ] **Step 2: Подобрать accession COI.** Отправить `Achl_ACHLO006-09_1_F` в BLAST (`chrom_core.blast.build_request`/`submit`, результат открыть вручную) и выбрать запись GenBank *Allolobophora chlorotica* с аннотированной CDS, покрывающую рид целиком с идентичностью ≥ 98%. Запустить `python scripts/fetch_test_data.py --coi <ACC>`.
- [ ] **Step 3: Написать `tests/data/README.md`** — источник каждого файла, commit/accession, лицензия (sangeranalyseR — MIT; BioPython — Biopython License; NCBI — public domain), пометка «полусинтетические данные генерируются в `tests/clone_fixtures.py`, на диске не лежат».
- [ ] **Step 4: Написать падающий тест**

```python
# tests/test_data_corpus.py
from clone_fixtures import REAL_READS, COI_GB, PUC19_GB, SNAPGENE_DNA, real_read, coi_reference
from chrom_core.align import align_read

def test_corpus_files_present():
    assert len(REAL_READS) == 8
    assert COI_GB.exists() and PUC19_GB.exists() and SNAPGENE_DNA.exists()

def test_every_real_read_aligns_to_coi_with_high_identity():
    ref = coi_reference()
    for path in REAL_READS:
        aln = align_read(real_read(path.name).seq, ref)
        assert aln.identity >= 0.97, path.name
        assert aln.strand == ("+" if "_F" in path.name else "-")
```

- [ ] **Step 5: Запустить** `python -m pytest tests/test_data_corpus.py -v` → FAIL (`ModuleNotFoundError: clone_fixtures`).
- [ ] **Step 6: Реализовать `tests/clone_fixtures.py`** (константы и функции из Interfaces; `coi_reference` читает `COI_GB` через `Bio.SeqIO`). Добавить `pythonpath = ["tests"]` в `[tool.pytest.ini_options]` в `pyproject.toml`.
- [ ] **Step 7: Запустить** тест → PASS. Если идентичность ниже порога у отдельных ридов — сменить accession (Step 2), а не порог.
- [ ] **Step 8: Commit** `git add scripts/fetch_test_data.py tests/data tests/clone_fixtures.py tests/test_data_corpus.py pyproject.toml && git commit -m "Test corpus: real sangeranalyseR reads, COI and pUC19 references"`

### Task 2: `peaks.py` — нормированный анализ вторичных пиков

**Files:**
- Create: `chrom_core/peaks.py`
- Modify: `chrom_core/snps.py` (использует `peaks`), `chrom_core/verify.py:166-181` (то же)
- Test: `tests/test_peaks.py`; регрессия — `tests/test_snps.py`, `tests/test_verify.py`

**Interfaces:**
- Consumes: `Chromatogram` (`chrom_core/abif.py`), `ReadAlignment` (`chrom_core/align.py`)
- Produces (`chrom_core/peaks.py`):
  - `BASELINE_SPAN = 20`
  - `peak_window(chrom: Chromatogram, i: int) -> int` — половина расстояния до ближайшего соседнего пика, не меньше 2
  - `channel_baseline(chrom: Chromatogram, i: int, span: int = BASELINE_SPAN) -> dict[str, float]` — для каждого канала медиана высот пиков тех из `span` ближайших (по индексу) вызванных баз, у которых этот канал — вызванная база; если таких нет — медиана высот всех `span` пиков по всем каналам
  - `@dataclass(frozen=True) SecondaryPeak(primary_base: str, base: str, ratio: float | None)`
  - `secondary_peak(chrom: Chromatogram, i: int) -> SecondaryPeak` — высоты в окне `peak_window`, делённые на `channel_baseline`; `ratio = норм. второй вне вызванных / норм. максимальный среди вызванных`
  - `scan_secondary_peaks(chrom: Chromatogram, alignment: ReadAlignment, min_ratio: float, min_q: int) -> list[tuple[int, SecondaryPeak]]` — по всем `ref_pos` из `alignment.ref_to_read` (в ориентации `chrom`, как в `call_snps`), где `quality ≥ min_q` и `ratio ≥ min_ratio`; список `(ref_pos, SecondaryPeak)` по возрастанию `ref_pos`
  - `peak_heights` остаётся в `snps.py` с прежней сигнатурой (обратная совместимость, используется UI-тултипами)

- [ ] **Step 1: Написать падающие тесты** (синтетика через `make_chrom` из `tests/test_snps.py` — перенести её в `tests/clone_fixtures.py` и импортировать в обоих местах):

```python
def test_weak_channel_normalisation_removes_false_het():
    # G-канал систематически вдвое слабее: «второй пик» G на 45% сырой высоты
    # при норм. высоте ≈ 0.9·baseline A → ratio по сырым 0.45 (ложный het),
    # после нормировки — сам G-пик нормален, а под A-пиком G 45/500 → < 0.40
    ...
    assert secondary_peak(chrom, 10).ratio < 0.40

def test_true_het_detected_after_normalisation():
    ...  # второй пик C на 50% → ratio ≈ 0.5
    assert secondary_peak(chrom, 10).base == "C"
    assert 0.45 <= secondary_peak(chrom, 10).ratio <= 0.55

def test_peak_window_tracks_spacing():
    assert peak_window(chrom_spacing_10, 5) == 5
    assert peak_window(chrom_spacing_3, 5) == 2

def test_scan_finds_het_where_reference_base_was_called():
    # рид = референс, но под позицией 30 второй пик T на 50%
    hits = scan_secondary_peaks(chrom, align_read(chrom.seq, ref), 0.40, 20)
    assert [pos for pos, _ in hits] == [30]
```

  Конкретные значения трасс — в теле тестов; утверждения — как выше.
- [ ] **Step 2: Запустить** `python -m pytest tests/test_peaks.py -v` → FAIL (`ModuleNotFoundError: chrom_core.peaks`).
- [ ] **Step 3: Реализовать `chrom_core/peaks.py`** по Interfaces.
- [ ] **Step 4: Перевести `call_snps` и `_verify_one` на `secondary_peak`.** Логика вердиктов не меняется, меняется только источник `ratio`/`secondary_base`.
- [ ] **Step 5: Добавить регрессию на реальных ридах** в `tests/test_peaks.py`:

```python
def test_real_reads_have_no_mass_false_hets():
    ref = coi_reference()
    for path in REAL_READS:
        chrom = real_read(path.name)
        aln = align_read(chrom.seq, ref)
        hits = scan_secondary_peaks(chrom, aln, 0.40, 20)
        assert len(hits) <= 2, (path.name, hits)
```

- [ ] **Step 6: Запустить всё** `python -m pytest -q` → PASS (включая `test_snps.py`, `test_verify.py`). Если старый тест SNP падает из-за нормировки, разобрать конкретный случай; менять ожидание теста можно только если новое значение ближе к истине — с комментарием в тесте.
- [ ] **Step 7: Commit** `git commit -m "peaks: normalised secondary-peak analysis at any read position"`

### Task 3: Исправить две ошибки интерфейса

**Files:**
- Modify: `chrom_core/verify.py` (`VariantVerdict`), `chrom_core/batch.py` (`FileReport`, `scan_folder`), `ui/main_window.py:722-757` (`_fill_verify_table`), `ui/batch_panel.py:76`, `i18n.py`, `translations_en.py`
- Test: `tests/test_verify.py`, `tests/test_batch.py`, `tests/test_smoke_gui.py`

**Interfaces:**
- Produces:
  - `VariantVerdict.note_code: str = ""` — коды: `only_snv`, `pos_outside_ref`, `ref_differs`, `deleted_in_read`, `outside_coverage`, `low_q` (с `note_args={"q": ..., "min_q": ...}`); поле `note_args: dict = field(default_factory=dict)`. Поле `note` (английский текст) остаётся — его пишут CSV/XLSX-отчёты.
  - i18n-ключи `verify.note.<code>` в RU и EN.
  - `FileReport.mean_q: float | None` — `None`, если у рида все Q = 0.

- [ ] **Step 1: Падающие тесты**

```python
# test_verify.py
def test_ref_mismatch_carries_note_code(): ...
    assert verdict.note_code == "ref_differs"

# test_batch.py
def test_read_without_quality_reports_none_mean_q(folder):
    by_name = {r.path.name: r for r in batch.scan_folder(folder)}
    assert by_name["310.ab1"].mean_q is None

# test_smoke_gui.py
def test_verify_table_note_is_localised(qapp, ...):
    # русский UI, вариант с неверной ref-базой → ячейка «Примечание»
    assert "table ref base" not in note_cell_text
    assert note_cell_text == tr("verify.note.ref_differs")
```

- [ ] **Step 2: Запустить** → FAIL.
- [ ] **Step 3: Реализовать:** `_verify_one` заполняет `note_code`/`note_args`; `_fill_verify_table` показывает `tr(f"verify.note.{code}", **args)` (если код пуст — `note`); `scan_folder` ставит `mean_q=None` при `not any(chrom.quality)`; `BatchPanel` показывает `"—"` для `None`; `status` не меняется (рид без качества по-прежнему оценивается по длине обрезки). Проверить, что `summary_csv` пишет пустую ячейку для `None`.
- [ ] **Step 4: Запустить** `python -m pytest -q` → PASS.
- [ ] **Step 5: Commit** `git commit -m "Fix untranslated verify notes and Q=0 for reads without quality"`

---

## Фаза 2 — ядро проверки клонов

### Task 4: `reference.py` — загрузка конструкта

**Files:**
- Create: `chrom_core/reference.py`
- Test: `tests/test_reference.py`

**Interfaces:**
- Produces:
  - `@dataclass(frozen=True) Feature(name: str, kind: str, start: int, end: int, strand: int)` — 0-based, полуоткрытый; `start > end` — фича через начало координат; `codon_start: int = 1` (из квалификатора `codon_start` для CDS)
  - `@dataclass(frozen=True) Construct(name: str, seq: str, circular: bool, features: tuple[Feature, ...] = ())` с методом `with_circular(flag: bool) -> Construct`
  - `load_construct(path: str | Path) -> Construct` — `.fa .fasta .fna .txt` → `fasta`, `.gb .gbk .genbank` → `genbank`, `.dna` → `snapgene`; неизвестное расширение — пробовать по первому символу (`>` → fasta, `LOCUS` → genbank), иначе `ValueError`. Первая запись; `seq.upper()`; недопустимые символы → `ValueError` (как `align.parse_sequence`). `circular` — из `record.annotations["topology"] == "circular"`; для FASTA — `True`. Имя фичи — первый из квалификаторов `label`, `gene`, `product`, `note`, иначе `kind`. Составные локации (`CompoundLocation`) через начало: `start` первой части, `end` последней.
  - `construct_from_text(text: str, name: str = "pasted") -> Construct` — через `align.parse_sequence`, `circular=True`

- [ ] **Step 1: Падающие тесты**

```python
def test_genbank_puc19():
    c = load_construct(PUC19_GB)
    assert len(c.seq) == 2686 and c.circular
    assert any(f.kind == "CDS" for f in c.features)

def test_snapgene_sample(): c = load_construct(SNAPGENE_DNA); assert c.seq and c.features

def test_fasta_multi_record_lowercase(tmp_path):
    p = tmp_path / "x.fasta"; p.write_text(">a\nacgtn\n>b\nTTTT\n")
    c = load_construct(p)
    assert c.seq == "ACGTN" and c.circular and c.features == ()

def test_unknown_extension_sniffed(tmp_path): ...  # "x.seq" с ">"-заголовком → FASTA

def test_garbage_rejected(tmp_path): ...  # pytest.raises(ValueError)

def test_with_circular_toggle(): assert not construct_from_text("ACGT").with_circular(False).circular
```

- [ ] **Step 2: Запустить** `python -m pytest tests/test_reference.py -v` → FAIL.
- [ ] **Step 3: Реализовать** по Interfaces.
- [ ] **Step 4: Запустить** → PASS.
- [ ] **Step 5: Commit** `git commit -m "reference: load constructs from FASTA, GenBank and SnapGene"`

### Task 5: Кольцевое выравнивание

**Files:**
- Modify: `chrom_core/align.py` (новые функции; `align_read` не меняется)
- Test: `tests/test_align_circular.py`

**Interfaces:**
- Consumes: `Construct` (Task 4)
- Produces (`chrom_core/align.py`):
  - `CIRCULAR_MARGIN = 50`
  - `ReadAlignment` получает поле `ref_length: int = 0` и методы `covers(pos: int) -> bool`, `covered_positions() -> list[int]` (для обёрнутых — `ref_start..L-1, 0..ref_end-1`); `ref_end < ref_start` означает обёртку
  - `align_to_construct(read: str, construct: Construct) -> ReadAlignment` — линейный: `align_read` + `ref_length`; кольцевой: `align_read(read, seq + seq[:len(read) + CIRCULAR_MARGIN])` (лимит проверяется по `len(construct.seq)`), затем все координаты `ref_*`, ключи `ref_to_read`, `Difference.ref_pos` — по модулю `L`; для `ref_end` используется `((ref_end - 1) % L) + 1`
  - `gapped_sequences` для обёрнутого выравнивания работает по развёрнутому референсу (`seq + seq[:…]`) — передавать ему `construct.seq * 2`

- [ ] **Step 1: Падающие тесты** (на реальном риде и COI):

```python
def test_linear_matches_align_read(): ...  # те же differences/identity, ref_length == len(ref)

def test_read_across_origin_of_rotated_reference():
    ref = coi_reference(); aln0 = align_read(read.seq, ref)
    mid = (aln0.ref_start + aln0.ref_end) // 2
    rotated = ref[mid:] + ref[:mid]
    aln = align_to_construct(read.seq, Construct("rot", rotated, True))
    assert aln.ref_end < aln.ref_start            # обёрнуто
    assert len(aln.covered_positions()) == aln0.covered
    assert [ (d.ref_pos - (len(ref) - mid)) % len(ref) for d in aln.differences ] \
        == [ d.ref_pos for d in aln0.differences ]
    assert aln.covers(0) and aln.covers(len(ref) - 1)

def test_same_rotation_linear_breaks_read():
    aln = align_to_construct(read.seq, Construct("rot", rotated, False))
    assert aln.covered < aln0.covered            # без кольца рид рвётся

def test_limit_applies_to_original_length(): ...  # 100_000 п.н. кольцо проходит; 100_001 — ValueError
```

- [ ] **Step 2: Запустить** → FAIL.
- [ ] **Step 3: Реализовать** по Interfaces.
- [ ] **Step 4: Запустить** `python -m pytest -q` → PASS (включая старые `test_align.py`).
- [ ] **Step 5: Commit** `git commit -m "align: circular constructs with origin-spanning reads"`

### Task 6: Группировка ридов в клоны

**Files:**
- Create: `chrom_core/clones.py` (раздел группировки)
- Modify: `chrom_core/config.py` (`AppConfig.primer_tokens: list[str]`, `AppConfig.group_pattern: str = ""`)
- Test: `tests/test_grouping.py`

**Interfaces:**
- Produces:
  - `DEFAULT_PRIMER_TOKENS = ("F", "R", "FW", "FWD", "REV", "RV", "T7", "T7TERM", "SP6", "T3", "M13F", "M13R", "M13", "PUC-F", "PUC-R", "BGH", "CMV")`
  - `@dataclass(frozen=True) ReadRef(path: Path, clone: str, primer: str)`
  - `@dataclass Grouping(reads: list[ReadRef], unassigned: list[Path])` с методами `clones() -> dict[str, list[ReadRef]]` (порядок клонов — натуральная сортировка: `cl2` < `cl10`), `assign(path: Path, clone: str, primer: str) -> None`, `to_dict() -> dict`, `from_dict(d: dict) -> Grouping`
  - `parse_name(stem: str, primer_tokens: Iterable[str] = DEFAULT_PRIMER_TOKENS) -> tuple[str, str] | None` — алгоритм из спецификации §2.4 (разделители `_ - .` и пробел; токен праймера — последний совпавший со словарём без учёта регистра, иначе последний токен; клон — исходная подстрока до токена праймера без хвостовых разделителей, без одиночной цифры прямо перед праймером и без ведущей лунки `^[A-H]\d{1,2}$`; хвостовой суффикс повтора `_2`/`-2` **после** токена праймера отбрасывается). Пустой клон → `None`.
  - `group_reads(paths: Iterable[Path], pattern: str | None = None, primer_tokens=DEFAULT_PRIMER_TOKENS) -> Grouping` — при `pattern` используется `re.fullmatch` с группами `clone`, `primer`; несовпавшие → `unassigned`. Невалидный regex → `ValueError`.

- [ ] **Step 1: Падающие тесты**

```python
@pytest.mark.parametrize("stem,expected", [
    ("Achl_ACHLO006-09_1_F", ("Achl_ACHLO006-09", "F")),
    ("Achl_ACHLO006-09_2_R", ("Achl_ACHLO006-09", "R")),
    ("cl3_T7", ("cl3", "T7")),
    ("cl3_T7_2", ("cl3", "T7")),          # повторное секвенирование
    ("A05_K3-F", ("K3", "F")),
    ("clone12.rev", ("clone12", "rev")),
    ("клон 14_SP6", ("клон 14", "SP6")),
    ("T7", None),
])
def test_parse_name(stem, expected): assert parse_name(stem) == expected

def test_group_real_reads_into_four_clones():
    g = group_reads(REAL_READS)
    assert len(g.clones()) == 4 and all(len(v) == 2 for v in g.clones().values())

def test_repeat_reads_both_kept(): ...  # cl3_T7.ab1 + cl3_T7_2.ab1 → 2 рида в cl3

def test_custom_pattern_and_unassigned(): ...  # r"(?P<clone>S\d+)-(?P<primer>\w+)" ; "junk" → unassigned

def test_bad_regex_raises(): ...

def test_natural_sort(): ...  # cl2, cl10

def test_grouping_roundtrip_and_manual_assign(): ...
```

- [ ] **Step 2: Запустить** → FAIL.
- [ ] **Step 3: Реализовать** по Interfaces; `AppConfig` — новые поля со значениями по умолчанию (`list(DEFAULT_PRIMER_TOKENS)`, `""`).
- [ ] **Step 4: Запустить** `python -m pytest -q` → PASS.
- [ ] **Step 5: Commit** `git commit -m "clones: group reads into clones by file name"`

### Task 7: Проверка клона и вердикт

**Files:**
- Modify: `chrom_core/clones.py` (раздел проверки)
- Modify: `tests/clone_fixtures.py` (полусинтетика)
- Test: `tests/test_clone_verify.py`

**Interfaces:**
- Consumes: `Construct` (Task 4), `align_to_construct`, `ReadAlignment.covers/covered_positions` (Task 5), `secondary_peak`, `scan_secondary_peaks` (Task 2), `ReadRef` (Task 6), `sessions.load_session/session_path`, `EditableRead.restore_state/current`
- Produces (`chrom_core/clones.py`):
  - `@dataclass(frozen=True) VerifyParams(min_q: int = 20, het_ratio: float = 0.40, mixed_min_positions: int = 3)`
  - `Region = tuple[int, int]` — 0-based полуоткрытый; `start > end` — через начало (только для кольцевого); `region_positions(region, length) -> list[int]`; `None` = весь конструкт
  - `@dataclass(frozen=True) ReadResult(ref: ReadRef, alignment: ReadAlignment | None, error: str = "", ambiguous: bool = False)`
  - `@dataclass(frozen=True) CloneDifference(ref_pos: int, kind: str, ref_base: str, read_base: str, support: int, against: int, low_quality: int, status: str, in_region: bool, effect: str = "")` — `status ∈ {"confirmed", "disputed", "low_quality"}`
  - `@dataclass(frozen=True) CloneReport(clone: str, reads: tuple[ReadResult, ...], coverage: dict[int, tuple[int, int]], differences: tuple[CloneDifference, ...], mixed_positions: tuple[int, ...], verdict: str)` + свойства `region_coverage_pct: float`, `both_strands_pct: float`, `region_differences` (только `in_region`). `coverage[pos] = (ридов на «+», ридов на «−»)` с Q ≥ `min_q`, только для позиций региона
  - `load_read(path: Path) -> Chromatogram` — с учётом сессии правок (как `batch._edited_read`, вынести её в общую функцию и переиспользовать)
  - `verify_clone(construct: Construct, region: Region | None, reads: list[ReadRef], params: VerifyParams = VerifyParams(), annotate=None) -> CloneReport` — `annotate` — колбэк `(construct, CloneDifference) -> str` (эффект; подключается в Task 8)
  - `verify_all(construct, region, grouping: Grouping, params, progress=None, cancelled=None) -> list[CloneReport]` — `progress(done: int, total: int)`, `cancelled() -> bool` (при `True` — `CancelledError` из `concurrent.futures`)
- Правила (спецификация §2.4): `N`/IUPAC в **референсе** — не расхождение; расхождение «подтверждает» рид с Q ≥ `min_q` в этой позиции (для вставки — Q соседних баз), «против» — рид, покрывающий позицию с референсной базой и Q ≥ `min_q`; `alternatives > 1` → `ambiguous=True`, рид исключается из покрытия и расхождений. Вердикт — первое сработавшее: `no_data`, `mixed`, `errors` (есть `confirmed` или `disputed` в регионе), `incomplete`, `ok`.
- Produces (`tests/clone_fixtures.py`):
  - `mutate(seq: str, pos: int, kind: str, base: str = "") -> str` — `kind ∈ {"sub", "ins", "del"}`
  - `with_secondary_peak(chrom: Chromatogram, index: int, base: str, ratio: float) -> Chromatogram` — копирует форму пика вызванной базы в канал `base`, масштабируя к `ratio` × высота основного пика
  - `write_read(chrom: Chromatogram, folder: Path, name: str) -> Path` — сохраняет как сессию правок рядом с копией исходного `.ab1` (через `sessions.save_session`), чтобы `load_read` вернул изменённую хроматограмму

- [ ] **Step 1: Падающие тесты** — один набор клонов, собранный фикстурой из реальных ридов `ACHLO006-09` (F+R), регион = покрытая обоими ридами часть COI:

```python
def test_ok_clone(): assert report.verdict == "ok" and report.region_differences == ()

def test_substitution_found_with_support():
    ref = mutate(coi, P, "sub", other_base)        # P внутри региона
    r = verify_clone(Construct("c", ref, False), region, reads)
    d, = r.region_differences
    assert (d.ref_pos, d.kind, d.status, d.support, d.against) == (P, "substitution", "confirmed", 2, 0)
    assert r.verdict == "errors"

def test_insertion_and_deletion_found(): ...      # mutate "del" → read shows insertion; "ins" → deletion

def test_disputed_when_reads_disagree(): ...      # мутация в одном риде (через write_read) → support 1, against 1, "disputed", verdict "errors"

def test_outside_region_does_not_affect_verdict(): ...

def test_incomplete_when_read_trimmed(): ...      # только F, обрезан → "incomplete", region_coverage_pct < 100

def test_mixed_50_percent_flagged():               # 3 позиции с ratio 0.5
    assert r.verdict == "mixed" and len(r.mixed_positions) == 3

def test_mixed_30_percent_not_flagged(): ...      # ratio 0.3 → "ok"

def test_no_data_for_foreign_read(): ...          # рид pUC19-несовместимый (COI-рид против pUC19) → "no_data"

def test_reverse_complemented_session_same_result(): ...   # Review Focus 1

def test_region_across_origin(): ...              # повёрнутый кольцевой COI, регион start > end → те же расхождения  (Review Focus 2)

def test_read_without_quality_is_incomplete_not_errors(): ...  # все Q=0 → "incomplete" (Review Focus 3)

def test_reference_N_is_not_a_difference(): ...   # Review Focus 4

def test_unreadable_file_marked_not_fatal(): ...  # битый .ab1 → ReadResult.error, остальные обработаны

def test_ambiguous_read_excluded(): ...           # линейный конструкт = COI + COI → ReadResult.ambiguous, verdict "no_data"

def test_invalid_region_rejected(): ...           # пустой регион, вне конструкта, start > end на линейном → ValueError до выравнивания

def test_verify_all_progress_and_cancel(): ...

def test_performance_12_clones_x_4_reads():
    # конструкт ~10 т.п.н. = COI + случайный фланг (seed 1); 12 клонов × 4 рида из реальных
    t = time.perf_counter(); verify_all(...); elapsed = time.perf_counter() - t
    assert elapsed < LIMIT   # LIMIT = 2 × замер на машине разработчика, вписать числом
```

- [ ] **Step 2: Запустить** `python -m pytest tests/test_clone_verify.py -v` → FAIL.
- [ ] **Step 3: Реализовать** фикстуры в `tests/clone_fixtures.py` и проверку в `chrom_core/clones.py` по Interfaces.
- [ ] **Step 4: Запустить** → PASS; замерить `test_performance…`, вписать `LIMIT` числом.
- [ ] **Step 5: Запустить всё** `python -m pytest -q` → PASS.
- [ ] **Step 6: Commit** `git commit -m "clones: per-clone coverage, merged differences and verdicts"`

### Task 8: Эффект расхождения в CDS

**Files:**
- Create: `chrom_core/effects.py`
- Test: `tests/test_effects.py`

**Interfaces:**
- Consumes: `Construct`, `Feature` (Task 4), `CloneDifference` (Task 7)
- Produces:
  - `difference_effect(construct: Construct, diff: CloneDifference) -> str` — первая CDS, содержащая `diff.ref_pos` (с учётом обёртки и `strand`, рамка из `codon_start`). Замена: `"synonymous"` / `"<Ref><номер><Alt>"` (однобуквенные коды, например `"D45N"`) / `"<Ref><номер>*"` для стопа; вставка/делеция: `"frameshift"`, если длина не кратна 3, иначе `"inframe_ins"` / `"inframe_del"`. Вне CDS или без фич — `""`. Трансляция — `Bio.Seq.translate`, таблица 11.
  - Строка-код эффекта переводится в UI ключами `effect.synonymous`, `effect.frameshift`, `effect.inframe_ins`, `effect.inframe_del` (миссенс/стоп показываются как есть).

- [ ] **Step 1: Падающие тесты** — на COI из GenBank (с CDS) и на pUC19 (`bla` на минус-цепи):

```python
def test_synonymous_third_codon_position(): ...   # замена 3-го нуклеотида кодона на синонимичный
def test_missense_label(): assert difference_effect(coi_c, d) == f"{aa_ref}{n}{aa_alt}"
def test_stop_gained(): assert difference_effect(...).endswith("*")
def test_frameshift_and_inframe(): ...
def test_minus_strand_cds_puc19(): ...            # замена в bla даёт верный номер аминокислоты
def test_no_features_empty(): assert difference_effect(Construct("x", seq, False), d) == ""
```

  Ожидаемые значения вычисляются в тесте из известного кодона, а не хардкодятся «на глаз».
- [ ] **Step 2: Запустить** → FAIL.
- [ ] **Step 3: Реализовать**; подключить `difference_effect` как `annotate` по умолчанию в `verify_clone`.
- [ ] **Step 4: Запустить** `python -m pytest -q` → PASS.
- [ ] **Step 5: Commit** `git commit -m "effects: amino-acid consequences of clone differences in CDS"`

### Task 9: Отчёты по клонам

**Files:**
- Create: `chrom_core/clone_report.py`
- Test: `tests/test_clone_report.py`

**Interfaces:**
- Consumes: `CloneReport`, `CloneDifference` (Task 7)
- Produces:
  - `clones_csv(reports: list[CloneReport]) -> str` — колонки `clone, verdict, region_coverage_pct, both_strands_pct, differences, reads`
  - `differences_csv(reports: list[CloneReport]) -> str` — колонки `clone, position, kind, ref, read, effect, support, against, status, note`; `position` 1-based
  - `STATUS_HEX: dict[str, tuple[str, str]]` — `verdict -> (светлая, тёмная)` hex-цвета для `ok, errors, incomplete, mixed, no_data` (единый источник палитры; UI импортирует его в Task 13)
  - `clone_report_xlsx(path, reports: list[CloneReport], construct_name: str) -> None` — листы `Clones` и `Differences` с теми же колонками; ячейка вердикта залита светлым цветом из `STATUS_HEX`
  - CSV в том же стиле, что `verify.verdicts_csv` (разделитель `,`, экранирование как `batch._csv_field`)

- [ ] **Step 1: Падающие тесты** — строки заголовков, число строк, 1-based позиция, цвет заливки ячейки `errors` в XLSX (`openpyxl.load_workbook`).
- [ ] **Step 2: Запустить** → FAIL.
- [ ] **Step 3: Реализовать.**
- [ ] **Step 4: Запустить** → PASS.
- [ ] **Step 5: Commit** `git commit -m "clone_report: CSV and Excel reports for clone verification"`

### Task 10: Проект `.ab1proj`

**Files:**
- Create: `chrom_core/project.py`
- Test: `tests/test_project.py`

**Interfaces:**
- Consumes: `Construct`, `load_construct` (Task 4), `Grouping` (Task 6), `VerifyParams`, `Region` (Task 7)
- Produces:
  - `PROJECT_VERSION = 1`
  - `@dataclass Project(construct_path: Path, circular: bool, region: Region | None, grouping: Grouping, params: VerifyParams)`
  - `save_project(path: Path, project: Project) -> None` — JSON UTF-8 (`ensure_ascii=False`); пути к конструкту и ридам относительные к папке проекта, если лежат внутри неё, иначе абсолютные
  - `load_project(path: Path) -> Project` — `ValueError` с понятным сообщением при битом JSON, чужом файле (`app != "ab1kit"`) или `version > PROJECT_VERSION`

- [ ] **Step 1: Падающие тесты:** round-trip; относительные пути при переносе папки проекта; путь `tmp_path / "Данные" / "клоны 1"` (Review Focus 5); битый файл и будущая версия → `ValueError`.
- [ ] **Step 2: Запустить** → FAIL.
- [ ] **Step 3: Реализовать.**
- [ ] **Step 4: Запустить** `python -m pytest -q` → PASS.
- [ ] **Step 5: Commit** `git commit -m "project: .ab1proj files for clone verification sessions"`

---

## Фаза 3 — интерфейс

### Task 11: Фоновые задачи

**Files:**
- Create: `ui/tasks.py`
- Modify: `ui/main_window.py` (`_scan_folder`, `_blast_submit` — через задачи)
- Test: `tests/test_ui_tasks.py`

**Interfaces:**
- Produces:
  - `class Task(QObject)` с сигналами `progress(int, int)`, `finished(object)`, `failed(str)`, `cancelled()`; методы `cancel()`, `is_cancelled() -> bool`
  - `run_task(fn: Callable[[Callable[[int, int], None], Callable[[], bool]], object], parent: QObject) -> Task` — `fn(progress, is_cancelled)` выполняется в `QThreadPool.globalInstance()`; сигналы приходят в GUI-поток; исключение → `failed(str(exc))`; `CancelledError` → `cancelled`
  - В строке состояния `MainWindow` — `QProgressBar` и кнопка «Отмена» (ключи `task.cancel`, `task.running`), видимые во время задачи

- [ ] **Step 1: Падающие тесты** (offscreen, `qtbot`-free: крутить `QApplication.processEvents()` до сигнала с таймаутом 5 с): результат приходит в `finished`; исключение — в `failed`; отмена — в `cancelled`; GUI-поток не блокируется (таймер в GUI-потоке тикает, пока задача спит 0,3 с).
- [ ] **Step 2: Запустить** → FAIL.
- [ ] **Step 3: Реализовать** `ui/tasks.py`; перевести сканирование папки и BLAST на `run_task`.
- [ ] **Step 4: Запустить** `python -m pytest -q` → PASS.
- [ ] **Step 5: Commit** `git commit -m "ui: background tasks with progress and cancel"`

### Task 12: Разбить `main_window.py`

**Files:**
- Create: `ui/read_view/__init__.py`, `ui/read_view/dialogs.py` (`SearchDialog`, `TrimDialog`), `ui/read_view/analysis_dock.py`, `ui/read_view/verify_dock.py`, `ui/read_view/batch_dock.py` (+ перенос `ui/batch_panel.py`), `ui/read_view/read_view.py` (`ReadView(QMainWindow)` — всё, что сейчас относится к одному риду: вью, миникарта, доки, правка, обрезка, экспорт, поиск)
- Modify: `ui/main_window.py` — остаётся оболочка: меню «Файл/Вид/Справка», язык, тема, недавние файлы, `run()`
- Test: существующий `tests/test_smoke_gui.py` (без изменения утверждений; правятся только пути monkeypatch)

**Interfaces:**
- Produces:
  - `ReadView.open_path(path) -> bool`, `ReadView.refresh_read()`, `ReadView.current_model() -> EditableRead | None`, сигнал `ReadView.read_edited(Path)` (испускается после правки, undo/redo, обрезки и сохранения сессии)
  - `MainWindow.read_view: ReadView`; публичные методы `MainWindow`, которые используют тесты (`open_path`, `apply_base_*`, `undo_read`, …), остаются как делегаты, чтобы `test_smoke_gui.py` не менялся по смыслу

- [ ] **Step 1: Перенести код** без изменения поведения, по одному блоку (диалоги → доки → ReadView), после каждого — `python -m pytest tests/test_smoke_gui.py -q` → PASS.
- [ ] **Step 2: Проверить размер:** ни один файл в `ui/` не длиннее ~450 строк (`wc -l ui/*.py ui/read_view/*.py`).
- [ ] **Step 3: Запустить всё** `python -m pytest -q` → PASS.
- [ ] **Step 4: Commit** `git commit -m "ui: split main window into read_view modules (no behaviour change)"`

### Task 13: Тема, иконки, палитра статусов

**Files:**
- Create: `ui/style.qss` (светлая), `ui/style_dark.qss`, `ui/icons/*.svg` + `ui/icons/LICENSE` (набор Lucide, ISC/MIT — проверить лицензию при скачивании; иначе Tabler Icons, MIT), `ui/status.py`
- Modify: `ui/theme.py` (загрузка QSS), `ui/read_view/*` (иконки на тулбаре, растяжение колонок, тултипы длинных ячеек), `ui/chrom_view.py` (цвет полос качества по Q), `ab1kit.spec` (включить `ui/style*.qss`, `ui/icons`), `pyproject.toml` (`package-data` для `ui`)
- Test: `tests/test_smoke_gui.py`

**Interfaces:**
- Produces (`ui/status.py`):
  - `status_color(verdict: str, theme: str) -> QColor` — из `chrom_core.clone_report.STATUS_HEX` (Task 9)
  - `quality_color(q: int, theme: str) -> QColor` — `<20` красный, `20–29` жёлтый, `≥30` зелёный (те же пороги, что миникарта)
  - `status_chip(verdict: str, text: str) -> QLabel` — скруглённая плашка
  - `icon(name: str) -> QIcon`

- [ ] **Step 1: Падающие тесты:** у каждого действия тулбара есть иконка; `quality_color(15) != quality_color(35)`; таблицы растягивают последнюю колонку и имеют тултип с полным текстом примечания; QSS применяется (`app.styleSheet()` не пуст) в обеих темах.
- [ ] **Step 2: Запустить** → FAIL.
- [ ] **Step 3: Реализовать.** Пакет иконок: только используемые (~20 файлов).
- [ ] **Step 4: Запустить** `python -m pytest -q` → PASS; сделать скриншоты обеих тем через `scripts/render_preview.py` и посмотреть глазами.
- [ ] **Step 5: Commit** `git commit -m "ui: QSS theme, icons, unified status palette, coloured quality bars"`

### Task 14: Строка-сводка и матрица клонов

**Files:**
- Create: `ui/construct/__init__.py`, `ui/construct/summary_bar.py`, `ui/construct/matrix.py`
- Test: `tests/test_construct_ui.py`

**Interfaces:**
- Consumes: `Construct`, `Region`, `CloneReport` (Tasks 4, 7), `status_color`, `status_chip` (Task 13)
- Produces:
  - `SummaryBar(QWidget)`: `set_data(construct: Construct, region: Region | None, reports: list[CloneReport])`; сигнал `filter_changed(set[str])` — множество показываемых вердиктов (клик по счётчику переключает)
  - `CloneMatrix(QWidget)` (собственная отрисовка `paintEvent`): `set_data(construct, region, reports)`, `set_filter(verdicts: set[str])`; дорожка референса (фичи прямоугольниками с подписью, рамка региона), строка на клон: покрытие (две цепи — насыщенный, одна — светлый, Q<20/нет — фон), метки расхождений цветом статуса, плашка вердикта справа. Колёсико — масштаб общей оси, перетаскивание — сдвиг. Протяжка по дорожке референса с зажатым Shift — задать регион.
  - Сигналы: `difference_clicked(str, int)` (клон, `ref_pos`), `clone_clicked(str)`, `region_changed(object)` (`Region`)
  - `CloneMatrix.hit_test(pos: QPoint) -> tuple[str, int] | None` — для тестов и кликов

- [ ] **Step 1: Падающие тесты** (offscreen, данные — `verify_all` на наборе из Task 7): сводка показывает счётчики по вердиктам; фильтр скрывает строки; `hit_test` по координате метки возвращает `(clone, ref_pos)`; Shift-протяжка испускает `region_changed`; виджет рендерится в `QImage` без исключений в обеих темах; 12 строк помещаются в 600 px по высоте.
- [ ] **Step 2: Запустить** → FAIL.
- [ ] **Step 3: Реализовать.**
- [ ] **Step 4: Запустить** → PASS; отрендерить PNG матрицы и проверить глазами против макета B.
- [ ] **Step 5: Commit** `git commit -m "ui: construct summary bar and clone matrix"`

### Task 15: Панель деталей расхождения

**Files:**
- Create: `ui/construct/detail_pane.py`
- Test: `tests/test_construct_ui.py`

**Interfaces:**
- Consumes: `CloneReport`, `ReadResult`, `load_read` (Task 7), `ChromatogramWidget` (`ui/chrom_view.py`: `set_chromatogram`, `set_reference`, `goto_base`, `set_highlights`)
- Produces:
  - `DetailPane(QWidget)`: `show_clone(report: CloneReport, construct: Construct)`, `show_position(ref_pos: int)`; слева — стопка `ChromatogramWidget` по одному на рид клона (компактная высота, общий масштаб), каждый центрирован на базе, соответствующей `ref_pos` (через `alignment.ref_to_read` и ориентацию), с подсветкой; над стопкой — строка референса ±10 п.н.; справа — таблица расхождений клона (позиция 1-based, `ref→read`, эффект через `tr`, за/против, статус); выбор строки → `show_position`
  - Сигнал `open_read(Path)` — по двойному клику на хроматограмму

- [ ] **Step 1: Падающие тесты:** число хроматограмм = числу ридов клона с выравниванием; после `show_position(P)` центральная база каждого вью — та, что выровнена на `P`; для рида на минус-цепи — тоже; выбор строки таблицы вызывает `show_position`; двойной клик испускает `open_read`.
- [ ] **Step 2: Запустить** → FAIL.
- [ ] **Step 3: Реализовать.**
- [ ] **Step 4: Запустить** → PASS.
- [ ] **Step 5: Commit** `git commit -m "ui: difference detail pane with stacked chromatograms"`

### Task 16: Экран «Проверка конструкта», группировка, связка с просмотрщиком

**Files:**
- Create: `ui/construct/grouping_dialog.py`, `ui/construct/workspace.py`
- Modify: `ui/main_window.py` (вкладки «Проверка конструкта» / «Просмотр рида», меню «Проект»: новый/открыть/сохранить, экспорт отчёта), `i18n.py`, `translations_en.py`
- Test: `tests/test_construct_ui.py`, `tests/test_smoke_gui.py`

**Interfaces:**
- Consumes: всё из Tasks 4–15; `run_task` (Task 11); `ReadView.read_edited` (Task 12)
- Produces:
  - `GroupingDialog(QDialog)`: `__init__(grouping: Grouping, pattern: str, parent)`, `result_grouping() -> Grouping`, `result_pattern() -> str`; таблица «файл / клон / праймер» с правкой; поле шаблона с предпросмотром (невалидный regex — красная рамка, кнопка OK неактивна); строки «без клона» выделены
  - `ConstructWorkspace(QWidget)`: `set_construct(path)`, `add_reads(paths)`, `set_region(region)`, `load_project(path)`, `save_project(path)`, `export_report(path)`; пустое состояние с двумя зонами перетаскивания (конструкт; папка/файлы ридов); любое изменение (конструкт, риды, группировка, регион, параметры, `read_edited`) → перезапуск `verify_all` через `run_task` с отменой предыдущего; сигнал `open_read(Path)`
  - Поля параметров `min_q`, `het_ratio`, `mixed_min_positions` и галочка «кольцевая» — в строке над матрицей; сохраняются в проекте
  - `MainWindow`: вкладка «Проверка конструкта» первая и активна при старте без аргументов; `open_read` → переключение на «Просмотр рида» + `open_path`; `ab1kit file.ab1proj` открывает проект

- [ ] **Step 1: Падающие тесты:** перетаскивание COI `.gb` и 8 реальных ридов → после завершения задачи 4 строки в матрице и вердикты как в Task 7 (для неизменённого референса — все `ok` или `incomplete`); правка рида в «Просмотре рида» (замена базы в позиции расхождения) → пересчёт и смена вердикта; сохранение/открытие проекта восстанавливает регион и группировку; битый `.ab1proj` → сообщение и пустое рабочее пространство без исключения; недопустимый регион → сообщение, проверка не запускается; отмена пересчёта оставляет прежние вердикты в матрице; `GroupingDialog` с невалидным regex блокирует OK; пустое состояние видно без данных; все новые ключи есть в RU и EN.
- [ ] **Step 2: Запустить** → FAIL.
- [ ] **Step 3: Реализовать.**
- [ ] **Step 4: Запустить** `python -m pytest -q` → PASS.
- [ ] **Step 5: Ручная проверка:** запустить `ab1kit`, пройти сценарий из критерия успеха спецификации (перетащить конструкт и папку, задать регион, открыть расхождение, перейти в рид, поправить, вернуться) в светлой и тёмной теме, RU и EN. Обновить скриншоты `docs/screenshots/` (`scripts/render_preview.py`), добавить скриншот экрана проверки.
- [ ] **Step 6: Commit** `git commit -m "ui: construct verification workspace with grouping and project files"`

---

## Фаза 4 — CLI и выпуск

### Task 17: Подкоманды CLI, документация, версия

**Files:**
- Modify: `chrom_cli.py`, `README.md`, `CHANGELOG.md`, `pyproject.toml` (`version = "1.1.0"`)
- Test: `tests/test_cli.py`

**Interfaces:**
- Consumes: `load_construct`, `group_reads`, `verify_all`, `VerifyParams`, `clones_csv`/`differences_csv`/`clone_report_xlsx`, `batch.scan_folder`/`summary_csv`/`batch_export`
- Produces:
  - `main(argv) -> int`: если `argv[0] in {"verify-clones", "qc", "export"}` — `argparse`-подкоманда, иначе GUI (как сейчас)
  - `verify-clones --construct FILE [--region A-B] [--linear] [--pattern REGEX] [--min-q N] [--het-ratio X] READS_DIR -o REPORT` — `A-B` 1-based включительно → `Region(A-1, B)`; `.xlsx` → `clone_report_xlsx`, `.csv` → `clones_csv` в `REPORT` и `differences_csv` в `REPORT` с суффиксом `_differences.csv`; печать короткой таблицы клонов в stdout; код возврата 0 — есть хотя бы один `ok`, 1 — нет, 2 — ошибка входных данных (сообщение в stderr)
  - `qc READS_DIR [-o CSV]` — `scan_folder` + `summary_csv` (stdout без `-o`)
  - `export READS_DIR --trim [--fasta|--fastq] -o DIR` — `batch_export`

- [ ] **Step 1: Падающие тесты:** `verify-clones` на `tmp_path / "Данные клоны"` с копиями реальных ридов и COI → код 0, XLSX с листами `Clones`/`Differences`; на изменённом референсе — строка с `errors`; без ридов/битый конструкт → 2; `qc` печатает CSV; `export` создаёт файлы; `main([])` и `main(["x.ab1"])` по-прежнему вызывают GUI (`monkeypatch` на `ui.main_window.run`).
- [ ] **Step 2: Запустить** → FAIL.
- [ ] **Step 3: Реализовать.**
- [ ] **Step 4: Документация:** в README — раздел «Clone verification» (сценарий, формат имён, шаблон, вердикты, CLI-пример) и новый скриншот; CHANGELOG — `[1.1.0]` с Added/Changed/Fixed по задачам 2–17.
- [ ] **Step 5: Запустить всё** `python -m pytest -q` → PASS; собрать exe `pyinstaller ab1kit.spec` и запустить — экран проверки открывается, иконки и тема на месте.
- [ ] **Step 6: Commit** `git commit -m "CLI subcommands verify-clones/qc/export; docs and v1.1.0"`
