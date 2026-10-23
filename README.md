# ab1kit

Open-source desktop tool for viewing and analyzing Sanger sequencing chromatograms
(`.ab1` files) — a free, offline alternative to FinchTV / Chromas for everyday
molecular biology lab work. English/Russian interface.

## Features

**MVP (in development)**

- Parse ABIF (`.ab1`): traces, phred quality, peak positions, base calls
- Chromatogram viewer: smooth zoom/scroll, quality bar, position ruler,
  subsequence search
- Read editing: manual base-call corrections, reverse-complement,
  quality trimming (sliding window, Chromas-style)
- Export: FASTA / FASTQ (with quality), copy sequence, save trimmed reads
- Batch mode: open a folder, per-file QC, batch trim & export

**Planned (v1.1+)**

- Reference alignment (built-in Smith-Waterman, optional MAFFT)
- SNP calling, including heterozygous double peaks (IUPAC codes)
- Expected-SNP verification: `.ab1` + list of expected variants →
  confirmed genotypes and mismatches highlighted
- Excel/CSV reports

## Stack

Python 3.11+, PySide6 (Qt) GUI. Core logic is UI-independent and covered by
pytest. MIT license.

## Getting started

Requires Python 3.11+:

```
python -m venv .venv
.venv/Scripts/pip install -e ".[dev]"    # Windows; on POSIX: .venv/bin/pip
.venv/Scripts/python -m pytest           # run the test suite
.venv/Scripts/ab1kit                     # launch the GUI
```

## Status

Concept stage — the full concept and roadmap (in Russian) are in
[concept.ru.md](concept.ru.md).
