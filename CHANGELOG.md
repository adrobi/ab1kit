# Changelog

All notable changes to ab1kit are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[SemVer](https://semver.org/).

## [1.0.0] — 2026-10-24

First feature-complete release: the full cycle «open a chromatogram → inspect →
edit → trim → export → align to a reference → verify expected genotypes».

### Added
- ABIF (`.ab1`) parsing: traces, phred quality, peak positions, metadata;
  dye-channel mapping verified empirically on ABI 310/3730 data
- Chromatogram viewer: zoom/pan, adaptive ruler, minimap with quality map and
  optional quality profile, dark theme, per-channel toggles, smoothing,
  IUPAC-aware subsequence search, go-to-base, drag&drop, recent files, file
  info panel
- Base editing (change/insert/delete) with undo/redo, edited-base markers,
  reverse-complement as one undoable operation
- Quality trimming: sliding window with live preview, draggable handles,
  editable/undoable application
- Export: FASTA, FASTQ (phred33), PNG, PDF; editing sessions autosaved as
  `<name>.ab1kit.json` and restored on reopen
- Batch mode: folder QC scan table, batch trim+export honoring sessions,
  CSV summary
- Reference alignment: exact semiglobal DP with affine gaps via
  Bio.Align.PairwiseAligner, both orientations, ambiguity detection
- SNP calling: homozygous / heterozygous (double-peak ratio + IUPAC) /
  low-quality verdicts
- Expected-SNP verification: VCF/CSV/TSV input, per-position verdicts with
  input-error guards, chromatogram marks, CSV and Excel reports
- NCBI BLAST submission (browser results)
- RU/EN interface, settings persisted in JSON, CI on Python 3.11–3.13,
  one-file Windows executable via PyInstaller

### Performance
- Min/max pixel-bucket decimation of traces at overview zoom: ~7× faster
  rendering on 3,000+ base reads (0.5 s → 73 ms per frame)

### Notes
- `.scf` support intentionally deferred: modern BioPython dropped its SCF
  reader and no reference sample files are available to validate a custom
  parser against
