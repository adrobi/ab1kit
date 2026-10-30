# Changelog

All notable changes to ab1kit are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versioning is
[SemVer](https://semver.org/).

## [1.1.0] — 2026-10-30

Clone verification: check a plate of clones against a construct, in the GUI or
from the command line.

### Added
- Reference loading from FASTA, GenBank and SnapGene files, with circular
  constructs and reads spanning the origin
- Normalised secondary-peak analysis at any read position
- Grouping of reads into clones by file name (default rules or a custom regex),
  per-clone coverage by strand, merged differences and verdicts
  (`ok` / `errors` / `incomplete` / `mixed` / `no_data`)
- Amino-acid consequences of differences inside CDS features
- Clone report as CSV and Excel (sheets `Clones` and `Differences`);
  `.ab1proj` project files for verification sessions
- Clone-verification workspace: summary bar, clone matrix, difference detail pane with stacked
  chromatograms, grouping editor, drag-and-drop of reads and projects
- Background tasks with progress and cancel
- Command-line subcommands `verify-clones`, `qc` and `export` (no Qt import)
- Real-data test corpus (sangeranalyseR reads, COI and pUC19 references)

### Changed
- Refreshed UI: QSS theme, icons, one status palette everywhere
- Quality colour thresholds are now <20 / 20-29 / >=30 everywhere
- Main window split into the `ui.read_view` package (no behaviour change)
- Alignment counts distinct placements; repeated hits are flagged as ambiguous

### Fixed
- Reverse-strand reads in expected-SNP verification and the reference overlay
  were indexed in the wrong orientation (now correct)
- Reverse-complementing a read misplaced peak positions; editing sessions saved
  after a reverse-complement with 1.0.0 or earlier hold wrong peak positions —
  re-open the original file and reverse-complement again
- Untranslated verification notes; mean quality for reads without quality data

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
