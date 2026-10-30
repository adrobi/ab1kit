# ab1kit

Open-source desktop tool for viewing, editing and analyzing Sanger sequencing
chromatograms (`.ab1` files) — a free, offline alternative to FinchTV / Chromas
for everyday molecular biology lab work. English/Russian interface.

![viewer](docs/screenshots/viewer.png)
![construct check](docs/screenshots/construct.png)

## Features

**Viewer**

- Parse ABIF (`.ab1`): traces, phred quality, peak positions, base calls
- Smooth chromatogram viewer: zoom/scroll (wheel, keyboard, Ctrl+drag to zoom
  into a selection), position ruler, quality bars, per-base tooltips
- Whole-read minimap with quality map and optional quality profile
- Dark theme, per-channel (A/C/G/T) visibility, trace smoothing
- Find subsequence (Ctrl+F, IUPAC-aware, both strands), go to base (Ctrl+G)

**Editing & trimming**

- Manual base editing (change / insert / delete) with unlimited undo/redo
- Reverse-complement as a single undoable operation
- Quality trimming (sliding window, Chromas-style) with live preview,
  draggable trim handles
- FASTA / FASTQ (phred33) / PNG / PDF export; editing sessions are saved and
  restored automatically

**Batch mode**

- Scan a folder of `.ab1` files: per-file QC (length, mean Q, Q20, trim
  length), status marks, sorting, hiding bad reads
- Batch trim + export respecting per-file editing sessions; CSV summary

**Analysis**

- Reference alignment: exact semiglobal DP with affine gaps
  (Bio.Align.PairwiseAligner, free end gaps — read ends are never dropped),
  both read orientations, ambiguity flag for repetitive references
- SNP calling: homozygous vs heterozygous (double-peak detection by
  secondary/primary height ratio, IUPAC codes), low-quality flagging
- **Expected-SNP verification**: load a VCF/CSV/TSV of expected variants and
  get a verdict per position (confirmed / heterozygous / not found / no data /
  low quality / input-error guards), with marks on the chromatogram and
  CSV / Excel reports
- Send the current read to NCBI BLAST (opens results in the browser)

![dark theme](docs/screenshots/viewer-dark.png)
![batch mode](docs/screenshots/batch.png)

## Clone verification

Check a set of plasmid or amplicon clones against the expected construct from
their Sanger reads: open the clone-verification workspace, load the reference (FASTA,
GenBank or SnapGene), drop the `.ab1` reads of the whole plate, optionally set
the region of interest, and read the verdict per clone off the colour matrix.

![clone verification](docs/screenshots/construct.png)

**Scenario.** Load the construct, add the reads (files or a folder), check the
automatic grouping into clones, choose the region; every clone gets its verdict right away. Click a
difference in the matrix to see the stacked chromatograms of every read at that
position. Export the result as Excel/CSV, or save the whole session as an
`.ab1proj` file. FASTA constructs are treated as circular (reads may span the
origin); GenBank/SnapGene files take the topology from the file; the *Circular*
checkbox (`--linear` on the command line) overrides the topology.

**File names.** Reads are grouped by name as `<clone>_<primer>`: the primer is
the last known primer token (`F`, `R`, `T7`, `SP6`, `M13F`, `pUC-F`, ...); if
the name has none, the last token is taken as the primer. The clone is what
precedes the primer, minus a single-digit order token and a leading plate well
(`A05`), the well being stripped only when another clone token follows it
(`B2_F` is clone `B2`). `A05_clone3_2_R.ab1` is clone `clone3`, primer `R`;
`Achl_ACHLO006-09_1_F.ab1` is clone `Achl_ACHLO006-09`. Files whose name yields
no clone stay unassigned and are skipped. For other conventions give a regular
expression with named groups `clone` and `primer` that matches the whole file
name without the extension, e.g. `(?P<clone>\d+)-(?P<primer>[A-Za-z0-9]+)`.

**Verdicts.** Checked in this order, the first match wins:

| Verdict | Meaning |
| --- | --- |
| `no_data` (grey) | no read of the clone aligned to the construct |
| `mixed` (purple) | at least 3 positions (default) with a double peak: probably a mixed colony |
| `errors` (red) | at least one difference in the region that is confirmed or disputed |
| `incomplete` (yellow) | no such difference, but part of the region is not covered by any read |
| `ok` (green) | the region is fully covered and shows no reliable difference |

A difference is *confirmed* (only supporting reads, Q >= min Q), *disputed*
(another read shows the reference there) or *low quality* (seen only below
min Q). Both confirmed and disputed differences give `errors`; only low-quality
ones do not change the verdict. Defaults: min Q 20, double-peak ratio 0.40.

**Command line** (no GUI, no Qt needed):

```
ab1kit verify-clones --construct pET28-gene.gb --region 120-1850 reads/ -o report.xlsx
ab1kit qc reads/ -o qc.csv
ab1kit export reads/ --trim --fasta -o trimmed/
```

`verify-clones` prints a table of clones and exits with 0 if at least one clone
is `ok`, 1 if none is, 2 on invalid input; `-o report.csv` writes `report.csv`
and `report_differences.csv`. Options: `--region A-B` (1-based, inclusive),
`--linear`, `--pattern REGEX`, `--min-q N`, `--het-ratio X`.

The exit-code contract applies to the `ab1kit` command installed with pip. The
Windows release `ab1kit.exe` is a windowed program: the shell returns to the
prompt immediately, so exit codes and output ordering are unreliable and output
may appear after the prompt. Wait for it explicitly: `start /wait ab1kit.exe
verify-clones ...` (cmd), or `Start-Process -Wait -NoNewWindow` /
`(Start-Process ab1kit.exe -ArgumentList ... -Wait -PassThru).ExitCode`
(PowerShell).

## Getting started

Requires Python 3.11+:

```
python -m venv .venv
.venv/Scripts/pip install -e ".[dev]"    # Windows; on POSIX: .venv/bin/pip
.venv/Scripts/python -m pytest           # run the test suite
.venv/Scripts/ab1kit                     # launch the GUI
```

`ab1kit` also accepts a file or folder argument. On Windows, a prebuilt
executable is available on the [Releases](https://github.com/adrobi/ab1kit/releases)
page.

## Keyboard shortcuts

| Shortcut | Action |
| --- | --- |
| Ctrl+O / Ctrl+Shift+O | open file / folder |
| wheel, Ctrl+drag | zoom / zoom to selection |
| ←→ PgUp/PgDn Home/End | navigate |
| + / − | zoom in / out |
| Ctrl+G | go to base |
| Ctrl+F / F3 | find subsequence / next match |
| Ctrl+Z / Ctrl+Y | undo / redo |
| Ctrl+R | reverse-complement |
| Ctrl+T | quality trimming |
| Ctrl+S | save editing session |
| Ctrl+L | load reference |
| Ctrl+Shift+V | check expected SNPs |

## Stack

Python 3.11+, PySide6 (Qt) GUI, BioPython for ABIF parsing and alignment.
Core logic (`chrom_core`) is UI-independent and covered by pytest. MIT
license — see [LICENSE](LICENSE).

## Contributing

Issues and pull requests are welcome at
[github.com/adrobi/ab1kit](https://github.com/adrobi/ab1kit).
Example chromatograms in `examples/` come from the
[BioPython test suite](https://github.com/biopython/biopython/tree/master/Tests/Abi)
(Biopython License).

## Status

v1.1.0 — adds clone verification, a command-line interface, a refreshed UI
and more accurate peak analysis on top of the 1.0 viewer, editing, batch mode,
alignment, SNP calling and expected-SNP verification. See [CHANGELOG.md](CHANGELOG.md).
