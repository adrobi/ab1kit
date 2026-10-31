# Test data corpus

Downloaded once by `scripts/fetch_test_data.py --coi HM417896 --coi KC569719 --coi HM417923` and committed.
Tests never use the network. Semi-synthetic data is generated in
`tests/clone_fixtures.py` and is not stored on disk.

| File | Source | Pin | License |
|---|---|---|---|
| `real/Achl_ACHLO{006,007,040,041}-09_{1_F,2_R}.ab1` (8 files) | [sangeranalyseR](https://github.com/roblanf/sangeranalyseR) `inst/extdata/Allolobophora_chlorotica/ACHLO/` | commit `4c5ebc9298a873c209819f12b7dfd85e126f1406` | MIT |
| `refs/pUC19_L09137.gb` | NCBI nuccore `L09137` (efetch, `gbwithparts`) | accession | public domain |
| `refs/coi_{HM417896,KC569719,HM417923}.gb` | NCBI nuccore (efetch, `gbwithparts`) | accessions | public domain |
| `refs/sample-d.dna` | [Biopython](https://github.com/biopython/biopython) `Tests/SnapGene/sample-d.dna` | commit `2fbe8755bb44731e60b7bd1e2cc55cf34d94e15b` | Biopython License |

## COI references (per-lineage)

*Allolobophora chlorotica* is a cryptic-species aggregate: its mitochondrial COI
lineages differ by up to ~16%, so no single GenBank record is >= 97% identical to
all eight reads. Candidates were found with NCBI E-utilities (esearch, ~930
records) and checked locally with `chrom_core.align.align_read`. Each record
below is an iBOL/BOLD mirror with an annotated CDS; `SAMPLE_REFERENCE` in
`tests/clone_fixtures.py` maps each sample to its reference.
`coi_reference()` is HM417896, the primary reference for ACHLO006/007.

| Accession | Samples | Identity F / R |
|---|---|---|
| `HM417896.1` (voucher AC2, 658 bp) | ACHLO006-09 | 0.995 / 0.973 |
| | ACHLO007-09 | 0.995 / 0.992 |
| `KC569719.1` (voucher Acy224, 625 bp) | ACHLO040-09 | 0.986 / 0.977 |
| `HM417923.1` (voucher B2, 658 bp) | ACHLO041-09 | 0.997 / 0.989 |

Cross-lineage: ACHLO040-09 F against HM417896 is only 0.822 identical (041:
0.836; 006/007 against KC569719/HM417923: 0.82-0.84), which makes it a real
"different lineage" read for later tests.
