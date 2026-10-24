# Test data corpus

Downloaded once by `scripts/fetch_test_data.py --coi HM417896` and committed.
Tests never use the network. Semi-synthetic data is generated in
`tests/clone_fixtures.py` and is not stored on disk.

| File | Source | Pin | License |
|---|---|---|---|
| `real/Achl_ACHLO{006,007,040,041}-09_{1_F,2_R}.ab1` (8 files) | [sangeranalyseR](https://github.com/roblanf/sangeranalyseR) `inst/extdata/Allolobophora_chlorotica/ACHLO/` | commit `4c5ebc9298a873c209819f12b7dfd85e126f1406` | MIT |
| `refs/pUC19_L09137.gb` | NCBI nuccore `L09137` (efetch, `gbwithparts`) | accession | public domain |
| `refs/coi_HM417896.gb` | NCBI nuccore `HM417896.1` (efetch, `gbwithparts`) | accession | public domain |
| `refs/sample-d.dna` | [Biopython](https://github.com/biopython/biopython) `Tests/SnapGene/sample-d.dna` | commit `2fbe8755bb44731e60b7bd1e2cc55cf34d94e15b` | Biopython License |

## Why HM417896 as the COI reference

*Allolobophora chlorotica* is a cryptic-species aggregate: its mitochondrial COI
lineages differ by up to ~16%. The four samples fall into at least three
lineages, so no single GenBank record is >= 97% identical to all eight reads.

Candidates were found with NCBI E-utilities (esearch on `Allolobophora chlorotica`
COI/COX1, ~930 records) and checked locally with `chrom_core.align.align_read`.
HM417896.1 (voucher AC2, iBOL/BOLD mirror, 658 bp, annotated CDS) is the record
that fits the most reads and is the one that matches the ACHLO006-09 / 007-09
lineage:

| Read | strand | identity vs HM417896 |
|---|---|---|
| ACHLO006-09 F / R | + / - | 0.995 / 0.973 |
| ACHLO007-09 F / R | + / - | 0.995 / 0.992 |
| ACHLO040-09 F / R | + / - | 0.822 / 0.821 |
| ACHLO041-09 F / R | + / - | 0.836 / 0.818 |

ACHLO040-09 and ACHLO041-09 are other lineages: their best records in GenBank are
KC569719.1 (040 F: 0.986) and HM417923.1 (041 F: 0.997). They are kept in the
corpus on purpose as real, divergent-from-reference reads (useful as a
"many real differences" case); the corpus test treats them separately.
