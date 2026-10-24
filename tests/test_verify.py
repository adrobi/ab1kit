import random
from pathlib import Path

import pytest

from chrom_core import read_ab1
from chrom_core.align import align_read
from chrom_core.verify import (
    ExpectedVariant,
    parse_expected,
    verify_expected,
    verdicts_csv,
)
from clone_fixtures import make_chrom

EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


# ------------------------------------------------------------------- parsing
def test_parse_vcf() -> None:
    text = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t101\trs1\tC\tT\t.\t.\t.\n"
        "chr1\t205\trs2\tG\tA\t.\t.\t.\n"
    )
    variants, warnings = parse_expected(text)
    assert warnings == []
    assert variants == [
        ExpectedVariant(101, "C", "T", "rs1"),
        ExpectedVariant(205, "G", "A", "rs2"),
    ]


def test_parse_csv_with_header() -> None:
    text = "name;pos;ref;alt\nsnpA;50;C;T\nsnpB;80;G;A\n"
    variants, warnings = parse_expected(text)
    assert warnings == []
    assert variants[0].position == 50 and variants[0].name == "snpA"
    assert variants[0].ref == "C" and variants[0].alt == "T"


def test_parse_tsv_without_header() -> None:
    variants, warnings = parse_expected("101\tC\tT\n205\tG\tA\n")
    assert warnings == []
    assert [v.position for v in variants] == [101, 205]


def test_parse_two_columns_position_alt() -> None:
    variants, _ = parse_expected("101,T\n")
    assert variants == [ExpectedVariant(101, ".", "T", "")]


def test_parse_bad_rows_warn() -> None:
    variants, warnings = parse_expected("101,C,T\nnotanumber,C,G\n")
    assert len(variants) == 1
    assert len(warnings) == 1


def test_name_first_column_with_header() -> None:
    # regression: a header line must be detected even when column 0 holds
    # names, and the name column must be honored
    variants, warnings = parse_expected(
        "name;pos;ref;alt\nsnp1;51;A;T\nsnp2;61;C;G\n"
    )
    assert warnings == []
    assert [v.position for v in variants] == [51, 61]
    assert variants[0].name == "snp1"
    assert variants[0].ref == "A" and variants[0].alt == "T"


def test_name_first_column_without_header_warns_not_swallows() -> None:
    # without a header the parser assumes (pos, ref, alt, ...) — a names-first
    # row must produce a warning, not silently vanish
    variants, warnings = parse_expected("snp1;51;A;T\n201;C;G\n")
    assert [v.position for v in variants] == [201]
    assert len(warnings) == 1


def test_parse_empty_raises() -> None:
    with pytest.raises(ValueError):
        parse_expected("   \n")


# --------------------------------------------------------------- verification
@pytest.fixture
def setup():
    rng = random.Random(21)
    reference = "".join(rng.choice("ACGT") for _ in range(120))
    return reference


def confirmed_case(reference, index, alt):
    read = list(reference)
    read[index] = alt
    return "".join(read)


def _verify(reference, read, variants, chrom=None, min_q=20):
    chrom = chrom or make_chrom(read)
    alignment = align_read(read, reference)
    return verify_expected(alignment, chrom, reference, variants, min_q=min_q), chrom


def test_confirmed_homozygous(setup) -> None:
    reference = setup
    read = confirmed_case(reference, 50, "T" if reference[50] != "T" else "A")
    verdicts, _ = _verify(reference, read, [ExpectedVariant(51, reference[50], read[50], "snp1")])
    assert verdicts[0].verdict == "confirmed"
    assert verdicts[0].observed == read[50]
    assert verdicts[0].read_pos == 50
    assert verdicts[0].quality == 40


def test_confirmed_heterozygous(setup) -> None:
    reference = setup
    alt = "T" if reference[50] != "T" else "A"
    read = confirmed_case(reference, 50, alt)
    chrom = make_chrom(read, heights={(50, "C"): 600})  # second peak 60%
    verdicts, _ = _verify(reference, read, [ExpectedVariant(51, reference[50], alt)], chrom=chrom)
    assert verdicts[0].verdict == "heterozygous"
    assert verdicts[0].observed == "Y" or verdicts[0].observed == "R" or len(verdicts[0].observed) == 2 or verdicts[0].observed in "RYSWMK"


def test_expected_het_iupac_alt_confirmed(setup) -> None:
    reference = setup
    alt = "T" if reference[50] != "T" else "A"
    read = confirmed_case(reference, 50, alt)
    chrom = make_chrom(read, heights={(50, "C"): 600})
    # expecting the mixed genotype itself (e.g. Y = C/T)
    expected_code = "Y" if {alt, "C"} == {"C", "T"} else None
    if expected_code is None:  # adapt to whatever pair alt+C forms
        from chrom_core.snps import alleles_to_iupac

        expected_code = alleles_to_iupac({alt, "C"})
    verdicts, _ = _verify(
        reference, read, [ExpectedVariant(51, reference[50], expected_code)], chrom=chrom
    )
    assert verdicts[0].verdict == "heterozygous"


def test_not_found_when_read_shows_reference(setup) -> None:
    reference = setup
    read = reference  # no mutation
    verdicts, _ = _verify(reference, read, [ExpectedVariant(51, reference[50], "T" if reference[50] != "T" else "A")])
    assert verdicts[0].verdict == "not_found"


def test_no_data_outside_coverage(setup) -> None:
    reference = setup
    read = reference[10:60]
    verdicts, _ = _verify(reference, read, [ExpectedVariant(5, reference[4], "A")])
    assert verdicts[0].verdict == "no_data"
    assert "coverage" in verdicts[0].note


def test_no_data_when_deleted(setup) -> None:
    reference = setup
    read = reference[:50] + reference[51:]  # deletion at 50
    verdicts, _ = _verify(reference, read, [ExpectedVariant(51, reference[50], "A")])
    assert verdicts[0].verdict == "no_data"
    assert "deleted" in verdicts[0].note


def test_ref_mismatch(setup) -> None:
    reference = setup
    wrong_ref = "A" if reference[50] != "A" else "C"
    verdicts, _ = _verify(reference, reference, [ExpectedVariant(51, wrong_ref, "T")])
    assert verdicts[0].verdict == "ref_mismatch"
    assert verdicts[0].observed == reference[50]


def test_ref_mismatch_carries_note_code(setup) -> None:
    reference = setup
    wrong_ref = "A" if reference[50] != "A" else "C"
    verdicts, _ = _verify(reference, reference, [ExpectedVariant(51, wrong_ref, "T")])
    assert verdicts[0].note_code == "ref_differs"
    assert verdicts[0].note  # English text stays for CSV/XLSX


def test_other_note_codes(setup) -> None:
    reference = setup
    verdicts, _ = _verify(reference, reference, [ExpectedVariant(51, "AT", "GC")])
    assert verdicts[0].note_code == "only_snv"
    verdicts, _ = _verify(reference, reference, [ExpectedVariant(5000, ".", "A")])
    assert verdicts[0].note_code == "pos_outside_ref"
    read = reference[:40] + reference[60:]
    verdicts, _ = _verify(reference, read, [ExpectedVariant(51, reference[50], "A")])
    assert verdicts[0].note_code == "deleted_in_read"
    alt = "T" if reference[50] != "T" else "A"
    chrom = make_chrom(confirmed_case(reference, 50, alt), quality=10)
    verdicts, _ = _verify(reference, confirmed_case(reference, 50, alt),
                          [ExpectedVariant(51, reference[50], alt)], chrom=chrom)
    assert verdicts[0].note_code == "low_q"
    assert verdicts[0].note_args == {"q": 10, "min_q": 20}


_COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


@pytest.mark.parametrize("name", ["Achl_ACHLO006-09_1_F.ab1", "Achl_ACHLO006-09_2_R.ab1"])
def test_real_read_confirmed_on_both_strands(name) -> None:
    from clone_fixtures import real_read, reference_for

    chrom = real_read(name)
    reference = reference_for(name)
    alignment = align_read("".join(chrom.seq), reference)
    assert alignment.strand == ("-" if name.endswith("_R.ab1") else "+")
    n = len(chrom.seq)
    checked = 0
    for pos, oriented_idx in sorted(alignment.ref_to_read.items()):
        if pos < 100:
            continue
        idx = oriented_idx if alignment.strand == "+" else n - 1 - oriented_idx
        if chrom.quality[idx] < 30 or chrom.seq[idx].upper() not in "ACGT":
            continue
        base = chrom.seq[idx].upper()
        ref_strand = base if alignment.strand == "+" else _COMP[base]
        verdict = verify_expected(
            alignment, chrom, reference, [ExpectedVariant(pos + 1, ".", ref_strand)]
        )[0]
        assert verdict.verdict == "confirmed", (pos, verdict)
        assert verdict.observed == ref_strand
        assert verdict.read_pos == idx
        checked += 1
        if checked >= 40:
            break
    assert checked == 40


def test_low_quality(setup) -> None:
    reference = setup
    read = confirmed_case(reference, 50, "T" if reference[50] != "T" else "A")
    chrom = make_chrom(read, quality=10)
    verdicts, _ = _verify(reference, read, [ExpectedVariant(51, reference[50], read[50])], chrom=chrom)
    assert verdicts[0].verdict == "low_quality"


def test_unsupported_non_snv(setup) -> None:
    reference = setup
    verdicts, _ = _verify(reference, reference, [ExpectedVariant(51, "AT", "GC")])
    assert verdicts[0].verdict == "unsupported"


# -------------------------------------------------------------------- reports
def test_verdicts_csv(setup) -> None:
    reference = setup
    read = confirmed_case(reference, 50, "T" if reference[50] != "T" else "A")
    verdicts, _ = _verify(reference, read, [ExpectedVariant(51, reference[50], read[50], "snp1")])
    text = verdicts_csv(verdicts)
    lines = text.strip().splitlines()
    assert lines[0].startswith("name;position;ref;alt;verdict")
    assert "snp1" in lines[1] and "confirmed" in lines[1]


def test_verdicts_xlsx(setup, tmp_path) -> None:
    openpyxl = pytest.importorskip("openpyxl")
    reference = setup
    read = confirmed_case(reference, 50, "T" if reference[50] != "T" else "A")
    verdicts, _ = _verify(reference, read, [ExpectedVariant(51, reference[50], read[50], "snp1")])
    path = tmp_path / "report.xlsx"
    from chrom_core.verify import verdicts_xlsx

    verdicts_xlsx(path, verdicts, sample_id="sample-1")
    book = openpyxl.load_workbook(path)
    sheet = book.active
    assert sheet.cell(row=1, column=1).value == "sample"
    assert sheet.cell(row=2, column=2).value == "snp1"
    assert sheet.cell(row=2, column=6).value == "confirmed"


# ---------------------------------------------------------------- real reads
@pytest.mark.skipif(not (EXAMPLES / "3730.ab1").exists(), reason="no example")
def test_real_read_end_to_end() -> None:
    chrom = read_ab1(EXAMPLES / "3730.ab1")
    reference = list(chrom.seq)
    alt = "A" if reference[500] != "A" else "G"
    reference[500] = alt  # reference carries the alt; read shows the ref base
    reference = "".join(reference)
    alignment = align_read(chrom.seq, reference)
    variants, warnings = parse_expected("501\t.\tA\n")  # alt=alt base
    # the read carries the ORIGINAL base, so the alt is not found
    variants = [ExpectedVariant(501, reference[500], alt)]
    verdicts = verify_expected(alignment, chrom, reference, variants)
    assert verdicts[0].verdict in ("not_found", "low_quality", "heterozygous")
    # now flip the expectation to the base the read actually carries
    actual = chrom.seq[500]
    verdicts = verify_expected(
        alignment, chrom, reference, [ExpectedVariant(501, reference[500], actual)]
    )
    assert verdicts[0].verdict in ("confirmed", "low_quality")
