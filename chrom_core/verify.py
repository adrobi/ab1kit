"""Expected-SNP verification: confirm designed genotypes from a variant list.

Input: a VCF or CSV/TSV table of expected substitutions (position in the
reference, ref/alt base, optional name) plus a read alignment. Every
expected variant gets a verdict:

- ``confirmed``     — read carries the alt allele (clean single peak);
- ``heterozygous``  — alt present together with a second allele (double
                      peak or mixed basecaller call);
- ``not_found``     — read shows something else (reference or third base);
- ``no_data``       — position outside read coverage or deleted in the read;
- ``low_quality``   — phred below the threshold, verdict untrusted;
- ``ref_mismatch``  — the ref base in the table disagrees with the actual
                      reference (input error guard);
- ``unsupported``   — not a single-nucleotide variant.
"""

from __future__ import annotations

from dataclasses import dataclass

from chrom_core.abif import Chromatogram
from chrom_core.align import ReadAlignment
from chrom_core.seqops import IUPAC
from chrom_core.peaks import secondary_peak
from chrom_core.snps import DEFAULT_HET_RATIO, DEFAULT_MIN_Q, alleles_to_iupac


@dataclass(frozen=True)
class ExpectedVariant:
    position: int      # 1-based position in the reference
    ref: str           # expected reference base ('.' / '' / '-' when unknown)
    alt: str
    name: str = ""


@dataclass(frozen=True)
class VariantVerdict:
    variant: ExpectedVariant
    verdict: str
    observed: str      # what the read shows (base or IUPAC code), '' when no data
    quality: int | None = None
    ratio: float | None = None
    read_pos: int | None = None  # 0-based index in the displayed read
    note: str = ""


# ------------------------------------------------------------------- parsing
def parse_expected(text: str) -> tuple[list[ExpectedVariant], list[str]]:
    """Parse a VCF or CSV/TSV table; bad lines are skipped with warnings."""
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if not lines:
        raise ValueError("empty variant list")

    variants: list[ExpectedVariant] = []
    warnings: list[str] = []

    if any(line.startswith("#") for line in lines):  # VCF
        for line in lines:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 5:
                warnings.append(f"VCF row skipped (too few columns): {line[:50]}")
                continue
            try:
                position = int(fields[1])
            except ValueError:
                warnings.append(f"VCF row skipped (bad POS): {fields[1][:20]}")
                continue
            name = fields[2] if fields[2] not in (".", "") else ""
            variants.append(ExpectedVariant(position, fields[3].upper(), fields[4].split(",")[0].upper(), name))
        return variants, warnings

    # CSV/TSV: sniff the delimiter from the first line
    sample = lines[0]
    delimiter = "\t" if "\t" in sample else ";" if ";" in sample else ","
    columns = {"pos": 0, "ref": 1, "alt": 2, "name": 3}
    header_used = False
    if not _looks_like_number(sample.split(delimiter)[0]) and any(
        keyword in sample.lower()
        for keyword in ("pos", "position", "позиция", "ref", "alt", "chr")
    ):
        header_used = True
        header = [h.strip().lower() for h in sample.split(delimiter)]
        for key, default in columns.items():
            for index, hname in enumerate(header):
                if key in hname:
                    columns[key] = index
                    break
        lines = lines[1:]

    for line in lines:
        fields = [f.strip() for f in line.split(delimiter)]
        try:
            position = int(fields[columns["pos"]])
        except (ValueError, IndexError):
            warnings.append(f"row skipped (bad position): {line[:50]}")
            continue
        if not header_used and len(fields) == 2:
            # two-column shorthand: (position, alt), ref unknown
            variants.append(ExpectedVariant(position, ".", fields[1].upper(), ""))
            continue
        if len(fields) <= columns["alt"]:
            warnings.append(f"row skipped (no alt column): {line[:50]}")
            continue
        ref = fields[columns["ref"]] if len(fields) > columns["ref"] else "."
        name = fields[columns["name"]] if len(fields) > columns["name"] else ""
        variants.append(
            ExpectedVariant(position, (ref or ".").upper(), fields[columns["alt"]].upper(), name)
        )
    return variants, warnings


def _looks_like_number(text: str) -> bool:
    try:
        int(text)
        return True
    except ValueError:
        return False


# --------------------------------------------------------------- verification
def verify_expected(
    alignment: ReadAlignment,
    chrom: Chromatogram,
    reference: str,
    variants: list[ExpectedVariant],
    min_q: int = DEFAULT_MIN_Q,
    het_ratio: float = DEFAULT_HET_RATIO,
) -> list[VariantVerdict]:
    reference = reference.upper()
    return [
        _verify_one(v, alignment, chrom, reference, min_q, het_ratio)
        for v in variants
    ]


def _verify_one(
    variant: ExpectedVariant,
    alignment: ReadAlignment,
    chrom: Chromatogram,
    reference: str,
    min_q: int,
    het_ratio: float,
) -> VariantVerdict:
    pos = variant.position - 1
    is_snv = (
        len(variant.alt) == 1
        and variant.alt in IUPAC
        and (len(variant.ref) == 1 or variant.ref in (".", "", "-"))
    )
    if not is_snv:
        return VariantVerdict(variant, "unsupported", "", note="only SNV rows are supported")
    if not 0 <= pos < len(reference):
        return VariantVerdict(variant, "unsupported", "", note="position outside the reference")

    if variant.ref not in (".", "", "-") and variant.ref.upper() != reference[pos]:
        return VariantVerdict(
            variant, "ref_mismatch", reference[pos], note="table ref base differs from the reference"
        )

    if pos not in alignment.ref_to_read:
        note = (
            "deleted in the read"
            if alignment.ref_start <= pos < alignment.ref_end
            else "outside read coverage"
        )
        return VariantVerdict(variant, "no_data", "", note=note)

    read_idx = alignment.ref_to_read[pos]
    quality = chrom.quality[read_idx]
    called = chrom.seq[read_idx].upper()
    called_alleles = set(IUPAC.get(called, "ACGT"))
    peak = secondary_peak(chrom, read_idx)
    secondary_base, ratio = peak.base, peak.ratio

    alleles = set(called_alleles)
    if ratio is not None and ratio >= het_ratio and secondary_base:
        alleles.add(secondary_base)
    observed = alleles_to_iupac(alleles)

    if quality < min_q:
        return VariantVerdict(
            variant, "low_quality", observed, quality, ratio, read_idx, f"Q{quality} < {min_q}"
        )

    alt_alleles = set(IUPAC.get(variant.alt, variant.alt))
    if alt_alleles <= alleles:
        verdict = "heterozygous" if len(alleles) > 1 else "confirmed"
    else:
        verdict = "not_found"
    return VariantVerdict(variant, verdict, observed, quality, ratio, read_idx)


# -------------------------------------------------------------------- reports
_CSV_HEADER = ("name", "position", "ref", "alt", "verdict", "observed", "quality", "ratio", "note")


def verdicts_csv(verdicts: list[VariantVerdict]) -> str:
    """Semicolon-separated report (RU-Excel friendly, like batch summary)."""
    lines = [";".join(_CSV_HEADER)]
    for v in verdicts:
        ratio = f"{v.ratio:.2f}" if v.ratio is not None else ""
        row = (
            v.variant.name,
            v.variant.position,
            v.variant.ref,
            v.variant.alt,
            v.verdict,
            v.observed,
            v.quality if v.quality is not None else "",
            ratio,
            v.note,
        )
        lines.append(";".join(str(f).replace(";", ",") for f in row))
    return "\n".join(lines) + "\n"


def verdicts_xlsx(path, verdicts: list[VariantVerdict], sample_id: str = "") -> None:
    """Excel report via openpyxl (imported lazily; optional dependency)."""
    try:
        from openpyxl import Workbook
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise RuntimeError("openpyxl is required for Excel reports: pip install openpyxl") from exc

    book = Workbook()
    sheet = book.active
    sheet.title = "expected SNPs"
    header = ("sample",) + _CSV_HEADER
    sheet.append(header)
    for cell in sheet[1]:
        cell.font = cell.font.copy(bold=True)
    for v in verdicts:
        ratio = round(v.ratio, 2) if v.ratio is not None else None
        sheet.append(
            (
                sample_id,
                v.variant.name,
                v.variant.position,
                v.variant.ref,
                v.variant.alt,
                v.verdict,
                v.observed,
                v.quality,
                ratio,
                v.note,
            )
        )
    widths = {"A": 28, "B": 14, "C": 10, "D": 8, "E": 8, "F": 14, "G": 10, "H": 9, "I": 8, "J": 30}
    for column, width in widths.items():
        sheet.column_dimensions[column].width = width
    book.save(path)
