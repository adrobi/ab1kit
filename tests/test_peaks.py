import random
from dataclasses import replace

from chrom_core.align import align_read
from chrom_core.peaks import (
    SecondaryPeak,
    channel_baseline,
    peak_window,
    scan_secondary_peaks,
    secondary_peak,
)
from chrom_core.snps import peak_heights
from clone_fixtures import REAL_READS, make_chrom, real_read, reference_for


def _random_seq(seed: int, n: int = 120) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def test_weak_channel_normalisation_removes_false_het() -> None:
    # The G dye is systematically half as bright: every called G peak is 500
    # high, every other called peak 1000. At position 10 (called G) the A
    # channel shows 225 — 45% of the raw G height, a false het by raw ratio.
    # Normalised, the G peak itself is at its baseline (500/500 = 1.0) while
    # A is 225/1000 = 0.225 of its own baseline -> ratio well below 0.40.
    seq = list(_random_seq(1, 40))
    seq[10] = "G"
    seq = "".join(seq)
    chrom = make_chrom(seq, primary={"G": 500}, heights={(10, "A"): 225})
    raw = peak_heights(chrom, 10)
    assert raw["A"] / raw["G"] >= 0.40  # the raw ratio would call a het
    peak = secondary_peak(chrom, 10)
    assert peak.primary_base == "G"
    assert peak.base == "A"
    assert peak.ratio < 0.40


def test_true_het_detected_after_normalisation() -> None:
    seq = list(_random_seq(2, 40))
    seq[10] = "A"
    seq = "".join(seq)
    chrom = make_chrom(seq, heights={(10, "C"): 500})  # second C peak at 50%
    peak = secondary_peak(chrom, 10)
    assert peak.base == "C"
    assert 0.45 <= peak.ratio <= 0.55


def test_peak_window_tracks_spacing() -> None:
    chrom_spacing_10 = make_chrom("ACGTACGTAC", spacing=10)
    chrom_spacing_3 = make_chrom("ACGTACGTAC", spacing=3)
    assert peak_window(chrom_spacing_10, 5) == 5
    assert peak_window(chrom_spacing_3, 5) == 2


def test_neighbour_shoulder_is_not_a_second_peak() -> None:
    # Position 11 is a broad A peak (spacing 10): its rising flank reaches
    # 600 at the edge of position 10's window. A plain window maximum would
    # read that as a 60% second peak; it is a slope, not a peak.
    seq = "CGTACGTACGCAGTCGATCGAT"
    chrom = make_chrom(seq)
    assert seq[10] == "C" and seq[11] == "A"
    a = list(chrom.traces["A"])
    for k, v in zip(range(111, 120), (300, 350, 400, 500, 600, 700, 800, 900, 950)):
        a[k] = v
    a[105:111] = [100, 120, 150, 180, 220, 260]
    chrom = replace(chrom, traces={**chrom.traces, "A": tuple(a)})
    assert peak_heights(chrom, 10)["A"] >= 600  # raw window max sees the shoulder
    peak = secondary_peak(chrom, 10)
    assert peak.ratio < 0.40


def test_channel_baseline_uses_called_peaks_per_channel() -> None:
    seq = "ACGT" * 10
    chrom = make_chrom(seq, primary={"G": 500})
    baseline = channel_baseline(chrom, 20)
    assert baseline == {"A": 1000.0, "C": 1000.0, "G": 500.0, "T": 1000.0}


def test_channel_baseline_falls_back_when_channel_never_called() -> None:
    # no T is called anywhere: T falls back to the median of all called peaks
    chrom = make_chrom("ACG" * 15, primary={"G": 400})
    baseline = channel_baseline(chrom, 20)
    assert baseline["G"] == 400.0
    assert baseline["T"] == 1000.0  # median of 1000/1000/400 heights


def test_secondary_peak_without_outside_channel_is_unmeasurable() -> None:
    chrom = make_chrom("ACGTNACGTA")
    peak = secondary_peak(chrom, 4)
    assert peak == SecondaryPeak(primary_base=peak.primary_base, base="", ratio=None)


def test_scan_finds_het_where_reference_base_was_called() -> None:
    ref = list(_random_seq(3))
    if ref[30] == "T":
        ref[30] = "A"
    ref = "".join(ref)
    # the read equals the reference, but under position 30 a T peak at 50%
    chrom = make_chrom(ref, heights={(30, "T"): 500})
    hits = scan_secondary_peaks(chrom, align_read(chrom.seq, ref), 0.40, 20)
    assert [pos for pos, _ in hits] == [30]
    assert hits[0][1].base == "T"


def test_scan_maps_reverse_strand_to_displayed_read() -> None:
    from chrom_core.seqops import reverse_complement

    forward = list(_random_seq(4))
    if forward[30] == "T":
        forward[30] = "A"
    forward = "".join(forward)
    chrom = make_chrom(forward, heights={(30, "T"): 500})
    ref = reverse_complement(forward)
    aln = align_read(chrom.seq, ref)
    assert aln.strand == "-"
    hits = scan_secondary_peaks(chrom, aln, 0.40, 20)
    assert [pos for pos, _ in hits] == [len(forward) - 1 - 30]


def test_scan_skips_low_quality() -> None:
    ref = list(_random_seq(5))
    if ref[30] == "T":
        ref[30] = "A"
    ref = "".join(ref)
    chrom = make_chrom(ref, quality=10, heights={(30, "T"): 500})
    assert scan_secondary_peaks(chrom, align_read(chrom.seq, ref), 0.40, 20) == []


def test_real_reads_have_no_mass_false_hets() -> None:
    assert REAL_READS
    for path in REAL_READS:
        chrom = real_read(path.name)
        aln = align_read(chrom.seq, reference_for(path.name))
        hits = scan_secondary_peaks(chrom, aln, 0.40, 20)
        assert len(hits) <= 2, (path.name, hits)
