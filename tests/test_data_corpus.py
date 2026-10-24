from clone_fixtures import (
    REAL_READS, COI_GB, PUC19_GB, SNAPGENE_DNA, real_read, coi_reference, reference_for,
)
from chrom_core.align import align_read


def _expected_strand(name: str) -> str:
    return "+" if "_F" in name else "-"


def test_corpus_files_present():
    assert len(REAL_READS) == 8
    assert COI_GB.exists() and PUC19_GB.exists() and SNAPGENE_DNA.exists()


def test_every_real_read_aligns_to_coi_with_high_identity():
    for path in REAL_READS:
        aln = align_read(real_read(path.name).seq, reference_for(path.name))
        assert aln.identity >= 0.97, path.name
        assert aln.strand == _expected_strand(path.name), path.name


def test_achlo040_is_a_different_lineage_from_the_primary_reference():
    aln = align_read(real_read("Achl_ACHLO040-09_1_F.ab1").seq, coi_reference())
    assert aln.identity < 0.90
