from clone_fixtures import REAL_READS, COI_GB, PUC19_GB, SNAPGENE_DNA, real_read, coi_reference
from chrom_core.align import align_read

# ACHLO040-09 and ACHLO041-09 belong to other, deeply divergent lineages of the
# Allolobophora chlorotica aggregate than the single GenBank COI record we ship
# (HM417896, matches ACHLO006-09/007-09). See tests/data/README.md.
DIVERGENT_SAMPLES = ("ACHLO040-09", "ACHLO041-09")


def _expected_strand(name: str) -> str:
    return "+" if "_F" in name else "-"


def test_corpus_files_present():
    assert len(REAL_READS) == 8
    assert COI_GB.exists() and PUC19_GB.exists() and SNAPGENE_DNA.exists()


def test_every_real_read_aligns_to_coi_with_high_identity():
    ref = coi_reference()
    for path in REAL_READS:
        if any(s in path.name for s in DIVERGENT_SAMPLES):
            continue
        aln = align_read(real_read(path.name).seq, ref)
        assert aln.identity >= 0.97, path.name
        assert aln.strand == _expected_strand(path.name)


def test_divergent_lineage_reads_still_align_on_the_right_strand():
    ref = coi_reference()
    names = [p.name for p in REAL_READS if any(s in p.name for s in DIVERGENT_SAMPLES)]
    assert len(names) == 4
    for name in names:
        aln = align_read(real_read(name).seq, ref)
        assert aln.strand == _expected_strand(name), name
        assert 0.75 <= aln.identity < 0.97, name
