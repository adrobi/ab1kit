import pytest
from Bio import SeqIO

from clone_fixtures import COI_GB, PUC19_GB, SNAPGENE_DNA
from chrom_core.reference import Construct, Feature, construct_from_text, load_construct


def test_genbank_puc19():
    c = load_construct(PUC19_GB)
    assert len(c.seq) == len(SeqIO.read(PUC19_GB, "genbank").seq) and c.circular
    # The committed L09137 record annotates only "source" (no CDS); CDS is covered by COI below.
    assert [f.kind for f in c.features] == ["source"]


def test_coi_has_cds():
    assert any(f.kind == "CDS" for f in load_construct(COI_GB).features)


def test_puc19_length_is_2686():
    assert len(load_construct(PUC19_GB).seq) == 2686


def test_snapgene_sample():
    c = load_construct(SNAPGENE_DNA)
    assert c.seq and c.features


def test_coi_cds_codon_start_and_name():
    c = load_construct(COI_GB)
    cds = [f for f in c.features if f.kind == "CDS"]
    assert cds and all(f.codon_start in (1, 2, 3) for f in cds)
    rec = SeqIO.read(COI_GB, "genbank")
    want = int(next(f for f in rec.features if f.type == "CDS").qualifiers.get("codon_start", ["1"])[0])
    assert cds[0].codon_start == want
    assert all(f.name for f in c.features)


def test_fasta_multi_record_lowercase(tmp_path):
    p = tmp_path / "x.fasta"
    p.write_text(">a\nacgtn\n>b\nTTTT\n")
    c = load_construct(p)
    assert c.seq == "ACGTN" and c.circular and c.features == ()


def test_unknown_extension_sniffed_fasta(tmp_path):
    p = tmp_path / "x.seq"
    p.write_text(">a\nACGT\n")
    assert load_construct(p).seq == "ACGT"


def test_unknown_extension_sniffed_genbank(tmp_path):
    p = tmp_path / "x.seq"
    p.write_bytes(PUC19_GB.read_bytes())
    assert len(load_construct(p).seq) == 2686


def test_garbage_rejected(tmp_path):
    p = tmp_path / "x.seq"
    p.write_text("hello world")
    with pytest.raises(ValueError):
        load_construct(p)


def test_invalid_chars_rejected(tmp_path):
    p = tmp_path / "x.fa"
    p.write_text(">a\nAC!GT\n")
    with pytest.raises(ValueError):
        load_construct(p)


def test_with_circular_toggle():
    assert not construct_from_text("ACGT").with_circular(False).circular


def test_construct_from_text_defaults():
    c = construct_from_text(">x\nacgt\n")
    assert c == Construct("pasted", "ACGT", True, ())


def test_origin_spanning_feature(tmp_path):
    p = tmp_path / "o.gb"
    p.write_text(
        "LOCUS       t                         12 bp    DNA     circular SYN 01-JAN-2000\n"
        "FEATURES             Location/Qualifiers\n"
        "     misc_feature    join(9..12,1..3)\n"
        '                     /label="wrap"\n'
        "ORIGIN\n        1 acgtacgtac gt\n//\n"
    )
    f = load_construct(p).features[0]
    assert (f.name, f.start, f.end, f.strand, f.codon_start) == ("wrap", 8, 3, 1, 1)
