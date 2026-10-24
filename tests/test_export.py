from chrom_core.export import to_fasta, to_fastq, write_text


def test_fasta_header_and_wrap() -> None:
    text = to_fasta("sample1", "A" * 75)
    lines = text.splitlines()
    assert lines[0] == ">sample1"
    assert [len(line) for line in lines[1:]] == [60, 15]
    assert text.endswith("\n")


def test_fastq_phred33() -> None:
    text = to_fastq("s", "ACGT", (10, 20, 30, 40))
    lines = text.splitlines()
    assert lines[0] == "@s"
    assert lines[1] == "ACGT"
    assert lines[2] == "+"
    assert lines[3] == chr(43) + chr(53) + chr(63) + chr(73)


def test_fastq_clamps_high_quality() -> None:
    text = to_fastq("s", "A", (99,))
    assert text.splitlines()[3] == chr(126)  # 33 + 93


def test_fastq_length_mismatch_raises() -> None:
    import pytest

    with pytest.raises(ValueError):
        to_fastq("s", "ACG", (1, 2))


def test_write_text_lf_only(tmp_path) -> None:
    path = tmp_path / "out.fasta"
    write_text(path, "X\nY\n")
    raw = path.read_bytes()
    assert raw == b"X\nY\n"
    assert b"\r" not in raw
