import io

import pytest

from chrom_core import blast


def test_build_request_fields() -> None:
    params = blast.build_request("ACGTACGT", program="blastn", database="nt")
    assert params["CMD"] == "Put"
    assert params["PROGRAM"] == "blastn"
    assert params["DATABASE"] == "nt"
    assert params["QUERY"] == "ACGTACGT"
    assert params["MEGABLAST"] == "on"


def test_build_request_without_megablast() -> None:
    params = blast.build_request("ACGT", megablast=False)
    assert "MEGABLAST" not in params
    assert params["FILTER"] == "F"


def test_build_request_empty_sequence_raises() -> None:
    with pytest.raises(ValueError):
        blast.build_request("   ")


def test_submit_parses_rid(monkeypatch) -> None:
    class FakeResponse(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    sent = {}

    def fake_urlopen(request, timeout=None):
        sent["data"] = request.data.decode()
        sent["timeout"] = timeout
        return FakeResponse(b"<html>Query submitted. RID = 4VBXA123001 </html>")

    monkeypatch.setattr(blast.urllib.request, "urlopen", fake_urlopen)
    rid = blast.submit(blast.build_request("ACGT" * 20))
    assert rid == "4VBXA123001"
    assert "CMD=Put" in sent["data"] and "QUERY=ACGT" in sent["data"]


def test_submit_without_rid_raises_with_server_message(monkeypatch) -> None:
    class FakeResponse(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    def fake_urlopen(request, timeout=None):
        return FakeResponse(
            b"<html>Message from server<br> Invalid sequence </html>"
        )

    monkeypatch.setattr(blast.urllib.request, "urlopen", fake_urlopen)
    with pytest.raises(RuntimeError, match="Invalid sequence"):
        blast.submit(blast.build_request("ACGT"))


def test_result_url() -> None:
    assert blast.result_url("AB12") == f"{blast.BLAST_URL}?CMD=Get&RID=AB12"
