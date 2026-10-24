"""One-time download of the real-data test corpus into ``tests/data``.

Usage::

    python scripts/fetch_test_data.py --coi ACCESSION [--force]

Uses only ``urllib``. Idempotent: files that already exist are skipped
unless ``--force`` is given. Tests never touch the network; the downloaded
files are committed to the repository.
"""

from __future__ import annotations

import argparse
import time
import urllib.parse
import urllib.request
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "tests" / "data"

SANGER_COMMIT = "4c5ebc9298a873c209819f12b7dfd85e126f1406"
BIOPYTHON_COMMIT = "2fbe8755bb44731e60b7bd1e2cc55cf34d94e15b"

SANGER_BASE = (
    "https://raw.githubusercontent.com/roblanf/sangeranalyseR/"
    f"{SANGER_COMMIT}/inst/extdata/Allolobophora_chlorotica/ACHLO/"
)
SNAPGENE_URL = (
    "https://raw.githubusercontent.com/biopython/biopython/"
    f"{BIOPYTHON_COMMIT}/Tests/SnapGene/sample-d.dna"
)
EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

READS = [
    f"Achl_{sample}_{suffix}.ab1"
    for sample in ("ACHLO006-09", "ACHLO007-09", "ACHLO040-09", "ACHLO041-09")
    for suffix in ("1_F", "2_R")
]

_last_ncbi_request = 0.0


def _download(url: str) -> bytes:
    request = urllib.request.Request(url, headers={"User-Agent": "ab1kit-fetch-test-data"})
    with urllib.request.urlopen(request, timeout=60) as response:
        return response.read()


def _genbank_url(accession: str) -> str:
    query = urllib.parse.urlencode(
        {"db": "nuccore", "id": accession, "rettype": "gbwithparts", "retmode": "text"}
    )
    return f"{EFETCH}?{query}"


def _fetch(url: str, target: Path, force: bool, ncbi: bool = False) -> None:
    global _last_ncbi_request
    if target.exists() and not force:
        print(f"skip   {target.relative_to(DATA)}")
        return
    if ncbi:  # stay polite: at most ~3 requests per second
        wait = 0.4 - (time.monotonic() - _last_ncbi_request)
        if wait > 0:
            time.sleep(wait)
        _last_ncbi_request = time.monotonic()
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(_download(url))
    print(f"fetched {target.relative_to(DATA)}")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coi", metavar="ACC", help="GenBank accession of the COI reference")
    parser.add_argument("--force", action="store_true", help="re-download existing files")
    args = parser.parse_args(argv)

    for name in READS:
        _fetch(SANGER_BASE + name, DATA / "real" / name, args.force)
    _fetch(_genbank_url("L09137"), DATA / "refs" / "pUC19_L09137.gb", args.force, ncbi=True)
    if args.coi:
        _fetch(_genbank_url(args.coi), DATA / "refs" / f"coi_{args.coi}.gb", args.force, ncbi=True)
    _fetch(SNAPGENE_URL, DATA / "refs" / "sample-d.dna", args.force)


if __name__ == "__main__":
    main()
