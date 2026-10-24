"""Submitting reads to NCBI BLAST via the www-blast URL API.

Only the submission (``CMD=Put``) is done by ab1kit; the RID is then opened
in the user's browser, where NCBI itself shows progress and results. No
polling, no API key required.
"""

from __future__ import annotations

import re
import urllib.parse
import urllib.request

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
_RID_RE = re.compile(r"RID\s*=\s*(\S+)")
_ERROR_RE = re.compile(r"Message from server.*?<br>\s*(.*?)\s*<", re.DOTALL)


def build_request(
    seq: str,
    program: str = "blastn",
    database: str = "nt",
    megablast: bool = True,
    hitlist: int = 50,
) -> dict[str, str]:
    """POST payload for ``CMD=Put``."""
    if not seq.strip():
        raise ValueError("empty sequence")
    params = {
        "CMD": "Put",
        "PROGRAM": program,
        "DATABASE": database,
        "QUERY": seq,
        "HITLIST_SIZE": str(hitlist),
        "FILTER": "L" if megablast else "F",  # low-complexity filter, megablast default
    }
    if program == "blastn" and megablast:
        params["MEGABLAST"] = "on"
    return params


def submit(params: dict[str, str], timeout: float = 20.0) -> str:
    """Submit the request and return the RID assigned by NCBI."""
    data = urllib.parse.urlencode(params).encode("ascii")
    request = urllib.request.Request(
        BLAST_URL, data=data, headers={"User-Agent": "ab1kit/1.0"}
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        page = response.read().decode("utf-8", errors="replace")
    match = _RID_RE.search(page)
    if match:
        return match.group(1)
    error = _ERROR_RE.search(page)
    detail = error.group(1) if error else page[:200]
    raise RuntimeError(f"NCBI BLAST did not return a RID: {detail}")


def result_url(rid: str) -> str:
    """Status/results page for a submitted RID."""
    return f"{BLAST_URL}?CMD=Get&RID={urllib.parse.quote(rid)}"
