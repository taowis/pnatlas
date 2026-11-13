"""Create or update a Zenodo draft record for the data archive and upload
the bundled files into it. The record is left as a draft: publishing, which
registers the DOI for good, is done in the browser after a look at the draft.

Usage: ZENODO_TOKEN=... python workflow/release/zenodo_upload.py --dir <bundles>
       [--metadata workflow/release/zenodo_metadata.json] [--record-id N]
       [--sandbox] [--publication-date YYYY-MM-DD]
The token is read from ZENODO_TOKEN or from ~/.config/zenodo/token and needs
the deposit:write scope. Files already in the draft with the same MD5 are
skipped, so the command can be rerun after an interrupted upload.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from datetime import date
from pathlib import Path

import requests


def token() -> str:
    t = os.environ.get("ZENODO_TOKEN") or ""
    p = Path.home() / ".config" / "zenodo" / "token"
    if not t and p.exists():
        t = p.read_text().strip()
    if not t:
        sys.exit("no token: set ZENODO_TOKEN or write ~/.config/zenodo/token")
    return t


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, type=Path)
    ap.add_argument("--metadata", default=Path("workflow/release/zenodo_metadata.json"), type=Path)
    ap.add_argument("--record-id", type=int)
    ap.add_argument("--sandbox", action="store_true")
    ap.add_argument("--publication-date", default=date.today().isoformat())
    a = ap.parse_args()
    base = "https://sandbox.zenodo.org" if a.sandbox else "https://zenodo.org"
    s = requests.Session()
    s.headers["Authorization"] = f"Bearer {token()}"

    if a.record_id:
        r = s.get(f"{base}/api/deposit/depositions/{a.record_id}", timeout=60)
    else:
        r = s.post(f"{base}/api/deposit/depositions", json={}, timeout=60)
    r.raise_for_status()
    dep = r.json()
    dep_id, bucket = dep["id"], dep["links"]["bucket"]

    meta = json.loads(a.metadata.read_text())
    meta["publication_date"] = a.publication_date
    r = s.put(f"{base}/api/deposit/depositions/{dep_id}", json={"metadata": meta}, timeout=60)
    if r.status_code >= 400:
        sys.exit(f"metadata rejected: {r.status_code} {r.text}")
    dep = r.json()
    doi = dep["metadata"].get("prereserve_doi", {}).get("doi", "")
    print(f"record {dep_id} draft at {dep['links']['html']}")
    print(f"reserved DOI {doi}")

    have = {f["filename"]: f["checksum"] for f in s.get(f"{base}/api/deposit/depositions/{dep_id}/files", timeout=60).json()}
    for path in sorted(p for p in a.dir.iterdir() if p.is_file()):
        local = md5(path)
        if have.get(path.name) == local:
            print(f"present  {path.name}")
            continue
        with path.open("rb") as fh:
            r = s.put(f"{bucket}/{path.name}", data=fh, timeout=7200)
        if r.status_code >= 400:
            sys.exit(f"upload of {path.name} failed: {r.status_code} {r.text}")
        remote = r.json()["checksum"].removeprefix("md5:")
        state = "uploaded" if remote == local else "MISMATCH"
        print(f"{state} {path.name} {path.stat().st_size} bytes")
        if state == "MISMATCH":
            sys.exit(1)
    print(f"DRAFT READY {dep_id} {doi}")


if __name__ == "__main__":
    main()
