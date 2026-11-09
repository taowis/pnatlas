"""Verify files downloaded by hand from a Dryad dataset and record them in the manifest.

Dryad serves file downloads only to a browser session or with an API token
(the API download link answers 401 without one and the web download link 403
to a script), so the files of a Dryad dataset are saved from the dataset page
in a browser and moved into <data_root>/raw/<dataset>/<subdir>/. This script
reads the file list of the dataset version from the Dryad API (name, size,
sha-256), checks every downloaded file against it, and appends one manifest
row per verified file (url = the web download link of the file, md5 computed
here, sha-256 in the notes). Files not yet present are reported and skipped,
so the script can be rerun as the downloads come in.

Usage:
  python workflow/scripts/dryad_verify.py --list
  python workflow/scripts/dryad_verify.py --data-root data --dir raw/yu2024_hdrg_xenium/dryad \
      --dataset yu2024_hdrg_xenium --manifest data/MANIFEST.tsv
"""
import argparse
import csv
import datetime as dt
import hashlib
import json
import pathlib
import sys
import urllib.request

DOI = "10.5061/dryad.gf1vhhmxq"
VERSION_ID = 324492  # version 3 of the dataset
API = "https://datadryad.org/api/v2"
WEB = "https://datadryad.org/downloads/file_stream"
MANIFEST_COLS = ["path", "dataset", "category", "source", "accession", "url", "version",
                 "access_date", "md5", "size_bytes", "fetch", "access", "notes"]


def api(path):
    with urllib.request.urlopen(f"{API}{path}", timeout=60) as r:
        return json.load(r)


def dryad_files(version_id):
    out = []
    page = f"/versions/{version_id}/files"
    while page:
        d = api(page)
        for f in d["_embedded"]["stash:files"]:
            fid = f["_links"]["stash:download"]["href"].rstrip("/").split("/")[-2]
            out.append({"name": f["path"], "size": int(f["size"]), "digest": f.get("digest", ""),
                        "digest_type": f.get("digestType", ""), "id": fid})
        page = d["_links"].get("next", {}).get("href")
    return out


def digests(path):
    md5, sha = hashlib.md5(), hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            md5.update(chunk)
            sha.update(chunk)
    return md5.hexdigest(), sha.hexdigest()


def read_manifest(path):
    if not path.exists():
        return {}
    with open(path, newline="") as f:
        return {r["path"]: r for r in csv.DictReader(f, delimiter="\t")}


def write_manifest(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=MANIFEST_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for k in sorted(rows):
            w.writerow(rows[k])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true", help="print the dataset's file list and exit")
    ap.add_argument("--data-root")
    ap.add_argument("--dir", help="directory under data root holding the downloaded files")
    ap.add_argument("--dataset")
    ap.add_argument("--manifest", default="data/MANIFEST.tsv")
    ap.add_argument("--version-id", type=int, default=VERSION_ID)
    ap.add_argument("--doi", default=DOI)
    a = ap.parse_args()
    files = dryad_files(a.version_id)
    version = api(f"/datasets/doi%3A{a.doi.replace('/', '%2F')}").get("versionNumber", "")
    if a.list:
        print("name\tsize_bytes\tdigest_type\tdigest\tfile_id")
        for f in files:
            print(f"{f['name']}\t{f['size']}\t{f['digest_type']}\t{f['digest']}\t{f['id']}")
        return 0
    if not (a.data_root and a.dir and a.dataset):
        ap.error("--data-root, --dir and --dataset are required unless --list")
    root, manifest_path = pathlib.Path(a.data_root), pathlib.Path(a.manifest)
    manifest = read_manifest(manifest_path)
    bad = 0
    for f in files:
        rel = f"{a.dir}/{f['name']}"
        p = root / rel
        if not p.exists():
            print(f"missing\t{rel}", flush=True)
            continue
        size = p.stat().st_size
        if size != f["size"]:
            print(f"SIZE_MISMATCH\t{rel}\t{size}\t{f['size']}", flush=True)
            bad += 1
            continue
        md5, sha = digests(p)
        if f["digest_type"] == "sha-256" and sha != f["digest"]:
            print(f"SHA_MISMATCH\t{rel}\t{sha}\t{f['digest']}", flush=True)
            bad += 1
            continue
        manifest[rel] = {
            "path": rel, "dataset": a.dataset, "category": "raw", "source": "Dryad",
            "accession": f"doi:{a.doi}", "url": f"{WEB}/{f['id']}", "version": f"version {version}",
            "access_date": manifest.get(rel, {}).get("access_date") or dt.date.today().isoformat(),
            "md5": md5, "size_bytes": str(size),
            "fetch": "manual: saved from the dataset page in a browser; workflow/scripts/dryad_verify.py",
            "access": "open", "notes": f"sha-256 {sha} matches Dryad",
        }
        print(f"ok\t{rel}\t{md5}\t{size}", flush=True)
    write_manifest(manifest_path, manifest)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
