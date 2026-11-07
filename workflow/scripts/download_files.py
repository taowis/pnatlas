"""Download the files listed in a stage file list and record them in the manifest.

Reads a TSV with columns dataset, accession, file, url, size_bytes, access
(as written by stage2_list_files.py), fetches each file with curl into
<data_root>/raw/<dataset>/, computes its md5 and size, and appends one row per
file to data/MANIFEST.tsv. Files already present with the recorded md5 are
skipped, so the script can be rerun after an interruption.

Usage:
  python workflow/scripts/download_files.py data/metadata/stage2_files.tsv \
      --data-root "$DATA_ROOT" --manifest data/MANIFEST.tsv --fetch "workflow/scripts/download_files.py data/metadata/stage2_files.tsv"
  python workflow/scripts/download_files.py data/metadata/stage2_files.tsv --data-root "$DATA_ROOT" --verify
"""
import argparse
import csv
import datetime as dt
import hashlib
import pathlib
import subprocess
import sys

MANIFEST_COLS = ["path", "dataset", "category", "source", "accession", "url", "version",
                 "access_date", "md5", "size_bytes", "fetch", "access", "notes"]


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_manifest(path):
    if not path.exists():
        return {}
    with open(path, newline="") as f:
        return {r["path"]: r for r in csv.DictReader(f, delimiter="\t")}


def write_manifest(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=MANIFEST_COLS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in sorted(rows.values(), key=lambda r: r["path"]):
            w.writerow(r)


def source_of(accession):
    if accession.startswith("GSE"):
        return "GEO"
    if accession.startswith("SPARC"):
        return "SPARC"
    if accession.startswith("Box"):
        return "Box (UT Dallas)"
    return "painseq portal (Dropbox)"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files_tsv")
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--manifest", default="data/MANIFEST.tsv")
    ap.add_argument("--fetch", default="workflow/scripts/download_files.py")
    ap.add_argument("--verify", action="store_true", help="only check md5 of present files")
    ap.add_argument("--only", help="comma-separated dataset ids to process")
    a = ap.parse_args()

    root = pathlib.Path(a.data_root)
    manifest_path = pathlib.Path(a.manifest)
    manifest = read_manifest(manifest_path)
    only = set(a.only.split(",")) if a.only else None
    failed = 0

    with open(a.files_tsv, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    for r in rows:
        if only and r["dataset"] not in only:
            continue
        rel = f"raw/{r['dataset']}/{pathlib.Path(r['file']).name}"
        dest = root / rel
        known = manifest.get(rel)
        if a.verify:
            if not dest.exists():
                print(f"MISSING\t{rel}")
                failed += 1
            elif known and md5sum(dest) == known["md5"]:
                print(f"OK\t{rel}")
            else:
                print(f"MISMATCH\t{rel}")
                failed += 1
            continue
        if known and dest.exists() and dest.stat().st_size == int(known["size_bytes"]):
            print(f"present\t{rel}", flush=True)
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        print(f"fetch\t{rel}", flush=True)
        cmd = ["curl", "-L", "--fail", "--retry", "5", "--retry-delay", "30", "-C", "-",
               "-A", "pnatlas-stage2", "-o", str(dest), r["url"]]
        if subprocess.run(cmd, check=False).returncode != 0:
            print(f"FAILED\t{rel}", flush=True)
            failed += 1
            continue
        manifest[rel] = {
            "path": rel, "dataset": r["dataset"], "category": "raw",
            "source": source_of(r["accession"]), "accession": r["accession"], "url": r["url"],
            "version": r["accession"] if r["accession"].startswith("SPARC") else "",
            "access_date": dt.date.today().isoformat(), "md5": md5sum(dest),
            "size_bytes": dest.stat().st_size, "fetch": a.fetch, "access": r["access"],
            "notes": "",
        }
        write_manifest(manifest_path, manifest)
        print(f"done\t{rel}\t{manifest[rel]['md5']}\t{manifest[rel]['size_bytes']}", flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
