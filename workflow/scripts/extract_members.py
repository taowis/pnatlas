"""Extract chosen members of an archive (rar, zip, tar, 7z) into a directory
with libarchive, so that no separate rar tool is needed. Members are
chosen by shell-style patterns matched against their path inside the
archive, and the leading path components are dropped as with
`tar --strip-components`. Used for the Dryad bundles of Yu 2024, whose
Xenium output folders are packed as rar archives of 2 to 5 GB each, of
which the analysis needs a few files.

Usage: python workflow/scripts/extract_members.py <archive> <out-dir> --strip 1 '<pattern>' ...
Prints one line per extracted member (path in the archive, bytes written).
"""
import argparse
import fnmatch
import pathlib
import sys

import libarchive


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("archive")
    ap.add_argument("out_dir")
    ap.add_argument("patterns", nargs="+", help="shell patterns matched against the member path")
    ap.add_argument("--strip", type=int, default=0, help="leading path components to drop")
    a = ap.parse_args()
    out = pathlib.Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    n = 0
    with libarchive.file_reader(a.archive) as arc:
        for entry in arc:
            path = entry.pathname
            if entry.isdir or not any(fnmatch.fnmatch(path, p) for p in a.patterns):
                continue
            parts = pathlib.PurePosixPath(path).parts[a.strip:]
            if not parts or ".." in parts or parts[0] == "/":
                raise SystemExit(f"refusing member path {path}")
            dest = out.joinpath(*parts)
            dest.parent.mkdir(parents=True, exist_ok=True)
            written = 0
            with open(dest, "wb") as f:
                for block in entry.get_blocks():
                    f.write(block)
                    written += len(block)
            if entry.size and written != entry.size:
                raise SystemExit(f"{path}: wrote {written} bytes, archive says {entry.size}")
            print(f"{path}\t{written}", flush=True)
            n += 1
    print(f"extracted {n} members to {out}", flush=True)
    return 0 if n else 1


if __name__ == "__main__":
    sys.exit(main())
