"""Build the download list for the stage 2 reference layer.

Lists the supplementary files of each GEO series from the GEO FTP directory
(over HTTPS), the files of SPARC dataset 465 from the Pennsieve discover API,
and the painseq Seurat object, keeps the files stage 2 needs (human count
matrices, author metadata, SDS sheets) and writes one row per file.

Usage: python workflow/scripts/stage2_list_files.py > data/metadata/stage2_files.tsv
"""
import json
import re
import sys
import urllib.request

UA = {"User-Agent": "pnatlas-stage2"}
GEO_HTTPS = "https://ftp.ncbi.nlm.nih.gov/geo/series"
S3 = "https://sparc-prod-aod-discover-publish50-use1.s3.amazonaws.com/465/"

# dataset id (config/datasets.yaml), GEO series, regex of file names to keep
GEO = [
    ("nguyen2021_hdrg", "GSE168243", r"_RAW\.tar$"),
    ("jung2023_hdrg", "GSE201586", r"_RAW\.tar$|humanallcelltypes\.rds\.gz$"),
    ("yu2024_singlesoma", "GSE249746", r"\.csv\.gz$"),
    ("zhang_hdrg_glia", "GSE189501", r"_DRG_.*\.txt\.gz$"),
    ("yang2022_htg", "GSE197289", r"snRNA-seq_human_.*\.gz$"),
    ("heming2025_sural_snrna", "GSE285983", r"_RAW\.tar$|metadata.*\.csv\.gz$"),
]
SPARC_KEEP = re.compile(r"^files/(primary/filtered_feature_bc_matrix/.*|primary/meta\.data\.csv|"
                        r"subjects\.xlsx|samples\.xlsx|README\.txt|dataset_description\.xlsx)$")
PAINSEQ = ("bhuiyan2025_hdrg", "humandrg.painseq.com", "seurat_object_preprint_release.Rds",
           ("https://www.dropbox.com/scl/fi/wx6472hewabt1vvafxet5/seurat_object_preprint_release.Rds"
            "?rlkey=struflwt1ql68okohdcvzosq1&dl=1"), "portal")

COLS = ["dataset", "accession", "file", "url", "size_bytes", "access"]


def fetch(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode("utf-8", "replace")


def size_bytes(s):
    units = {"K": 1 << 10, "M": 1 << 20, "G": 1 << 30}
    return int(float(s[:-1]) * units[s[-1]]) if s[-1] in units else int(s)


def geo_rows():
    for ds, gse, keep in GEO:
        base = f"{GEO_HTTPS}/{gse[:-3]}nnn/{gse}/suppl/"
        for m in re.finditer(r'href="([^"?/]+)">[^<]*</a>\s+\S+ \S+\s+(\S+)', fetch(base)):
            name, size = m.group(1), m.group(2)
            if re.search(keep, name):
                yield [ds, gse, name, base + name, size_bytes(size), "open"]


def sparc_rows():
    d = json.loads(fetch("https://api.pennsieve.io/discover/datasets/465"))
    v = d["version"]

    def walk(path=""):
        r = json.loads(fetch(f"https://api.pennsieve.io/discover/datasets/465/versions/{v}"
                             f"/files/browse?path={path}&limit=200"))
        for f in r.get("files", []):
            if f["type"] == "Directory":
                yield from walk(f["path"])
            elif SPARC_KEEP.match(f["path"]):
                yield ["yang2025_hsg_multiome", f"SPARC 465 v{v}", f["path"], S3 + f["path"],
                       f["size"], "open"]
    yield from walk()


def painseq_rows():
    ds, acc, name, url, access = PAINSEQ
    req = urllib.request.Request(url, headers=UA, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        size = r.headers.get("original-content-length") or r.headers.get("content-length", "")
    yield [ds, acc, name, url, size, access]


def main():
    print("\t".join(COLS))
    for gen in (geo_rows, sparc_rows, painseq_rows):
        for row in gen():
            print("\t".join(str(x) for x in row))


if __name__ == "__main__":
    sys.exit(main())
