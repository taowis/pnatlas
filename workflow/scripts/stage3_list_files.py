"""Build the download list for the stage 3 spatial layer.

Lists the supplementary files of each GEO series from the GEO FTP directory
(over HTTPS) and the files of the two public Box folders of the UT Dallas
apps (parsed from the folder pages; each file is fetched through the
shared-link download endpoint), keeps what stage 3 needs and writes one
row per file. From the Tavares-Ferreira folder only the Space Ranger H5
files, the sample sheet and the combined Seurat object are taken (the
H and E images and Loupe files, 37 GB, are not); from the DPN nerve folder
the annotated Seurat object and the sample sheet (the Loupe files are not).

Usage: python workflow/scripts/stage3_list_files.py > data/metadata/stage3_files.tsv
"""
import html
import re
import sys
import urllib.request

UA = {"User-Agent": "Mozilla/5.0"}
GEO_HTTPS = "https://ftp.ncbi.nlm.nih.gov/geo/series"
BOX = "https://utdallas.app.box.com"

# dataset id (config/datasets.yaml), GEO series, regex of file names to keep
GEO = [
    ("heming2025_sural_xenium", "GSE285984", r"_RAW\.tar$"),
    ("yu2024_hdrg_xenium", "GSE273557", r"_RAW\.tar$|matrix.*\.csv\.gz$"),
    ("zhang2024_hdrg_visium", "GSE243073", r"_RAW\.tar$"),
    ("vroman2023_hdrg_visium", "GSE215994", r"_RAW\.tar$"),
    ("shiers2025_dpn_drg", "GSE295206", r"_RAW\.tar$|metadata.*\.docx\.gz$"),
]
# dataset id, shared link name, folder id within the link (None for the root), regex to keep
BOX_FOLDERS = [
    ("tavaresferreira2022_hdrg_visium", "dgs8a7vo9uu9cx6l81somaxb468hpgnz", None, r"\.rds$|\.xlsx$"),
    ("tavaresferreira2022_hdrg_visium", "dgs8a7vo9uu9cx6l81somaxb468hpgnz", "156483046434", r"\.h5$"),
    ("utdallas_dpn_nerves", "apafnhru6dihoi2va20wxexi40lia9b0", None, r"\.rds$|\.csv$"),
]

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


def box_rows():
    for ds, shared, folder, keep in BOX_FOLDERS:
        page = f"{BOX}/s/{shared}" + (f"/folder/{folder}" if folder else "")
        items = re.finditer(r'"typedID":"(f_\d+)","type":"file","id":\d+.*?"name":"([^"]*)","itemSize":(\d+)', fetch(page))
        for m in items:
            name = html.unescape(m.group(2))
            if re.search(keep, name):
                url = f"{BOX}/index.php?rm=box_download_shared_file&shared_name={shared}&file_id={m.group(1)}"
                yield [ds, f"Box {shared}", name, url, int(m.group(3)), "open"]


def main():
    print("\t".join(COLS))
    for row in list(geo_rows()) + list(box_rows()):
        print("\t".join(str(x) for x in row))


if __name__ == "__main__":
    sys.exit(main())
