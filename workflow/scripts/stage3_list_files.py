"""Build the download list for the stage 3 spatial layer.

Lists the supplementary files of each GEO series from the GEO FTP directory
(over HTTPS) and the files of the two public Box folders of the UT Dallas
apps (parsed from the folder pages; each file is fetched through the
shared-link download endpoint), keeps what stage 3 needs and writes one
row per file. From the Tavares-Ferreira folder only the Space Ranger H5
files, the sample sheet and the combined Seurat object are taken (the
H and E images and Loupe files, 37 GB, are not); from the DPN nerve folder
the annotated Seurat object and the sample sheet (the Loupe files are not).
With --images the list holds instead the files for the anatomy check: the
16 H and E images (TIFF, 26.5 GB) of the Tavares-Ferreira Box folder, and
the Space Ranger high-resolution tissue image and scale factors of each
Shiers 2025 section from the SPARC copy of that dataset (Pennsieve
discover dataset 490, version 1, DOI 10.26275/vz0c-uet7, CC BY), read
from the dataset's file list and fetched from its public S3 bucket; the
Shiers files are renamed <section>_<name> so that the 16 sections do not
collide.

Usage: python workflow/scripts/stage3_list_files.py > data/metadata/stage3_files.tsv
       python workflow/scripts/stage3_list_files.py --images > data/metadata/stage3_images.tsv
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
BOX_IMAGES = [
    ("tavaresferreira2022_hdrg_visium", "dgs8a7vo9uu9cx6l81somaxb468hpgnz", "206669782611", r"\.tif$"),
]
# SPARC copy of the Shiers 2025 dataset: discover dataset id, version, S3 bucket,
# files to take (regex on the path) and the sample folder to config section map
SPARC_ID, SPARC_VERSION = 490, 1
SPARC_S3 = "https://sparc-prod-aod-discover-publish50-use1.s3.amazonaws.com"
SPARC_KEEP = r"/spatial/(tissue_hires_image[^/]*\.png|scalefactors[^/]*\.json)$"
SPARC_SECTIONS = {
    "sam-UTD-DN0003-DRG-4": "M2-1", "sam-UTD-DN0003-DRG-8": "M2-2",
    "sam-UTD-180714DHB-DRGL4-1": "180714-L4-1", "sam-UTD-180714DHB-DRGL4-5": "180714-L4-5",
    "sam-UTD-181030DHB-1L4-3": "181030-1L4-3", "sam-UTD-181030DHB-1L4-7": "181030-1L4-7",
    "sam-UTD-181030DHB-2L4-2": "181030-2L4-2", "sam-UTD-181030DHB-2L4-6": "181030-2L4-6",
    "sam-UTD-DN0134-DRGL3-S1-1": "DN0134-S1-1", "sam-UTD-DN0134-DRGL3-S1-2": "DN0134-S1-2",
    "sam-UTD-DN0134-DRGL3-S1-3": "DN0134-S1-3",
    "sam-UTD-DN0095-DRGL3-S2-1": "DN0095-S2-1", "sam-UTD-DN0095-DRGL3-S2-2": "DN0095-S2-2",
    "sam-UTD-DN0095-DRGL3-S2-3": "DN0095-S2-3",
    "sam-UTD-DN0106-DRGL-S5-1": "DN0106-S5-1", "sam-UTD-DN0106-DRGL-S5-2": "DN0106-S5-2",
}

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


def box_rows(folders):
    for ds, shared, folder, keep in folders:
        page = f"{BOX}/s/{shared}" + (f"/folder/{folder}" if folder else "")
        items = re.finditer(r'"typedID":"(f_\d+)","type":"file","id":\d+.*?"name":"([^"]*)","itemSize":(\d+)', fetch(page))
        for m in items:
            name = html.unescape(m.group(2))
            if re.search(keep, name):
                url = f"{BOX}/index.php?rm=box_download_shared_file&shared_name={shared}&file_id={m.group(1)}"
                yield [ds, f"Box {shared}", name, url, int(m.group(3)), "open"]


def sparc_rows():
    import json
    meta = json.loads(fetch(f"https://api.pennsieve.io/discover/datasets/{SPARC_ID}/versions/{SPARC_VERSION}/metadata"))
    for f in meta["files"]:
        path = f["path"]
        if not re.search(SPARC_KEEP, path):
            continue
        section = SPARC_SECTIONS[path.split("/")[-3]]
        name = f"{section}_{path.rsplit('/', 1)[1]}"
        yield ["shiers2025_dpn_drg", f"SPARC {SPARC_ID} v{SPARC_VERSION}", name,
               f"{SPARC_S3}/{SPARC_ID}/{path}", int(f["size"]), "open"]


def main():
    images = "--images" in sys.argv[1:]
    print("\t".join(COLS))
    rows = list(box_rows(BOX_IMAGES)) + list(sparc_rows()) if images else list(geo_rows()) + list(box_rows(BOX_FOLDERS))
    for row in rows:
        print("\t".join(str(x) for x in row))


if __name__ == "__main__":
    sys.exit(main())
