"""Check access terms of the non-GEO candidate sources (SPARC, Dryad, painseq/Dropbox).

Queries DataCite for each DOI, the Pennsieve discover API for SPARC datasets, the
Dryad API for Dryad datasets, and sends a HEAD request to each Dropbox link to get
the file size. Nothing is downloaded except response headers.

Usage: python workflow/scripts/stage1_verify_nongeo.py > data/metadata/nongeo_access.tsv
"""
import json
import sys
import urllib.parse
import urllib.request

UA = {"User-Agent": "pnatlas-stage1-check"}

SPARC = [
    ("sparc_sg_multiome", "10.26275/oqxx-si9h", 465),
    ("sparc_drg_multiome", "10.26275/sojr-qhvm", 464),
]
DRYAD = [("dryad_yu_xenium", "10.5061/dryad.gf1vhhmxq")]
DROPBOX = [
    ("painseq_humandrg_seurat",
     ("https://www.dropbox.com/scl/fi/wx6472hewabt1vvafxet5/seurat_object_preprint_release.Rds"
      "?rlkey=struflwt1ql68okohdcvzosq1&dl=1"),
     "https://humandrg.painseq.com"),
    ("painseq_harmonized_folder",
     ("https://www.dropbox.com/scl/fo/x9c173p8p0wxze15c29ln/AIpoBvGDp28rS1aDLUnzdCk"
      "?rlkey=gsg6vcp5sm8e88s36ffml4uur&dl=1"),
     "https://harmonized.painseq.com"),
]
COLS = ["id", "source", "doi", "url", "title", "licence", "embargo", "embargo_release",
        "size_bytes", "n_files", "published", "related"]


def get(url):
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def head_size(url):
    req = urllib.request.Request(url, headers=UA, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.headers.get("original-content-length") or r.headers.get("content-length", "")


def datacite(doi):
    a = get(f"https://api.datacite.org/dois/{doi}")["data"]["attributes"]
    rights = ";".join(r.get("rightsIdentifier") or r.get("rights", "") for r in a.get("rightsList", []))
    rel = ";".join(f"{r.get('relationType')}:{r.get('relatedIdentifier')}"
                   for r in a.get("relatedIdentifiers", []))
    return a["titles"][0]["title"], a.get("url", ""), rights, rel


def sparc_rows():
    for rid, doi, ds in SPARC:
        title, url, rights, rel = datacite(doi)
        d = get(f"https://api.pennsieve.io/discover/datasets/{ds}")
        ext = ";".join(f"{e['relationshipType']}:{e['doi']}" for e in d.get("externalPublications", []))
        yield {"id": rid, "source": "SPARC", "doi": f"{doi} (current: {d.get('doi')})", "url": url,
               "title": title, "licence": d.get("license") or rights,
               "embargo": d.get("embargo"), "embargo_release": d.get("embargoReleaseDate") or "",
               "size_bytes": d.get("size"), "n_files": d.get("fileCount"),
               "published": d.get("versionPublishedAt", ""), "related": ext or rel}


def dryad_rows():
    for rid, doi in DRYAD:
        title, url, rights, rel = datacite(doi)
        d = get("https://datadryad.org/api/v2/datasets/" + urllib.parse.quote(f"doi:{doi}", safe=""))
        v = d["_links"]["stash:version"]["href"]
        files = get(f"https://datadryad.org{v}/files")["_embedded"]["stash:files"]
        rel = ";".join(f"{r['relationship']}:{r['identifier']}" for r in d.get("relatedWorks", []))
        yield {"id": rid, "source": "Dryad", "doi": doi, "url": url, "title": title,
               "licence": d.get("license", rights), "embargo": False, "embargo_release": "",
               "size_bytes": d.get("storageSize"), "n_files": len(files),
               "published": d.get("publicationDate", ""), "related": rel}


def dropbox_rows():
    for rid, link, portal in DROPBOX:
        yield {"id": rid, "source": "painseq portal (Dropbox)", "doi": "", "url": portal,
               "title": link.split("?")[0].rsplit("/", 1)[-1], "licence": "none stated on portal",
               "embargo": False, "embargo_release": "", "size_bytes": head_size(link),
               "n_files": "", "published": "", "related": link.split("?")[0]}


def main():
    print("\t".join(COLS))
    for gen in (sparc_rows, dryad_rows, dropbox_rows):
        for row in gen():
            print("\t".join(str(row.get(c, "")).replace("\t", " ") for c in COLS))


if __name__ == "__main__":
    sys.exit(main())
