"""Look up candidate GEO series with NCBI E-utilities and write one row per series.

Usage: python workflow/scripts/stage1_verify_geo.py GSE1 GSE2 ... > data/metadata/geo_verification.tsv
"""
import json
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def get(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.load(r)


def lookup(acc):
    q = urllib.parse.urlencode({"db": "gds", "term": f"{acc}[ACCN] AND gse[ETYP]", "retmode": "json"})
    ids = get(f"{EUTILS}/esearch.fcgi?{q}")["esearchresult"]["idlist"]
    if not ids:
        return {"accession": acc, "status": "not found"}
    s = get(f"{EUTILS}/esummary.fcgi?db=gds&id={ids[0]}&retmode=json")["result"][ids[0]]
    return {
        "accession": acc,
        "status": "found",
        "title": s.get("title", ""),
        "taxon": s.get("taxon", ""),
        "gdstype": s.get("gdstype", ""),
        "n_samples": s.get("n_samples", ""),
        "pubmed_ids": ";".join(str(p) for p in s.get("pubmedids", [])),
        "suppfile": s.get("suppfile", ""),
        "public_date": s.get("pdat", ""),
        "ftp": s.get("ftplink", ""),
    }


def main(accs):
    cols = ["accession", "status", "title", "taxon", "gdstype", "n_samples",
            "pubmed_ids", "suppfile", "public_date", "ftp"]
    print("\t".join(cols))
    for acc in accs:
        row = lookup(acc)
        print("\t".join(str(row.get(c, "")).replace("\t", " ") for c in cols))
        time.sleep(0.4)


if __name__ == "__main__":
    main(sys.argv[1:])
