"""Verify the genetics-layer sources (stage 5 inputs) against live APIs and
print one TSV row per source.

Queries: GWAS Catalog REST API v2 (study records), the eQTL Catalogue dataset
table on GitHub (the v2 API has been retired), the GTEx public bucket listing
(Google Cloud Storage JSON API) and the GTEx portal API v2 (tissue counts),
and PanelApp (Genomics England and Australia) panel records.

Usage: python workflow/scripts/stage1_verify_genetics.py > data/metadata/genetics_sources.tsv
"""
import csv
import io
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

UA = {"User-Agent": "pnatlas-stage1-check", "Accept": "application/json"}

GWAS = ["GCST90129436", "GCST90475676", "GCST90475844", "GCST90859373", "GCST007581",
        "GCST90018813", "GCST90104786", "GCST008512", "GCST008223", "GCST008530"]
EQTL_TABLE = ("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/"
              "master/tabix/tabix_ftp_paths.tsv")
GTEX_PREFIXES = [
    "bulk-qtl/v10/single-tissue-cis-qtl/GTEx_Analysis_v10_eQTL",
    "bulk-qtl/v10/single-tissue-cis-qtl/GTEx_Analysis_v10_sQTL",
    "bulk-gex/v10/rna-seq/counts-by-tissue/gene_reads_v10_nerve",
    "bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene",
    "annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_S",
]
PANELS = [("GE", "https://panelapp.genomicsengland.co.uk/api/v1/panels/", [85, 846, 79, 165]),
          ("AU", "https://panelapp-aus.org/api/v1/panels/", [3070, 3439])]
COLS = ["source", "id", "name", "detail", "n", "access", "url"]


def get(url, tries=3):
    req = urllib.request.Request(url, headers=UA)
    for i in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            if e.code < 500 or i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


def gwas_rows():
    for acc in GWAS:
        d = json.loads(get(f"https://www.ebi.ac.uk/gwas/rest/api/v2/studies/{acc}"))
        full = d.get("full_summary_stats") or ""
        yield {"source": "GWAS Catalog", "id": acc, "name": d.get("disease_trait", ""),
               "detail": f"pmid {d.get('pubmed_id')}; {d.get('initial_sample_size', '')}",
               "n": d.get("snp_count", ""),
               "access": "open, full summary statistics" if full and full != "NA" else "open, top associations only",
               "url": full if full and full != "NA" else f"https://www.ebi.ac.uk/gwas/studies/{acc}"}


def eqtl_rows():
    rows = csv.DictReader(io.StringIO(get(EQTL_TABLE)), delimiter="\t")
    for r in rows:
        if r["tissue_label"] == "tibial nerve":
            yield {"source": "eQTL Catalogue", "id": r["dataset_id"],
                   "name": f"{r['study_label']} {r['tissue_label']} {r['quant_method']}",
                   "detail": f"study {r['study_id']}; condition {r['condition_label']}",
                   "n": r["sample_size"], "access": "open", "url": r["ftp_path"]}


def gtex_rows():
    for tissue in json.loads(get("https://gtexportal.org/api/v2/dataset/tissueSiteDetail?datasetId=gtex_v10&itemsPerPage=100"))["data"]:
        if tissue["tissueSiteDetailId"] == "Nerve_Tibial":
            q = urllib.parse.urlencode({"datasetId": "gtex_v10", "tissueSiteDetailId": "Nerve_Tibial",
                                        "dataType": "RNASEQ", "itemsPerPage": 1})
            n_rna = json.loads(get(f"https://gtexportal.org/api/v2/dataset/sample?{q}"))["paging_info"]["totalNumberOfItems"]
            yield {"source": "GTEx portal", "id": "gtex_v10 Nerve_Tibial",
                   "name": "tibial nerve, GTEx v10",
                   "detail": f"eGenes {tissue['eGeneCount']}; sGenes {tissue['sGeneCount']}",
                   "n": n_rna, "access": "open (expression, QTL); donor phenotypes controlled (dbGaP phs000424)",
                   "url": "https://gtexportal.org/home/tissue/Nerve_Tibial"}
    for prefix in GTEX_PREFIXES:
        q = urllib.parse.urlencode({"prefix": prefix, "fields": "items(name,size)"})
        for item in json.loads(get(f"https://storage.googleapis.com/storage/v1/b/adult-gtex/o?{q}")).get("items", []):
            yield {"source": "GTEx bucket", "id": item["name"].rsplit("/", 1)[-1], "name": item["name"],
                   "detail": "gs://adult-gtex", "n": item["size"], "access": "open",
                   "url": f"https://storage.googleapis.com/adult-gtex/{item['name']}"}


def panel_rows():
    for label, base, ids in PANELS:
        for pid in ids:
            d = json.loads(get(f"{base}{pid}/"))
            genes = d.get("genes", [])
            green = sum(1 for g in genes if g.get("confidence_level") == "3")
            amber = sum(1 for g in genes if g.get("confidence_level") == "2")
            yield {"source": f"PanelApp {label}", "id": str(pid), "name": d.get("name", ""),
                   "detail": f"version {d.get('version')}; green {green}; amber {amber}",
                   "n": len(genes), "access": "open", "url": f"{base}{pid}/"}


def main():
    w = csv.DictWriter(sys.stdout, fieldnames=COLS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    for gen in (gwas_rows, eqtl_rows, gtex_rows, panel_rows):
        for row in gen():
            w.writerow({k: str(row.get(k, "")).replace("\t", " ") for k in COLS})


if __name__ == "__main__":
    main()
