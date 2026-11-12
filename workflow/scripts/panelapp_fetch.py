"""Fetch the gene lists of the PanelApp panels the plan names, at the panel
versions the plan fixes, and write one small table per panel under
data/metadata/ (a snapshot: the API serves a given version on request,
but the record is kept in the repository). Columns: panel, version,
gene_symbol, ensembl_gene_id (GRCh38), confidence (green, amber, red),
mode_of_inheritance, phenotypes.

Usage: python workflow/scripts/panelapp_fetch.py --out-dir data/metadata
"""
import argparse
import json
import pathlib
import urllib.request

PANELS = [  # (site, panel id, version), plan section 6
    ("https://panelapp.genomicsengland.co.uk/api/v1", 846, "8.38"),
    ("https://panelapp.genomicsengland.co.uk/api/v1", 85, "1.513"),
    ("https://panelapp-aus.org/api/v1", 3070, "3.1"),
    ("https://panelapp-aus.org/api/v1", 3439, "2.0"),
]
CONF = {"3": "green", "2": "amber", "1": "red", "0": "none"}


def get(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "pnatlas-stage5"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="data/metadata")
    args = ap.parse_args()
    out = pathlib.Path(args.out_dir)
    for site, pid, version in PANELS:
        panel = get(f"{site}/panels/{pid}/?version={version}")
        if str(panel["version"]) != version:
            raise SystemExit(f"panel {pid}: asked for version {version}, got {panel['version']}")
        host = "ge" if "genomicsengland" in site else "aus"
        rows = []
        for g in panel["genes"]:
            gd = g["gene_data"]
            ens = gd.get("ensembl_genes", {}).get("GRch38", {})
            ens_id = next(iter(ens.values()))["ensembl_id"] if ens else ""
            rows.append([f"{host}:{pid}", version, gd["gene_symbol"], ens_id, CONF.get(str(g["confidence_level"]), ""),
                         g.get("mode_of_inheritance", ""), "; ".join(g.get("phenotypes", []))])
        path = out / f"panelapp_{host}_{pid}_v{version}.tsv"
        with open(path, "w") as f:
            f.write("panel\tversion\tgene_symbol\tensembl_gene_id\tconfidence\tmode_of_inheritance\tphenotypes\n")
            for r in rows:
                f.write("\t".join(str(x).replace("\t", " ").replace("\n", " ") for x in r) + "\n")
        green = sum(r[4] == "green" for r in rows)
        print(f"{path}: {panel['name']} v{panel['version']}, {len(rows)} genes, {green} green", flush=True)


if __name__ == "__main__":
    main()
