"""Build the static gene lookup site from the result tables.

Writes <out>/index.html (a copy of portal/index.html), <out>/data/ (the JSON the
page reads: a gene index, one shard of gene records per two-letter prefix, the
state table and the site metadata) and <out>/tables/ (copies of the result
tables offered for download). Everything on the page is precomputed here, so
the site needs no server beyond a file host.
"""
import argparse
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]

PANEL_LABELS = {
    "truth:CMT": "PNPrank truth set, Charcot-Marie-Tooth",
    "truth:dHMN": "PNPrank truth set, distal hereditary motor neuropathy",
    "truth:HSAN": "PNPrank truth set, hereditary sensory and autonomic neuropathy",
    "truth:HSP": "PNPrank truth set, hereditary spastic paraplegia with neuropathy",
    "truth:SYND": "PNPrank truth set, syndromic neuropathy",
    "truth:all": "PNPrank truth set, all categories",
    "panelapp:ge_846_v8.38": "PanelApp England 846, hereditary neuropathy or pain disorder (v8.38, green)",
    "panelapp:ge_85_v1.513": "PanelApp England 85, hereditary neuropathy (v1.513, green)",
    "panelapp:aus_3070_v3.1": "PanelApp Australia 3070, hereditary neuropathy (v3.1, green)",
    "panelapp:aus_3439_v2.0": "PanelApp Australia 3439, autonomic neuropathy (v2.0, green)",
}

DOWNLOAD_TABLES = [
    ("results/backtest/pnatlas_gene_features.tsv.gz", "Gene by state features used in the PNPrank back-test (17,154 genes by 83 columns)"),
    ("results/reference/state_rule.tsv", "Reference states with cell and donor counts"),
    ("results/reference/final_state_by_dataset.tsv", "Cells per state and dataset"),
    ("results/spatial/state_mapping_status.tsv", "Spatial mapping status per state and compartment"),
    ("results/deconvolution/gtex_proportions.tsv", "GTEx tibial nerve state fractions (MuSiC, primary) with donor covariates"),
    ("results/deconvolution/gtex_bayesprism.tsv", "GTEx tibial nerve state fractions (BayesPrism)"),
    ("results/deconvolution/gtex_covariate_models.tsv", "State fraction against donor covariates"),
    ("results/genetics/gene_classes.tsv", "Neuropathy gene classes with each gene's most specific state"),
    ("results/genetics/gene_class_enrichment.tsv", "Enrichment of gene classes in state-specific genes"),
    ("results/genetics/coloc_all_pairs.tsv", "All GWAS region by eQTL gene colocalisation tests"),
    ("results/genetics/coloc_colocalised.tsv", "Colocalised pairs (PP.H4 above 0.8)"),
    ("results/backtest/backtest_summary.tsv", "PNPrank back-test with and without the pnatlas features"),
    ("results/backtest/backtest_by_cutoff.tsv", "PNPrank back-test per cutoff year"),
]


def rnd(x, d):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return None
    return round(float(x), d)


def shard_key(symbol):
    k = re.sub(r"[^A-Z0-9]", "_", symbol.upper()[:2])
    return k if len(k) == 2 else (k + "_")[:2]


def gwas_labels(datasets_yaml):
    labels = {}
    with open(datasets_yaml) as fh:
        entries = yaml.safe_load(fh)["datasets"]
    for d in entries:
        if str(d["id"]).startswith("gwas_"):
            name = d["id"][5:].replace("_", " ")
            labels[d["id"]] = f"{name} ({d['accession']})"
    return labels


def build(args):
    out = Path(args.out)
    (out / "data" / "genes").mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)

    feats = pd.read_csv(args.features, sep="\t", index_col=0)
    spec_cols = [c for c in feats.columns if c.startswith("pnatlas|spec|")]
    states = [c.split("|", 2)[2] for c in spec_cols]
    deconv_cols = [c for c in feats.columns if c.startswith("pnatlas|deconv|")]
    deconv_states = [c.split("|", 2)[2] for c in deconv_cols]

    genes = pd.read_csv(args.genes, sep="\t", index_col=0)
    classes = pd.read_csv(args.classes, sep="\t")
    classes = classes[classes["class"] != "truth:all"]
    class_by_symbol = classes.groupby("gene")["class"].apply(list).to_dict()
    truth = pd.read_csv(args.truth, sep="\t").set_index("ensembl_gene_id")
    coloc = pd.read_csv(args.coloc, sep="\t")
    labels = gwas_labels(args.datasets)
    tested = coloc.groupby("gene").size().to_dict()
    strong = coloc[coloc["pp_h4"] >= 0.5].sort_values("pp_h4", ascending=False)
    pairs = {}
    for g, grp in strong.groupby("gene"):
        pairs[g] = [
            {"gwas": labels.get(r.gwas, r.gwas), "region": r.region, "pp_h4": rnd(r.pp_h4, 3),
             "method": r.method, "lead": r.lead_variant, "lead_p": f"{r.lead_p:.1e}"}
            for r in grp.itertuples()
        ]

    index = []
    shards = {}
    for gid, row in feats.iterrows():
        symbol = genes["gene_name"].get(gid) if gid in genes.index else None
        if not isinstance(symbol, str) or not symbol:
            symbol = gid
        cls = class_by_symbol.get(symbol, [])
        key = shard_key(symbol)
        rec = {
            "symbol": symbol, "id": gid,
            "biotype": genes["gene_biotype"].get(gid, "") if gid in genes.index else "",
            "chromosome": str(genes["chromosome"].get(gid, "")) if gid in genes.index else "",
            "spec": [rnd(row[f"pnatlas|spec|{s}"], 4) for s in states],
            "logcpm": [rnd(row[f"pnatlas|logcpm|{s}"], 2) for s in states],
            "compartment": {"ganglion": rnd(row["pnatlas|compartment|ganglion"], 3),
                            "nerve": rnd(row["pnatlas|compartment|nerve"], 3)},
            "deconv": [rnd(row[f"pnatlas|deconv|{s}"], 3) for s in deconv_states],
            "coloc": {"tested": int(tested.get(gid, 0)), "max_pp_h4": rnd(row["pnatlas|coloc|max_pp_h4"], 3),
                      "pairs": pairs.get(gid, [])},
            "classes": cls,
        }
        if gid in truth.index:
            t = truth.loc[gid]
            rec["truth"] = {"category": t["disease_category"], "tier": int(t["tier"]),
                            "inheritance": t["inheritance"], "year": int(t["discovery_year"]),
                            "pmid": int(t["pmid"]) if pd.notna(t["pmid"]) else None}
        shards.setdefault(key, {})[gid] = rec
        index.append({"s": symbol, "e": gid, "k": key, "c": len(cls), "t": "truth" in rec})
    for key, recs in shards.items():
        with open(out / "data" / "genes" / f"{key}.json", "w") as fh:
            json.dump(recs, fh, separators=(",", ":"))
    index.sort(key=lambda r: r["s"].upper())
    with open(out / "data" / "index.json", "w") as fh:
        json.dump(index, fh, separators=(",", ":"))

    rule = pd.read_csv(args.state_rule, sep="\t").set_index("state_final")
    by_ds = pd.read_csv(args.by_dataset, sep="\t").set_index("state_final")
    mapping = pd.read_csv(args.mapping, sep="\t")
    music = pd.read_csv(args.music, sep="\t", index_col=0)
    bp = pd.read_csv(args.bayesprism, sep="\t", index_col=0)
    enr = pd.read_csv(args.enrichment, sep="\t")
    enr = enr[(enr["enriched"]) & (enr["class"] != "truth:all")]
    state_rows = []
    for s in states:
        r = rule.loc[s]
        m = {}
        for mr in mapping[mapping["state"] == s].itertuples():
            m[mr.compartment] = {"status": mr.status, "sections": int(mr.sections_checked),
                                 "auroc": rnd(mr.median_auroc, 3), "spearman": rnd(mr.median_spearman, 3)}
        gtex = None
        if s in music.columns:
            gtex = {"music_mean": rnd(music[s].mean(), 4), "music_zero": rnd((music[s] == 0).mean(), 3),
                    "bayesprism_mean": rnd(bp[s].mean(), 4)}
        state_rows.append({
            "state": s, "class": r["class"], "cells": int(r["cells"]), "donors": int(r["known_donors"]),
            "datasets": [d for d in by_ds.columns if s in by_ds.index and by_ds.loc[s, d] > 0],
            "mapping": m, "gtex": gtex,
            "enriched": [{"class": e["class"], "fold": rnd(e["fold"], 2), "q": f"{e['q']:.1e}"}
                         for _, e in enr[enr["state"] == s].iterrows()],
        })
    with open(out / "data" / "states.json", "w") as fh:
        json.dump({"states": state_rows, "deconv_states": deconv_states, "labels": PANEL_LABELS,
                   "datasets": list(by_ds.columns)}, fh, separators=(",", ":"))

    summary = pd.read_csv(args.backtest_summary, sep="\t")
    primary = summary[summary["primary"]].iloc[0]
    commit = args.commit
    if not commit:
        try:
            commit = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"],
                                    capture_output=True, text=True, check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            commit = "unknown"
    meta = {
        "commit": commit, "genes": len(index), "states": len(states),
        "cells": int(rule["cells"].sum()), "gtex_samples": len(music),
        "coloc_pairs": len(coloc), "coloc_colocalised": int((coloc["pp_h4"] > 0.8).sum()),
        "backtest": {"mean_A": rnd(primary["mean_A"], 3), "mean_B": rnd(primary["mean_B"], 3),
                     "p": rnd(primary["p_sign_flip"], 3), "improvement": bool(primary["improvement"])},
        "tables": [],
    }
    for rel, desc in DOWNLOAD_TABLES:
        src = ROOT / rel
        shutil.copy2(src, out / "tables" / src.name)
        meta["tables"].append({"file": src.name, "description": desc, "bytes": src.stat().st_size})
    with open(out / "data" / "meta.json", "w") as fh:
        json.dump(meta, fh, indent=1)
    shutil.copy2(Path(__file__).with_name("index.html"), out / "index.html")
    (out / ".nojekyll").touch()
    print(f"{len(index)} genes in {len(shards)} shards, {len(states)} states, written to {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    r = ROOT
    ap.add_argument("--out", default=str(r / "portal" / "site"))
    ap.add_argument("--features", default=str(r / "results/backtest/pnatlas_gene_features.tsv.gz"))
    ap.add_argument("--genes", default=str(r / "data/metadata/ensembl105_gene_length_gc.tsv.gz"))
    ap.add_argument("--classes", default=str(r / "results/genetics/gene_classes.tsv"))
    ap.add_argument("--enrichment", default=str(r / "results/genetics/gene_class_enrichment.tsv"))
    ap.add_argument("--truth", default=str(r / "data/metadata/pnprank_truth_set.tsv"))
    ap.add_argument("--coloc", default=str(r / "results/genetics/coloc_all_pairs.tsv"))
    ap.add_argument("--datasets", default=str(r / "config/datasets.yaml"))
    ap.add_argument("--state-rule", default=str(r / "results/reference/state_rule.tsv"))
    ap.add_argument("--by-dataset", default=str(r / "results/reference/final_state_by_dataset.tsv"))
    ap.add_argument("--mapping", default=str(r / "results/spatial/state_mapping_status.tsv"))
    ap.add_argument("--music", default=str(r / "results/deconvolution/gtex_music.tsv"))
    ap.add_argument("--bayesprism", default=str(r / "results/deconvolution/gtex_bayesprism.tsv"))
    ap.add_argument("--backtest-summary", default=str(r / "results/backtest/backtest_summary.tsv"))
    ap.add_argument("--commit", default="")
    build(ap.parse_args())


if __name__ == "__main__":
    main()
