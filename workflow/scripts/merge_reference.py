"""Merge the QC'd reference datasets into one AnnData of raw counts.

Applies config/celltype_map.tsv (every author label must have a row),
removes the painseq atlas cells that duplicate datasets read separately,
keeps the genes present in every dataset, and adds the columns used
downstream: tissue, compartment, protocol, study and batch (donor, or
library where the donor is unknown). Only the shared obs columns are kept;
the authors' extra columns stay in the per-dataset QC files.

Usage: python workflow/scripts/merge_reference.py --qc-dir "$DATA_ROOT/interim/reference_qc" \
           --out "$DATA_ROOT/interim/reference_merged.h5ad" --summary-dir results/reference
"""
import argparse
import pathlib

import anndata as ad
import pandas as pd
import yaml

from pnatlas.labels import LEVELS, UNLABELLED, apply_celltype_map, read_celltype_map

KEEP_OBS = ["dataset", "sample", "donor", "author_label", "author_label_fine", "n_counts", "n_genes",
            "pct_mt", "pct_ribo", "flag_low_genes", "flag_high_mt", "doublet_score", "flag_doublet",
            "qc_pass"]


def load_one(path, ds, cfg, cmap):
    a = ad.read_h5ad(path)
    n0 = a.n_obs
    if ds == "bhuiyan2025_hdrg":
        # amendment 1: drop the atlas cells that come from datasets read here
        a = a[~a.obs["author_Dataset"].isin(cfg["exclude_painseq_sources"])].copy()
        # cells without a library name are grouped by their atlas source study
        no_lib = a.obs["sample"].astype(str) == "SeuratProject"
        a.obs["sample"] = a.obs["sample"].astype(str).where(~no_lib, a.obs["author_Dataset"].astype(str))
        a.obs["donor"] = "unknown"
    obs = a.obs[KEEP_OBS].copy()
    for c in ["sample", "donor", "author_label", "author_label_fine"]:
        obs[c] = obs[c].astype(str)
    obs[LEVELS] = apply_celltype_map(obs, cmap)[LEVELS]
    if ds in cfg.get("labels_unusable", []):
        # the authors' labels are kept in author_label* but do not seed anything
        obs[LEVELS] = UNLABELLED
    info = cfg["datasets"][ds]
    obs["tissue"], obs["compartment"], obs["protocol"] = info["tissue"], info["compartment"], info["protocol"]
    known = obs["donor"] != "unknown"
    obs["batch"] = ds + ":" + obs["donor"].where(known, obs["sample"])
    a.obs = obs
    a.var = a.var[[]]
    a.obs_names = ds + ":" + a.obs_names.astype(str)
    return a, n0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--summary-dir", default="results/reference")
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--map", default="config/celltype_map.tsv")
    a = ap.parse_args()
    with open(a.config) as f:
        cfg = yaml.safe_load(f)["integration"]
    cmap = read_celltype_map(a.map)

    parts, rows = {}, []
    for ds in cfg["datasets"]:
        adata, n0 = load_one(pathlib.Path(a.qc_dir) / f"{ds}.h5ad", ds, cfg, cmap)
        parts[ds] = adata
        rows.append({"dataset": ds, "cells_in": n0, "cells_kept": adata.n_obs, "genes": adata.n_vars,
                     "donors": adata.obs.loc[adata.obs["donor"] != "unknown", "donor"].nunique(),
                     "batches": adata.obs["batch"].nunique(),
                     "labelled": int((adata.obs["state"] != "unlabelled").sum())})
    genes = sorted(set.intersection(*(set(p.var_names) for p in parts.values())))
    for ds in parts:
        parts[ds] = parts[ds][:, genes].copy()
    merged = ad.concat(list(parts.values()), join="inner", merge="same")
    for c in ["dataset", "sample", "donor", "batch", "tissue", "compartment", "protocol"] + LEVELS:
        merged.obs[c] = merged.obs[c].astype("category")
    merged.uns["merge"] = {"genes": cfg["genes"], "n_genes": len(genes),
                           "exclude_painseq_sources": cfg["exclude_painseq_sources"],
                           "labels_unusable": cfg.get("labels_unusable", [])}
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    merged.write_h5ad(a.out, compression="gzip")

    out = pathlib.Path(a.summary_dir)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out / "merged_summary.tsv", sep="\t", index=False)
    pd.crosstab(merged.obs["state"], merged.obs["dataset"]).to_csv(out / "merged_state_by_dataset.tsv", sep="\t")
    pd.crosstab(merged.obs["class"], merged.obs["dataset"]).to_csv(out / "merged_class_by_dataset.tsv", sep="\t")
    print(f"merged\t{merged.n_obs}\t{merged.n_vars}\t{merged.obs['batch'].nunique()} batches\t"
          f"{merged.obs['donor'].nunique() - 1} known donors")


if __name__ == "__main__":
    main()
