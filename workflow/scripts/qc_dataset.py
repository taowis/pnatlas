"""Per-cell QC of one reference dataset: counts, genes, mitochondrial and
ribosomal fractions, a scrublet doublet call per sample, and the QC flags
from the thresholds in config/config.yaml. Writes the dataset back as h5ad
with the metrics and flags in obs and a `qc_pass` column, and prints one
summary line.

Author-annotated datasets (config qc.keep_author_cells) keep every cell and
only carry the flags; the others keep the cells with qc_pass. Whole-cell
datasets use the whole-cell mitochondrial threshold and get no doublet call.

Usage: python workflow/scripts/qc_dataset.py <dataset_id> --in "$DATA_ROOT/interim/reference/<id>.h5ad" \
           --out "$DATA_ROOT/interim/reference_qc/<id>.h5ad" [--config config/config.yaml]
"""
import argparse
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml


def doublets(adata, expected_rate, min_cells, seed):
    import scrublet as scr

    score = np.full(adata.n_obs, np.nan)
    call = np.zeros(adata.n_obs, dtype=bool)
    for s in adata.obs["sample"].unique():
        idx = np.where(adata.obs["sample"].values == s)[0]
        if len(idx) < min_cells:
            continue
        x = adata.X[idx]
        sub = scr.Scrublet(x, expected_doublet_rate=expected_rate, random_state=seed)
        sc_, pred = sub.scrub_doublets(min_counts=2, min_cells=3, n_prin_comps=30, verbose=False)
        if pred is None:  # threshold could not be set automatically
            pred = sc_ > 0.25
        score[idx], call[idx] = sc_, pred
    return score, call


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    a = ap.parse_args()
    with open(a.config) as f:
        cfg = yaml.safe_load(f)
    q, seed = cfg["qc"], cfg["seed"]

    adata = ad.read_h5ad(a.inp)
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.upper().str.match(r"^RP[SL]\d")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, inplace=True, log1p=False)
    adata.obs = adata.obs.rename(columns={"total_counts": "n_counts", "n_genes_by_counts": "n_genes",
                                          "pct_counts_mt": "pct_mt", "pct_counts_ribo": "pct_ribo"})
    whole_cell = a.dataset in q.get("whole_cell", [])
    max_mt = q["max_pct_mt_cells"] if whole_cell else q["max_pct_mt_nuclei"]
    adata.obs["flag_low_genes"] = adata.obs["n_genes"] < q["min_genes"]
    adata.obs["flag_high_mt"] = adata.obs["pct_mt"] > max_mt
    if a.dataset in q.get("no_doublet_call", []):
        adata.obs["doublet_score"], adata.obs["flag_doublet"] = np.nan, False
    else:
        adata.obs["doublet_score"], adata.obs["flag_doublet"] = doublets(
            adata, q["scrublet_expected_rate"], q["scrublet_min_cells_per_sample"], seed)
    flags = adata.obs[["flag_low_genes", "flag_high_mt", "flag_doublet"]]
    adata.obs["qc_pass"] = ~flags.any(axis=1)
    keep_all = a.dataset in q.get("keep_author_cells", [])
    n_in = adata.n_obs
    out = adata if keep_all else adata[adata.obs["qc_pass"]].copy()
    out.uns["qc"] = {"min_genes": q["min_genes"], "max_pct_mt": max_mt, "author_cells_kept": keep_all,
                     "scrublet_expected_rate": q["scrublet_expected_rate"], "seed": seed}
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(a.out, compression="gzip")
    med = adata.obs[["n_counts", "n_genes", "pct_mt"]].median()
    print("\t".join(str(x) for x in [
        a.dataset, n_in, out.n_obs, int(flags["flag_low_genes"].sum()), int(flags["flag_high_mt"].sum()),
        int(flags["flag_doublet"].sum()), int(med["n_counts"]), int(med["n_genes"]), round(float(med["pct_mt"]), 2),
        "kept_all" if keep_all else "filtered"]))


if __name__ == "__main__":
    main()
