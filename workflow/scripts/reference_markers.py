"""Top reference markers per state, used by check (a) of the spatial
mapping. One-vs-rest Wilcoxon test (scanpy rank_genes_groups) on a
stratified subsample of the joint reference (up to N cells per state,
counts normalised per cell and log-transformed); a marker is a gene with
adjusted p below the cut-off that is detected in at least the given
fraction of the state's cells, ranked by log fold change; the top n are
kept per state.

Usage: python workflow/scripts/reference_markers.py --counts <reference_merged.h5ad> \
           --annotation <reference_annotation.tsv.gz> --out results/spatial/reference_markers.tsv
"""
import argparse

import numpy as np
import pandas as pd
import scanpy as sc
import yaml

from pnatlas.io import read_counts_float32

UNASSIGNED = "unassigned"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    p = cfg["spatial_checks"]
    rng = np.random.default_rng(cfg["seed"])

    ann = pd.read_csv(args.annotation, sep="\t", index_col=0)
    ann = ann[~ann["state_final"].str.endswith(UNASSIGNED)]
    pick = []
    for state, g in ann.groupby("state_final"):
        n = min(len(g), p["markers_cells_per_state"])
        pick.extend(rng.choice(g.index.values, n, replace=False))
    a = read_counts_float32(args.counts, pick)
    a.obs["state_final"] = ann.loc[a.obs_names, "state_final"].values
    print("cells", a.n_obs, "states", a.obs["state_final"].nunique(), flush=True)

    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.tl.rank_genes_groups(a, "state_final", method="wilcoxon", pts=True)
    rows = []
    for state in sorted(a.obs["state_final"].unique()):
        df = sc.get.rank_genes_groups_df(a, state)
        pct = a.uns["rank_genes_groups"]["pts"][state]
        df["pct"] = pct.loc[df["names"]].values
        df = df[(df["pvals_adj"] < p["markers_padj"]) & (df["pct"] >= p["markers_min_pct"]) & (df["logfoldchanges"] > 0)]
        df = df.sort_values("logfoldchanges", ascending=False).head(p["markers_top_n"])
        for rank, r in enumerate(df.itertuples(), 1):
            rows.append([state, rank, r.names, round(float(r.logfoldchanges), 3), round(float(r.pct), 3), float(r.pvals_adj)])
    out = pd.DataFrame(rows, columns=["state", "rank", "gene", "log2fc", "pct_in_state", "padj"])
    out.to_csv(args.out, sep="\t", index=False)
    print(out.groupby("state").size().to_string())


if __name__ == "__main__":
    main()
