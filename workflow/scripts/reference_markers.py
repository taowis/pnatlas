"""Top reference markers per state, used by check (a) of the spatial
mapping. One-vs-rest Wilcoxon test (scanpy rank_genes_groups) on a
stratified subsample of the joint reference (up to N cells per state,
counts normalised per cell and log-transformed); a marker is a gene with
adjusted p below the cut-off that is detected in at least the given
fraction of the state's cells, ranked by log fold change; the top n are
kept per state. With --genes the test is restricted to the listed genes
(one per line), for a targeted panel.

Usage: python workflow/scripts/reference_markers.py --counts <reference_merged.h5ad> \
           --annotation <reference_annotation.tsv.gz> --out results/spatial/reference_markers.tsv
"""
import argparse

import numpy as np
import pandas as pd
import yaml

from pnatlas.io import read_counts_float32
from pnatlas.markers import subsample_per_label, top_markers

UNASSIGNED = "unassigned"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--genes", default=None, help="file with one gene per line; restrict the test to these")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    p = cfg["spatial_checks"]
    rng = np.random.default_rng(cfg["seed"])

    ann = pd.read_csv(args.annotation, sep="\t", index_col=0)
    ann = ann[~ann["state_final"].str.endswith(UNASSIGNED)]
    pick = subsample_per_label(ann["state_final"], p["markers_cells_per_state"], rng)
    a = read_counts_float32(args.counts, pick)
    a.obs["state_final"] = ann.loc[a.obs_names, "state_final"].values
    if args.genes:
        genes = [g.strip() for g in open(args.genes) if g.strip()]
        a = a[:, a.var_names.isin(genes)].copy()
    print("cells", a.n_obs, "genes", a.n_vars, "states", a.obs["state_final"].nunique(), flush=True)

    out = top_markers(a, "state_final", p["markers_padj"], p["markers_min_pct"], p["markers_top_n"])
    out.to_csv(args.out, sep="\t", index=False)
    print(out.groupby("state").size().to_string())


if __name__ == "__main__":
    main()
