"""Cell state specificity of every gene from donor pseudobulk of the
reference (plan, sections 2 and 6).

For each cell with a named state the raw counts are summed per donor and
state; donor by state pseudobulks with fewer than specificity_min_cells
cells are dropped; each pseudobulk is scaled to counts per million; the
state's expression is the mean CPM over its donors (so every donor weighs
the same whatever its cell count); specificity is the state's expression
divided by the sum over states, so the specificities of a gene sum to
one across states. Cells with an unknown donor form one pseudobulk per
dataset and state.

Writes to --out-dir: specificity.tsv.gz (genes x states), expression_cpm.tsv.gz
(genes x states, mean CPM), pseudobulk_cells.tsv (donor, state, cells).
Usage: python workflow/scripts/reference_specificity.py --counts <reference_merged.h5ad> \
           --annotation <reference_annotation.tsv.gz> --out-dir <dir>
"""
import argparse
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

UNASSIGNED = "unassigned"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()
    gcfg = yaml.safe_load(open(args.config))["genetics"]
    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    ann = pd.read_csv(args.annotation, sep="\t", index_col=0, low_memory=False)
    ann = ann[~ann["state_final"].str.endswith(UNASSIGNED)]
    donor = ann["donor"].astype(str)
    donor = np.where(donor.isin(["unknown", "nan", ""]), ann["dataset"].astype(str) + ":unknown", donor)
    unit = pd.Series(donor, index=ann.index) + "|" + ann["state_final"].astype(str)
    codes, uniques = pd.factorize(unit)
    print(f"{len(ann)} cells, {len(uniques)} donor by state pseudobulks", flush=True)

    a = ad.read_h5ad(args.counts, backed="r")
    rows = a.obs_names.get_indexer(ann.index)
    order = np.argsort(rows)
    sums = np.zeros((len(uniques), a.n_vars))
    ncell = np.bincount(codes, minlength=len(uniques))
    chunk = 50000
    for i in range(0, len(order), chunk):
        sel = order[i:i + chunk]
        x = sp.csr_matrix(a[rows[sel]].to_memory().X)
        ind = sp.csr_matrix((np.ones(len(sel)), (codes[sel], np.arange(len(sel)))), shape=(len(uniques), len(sel)))
        sums += (ind @ x).toarray()
        print(f"  {min(i + chunk, len(order))} cells summed", flush=True)
    pb = pd.DataFrame(sums, index=uniques, columns=a.var_names)
    keep = ncell >= gcfg["specificity_min_cells"]
    pb = pb[keep]
    cells = pd.DataFrame({"donor": [u.split("|")[0] for u in pb.index], "state": [u.split("|")[1] for u in pb.index],
                          "cells": ncell[keep]})
    cells.to_csv(out / "pseudobulk_cells.tsv", sep="\t", index=False)
    cpm = pb.div(pb.sum(axis=1), axis=0) * 1e6
    expr = cpm.groupby(cells["state"].values).mean().T
    spec = expr.div(expr.sum(axis=1), axis=0).fillna(0.0)
    expr.to_csv(out / "expression_cpm.tsv.gz", sep="\t")
    spec.to_csv(out / "specificity.tsv.gz", sep="\t")
    print(f"{keep.sum()} pseudobulks kept, {expr.shape[1]} states, {expr.shape[0]} genes; written to {out}", flush=True)


if __name__ == "__main__":
    main()
