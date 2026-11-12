"""Prepare the inputs of the stage 4 deconvolution (plan, section 5).

From the merged reference and its annotation, take the cells of the nerve
compartment with a named state, keep the states with enough cells, and
hold a fixed number of donors out. The remaining donors give the
reference handed to the methods: up to max_cells_per_state cells per
state, drawn at random. The held-out donors give validation (a): each
pseudo-bulk mixture draws its state proportions from a flat Dirichlet,
picks one held-out donor, draws cells_per_mixture cells of that donor
with those proportions (a state the donor lacks is drawn from the other
held-out donors) and sums their counts; the realised cell fractions are
the truth. The GTEx tibial nerve read counts are read from the GCT file,
mapped to gene symbols (the Description column; duplicates summed) and
restricted, with everything else, to the genes shared with the reference.

Writes to --out-dir: genes.txt, reference_counts.mtx (cells x genes,
integer), reference_cells.tsv (cell, state, donor), donors.tsv (donor,
cells, role), mixtures_counts.tsv.gz (genes x mixtures),
mixtures_truth.tsv (mixture, donor, one column per state),
gtex_counts.tsv.gz (genes x samples), prepare_summary.tsv.
Usage: python workflow/scripts/deconv_prepare.py --counts <reference_merged.h5ad> \
           --annotation <reference_annotation.tsv.gz> --gtex-counts <gene_reads_v10_nerve_tibial.gct.gz> --out-dir <dir>
"""
import argparse
import gzip
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import yaml

UNASSIGNED = "unassigned"


def read_gct(path):
    with gzip.open(path, "rt") as f:
        f.readline()
        f.readline()
        df = pd.read_csv(f, sep="\t")
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--gtex-counts", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--n-mixtures", type=int, default=None, help="override, for a probe run")
    ap.add_argument("--max-cells-per-state", type=int, default=None, help="override, for a probe run")
    ap.add_argument("--gtex-samples", type=int, default=0, help="keep only the first N GTEx samples, for a probe run")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    dc = cfg["deconvolution"]
    rng = np.random.default_rng(cfg["seed"])
    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    n_mix = args.n_mixtures or dc["n_mixtures"]
    max_cells = args.max_cells_per_state or dc["max_cells_per_state"]

    ann = pd.read_csv(args.annotation, sep="\t", index_col=0, low_memory=False)
    ann = ann[(ann["compartment"] == dc["compartment"]) & ~ann["state_final"].str.endswith(UNASSIGNED)]
    sizes = ann["state_final"].value_counts()
    states = sorted(sizes[sizes >= dc["min_cells_per_state"]].index)
    ann = ann[ann["state_final"].isin(states)]
    ann["donor"] = ann["donor"].astype(str)
    donors = sorted(ann["donor"].unique())
    held = sorted(str(d) for d in rng.choice(donors, size=dc["holdout_donors"], replace=False))
    print(f"nerve cells {len(ann)}, states {len(states)}, donors {len(donors)}, held out {held}", flush=True)

    # reference cells: up to max_cells per state from the reference donors
    ref_ann = ann[~ann["donor"].isin(held)]
    keep = []
    for s in states:
        idx = ref_ann.index[ref_ann["state_final"] == s].to_numpy()
        keep.extend(rng.choice(idx, size=min(max_cells, len(idx)), replace=False))
    keep = sorted(keep)
    held_ann = ann[ann["donor"].isin(held)]

    a = ad.read_h5ad(args.counts, backed="r")
    genes_ref = a.var_names.to_numpy()
    gtex = read_gct(args.gtex_counts)
    sym = gtex["Description"].astype(str)
    gtex = gtex.drop(columns=["Name", "Description"]).groupby(sym.values).sum()
    if args.gtex_samples:
        gtex = gtex.iloc[:, :args.gtex_samples]
    common = sorted(set(genes_ref) & set(gtex.index))
    gpos = pd.Index(genes_ref).get_indexer(common)
    print(f"genes: reference {len(genes_ref)}, GTEx symbols {len(gtex)}, shared {len(common)}", flush=True)

    def load(cells):
        rows = a.obs_names.get_indexer(cells)
        order = np.argsort(rows)
        x = a[np.sort(rows)].to_memory().X
        x = sp.csr_matrix(x)[:, gpos]
        inv = np.empty_like(order)
        inv[order] = np.arange(len(order))
        return x[inv]

    ref_x = load(keep).astype(np.int64)
    scipy.io.mmwrite(str(out / "reference_counts.mtx"), ref_x, field="integer")
    pd.DataFrame({"cell": keep, "state": ref_ann.loc[keep, "state_final"].values,
                  "donor": ref_ann.loc[keep, "donor"].values}).to_csv(out / "reference_cells.tsv", sep="\t", index=False)
    pd.Series(common).to_csv(out / "genes.txt", index=False, header=False)
    print(f"reference written: {ref_x.shape[0]} cells", flush=True)

    # validation mixtures from the held-out donors
    held_x = load(held_ann.index.to_numpy()).tocsr()
    hs = held_ann["state_final"].to_numpy()
    hd = held_ann["donor"].to_numpy()
    mix, truth = [], []
    for m in range(n_mix):
        p = rng.dirichlet(np.full(len(states), dc["dirichlet_alpha"]))
        donor = held[m % len(held)]
        n_per = rng.multinomial(dc["cells_per_mixture"], p)
        picked = []
        for s, n in zip(states, n_per):
            if n == 0:
                continue
            pool = np.where((hs == s) & (hd == donor))[0]
            if len(pool) == 0:
                pool = np.where(hs == s)[0]
            picked.append(rng.choice(pool, size=n, replace=True))
        picked = np.concatenate(picked)
        mix.append(np.asarray(held_x[picked].sum(axis=0)).ravel())
        truth.append([f"mix{m + 1:03d}", donor] + list(n_per / n_per.sum()))
    pd.DataFrame(np.array(mix).T, index=common, columns=[t[0] for t in truth]).to_csv(
        out / "mixtures_counts.tsv.gz", sep="\t")
    pd.DataFrame(truth, columns=["mixture", "donor"] + states).to_csv(out / "mixtures_truth.tsv", sep="\t", index=False)
    gtex.loc[common].astype(np.int64).to_csv(out / "gtex_counts.tsv.gz", sep="\t")

    roles = pd.DataFrame({"donor": donors, "cells": [int((ann["donor"] == d).sum()) for d in donors],
                          "role": ["held out" if d in held else "reference" for d in donors]})
    roles.to_csv(out / "donors.tsv", sep="\t", index=False)
    summ = pd.DataFrame({"state": states, "cells_nerve": [int(sizes[s]) for s in states],
                         "cells_reference": [int((ref_ann.loc[keep, "state_final"] == s).sum()) for s in states],
                         "cells_held_out": [int((hs == s).sum()) for s in states]})
    summ.to_csv(out / "prepare_summary.tsv", sep="\t", index=False)
    print(f"mixtures {n_mix}, GTEx samples {gtex.shape[1]}; written to {out}", flush=True)


if __name__ == "__main__":
    main()
