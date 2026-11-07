"""Check that the authors' labels of a dataset follow the expression of the
cells they are attached to. Two numbers per dataset, computed on the
dataset alone (2,000 highly variable genes, 30 PCs, at most 50,000 cells):
the adjusted Rand index between k-means clusters (k = number of author
labels) and the author labels, and the fraction of a cell's 15 nearest
neighbours that carry the same author label. Labels that match the cells
give values well above zero; labels attached to the wrong cells give an
ARI near zero and a neighbour agreement equal to what chance gives.

Usage: python workflow/scripts/check_label_alignment.py --qc-dir "$DATA_ROOT/interim/reference_qc" \
           --out results/reference/label_alignment.tsv <dataset_id> [<dataset_id> ...]
"""
import argparse
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors


def check(path, seed, max_cells=50000):
    a = ad.read_h5ad(path)
    lab = a.obs["author_label_fine"].astype(str)
    a = a[lab != ""].copy()
    lab = a.obs["author_label_fine"].astype(str).values
    if a.n_obs > max_cells:
        idx = np.random.default_rng(seed).choice(a.n_obs, max_cells, replace=False)
        a, lab = a[idx].copy(), lab[idx]
    sc.pp.highly_variable_genes(a, flavor="seurat_v3", n_top_genes=2000)
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.pca(a, n_comps=30, mask_var="highly_variable", random_state=seed)
    k = len(np.unique(lab))
    km = KMeans(k, random_state=seed, n_init=3).fit_predict(a.obsm["X_pca"])
    nn = NearestNeighbors(n_neighbors=16).fit(a.obsm["X_pca"])
    _, nbr = nn.kneighbors(a.obsm["X_pca"])
    same = (lab[nbr[:, 1:]] == lab[:, None]).mean()
    p = pd.Series(lab).value_counts(normalize=True).values
    return {"cells": a.n_obs, "labels": k, "ari_kmeans": round(float(adjusted_rand_score(km, lab)), 3),
            "knn_same_label": round(float(same), 3), "knn_same_label_chance": round(float((p ** 2).sum()), 3)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("datasets", nargs="+")
    ap.add_argument("--qc-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seed", type=int, default=20261102)
    a = ap.parse_args()
    rows = []
    for ds in a.datasets:
        r = check(pathlib.Path(a.qc_dir) / f"{ds}.h5ad", a.seed)
        rows.append({"dataset": ds, **r})
        print("\t".join(str(v) for v in rows[-1].values()), flush=True)
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
