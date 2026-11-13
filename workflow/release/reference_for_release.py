"""Write the reference object for the data archive.

Takes the merged counts (interim/reference_merged.h5ad), the integration
embeddings (interim/reference_embeddings.h5ad) and the final annotation
(interim/reference_annotation.tsv.gz), drops the datasets whose cell-level
data may not be redistributed (--exclude, by dataset id), keeps the raw counts
and a short set of cell columns, and writes one compressed h5ad plus the cell
table. Cells are read in backed mode so that only the kept cells are loaded.
"""
import argparse

import anndata as ad
import h5py
import numpy as np
import pandas as pd

OBS_KEEP = ["dataset", "sample", "donor", "batch", "tissue", "compartment", "protocol",
            "author_label", "class_final", "state_final", "label_source", "transfer_conf",
            "n_counts", "n_genes", "pct_mt"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--merged", required=True)
    ap.add_argument("--embeddings", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--exclude", nargs="*", default=[], help="dataset ids to leave out")
    ap.add_argument("--out", required=True, help="output .h5ad")
    ap.add_argument("--cells-out", required=True, help="output cell table .tsv.gz")
    a = ap.parse_args()

    ann = pd.read_csv(a.annotation, sep="\t", index_col=0)
    merged = ad.read_h5ad(a.merged, backed="r")
    keep = merged.obs.index.isin(ann.index) & ~merged.obs["dataset"].isin(a.exclude).values
    print(f"{keep.sum()} of {merged.n_obs} cells kept; excluded datasets: {a.exclude or 'none'}")
    sub = merged[np.where(keep)[0]].to_memory()
    merged.file.close()
    obs = sub.obs.join(ann[["class_final", "state_final", "label_source", "transfer_conf"]], how="left")
    sub.obs = obs[[c for c in OBS_KEEP if c in obs.columns]].copy()
    for c in sub.obs.columns:
        if sub.obs[c].dtype == object:
            sub.obs[c] = sub.obs[c].astype(str)
    with h5py.File(a.embeddings) as h:
        idx = h["obs"]["_index"][:].astype(str)
        pos = pd.Series(np.arange(len(idx)), index=idx).reindex(sub.obs.index)
        assert pos.notna().all(), "cells missing from the embeddings file"
        order = pos.to_numpy().astype(int)
        for k in ["X_scanvi", "X_scvi"]:
            sub.obsm[k] = h["obsm"][k][:][order]
    sub.uns = {"source": "pnatlas reference, raw counts after quality control; "
                         "state_final is the reference state, label_source how it was assigned"}
    sub.write_h5ad(a.out, compression="gzip")
    sub.obs.to_csv(a.cells_out, sep="\t", compression="gzip")
    print(f"written {sub.n_obs} cells x {sub.n_vars} genes to {a.out}")


if __name__ == "__main__":
    main()
