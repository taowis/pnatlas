"""Estimate the reference cell-state signatures for cell2location from the
joint ganglia and nerve reference: a negative binomial regression of the
raw counts on the cell state (state_final of the stage 2 annotation), with
the donor (batch) and the study (dataset) as batch covariates
(cell2location RegressionModel). Cells without a named state ("unassigned"
and "<class> unassigned") are left out. Genes are filtered with the
cell2location defaults before the fit.

Writes to --out-dir:
  reference_signatures.csv.gz   genes x states, the posterior mean of the
                                per-state expression (means_per_cluster_mu_fg)
  signature_cells.tsv           cells per state, per compartment, used in the fit
  training_history.tsv          ELBO per epoch
  model/                        the saved scvi model
Usage: python workflow/scripts/c2l_reference_signatures.py --counts <reference_merged.h5ad> \
           --annotation <reference_annotation.tsv.gz> --out-dir <dir> [--max-epochs 250] [--subsample N]
"""
import argparse
import pathlib

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

from pnatlas.io import read_counts_float32

UNASSIGNED = "unassigned"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--max-epochs", type=int, default=None)
    ap.add_argument("--subsample", type=int, default=0, help="cells, for a probe run; 0 = all")
    ap.add_argument("--seed", type=int, default=None)
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    c2l = cfg["cell2location"]
    seed = args.seed if args.seed is not None else cfg["seed"]
    max_epochs = args.max_epochs or c2l["signature_max_epochs"]
    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    import scvi
    import cell2location
    from cell2location.models import RegressionModel
    scvi.settings.seed = seed

    ann = pd.read_csv(args.annotation, sep="\t", index_col=0)
    keep = ~ann["state_final"].str.endswith(UNASSIGNED)
    ann = ann[keep]
    if args.subsample:
        ann = ann.sample(n=min(args.subsample, len(ann)), random_state=seed)
    print("cells with a named state", keep.sum(), "used", len(ann), flush=True)

    a = read_counts_float32(args.counts, ann.index)
    a.obs = a.obs.join(ann[["state_final", "compartment"]].rename(columns={"compartment": "ann_compartment"}))
    a.obs["state_final"] = a.obs["state_final"].astype(str)
    a.obs["batch"] = a.obs["batch"].astype(str)
    a.obs["dataset"] = a.obs["dataset"].astype(str)
    print("loaded", a.shape, flush=True)

    sel = cell2location.utils.filtering.filter_genes(
        a, cell_count_cutoff=c2l["gene_cell_count_cutoff"], cell_percentage_cutoff2=c2l["gene_cell_percentage_cutoff2"],
        nonz_mean_cutoff=c2l["gene_nonz_mean_cutoff"])
    a = a[:, sel].copy()
    print("genes kept", a.n_vars, flush=True)

    counts = a.obs.groupby(["state_final", "ann_compartment"], observed=True).size().rename("cells").reset_index()
    counts.to_csv(out / "signature_cells.tsv", sep="\t", index=False)

    RegressionModel.setup_anndata(adata=a, batch_key="batch", labels_key="state_final",
                                  categorical_covariate_keys=["dataset"])
    mod = RegressionModel(a)
    mod.train(max_epochs=max_epochs, batch_size=c2l["signature_batch_size"], lr=c2l["signature_lr"],
              accelerator="auto")
    hist = mod.history["elbo_train"]
    hist.to_csv(out / "training_history.tsv", sep="\t")
    print("elbo first/last", float(hist.iloc[0, 0]), float(hist.iloc[-1, 0]), flush=True)

    a = mod.export_posterior(a, sample_kwargs={"num_samples": c2l["posterior_samples"],
                                               "batch_size": c2l["signature_batch_size"]})
    key = "means_per_cluster_mu_fg"
    cols = [c for c in a.varm[key].columns] if hasattr(a.varm[key], "columns") else None
    sig = pd.DataFrame(np.asarray(a.varm[key]), index=a.var_names,
                       columns=[c.replace(f"{key}_", "") for c in cols] if cols else a.uns["mod"]["factor_names"])
    sig.to_csv(out / "reference_signatures.csv.gz")
    mod.save(str(out / "model"), overwrite=True)
    print("signatures", sig.shape, "written to", out, flush=True)


if __name__ == "__main__":
    main()
