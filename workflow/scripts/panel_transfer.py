"""Label transfer of a targeted-panel spatial dataset (Xenium) onto the
reference states with scANVI, on the genes shared between the panel and
the reference (plan, section 4).

The reference is restricted to the dataset's compartment (and class, for
the Yu 2024 neurons) and to the shared genes minus a held-out tenth; scVI
and scANVI are trained on it (donor key as batch, study as covariate,
state_final as label), the query cells of every section are mapped with
scArches (query batch = section) and the scANVI soft prediction gives the
state and its probability; below min_prob the cell is "unassigned".

Check (a): markers are recomputed within the shared genes (top n per
state by the same one-vs-rest test on a reference subsample), and the
AUROC of the marker score for the cells assigned to a state against the
rest is reported per section. Check (b): the Spearman correlation across
the held-out genes between their mean normalised expression in the
assigned cells and in the reference cells of that state.

Writes to --out-dir: <section>_cells.tsv.gz (state, probability per cell),
training_history.tsv; to --results-dir: panel_markers.tsv,
transfer_checks.tsv, transfer_summary.tsv.
Usage: python workflow/scripts/panel_transfer.py <dataset_id> --counts <reference_merged.h5ad> \
           --annotation <reference_annotation.tsv.gz> --sections-dir <dir of section h5ads> \
           --out-dir <dir> --results-dir <dir>
"""
import argparse
import glob
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

from pnatlas.io import read_counts_float32
from pnatlas.markers import subsample_per_label, top_markers

UNASSIGNED = "unassigned"


def mean_cp10k(adata, genes):
    """Mean counts-per-10k expression of genes (normalised over all genes of adata)."""
    a = adata.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    return np.asarray(a[:, genes].X.mean(axis=0)).ravel()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("--counts", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--sections-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    pt, chk = cfg["panel_transfer"], cfg["spatial_checks"]
    dcfg = pt["datasets"][args.dataset]
    seed = cfg["seed"]
    rng = np.random.default_rng(seed)
    out, res = pathlib.Path(args.out_dir), pathlib.Path(args.results_dir)
    out.mkdir(parents=True, exist_ok=True)
    res.mkdir(parents=True, exist_ok=True)

    import scvi
    scvi.settings.seed = seed

    # query sections
    secs = [ad.read_h5ad(p) for p in sorted(glob.glob(str(pathlib.Path(args.sections_dir) / "*.h5ad")))]
    query = ad.concat(secs, join="inner", index_unique="-")
    query.X = sp.csr_matrix(query.X).astype(np.float32)
    query.obs["section"] = query.obs["section"].astype(str)
    panel = list(query.var_names)

    # reference cells of the compartment (and class)
    ann = pd.read_csv(args.annotation, sep="\t", index_col=0)
    ann = ann[(ann["compartment"] == dcfg["compartment"]) & ~ann["state_final"].str.endswith(UNASSIGNED)]
    if dcfg.get("class"):
        ann = ann[ann["class_final"] == dcfg["class"]]
    ref = read_counts_float32(args.counts, ann.index)
    ref.obs["state_final"] = ann.loc[ref.obs_names, "state_final"].astype(str).values
    shared = [g for g in ref.var_names if g in set(panel)]
    if len(shared) < pt["min_shared_genes"]:
        raise SystemExit(f"only {len(shared)} shared genes, below {pt['min_shared_genes']}")
    held = sorted(rng.choice(shared, int(round(len(shared) * cfg["cell2location"]["heldout_gene_fraction"])),
                             replace=False).tolist())
    fit_genes = [g for g in shared if g not in set(held)]
    print(f"{args.dataset}: {query.n_obs} query cells in {len(secs)} sections, {len(panel)} panel genes, "
          f"{len(shared)} shared, {len(held)} held out; reference {ref.n_obs} cells, "
          f"{ref.obs['state_final'].nunique()} states", flush=True)

    # markers within the shared genes, and reference means for check (b)
    ref_shared = ref[:, shared].copy()
    pick = subsample_per_label(ref_shared.obs["state_final"], chk["markers_cells_per_state"], rng)
    markers = top_markers(ref_shared[pick].copy(), "state_final", chk["markers_padj"], chk["markers_min_pct"],
                          chk["markers_top_n"])
    markers.to_csv(res / "panel_markers.tsv", sep="\t", index=False)
    ref_means = {}
    for state in sorted(ref_shared.obs["state_final"].unique()):
        ref_means[state] = mean_cp10k(ref_shared[ref_shared.obs["state_final"] == state], held)

    # scVI and scANVI on the reference, scArches on the query
    r = ref[:, fit_genes].copy()
    r.obs["batch"] = r.obs["batch"].astype(str)
    r.obs["dataset"] = r.obs["dataset"].astype(str)
    scvi.model.SCVI.setup_anndata(r, batch_key="batch", categorical_covariate_keys=["dataset"], labels_key="state_final")
    vae = scvi.model.SCVI(r, n_latent=pt["n_latent"], n_layers=pt["n_layers"], n_hidden=pt["n_hidden"],
                          gene_likelihood="nb")
    vae.train(max_epochs=pt["scvi_max_epochs"], batch_size=pt["batch_size"], early_stopping=True)
    lvae = scvi.model.SCANVI.from_scvi_model(vae, unlabeled_category="unknown", labels_key="state_final")
    lvae.train(max_epochs=pt["scanvi_max_epochs"], n_samples_per_label=pt["n_samples_per_label"],
               batch_size=pt["batch_size"])
    hist = pd.concat({"scvi": vae.history["elbo_train"], "scanvi": lvae.history["elbo_train"]}, axis=1)
    hist.to_csv(out / "training_history.tsv", sep="\t")

    q = query[:, fit_genes].copy()
    q.obs["batch"] = ("query:" + q.obs["section"]).values
    q.obs["dataset"] = args.dataset
    q.obs["state_final"] = "unknown"
    scvi.model.SCANVI.prepare_query_anndata(q, lvae)
    qm = scvi.model.SCANVI.load_query_data(q, lvae)
    qm.train(max_epochs=pt["query_max_epochs"], plan_kwargs={"weight_decay": 0.0}, batch_size=pt["batch_size"],
             check_val_every_n_epoch=None)
    probs = qm.predict(soft=True)
    pred = probs.idxmax(axis=1)
    pmax = probs.max(axis=1)
    cells = pd.DataFrame({"section": query.obs["section"].values, "pred_state": pred.values, "prob": pmax.values},
                         index=query.obs_names)
    cells["state"] = np.where(cells["prob"] >= pt["min_prob"], cells["pred_state"], UNASSIGNED)

    # normalised query expression for the checks
    qn = query.copy()
    sc.pp.normalize_total(qn, target_sum=1e4)
    cp10k = qn.X.copy()
    sc.pp.log1p(qn)
    genes = pd.Index(qn.var_names)
    rows, summary = [], []
    for sec in sorted(cells["section"].unique()):
        m = (cells["section"] == sec).values
        c = cells[m]
        c.to_csv(out / f"{sec}_cells.tsv.gz", sep="\t")
        n_assigned = int((c["state"] != UNASSIGNED).sum())
        summary.append([sec, len(c), n_assigned, round(1 - n_assigned / len(c), 3), c["state"].nunique() - 1,
                        round(float(c["prob"].median()), 3)])
        for state in sorted(ref_means):
            assigned = (c["state"] == state).values
            n = int(assigned.sum())
            mk = [g for g in markers.loc[markers["state"] == state, "gene"]]
            auroc, rho = np.nan, np.nan
            if n >= chk["min_assigned_spots"] and (~assigned).sum() >= chk["min_assigned_spots"]:
                if mk:
                    score = np.asarray(qn[m][:, mk].X.mean(axis=1)).ravel()
                    auroc = roc_auc_score(assigned, score)
                obs = np.asarray(cp10k[m][assigned][:, genes.get_indexer(held)].mean(axis=0)).ravel()
                if np.std(obs) > 0 and np.std(ref_means[state]) > 0:
                    rho = spearmanr(obs, ref_means[state]).correlation
            rows.append([sec, state, n, len(mk), auroc, len(held), rho,
                         bool(auroc >= chk["auroc_min"]) if not np.isnan(auroc) else "",
                         bool(rho >= chk["spearman_min"]) if not np.isnan(rho) else ""])
    pd.DataFrame(rows, columns=["section", "state", "assigned_cells", "markers_present", "auroc_markers",
                                "heldout_genes", "spearman_heldout", "pass_a", "pass_b"]
                 ).to_csv(res / "transfer_checks.tsv", sep="\t", index=False)
    s = pd.DataFrame(summary, columns=["section", "cells", "assigned", "unassigned_fraction", "states",
                                       "median_prob"])
    s.to_csv(res / "transfer_summary.tsv", sep="\t", index=False)
    print(s.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
