"""Label transfer of a targeted-panel spatial dataset (Xenium) onto the
reference states, on the genes shared between the panel and the reference
(plan, section 4), with the cell2location reference signatures as the
reference.

Each cell's counts over the shared genes are scored against every
candidate state with a multinomial likelihood whose gene probabilities
are the state's reference signature, normalised over those genes; the
posterior over states with a flat prior gives the label and its
probability. A cell below min_prob, or with fewer than min_counts
transcripts over the shared genes, is "unassigned" (the reason is kept).
Candidate states are those of the dataset's compartment (and class, for
the Yu 2024 neurons). The scANVI and scArches transfer tried first is
recorded in the analysis log (3.5).

Check (a): the AUROC of the marker score (markers computed within the
panel genes by reference_markers.py --genes) for the cells assigned to a
state against the rest, per section. Check (b), held-out genes: the
shared genes are split into ten folds (plan amendment 2); for each fold
the cells are classified without the fold's genes and the mean
normalised expression of those genes in the cells assigned to each state
is recorded, so every gene is held out once; the Spearman correlation
across all shared genes between these held-out means and the state's
signature is the check.

Writes to --out-dir: <section>_cells.tsv.gz (state, probability, counts,
reason per cell); to --results-dir: transfer_checks.tsv, transfer_summary.tsv.
Usage: python workflow/scripts/panel_transfer.py <dataset_id> --signatures <reference_signatures.csv.gz> \
           --annotation <reference_annotation.tsv.gz> --markers <panel_markers.tsv> \
           --sections-dir <dir of section h5ads> --out-dir <dir> --results-dir <dir>
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
from scipy.special import logsumexp
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

UNASSIGNED = "unassigned"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("--signatures", required=True)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--sections-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    pt, chk = cfg["panel_transfer"], cfg["spatial_checks"]
    dcfg = pt["datasets"][args.dataset]
    rng = np.random.default_rng(cfg["seed"])
    out, res = pathlib.Path(args.out_dir), pathlib.Path(args.results_dir)
    out.mkdir(parents=True, exist_ok=True)
    res.mkdir(parents=True, exist_ok=True)

    secs = [ad.read_h5ad(p) for p in sorted(glob.glob(str(pathlib.Path(args.sections_dir) / "*.h5ad")))]
    query = ad.concat(secs, join="inner", index_unique="-")
    query.X = sp.csr_matrix(query.X).astype(np.float32)
    query.obs["section"] = query.obs["section"].astype(str)

    # candidate states: those of the compartment (and class) in the annotation
    ann = pd.read_csv(args.annotation, sep="\t", index_col=0)
    ann = ann[(ann["compartment"] == dcfg["compartment"]) & ~ann["state_final"].str.endswith(UNASSIGNED)]
    if dcfg.get("class"):
        ann = ann[ann["class_final"] == dcfg["class"]]
    states = sorted(ann["state_final"].unique())
    sig = pd.read_csv(args.signatures, index_col=0)[states]
    markers = pd.read_csv(args.markers, sep="\t")

    shared = [g for g in sig.index if g in set(query.var_names)]
    if len(shared) < pt["min_shared_genes"]:
        raise SystemExit(f"only {len(shared)} shared genes, below {pt['min_shared_genes']}")
    n_folds = int(round(1 / cfg["cell2location"]["heldout_gene_fraction"]))
    folds = np.array_split(rng.permutation(shared), n_folds)
    print(f"{args.dataset}: {query.n_obs} cells in {len(secs)} sections, {query.n_vars} panel genes, "
          f"{len(shared)} shared in {n_folds} held-out folds; {len(states)} candidate states", flush=True)

    def classify(genes_used):
        # multinomial likelihood of each cell under each state, flat prior
        x = query[:, genes_used].X
        p = sig.loc[genes_used].values + pt["pseudocount"]
        logp = np.log(p / p.sum(axis=0, keepdims=True))
        ll = np.asarray(x @ logp)
        post = np.exp(ll - logsumexp(ll, axis=1, keepdims=True))
        counts = np.asarray(x.sum(axis=1)).ravel()
        pred = np.array(states)[post.argmax(axis=1)]
        prob = post.max(axis=1)
        reason = np.where(counts < pt["min_counts"], "low counts", np.where(prob < pt["min_prob"], "ambiguous", ""))
        return pred, prob, counts, reason, np.where(reason == "", pred, UNASSIGNED)

    pred, prob, counts, reason, state = classify(shared)
    cells = pd.DataFrame({"section": query.obs["section"].values, "pred_state": pred, "prob": prob,
                          "counts_shared": counts, "reason": reason, "state": state}, index=query.obs_names)

    qn = query.copy()
    sc.pp.normalize_total(qn, target_sum=1e4)
    cp10k = qn.X.copy()
    sc.pp.log1p(qn)
    genes = pd.Index(qn.var_names)
    sig_norm = sig.loc[shared] / sig.loc[shared].sum(axis=0)
    # held-out means: gene x state x section, each gene from the fold where it was held out
    sections = sorted(cells["section"].unique())
    heldout = {sec: pd.DataFrame(np.nan, index=shared, columns=states) for sec in sections}
    for fold in folds:
        used = [g for g in shared if g not in set(fold)]
        _, _, _, _, st = classify(used)
        gi = genes.get_indexer(fold)
        for sec in sections:
            m = (cells["section"] == sec).values
            for state in states:
                a = m & (st == state)
                if a.sum() >= chk["min_assigned_spots"]:
                    heldout[sec].loc[fold, state] = np.asarray(cp10k[a][:, gi].mean(axis=0)).ravel()
    rows, summary = [], []
    for sec in sorted(cells["section"].unique()):
        m = (cells["section"] == sec).values
        c = cells[m]
        c.to_csv(out / f"{sec}_cells.tsv.gz", sep="\t")
        n_assigned = int((c["state"] != UNASSIGNED).sum())
        summary.append([sec, len(c), n_assigned, round(1 - n_assigned / len(c), 3),
                        int((c["reason"] == "low counts").sum()), int((c["reason"] == "ambiguous").sum()),
                        c.loc[c["state"] != UNASSIGNED, "state"].nunique(), round(float(c["prob"].median()), 3)])
        for state in states:
            assigned = (c["state"] == state).values
            n = int(assigned.sum())
            mk = [g for g in markers.loc[markers["state"] == state, "gene"] if g in set(genes)]
            auroc, rho, n_held = np.nan, np.nan, 0
            if n >= chk["min_assigned_spots"] and (~assigned).sum() >= chk["min_assigned_spots"]:
                if mk:
                    score = np.asarray(qn[m][:, mk].X.mean(axis=1)).ravel()
                    auroc = roc_auc_score(assigned, score)
                obs = heldout[sec][state]
                ok = obs.notna().values
                ref = sig_norm.loc[shared, state].values
                if ok.sum() >= 10 and np.std(obs[ok]) > 0 and np.std(ref[ok]) > 0:
                    rho = spearmanr(obs[ok], ref[ok]).correlation
                n_held = int(ok.sum())
            rows.append([sec, state, n, len(mk), auroc, n_held, rho,
                         bool(auroc >= chk["auroc_min"]) if not np.isnan(auroc) else "",
                         bool(rho >= chk["spearman_min"]) if not np.isnan(rho) else ""])
    pd.DataFrame(rows, columns=["section", "state", "assigned_cells", "markers_present", "auroc_markers",
                                "heldout_genes", "spearman_heldout", "pass_a", "pass_b"]
                 ).to_csv(res / "transfer_checks.tsv", sep="\t", index=False)
    s = pd.DataFrame(summary, columns=["section", "cells", "assigned", "unassigned_fraction", "low_counts",
                                       "ambiguous", "states", "median_prob"])
    s.to_csv(res / "transfer_summary.tsv", sep="\t", index=False)
    print(s.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
