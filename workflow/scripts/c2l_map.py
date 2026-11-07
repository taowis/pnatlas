"""Map one Visium section onto the reference states with cell2location and
run checks (a) and (b) of the plan (section 4).

Genes shared between the section and the reference signatures are split
once per section (seeded): a held-out tenth is left out of the mapping, the
rest is fitted with the Cell2location model (default priors, the expected
cells per spot and detection alpha from the config). The posterior 5th
percentile of the cell abundance per state is the mapped abundance; the
state with the largest abundance in a spot is the spot's dominant state.

Check (a): for every state with enough dominant spots, the AUROC of the
score of its top reference markers (mean log-normalised expression over
the markers present in the section) for those spots against the rest.
Check (b): for every such state, the Spearman correlation across the
held-out genes between their mean normalised expression in the state's
dominant spots and the reference signature.

Writes to --out-dir: abundance.tsv.gz (spots x states, q05), spots.tsv.gz
(dominant state and its share per spot), checks.tsv (one row per state),
training_history.tsv, and prints one summary line.
Usage: python workflow/scripts/c2l_map.py --section <h5ad> --signatures <csv.gz> \
           --markers <tsv> --out-dir <dir> [--max-epochs N]
"""
import argparse
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--section", required=True)
    ap.add_argument("--signatures", required=True)
    ap.add_argument("--markers", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--max-epochs", type=int, default=None)
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    c2l, chk = cfg["cell2location"], cfg["spatial_checks"]
    seed = cfg["seed"]
    out = pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    import scvi
    from cell2location.models import Cell2location
    scvi.settings.seed = seed

    sec = ad.read_h5ad(args.section)
    sec.var_names_make_unique()
    sec.X = sp.csr_matrix(sec.X).astype(np.float32)
    sig = pd.read_csv(args.signatures, index_col=0)
    markers = pd.read_csv(args.markers, sep="\t")
    shared = [g for g in sig.index if g in set(sec.var_names)]
    rng = np.random.default_rng(seed)
    n_hold = int(round(len(shared) * c2l["heldout_gene_fraction"]))
    held = sorted(rng.choice(shared, n_hold, replace=False).tolist())
    fit_genes = [g for g in shared if g not in set(held)]
    name = sec.obs["section"].iloc[0]
    print(f"{name}: {sec.n_obs} spots, {len(shared)} shared genes, {len(held)} held out", flush=True)

    adata = sec[:, fit_genes].copy()
    keep = np.asarray(adata.X.sum(axis=1)).ravel() > 0
    adata = adata[keep].copy()
    Cell2location.setup_anndata(adata=adata)
    mod = Cell2location(adata, cell_state_df=sig.loc[fit_genes], N_cells_per_location=c2l["n_cells_per_location"],
                        detection_alpha=c2l["detection_alpha"])
    mod.train(max_epochs=args.max_epochs or c2l["map_max_epochs"], batch_size=c2l["map_batch_size"],
              train_size=1, accelerator="auto")
    mod.history["elbo_train"].to_csv(out / "training_history.tsv", sep="\t")
    adata = mod.export_posterior(adata, sample_kwargs={"num_samples": c2l["posterior_samples"],
                                                       "batch_size": adata.n_obs})
    ab = adata.obsm["q05_cell_abundance_w_sf"].copy()
    ab.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in ab.columns]
    ab.to_csv(out / "abundance.tsv.gz", sep="\t")
    share = ab.div(ab.sum(axis=1), axis=0)
    spots = pd.DataFrame({"dominant_state": share.idxmax(axis=1), "dominant_share": share.max(axis=1),
                          "total_abundance": ab.sum(axis=1)}, index=ab.index)
    spots.to_csv(out / "spots.tsv.gz", sep="\t")

    # normalised expression of the whole section for the checks
    norm = sec[adata.obs_names].copy()
    sc.pp.normalize_total(norm, target_sum=1e4)
    cp10k = norm.X.copy()
    sc.pp.log1p(norm)
    genes = pd.Index(norm.var_names)
    rows = []
    for state in sig.columns:
        assigned = (spots["dominant_state"] == state).values
        n = int(assigned.sum())
        mk = [g for g in markers.loc[markers["state"] == state, "gene"] if g in set(genes)]
        auroc, rho, n_held = np.nan, np.nan, 0
        if n >= chk["min_assigned_spots"] and (~assigned).sum() >= chk["min_assigned_spots"]:
            if mk:
                score = np.asarray(norm[:, mk].X.mean(axis=1)).ravel()
                auroc = roc_auc_score(assigned, score)
            hg = [g for g in held if sig.loc[g, state] > 0]
            if len(hg) >= 10:
                obs = np.asarray(cp10k[assigned][:, genes.get_indexer(hg)].mean(axis=0)).ravel()
                rho = spearmanr(obs, sig.loc[hg, state].values).correlation
                n_held = len(hg)
        rows.append([name, state, n, len(mk), auroc, n_held, rho,
                     bool(auroc >= chk["auroc_min"]) if not np.isnan(auroc) else "",
                     bool(rho >= chk["spearman_min"]) if not np.isnan(rho) else ""])
    checks = pd.DataFrame(rows, columns=["section", "state", "dominant_spots", "markers_present", "auroc_markers",
                                         "heldout_genes", "spearman_heldout", "pass_a", "pass_b"])
    checks.to_csv(out / "checks.tsv", sep="\t", index=False)
    done = checks[checks["dominant_spots"] >= chk["min_assigned_spots"]]
    print(f"{name}: {len(done)} states with dominant spots; (a) pass {int((done['pass_a'] == True).sum())}, "
          f"(b) pass {int((done['pass_b'] == True).sum())}; elbo last {float(mod.history['elbo_train'].iloc[-1, 0]):.0f}",
          flush=True)


if __name__ == "__main__":
    main()
