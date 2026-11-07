"""scIB metrics for the integration embeddings on the ganglia and the nerve
subsets, and the method choice by the fixed rule of the analysis plan
(section 3): the candidate with the higher mean biological conservation is
kept if its batch correction score is within the margin of the other
candidate's; otherwise the higher total score wins. Metrics use the mapped
state labels on labelled cells, subsampled to eval_max_cells per subset.

Usage: python workflow/scripts/evaluate_integration.py --in "$DATA_ROOT/interim/reference_embeddings.h5ad" \
           --out-dir results/reference
"""
import argparse
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import yaml

EMBEDDINGS = ["X_scvi", "X_scanvi", "X_harmony"]
CANDIDATES = {"X_scanvi": "scanvi", "X_harmony": "harmony"}


def metrics(adata, cfg, seed, n_jobs):
    from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation

    bm = Benchmarker(
        adata, batch_key=cfg["batch_key"], label_key=cfg["label_key"], embedding_obsm_keys=EMBEDDINGS,
        bio_conservation_metrics=BioConservation(isolated_labels=False, nmi_ari_cluster_labels_kmeans=True,
                                                 nmi_ari_cluster_labels_leiden=False, silhouette_label=True,
                                                 clisi_knn=False),
        # bras is this library's batch silhouette (batch ASW in the plan)
        batch_correction_metrics=BatchCorrection(bras=True, ilisi_knn=True, kbet_per_label=False,
                                                 graph_connectivity=True, pcr_comparison=False),
        n_jobs=n_jobs)
    bm.benchmark()
    return bm.get_results(min_max_scale=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out-dir", default="results/reference")
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--n-jobs", type=int, default=8)
    a = ap.parse_args()
    with open(a.config) as f:
        full = yaml.safe_load(f)
    cfg, seed = full["integration"], full["seed"]
    rng = np.random.default_rng(seed)

    adata = ad.read_h5ad(a.inp)
    st = adata.obs[cfg["label_key"]].astype(str)
    labelled = ~((st == "unlabelled") | st.str.endswith("unassigned"))
    tables = []
    for subset in ["ganglion", "nerve"]:
        idx = np.where(labelled.values & (adata.obs["compartment"].astype(str) == subset).values)[0]
        if len(idx) > cfg["eval_max_cells"]:
            idx = np.sort(rng.choice(idx, cfg["eval_max_cells"], replace=False))
        sub = adata[idx].copy()
        for c in [cfg["batch_key"], cfg["label_key"]]:
            sub.obs[c] = sub.obs[c].astype(str).astype("category")
        res = metrics(sub, cfg, seed, a.n_jobs)
        res = res.loc[EMBEDDINGS].copy()
        res.insert(0, "n_cells", len(idx))
        res.insert(0, "embedding", res.index)
        res.insert(0, "subset", subset)
        tables.append(res)
    out = pathlib.Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    table = pd.concat(tables).reset_index(drop=True)
    table.to_csv(out / "scib_metrics.tsv", sep="\t", index=False)

    # the rule, applied to the mean over the two subsets
    mean = table.groupby("embedding")[["Batch correction", "Bio conservation", "Total"]].mean()
    c1, c2 = "X_scanvi", "X_harmony"
    best_bio = c1 if mean.loc[c1, "Bio conservation"] >= mean.loc[c2, "Bio conservation"] else c2
    other = c2 if best_bio == c1 else c1
    within = mean.loc[best_bio, "Batch correction"] >= mean.loc[other, "Batch correction"] - cfg["rule_margin"]
    if within:
        chosen, why = best_bio, "higher mean bio conservation, batch correction within the margin"
    else:
        chosen = c1 if mean.loc[c1, "Total"] >= mean.loc[c2, "Total"] else c2
        why = "batch correction outside the margin, higher total score"
    lines = [f"chosen\t{CANDIDATES[chosen]}", f"reason\t{why}", f"margin\t{cfg['rule_margin']}", ""]
    lines += ["embedding\tbatch_correction\tbio_conservation\ttotal"]
    lines += [f"{e}\t{r['Batch correction']:.4f}\t{r['Bio conservation']:.4f}\t{r['Total']:.4f}"
              for e, r in mean.iterrows()]
    (out / "integration_choice.tsv").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
