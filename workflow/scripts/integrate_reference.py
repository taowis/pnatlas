"""Integrate the merged reference: scVI (donor-level batch, study covariate),
scANVI seeded with the mapped author labels, and Harmony on the PCA of the
same highly variable genes. Parameters come from the integration block of
config/config.yaml. Writes an h5ad with obs and the embeddings (no counts;
the merged file keeps those), the models, and a training summary.

Usage: python workflow/scripts/integrate_reference.py --in "$DATA_ROOT/interim/reference_merged.h5ad" \
           --out "$DATA_ROOT/interim/reference_embeddings.h5ad" --model-dir "$DATA_ROOT/interim/models" \
           --summary results/reference/integration_training.tsv
"""
import argparse
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import yaml


def harmony(sub, cfg, seed):
    import harmonypy as hm

    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)
    sc.pp.scale(sub, max_value=10)
    sc.tl.pca(sub, n_comps=cfg["n_pcs"], svd_solver="arpack", random_state=seed)
    ho = hm.run_harmony(sub.obsm["X_pca"], sub.obs, cfg["harmony_vars"], random_state=seed, verbose=False)
    z = np.asarray(ho.Z_corr)
    if z.shape[0] != sub.n_obs:
        z = z.T
    return sub.obsm["X_pca"].copy(), z


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--model-dir", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    a = ap.parse_args()
    with open(a.config) as f:
        full = yaml.safe_load(f)
    cfg, seed = full["integration"], full["seed"]

    import scvi
    scvi.settings.seed = seed
    adata = ad.read_h5ad(a.inp)
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=cfg["n_hvg"],
                                batch_key=cfg["hvg_batch_key"], subset=False)
    hvg = adata.var_names[adata.var["highly_variable"]]
    sub = adata[:, hvg].copy()
    del adata

    # scANVI seed labels: mapped states, with unassigned and unlabelled cells unknown
    st = sub.obs[cfg["label_key"]].astype(str)
    unknown = (st == "unlabelled") | st.str.endswith("unassigned")
    sub.obs["scanvi_label"] = st.where(~unknown, "unknown").astype("category")

    mdir = pathlib.Path(a.model_dir)
    scvi.model.SCVI.setup_anndata(sub, batch_key=cfg["batch_key"], categorical_covariate_keys=[cfg["study_key"]])
    m = scvi.model.SCVI(sub, n_latent=cfg["n_latent"], n_layers=cfg["n_layers"], n_hidden=cfg["n_hidden"],
                        gene_likelihood=cfg["gene_likelihood"])
    m.train(max_epochs=cfg["scvi_max_epochs"], batch_size=cfg["batch_size"], early_stopping=True)
    m.save(mdir / "scvi", overwrite=True)
    sub.obsm["X_scvi"] = m.get_latent_representation()
    rows = [{"model": "scvi", "epochs": len(m.history["elbo_train"]),
             "elbo_train": float(m.history["elbo_train"].iloc[-1, 0]),
             "elbo_validation": float(m.history["elbo_validation"].iloc[-1, 0])}]

    s = scvi.model.SCANVI.from_scvi_model(m, unlabeled_category="unknown", labels_key="scanvi_label")
    s.train(max_epochs=cfg["scanvi_max_epochs"], n_samples_per_label=cfg["scanvi_n_samples_per_label"],
            batch_size=cfg["batch_size"], early_stopping=True)
    s.save(mdir / "scanvi", overwrite=True)
    sub.obsm["X_scanvi"] = s.get_latent_representation()
    sub.obs["scanvi_pred"] = s.predict()
    rows.append({"model": "scanvi", "epochs": len(s.history["elbo_train"]),
                 "elbo_train": float(s.history["elbo_train"].iloc[-1, 0]),
                 "elbo_validation": float(s.history["elbo_validation"].iloc[-1, 0])})

    sub.obsm["X_pca"], sub.obsm["X_harmony"] = harmony(sub, cfg, seed)
    rows.append({"model": "harmony", "epochs": np.nan, "elbo_train": np.nan, "elbo_validation": np.nan})

    out = ad.AnnData(obs=sub.obs, obsm=dict(sub.obsm))
    out.uns["integration"] = {k: v for k, v in cfg.items() if k != "datasets"}
    out.uns["integration"]["hvg"] = list(hvg)
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.write_h5ad(a.out)
    pathlib.Path(a.summary).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(a.summary, sep="\t", index=False)
    print(f"integrated\t{out.n_obs} cells\t{len(hvg)} hvg\t" + "\t".join(f"{r['model']}:{r['epochs']}" for r in rows))


if __name__ == "__main__":
    main()
