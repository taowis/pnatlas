"""Read one reference-layer dataset into an AnnData of raw counts and write it
as h5ad with a fixed set of obs columns.

The per-dataset inputs (format, files, which columns give sample, donor and
the authors' labels) are in config/reference_inputs.yaml. Formats:
  csv_genes_x_cells   one gzipped CSV/TSV per sample, genes in rows
  tenx_h5             one Cell Ranger h5 per sample
  mtx_dir             a 10x-style directory (matrix.mtx.gz, features, barcodes),
                      as written by Cell Ranger or by rds_to_mtx.R
An optional metadata table (CSV) is joined to the cells by barcode.

Obs columns written for every dataset: dataset, sample, donor, author_label,
author_label_fine (empty string when the authors gave none), plus every
column of the authors' metadata prefixed with "author_".

Usage: python workflow/scripts/build_anndata.py <dataset_id> --data-root "$DATA_ROOT" \
           --config config/reference_inputs.yaml --out "$DATA_ROOT/interim/reference/<dataset_id>.h5ad"
"""
import argparse
import gzip
import pathlib
import re

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml


def read_csv_genes_x_cells(path):
    sep = "\t" if path.name.endswith((".txt.gz", ".tsv.gz", ".txt", ".tsv")) else ","
    df = pd.read_csv(path, sep=sep, index_col=0)
    x = sp.csr_matrix(df.values.T.astype(np.float32))
    return ad.AnnData(X=x, obs=pd.DataFrame(index=df.columns.astype(str)),
                      var=pd.DataFrame(index=df.index.astype(str)))


def read_tenx_h5(path):
    a = sc.read_10x_h5(path)
    a.var_names_make_unique()
    return a


def read_mtx_dir(path):
    a = sc.read_10x_mtx(path, var_names="gene_symbols", cache=False) if (path / "features.tsv.gz").exists() \
        and _is_tenx_features(path / "features.tsv.gz") else _read_plain_mtx(path)
    a.var_names_make_unique()
    return a


def _is_tenx_features(p):
    with gzip.open(p, "rt") as f:
        return len(f.readline().rstrip("\n").split("\t")) >= 2


def _read_plain_mtx(path):
    x = sc.read_mtx(path / "matrix.mtx.gz").X.T.tocsr()
    genes = pd.read_csv(path / "features.tsv.gz", header=None, sep="\t")[0].astype(str)
    cells = pd.read_csv(path / "barcodes.tsv.gz", header=None, sep="\t")[0].astype(str)
    return ad.AnnData(X=x, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))


READERS = {"csv_genes_x_cells": read_csv_genes_x_cells, "tenx_h5": read_tenx_h5, "mtx_dir": read_mtx_dir}


def sample_from_name(name, pattern):
    m = re.search(pattern, name)
    return m.group(1) if m else name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--config", default="config/reference_inputs.yaml")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    with open(a.config) as f:
        cfg = yaml.safe_load(f)[a.dataset]
    root = pathlib.Path(a.data_root)
    reader = READERS[cfg["format"]]
    paths = []
    for f in cfg["files"]:
        paths += sorted(root.glob(f)) if any(ch in f for ch in "*?[") else [root / f]
    parts = []
    for p in paths:
        x = reader(p)
        sample = sample_from_name(p.name if cfg["format"] != "mtx_dir" else p.parent.name,
                                  cfg.get("sample_pattern", r"^(.*)$")) \
            if cfg.get("sample_from") == "file" else None
        if sample is not None:
            x.obs["sample"] = sample
            if cfg.get("prefix_barcodes", True):
                x.obs_names = [f"{sample}_{b}" for b in x.obs_names]
        parts.append(x)
    adata = ad.concat(parts, join="outer", fill_value=0) if len(parts) > 1 else parts[0]
    adata.obs_names_make_unique()
    if "feature_types" in adata.var:
        # Cell Ranger ARC matrices carry ATAC peaks as well; keep the RNA features
        adata = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()

    meta = cfg.get("metadata")
    if meta:
        m = pd.read_csv(root / meta["file"])
        key = meta.get("barcode_column", m.columns[0])
        m[key] = m[key].astype(str)
        if meta.get("barcode_transform"):
            m[key] = m[key].str.replace(*meta["barcode_transform"], regex=True)
        m = m.set_index(key)
        if meta.get("keep_only_annotated", False):
            adata = adata[adata.obs_names.isin(m.index)].copy()
        m = m.reindex(adata.obs_names)
        m.columns = ["author_" + c for c in m.columns]
        for c in m.columns:
            adata.obs[c] = m[c].values
        for ours, theirs in meta.get("columns", {}).items():
            adata.obs[ours] = m["author_" + theirs].astype(str).where(m["author_" + theirs].notna(), "")
    if "sample" not in adata.obs:
        adata.obs["sample"] = adata.obs_names.str.extract(cfg["sample_pattern"])[0].values
    if "donor" not in adata.obs:
        if cfg.get("donor_pattern"):
            adata.obs["donor"] = adata.obs["sample"].astype(str).str.extract(cfg["donor_pattern"])[0].values
        elif cfg.get("donor_map"):
            adata.obs["donor"] = adata.obs["sample"].map(cfg["donor_map"]).astype(str)
        else:
            adata.obs["donor"] = adata.obs["sample"].astype(str)
    for c in ("author_label", "author_label_fine"):
        if c not in adata.obs:
            adata.obs[c] = ""
    adata.obs["dataset"] = a.dataset
    adata.obs = adata.obs.astype({c: str for c in ("dataset", "sample", "donor", "author_label", "author_label_fine")})
    adata.var["gene_symbol"] = adata.var_names
    front = ["dataset", "sample", "donor", "author_label", "author_label_fine"]
    adata.obs = adata.obs[front + [c for c in adata.obs.columns if c not in front]]

    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(out, compression="gzip")
    print(f"{a.dataset}\tcells\t{adata.n_obs}\tgenes\t{adata.n_vars}\tsamples\t{adata.obs['sample'].nunique()}"
          f"\tdonors\t{adata.obs['donor'].nunique()}\tlabelled\t{(adata.obs['author_label'] != '').sum()}")


if __name__ == "__main__":
    main()
