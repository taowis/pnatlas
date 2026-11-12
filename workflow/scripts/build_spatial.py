"""Read one spatial dataset into one AnnData of raw counts per section and
write them as h5ad with a fixed set of obs columns and, where the authors
give positions, the coordinates in obsm["spatial"].

The per-dataset inputs are in config/spatial_inputs.yaml. Formats:
  xenium         Xenium cell_feature_matrix.h5 (gene probes only) plus the
                 cells table (cell centroids in microns, areas, counts)
  visium_h5      Space Ranger filtered h5, optionally the tissue positions
                 CSV (with or without header) and a Loupe selections CSV
  visium_mtx     a filtered_feature_bc_matrix directory, no positions
  seurat_export  the directory written by rds_spatial_to_mtx.R; split into
                 sections by the image column of coordinates.csv.gz; the
                 exported x,y are read as row,column when position_order
                 is "yx" (Seurat 5 objects)
  neuron_csv     genes x cells CSV with a metadata CSV whose cell_id holds
                 the position as "x-y" (or "y-x": position_order) in pixels
                 (pixel_size_um converts them to microns) and whose section
                 column names the slide

Obs columns written for every section: dataset, section, donor, condition,
batch, technology, tissue, compartment, plus the authors' per-cell columns
prefixed with "author_". uns["spatial_units"] says what the coordinates are.

Usage: python workflow/scripts/build_spatial.py <dataset_id> --data-root "$DATA_ROOT" \
           --config config/spatial_inputs.yaml --out-dir "$DATA_ROOT/interim/spatial/<dataset_id>"
Prints one summary line per section (tab separated).
"""
import argparse
import pathlib
import re

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml

OBS = ["dataset", "section", "donor", "condition", "batch", "technology", "tissue", "compartment"]
POS_COLS = ["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"]


def read_positions(path):
    first = open(path).readline().split(",")[0]
    return pd.read_csv(path, header=0 if first == "barcode" else None, names=POS_COLS, index_col=0)


def read_selections(path):
    sel = pd.read_csv(path, header=None, names=["barcode", "label"], index_col=0)
    return sel["label"].astype(str)


def read_xenium(root, sec):
    a = sc.read_10x_h5(root / sec["h5"], gex_only=True)
    a.var_names_make_unique()
    cells = pd.read_parquet(root / sec["cells"]).set_index("cell_id")
    cells = cells.loc[a.obs_names]
    for c in ["transcript_counts", "cell_area", "nucleus_area"]:
        a.obs["author_" + c] = cells[c].values
    a.obsm["spatial"] = cells[["x_centroid", "y_centroid"]].values.astype(np.float64)
    a.uns["spatial_units"] = "microns"
    return a


def read_visium_h5(root, sec):
    a = sc.read_10x_h5(root / sec["h5"])
    a.var_names_make_unique()
    if sec.get("positions"):
        pos = read_positions(root / sec["positions"]).loc[a.obs_names]
        for c in ["in_tissue", "array_row", "array_col"]:
            a.obs["author_" + c] = pos[c].values
        a.obsm["spatial"] = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].values.astype(np.float64)
        a.uns["spatial_units"] = "pixels of the full-resolution image (x = column, y = row)"
    if sec.get("selections"):
        sel = read_selections(root / sec["selections"])
        a.obs["author_selection"] = sel.reindex(a.obs_names).fillna("").values
    return a


def read_visium_mtx(root, sec):
    a = sc.read_10x_mtx(root / sec["dir"], var_names="gene_symbols", cache=False)
    a.var_names_make_unique()
    return a


def read_seurat_export(root, ds):
    d = root / ds["dir"]
    x = sc.read_mtx(d / "matrix.mtx.gz").X.T.tocsr()
    genes = pd.read_csv(d / "features.tsv.gz", header=None)[0].astype(str).values
    cells = pd.read_csv(d / "barcodes.tsv.gz", header=None)[0].astype(str).values
    a = ad.AnnData(X=x, obs=pd.DataFrame(index=cells), var=pd.DataFrame(index=genes))
    a.var_names_make_unique()
    meta = pd.read_csv(d / "metadata.csv.gz", index_col="barcode", low_memory=False).loc[cells]
    for c in meta.columns:
        a.obs["author_" + c] = meta[c].values
    coords = pd.read_csv(d / "coordinates.csv.gz").drop_duplicates("barcode").set_index("barcode")
    coords = coords.reindex(cells)
    a.obs["author_image"] = coords["image"].fillna("").values
    xy = coords[["x", "y"]].values.astype(np.float64)
    if ds.get("position_order", "xy") == "yx":  # Seurat 5 GetTissueCoordinates: x is the image row
        xy = xy[:, ::-1]
    a.obsm["spatial"] = xy
    a.uns["spatial_units"] = ds.get("spatial_units", "pixels of the full-resolution image (x = column, y = row)")
    out = {}
    for img, sec in ds["sections"].items():
        sub = a[a.obs["author_image"] == img].copy()
        if sub.n_obs == 0:
            raise SystemExit(f"no cells for image {img}")
        sub = sub[~np.isnan(sub.obsm["spatial"][:, 0])].copy()
        out[sec["id"]] = (sub, sec)
    return out


def read_neuron_csv(root, ds):
    df = pd.read_csv(root / ds["matrix"], index_col=0)
    meta = pd.read_csv(root / ds["metadata"], index_col=0).loc[df.columns]
    a = ad.AnnData(X=sp.csr_matrix(df.values.T.astype(np.float32)),
                   obs=pd.DataFrame(index=df.columns.astype(str)), var=pd.DataFrame(index=df.index.astype(str)))
    for c in meta.columns:
        a.obs["author_" + c] = meta[c].values
    # a cell id that is not "x-y" (a few neurons carry a file name instead) gets no position
    xy = meta[ds["position_column"]].astype(str).str.extract(r"^(\d+)-(\d+)$").apply(pd.to_numeric)
    xy = xy.values.astype(np.float64)
    if ds.get("position_order", "xy") == "yx":  # ImageJ names a ROI row-column
        xy = xy[:, ::-1]
    a.obsm["spatial"] = xy * float(ds.get("pixel_size_um", 1.0))
    a.obs["author_has_position"] = ~np.isnan(a.obsm["spatial"][:, 0])
    a.uns["spatial_units"] = ds.get("spatial_units", "image pixels")
    which = meta[ds["section_column"]].astype(str).str.extract(ds["section_pattern"])[0]
    out = {}
    for sid, sec in ds["sections"].items():
        sub = a[(which == sid).values].copy()
        out[sid] = (sub, sec)
    return out


PER_SECTION = {"xenium": read_xenium, "visium_h5": read_visium_h5, "visium_mtx": read_visium_mtx}
PER_DATASET = {"seurat_export": read_seurat_export, "neuron_csv": read_neuron_csv}


def finish(a, dataset, sid, sec, ds):
    a.obs["dataset"] = dataset
    a.obs["section"] = sid
    a.obs["donor"] = str(sec.get("donor", sid))
    a.obs["condition"] = str(sec.get("condition", ds.get("condition", "")))
    a.obs["batch"] = str(sec.get("batch", ds.get("batch", sid)))
    for c in ["technology", "tissue", "compartment"]:
        a.obs[c] = ds[c]
    a.obs = a.obs[OBS + [c for c in a.obs.columns if c not in OBS]]
    if not sp.issparse(a.X):
        a.X = sp.csr_matrix(a.X)
    a.X = sp.csr_matrix(a.X).astype(np.float32)
    counts = np.asarray(a.X.sum(axis=1)).ravel()
    return [dataset, sid, a.obs["donor"].iloc[0], a.obs["condition"].iloc[0], ds["technology"],
            a.n_obs, a.n_vars, int(np.median(counts)), "spatial" in a.obsm]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--config", default="config/spatial_inputs.yaml")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()
    root = pathlib.Path(args.data_root)
    ds = yaml.safe_load(open(args.config))[args.dataset]
    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if ds["format"] in PER_SECTION:
        items = {sid: (PER_SECTION[ds["format"]](root, sec), sec) for sid, sec in ds["sections"].items()}
    else:
        items = PER_DATASET[ds["format"]](root, ds)
    print("\t".join(["dataset", "section", "donor", "condition", "technology", "cells", "genes",
                     "median_counts", "coordinates"]))
    for sid, (a, sec) in items.items():
        row = finish(a, args.dataset, sid, sec, ds)
        a.write_h5ad(out_dir / f"{sid}.h5ad")
        print("\t".join(str(x) for x in row))


if __name__ == "__main__":
    main()
