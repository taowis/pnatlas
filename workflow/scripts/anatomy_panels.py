"""Panels for the anatomy check (plan, section 4, check (c)): for every
section of a dataset that has an image and coordinates, the tissue image
alone and the same image with the spots or cells of the perineurial,
epineurial fibroblast and endoneurial fibroblast states drawn on it, so
that a reader can say whether those states lie where the image shows
the perineurium, the epineurium and the endoneurium. Other assigned
spots or cells are drawn in grey; unassigned ones are not drawn.

The image source per dataset is the `anatomy` block of the dataset in
config/spatial_inputs.yaml: a Space Ranger tissue image with its scale
factor (Visium from H5), the low-resolution image exported from a Seurat
object with the scale factor from the exported coordinates (Visium from
a Seurat export), or the Xenium morphology OME-TIFF read at a pyramid
level with the pixel size from its metadata (cells in microns). The
labels are the dominant state per spot from the cell2location run or
the state per cell from the panel transfer.

Writes <out-dir>/<section>.png per section and appends the sections to
<out-dir>/../anatomy_check.tsv, the sheet the reader fills in (one row
per section and state: yes, no or not present).
Usage: python workflow/scripts/anatomy_panels.py <dataset> --data-root data --out-dir results/spatial/anatomy/<dataset>
"""
import argparse
import glob
import json
import pathlib
import re

import anndata as ad
import matplotlib
import numpy as np
import pandas as pd
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt

STATES = {"perineurial": "#d62728", "epineurial fibroblast": "#1f77b4", "endoneurial fibroblast": "#2ca02c"}
UNASSIGNED = "unassigned"


def read_image(spec, root, section, image_name):
    """Return (rgb array, scale from the coordinate units to image pixels)."""
    kind = spec["kind"]
    if kind == "spaceranger":
        png = glob.glob(str(root / spec["image"].format(section=section)))
        js = glob.glob(str(root / spec["scale_json"].format(section=section)))
        if not png or not js:
            return None, None
        with open(js[0]) as f:
            scale = json.load(f)[spec["scale_key"]]
        return plt.imread(png[0]), scale
    if kind == "seurat_export":
        png = root / spec["image"].format(image=image_name)
        if not png.exists():
            return None, None
        coords = pd.read_csv(root / spec["coordinates"])
        scale = float(coords.loc[coords["image"] == image_name, spec["scale_column"]].iloc[0])
        return plt.imread(png), scale
    if kind == "ome_tiff":
        import tifffile
        paths = glob.glob(str(root / spec["image"].format(section=section)))
        if not paths:
            return None, None
        with tifffile.TiffFile(paths[0]) as t:
            series = t.series[0]
            level = min(spec.get("level", 3), len(series.levels) - 1)
            arr = series.levels[level].asarray()
            down = series.levels[0].shape[-1] / series.levels[level].shape[-1]
            m = re.search(r'PhysicalSizeX="([0-9.]+)"', t.ome_metadata or "")
            px = float(m.group(1)) if m else spec.get("pixel_size_um", 0.2125)
        arr = np.squeeze(arr)
        if arr.ndim == 3:  # channels first: take the first (DAPI)
            arr = arr[0]
        lo, hi = np.percentile(arr, (1, 99.8))
        gray = np.clip((arr.astype(np.float32) - lo) / max(hi - lo, 1), 0, 1)
        return np.dstack([gray] * 3), 1.0 / (px * down)
    raise SystemExit(f"unknown image kind {kind}")


def read_labels(root, dataset, section, obs_names):
    c2l = root / "interim/spatial/c2l" / dataset / section / "spots.tsv.gz"
    tr = root / "interim/spatial/transfer" / dataset / f"{section}_cells.tsv.gz"
    if c2l.exists():
        lab = pd.read_csv(c2l, sep="\t", index_col=0)["dominant_state"]
    elif tr.exists():
        lab = pd.read_csv(tr, sep="\t", index_col=0)["state"]
    else:
        return None
    if not lab.index.isin(obs_names).all():
        lab.index = lab.index.str.replace(r"-\d+$", "", regex=True)
    return lab.reindex(obs_names)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--config", default="config/spatial_inputs.yaml")
    ap.add_argument("--width-px", type=int, default=1400, help="width of each panel in pixels")
    args = ap.parse_args()
    root, out = pathlib.Path(args.data_root), pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with open(args.config) as f:
        spec = yaml.safe_load(f)[args.dataset].get("anatomy")
    if spec is None:
        raise SystemExit(f"{args.dataset}: no anatomy block in the config")
    sheet_path = out.parent / "anatomy_check.tsv"
    sheet = pd.read_csv(sheet_path, sep="\t") if sheet_path.exists() else \
        pd.DataFrame(columns=["dataset", "section", "state", "spots_or_cells", "placed_correctly", "note"])
    rows = []
    for p in sorted(glob.glob(str(root / "interim/spatial" / args.dataset / "*.h5ad"))):
        a = ad.read_h5ad(p)
        section = pathlib.Path(p).stem
        image_name = str(a.obs["author_image"].iloc[0]) if "author_image" in a.obs else section
        img, scale = read_image(spec, root, section, image_name)
        if img is None:
            print(f"{section}: no image", flush=True)
            continue
        lab = read_labels(root, args.dataset, section, a.obs_names)
        if lab is None:
            print(f"{section}: no labels", flush=True)
            continue
        xy = a.obsm["spatial"] * scale
        ok = ~np.isnan(xy[:, 0]) & lab.notna().values & (lab.values != UNASSIGNED)
        h, w = img.shape[:2]
        fig_w = 2 * args.width_px / 100
        fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_w / 2 * h / w), dpi=100)
        for ax in axes:
            ax.imshow(img)
            ax.set_xlim(0, w)
            ax.set_ylim(h, 0)
            ax.axis("off")
        size = (spec.get("marker_px", 6) * args.width_px / w) ** 2
        other = ok & ~lab.isin(list(STATES)).values
        axes[1].scatter(xy[other, 0], xy[other, 1], s=size, c="#bbbbbb", alpha=0.5, linewidths=0)
        for state, colour in STATES.items():
            m = ok & (lab.values == state)
            axes[1].scatter(xy[m, 0], xy[m, 1], s=size, c=colour, alpha=0.9, linewidths=0, label=f"{state} ({m.sum()})")
            rows.append([args.dataset, section, state, int(m.sum()), "" if m.sum() else "not present", ""])
        axes[1].legend(loc="lower right", fontsize=8, markerscale=2, framealpha=0.8)
        axes[0].set_title(f"{args.dataset} {section}", fontsize=10, loc="left")
        axes[1].set_title(f"{int(ok.sum())} assigned of {a.n_obs}; grey = other states", fontsize=10, loc="left")
        fig.tight_layout(pad=0.2)
        fig.savefig(out / f"{section}.png")
        plt.close(fig)
        print(f"{section}: {int(ok.sum())} drawn, image {w}x{h}", flush=True)
    new = pd.DataFrame(rows, columns=sheet.columns)
    sheet = sheet[~sheet["section"].isin(new["section"]) | (sheet["dataset"] != args.dataset)]
    pd.concat([sheet, new]).sort_values(["dataset", "section", "state"]).to_csv(sheet_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
