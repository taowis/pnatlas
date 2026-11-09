"""Panels for the anatomy check (plan, section 4, check (c)): for every
section of a dataset that has an image and coordinates, the tissue image
and, next to it, where the perineurial, epineurial fibroblast and
endoneurial fibroblast states were placed, so that a reader can say
whether they lie where the image shows the perineurium, the epineurium
and the endoneurium.

For a Visium section mapped with cell2location the panel of each state
colours the spots by the estimated number of cells of that state in the
spot (the q05 abundance), spots with fewer than a third of a cell left
blank, so the reader sees the state's abundance over the section and not
only the spots it dominates (a nerve state rarely dominates a 30 cell
spot of a ganglion). For a targeted panel section the cells assigned to
the three states are drawn in colour and the other assigned cells in
grey.

The image source per dataset is the `anatomy` block of the dataset in
config/spatial_inputs.yaml: a Space Ranger tissue image with its scale
factor (Visium from H5), the low-resolution image exported from a Seurat
object with the scale factor from the exported coordinates (Visium from
a Seurat export), the Xenium morphology OME-TIFF read at a pyramid
level with the pixel size from its metadata (cells in microns), or an
H and E image aligned to the Xenium run in Xenium Explorer, whose
alignment matrix (3 by 3, H and E pixels to morphology pixels) is
inverted to place the cells (microns over the morphology pixel size)
on the H and E. When a section has several scale factor files (two
Space Ranger runs), the one whose scaled spot coordinates fit the image
is used.

Writes <out-dir>/<section>.jpg per section and the reader's sheet
<out-dir>/../anatomy_check.tsv (one row per section and state with the
count drawn and an empty `placed_correctly` column: yes, no or not present).
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

STATES = {"perineurial": ("#d62728", "Reds"), "epineurial fibroblast": ("#1f77b4", "Blues"),
          "endoneurial fibroblast": ("#2ca02c", "Greens")}
UNASSIGNED = "unassigned"
MIN_CELLS = 1 / 3  # estimated cells of a state below which a spot is left blank


def read_image(spec, root, section, image_name, xy):
    """Return (rgb array, scale from the coordinate units to image pixels);
    the scale is a number or, for an aligned image, a function of the coordinates."""
    kind = spec["kind"]
    if kind == "spaceranger":
        png = glob.glob(str(root / spec["image"].format(section=section)))
        jss = sorted(glob.glob(str(root / spec["scale_json"].format(section=section))))
        if not png or not jss:
            return None, None
        img = plt.imread(png[0])
        scales = []
        for js in jss:
            with open(js) as f:
                scales.append(json.load(f)[spec["scale_key"]])
        h, w = img.shape[:2]
        fit = [s for s in scales if (np.nanmax(xy[:, 0]) * s <= w) and (np.nanmax(xy[:, 1]) * s <= h)]
        return img, max(fit) if fit else min(scales)
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
        if arr.ndim == 3:  # a DAPI z-stack (Xenium morphology.ome.tif): maximum projection
            arr = arr.max(axis=0)
        lo, hi = np.percentile(arr[arr > 0], (1, 99.8))
        gray = np.clip((arr.astype(np.float32) - lo) / max(hi - lo, 1), 0, 1)
        return np.dstack([gray] * 3), 1.0 / (px * down)
    if kind == "he_aligned":
        import tifffile
        paths = glob.glob(str(root / spec["image"].format(section=section)))
        mats = glob.glob(str(root / spec["matrix"].format(section=section)))
        if not paths or not mats:
            return None, None
        with tifffile.TiffFile(paths[0]) as t:
            series = t.series[0]
            level = min(spec.get("level", 1), len(series.levels) - 1)
            arr = np.squeeze(series.levels[level].asarray())
            down = series.levels[0].shape[1] / series.levels[level].shape[1]
        # the matrix maps H and E pixels to morphology pixels; cells are in microns
        inv = np.linalg.inv(np.loadtxt(mats[0], delimiter=","))
        px = float(spec.get("pixel_size_um", 0.2125))

        def to_image(xy_um):
            pts = np.c_[xy_um / px, np.ones(len(xy_um))] @ inv.T
            return pts[:, :2] / down

        return arr, to_image
    raise SystemExit(f"unknown image kind {kind}")


def align(df, obs_names):
    if not df.index.isin(obs_names).all():
        df.index = df.index.str.replace(r"-\d+$", "", regex=True)
    return df.reindex(obs_names)


def read_mapping(root, dataset, section, obs_names):
    """(abundance table, None) for a cell2location section, (None, labels) for a transfer."""
    c2l = root / "interim/spatial/c2l" / dataset / section / "abundance.tsv.gz"
    tr = root / "interim/spatial/transfer" / dataset / f"{section}_cells.tsv.gz"
    if c2l.exists():
        return align(pd.read_csv(c2l, sep="\t", index_col=0), obs_names), None
    if tr.exists():
        return None, align(pd.read_csv(tr, sep="\t", index_col=0), obs_names)["state"]
    return None, None


def draw_image(ax, img):
    h, w = img.shape[:2]
    ax.imshow(img)
    ax.set_xlim(0, w)
    ax.set_ylim(h, 0)
    ax.axis("off")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dataset")
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--config", default="config/spatial_inputs.yaml")
    ap.add_argument("--width-px", type=int, default=1100, help="width of each panel in pixels")
    args = ap.parse_args()
    root, out = pathlib.Path(args.data_root), pathlib.Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    with open(args.config) as f:
        spec = yaml.safe_load(f)[args.dataset].get("anatomy")
    if spec is None:
        raise SystemExit(f"{args.dataset}: no anatomy block in the config")
    sheet_path = out.parent / "anatomy_check.tsv"
    cols = ["dataset", "section", "state", "spots_or_cells", "placed_correctly", "note"]
    sheet = pd.read_csv(sheet_path, sep="\t", dtype=str).fillna("") if sheet_path.exists() else pd.DataFrame(columns=cols)
    rows = []
    for p in sorted(glob.glob(str(root / "interim/spatial" / args.dataset / "*.h5ad"))):
        a = ad.read_h5ad(p)
        section = pathlib.Path(p).stem
        image_name = str(a.obs["author_image"].iloc[0]) if "author_image" in a.obs else section
        img, scale = read_image(spec, root, section, image_name, a.obsm["spatial"])
        if img is None:
            print(f"{section}: no image", flush=True)
            continue
        abundance, labels = read_mapping(root, args.dataset, section, a.obs_names)
        if abundance is None and labels is None:
            print(f"{section}: no mapping", flush=True)
            continue
        xy = scale(a.obsm["spatial"]) if callable(scale) else a.obsm["spatial"] * scale
        has_xy = ~np.isnan(xy[:, 0])
        h, w = img.shape[:2]
        n_panels = 4 if abundance is not None else 2
        fig_w = n_panels * args.width_px / 100
        fig, axes = plt.subplots(1, n_panels, figsize=(fig_w, fig_w / n_panels * h / w + 0.4), dpi=100)
        for ax in axes:
            draw_image(ax, img)
        size = (spec.get("marker_px", 6) * args.width_px / w) ** 2
        axes[0].set_title(f"{args.dataset} {section}", fontsize=10, loc="left")
        if abundance is not None:
            for ax, (state, (_, cmap)) in zip(axes[1:], STATES.items()):
                v = abundance[state].values
                m = has_xy & (v >= MIN_CELLS)
                vmax = max(np.percentile(v[m], 98), 1.0) if m.sum() else 1.0
                sc = ax.scatter(xy[m, 0], xy[m, 1], s=size, c=v[m], cmap=cmap, vmin=0, vmax=vmax, alpha=0.9, linewidths=0)
                fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.01).set_label("estimated cells in spot", fontsize=8)
                ax.set_title(f"{state}: {int(m.sum())} spots with {MIN_CELLS:.2f}+ cells, "
                             f"{v.sum() / abundance.values.sum() * 100:.1f}% of all cells", fontsize=9, loc="left")
                rows.append([args.dataset, section, state, int(m.sum()), "" if m.sum() else "not present", ""])
        else:
            ok = has_xy & labels.notna().values & (labels.values != UNASSIGNED)
            other = ok & ~labels.isin(list(STATES)).values
            axes[1].scatter(xy[other, 0], xy[other, 1], s=size, c="#bbbbbb", alpha=0.5, linewidths=0)
            for state, (colour, _) in STATES.items():
                m = ok & (labels.values == state)
                axes[1].scatter(xy[m, 0], xy[m, 1], s=size, c=colour, alpha=0.9, linewidths=0, label=f"{state} ({m.sum()})")
                rows.append([args.dataset, section, state, int(m.sum()), "" if m.sum() else "not present", ""])
            axes[1].legend(loc="lower right", fontsize=8, markerscale=2, framealpha=0.8)
            axes[1].set_title(f"{int(ok.sum())} assigned of {a.n_obs} cells; grey = other states", fontsize=9, loc="left")
        fig.tight_layout(pad=0.2)
        fig.savefig(out / f"{section}.jpg", pil_kwargs={"quality": 88})
        plt.close(fig)
        print(f"{section}: image {w}x{h}, scale {'affine' if callable(scale) else format(scale, '.4g')}", flush=True)
    new = pd.DataFrame(rows, columns=cols)
    keep = ~((sheet["dataset"] == args.dataset) & sheet["section"].isin(new["section"]))
    pd.concat([sheet[keep], new]).sort_values(["dataset", "section", "state"]).to_csv(sheet_path, sep="\t", index=False)


if __name__ == "__main__":
    main()
