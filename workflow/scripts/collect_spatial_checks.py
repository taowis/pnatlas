"""Collect the per-section checks of the spatial mapping into one table and
apply the plan's rule: a state that fails check (a) or (b) in more than
half of the sections where it could be checked is "not mapped".

Inputs: the checks.tsv files written by c2l_map.py (one per Visium
section, under <c2l_dir>/<dataset>/<section>/) and the transfer_checks.tsv
files of the panel datasets (<results_dir>/<dataset>/transfer_checks.tsv).
Writes <results_dir>/section_checks.tsv (every row, with dataset and
technology) and <results_dir>/state_mapping_status.tsv (per state: sections
checked, failures of (a) and (b), mapped or not, by compartment of the
sections).

Usage: python workflow/scripts/collect_spatial_checks.py --c2l-dir <dir> --results-dir results/spatial \
           --sections results/spatial/section_summary.tsv
"""
import argparse
import glob
import pathlib

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--c2l-dir", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--sections", required=True)
    args = ap.parse_args()
    res = pathlib.Path(args.results_dir)
    secs = pd.read_csv(args.sections, sep="\t")
    comp = {ds: ("nerve" if "sural" in t else "ganglion")
            for ds, t in zip(secs["dataset"], secs["tissue"])} if "tissue" in secs else \
        {ds: ("nerve" if ("sural" in ds or "nerve" in ds) else "ganglion") for ds in secs["dataset"].unique()}

    frames = []
    for p in sorted(glob.glob(str(pathlib.Path(args.c2l_dir) / "*" / "*" / "checks.tsv"))):
        df = pd.read_csv(p, sep="\t")
        df.insert(0, "dataset", pathlib.Path(p).parts[-3])
        df = df.rename(columns={"dominant_spots": "assigned"})
        df["method"] = "cell2location"
        frames.append(df)
    for p in sorted(glob.glob(str(res / "*" / "transfer_checks.tsv"))):
        df = pd.read_csv(p, sep="\t")
        df.insert(0, "dataset", pathlib.Path(p).parent.name)
        df = df.rename(columns={"assigned_cells": "assigned"})
        df["method"] = "scANVI transfer"
        frames.append(df)
    allc = pd.concat(frames, ignore_index=True)
    allc["compartment"] = allc["dataset"].map(comp)
    allc.to_csv(res / "section_checks.tsv", sep="\t", index=False)

    checked = allc[allc["pass_a"].notna() | allc["pass_b"].notna()].copy()
    for c in ["pass_a", "pass_b"]:
        checked[c] = checked[c].map({True: True, False: False, "True": True, "False": False})
    rows = []
    for (state, compartment), g in checked.groupby(["state", "compartment"]):
        n = len(g)
        fa = int((g["pass_a"] == False).sum())
        fb = int((g["pass_b"] == False).sum())
        na, nb = int(g["pass_a"].notna().sum()), int(g["pass_b"].notna().sum())
        not_mapped = (na and fa > na / 2) or (nb and fb > nb / 2)
        rows.append([state, compartment, n, na, fa, nb, fb, round(float(g["auroc_markers"].median()), 3),
                     round(float(g["spearman_heldout"].median()), 3), "not mapped" if not_mapped else "mapped"])
    status = pd.DataFrame(rows, columns=["state", "compartment", "sections_checked", "checked_a", "failed_a",
                                         "checked_b", "failed_b", "median_auroc", "median_spearman", "status"])
    status.to_csv(res / "state_mapping_status.tsv", sep="\t", index=False)
    print(status.to_string(index=False))


if __name__ == "__main__":
    main()
