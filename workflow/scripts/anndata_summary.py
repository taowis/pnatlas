"""One row per h5ad in a directory: cells, genes, samples, donors, labelled
cells and file size. Reads the files in backed mode so it is cheap.

Usage: python workflow/scripts/anndata_summary.py "$DATA_ROOT/interim/reference" > results/reference/anndata_summary.tsv
"""
import glob
import os
import sys

import anndata as ad


def main(d):
    print("dataset\tcells\tgenes\tsamples\tdonors\tlabelled\tfile_bytes")
    for f in sorted(glob.glob(os.path.join(d, "*.h5ad"))):
        a = ad.read_h5ad(f, backed="r")
        o = a.obs
        print("\t".join(str(x) for x in [
            o["dataset"].iloc[0], a.n_obs, a.n_vars, o["sample"].nunique(), o["donor"].nunique(),
            int((o["author_label"] != "").sum()), os.path.getsize(f)]))


if __name__ == "__main__":
    main(sys.argv[1])
