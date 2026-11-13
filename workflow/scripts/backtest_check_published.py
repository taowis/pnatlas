"""Check that feature set A at PNPrank's global seed reproduces its published back-test.

Compares the seed rows of backtest_seeds.tsv (feature set A) with the rows of
PNPrank's backtest_performance.tsv for the primary evaluation (all post-cutoff
positives, within-fold percentile, full feature set): AUROC to three decimals,
average precision to four, and the positive counts. Exits 1 on any difference.
"""
import argparse
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--published", required=True, help="backtest_performance.tsv of PNPrank at the tag")
    ap.add_argument("--seeds", required=True, help="backtest_seeds.tsv of this project")
    ap.add_argument("--seed", type=int, default=20261006)
    a = ap.parse_args()
    pub = pd.read_csv(a.published, sep="\t")
    pub = pub[(pub["eval"] == "all_after") & (pub["transform"] == "fold_rank") & (pub["score"] == "full")]
    ours = pd.read_csv(a.seeds, sep="\t")
    ours = ours[(ours["seed"] == a.seed) & (ours["features"] == "A")]
    m = pub.merge(ours, on=["positives", "cutoff"], suffixes=("_pub", "_ours"), how="left")
    m["auroc_ours"] = m["auroc_ours"].round(3)
    m["avg_precision_ours"] = m["avg_precision_ours"].round(4)
    cols = ["positives", "cutoff", "n_pos_pub", "n_pos_ours", "auroc_pub", "auroc_ours",
            "avg_precision_pub", "avg_precision_ours"]
    print(m[cols].to_string(index=False))
    same = ((m["auroc_pub"] == m["auroc_ours"]) & (m["avg_precision_pub"] == m["avg_precision_ours"])
            & (m["n_pos_pub"] == m["n_pos_ours"]))
    print(f"{int(same.sum())} of {len(m)} rows equal the published values")
    sys.exit(0 if same.all() else 1)


if __name__ == "__main__":
    main()
