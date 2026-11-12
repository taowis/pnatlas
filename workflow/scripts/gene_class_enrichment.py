"""EWCE-style enrichment of neuropathy gene classes in the cell states of the
reference (plan, section 6).

For each gene class (the disease categories of the PNPrank truth set,
tiers 1 and 2, and the green genes of each PanelApp panel) and each state,
the statistic is the mean specificity of the class genes in that state.
Its null is the same mean over random gene sets of the same size drawn
from the background (the genes with a specificity and an annotation),
matched on exon length and GC content: the background is cut into
quantile bins of length and of GC, and each bootstrap draws, for every
class gene, one background gene from the same length by GC bin. The
p value is the fraction of bootstraps whose mean is at least the
observed one (with the +1 correction); the fold change is the observed
mean over the bootstrap mean; q values are Benjamini-Hochberg across
states within a class. A class is enriched in a state at q below
enrichment_q.

Writes to --results-dir: gene_class_enrichment.tsv (class, state, genes,
observed, expected, fold, p, q, enriched), gene_classes.tsv (the genes of
each class that are in the background, with their state of highest
specificity).
Usage: python workflow/scripts/gene_class_enrichment.py --specificity <specificity.tsv.gz> \
           --length-gc <ensembl105_gene_length_gc.tsv.gz> --truth data/metadata/pnprank_truth_set.tsv \
           --panels-dir data/metadata --results-dir results/genetics
"""
import argparse
import pathlib

import numpy as np
import pandas as pd
import yaml


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--specificity", required=True)
    ap.add_argument("--length-gc", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--panels-dir", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--bootstraps", type=int, default=None, help="override, for a probe run")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    g = cfg["genetics"]
    rng = np.random.default_rng(cfg["seed"])
    res = pathlib.Path(args.results_dir)
    res.mkdir(parents=True, exist_ok=True)
    n_boot = args.bootstraps or g["bootstraps"]

    spec = pd.read_csv(args.specificity, sep="\t", index_col=0)
    lg = pd.read_csv(args.length_gc, sep="\t")
    lg = lg[(lg["gene_biotype"] == "protein_coding") & lg["gc_percent"].notna() & (lg["exon_length"] > 0)]
    lg = lg.drop_duplicates("gene_name").set_index("gene_name")
    bg = spec.index.intersection(lg.index)
    spec = spec.loc[bg]
    lg = lg.loc[bg]
    nb = g["background_bins"]
    lbin = pd.qcut(lg["exon_length"].rank(method="first"), nb, labels=False)
    gbin = pd.qcut(lg["gc_percent"].rank(method="first"), nb, labels=False)
    cell = (lbin * nb + gbin).to_numpy()
    members = {c: np.where(cell == c)[0] for c in np.unique(cell)}
    print(f"background {len(bg)} protein coding genes with specificity, length and GC; {len(members)} bins", flush=True)

    classes = {}
    truth = pd.read_csv(args.truth, sep="\t")
    for c in g["classes"]:
        classes[f"truth:{c}"] = sorted(set(truth.loc[truth["disease_category"] == c, "hgnc_symbol"]))
    classes["truth:all"] = sorted(set(truth["hgnc_symbol"]))
    for p in g["panels"]:
        t = pd.read_csv(pathlib.Path(args.panels_dir) / f"panelapp_{p}.tsv", sep="\t")
        classes[f"panelapp:{p}"] = sorted(set(t.loc[t["confidence"] == "green", "gene_symbol"]))

    S = spec.to_numpy()
    pos = {gname: i for i, gname in enumerate(spec.index)}
    rows, members_rows = [], []
    for cname, genes in classes.items():
        idx = np.array([pos[x] for x in genes if x in pos])
        if len(idx) < 5:
            print(f"{cname}: {len(idx)} genes in the background, skipped", flush=True)
            continue
        obs = S[idx].mean(axis=0)
        cells_of = cell[idx]
        boots = np.empty((n_boot, S.shape[1]))
        for b in range(n_boot):
            draw = np.array([rng.choice(members[c]) for c in cells_of])
            boots[b] = S[draw].mean(axis=0)
        p = (1 + (boots >= obs).sum(axis=0)) / (n_boot + 1)
        exp = boots.mean(axis=0)
        q = bh(p)
        for j, state in enumerate(spec.columns):
            rows.append([cname, state, len(idx), round(float(obs[j]), 5), round(float(exp[j]), 5),
                         round(float(obs[j] / exp[j]), 3) if exp[j] > 0 else np.nan, float(p[j]), float(q[j]),
                         bool(q[j] < g["enrichment_q"])])
        top = spec.iloc[idx].idxmax(axis=1)
        for gname, st in top.items():
            members_rows.append([cname, gname, st, round(float(spec.loc[gname, st]), 4)])
        print(f"{cname}: {len(idx)} genes, enriched in {int((q < g['enrichment_q']).sum())} states", flush=True)
    out = pd.DataFrame(rows, columns=["class", "state", "genes", "observed", "expected", "fold", "p", "q", "enriched"])
    out.to_csv(res / "gene_class_enrichment.tsv", sep="\t", index=False)
    pd.DataFrame(members_rows, columns=["class", "gene", "top_state", "top_specificity"]).to_csv(
        res / "gene_classes.tsv", sep="\t", index=False)
    print(out[out["enriched"]].sort_values(["class", "q"]).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
