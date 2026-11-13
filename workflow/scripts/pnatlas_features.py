"""Gene features of the pnatlas reference for the PNPrank v2 back-test (plan,
section 7).

One row per Ensembl gene id of the PNPrank universe (the index of its
static feature matrix) that the reference knows; columns:

- pnatlas|spec|<state> and pnatlas|logcpm|<state>: the specificity and
  log1p CPM of the gene in every state of the reference (stage 5
  pseudobulk tables; reference symbols mapped to Ensembl 105 ids, an id
  of the universe preferred when a symbol has several, then a protein
  coding one).
- pnatlas|compartment|<compartment>: the sum of the gene's specificity
  over the states mapped in that compartment in stage 3 (status table,
  status "mapped").
- pnatlas|deconv|<nerve state>: the Spearman correlation across the GTEx
  tibial nerve samples between the gene's log CPM and the estimated
  fraction of the state (stage 4; the fractions of the method named by
  --fractions, see decisions.md).
- pnatlas|coloc|tested, pnatlas|coloc|max_pp_h4, pnatlas|coloc|colocalised:
  whether the gene was in any colocalisation pair of stage 5, its largest
  PP.H4 over pairs (0 when not tested) and whether that exceeds the
  threshold of the config.

Usage: python workflow/scripts/pnatlas_features.py --specificity <specificity.tsv.gz> --expression <expression_cpm.tsv.gz> \
           --length-gc data/metadata/ensembl105_gene_length_gc.tsv.gz --status results/spatial/state_mapping_status.tsv \
           --gtex-counts <interim/deconv/gtex_counts.tsv.gz> --fractions results/deconvolution/gtex_bayesprism.tsv \
           --coloc results/genetics/coloc_all_pairs.tsv --universe <features_static.tsv.gz> --out <pnatlas_gene_features.tsv.gz>
"""
import argparse
import gzip

import numpy as np
import pandas as pd
import yaml


def symbol_to_id(symbols, lg, universe):
    """Map symbols to Ensembl ids: an id of the universe first, then protein coding, then the first id."""
    groups = lg.groupby("gene_name")
    out = {}
    for s in symbols:
        if s not in groups.groups:
            continue
        g = lg.loc[groups.groups[s]].sort_values("ensembl_gene_id")
        inu = g[g["ensembl_gene_id"].isin(universe)]
        pc = g[g["gene_biotype"] == "protein_coding"]
        pick = inu if len(inu) else (pc if len(pc) else g)
        out[s] = pick["ensembl_gene_id"].iloc[0]
    return out


def rank_rows(a):
    """Average ranks along axis 1 (for a Spearman correlation as Pearson of ranks)."""
    return pd.DataFrame(a).rank(axis=1).to_numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--specificity", required=True)
    ap.add_argument("--expression", required=True)
    ap.add_argument("--length-gc", required=True)
    ap.add_argument("--status", required=True)
    ap.add_argument("--gtex-counts", required=True)
    ap.add_argument("--fractions", required=True)
    ap.add_argument("--coloc", required=True)
    ap.add_argument("--universe", required=True, help="PNPrank static feature matrix (its index is the universe)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()
    cfg = yaml.safe_load(open(args.config))
    pp_h4 = cfg["genetics"]["coloc"]["pp_h4"]

    with gzip.open(args.universe, "rt") as f:
        universe = [line.split("\t", 1)[0] for line in f][1:]
    universe_set = set(universe)
    lg = pd.read_csv(args.length_gc, sep="\t")
    spec = pd.read_csv(args.specificity, sep="\t", index_col=0)
    expr = pd.read_csv(args.expression, sep="\t", index_col=0).reindex(spec.index)
    ids = symbol_to_id(spec.index, lg, universe_set)
    keep = [s for s in spec.index if s in ids]
    print(f"reference genes {len(spec)}, with an Ensembl id {len(keep)}, in the universe "
          f"{sum(ids[s] in universe_set for s in keep)}", flush=True)
    spec, expr = spec.loc[keep], expr.loc[keep]
    index = pd.Index([ids[s] for s in keep], name="ensembl_gene_id")
    feats = {}
    for st in spec.columns:
        feats[f"pnatlas|spec|{st}"] = spec[st].to_numpy()
        feats[f"pnatlas|logcpm|{st}"] = np.log1p(expr[st].to_numpy())
    F = pd.DataFrame(feats, index=index)
    F = F[~F.index.duplicated()]

    status = pd.read_csv(args.status, sep="\t")
    for comp, states in status[status["status"] == "mapped"].groupby("compartment")["state"]:
        cols = [f"pnatlas|spec|{s}" for s in states if f"pnatlas|spec|{s}" in F]
        F[f"pnatlas|compartment|{comp}"] = F[cols].sum(axis=1)
        print(f"compartment {comp}: {len(cols)} mapped states", flush=True)

    counts = pd.read_csv(args.gtex_counts, sep="\t", index_col=0)
    frac = pd.read_csv(args.fractions, sep="\t", index_col=0)
    samples = counts.columns.intersection(frac.index)
    logcpm = np.log1p(counts[samples] / counts[samples].sum(axis=0) * 1e6)
    gids = symbol_to_id(logcpm.index, lg, universe_set)
    gkeep = [s for s in logcpm.index if s in gids]
    R = rank_rows(logcpm.loc[gkeep].to_numpy())
    R = (R - R.mean(axis=1, keepdims=True)) / R.std(axis=1, keepdims=True)
    Q = rank_rows(frac.loc[samples].T.to_numpy())
    Q = (Q - Q.mean(axis=1, keepdims=True)) / Q.std(axis=1, keepdims=True)
    rho = pd.DataFrame(R @ Q.T / R.shape[1], index=[gids[s] for s in gkeep], columns=frac.columns)
    rho = rho[~rho.index.duplicated()]
    for st in frac.columns:
        F[f"pnatlas|deconv|{st}"] = rho[st].reindex(F.index)
    print(f"deconvolution correlations over {len(samples)} samples for {len(rho)} genes", flush=True)

    coloc = pd.read_csv(args.coloc, sep="\t")
    best = coloc.groupby("gene")["pp_h4"].max()
    F["pnatlas|coloc|tested"] = F.index.isin(best.index).astype(float)
    F["pnatlas|coloc|max_pp_h4"] = best.reindex(F.index).fillna(0.0)
    F["pnatlas|coloc|colocalised"] = (F["pnatlas|coloc|max_pp_h4"] > pp_h4).astype(float)
    print(f"colocalisation: {int(F['pnatlas|coloc|tested'].sum())} genes tested, "
          f"{int(F['pnatlas|coloc|colocalised'].sum())} colocalised", flush=True)

    F = F.loc[F.index.isin(universe_set)].round(5)
    F.to_csv(args.out, sep="\t", compression="gzip" if args.out.endswith(".gz") else None)
    print(f"{F.shape[0]} genes x {F.shape[1]} features written to {args.out}", flush=True)


if __name__ == "__main__":
    main()
