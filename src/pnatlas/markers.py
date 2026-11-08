"""Top markers per label by a one-vs-rest Wilcoxon test on log-normalised
counts (scanpy rank_genes_groups): genes with adjusted p below a cut-off,
detected in at least a fraction of the label's cells, ranked by log fold
change."""
import pandas as pd
import scanpy as sc


def subsample_per_label(labels, n_per_label, rng):
    """Index values of up to n_per_label cells per label."""
    pick = []
    for _, g in pd.Series(labels).groupby(labels):
        n = min(len(g), n_per_label)
        pick.extend(rng.choice(g.index.values, n, replace=False))
    return pick


def top_markers(adata, label_key, padj, min_pct, top_n):
    """adata holds raw counts; returns a DataFrame with label, rank, gene,
    log2fc, pct_in_label, padj."""
    a = adata.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.tl.rank_genes_groups(a, label_key, method="wilcoxon", pts=True)
    rows = []
    for label in sorted(a.obs[label_key].unique()):
        df = sc.get.rank_genes_groups_df(a, label)
        pct = a.uns["rank_genes_groups"]["pts"][label]
        df["pct"] = pct.loc[df["names"]].values
        df = df[(df["pvals_adj"] < padj) & (df["pct"] >= min_pct) & (df["logfoldchanges"] > 0)]
        df = df.sort_values("logfoldchanges", ascending=False).head(top_n)
        for rank, r in enumerate(df.itertuples(), 1):
            rows.append([label, rank, r.names, round(float(r.logfoldchanges), 3), round(float(r.pct), 3),
                         float(r.pvals_adj)])
    return pd.DataFrame(rows, columns=["state", "rank", "gene", "log2fc", "pct_in_state", "padj"])
