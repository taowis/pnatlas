"""Collect the colocalisation results of every GWAS and assign the
colocalised eGenes to cell states (plan, section 6).

A pair is colocalised at PP.H4 above pp_h4. For each colocalised gene
the state of highest specificity in the reference and that specificity
are given, with the gene's specificity across all states kept in a
second table, and the spatial compartment where the top state is mapped
(from the stage 3 status table: the compartments in which the state is
"mapped"). Gene ids are mapped to symbols with the Ensembl table.

Writes to --results-dir: coloc_all_pairs.tsv (every pair of every GWAS),
coloc_colocalised.tsv (PP.H4 above the threshold, with the state
assignment), coloc_summary.tsv (per GWAS: regions, pairs, colocalised
genes).
Usage: python workflow/scripts/coloc_collect.py --coloc-dir <dir with <gwas id>/coloc.tsv and regions.tsv> \
           --specificity <specificity.tsv.gz> --length-gc <ensembl105_gene_length_gc.tsv.gz> \
           --status results/spatial/state_mapping_status.tsv --results-dir results/genetics
"""
import argparse
import pathlib

import pandas as pd
import yaml


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coloc-dir", required=True)
    ap.add_argument("--specificity", required=True)
    ap.add_argument("--length-gc", required=True)
    ap.add_argument("--status", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()
    cc = yaml.safe_load(open(args.config))["genetics"]["coloc"]
    res = pathlib.Path(args.results_dir)
    res.mkdir(parents=True, exist_ok=True)
    spec = pd.read_csv(args.specificity, sep="\t", index_col=0)
    ens = pd.read_csv(args.length_gc, sep="\t").set_index("ensembl_gene_id")["gene_name"]
    status = pd.read_csv(args.status, sep="\t")
    mapped = status[status["status"] == "mapped"].groupby("state")["compartment"].apply(lambda s: ", ".join(sorted(s)))

    frames, summary = [], []
    for d in sorted(pathlib.Path(args.coloc_dir).iterdir()):
        p = d / "coloc.tsv"
        if not p.exists():
            continue
        regions = pd.read_csv(d / "regions.tsv", sep="\t")
        pairs = pd.read_csv(d / "pairs.tsv", sep="\t")
        df = pd.read_csv(p, sep="\t") if p.stat().st_size > 1 else pd.DataFrame()
        n_coloc = 0
        if len(df):
            df.insert(0, "gwas", d.name)
            df = df.merge(regions[["region", "lead_variant", "lead_p"]], on="region", how="left")
            df["gene_symbol"] = df["gene"].map(ens).fillna("")
            frames.append(df)
            n_coloc = int((df["pp_h4"] > cc["pp_h4"]).sum())
        summary.append([d.name, len(regions), len(pairs), len(df), n_coloc])
    allp = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    allp.to_csv(res / "coloc_all_pairs.tsv", sep="\t", index=False)
    pd.DataFrame(summary, columns=["gwas", "regions", "pairs", "tested", "colocalised"]).to_csv(
        res / "coloc_summary.tsv", sep="\t", index=False)
    if len(allp):
        col = allp[allp["pp_h4"] > cc["pp_h4"]].copy()
        tops, specs, comps = [], [], []
        for sym in col["gene_symbol"]:
            if sym in spec.index:
                s = spec.loc[sym]
                tops.append(s.idxmax())
                specs.append(round(float(s.max()), 4))
                comps.append(mapped.get(s.idxmax(), "not mapped"))
            else:
                tops.append("")
                specs.append(float("nan"))
                comps.append("")
        col["top_state"] = tops
        col["top_specificity"] = specs
        col["top_state_mapped_in"] = comps
        col = col.sort_values(["gwas", "pp_h4"], ascending=[True, False])
        col.to_csv(res / "coloc_colocalised.tsv", sep="\t", index=False)
        print(col[["gwas", "region", "gene_symbol", "pp_h4", "method", "top_state", "top_specificity",
                   "top_state_mapped_in"]].to_string(index=False))
        spec.loc[[g for g in col["gene_symbol"].unique() if g in spec.index]].round(4).to_csv(
            res / "coloc_gene_specificity.tsv", sep="\t")
    print(pd.DataFrame(summary, columns=["gwas", "regions", "pairs", "tested", "colocalised"]).to_string(index=False))


if __name__ == "__main__":
    main()
