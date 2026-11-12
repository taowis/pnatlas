"""Prepare the colocalisation of one GWAS with the tibial nerve eQTL (plan,
section 6).

Lead variants: every row of the harmonised GWAS file below gwas_p, clumped
by distance (the strongest first, then any other lead at least window_bp
away). Region: lead position plus and minus window_bp. In each region the
GWAS rows (tabix) and the eQTL Catalogue rows (tabix, gene expression
dataset) are read; variants are matched on chromosome, position and
allele pair in either orientation, which is enough because coloc uses
squared effects; a gene is a candidate when it shares at least
min_shared variants with the GWAS. Genes with a SuSiE credible set in
the region are flagged and their per-variant log Bayes factors are
pulled from the lbf file in one pass over it.

Writes to <out-dir>/<gwas id>/: regions.tsv, pairs.tsv, and per region
gwas_<region>.tsv (variant, chromosome, position, beta, se, p, maf) and
per pair eqtl_<region>_<gene>.tsv (variant, beta, se, maf, an) and, where
a credible set exists, lbf_<region>_<gene>.tsv (variant, lbf columns).
The GWAS meta yaml gives the sample size and, for a case-control study,
the case count, written to study.tsv.
Usage: python workflow/scripts/coloc_prepare.py --gwas-id <id> --gwas <GCST.h.tsv.gz> --meta <GCST.h.tsv.gz-meta.yaml> \
           --eqtl <QTD000286.all.tsv.gz> --credible-sets <QTD000286.credible_sets.tsv.gz> \
           --lbf <QTD000286.lbf_variable.txt.gz> --out-dir <dir>
"""
import argparse
import gzip
import pathlib

import numpy as np
import pandas as pd
import pysam
import yaml

MIN_SHARED = 50


def read_leads(path, p_max):
    cols = ["hm_chrom", "hm_pos", "hm_variant_id", "p_value"]
    leads = []
    for chunk in pd.read_csv(path, sep="\t", usecols=cols, chunksize=2_000_000, low_memory=False):
        chunk = chunk[pd.to_numeric(chunk["p_value"], errors="coerce") < p_max]
        leads.append(chunk)
    leads = pd.concat(leads, ignore_index=True)
    leads["p_value"] = leads["p_value"].astype(float)
    leads = leads.dropna(subset=["hm_chrom", "hm_pos"])
    leads["hm_chrom"] = leads["hm_chrom"].astype(str)
    leads["hm_pos"] = leads["hm_pos"].astype(int)
    return leads.sort_values("p_value")


def clump(leads, window):
    kept = []
    for _, r in leads.iterrows():
        if all(not (k["hm_chrom"] == r["hm_chrom"] and abs(k["hm_pos"] - r["hm_pos"]) < window) for k in kept):
            kept.append(r)
    return pd.DataFrame(kept)


def fetch(tbx, chrom, start, end, header):
    rows = []
    for c in ({chrom, f"chr{chrom}", chrom.replace("chr", "")}):
        try:
            rows = [line.split("\t") for line in tbx.fetch(c, max(start, 1), end)]
        except ValueError:
            continue
        if rows:
            break
    return pd.DataFrame(rows, columns=header) if rows else pd.DataFrame(columns=header)


def key(chrom, pos, a1, a2):
    a, b = sorted([str(a1).upper(), str(a2).upper()])
    return f"{str(chrom).replace('chr', '')}:{int(pos)}:{a}:{b}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gwas-id", required=True)
    ap.add_argument("--gwas", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--eqtl", required=True)
    ap.add_argument("--credible-sets", required=True)
    ap.add_argument("--lbf", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--max-regions", type=int, default=0, help="probe: stop after this many regions")
    args = ap.parse_args()
    cc = yaml.safe_load(open(args.config))["genetics"]["coloc"]
    out = pathlib.Path(args.out_dir) / args.gwas_id
    out.mkdir(parents=True, exist_ok=True)

    meta = yaml.safe_load(open(args.meta))
    samples = meta.get("samples", [{}])
    n_total = sum(int(s.get("sample_size", 0) or 0) for s in samples) or None
    n_case = sum(int(s.get("case_count", 0) or 0) for s in samples) or None
    is_cc = bool(meta.get("case_control_study", False)) or n_case is not None
    pd.DataFrame([{"gwas_id": args.gwas_id, "sample_size": n_total, "case_count": n_case,
                   "case_control": is_cc, "build": meta.get("genome_assembly", ""),
                   "trait": meta.get("trait_description", "")}]).to_csv(out / "study.tsv", sep="\t", index=False)

    leads = clump(read_leads(args.gwas, cc["gwas_p"]), cc["window_bp"])
    if args.max_regions:
        leads = leads.head(args.max_regions)
    print(f"{args.gwas_id}: {len(leads)} lead regions", flush=True)
    if len(leads) == 0:
        pd.DataFrame(columns=["region", "chromosome", "start", "end", "lead_variant", "lead_p"]).to_csv(
            out / "regions.tsv", sep="\t", index=False)
        pd.DataFrame(columns=["region", "gene", "n_shared", "has_cs"]).to_csv(out / "pairs.tsv", sep="\t", index=False)
        return

    with gzip.open(args.gwas, "rt") as f:
        gwas_header = f.readline().rstrip("\n").split("\t")
    with gzip.open(args.eqtl, "rt") as f:
        eqtl_header = f.readline().rstrip("\n").split("\t")
    gtbx, etbx = pysam.TabixFile(args.gwas), pysam.TabixFile(args.eqtl)
    cs = pd.read_csv(args.credible_sets, sep="\t")
    cs_genes_by_region = {}
    regions, pairs, lbf_wanted = [], [], {}
    for i, (_, lead) in enumerate(leads.iterrows(), 1):
        chrom, pos = lead["hm_chrom"], int(lead["hm_pos"])
        start, end = max(pos - cc["window_bp"], 1), pos + cc["window_bp"]
        rid = f"r{i:03d}_{chrom}_{pos}"
        g = fetch(gtbx, chrom, start, end, gwas_header)
        if g.empty:
            continue
        g = g[["hm_chrom", "hm_pos", "hm_other_allele", "hm_effect_allele", "hm_beta", "standard_error", "p_value",
               "hm_effect_allele_frequency"]].copy()
        for c in ["hm_pos", "hm_beta", "standard_error", "p_value", "hm_effect_allele_frequency"]:
            g[c] = pd.to_numeric(g[c], errors="coerce")
        g = g.dropna(subset=["hm_pos", "hm_beta", "standard_error"])
        g["variant"] = [key(c, p, a, b) for c, p, a, b in zip(g["hm_chrom"], g["hm_pos"], g["hm_other_allele"],
                                                                  g["hm_effect_allele"])]
        g = g.drop_duplicates("variant")
        eaf = g["hm_effect_allele_frequency"]
        g["maf"] = np.where(eaf > 0.5, 1 - eaf, eaf)
        gw = pd.DataFrame({"variant": g["variant"], "chromosome": g["hm_chrom"], "position": g["hm_pos"].astype(int),
                           "beta": g["hm_beta"], "se": g["standard_error"], "p": g["p_value"], "maf": g["maf"]})
        gw.to_csv(out / f"gwas_{rid}.tsv", sep="\t", index=False)
        e = fetch(etbx, chrom, start, end, eqtl_header)
        if e.empty:
            regions.append([rid, chrom, start, end, lead["hm_variant_id"], lead["p_value"], 0])
            continue
        e = e[["molecular_trait_id", "gene_id", "chromosome", "position", "ref", "alt", "beta", "se", "maf", "an"]].copy()
        for c in ["position", "beta", "se", "maf", "an"]:
            e[c] = pd.to_numeric(e[c], errors="coerce")
        e = e.dropna(subset=["beta", "se"])
        e["variant"] = [key(c, p, a, b) for c, p, a, b in zip(e["chromosome"], e["position"], e["ref"], e["alt"])]
        shared = set(gw["variant"])
        e = e[e["variant"].isin(shared)]
        n_genes = 0
        cs_here = cs[(cs["variant"].str.startswith(f"chr{str(chrom).replace('chr', '')}_"))] if "variant" in cs else cs
        for gene, eg in e.groupby("gene_id"):
            eg = eg.drop_duplicates("variant")
            if len(eg) < MIN_SHARED:
                continue
            has_cs = False
            if "gene_id" in cs_here:
                c = cs_here[cs_here["gene_id"] == gene]
                if len(c):
                    cpos = c["variant"].str.split("_").str[1].astype(int)
                    has_cs = bool(((cpos >= start) & (cpos <= end)).any())
            eg[["variant", "beta", "se", "maf", "an"]].to_csv(out / f"eqtl_{rid}_{gene}.tsv", sep="\t", index=False)
            pairs.append([rid, gene, len(eg), has_cs])
            if has_cs:
                lbf_wanted.setdefault(gene, []).append((rid, set(eg["variant"])))
            n_genes += 1
        regions.append([rid, chrom, start, end, lead["hm_variant_id"], lead["p_value"], n_genes])
        print(f"  {rid}: {len(gw)} GWAS variants, {n_genes} candidate genes", flush=True)
    pd.DataFrame(regions, columns=["region", "chromosome", "start", "end", "lead_variant", "lead_p", "genes"]).to_csv(
        out / "regions.tsv", sep="\t", index=False)
    pairs = pd.DataFrame(pairs, columns=["region", "gene", "n_shared", "has_cs"])
    pairs.to_csv(out / "pairs.tsv", sep="\t", index=False)
    print(f"{len(pairs)} region by gene pairs, {int(pairs['has_cs'].sum())} with a credible set", flush=True)

    if lbf_wanted:
        wanted = set(lbf_wanted)
        found = {}
        with gzip.open(args.lbf, "rt") as f:
            header = f.readline().rstrip("\n").split("\t")
            gi = header.index("molecular_trait_id")
            for line in f:
                p = line.rstrip("\n").split("\t")
                if p[gi] in wanted:
                    found.setdefault(p[gi], []).append(p)
        lcols = [h for h in header if h.startswith("lbf_variable")]
        n_written = 0
        for gene, rows in found.items():
            df = pd.DataFrame(rows, columns=header)
            df["variant"] = [key(c, p, r, a) for c, p, r, a in zip(df["chromosome"], df["position"],
                                                                      df["variant"].str.split("_").str[2],
                                                                      df["variant"].str.split("_").str[3])]
            for rid, variants in lbf_wanted[gene]:
                sub = df[df["variant"].isin(variants)][["variant"] + lcols]
                if len(sub):
                    sub.to_csv(out / f"lbf_{rid}_{gene}.tsv", sep="\t", index=False)
                    n_written += 1
        print(f"log Bayes factors written for {n_written} pairs", flush=True)


if __name__ == "__main__":
    main()
