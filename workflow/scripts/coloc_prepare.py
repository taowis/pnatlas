"""Prepare the colocalisation of one GWAS with the tibial nerve eQTL (plan,
section 6).

The GWAS files are the GWAS Catalog harmonised files in the GWAS-SSF
layout (chromosome, base_pair_location, effect_allele, other_allele,
beta or odds_ratio, standard_error, effect_allele_frequency, p_value,
rsid; some studies add n, num_cases, num_controls and the confidence
interval). The effect is beta where given and log odds ratio otherwise;
the standard error is standard_error where given and otherwise derived
from the 95 percent confidence interval of the odds ratio. The
harmonisation code is not used: alleles are matched as an unordered
pair and coloc uses squared effects, so strand and orientation do not
matter.

Lead variants: every row below gwas_p, clumped by distance (the strongest
first, then any other lead at least window_bp away). Region: lead
position plus and minus window_bp. In each region the GWAS rows (tabix)
and the eQTL Catalogue rows (tabix, gene expression dataset) are read
and matched on chromosome, position and allele pair; a gene is a
candidate when it shares at least min_shared variants with the GWAS.
Genes with a SuSiE credible set in the region are flagged and their
per-variant log Bayes factors are pulled from the lbf file in one pass
over it.

Writes to <out-dir>/<gwas id>/: regions.tsv, pairs.tsv, and per region
gwas_<region>.tsv (variant, chromosome, position, beta, se, p, maf) and
per pair eqtl_<region>_<gene>.tsv (variant, beta, se, maf, an) and, where
a credible set exists, lbf_<region>_<gene>.tsv (variant, lbf columns).
study.tsv holds the sample size (the meta yaml, identical sample
entries counted once) and the case count: the meta yaml where it gives
one, otherwise the median of the per-variant num_cases column where the
file has one; a study flagged case-control without a case count is
treated as quantitative by coloc_run.R.
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
    """Rows below p_max and, from the first chunk, the medians of the per-variant n, num_cases and num_controls."""
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    idcol = "rsid" if "rsid" in header else "variant_id"
    extra = [c for c in ["n", "num_cases", "num_controls"] if c in header]
    cols = ["chromosome", "base_pair_location", idcol, "p_value"] + extra
    leads, medians = [], {}
    for i, chunk in enumerate(pd.read_csv(path, sep="\t", usecols=cols, chunksize=2_000_000, low_memory=False)):
        if i == 0:
            medians = {c: float(pd.to_numeric(chunk[c], errors="coerce").median()) for c in extra}
        pv = pd.to_numeric(chunk["p_value"], errors="coerce")
        leads.append(chunk[pv < p_max].assign(p_value=pv[pv < p_max]))
    leads = pd.concat(leads, ignore_index=True)
    leads = leads.dropna(subset=["chromosome", "base_pair_location"])
    leads["chromosome"] = leads["chromosome"].astype(str)
    leads["base_pair_location"] = leads["base_pair_location"].astype(int)
    leads["lead_id"] = leads[idcol].astype(str)
    return leads.sort_values("p_value"), medians


def gwas_effects(g):
    """beta and se columns from a GWAS-SSF frame: beta or log odds ratio; standard error or from the CI."""
    beta = pd.to_numeric(g["beta"], errors="coerce") if "beta" in g else pd.Series(np.nan, index=g.index)
    if "odds_ratio" in g:
        orr = pd.to_numeric(g["odds_ratio"], errors="coerce")
        beta = beta.where(beta.notna(), np.log(orr.where(orr > 0)))
    se = pd.to_numeric(g["standard_error"], errors="coerce") if "standard_error" in g else pd.Series(np.nan, index=g.index)
    if "ci_upper" in g and "ci_lower" in g:
        up, lo = pd.to_numeric(g["ci_upper"], errors="coerce"), pd.to_numeric(g["ci_lower"], errors="coerce")
        se_ci = (np.log(up.where(up > 0)) - np.log(lo.where(lo > 0))) / (2 * 1.959964)
        se = se.where(se.notna() & (se > 0), se_ci)
    return beta, se


def clump(leads, window):
    kept = []
    for _, r in leads.iterrows():
        if all(not (k["chromosome"] == r["chromosome"] and abs(k["base_pair_location"] - r["base_pair_location"]) < window)
               for k in kept):
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
    samples = meta.get("samples", [{}]) or [{}]
    seen, uniq = set(), []
    for smp in samples:                      # some meta files list the same sample set twice
        k = (str(smp.get("sample_ancestry_category")), smp.get("sample_size"))
        if k not in seen:
            seen.add(k)
            uniq.append(smp)
    n_total = sum(int(smp.get("sample_size", 0) or 0) for smp in uniq) or None
    n_case = sum(int(smp.get("case_count", 0) or 0) for smp in uniq) or None
    is_cc = bool(meta.get("case_control_study", False)) or any(smp.get("case_control_study") for smp in uniq) \
        or n_case is not None
    trait = meta.get("trait_description", "")
    trait = "; ".join(trait) if isinstance(trait, list) else trait

    leads, medians = read_leads(args.gwas, cc["gwas_p"])
    if n_case is None and medians.get("num_cases"):
        n_case = int(round(medians["num_cases"]))
        if n_total is None and medians.get("n"):
            n_total = int(round(medians["n"]))
    pd.DataFrame([{"gwas_id": args.gwas_id, "sample_size": n_total, "case_count": n_case,
                   "case_control": is_cc, "build": meta.get("genome_assembly", ""), "trait": trait,
                   "median_n": medians.get("n"), "median_num_cases": medians.get("num_cases")}]).to_csv(
        out / "study.tsv", sep="\t", index=False)
    leads = clump(leads, cc["window_bp"])
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
        chrom, pos = lead["chromosome"], int(lead["base_pair_location"])
        start, end = max(pos - cc["window_bp"], 1), pos + cc["window_bp"]
        rid = f"r{i:03d}_{chrom}_{pos}"
        g = fetch(gtbx, chrom, start, end, gwas_header)
        if g.empty:
            continue
        g["beta"], g["se"] = gwas_effects(g)
        for c in ["base_pair_location", "p_value", "effect_allele_frequency"]:
            g[c] = pd.to_numeric(g[c], errors="coerce")
        g = g.dropna(subset=["base_pair_location", "beta", "se"])
        g = g[g["se"] > 0]
        g["variant"] = [key(c, p, a, b) for c, p, a, b in zip(g["chromosome"], g["base_pair_location"],
                                                                  g["other_allele"], g["effect_allele"])]
        g = g.drop_duplicates("variant")
        eaf = g["effect_allele_frequency"]
        g["maf"] = np.where(eaf > 0.5, 1 - eaf, eaf)
        gw = pd.DataFrame({"variant": g["variant"], "chromosome": g["chromosome"],
                           "position": g["base_pair_location"].astype(int), "beta": g["beta"], "se": g["se"],
                           "p": g["p_value"], "maf": g["maf"]})
        gw.to_csv(out / f"gwas_{rid}.tsv", sep="\t", index=False)
        e = fetch(etbx, chrom, start, end, eqtl_header)
        if e.empty:
            regions.append([rid, chrom, start, end, lead["lead_id"], lead["p_value"], 0])
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
        regions.append([rid, chrom, start, end, lead["lead_id"], lead["p_value"], n_genes])
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
