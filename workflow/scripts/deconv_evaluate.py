"""Evaluate the stage 4 deconvolution (plan, section 5) and run the GTEx
covariate models.

Validation (a): for each method, the estimated proportions of the
pseudo-bulk mixtures against the truth: Pearson correlation across states
within each mixture (its mean over mixtures is the plan's statistic), the
root mean squared error over all mixture by state entries, and per state
the Pearson correlation across mixtures. A method passes when the mean
Pearson is above pass_pearson and the RMSE below pass_rmse; BayesPrism is
the primary method and is kept unless it fails and MuSiC passes.

Validation (b): among the GTEx tibial nerve samples, those whose
pathology notes mention adipose tissue (adipos, fat) against the rest on
the adipocyte fraction, and those whose notes mention fibrosis (fibro)
against the rest on the fibroblast fraction (endoneurial plus epineurial
fibroblast), one-sided Mann-Whitney (Wilcoxon rank sum) in the expected
direction.

GTEx model: for each state, ordinary least squares of the proportion on
age (bracket midpoint), sex, ischaemic time (SMTSISCH, minutes) and the
Hardy scale (categorical, 0 to 4) together; coefficients, p values and
Benjamini-Hochberg q values across states for each term.

Writes to --results-dir: validation_pseudobulk.tsv, validation_per_state.tsv,
validation_histology.tsv, gtex_proportions.tsv (primary method, with the
covariates), gtex_covariate_models.tsv, method_choice.tsv.
Usage: python workflow/scripts/deconv_evaluate.py --prep-dir <dir> --results-dir results/deconvolution \
           --sample-attributes <SampleAttributesDS.txt> --subject-phenotypes <SubjectPhenotypesDS.txt>
"""
import argparse
import pathlib
import re

import numpy as np
import pandas as pd
import yaml
from scipy import stats

METHODS = ["bayesprism", "music"]


def bh(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    return np.minimum(q, 1)


def ols(y, X, names):
    X = np.column_stack([np.ones(len(y)), X])
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    s2 = resid @ resid / dof
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.diag(cov))
    t = beta / se
    p = 2 * stats.t.sf(np.abs(t), dof)
    return pd.DataFrame({"term": ["intercept"] + names, "coef": beta, "se": se, "p": p})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prep-dir", required=True)
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--sample-attributes", required=True)
    ap.add_argument("--subject-phenotypes", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    args = ap.parse_args()
    dc = yaml.safe_load(open(args.config))["deconvolution"]
    prep, res = pathlib.Path(args.prep_dir), pathlib.Path(args.results_dir)
    res.mkdir(parents=True, exist_ok=True)

    truth = pd.read_csv(prep / "mixtures_truth.tsv", sep="\t", index_col=0).drop(columns=["donor"])
    states = list(truth.columns)
    rows, per_state, available = [], [], {}
    for m in METHODS:
        p = res / f"validation_{m}.tsv"
        if not p.exists():
            continue
        est = pd.read_csv(p, sep="\t", index_col=0).reindex(index=truth.index, columns=states).fillna(0.0)
        available[m] = est
        r_mix = [stats.pearsonr(truth.loc[i], est.loc[i])[0] for i in truth.index]
        rmse = float(np.sqrt(((est.values - truth.values) ** 2).mean()))
        mean_r = float(np.mean(r_mix))
        rows.append([m, len(truth), round(mean_r, 3), round(float(np.median(r_mix)), 3), round(rmse, 4),
                     mean_r > dc["pass_pearson"] and rmse < dc["pass_rmse"]])
        for s in states:
            r = stats.pearsonr(truth[s], est[s])[0] if truth[s].std() > 0 else np.nan
            per_state.append([m, s, round(float(truth[s].mean()), 4), round(float(est[s].mean()), 4),
                              round(float(r), 3), round(float(np.sqrt(((est[s] - truth[s]) ** 2).mean())), 4)])
    val = pd.DataFrame(rows, columns=["method", "mixtures", "mean_pearson", "median_pearson", "rmse", "pass"])
    val.to_csv(res / "validation_pseudobulk.tsv", sep="\t", index=False)
    pd.DataFrame(per_state, columns=["method", "state", "true_mean", "estimated_mean", "pearson", "rmse"]).to_csv(
        res / "validation_per_state.tsv", sep="\t", index=False)
    print(val.to_string(index=False), flush=True)
    passed = {r.method: r["pass"] for _, r in val.iterrows()}
    primary = "bayesprism" if passed.get("bayesprism") else ("music" if passed.get("music") else None)
    pd.DataFrame([{"primary": primary or "none", "bayesprism_pass": passed.get("bayesprism"),
                   "music_pass": passed.get("music")}]).to_csv(res / "method_choice.tsv", sep="\t", index=False)
    if primary is None:
        print("no method passes validation (a); GTEx results are not interpreted", flush=True)
        return

    # GTEx samples with covariates
    attr = pd.read_csv(args.sample_attributes, sep="\t", low_memory=False)
    attr = attr[attr["SMTSD"] == dc["tissue"]].set_index("SAMPID")
    pheno = pd.read_csv(args.subject_phenotypes, sep="\t").set_index("SUBJID")
    est = pd.read_csv(res / f"gtex_{primary}.tsv", sep="\t", index_col=0)
    est = est[est.index.isin(attr.index)]
    subj = est.index.str.split("-").str[:2].str.join("-")
    meta = pd.DataFrame(index=est.index)
    meta["subject"] = subj
    meta["sex"] = pheno.reindex(subj)["SEX"].values          # 1 male, 2 female
    meta["age"] = pheno.reindex(subj)["AGE"].astype(str).map(
        lambda s: np.mean([float(x) for x in s.split("-")]) if "-" in s else np.nan).values
    meta["hardy"] = pheno.reindex(subj)["DTHHRDY"].values
    meta["ischaemic_min"] = attr.loc[est.index, "SMTSISCH"].values
    meta["pathology_notes"] = attr.loc[est.index, "SMPTHNTS"].fillna("").values
    gtex = pd.concat([meta, est], axis=1)
    gtex.to_csv(res / "gtex_proportions.tsv", sep="\t")
    print(f"GTEx {primary}: {len(gtex)} tibial nerve samples", flush=True)

    # validation (b): histology notes
    notes = gtex["pathology_notes"].str.lower()
    fib_states = [s for s in est.columns if "fibroblast" in s]
    tests = [("adipose", r"adipos|\bfat", "adipocyte", est["adipocyte"] if "adipocyte" in est else None),
             ("fibrosis", r"fibro", "+".join(fib_states), est[fib_states].sum(axis=1) if fib_states else None)]
    hrows = []
    for name, pat, state, frac in tests:
        if frac is None:
            continue
        hit = notes.str.contains(pat, regex=True)
        a, b = frac[hit], frac[~hit]
        u, p = stats.mannwhitneyu(a, b, alternative="greater") if len(a) and len(b) else (np.nan, np.nan)
        hrows.append([name, state, int(hit.sum()), int((~hit).sum()), round(float(a.median()), 4),
                      round(float(b.median()), 4), p, bool(p < 0.05) if p == p else False])
    hist = pd.DataFrame(hrows, columns=["note", "state", "samples_mentioned", "samples_other", "median_mentioned",
                                        "median_other", "p_greater", "pass"])
    hist.to_csv(res / "validation_histology.tsv", sep="\t", index=False)
    print(hist.to_string(index=False), flush=True)

    # covariate models
    ok = gtex[["sex", "age", "hardy", "ischaemic_min"]].notna().all(axis=1)
    g = gtex[ok]
    hardy = pd.get_dummies(g["hardy"].astype(int), prefix="hardy", drop_first=True).astype(float)
    X = np.column_stack([g["age"].values, (g["sex"] == 2).astype(float).values, g["ischaemic_min"].values / 60.0,
                         hardy.values])
    names = ["age_years", "female", "ischaemic_hours"] + list(hardy.columns)
    mrows = []
    for s in est.columns:
        d = ols(g[s].values.astype(float), X, names)
        d.insert(0, "state", s)
        mrows.append(d)
    models = pd.concat(mrows, ignore_index=True)
    models = models[models["term"] != "intercept"].copy()
    models["q"] = np.nan
    for term in names:
        m = models["term"] == term
        models.loc[m, "q"] = bh(models.loc[m, "p"])
    models["samples"] = len(g)
    models.to_csv(res / "gtex_covariate_models.tsv", sep="\t", index=False)
    sig = models[models["q"] < 0.05]
    print(f"covariate models on {len(g)} samples: {len(sig)} state by term pairs at q < 0.05", flush=True)
    print(sig.sort_values("q").head(20).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
