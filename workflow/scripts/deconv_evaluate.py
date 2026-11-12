"""Evaluate the stage 4 deconvolution (plan, section 5) and run the GTEx
covariate models.

Validation (a): for each method, the estimated proportions of the
pseudo-bulk mixtures against the truth: Pearson correlation across states
within each mixture (its mean over mixtures is the plan's statistic), the
root mean squared error over all mixture by state entries, and per state
the Pearson correlation across mixtures. A method passes when the mean
Pearson is above pass_pearson and the RMSE below pass_rmse; BayesPrism is
the primary method and is kept unless it fails and MuSiC passes.

Validation (b), for every method: among the GTEx tibial nerve samples,
the fat content the pathologist wrote in the notes (the largest
percentage near "fat" or "adipose", for example "20-30% fat" gives 30;
"no fat", "fat free" or "trace fat" give 0; notes that do not mention
fat are left out) against the adipocyte fraction, Spearman correlation
one-sided; and those whose
notes mention fibrosis (fibro) against the rest on the fibroblast
fraction (endoneurial plus epineurial fibroblast), one-sided Mann-Whitney
(Wilcoxon rank sum). A test passes at p below 0.05 in the expected
direction.

Marker concordance (when --gtex-counts and --specificity are given): for
each state the ten protein coding genes of highest specificity in the
reference among those expressed at 50 CPM or more in the state's
pseudobulk and present in the GTEx table, their mean log2 CPM in each GTEx sample as a marker
score, and the Spearman correlation of that score with the estimated
fraction of the state across samples, for every method. This checks that
a method tracks the signal that is in the bulk RNA, independently of the
pathology notes. The adipocyte marker score is also tested against the
fat percentage of the notes and the result added to the histology table
as method "markers": it is the ceiling of what any deconvolution can
reach on validation (b).

GTEx model, for every method: for each state, ordinary least squares of
the proportion on age (bracket midpoint), sex, ischaemic time (SMTSISCH,
minutes) and the Hardy scale (categorical, 0 to 4) together; coefficients,
p values and Benjamini-Hochberg q values across states for each term
within a method. The primary method's rows come first.

Writes to --results-dir: validation_pseudobulk.tsv, validation_per_state.tsv,
validation_histology.tsv, gtex_proportions.tsv (primary method, with the
covariates), gtex_covariate_models.tsv, method_choice.tsv, marker_concordance.tsv.
Usage: python workflow/scripts/deconv_evaluate.py --prep-dir <dir> --results-dir results/deconvolution \
           --sample-attributes <SampleAttributesDS.txt> --subject-phenotypes <SubjectPhenotypesDS.txt> \
           [--gtex-counts <prep dir>/gtex_counts.tsv.gz --specificity <interim>/genetics/specificity.tsv.gz \
            --expression <interim>/genetics/expression_cpm.tsv.gz --length-gc data/metadata/ensembl105_gene_length_gc.tsv.gz]
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
    q = np.full(len(p), np.nan)
    ok = np.where(np.isfinite(p))[0]      # a state with no variance has no p value
    n = len(ok)
    order = ok[np.argsort(p[ok])]
    q[order] = np.minimum(np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1], 1)
    return q


def fat_percent(note):
    """Fat content from a pathology note: the largest percentage within 40 characters before or 25 after each
    "fat" or "adipos" (a range gives its upper end, "<5%" gives 5); "no fat", "fat free" or "trace fat" without a
    percentage give 0; a note that does not mention fat gives NaN."""
    note = note.lower()
    best = np.nan
    for m in re.finditer(r"fat|adipos", note):
        window = note[max(0, m.start() - 40):m.end() + 25]
        cand = [float(x) for x in re.findall(r"(\d+(?:\.\d+)?)\s*%", window)]
        if not cand and re.search(r"no (significant )?(adherent |attached |external )?(fat|adipos)|fat[- ]free|trace", window):
            cand = [0.0]
        if cand:
            best = max(cand) if best != best else max(best, max(cand))
    return best


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
    ap.add_argument("--gtex-counts", default=None)
    ap.add_argument("--specificity", default=None)
    ap.add_argument("--expression", default=None, help="state pseudobulk CPM, same layout as the specificity")
    ap.add_argument("--length-gc", default=None, help="Ensembl table with gene_name and gene_biotype")
    ap.add_argument("--n-markers", type=int, default=10)
    ap.add_argument("--marker-min-cpm", type=float, default=50.0)
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

    # validation (b): histology notes, for every method
    notes = gtex["pathology_notes"]
    fat = notes.map(fat_percent)
    fibro_hit = notes.str.lower().str.contains("fibro", regex=False)
    hrows = []
    for m, path in [(mm, res / f"gtex_{mm}.tsv") for mm in METHODS]:
        if not path.exists():
            continue
        e = pd.read_csv(path, sep="\t", index_col=0).reindex(gtex.index)
        fib_states = [s for s in e.columns if "fibroblast" in s]
        if "adipocyte" in e:
            ok = fat.notna()
            rho, p = stats.spearmanr(fat[ok], e.loc[ok, "adipocyte"]) if ok.sum() > 10 else (np.nan, np.nan)
            p1 = p / 2 if rho > 0 else 1 - p / 2
            hrows.append([m, "fat percent", "adipocyte", "spearman", int(ok.sum()), int((~ok).sum()),
                          round(float(rho), 3), p1, bool(p1 < 0.05) if p1 == p1 else False])
        if fib_states:
            frac = e[fib_states].sum(axis=1)
            a, b = frac[fibro_hit], frac[~fibro_hit]
            u, p = stats.mannwhitneyu(a, b, alternative="greater") if len(a) and len(b) else (np.nan, np.nan)
            hrows.append([m, "fibrosis mentioned", "+".join(fib_states), "mann-whitney", int(len(a)), int(len(b)),
                          round(float(a.median() - b.median()), 4), p, bool(p < 0.05) if p == p else False])
    hist = pd.DataFrame(hrows, columns=["method", "note", "state", "test", "samples_with_note", "samples_other",
                                        "statistic", "p_one_sided", "pass"])
    hist.to_csv(res / "validation_histology.tsv", sep="\t", index=False)
    print(hist.to_string(index=False), flush=True)

    # marker concordance
    if args.gtex_counts and args.specificity:
        counts = pd.read_csv(args.gtex_counts, sep="\t", index_col=0)
        logcpm = np.log2(counts / counts.sum(axis=0) * 1e6 + 1)
        spec = pd.read_csv(args.specificity, sep="\t", index_col=0)
        spec = spec.loc[spec.index.intersection(logcpm.index)]
        if args.length_gc:
            lg = pd.read_csv(args.length_gc, sep="\t")
            spec = spec.loc[spec.index.isin(lg.loc[lg["gene_biotype"] == "protein_coding", "gene_name"])]
        expr = pd.read_csv(args.expression, sep="\t", index_col=0).reindex(spec.index) if args.expression else None
        crows = []
        for m in METHODS:
            path = res / f"gtex_{m}.tsv"
            if not path.exists():
                continue
            e = pd.read_csv(path, sep="\t", index_col=0)
            common = e.index.intersection(logcpm.columns)
            for s in e.columns:
                if s not in spec.columns:
                    continue
                cand = spec[s] if expr is None else spec.loc[expr[s] >= args.marker_min_cpm, s]
                markers = cand.sort_values(ascending=False).head(args.n_markers).index
                score = logcpm.loc[markers, common].mean(axis=0)
                frac = e.loc[common, s]
                rho, pv = stats.spearmanr(score, frac) if frac.std() > 0 else (np.nan, np.nan)
                crows.append([m, s, ";".join(markers), round(float(rho), 3) if rho == rho else np.nan, pv,
                              round(float(frac.mean()), 4), round(float((frac == 0).mean()), 3)])
                if s == "adipocyte" and m == METHODS[0]:
                    okf = fat.reindex(common).notna()
                    r2, p2 = stats.spearmanr(fat.reindex(common)[okf], score[okf])
                    p2 = p2 / 2 if r2 > 0 else 1 - p2 / 2
                    hist.loc[len(hist)] = ["markers", "fat percent", "adipocyte marker score", "spearman",
                                           int(okf.sum()), int((~okf).sum()), round(float(r2), 3), p2, bool(p2 < 0.05)]
                    hist.to_csv(res / "validation_histology.tsv", sep="\t", index=False)
        conc = pd.DataFrame(crows, columns=["method", "state", "markers", "spearman", "p", "mean_fraction",
                                            "fraction_zero"])
        conc.to_csv(res / "marker_concordance.tsv", sep="\t", index=False)
        print(conc.pivot(index="state", columns="method", values="spearman").to_string(), flush=True)

    # covariate models, for every method
    ok = gtex[["sex", "age", "hardy", "ischaemic_min"]].notna().all(axis=1)
    g = gtex[ok]
    hardy = pd.get_dummies(g["hardy"].astype(int), prefix="hardy", drop_first=True).astype(float)
    X = np.column_stack([g["age"].values, (g["sex"] == 2).astype(float).values, g["ischaemic_min"].values / 60.0,
                         hardy.values])
    names = ["age_years", "female", "ischaemic_hours"] + list(hardy.columns)
    varying = X.std(axis=0) > 0           # a Hardy level absent from the samples would make the design singular
    X, names = X[:, varying], [n for n, v in zip(names, varying) if v]
    tables = []
    for m in [primary] + [mm for mm in METHODS if mm != primary]:
        path = res / f"gtex_{m}.tsv"
        if not path.exists():
            continue
        e = pd.read_csv(path, sep="\t", index_col=0).reindex(g.index)
        mrows = []
        for s in e.columns:
            d = ols(e[s].values.astype(float), X, names)
            d.insert(0, "state", s)
            mrows.append(d)
        models = pd.concat(mrows, ignore_index=True)
        models = models[models["term"] != "intercept"].copy()
        models.insert(0, "method", m)
        models["q"] = np.nan
        for term in names:
            sel = models["term"] == term
            models.loc[sel, "q"] = bh(models.loc[sel, "p"])
        models["samples"] = len(g)
        tables.append(models)
    models = pd.concat(tables, ignore_index=True)
    models.to_csv(res / "gtex_covariate_models.tsv", sep="\t", index=False)
    for m in models["method"].unique():
        sig = models[(models["method"] == m) & (models["q"] < 0.05)]
        print(f"{m}: covariate models on {len(g)} samples, {len(sig)} state by term pairs at q < 0.05", flush=True)
        print(sig.sort_values("q").head(25).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
