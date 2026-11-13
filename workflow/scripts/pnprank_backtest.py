"""Temporal back-test of PNPrank with and without the pnatlas features (plan,
section 7), with PNPrank's own evaluation code at its submission tag.

Feature set A is PNPrank's static feature matrix with its co-expression
label features; feature set B is A plus the pnatlas gene features. For
every positive set (tier 1; tiers 1 and 2), cutoff year and fold seed, the
genes labelled 0 are scored out of fold by `pnprank.backtest.
oof_scores_at_cutoff` with the same folds for A and B (the folds depend on
the labels and the seed only, which is asserted), the scores become
percentiles within the fold model that scored them, and the post-cutoff
positives of the set are evaluated against the non-candidate genes:
AUROC, average precision and the recall of the positives in the top
fraction of the config. The seeds are PNPrank's global seed and the next
ones, as many as its repeated-CV setting (the ten seeds of its config).

Primary comparison (plan): tier-1 AUROC averaged over the cutoff years.
For each seed the mean over cutoffs is taken for A and for B; the noise
baseline is the range of A's means over the seeds; an improvement is
claimed when B's mean over seeds exceeds A's by more than that range and
a paired sign-flip permutation test over the cutoff years (the mean over
seeds of B minus A at each cutoff; exact, one-sided) gives p below 0.05.
The same statistics are given for the top-fraction recall and for the
tier 1 and 2 set.

The PNPrank checkout must be at the tag named in the config (its code,
configs and truth set unchanged against the tag, checked with git);
PNPrank's config gives the atlases, cutoffs, folds and seeds. Scores are
cached per job under <out>/scores/.

Writes to --out: backtest_seeds.tsv (one row per positives, cutoff, seed
and feature set), backtest_by_cutoff.tsv (means over seeds),
backtest_summary.tsv (the primary and secondary comparisons with the
verdict), model_params.tsv (the C chosen per fold).
Usage: python workflow/scripts/pnprank_backtest.py --pnprank-repo <checkout> --features-static <features_static.tsv.gz> \
           --pnatlas-features <pnatlas_gene_features.tsv.gz> --out results/backtest [--n-jobs N] [--probe]
"""
import argparse
import itertools
import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd
import yaml
from joblib import Parallel, delayed
from sklearn.metrics import average_precision_score

TIERS = {"tier1": {"1"}, "tier1+2": {"1", "2"}}


def check_repo(repo, tag, commit):
    head = subprocess.run(["git", "-C", repo, "rev-parse", tag + "^{commit}"], capture_output=True, text=True, check=True).stdout.strip()
    if not head.startswith(commit):
        raise SystemExit(f"{tag} of {repo} is {head[:7]}, not {commit}")
    diff = subprocess.run(["git", "-C", repo, "diff", "--stat", tag, "--", "src", "workflow", "config", "data/truth",
                           "pyproject.toml"], capture_output=True, text=True, check=True).stdout.strip()
    if diff:
        raise SystemExit(f"{repo} differs from {tag} in code, config or truth:\n{diff}")
    return head


def sign_flip_p(d):
    """Exact one-sided sign-flip test of mean(d) > 0 over all 2^n sign assignments."""
    d = np.asarray(d, dtype=float)
    obs = d.mean()
    n = 0
    total = 2 ** len(d)
    for signs in itertools.product([1, -1], repeat=len(d)):
        if (d * np.array(signs)).mean() >= obs - 1e-12:
            n += 1
    return n / total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pnprank-repo", required=True)
    ap.add_argument("--features-static", required=True)
    ap.add_argument("--pnatlas-features", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--n-jobs", type=int, default=1)
    ap.add_argument("--probe", action="store_true", help="one seed, the first cutoff, tier 1 only")
    args = ap.parse_args()
    bt = yaml.safe_load(open(args.config))["backtest"]
    repo = pathlib.Path(args.pnprank_repo)
    head = check_repo(str(repo), bt["pnprank_tag"], bt["pnprank_commit"])
    sys.path.insert(0, str(repo / "src"))
    from pnprank.backtest import (auroc, evaluation_labels, oof_scores_at_cutoff, training_labels,
                                  within_fold_percentile)
    from pnprank.features import label_features, profiles_from_static

    pcfg = yaml.safe_load((repo / "config" / "config.yaml").read_text())
    fc, mc, bc = pcfg["features"], pcfg["model"], pcfg["backtest"]
    base = pcfg["seeds"]["global"]
    seeds = [base + i for i in range(bc["repeated_cv_seeds"])]
    cutoffs = list(bc["cutoff_years"])
    sets = list(TIERS)
    if args.probe:
        seeds, cutoffs, sets = seeds[:1], cutoffs[:1], ["tier1"]
    out = pathlib.Path(args.out)
    (out / "scores").mkdir(parents=True, exist_ok=True)

    XA = pd.read_csv(args.features_static, sep="\t", index_col=0)
    universe = list(XA.index)
    P = pd.read_csv(args.pnatlas_features, sep="\t", index_col=0).reindex(universe)
    XB = pd.concat([XA, P], axis=1)
    profiles = profiles_from_static(XA, fc["atlases"])
    truth = pd.read_csv(repo / pcfg["truth"]["path"], sep="\t", dtype=str)
    cands = pd.read_csv(repo / pcfg["truth"]["candidates"], sep="\t", dtype=str)
    negatives = set(universe) - set(cands.ensembl_gene_id) - set(truth.ensembl_gene_id)
    print(f"PNPrank {bt['pnprank_tag']} {head[:7]}: {len(universe)} genes, A {XA.shape[1]} features, "
          f"B {XB.shape[1]} features ({P.notna().any(axis=1).sum()} genes with pnatlas features); "
          f"seeds {seeds}, cutoffs {cutoffs}", flush=True)

    def label_fn(known):
        return label_features(profiles, universe, known)

    def run_job(set_name, cutoff, seed, fs):
        f = out / "scores" / f"{set_name}_{cutoff}_{seed}_{fs}.tsv.gz"
        if f.exists():
            return None
        y = training_labels(universe, truth, cutoff, TIERS[set_name])
        X = XA if fs == "A" else XB
        scores, params, fold = oof_scores_at_cutoff(X, y, "logistic", seed, bc["scoring_folds"], mc["inner_folds"],
                                                    label_fn, 1)
        pd.DataFrame({"fold": fold, "score": scores}).to_csv(f, sep="\t", compression="gzip")
        return [{"positives": set_name, "cutoff": cutoff, "seed": seed, "features": fs, "fold": k,
                 "C": p["logisticregression__C"], "n_pos_train": int(y.sum())} for k, p in enumerate(params)]

    jobs = [(s, c, sd, fs) for s in sets for c in cutoffs for sd in seeds for fs in ("A", "B")]
    res = Parallel(n_jobs=args.n_jobs, verbose=5)(delayed(run_job)(*j) for j in jobs)
    params = [r for rr in res if rr for r in rr]
    pf = out / "model_params.tsv"
    if params:
        old = pd.read_csv(pf, sep="\t") if pf.exists() else pd.DataFrame()
        pd.concat([old, pd.DataFrame(params)], ignore_index=True).drop_duplicates().to_csv(pf, sep="\t", index=False)

    rows = []
    for s, c, sd in itertools.product(sets, cutoffs, seeds):
        dfs = {fs: pd.read_csv(out / "scores" / f"{s}_{c}_{sd}_{fs}.tsv.gz", sep="\t", index_col=0) for fs in ("A", "B")}
        assert dfs["A"]["fold"].equals(dfs["B"]["fold"]), f"folds differ between A and B at {s} {c} {sd}"
        ev = evaluation_labels(truth, negatives, c, TIERS[s], None)
        ev = ev[ev.index.isin(dfs["A"].index)]
        lab = ev.to_numpy().astype(bool)
        for fs, df in dfs.items():
            pct = within_fold_percentile(df["score"], df["fold"]).reindex(ev.index).to_numpy(float)
            ok = ~np.isnan(pct)
            rows.append({"positives": s, "cutoff": c, "seed": sd, "features": fs, "n_pos": int(lab.sum()),
                         "n_neg": int((~lab).sum()), "auroc": auroc(pct, lab),
                         "avg_precision": average_precision_score(lab[ok], pct[ok]),
                         "top_recall": float((pct[lab] >= 1 - bt["top_fraction"]).mean())})
    seeds_df = pd.DataFrame(rows)
    seeds_df.round(4).to_csv(out / "backtest_seeds.tsv", sep="\t", index=False)
    by_cut = seeds_df.groupby(["positives", "cutoff", "features"]).agg(
        n_pos=("n_pos", "first"), auroc=("auroc", "mean"), auroc_min=("auroc", "min"), auroc_max=("auroc", "max"),
        top_recall=("top_recall", "mean"), avg_precision=("avg_precision", "mean")).reset_index()
    by_cut.round(4).to_csv(out / "backtest_by_cutoff.tsv", sep="\t", index=False)

    summary = []
    for s in sets:
        for metric in ["auroc", "top_recall", "avg_precision"]:
            per_seed = seeds_df[seeds_df["positives"] == s].groupby(["features", "seed"])[metric].mean().unstack(0)
            mean_a, mean_b = per_seed["A"].mean(), per_seed["B"].mean()
            range_a = per_seed["A"].max() - per_seed["A"].min()
            per_cut = seeds_df[seeds_df["positives"] == s].groupby(["cutoff", "features"])[metric].mean().unstack(1)
            d = (per_cut["B"] - per_cut["A"]).to_numpy()
            p = sign_flip_p(d) if len(d) > 1 else np.nan
            summary.append({"positives": s, "metric": metric, "cutoffs": len(d), "seeds": per_seed.shape[0],
                            "mean_A": mean_a, "min_A": per_seed["A"].min(), "max_A": per_seed["A"].max(),
                            "range_A": range_a, "mean_B": mean_b, "diff": mean_b - mean_a,
                            "cutoffs_B_better": int((d > 0).sum()), "p_sign_flip": p,
                            "improvement": bool((mean_b - mean_a) > range_a and p == p and p < 0.05),
                            "primary": bool(s == "tier1" and metric == "auroc")})
    summary = pd.DataFrame(summary)
    summary.round(4).to_csv(out / "backtest_summary.tsv", sep="\t", index=False)
    print(by_cut.round(3).to_string(index=False), flush=True)
    print(summary.round(4).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
