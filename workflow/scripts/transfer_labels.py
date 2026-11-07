"""Give every reference cell a state: the mapped author state where the
authors annotated the cell, and otherwise the state of its nearest
neighbours in the chosen embedding, within the cell's own compartment
(ganglion or nerve) and, where the authors gave a class, within that
class. In the sympathetic ganglion dataset, cells that land among neurons
are sympathetic neurons (the tissue has no sensory neurons and the
reference has no sympathetic ones); the catecholamine markers of those
cells are tabulated as the check. The reliability of the transfer is
measured by holding out whole batches of labelled cells and predicting
them from the rest. Finally the state rule of the plan (section 3) is
applied: at least state_min_cells cells from at least state_min_donors
known donors.

Usage: python workflow/scripts/transfer_labels.py --emb "$DATA_ROOT/interim/reference_embeddings.h5ad" \
           --counts "$DATA_ROOT/interim/reference_merged.h5ad" \
           --out "$DATA_ROOT/interim/reference_annotation.tsv.gz" --out-dir results/reference
"""
import argparse
import pathlib

import anndata as ad
import numpy as np
import pandas as pd
import yaml

from pnatlas.labels import UNLABELLED
from pnatlas.transfer import UNASSIGNED, knn_transfer


def seeded(state):
    st = state.astype(str)
    return ~((st == UNLABELLED) | st.str.endswith(UNASSIGNED))


def transfer_queries(obs, emb, cfg, n_jobs):
    """Transfer within compartment; cells with an author class (the
    "<class> unassigned" states) use only the reference cells of that class."""
    is_ref = seeded(obs["state"]).values
    out = pd.DataFrame(index=obs.index, columns=["transfer_state", "transfer_conf", "state_final", "class_final"])
    comp, klass, state = obs["compartment"].astype(str).values, obs["class"].astype(str).values, obs["state"].astype(str).values
    for c in np.unique(comp):
        for q_class in np.unique(klass[(comp == c) & ~is_ref]):
            q = np.where((comp == c) & ~is_ref & (klass == q_class))[0]
            r = (comp == c) & is_ref & ((klass == q_class) if q_class != UNLABELLED else True)
            r = np.where(r)[0]
            if len(r) < cfg["k"]:
                raise ValueError(f"{c}/{q_class}: {len(r)} reference cells")
            res = knn_transfer(emb[r], state[r], klass[r], emb[q], cfg["k"], cfg["min_conf"], n_jobs)
            out.iloc[q] = res.values
            print(f"transfer\t{c}\t{q_class}\t{len(q)} cells from {len(r)} reference cells", flush=True)
    return out


def holdout(obs, emb, cfg, seed, n_jobs):
    """Predict held-out batches of labelled cells from the other batches of
    the same compartment, all classes together (the unlabelled case)."""
    rng = np.random.default_rng(seed)
    is_ref = seeded(obs["state"]).values
    comp, state, klass = obs["compartment"].astype(str).values, obs["state"].astype(str).values, obs["class"].astype(str).values
    batch, ds = obs["batch"].astype(str).values, obs["dataset"].astype(str).values
    pred = pd.Series(index=obs.index, dtype=object)
    conf = pd.Series(index=obs.index, dtype=float)
    for c in np.unique(comp):
        cells = np.where((comp == c) & is_ref)[0]
        b = np.unique(batch[cells])
        fold_of = dict(zip(b, rng.permutation(len(b)) % cfg["n_folds"]))
        folds = np.array([fold_of[x] for x in batch[cells]])
        for f in range(cfg["n_folds"]):
            q, r = cells[folds == f], cells[folds != f]
            if not len(q):
                continue
            res = knn_transfer(emb[r], state[r], klass[r], emb[q], cfg["k"], cfg["min_conf"], n_jobs)
            pred.iloc[q], conf.iloc[q] = res["state_final"].values, res["transfer_conf"].values
            print(f"holdout\t{c}\tfold {f}\t{len(q)} cells", flush=True)
    rows = []
    for key, name in [("state", "state"), ("dataset", "dataset")]:
        col = state if key == "state" else ds
        for c in np.unique(comp):
            sel = (comp == c) & is_ref
            for v in np.unique(col[sel]):
                m = sel & (col == v)
                p, t = pred.values[m], state[m]
                assigned = conf.values[m] >= cfg["min_conf"]
                wrong = pd.Series(p[assigned & (p != t)]).value_counts()
                rows.append({"compartment": c, "level": name, "group": v, "n": int(m.sum()),
                             "assigned": round(float(assigned.mean()), 3),
                             "accuracy_assigned": round(float((p[assigned] == t[assigned]).mean()), 3) if assigned.any() else np.nan,
                             "accuracy_all": round(float((p == t).mean()), 3),
                             "top_confusion": wrong.index[0] if len(wrong) else "",
                             "top_confusion_n": int(wrong.iloc[0]) if len(wrong) else 0})
    for c in np.unique(comp):
        m = (comp == c) & is_ref
        p, t = pred.values[m], state[m]
        assigned = conf.values[m] >= cfg["min_conf"]
        rows.append({"compartment": c, "level": "all", "group": "all", "n": int(m.sum()),
                     "assigned": round(float(assigned.mean()), 3),
                     "accuracy_assigned": round(float((p[assigned] == t[assigned]).mean()), 3),
                     "accuracy_all": round(float((p == t).mean()), 3), "top_confusion": "", "top_confusion_n": 0})
    return pd.DataFrame(rows)


def sympathetic(obs, counts_path, cfg):
    """Cells of the sympathetic ganglion dataset that land among neurons are
    sympathetic neurons; the catecholamine markers are counted as the check."""
    sym = (obs["tissue"].astype(str) == cfg["sympathetic_tissue"]).values
    idx = np.where(sym)[0]
    if not len(idx):
        return obs, pd.DataFrame()
    neuron = sym & (obs["class_final"].astype(str) == "neuron").values & (obs["label_source"] == "transfer").values
    obs.loc[neuron, "state_final"] = cfg["sympathetic_state"]
    obs.loc[neuron, "label_source"] = "transfer, tissue"
    m = ad.read_h5ad(counts_path, backed="r")
    genes = [g for g in cfg["sympathetic_markers"] if g in m.var_names]
    sub = m[idx, genes].to_memory()
    x = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    det = pd.DataFrame(x >= 1, columns=genes, index=sub.obs_names)
    n_det = det.sum(axis=1)
    tab = pd.crosstab(obs.loc[sub.obs_names, "state_final"].astype(str), n_det)
    tab.columns = [f"markers_detected_{c}" for c in tab.columns]
    tab.insert(0, "cells", tab.sum(axis=1))
    for g in genes:
        tab[f"{g}_detected"] = det.groupby(obs.loc[sub.obs_names, "state_final"].astype(str).values)[g].sum()
    return obs, tab.reset_index().rename(columns={"index": "state_final"})


def state_rule(obs, cfg):
    o = obs[~obs["state_final"].astype(str).str.endswith(UNASSIGNED)]
    known = o["donor"].astype(str) != "unknown"
    donor_key = o["dataset"].astype(str) + ":" + o["donor"].astype(str)
    g = o.groupby(o["state_final"].astype(str))
    t = pd.DataFrame({"cells": g.size(),
                      "known_donors": donor_key[known].groupby(o["state_final"].astype(str)[known]).nunique(),
                      "cells_donor_unknown": (~known).groupby(o["state_final"].astype(str)).sum(),
                      "batches": g["batch"].nunique(), "datasets": g["dataset"].nunique(),
                      "class": g["class_final"].first()}).fillna({"known_donors": 0})
    t["known_donors"] = t["known_donors"].astype(int)
    t["passes"] = (t["cells"] >= cfg["state_min_cells"]) & (t["known_donors"] >= cfg["state_min_donors"])
    return t.sort_values(["class", "cells"], ascending=[True, False]).reset_index().rename(columns={"index": "state_final"})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-dir", default="results/reference")
    ap.add_argument("--config", default="config/config.yaml")
    ap.add_argument("--n-jobs", type=int, default=8)
    a = ap.parse_args()
    with open(a.config) as f:
        full = yaml.safe_load(f)
    cfg, seed = full["transfer"], full["seed"]

    e = ad.read_h5ad(a.emb)
    obs, emb = e.obs.copy(), np.asarray(e.obsm[cfg["embedding"]], dtype=np.float32)
    del e

    res = transfer_queries(obs, emb, cfg, a.n_jobs)
    is_ref = seeded(obs["state"]).values
    for c in res.columns:
        obs[c] = res[c].values
    obs["transfer_conf"] = pd.to_numeric(obs["transfer_conf"])
    obs.loc[is_ref, "state_final"] = obs.loc[is_ref, "state"].astype(str)
    obs.loc[is_ref, "class_final"] = obs.loc[is_ref, "class"].astype(str)
    obs["label_source"] = np.where(is_ref, "author", "transfer")
    obs, sym_tab = sympathetic(obs, a.counts, cfg)

    ho = holdout(obs, emb, cfg, seed, a.n_jobs)
    rule = state_rule(obs, cfg)

    out = pathlib.Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    ho.to_csv(out / "transfer_holdout.tsv", sep="\t", index=False)
    rule.to_csv(out / "state_rule.tsv", sep="\t", index=False)
    if len(sym_tab):
        sym_tab.to_csv(out / "sympathetic_neuron_markers.tsv", sep="\t", index=False)
    pd.crosstab(obs["state_final"].astype(str), obs["dataset"].astype(str)).to_csv(out / "final_state_by_dataset.tsv", sep="\t")
    ds = obs.groupby(obs["dataset"].astype(str))
    summary = pd.DataFrame({
        "cells": ds.size(),
        "author": ds["label_source"].apply(lambda s: int((s == "author").sum())),
        "transferred": ds.apply(lambda d: int(((d["label_source"] != "author") & ~d["state_final"].astype(str).str.endswith(UNASSIGNED)).sum())),
        "unassigned": ds.apply(lambda d: int(d["state_final"].astype(str).str.endswith(UNASSIGNED).sum())),
        "median_conf": ds.apply(lambda d: round(float(pd.to_numeric(d.loc[d["label_source"] != "author", "transfer_conf"]).median()), 3)
                                if (d["label_source"] != "author").any() else np.nan)})
    summary.reset_index().to_csv(out / "label_transfer_summary.tsv", sep="\t", index=False)

    keep = ["dataset", "sample", "donor", "batch", "tissue", "compartment", "protocol", "author_label",
            "author_label_fine", "class", "state", "fine", "scanvi_pred", "transfer_state", "transfer_conf",
            "class_final", "state_final", "label_source"]
    pathlib.Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    obs[keep].to_csv(a.out, sep="\t", compression="gzip")
    print(summary.to_string())
    print(f"states passing the rule: {int(rule['passes'].sum())} of {len(rule)}")


if __name__ == "__main__":
    main()
