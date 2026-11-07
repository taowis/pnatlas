"""Apply config/celltype_map.tsv to the obs table of a reference dataset."""
import pandas as pd

LEVELS = ["class", "state", "fine"]
UNLABELLED = "unlabelled"


def read_celltype_map(path):
    m = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    dup = m.duplicated(["dataset", "author_label_fine"])
    if dup.any():
        raise ValueError(f"duplicate map rows: {m.loc[dup, ['dataset', 'author_label_fine']].values.tolist()}")
    return m


def apply_celltype_map(obs, cmap):
    """Return a frame with the three unified label levels for every row of obs.

    obs needs the columns dataset and author_label_fine. A dataset that has no
    rows in the map is unlabelled at every level; a dataset that has rows
    must have one for every label it carries, otherwise a ValueError lists
    the missing labels.
    """
    out = pd.DataFrame(index=obs.index)
    mapped = set(cmap["dataset"])
    for lv in LEVELS:
        out[lv] = UNLABELLED
    for ds in obs["dataset"].unique():
        if ds not in mapped:
            continue
        rows = obs["dataset"] == ds
        sub = cmap[cmap["dataset"] == ds].set_index("author_label_fine")
        labels = obs.loc[rows, "author_label_fine"].astype(str)
        missing = sorted(set(labels) - set(sub.index))
        if missing:
            raise ValueError(f"{ds}: no map row for {missing}")
        for lv in LEVELS:
            out.loc[rows, lv] = labels.map(sub[lv]).values
    return out
