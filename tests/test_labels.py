import pandas as pd
import pytest

from pnatlas.labels import UNLABELLED, apply_celltype_map, read_celltype_map

ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]


def test_map_file_is_consistent():
    m = read_celltype_map(ROOT / "config/celltype_map.tsv")
    assert set(m.columns) >= {"dataset", "author_label_fine", "author_label", "class", "state", "fine"}
    assert (m["state"] != "").all() and (m["class"] != "").all()
    # a fine label carries one state and one class
    assert m.groupby(["dataset", "fine"])["state"].nunique().max() == 1
    assert m.groupby("state")["class"].nunique().max() == 1


def test_apply_map_and_unlabelled():
    cmap = pd.DataFrame({
        "dataset": ["a", "a"], "author_label_fine": ["x", "y"], "author_label": ["x", "y"],
        "class": ["neuron", "glia"], "state": ["s1", "s2"], "fine": ["f1", "f2"]})
    obs = pd.DataFrame({"dataset": ["a", "a", "b"], "author_label_fine": ["y", "x", ""]},
                       index=["c1", "c2", "c3"])
    out = apply_celltype_map(obs, cmap)
    assert out["state"].tolist() == ["s2", "s1", UNLABELLED]
    assert out["class"].tolist() == ["glia", "neuron", UNLABELLED]


def test_missing_label_raises():
    cmap = pd.DataFrame({"dataset": ["a"], "author_label_fine": ["x"], "author_label": ["x"],
                         "class": ["neuron"], "state": ["s1"], "fine": ["f1"]})
    obs = pd.DataFrame({"dataset": ["a"], "author_label_fine": ["z"]}, index=["c1"])
    with pytest.raises(ValueError, match="no map row"):
        apply_celltype_map(obs, cmap)
