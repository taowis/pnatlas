import numpy as np

from pnatlas.transfer import UNASSIGNED, knn_transfer, vote


def test_vote_levels():
    k = 5
    nbr_state = np.array([
        ["C-NP"] * 4 + ["C-PEP"],                      # clear state
        ["C-NP", "C-NP", "C-PEP", "C-PEP", "A-PEP"],    # class clear, state not
        ["C-NP", "C-NP", "mural", "mural", "macrophage"],  # nothing clear
    ])
    nbr_class = np.where(np.isin(nbr_state, ["C-NP", "C-PEP", "A-PEP"]), "neuron",
                         np.where(nbr_state == "mural", "vascular", "immune"))
    out = vote(nbr_state, nbr_class, min_conf=0.5)
    assert list(out["state_final"]) == ["C-NP", "neuron " + UNASSIGNED, UNASSIGNED]
    assert list(out["class_final"]) == ["neuron", "neuron", UNASSIGNED]
    assert out["transfer_conf"].tolist() == [4 / k, 2 / k, 2 / k]
    assert list(out["transfer_state"]) == ["C-NP", "C-NP", "C-NP"]


def test_knn_transfer_recovers_clusters():
    rng = np.random.default_rng(0)
    ref = np.vstack([rng.normal(0, 0.1, (50, 3)), rng.normal(5, 0.1, (50, 3))])
    state = np.array(["a"] * 50 + ["b"] * 50)
    klass = np.array(["x"] * 50 + ["y"] * 50)
    query = np.array([[0, 0, 0], [5, 5, 5], [2.5, 2.5, 2.5]], dtype=float)
    out = knn_transfer(ref, state, klass, query, k=10, min_conf=0.5, n_jobs=1)
    assert list(out["state_final"][:2]) == ["a", "b"]
    assert list(out["class_final"][:2]) == ["x", "y"]
