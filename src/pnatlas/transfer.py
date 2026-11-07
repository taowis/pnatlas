"""Label transfer by majority vote of the nearest reference cells in an
integrated embedding. A cell takes the state of most of its k neighbours
when that state holds at least min_conf of them; otherwise it takes the
majority class with the suffix "unassigned", or "unassigned" alone when
no class reaches min_conf either."""
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

UNASSIGNED = "unassigned"


def vote_codes(codes, n_cat):
    """codes: integer array of shape (cells, k). Returns the winning code per
    cell and the fraction of neighbours that carry it. Ties go to the lowest
    code, which is deterministic."""
    n, k = codes.shape
    counts = np.zeros((n, n_cat), dtype=np.int32)
    np.add.at(counts, (np.repeat(np.arange(n), k), codes.ravel()), 1)
    win = counts.argmax(axis=1)
    return win, counts[np.arange(n), win] / k


def vote(nbr_state, nbr_class, min_conf):
    """nbr_state, nbr_class: string arrays of shape (cells, k), the states and
    classes of each cell's neighbours. Returns a DataFrame with transfer_state,
    transfer_conf, state_final and class_final."""
    nbr_state, nbr_class = np.asarray(nbr_state, dtype=object), np.asarray(nbr_class, dtype=object)
    states, s_codes = np.unique(nbr_state, return_inverse=True)
    classes, c_codes = np.unique(nbr_class, return_inverse=True)
    s_codes, c_codes = s_codes.reshape(nbr_state.shape), c_codes.reshape(nbr_class.shape)
    s_win, s_conf = vote_codes(s_codes, len(states))
    c_win, c_conf = vote_codes(c_codes, len(classes))
    # the class of the winning state, read off the neighbours that carry it
    state_class = {}
    for s, c in zip(nbr_state.ravel(), nbr_class.ravel()):
        state_class.setdefault(s, c)
    t_state = states[s_win]
    final_state = np.where(s_conf >= min_conf, t_state,
                           np.where(c_conf >= min_conf, classes[c_win].astype(object) + " " + UNASSIGNED, UNASSIGNED))
    final_class = np.where(s_conf >= min_conf, [state_class[s] for s in t_state],
                           np.where(c_conf >= min_conf, classes[c_win], UNASSIGNED))
    return pd.DataFrame({"transfer_state": t_state, "transfer_conf": s_conf,
                         "state_final": final_state, "class_final": final_class})


def knn_transfer(ref_emb, ref_state, ref_class, query_emb, k, min_conf, n_jobs=8):
    """Exact k nearest neighbours of each query cell among the reference cells
    (Euclidean distance in the embedding), then the vote."""
    nn = NearestNeighbors(n_neighbors=k, algorithm="brute", n_jobs=n_jobs).fit(ref_emb)
    _, idx = nn.kneighbors(query_emb)
    ref_state, ref_class = np.asarray(ref_state, dtype=object), np.asarray(ref_class, dtype=object)
    return vote(ref_state[idx], ref_class[idx], min_conf)
