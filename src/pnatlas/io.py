"""Readers for the large interim files."""
import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp
from anndata.io import read_elem


def read_counts_float32(path, cells=None):
    """Read a CSR h5ad of counts as float32 and keep only the given cells.
    The merged reference is a float64 CSR of 1.1e9 non-zeros; reading it
    through anndata takes twice that in memory, so the arrays are taken
    with h5py and the subset made before the object is built."""
    with h5py.File(path) as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        indptr = f["X/indptr"][:]
        data = f["X/data"][:].astype(np.float32)
        indices = f["X/indices"][:]
    x = sp.csr_matrix((data, indices, indptr), shape=(len(obs), len(var)))
    del data, indices
    if cells is not None:
        rows = np.where(obs.index.isin(cells))[0]
        x = x[rows]
        obs = obs.iloc[rows].copy()
    return ad.AnnData(X=x, obs=obs, var=var)
