# utils.py
# Shared helper functions used across all project scripts.

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp


def ensure_dirs(*dirs):
    """Create all listed directories if they do not already exist."""
    for d in dirs:
        os.makedirs(d, exist_ok=True)


def print_section(title):
    """Print a clearly visible section header to stdout."""
    width = 60
    print()
    print("=" * width)
    print(f"  {title}")
    print("=" * width)


def find_genes(var_names, genes, label="", verbose=True):
    """
    Locate gene names in an AnnData .var_names index.

    Tries exact match first, then upper/lower/title-case variants so the
    code works regardless of whether gene names are ALL_CAPS (human) or
    Title_case (mouse).

    Parameters
    ----------
    var_names : array-like
        adata.var_names from the dataset.
    genes : list of str
        Gene names to look up.
    label : str
        Optional label printed in warning messages.
    verbose : bool
        Whether to print warnings for missing genes.

    Returns
    -------
    found : list of str
        Gene names exactly as they appear in var_names.
    missing : list of str
        Genes that could not be located under any case variant.
    """
    var_set = set(var_names)
    found, missing = [], []
    for gene in genes:
        matched = None
        for candidate in [gene, gene.upper(), gene.lower(), gene.capitalize()]:
            if candidate in var_set:
                matched = candidate
                break
        if matched:
            found.append(matched)
        else:
            missing.append(gene)

    if missing and verbose:
        prefix = f"[{label}] " if label else ""
        print(f"  {prefix}Warning: genes not found in dataset: {missing}")
        print(f"  This may reflect genuine absence or low expression (dropout).")

    return found, missing


def to_dense(mat):
    """
    Convert a matrix to a dense numpy float32 array.

    Handles scipy sparse matrices, numpy arrays, and numpy matrices.
    """
    if sp.issparse(mat):
        return mat.toarray().astype(np.float32)
    return np.asarray(mat, dtype=np.float32)


def get_layer(adata, layer="log_norm"):
    """
    Return the expression matrix as a dense float32 array.

    Parameters
    ----------
    adata : AnnData
    layer : str
        'log_norm' → adata.layers['log_norm']  (log-normalised, unscaled)
        'X'        → adata.X                   (usually scaled for PCA)

    Returns
    -------
    np.ndarray, shape (n_cells, n_genes)
    """
    if layer == "X":
        return to_dense(adata.X)
    if layer not in adata.layers:
        raise KeyError(
            f"Layer '{layer}' not found in adata. "
            "Available layers: " + str(list(adata.layers.keys()))
        )
    return to_dense(adata.layers[layer])


def save_fig(fig, path, dpi=150):
    """Save a matplotlib figure and close it."""
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


def summarise_adata(adata, label=""):
    """Print a brief summary of an AnnData object."""
    prefix = f"[{label}] " if label else ""
    print(f"  {prefix}{adata.n_obs:,} cells  x  {adata.n_vars:,} genes")
    for col in ["leiden", "cell_type", "astrocyte_cluster", "astrocyte_subtype"]:
        if col in adata.obs.columns:
            n = adata.obs[col].nunique()
            print(f"  {prefix}'{col}' — {n} unique values")
