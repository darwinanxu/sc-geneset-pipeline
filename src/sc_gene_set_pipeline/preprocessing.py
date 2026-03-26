import numpy as np
from anndata import AnnData
from scipy import sparse


def _sum_axis(X, axis: int) -> np.ndarray:
    values = X.sum(axis=axis)
    return np.asarray(values).ravel()


def _nnz_axis(X, axis: int) -> np.ndarray:
    if sparse.issparse(X):
        return np.asarray((X > 0).sum(axis=axis)).ravel()
    return np.count_nonzero(X > 0, axis=axis)


def basic_qc_filter(
    adata: AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
) -> AnnData:
    """
    Apply basic cell/gene filtering.

    Parameters
    ----------
    min_genes
        Minimum number of detected genes per cell.
    min_cells
        Minimum number of cells per gene.
    """
    adata = adata.copy()
    X = adata.X

    genes_per_cell = _nnz_axis(X, axis=1)
    cell_mask = genes_per_cell >= min_genes
    adata = adata[cell_mask].copy()
    if adata.n_obs == 0:
        raise ValueError("No cells remain after filtering. Lower `min_genes` or inspect the input data.")

    X = adata.X
    cells_per_gene = _nnz_axis(X, axis=0)
    gene_mask = cells_per_gene >= min_cells
    adata = adata[:, gene_mask].copy()
    if adata.n_vars == 0:
        raise ValueError("No genes remain after filtering. Lower `min_cells` or inspect the input data.")

    return adata


def normalize_log1p(
    adata: AnnData,
    target_sum: float = 1e4,
) -> AnnData:
    """
    Normalize counts per cell and log-transform.
    """
    adata = adata.copy()
    X = adata.X

    counts_per_cell = _sum_axis(X, axis=1)
    scale = np.ones_like(counts_per_cell, dtype=float)
    nonzero_mask = counts_per_cell > 0
    scale[nonzero_mask] = target_sum / counts_per_cell[nonzero_mask]

    if sparse.issparse(X):
        X = sparse.diags(scale) @ X
        X.data = np.log1p(X.data)
    else:
        X = np.asarray(X, dtype=float) * scale[:, None]
        X = np.log1p(X)

    adata.X = X
    return adata


def add_basic_qc_metrics(adata: AnnData) -> AnnData:
    """
    Add simple QC-like summary columns to adata.obs.

    Adds:
    - n_counts
    - n_genes
    - sparsity
    """
    adata = adata.copy()

    X = adata.X
    counts_per_cell = _sum_axis(X, axis=1)
    genes_per_cell = _nnz_axis(X, axis=1)

    n_total_genes = adata.n_vars
    sparsity_values = 1.0 - (genes_per_cell / n_total_genes)

    adata.obs["n_counts"] = counts_per_cell
    adata.obs["n_genes"] = genes_per_cell
    adata.obs["sparsity"] = sparsity_values

    return adata


def run_basic_preprocessing(
    adata: AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
    target_sum: float = 1e4,
) -> AnnData:
    """
    Convenience wrapper for the default preprocessing steps.
    """
    adata = basic_qc_filter(adata, min_genes=min_genes, min_cells=min_cells)
    adata = normalize_log1p(adata, target_sum=target_sum)
    adata = add_basic_qc_metrics(adata)
    return adata
