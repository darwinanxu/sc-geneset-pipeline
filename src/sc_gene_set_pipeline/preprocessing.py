from anndata import AnnData
import scanpy as sc

def basic_qc_filter(
    adata: AnnData,
    min_genes: int = 200,
    min_cells: int = 3,
) -> AnnData:
    adata = adata.copy()
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata

def normalize_log1p(adata: AnnData, target_sum: float = 1e4) -> AnnData:
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)
    return adata

def add_basic_qc_metrics(adata: AnnData) -> AnnData:
    adata = adata.copy()
    adata.obs["n_counts"] = adata.X.sum(axis=1).A1 if hasattr(adata.X, "A1") else adata.X.sum(axis=1)
    adata.obs["n_genes"] = (adata.X > 0).sum(axis=1).A1 if hasattr((adata.X > 0).sum(axis=1), "A1") else (adata.X > 0).sum(axis=1)
    return adata