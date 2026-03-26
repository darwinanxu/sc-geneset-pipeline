from __future__ import annotations

import pandas as pd


def summarize_method_performance(
    method_name: str,
    runtime_sec: float,
    auc_df: pd.DataFrame | None = None,
    qc_corr_df: pd.DataFrame | None = None,
    n_gene_sets_scored: int | None = None,
    n_gene_sets_requested: int | None = None,
    mean_genes_per_set: float | None = None,
) -> pd.DataFrame:
    summary = {
        "method": method_name,
        "runtime_sec": runtime_sec,
    }

    if n_gene_sets_scored is not None:
        summary["n_gene_sets_scored"] = n_gene_sets_scored

    if n_gene_sets_requested is not None:
        summary["n_gene_sets_requested"] = n_gene_sets_requested

    if mean_genes_per_set is not None:
        summary["mean_genes_per_set"] = mean_genes_per_set

    if auc_df is not None and not auc_df.empty:
        summary["mean_auroc"] = auc_df["auroc"].mean()
        summary["mean_auprc"] = auc_df["auprc"].mean()

    if qc_corr_df is not None and not qc_corr_df.empty:
        summary["mean_abs_qc_corr"] = qc_corr_df["spearman_rho"].abs().mean()

    return pd.DataFrame([summary])
