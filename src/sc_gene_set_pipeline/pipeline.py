from typing import Iterable

import pandas as pd

from .gene_sets import filter_gene_sets_to_var_names, gene_set_overlap_frame
from .scoring.registry import get_scorer
from .evaluation.confounding import score_qc_correlations
from .evaluation.summary import summarize_method_performance


def run_pipeline(
    adata,
    gene_sets: dict[str, list[str]],
    methods: Iterable[str],
    qc_columns: Iterable[str] = ("n_counts", "n_genes"),
    min_gene_set_overlap: int = 1,
):
    methods = list(methods)
    if not methods:
        raise ValueError("At least one scoring method must be provided.")

    overlap_summary = gene_set_overlap_frame(gene_sets, adata.var_names)
    filtered_gene_sets = filter_gene_sets_to_var_names(
        gene_sets,
        adata.var_names,
        min_overlap=min_gene_set_overlap,
    )

    if not filtered_gene_sets:
        raise ValueError("No gene sets remain after filtering to dataset genes.")

    all_scores = {}
    all_summaries = []
    all_qc = {}

    for method in methods:
        scorer = get_scorer(method)
        result = scorer.score(adata, filtered_gene_sets)

        if result.score_matrix.empty:
            raise ValueError(f"Method '{method}' did not score any gene sets.")

        all_scores[method] = result.score_matrix

        qc_df = score_qc_correlations(result.score_matrix, adata.obs, qc_columns=qc_columns)
        all_qc[method] = qc_df

        genes_per_set = result.metadata.get("genes_per_set", {})
        mean_genes_per_set = (
            sum(genes_per_set.values()) / len(genes_per_set)
            if genes_per_set
            else 0.0
        )

        summary_df = summarize_method_performance(
            method_name=method,
            runtime_sec=result.runtime_sec,
            qc_corr_df=qc_df,
            n_gene_sets_scored=result.metadata.get("n_gene_sets"),
            n_gene_sets_requested=len(filtered_gene_sets),
            mean_genes_per_set=mean_genes_per_set,
        )
        all_summaries.append(summary_df)

    summary = pd.concat(all_summaries, axis=0, ignore_index=True)
    return {
        "scores": all_scores,
        "qc": all_qc,
        "summary": summary,
        "gene_set_overlap": overlap_summary,
        "filtered_gene_sets": filtered_gene_sets,
    }
