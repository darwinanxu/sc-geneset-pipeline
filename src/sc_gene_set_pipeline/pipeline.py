from typing import List, Dict
import pandas as pd

from .gene_sets import filter_gene_sets_to_var_names
from .scoring.registry import get_scorer
from .evaluation.confounding import score_qc_correlations
from .evaluation.summary import summarize_method_performance

def run_pipeline(adata, gene_sets: Dict[str, list], methods: List[str]):
    gene_sets = filter_gene_sets_to_var_names(gene_sets, adata.var_names)

    all_scores = {}
    all_summaries = []
    all_qc = {}

    for method in methods:
        scorer = get_scorer(method)
        result = scorer.score(adata, gene_sets)

        all_scores[method] = result.score_matrix

        qc_df = score_qc_correlations(result.score_matrix, adata.obs)
        all_qc[method] = qc_df

        summary_df = summarize_method_performance(
            method_name=method,
            runtime_sec=result.runtime_sec,
            qc_corr_df=qc_df,
        )
        all_summaries.append(summary_df)

    summary = pd.concat(all_summaries, axis=0, ignore_index=True)
    return {
        "scores": all_scores,
        "qc": all_qc,
        "summary": summary,
    }