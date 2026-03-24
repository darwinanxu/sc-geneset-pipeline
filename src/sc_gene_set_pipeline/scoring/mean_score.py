import time
import pandas as pd
import numpy as np
from anndata import AnnData
from .base import BaseScorer
from ..types import ScoreResult

class MeanExpressionScorer(BaseScorer):
    name = "mean_score"

    def score(self, adata: AnnData, gene_sets: dict) -> ScoreResult:
        start = time.time()

        gene_to_idx = {gene: i for i, gene in enumerate(adata.var_names)}
        scores = {}

        for set_name, genes in gene_sets.items():
            valid_genes = [g for g in genes if g in gene_to_idx]
            if len(valid_genes) == 0:
                continue

            idx = [gene_to_idx[g] for g in valid_genes]
            subset = adata[:, idx].X
            if hasattr(subset, "mean"):
                score = np.asarray(subset.mean(axis=1)).ravel()
            else:
                score = subset.mean(axis=1)

            scores[set_name] = score

        score_df = pd.DataFrame(scores, index=adata.obs_names)
        runtime = time.time() - start

        return ScoreResult(
            method=self.name,
            score_matrix=score_df,
            runtime_sec=runtime,
            metadata={"n_gene_sets": len(score_df.columns)}
        )