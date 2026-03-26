import time

import numpy as np
import pandas as pd
from anndata import AnnData

from .base import BaseScorer
from ..types import ScoreResult


class ZScoreMeanScorer(BaseScorer):
    name = "zscore_mean"

    def score(self, adata: AnnData, gene_sets: dict) -> ScoreResult:
        start = time.time()

        X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
        X = np.asarray(X, dtype=float)
        gene_to_idx = {gene: i for i, gene in enumerate(adata.var_names)}

        gene_means = X.mean(axis=0)
        gene_stds = X.std(axis=0)
        gene_stds[gene_stds == 0] = 1.0
        z_matrix = (X - gene_means) / gene_stds

        scores = {}
        genes_per_set = {}

        for set_name, genes in gene_sets.items():
            valid_genes = [g for g in genes if g in gene_to_idx]
            if len(valid_genes) == 0:
                continue

            idx = [gene_to_idx[g] for g in valid_genes]
            scores[set_name] = z_matrix[:, idx].mean(axis=1)
            genes_per_set[set_name] = len(valid_genes)

        score_df = pd.DataFrame(scores, index=adata.obs_names)
        runtime = time.time() - start

        return ScoreResult(
            method=self.name,
            score_matrix=score_df,
            runtime_sec=runtime,
            metadata={
                "n_gene_sets": len(score_df.columns),
                "genes_per_set": genes_per_set,
            },
        )
