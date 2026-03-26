import time
import pandas as pd
import numpy as np
from scipy.stats import rankdata
from anndata import AnnData
from .base import BaseScorer
from ..types import ScoreResult

class RankBasedScorer(BaseScorer):
    name = "rank_score"

    def score(self, adata: AnnData, gene_sets: dict) -> ScoreResult:
        start = time.time()

        X = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
        gene_to_idx = {gene: i for i, gene in enumerate(adata.var_names)}

        ranked = np.apply_along_axis(rankdata, 1, X)
        scores = {}
        genes_per_set = {}

        for set_name, genes in gene_sets.items():
            valid_genes = [g for g in genes if g in gene_to_idx]
            if len(valid_genes) == 0:
                continue

            idx = [gene_to_idx[g] for g in valid_genes]
            scores[set_name] = ranked[:, idx].mean(axis=1)
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
