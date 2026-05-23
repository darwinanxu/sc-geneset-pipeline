from pathlib import Path
import sys

import numpy as np
import pandas as pd
from anndata import AnnData

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from sc_gene_set_pipeline.pipeline import run_pipeline
from sc_gene_set_pipeline.preprocessing import run_basic_preprocessing


def main():
    counts = np.array(
        [
            [5, 3, 0, 0, 1, 0],
            [4, 2, 0, 0, 2, 0],
            [0, 0, 6, 4, 0, 2],
            [0, 1, 5, 3, 0, 1],
        ],
        dtype=float,
    )
    obs = pd.DataFrame(index=["cell_1", "cell_2", "cell_3", "cell_4"])
    var = pd.DataFrame(index=["NKG7", "PRF1", "IFIT1", "ISG15", "IL7R", "MALAT1"])
    adata = AnnData(X=counts, obs=obs, var=var)

    gene_sets = {
        "cytotoxicity": ["NKG7", "PRF1"],
        "interferon_response": ["IFIT1", "ISG15"],
        "naive_t_cell": ["IL7R", "MALAT1"],
    }

    adata = run_basic_preprocessing(adata, min_genes=1, min_cells=1)
    outputs = run_pipeline(
        adata=adata,
        gene_sets=gene_sets,
        methods=["mean_score", "rank_score", "zscore_mean"],
        qc_columns=["n_counts", "n_genes", "sparsity"],
    )

    print("Summary")
    print(outputs["summary"])
    print("\nGene set diagnostics")
    print(outputs["gene_set_diagnostics"])
    print("\nCombined score table preview")
    print(outputs["combined_scores"].head())


if __name__ == "__main__":
    main()
