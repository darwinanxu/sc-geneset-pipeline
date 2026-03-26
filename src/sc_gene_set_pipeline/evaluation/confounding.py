from typing import Iterable
import pandas as pd
from scipy.stats import spearmanr


def score_qc_correlations(
    score_df: pd.DataFrame,
    obs_df: pd.DataFrame,
    qc_columns: Iterable[str] = ("n_counts", "n_genes"),
) -> pd.DataFrame:
    required_cols = list(qc_columns)
    for col in required_cols:
        if col not in obs_df.columns:
            raise ValueError(f"Missing QC column in obs: {col}")

    rows = []
    for score_name in score_df.columns:
        for qc_col in required_cols:
            rho, pval = spearmanr(score_df[score_name], obs_df[qc_col], nan_policy="omit")
            rows.append({
                "score_name": score_name,
                "qc_metric": qc_col,
                "spearman_rho": rho,
                "p_value": pval,
            })

    return pd.DataFrame(rows)
