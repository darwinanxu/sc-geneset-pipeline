import pandas as pd

def summarize_method_performance(
    method_name: str,
    runtime_sec: float,
    auc_df: pd.DataFrame | None = None,
    qc_corr_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    summary = {
        "method": method_name,
        "runtime_sec": runtime_sec,
    }

    if auc_df is not None and not auc_df.empty:
        summary["mean_auroc"] = auc_df["auroc"].mean()
        summary["mean_auprc"] = auc_df["auprc"].mean()

    if qc_corr_df is not None and not qc_corr_df.empty:
        summary["mean_abs_qc_corr"] = qc_corr_df["spearman_rho"].abs().mean()

    return pd.DataFrame([summary])