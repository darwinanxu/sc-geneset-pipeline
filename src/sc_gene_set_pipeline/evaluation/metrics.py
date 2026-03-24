from typing import Dict
import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score

def one_vs_rest_auc(score_vector, labels, positive_label) -> Dict[str, float]:
    y_true = (labels == positive_label).astype(int)
    if len(np.unique(y_true)) < 2:
        return {"auroc": np.nan, "auprc": np.nan}

    return {
        "auroc": roc_auc_score(y_true, score_vector),
        "auprc": average_precision_score(y_true, score_vector),
    }