from .mean_score import MeanExpressionScorer
from .rank_score import RankBasedScorer
from .zscore_score import ZScoreMeanScorer

SCORER_REGISTRY = {
    "mean_score": MeanExpressionScorer,
    "rank_score": RankBasedScorer,
    "zscore_mean": ZScoreMeanScorer,
}

def get_scorer(method_name: str):
    if method_name not in SCORER_REGISTRY:
        raise ValueError(f"Unknown method: {method_name}")
    return SCORER_REGISTRY[method_name]()
