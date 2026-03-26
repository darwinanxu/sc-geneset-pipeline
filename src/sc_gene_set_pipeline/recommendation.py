def recommend_method(n_cells: int, sparsity: float, prioritize_speed: bool = False) -> str:
    if prioritize_speed and n_cells > 100000:
        return "rank_score"
    if sparsity > 0.9:
        return "rank_score"
    if 10000 <= n_cells <= 100000:
        return "zscore_mean"
    return "mean_score"
