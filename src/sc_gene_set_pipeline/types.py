from dataclasses import dataclass, field
from typing import Dict
import pandas as pd


@dataclass
class ScoreResult:
    method: str
    score_matrix: pd.DataFrame
    runtime_sec: float
    metadata: Dict[str, object] = field(default_factory=dict)


@dataclass
class EvaluationResult:
    method: str
    metrics: Dict[str, float]
    metadata: Dict[str, object] = field(default_factory=dict)
