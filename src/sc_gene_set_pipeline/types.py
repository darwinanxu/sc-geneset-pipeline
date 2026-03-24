from dataclasses import dataclass, field
from typing import Dict, List, Optional
import pandas as pd

@dataclass
class ScoreResult:
    method: str
    score_matrix: pd.DataFrame
    runtime_sec: float
    metadata: Dict = field(default_factory=dict)

@dataclass
class EvaluationResult:
    method: str
    metrics: Dict[str, float]
    metadata: Dict = field(default_factory=dict)