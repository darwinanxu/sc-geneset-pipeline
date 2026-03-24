from abc import ABC, abstractmethod
from anndata import AnnData
from ..types import ScoreResult

class BaseScorer(ABC):
    name: str

    @abstractmethod
    def score(self, adata: AnnData, gene_sets: dict) -> ScoreResult:
        pass