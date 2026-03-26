from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class PreprocessingConfig:
    min_genes: int = 200
    min_cells: int = 3
    target_sum: float = 1e4


@dataclass
class PipelineConfig:
    methods: list[str] = field(default_factory=lambda: ["mean_score", "rank_score", "zscore_mean"])
    qc_columns: list[str] = field(default_factory=lambda: ["n_counts", "n_genes", "sparsity"])
    min_gene_set_overlap: int = 1


@dataclass
class OutputConfig:
    save_filtered_gene_sets: bool = True


@dataclass
class ProjectConfig:
    preprocessing: PreprocessingConfig = field(default_factory=PreprocessingConfig)
    pipeline: PipelineConfig = field(default_factory=PipelineConfig)
    output: OutputConfig = field(default_factory=OutputConfig)


def _merge_dicts(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge_dicts(merged[key], value)
        else:
            merged[key] = value
    return merged


def _to_dataclass(config_dict: dict[str, Any]) -> ProjectConfig:
    preprocessing = PreprocessingConfig(**config_dict.get("preprocessing", {}))
    pipeline = PipelineConfig(**config_dict.get("pipeline", {}))
    output = OutputConfig(**config_dict.get("output", {}))
    return ProjectConfig(
        preprocessing=preprocessing,
        pipeline=pipeline,
        output=output,
    )


def default_config() -> ProjectConfig:
    return ProjectConfig()


def load_config(path: str | Path | None = None) -> ProjectConfig:
    config = default_config()
    if path is None:
        return config

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle) or {}

    base_dict = {
        "preprocessing": vars(config.preprocessing),
        "pipeline": vars(config.pipeline),
        "output": vars(config.output),
    }
    merged = _merge_dicts(base_dict, loaded)
    return _to_dataclass(merged)
