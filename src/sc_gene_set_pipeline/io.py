from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from anndata import AnnData, read_h5ad


def load_anndata(path: str | Path) -> AnnData:
    """
    Load an AnnData object from disk.

    Currently supports:
    - .h5ad
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"AnnData file not found: {path}")

    if path.suffix == ".h5ad":
        return read_h5ad(path)

    raise ValueError(f"Unsupported file format: {path.suffix}")


def save_dataframe(df: pd.DataFrame, path: str | Path) -> None:
    """
    Save a pandas DataFrame to CSV.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=True)


def save_json(data, path: str | Path) -> None:
    """
    Save JSON-serializable data to disk.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, sort_keys=True)


def save_anndata(adata: AnnData, path: str | Path) -> None:
    """
    Save an AnnData object to disk as .h5ad.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    if path.suffix != ".h5ad":
        raise ValueError("AnnData output path must end with .h5ad")

    adata.write(path)


def ensure_results_dirs(base_dir: str | Path) -> None:
    """
    Create standard results subdirectories.
    """
    base = Path(base_dir)
    for subdir in ["figures", "tables", "logs"]:
        (base / subdir).mkdir(parents=True, exist_ok=True)
