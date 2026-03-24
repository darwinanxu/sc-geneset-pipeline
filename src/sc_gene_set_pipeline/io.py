from pathlib import Path
import scanpy as sc
from anndata import AnnData

def load_anndata(path: str | Path) -> AnnData:
    path = Path(path)
    if path.suffix == ".h5ad":
        return sc.read_h5ad(path)
    raise ValueError(f"Unsupported file format: {path.suffix}")

def save_dataframe(df, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=True)

def ensure_results_dirs(base_dir: str | Path) -> None:
    base = Path(base_dir)
    for sub in ["figures", "tables", "logs"]:
        (base / sub).mkdir(parents=True, exist_ok=True)