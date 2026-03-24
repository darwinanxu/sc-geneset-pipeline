import json
from pathlib import Path
from typing import Dict, List

GeneSets = Dict[str, List[str]]

def load_gene_sets(path: str | Path) -> GeneSets:
    path = Path(path)
    if path.suffix == ".json":
        with open(path, "r") as f:
            gene_sets = json.load(f)
    else:
        raise ValueError(f"Unsupported gene set format: {path.suffix}")

    validate_gene_sets(gene_sets)
    return gene_sets

def validate_gene_sets(gene_sets: GeneSets) -> None:
    if not isinstance(gene_sets, dict):
        raise TypeError("Gene sets must be a dictionary: {set_name: [genes]}")

    for name, genes in gene_sets.items():
        if not isinstance(name, str):
            raise TypeError("Gene set name must be a string")
        if not isinstance(genes, list) or len(genes) == 0:
            raise ValueError(f"Gene set '{name}' must be a non-empty list")
        if not all(isinstance(g, str) for g in genes):
            raise TypeError(f"All genes in '{name}' must be strings")

def filter_gene_sets_to_var_names(gene_sets: GeneSets, var_names) -> GeneSets:
    var_name_set = set(map(str, var_names))
    filtered = {}
    for set_name, genes in gene_sets.items():
        kept = [g for g in genes if g in var_name_set]
        if len(kept) > 0:
            filtered[set_name] = kept
    return filtered