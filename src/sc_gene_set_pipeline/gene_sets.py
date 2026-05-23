from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Dict, List

GeneSets = Dict[str, List[str]]


def load_gene_sets(path: str | Path) -> GeneSets:
    """
    Load gene sets from a JSON file.

    Expected format:
    {
        "set_name_1": ["GENE_A", "GENE_B"],
        "set_name_2": ["GENE_X", "GENE_Y"]
    }
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Gene set file not found: {path}")

    if path.suffix != ".json":
        raise ValueError(
            f"Unsupported gene set format: {path.suffix}. Only .json is supported for now."
        )

    with open(path, "r") as f:
        gene_sets = json.load(f)

    validate_gene_sets(gene_sets)
    return gene_sets


def validate_gene_sets(gene_sets: GeneSets) -> None:
    """
    Validate the gene set dictionary structure.
    """
    if not isinstance(gene_sets, dict):
        raise TypeError("Gene sets must be a dictionary of {set_name: [genes]}")

    if len(gene_sets) == 0:
        raise ValueError("Gene set dictionary is empty")

    for set_name, genes in gene_sets.items():
        if not isinstance(set_name, str):
            raise TypeError("Each gene set name must be a string")

        if not isinstance(genes, list):
            raise TypeError(f"Gene set '{set_name}' must map to a list of genes")

        if len(genes) == 0:
            raise ValueError(f"Gene set '{set_name}' is empty")

        if not all(isinstance(g, str) for g in genes):
            raise TypeError(f"All genes in gene set '{set_name}' must be strings")


def filter_gene_sets_to_var_names(
    gene_sets: GeneSets,
    var_names,
    min_overlap: int = 1,
) -> GeneSets:
    """
    Keep only genes that exist in adata.var_names.

    Parameters
    ----------
    gene_sets
        Dictionary of gene sets.
    var_names
        Typically adata.var_names.
    min_overlap
        Minimum number of genes required to keep a gene set.
    """
    var_name_set = set(map(str, var_names))
    filtered: GeneSets = {}

    for set_name, genes in gene_sets.items():
        kept_genes = [g for g in genes if g in var_name_set]
        if len(kept_genes) >= min_overlap:
            filtered[set_name] = kept_genes

    return filtered


def _unique_preserving_order(values: list[str]) -> list[str]:
    seen = set()
    unique_values = []
    for value in values:
        if value not in seen:
            unique_values.append(value)
            seen.add(value)
    return unique_values


def diagnose_gene_sets(
    gene_sets: GeneSets,
    var_names,
    min_overlap: int = 1,
) -> List[dict]:
    """
    Return detailed diagnostics for each gene set.

    The diagnostics include matched genes, missing genes, duplicate genes, and
    whether a gene set passes the requested minimum dataset overlap.
    """
    var_name_set = set(map(str, var_names))
    diagnostics = []

    for set_name, genes in gene_sets.items():
        unique_genes = _unique_preserving_order(genes)
        gene_counts = Counter(genes)
        duplicate_genes = sorted(gene for gene, count in gene_counts.items() if count > 1)
        matched_genes = [gene for gene in unique_genes if gene in var_name_set]
        missing_genes = [gene for gene in unique_genes if gene not in var_name_set]

        diagnostics.append(
            {
                "gene_set": set_name,
                "n_genes_input": len(genes),
                "n_unique_genes_input": len(unique_genes),
                "n_duplicate_genes": len(duplicate_genes),
                "n_genes_matched": len(matched_genes),
                "n_genes_missing": len(missing_genes),
                "match_fraction": len(matched_genes) / len(unique_genes) if unique_genes else 0.0,
                "passes_min_overlap": len(matched_genes) >= min_overlap,
                "matched_genes": ";".join(matched_genes),
                "missing_genes": ";".join(missing_genes),
                "duplicate_genes": ";".join(duplicate_genes),
            }
        )

    return diagnostics


def gene_set_diagnostics_frame(
    gene_sets: GeneSets,
    var_names,
    min_overlap: int = 1,
):
    """
    Return detailed gene set diagnostics as a DataFrame.
    """
    import pandas as pd

    return pd.DataFrame(
        diagnose_gene_sets(
            gene_sets,
            var_names,
            min_overlap=min_overlap,
        )
    )


def summarize_gene_set_overlap(gene_sets: GeneSets, var_names) -> List[dict]:
    """
    Return a simple summary of overlap between each gene set and the dataset genes.
    """
    var_name_set = set(map(str, var_names))
    summary = []

    for set_name, genes in gene_sets.items():
        overlap = [g for g in genes if g in var_name_set]
        summary.append(
            {
                "gene_set": set_name,
                "n_genes_input": len(genes),
                "n_genes_matched": len(overlap),
                "match_fraction": len(overlap) / len(genes) if len(genes) > 0 else 0.0,
            }
        )

    return summary


def gene_set_overlap_frame(gene_sets: GeneSets, var_names):
    """
    Return gene set overlap summary as a DataFrame.
    """
    import pandas as pd

    return pd.DataFrame(summarize_gene_set_overlap(gene_sets, var_names))
