# sc-gene-set-pipeline

`sc-gene-set-pipeline` is a small Python package for scoring curated gene sets in single-cell RNA-seq data, comparing scoring methods, and checking whether those scores are confounded by basic QC metrics.

The project is built around `AnnData` inputs and is meant to be usable both as:
- a command-line pipeline for `.h5ad` datasets
- a lightweight library for method benchmarking and exploratory analysis

## What It Does

Given:
- a single-cell dataset in `.h5ad` format
- a JSON file of named gene sets

the pipeline will:
1. apply basic preprocessing
2. filter gene sets to genes present in the dataset
3. score every cell with one or more gene-set scoring methods
4. compute correlations between scores and QC metrics such as `n_counts`, `n_genes`, and `sparsity`
5. write summary tables for downstream inspection

## Implemented Scoring Methods

- `mean_score`: average expression across genes in a set
- `rank_score`: average within-cell gene rank across genes in a set
- `zscore_mean`: average per-gene z-scored expression across genes in a set

## Installation

```bash
python3 -m pip install -e .
```

For development, including tests:

```bash
python3 -m pip install -e ".[dev]"
```

## Input Formats

### Expression Data

Input data must be an `AnnData` object saved as `.h5ad`.

Expected:
- `adata.X` contains raw counts or count-like expression values
- `adata.var_names` contains gene symbols or other identifiers that match the gene set file

The pipeline will add these QC columns to `adata.obs` during preprocessing:
- `n_counts`
- `n_genes`
- `sparsity`

These QC metrics are computed after filtering but before normalization/log transformation, so `n_counts` reflects the count-like input scale.

### Gene Sets

Gene sets are loaded from JSON with this structure:

```json
{
  "cytotoxicity": ["NKG7", "PRF1", "GZMB", "GNLY"],
  "interferon_response": ["IFIT1", "IFIT2", "ISG15", "MX1"]
}
```

A small example is provided at [data/gene_sets/immune_programs_example.json](/Users/darwin/projects/sc-gene-set-pipeline/data/gene_sets/immune_programs_example.json).

## Quick Start

Run the bundled synthetic example:

```bash
PYTHONPATH=src python3 examples/quickstart.py
```

Or use the CLI on a real dataset:

```bash
PYTHONPATH=src python3 scripts/run_pipeline.py \
  --data path/to/data.h5ad \
  --gene_sets data/gene_sets/immune_programs_example.json \
  --config configs/default.yaml \
  --outdir results/example_run
```

You can override methods directly from the command line:

```bash
PYTHONPATH=src python3 scripts/run_pipeline.py \
  --data path/to/data.h5ad \
  --gene_sets data/gene_sets/immune_programs_example.json \
  --methods mean_score zscore_mean \
  --outdir results/example_run
```

## Configuration

The CLI accepts an optional YAML config file. The default example is [configs/default.yaml](/Users/darwin/projects/sc-gene-set-pipeline/configs/default.yaml).

```yaml
preprocessing:
  min_genes: 1
  min_cells: 1
  target_sum: 10000

pipeline:
  methods:
    - mean_score
    - rank_score
    - zscore_mean
  qc_columns:
    - n_counts
    - n_genes
    - sparsity
  min_gene_set_overlap: 1

output:
  save_filtered_gene_sets: true
```

If no config is supplied, the package uses built-in defaults from [config.py](/Users/darwin/projects/sc-gene-set-pipeline/src/sc_gene_set_pipeline/config.py).

## Outputs

The CLI creates `figures/`, `tables/`, and `logs/` under the output directory and writes the main tables to `tables/`.

Expected outputs:
- `method_summary.csv`: one row per scoring method, including runtime and mean absolute QC correlation
- `combined_scores.csv`: long-form score table with `method`, `cell_id`, `gene_set`, and `score`
- `{method}_scores.csv`: per-cell gene set score matrix for each method
- `{method}_qc_correlations.csv`: Spearman correlations between scores and QC metrics
- `gene_set_overlap.csv`: overlap between input gene sets and dataset genes
- `gene_set_diagnostics.csv`: detailed gene set diagnostics, including matched genes, missing genes, duplicate genes, and whether each set passes the minimum overlap threshold
- `filtered_gene_sets.json`: gene sets retained after filtering to dataset genes

## Library Usage

```python
from anndata import AnnData
from sc_gene_set_pipeline.pipeline import run_pipeline
from sc_gene_set_pipeline.preprocessing import run_basic_preprocessing

adata = run_basic_preprocessing(adata, min_genes=1, min_cells=1)
outputs = run_pipeline(
    adata=adata,
    gene_sets={
        "cytotoxicity": ["NKG7", "PRF1"],
        "interferon_response": ["IFIT1", "ISG15"],
    },
    methods=["mean_score", "rank_score", "zscore_mean"],
    qc_columns=["n_counts", "n_genes", "sparsity"],
)
```

Returned keys:
- `scores`
- `combined_scores`
- `qc`
- `summary`
- `gene_set_overlap`
- `gene_set_diagnostics`
- `filtered_gene_sets`

## Project Layout

```text
sc-gene-set-pipeline/
├── configs/
├── data/
├── examples/
├── scripts/
├── src/sc_gene_set_pipeline/
└── tests/
```

Important modules:
- [pipeline.py](/Users/darwin/projects/sc-gene-set-pipeline/src/sc_gene_set_pipeline/pipeline.py): main orchestration
- [preprocessing.py](/Users/darwin/projects/sc-gene-set-pipeline/src/sc_gene_set_pipeline/preprocessing.py): default single-cell preprocessing steps
- [gene_sets.py](/Users/darwin/projects/sc-gene-set-pipeline/src/sc_gene_set_pipeline/gene_sets.py): gene set loading, validation, and overlap reporting
- [scoring/registry.py](/Users/darwin/projects/sc-gene-set-pipeline/src/sc_gene_set_pipeline/scoring/registry.py): scorer registry
- [evaluation/confounding.py](/Users/darwin/projects/sc-gene-set-pipeline/src/sc_gene_set_pipeline/evaluation/confounding.py): QC-confounding analysis

## Development

Useful commands:

```bash
make install-dev
make test
make example
make lint
```

## Current Scope

This project currently focuses on:
- simple, interpretable scoring baselines
- consistent preprocessing and reporting
- QC-confounding checks for benchmarking

It does not yet include:
- visualization/report generation
- large-scale workflow orchestration
- advanced enrichment models beyond the implemented baselines
