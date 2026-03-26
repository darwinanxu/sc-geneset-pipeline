import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from sc_gene_set_pipeline.config import load_config
from sc_gene_set_pipeline.io import (
    ensure_results_dirs,
    load_anndata,
    save_dataframe,
    save_json,
)
from sc_gene_set_pipeline.gene_sets import load_gene_sets
from sc_gene_set_pipeline.preprocessing import basic_qc_filter, normalize_log1p, add_basic_qc_metrics
from sc_gene_set_pipeline.pipeline import run_pipeline


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True, help="Path to .h5ad file")
    parser.add_argument("--gene_sets", required=True, help="Path to gene sets JSON")
    parser.add_argument("--config", help="Optional YAML config file")
    parser.add_argument("--methods", nargs="+", help="Override methods from config")
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    config = load_config(args.config)
    methods = args.methods or config.pipeline.methods

    ensure_results_dirs(args.outdir)

    adata = load_anndata(args.data)
    adata = basic_qc_filter(
        adata,
        min_genes=config.preprocessing.min_genes,
        min_cells=config.preprocessing.min_cells,
    )
    adata = normalize_log1p(adata, target_sum=config.preprocessing.target_sum)
    adata = add_basic_qc_metrics(adata)

    gene_sets = load_gene_sets(args.gene_sets)

    outputs = run_pipeline(
        adata=adata,
        gene_sets=gene_sets,
        methods=methods,
        qc_columns=config.pipeline.qc_columns,
        min_gene_set_overlap=config.pipeline.min_gene_set_overlap,
    )

    outdir = Path(args.outdir)

    save_dataframe(outputs["summary"], outdir / "tables" / "method_summary.csv")
    save_dataframe(outputs["gene_set_overlap"], outdir / "tables" / "gene_set_overlap.csv")

    for method, df in outputs["scores"].items():
        save_dataframe(df, outdir / "tables" / f"{method}_scores.csv")

    for method, df in outputs["qc"].items():
        save_dataframe(df, outdir / "tables" / f"{method}_qc_correlations.csv")

    if config.output.save_filtered_gene_sets:
        save_json(outputs["filtered_gene_sets"], outdir / "tables" / "filtered_gene_sets.json")

    print("Pipeline completed successfully.")
    print(outputs["summary"])

if __name__ == "__main__":
    main()
