import argparse

from sc_gene_set_pipeline.io import load_anndata, save_dataframe, ensure_results_dirs
from sc_gene_set_pipeline.gene_sets import load_gene_sets
from sc_gene_set_pipeline.preprocessing import basic_qc_filter, normalize_log1p, add_basic_qc_metrics
from sc_gene_set_pipeline.pipeline import run_pipeline

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True, help="Path to .h5ad file")
    parser.add_argument("--gene_sets", required=True, help="Path to gene sets JSON")
    parser.add_argument("--methods", nargs="+", default=["mean_score", "rank_score"])
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()

def main():
    args = parse_args()

    ensure_results_dirs(args.outdir)

    adata = load_anndata(args.data)
    adata = basic_qc_filter(adata)
    adata = normalize_log1p(adata)
    adata = add_basic_qc_metrics(adata)

    gene_sets = load_gene_sets(args.gene_sets)

    outputs = run_pipeline(
        adata=adata,
        gene_sets=gene_sets,
        methods=args.methods,
    )

    save_dataframe(outputs["summary"], f"{args.outdir}/tables/method_summary.csv")

    for method, df in outputs["scores"].items():
        save_dataframe(df, f"{args.outdir}/tables/{method}_scores.csv")

    for method, df in outputs["qc"].items():
        save_dataframe(df, f"{args.outdir}/tables/{method}_qc_correlations.csv")

    print("Pipeline completed successfully.")
    print(outputs["summary"])

if __name__ == "__main__":
    main()