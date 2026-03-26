import json

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from sc_gene_set_pipeline.config import default_config, load_config
from sc_gene_set_pipeline.gene_sets import (
    filter_gene_sets_to_var_names,
    gene_set_overlap_frame,
    load_gene_sets,
)
from sc_gene_set_pipeline.pipeline import run_pipeline
from sc_gene_set_pipeline.preprocessing import run_basic_preprocessing
from sc_gene_set_pipeline.scoring.registry import SCORER_REGISTRY, get_scorer


@pytest.fixture
def toy_adata():
    counts = np.array(
        [
            [5, 3, 0, 0, 1, 0],
            [4, 2, 0, 0, 2, 0],
            [0, 0, 6, 4, 0, 2],
            [0, 1, 5, 3, 0, 1],
        ],
        dtype=float,
    )
    obs = pd.DataFrame(index=["cell_1", "cell_2", "cell_3", "cell_4"])
    var = pd.DataFrame(index=["NKG7", "PRF1", "IFIT1", "ISG15", "IL7R", "MALAT1"])
    return AnnData(X=counts, obs=obs, var=var)


@pytest.fixture
def toy_gene_sets():
    return {
        "cytotoxicity": ["NKG7", "PRF1"],
        "interferon_response": ["IFIT1", "ISG15"],
        "naive_t_cell": ["IL7R", "MALAT1"],
        "missing": ["DOES_NOT_EXIST"],
    }


def test_load_gene_sets_json(tmp_path):
    gene_sets = {"program_a": ["GENE1", "GENE2"]}
    path = tmp_path / "gene_sets.json"
    path.write_text(json.dumps(gene_sets), encoding="utf-8")

    assert load_gene_sets(path) == gene_sets


def test_filter_gene_sets_to_var_names(toy_adata, toy_gene_sets):
    filtered = filter_gene_sets_to_var_names(toy_gene_sets, toy_adata.var_names)

    assert "missing" not in filtered
    assert set(filtered) == {"cytotoxicity", "interferon_response", "naive_t_cell"}


def test_gene_set_overlap_frame_reports_matches(toy_adata, toy_gene_sets):
    overlap_df = gene_set_overlap_frame(toy_gene_sets, toy_adata.var_names)

    assert set(overlap_df.columns) == {
        "gene_set",
        "n_genes_input",
        "n_genes_matched",
        "match_fraction",
    }
    missing_row = overlap_df.loc[overlap_df["gene_set"] == "missing"].iloc[0]
    assert missing_row["n_genes_matched"] == 0


@pytest.mark.parametrize("method_name", sorted(SCORER_REGISTRY))
def test_each_scorer_returns_expected_shape(toy_adata, toy_gene_sets, method_name):
    adata = run_basic_preprocessing(toy_adata, min_genes=1, min_cells=1)
    scorer = get_scorer(method_name)
    filtered = filter_gene_sets_to_var_names(toy_gene_sets, adata.var_names)

    result = scorer.score(adata, filtered)

    assert result.score_matrix.shape == (adata.n_obs, 3)
    assert list(result.score_matrix.index) == list(adata.obs_names)
    assert result.metadata["n_gene_sets"] == 3


def test_run_pipeline_returns_expected_outputs(toy_adata, toy_gene_sets):
    adata = run_basic_preprocessing(toy_adata, min_genes=1, min_cells=1)

    outputs = run_pipeline(
        adata=adata,
        gene_sets=toy_gene_sets,
        methods=["mean_score", "rank_score", "zscore_mean"],
        qc_columns=["n_counts", "n_genes", "sparsity"],
    )

    assert set(outputs) == {
        "scores",
        "qc",
        "summary",
        "gene_set_overlap",
        "filtered_gene_sets",
    }
    assert outputs["summary"].shape[0] == 3
    assert "mean_abs_qc_corr" in outputs["summary"].columns
    assert outputs["gene_set_overlap"].shape[0] == 4
    assert set(outputs["filtered_gene_sets"]) == {
        "cytotoxicity",
        "interferon_response",
        "naive_t_cell",
    }
    for qc_df in outputs["qc"].values():
        assert set(qc_df["qc_metric"]) == {"n_counts", "n_genes", "sparsity"}


def test_run_pipeline_requires_qc_columns(toy_adata, toy_gene_sets):
    with pytest.raises(ValueError, match="Missing QC column"):
        run_pipeline(
            adata=toy_adata,
            gene_sets=toy_gene_sets,
            methods=["mean_score"],
        )


def test_run_pipeline_fails_when_no_overlap(toy_adata):
    adata = run_basic_preprocessing(toy_adata, min_genes=1, min_cells=1)

    with pytest.raises(ValueError, match="No gene sets remain"):
        run_pipeline(
            adata=adata,
            gene_sets={"missing": ["X", "Y"]},
            methods=["mean_score"],
            qc_columns=["n_counts", "n_genes", "sparsity"],
        )


def test_run_pipeline_requires_at_least_one_method(toy_adata, toy_gene_sets):
    adata = run_basic_preprocessing(toy_adata, min_genes=1, min_cells=1)

    with pytest.raises(ValueError, match="At least one scoring method"):
        run_pipeline(
            adata=adata,
            gene_sets=toy_gene_sets,
            methods=[],
            qc_columns=["n_counts", "n_genes", "sparsity"],
        )


def test_basic_qc_filter_raises_when_no_cells_remain(toy_adata):
    with pytest.raises(ValueError, match="No cells remain after filtering"):
        run_basic_preprocessing(toy_adata, min_genes=100, min_cells=1)


def test_basic_qc_filter_raises_when_no_genes_remain(toy_adata):
    with pytest.raises(ValueError, match="No genes remain after filtering"):
        run_basic_preprocessing(toy_adata, min_genes=1, min_cells=100)


def test_default_config_values():
    config = default_config()

    assert config.pipeline.methods == ["mean_score", "rank_score", "zscore_mean"]
    assert config.pipeline.qc_columns == ["n_counts", "n_genes", "sparsity"]


def test_load_config_overrides_defaults(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text(
        "\n".join(
            [
                "preprocessing:",
                "  min_genes: 5",
                "pipeline:",
                "  methods:",
                "    - zscore_mean",
                "  qc_columns:",
                "    - n_counts",
                "output:",
                "  save_filtered_gene_sets: false",
            ]
        ),
        encoding="utf-8",
    )

    config = load_config(path)

    assert config.preprocessing.min_genes == 5
    assert config.pipeline.methods == ["zscore_mean"]
    assert config.pipeline.qc_columns == ["n_counts"]
    assert config.output.save_filtered_gene_sets is False
