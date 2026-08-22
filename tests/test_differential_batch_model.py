import json
import numpy as np
import pandas as pd
import pytest
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

from src import differential_abundance as da


def _metadata(groups, batches):
    sample_ids = [f"S{index + 1}" for index in range(len(groups))]
    return pd.DataFrame({
        "sample_id": sample_ids,
        "group": groups,
        "batch": batches,
        "original_column": sample_ids,
    })


def _data(rows, sample_ids):
    frame = pd.DataFrame({
        "Accession": [f"A{index}" for index in range(len(rows))],
        "Gene.name": [f"G{index}" for index in range(len(rows))],
        "Name": [f"Protein {index}" for index in range(len(rows))],
        "PeptideCounts": [2] * len(rows),
    })
    for position, sample_id in enumerate(sample_ids):
        frame[sample_id] = [row[position] for row in rows]
    frame.attrs["abundance_scale"] = "log2"
    return frame


def test_default_da_input_is_normalized_log2_csv():
    assert da.FINAL_DATA_PATH.name == "data_log2.csv"


def test_group_plus_batch_model_reports_group_effect_df_and_bh_values():
    metadata = _metadata(
        ["NDM", "NDM", "NDM", "NDM", "DM", "DM", "DM", "DM"],
        ["B1", "B1", "B2", "B2", "B1", "B1", "B2", "B2"],
    )
    data = _data(
        [
            [1.0, 1.2, 3.0, 3.2, 2.0, 2.2, 4.0, 4.2],
            [2.0, 2.2, 2.4, 2.6, 2.1, 2.3, 2.5, 2.7],
        ],
        metadata["sample_id"].tolist(),
    )

    results, untested = da.run_differential_abundance(
        data, metadata, min_observed_per_group=2,
    )

    assert untested.empty
    design = np.column_stack([
        np.ones(8),
        np.array([0, 0, 0, 0, 1, 1, 1, 1], dtype=float),
        np.array([0, 0, 1, 1, 0, 0, 1, 1], dtype=float),
    ])
    expected_models = [sm.OLS(np.asarray(row), design).fit() for row in [
        [1.0, 1.2, 3.0, 3.2, 2.0, 2.2, 4.0, 4.2],
        [2.0, 2.2, 2.4, 2.6, 2.1, 2.3, 2.5, 2.7],
    ]]
    by_accession = results.set_index("Accession")
    for index, model in enumerate(expected_models):
        result = by_accession.loc[f"A{index}"]
        assert result["group_coefficient"] == pytest.approx(model.params[1])
        assert result["log2FC"] == pytest.approx(model.params[1])
        assert result["standard_error"] == pytest.approx(model.bse[1])
        assert result["t_statistic"] == pytest.approx(model.tvalues[1])
        assert result["p_value"] == pytest.approx(model.pvalues[1])
        assert result["residual_degrees_of_freedom"] == 5
    expected_q = multipletests(
        [model.pvalues[1] for model in expected_models], method="fdr_bh"
    )[1]
    observed_q = by_accession.loc[["A0", "A1"], "q_value"].to_numpy()
    assert observed_q == pytest.approx(expected_q)


def test_one_batch_falls_back_to_group_only_model():
    metadata = _metadata(
        ["NDM", "NDM", "DM", "DM"], ["B1", "B1", "B1", "B1"],
    )
    data = _data([[1.0, 1.2, 2.0, 2.2]], metadata["sample_id"].tolist())
    results, untested = da.run_differential_abundance(
        data, metadata, min_observed_per_group=2,
    )
    assert untested.empty
    assert results.iloc[0]["log2FC"] == pytest.approx(1.0)
    assert results.iloc[0]["residual_degrees_of_freedom"] == 2


def test_perfect_group_batch_confounding_fails_clearly():
    metadata = _metadata(
        ["NDM", "NDM", "DM", "DM"], ["B1", "B1", "B2", "B2"],
    )
    data = _data([[1.0, 1.2, 2.0, 2.2]], metadata["sample_id"].tolist())
    with pytest.raises(ValueError, match="perfectly confounded"):
        da.run_differential_abundance(data, metadata, min_observed_per_group=2)


def test_analysis_metadata_documents_model_effect_and_cumulative_counts(
    tmp_path, monkeypatch,
):
    metadata = _metadata(
        ["NDM", "NDM", "NDM", "NDM", "DM", "DM", "DM", "DM"],
        ["B1", "B1", "B2", "B2", "B1", "B1", "B2", "B2"],
    )
    data = _data(
        [[1.0, 1.2, 3.0, 3.2, 2.0, 2.2, 4.0, 4.2]],
        metadata["sample_id"].tolist(),
    )
    monkeypatch.setattr(da, "RESULTS_DIR", tmp_path)
    monkeypatch.setattr(da, "FIGURES_DIR", tmp_path)
    monkeypatch.setattr(da, "create_volcano_plot", lambda *_args, **_kwargs: tmp_path / "primary.png")
    monkeypatch.setattr(
        da, "create_top_protein_heatmap",
        lambda *_args, **_kwargs: (None, "unchanged test rule", []),
    )
    monkeypatch.setattr(
        da, "create_exploratory_visualizations",
        lambda *_args, **_kwargs: {
            "level_4_count": 0, "level_5_count": 0, "total_count": 0,
            "volcano_path": None, "full_heatmap_path": None,
            "top_heatmap_path": None, "ranking_rule": "unchanged",
        },
    )

    results, *_ = da.run_differential_abundance_analysis(
        data, metadata, min_observed_per_group=2,
    )
    saved = json.loads(
        (tmp_path / "analysis_metadata_DM_vs_NDM.json").read_text()
    )

    assert results.iloc[0]["log2FC"] == pytest.approx(
        results.iloc[0]["group_coefficient"]
    )
    assert saved["model"] == "abundance ~ group + batch"
    assert saved["log2fc_definition"] == (
        "Batch-adjusted fitted group coefficient from abundance ~ group + batch"
    )
    assert saved["candidate_threshold_counts_are_cumulative"] is True
