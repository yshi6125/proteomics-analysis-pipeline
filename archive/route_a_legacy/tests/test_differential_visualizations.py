from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
import pytest

from src import differential_abundance as da


matplotlib.use("Agg")


def differential_results():
    return pd.DataFrame({
        "row_position": [0, 1, 2, 3, 4],
        "Accession": ["A0", "A1", "A2", "A3", "A4"],
        "Gene.name": ["G0", "G1", "G2", "G3", "G4"],
        "log2FC": [1.2, 0.2, 0.7, 0.2, 0.1],
        "p_value": [0.001, 0.01, 0.01, 0.02, 0.5],
        "q_value": [0.01, 0.04, 0.2, 0.3, 0.8],
    })


def matrix_and_metadata():
    data = pd.DataFrame({
        "Accession": ["A0", "A1", "A2", "A3", "A4"],
        "Gene.name": ["G0", "G1", "G2", "G3", "G4"],
        "S1": [1, 2, 3, 4, 5], "S2": [2, 3, 4, 5, 6],
        "S3": [3, 4, 6, 8, 7], "S4": [4, 5, 7, 9, 8],
    })
    metadata = pd.DataFrame({
        "sample_id": ["S1", "S2", "S3", "S4"],
        "group": ["NDM", "NDM", "DM", "DM"],
        "batch": ["B1", "B2", "B1", "B2"],
    })
    return data, metadata


def test_exploratory_set_contains_only_exclusive_levels_4_and_5():
    results = differential_results()
    _, by_level = da.report_candidate_thresholds(results)
    exploratory = da.build_exploratory_candidate_set(by_level)

    assert exploratory["row_position"].tolist() == [2, 3]
    assert exploratory["candidate_origin"].tolist() == ["Level 4", "Level 5"]
    assert not set(exploratory["row_position"]).intersection({0, 1})


def test_exploratory_volcano_uses_nominal_p_values(tmp_path, monkeypatch):
    monkeypatch.setattr(da, "EXPLORATORY_FIGURES_DIR", tmp_path)
    monkeypatch.setattr(da.plt, "show", lambda: None)
    monkeypatch.setattr(da.plt, "close", lambda _figure: None)
    results = differential_results()
    _, by_level = da.report_candidate_thresholds(results)
    exploratory = da.build_exploratory_candidate_set(by_level)

    da.create_exploratory_volcano_plot(results, exploratory, "NDM", "DM")
    axis = da.plt.gcf().axes[0]

    background_y = np.asarray(axis.collections[0].get_offsets()[:, 1], dtype=float)
    assert np.allclose(background_y, -np.log10(results["p_value"]))
    assert axis.get_ylabel() == "-log10(p-value)"
    horizontal_lines = [
        line for line in axis.lines
        if np.allclose(np.asarray(line.get_ydata(), dtype=float), -np.log10(0.05))
    ]
    assert len(horizontal_lines) == 1
    da.plt.close("all")


def test_primary_volcano_remains_q_value_based(tmp_path, monkeypatch):
    monkeypatch.setattr(da, "FIGURES_DIR", tmp_path)
    monkeypatch.setattr(da.plt, "show", lambda: None)
    monkeypatch.setattr(da.plt, "close", lambda _figure: None)
    results = differential_results()

    da.create_volcano_plot(results, "NDM", "DM")
    axis = da.plt.gcf().axes[0]

    plotted_y = np.concatenate([
        np.asarray(collection.get_offsets()[:, 1], dtype=float)
        for collection in axis.collections
        if len(collection.get_offsets())
    ])
    assert np.allclose(
        np.sort(plotted_y), np.sort(-np.log10(results["q_value"].to_numpy()))
    )
    assert axis.get_ylabel() == "-log10(q-value)"
    da.plt.close("all")


def test_candidate_threshold_definitions_and_cumulative_counts_are_unchanged():
    definitions = da._threshold_definitions()
    assert [item["level"] for item in definitions] == [1, 2, 3, 4, 5]
    assert [item["q_threshold"] for item in definitions[:3]] == [0.05] * 3
    assert all(np.isnan(item["p_threshold"]) for item in definitions[:3])
    assert [item["log2fc_threshold"] for item in definitions[:2]] == [1.0, 0.585]
    assert np.isnan(definitions[2]["log2fc_threshold"])
    assert [item["p_threshold"] for item in definitions[3:]] == [0.05, 0.05]
    assert all(np.isnan(item["q_threshold"]) for item in definitions[3:])
    assert definitions[3]["log2fc_threshold"] == 0.585
    assert np.isnan(definitions[4]["log2fc_threshold"])
    _, by_level = da.report_candidate_thresholds(differential_results())
    assert [len(by_level[level]) for level in range(1, 6)] == [1, 1, 2, 2, 4]


def test_exploratory_figures_are_generated_from_existing_matrix(tmp_path, monkeypatch):
    monkeypatch.setattr(da, "EXPLORATORY_FIGURES_DIR", tmp_path)
    monkeypatch.setattr(da.plt, "show", lambda: None)
    results = differential_results()
    _, by_level = da.report_candidate_thresholds(results)
    data, metadata = matrix_and_metadata()

    outputs = da.create_exploratory_visualizations(
        data, metadata, results, by_level, "NDM", "DM", top_n=1,
    )

    assert outputs["level_4_count"] == 1
    assert outputs["level_5_count"] == 1
    assert outputs["total_count"] == 2
    assert Path(outputs["volcano_path"]).exists()
    assert Path(outputs["full_heatmap_path"]).exists()
    assert Path(outputs["top_heatmap_path"]).exists()
    assert set(outputs["full_heatmap_row_positions"]) == {2, 3}
    assert outputs["top_heatmap_row_positions"] == [2]
    assert data.loc[2, ["S1", "S2", "S3", "S4"]].tolist() == [3, 4, 6, 7]


def test_exploratory_ranking_is_deterministic():
    candidates = pd.DataFrame({
        "row_position": [3, 2, 1, 0],
        "candidate_origin": ["Level 5", "Level 4", "Level 4", "Level 4"],
        "p_value": [0.001, 0.02, 0.01, 0.01],
        "log2FC": [2.0, 1.0, 0.6, 0.8],
        "Gene.name": ["D", "C", "B", "A"],
        "Accession": ["D", "C", "B", "A"],
    })
    ranked = da.rank_exploratory_candidates(candidates)
    assert ranked["row_position"].tolist() == [0, 1, 2, 3]


def test_empty_exploratory_set_is_handled_gracefully(tmp_path, monkeypatch):
    monkeypatch.setattr(da, "EXPLORATORY_FIGURES_DIR", tmp_path)
    monkeypatch.setattr(da.plt, "show", lambda: None)
    results = differential_results().assign(p_value=0.5, q_value=0.8)
    _, by_level = da.report_candidate_thresholds(results)
    data, metadata = matrix_and_metadata()

    outputs = da.create_exploratory_visualizations(
        data, metadata, results, by_level, "NDM", "DM",
    )
    assert outputs["total_count"] == 0
    assert Path(outputs["volcano_path"]).exists()
    assert outputs["full_heatmap_path"] is None
    assert outputs["top_heatmap_path"] is None
