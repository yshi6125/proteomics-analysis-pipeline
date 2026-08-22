import json
from pathlib import Path

import pandas as pd

from src.diagnosis import reserve_next_run_directory, run_diagnosis


def test_run_number_uses_maximum_existing_directory(tmp_path):
    (tmp_path / "diagnosis_001").mkdir()
    (tmp_path / "diagnosis_003").mkdir()
    (tmp_path / "unrelated").mkdir()

    run_id, run_dir = reserve_next_run_directory(tmp_path)

    assert run_id == "diagnosis_004"
    assert run_dir.is_dir()
    assert (tmp_path / "diagnosis_003").is_dir()


def diagnostic_result(prefix):
    return {
        "output_prefix": prefix,
        "abundance_scale": "linear",
        "distribution": (None, None, {
            "median_skewness": 1.2,
            "proportion_strongly_skewed": 0.5,
            "interpretation": "Review skewness.",
        }),
        "pca": (None, None, {
            "initial_features": 2, "features_used": 2,
            "temporarily_imputed_values": 1,
            "pc1_variance": 0.6, "pc2_variance": 0.3,
        }),
        "batch_assessment": (None, {
            "recommendation": "Repeat after log2 transformation",
            "rationale": "Linear-scale skew may affect PCA.",
            "weighted_unique_batch_r_squared_pc1_pc2_subspace": 0.1,
        }),
        "correlation": (None, {
            "median_pairwise_correlation": 0.9,
            "potential_low_correlation_samples": [],
        }),
        "hierarchical_clustering": (None, None, {
            "metric": "correlation", "linkage_method": "average",
        }),
        "outliers": (None, {
            "flagged_samples": [], "flag_rule": "At least 2 of 4 flags",
        }),
    }


def test_diagnosis_run_is_self_contained_and_updates_latest_pointer(tmp_path):
    data_path = tmp_path / "data.pkl"
    metadata_path = tmp_path / "metadata.csv"
    diagnostics_root = tmp_path / "diagnostics"
    figures_source = tmp_path / "shared_figures"
    tables_source = tmp_path / "shared_tables"
    data = pd.DataFrame({
        "gene": ["A", "B"], "S1": [1.0, 2.0], "S2": [2.0, None],
    })
    data.attrs["abundance_scale"] = "linear"
    data.to_pickle(data_path)
    pd.DataFrame({
        "sample_id": ["S1", "S2"], "group": ["DM", "NDM"],
        "batch": ["B1", "B2"], "original_column": ["raw1", "raw2"],
    }).to_csv(metadata_path, index=False)

    def runner(_data, _metadata, prefix):
        figures_source.mkdir()
        tables_source.mkdir()
        (figures_source / f"{prefix}_pca_samples.png").write_bytes(b"png")
        (tables_source / f"{prefix}_batch_assessment.csv").write_text("status\nok\n")
        return diagnostic_result(prefix)

    summary, report_path = run_diagnosis(
        comparison="DM_vs_NDM", data_path=data_path, metadata_path=metadata_path,
        diagnostics_root=diagnostics_root, runner=runner,
        figures_source=figures_source, tables_source=tables_source,
    )

    assert summary["run_id"] == "diagnosis_001"
    assert report_path.exists()
    assert (report_path.parent / "figures" / "pca_samples.png").read_bytes() == b"png"
    assert (report_path.parent / "tables" / "batch_assessment.csv").exists()
    assert not list(figures_source.glob("diagnosis_001_*"))
    report = report_path.read_text()
    assert "Missing abundance values: 1" in report
    assert "Repeat after log2 transformation" in report
    assert "figures/pca_samples.png" in report
    latest = json.loads((diagnostics_root / "latest_diagnosis.json").read_text())
    assert latest["run_id"] == "diagnosis_001"
    assert latest["report"] == "diagnosis_001/diagnosis_report_DM_vs_NDM.md"


def test_subsequent_diagnosis_never_overwrites_prior_report(tmp_path):
    first_id, first_dir = reserve_next_run_directory(tmp_path)
    first_report = first_dir / "report.md"
    first_report.write_text("first")
    second_id, second_dir = reserve_next_run_directory(tmp_path)

    assert (first_id, second_id) == ("diagnosis_001", "diagnosis_002")
    assert first_report.read_text() == "first"
    assert second_dir != first_dir


def test_diagnosis_accepts_csv_input_with_explicit_scale(tmp_path):
    data_path = tmp_path / "data.csv"
    metadata_path = tmp_path / "metadata.csv"
    pd.DataFrame({"S1": [1.0], "S2": [2.0]}).to_csv(data_path, index=False)
    pd.DataFrame({
        "sample_id": ["S1", "S2"], "group": ["A", "B"],
        "batch": ["X", "Y"], "original_column": ["raw1", "raw2"],
    }).to_csv(metadata_path, index=False)

    observed_scale = []

    def runner(data, _metadata, prefix):
        observed_scale.append(data.attrs["abundance_scale"])
        return diagnostic_result(prefix) | {"abundance_scale": data.attrs["abundance_scale"]}

    summary, _ = run_diagnosis(
        data_path=data_path, metadata_path=metadata_path,
        abundance_scale="log2", diagnostics_root=tmp_path / "diagnostics",
        runner=runner, figures_source=tmp_path / "figures",
        tables_source=tmp_path / "tables",
    )
    assert observed_scale == ["log2"]
    assert summary["preprocessing"]["log2_transformation_applied"] is True
