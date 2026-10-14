"""Versioned, deterministic QC and preprocessing-diagnosis runs."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import shutil
import sys
from typing import Any, Callable

import pandas as pd

if __package__ is None:  # Support ``python src/diagnosis.py``.
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src import normalization_diagnostics


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DIAGNOSTICS_ROOT = PROJECT_ROOT / "results" / "diagnostics"
RUN_PATTERN = re.compile(r"^diagnosis_(\d{3,})$")


def reserve_next_run_directory(root: Path = DIAGNOSTICS_ROOT) -> tuple[str, Path]:
    """Atomically reserve max(existing diagnosis numbers) + 1."""
    root.mkdir(parents=True, exist_ok=True)
    while True:
        existing = [
            int(match.group(1))
            for path in root.iterdir()
            if path.is_dir() and (match := RUN_PATTERN.fullmatch(path.name))
        ]
        run_number = max(existing, default=0) + 1
        run_id = f"diagnosis_{run_number:03d}"
        run_dir = root / run_id
        try:
            run_dir.mkdir()
        except FileExistsError:
            continue
        return run_id, run_dir


def collect_run_artifacts(
    run_id: str,
    run_dir: Path,
    *,
    figures_source: Path = normalization_diagnostics.FIGURES_DIR,
    tables_source: Path = normalization_diagnostics.RESULTS_DIR,
) -> dict[str, list[str]]:
    """Move only artifacts created with this run's unique prefix."""
    collected: dict[str, list[str]] = {"figures": [], "tables": []}
    for label, source, destination_name in (
        ("figures", figures_source, "figures"),
        ("tables", tables_source, "tables"),
    ):
        destination = run_dir / destination_name
        for artifact in sorted(source.glob(f"{run_id}_*")) if source.exists() else []:
            destination.mkdir(parents=True, exist_ok=True)
            target = destination / artifact.name.removeprefix(f"{run_id}_")
            if target.exists():
                raise FileExistsError(f"Refusing to overwrite diagnosis artifact: {target}")
            shutil.move(str(artifact), target)
            collected[label].append(str(target.relative_to(run_dir)))
    return collected


def _diagnostic_details(results: dict[str, Any], key: str, index: int = -1) -> dict[str, Any]:
    value = results[key]
    details = value[index] if isinstance(value, tuple) else value
    return details if isinstance(details, dict) else {}


def build_diagnosis_summary(
    *,
    run_id: str,
    comparison: str,
    generated_at: str,
    data_path: Path,
    metadata_path: Path,
    data: pd.DataFrame,
    metadata: pd.DataFrame,
    diagnostic_results: dict[str, Any],
    artifacts: dict[str, list[str]],
) -> dict[str, Any]:
    """Build the single structured source used by terminal and Markdown output."""
    distribution = _diagnostic_details(diagnostic_results, "distribution")
    pca = _diagnostic_details(diagnostic_results, "pca")
    batch = _diagnostic_details(diagnostic_results, "batch_assessment")
    correlation = _diagnostic_details(diagnostic_results, "correlation")
    clustering = _diagnostic_details(diagnostic_results, "hierarchical_clustering")
    outliers = _diagnostic_details(diagnostic_results, "outliers")
    abundance_columns = normalization_diagnostics.resolve_abundance_columns(data, metadata)
    abundance = normalization_diagnostics.copy_numeric_abundances(data, metadata)
    missing_values = int(abundance.isna().sum().sum())
    total_values = int(abundance.size)
    scale = str(diagnostic_results.get("abundance_scale", "unknown"))
    warnings: list[str] = []
    if scale.strip().lower() in {"linear", "linear_scaled"}:
        warnings.append(
            "Batch assessment on linear-scale data should be repeated after an "
            "approved log2 transformation before finalizing model specification."
        )
    if outliers.get("flagged_samples"):
        warnings.append("Potential outliers require scientific review; none are removed automatically.")
    return {
        "run_id": run_id,
        "comparison": comparison,
        "generated_at": generated_at,
        "inputs": {
            "data": str(data_path), "metadata": str(metadata_path),
            "rows": int(data.shape[0]), "columns": int(data.shape[1]),
            "samples": len(abundance_columns),
            "groups": metadata["group"].astype(str).value_counts().sort_index().to_dict(),
            "batches": metadata["batch"].astype(str).value_counts().sort_index().to_dict(),
        },
        "missing_values": {
            "count": missing_values,
            "percent": (100.0 * missing_values / total_values) if total_values else 0.0,
        },
        "preprocessing": {
            "abundance_scale": scale,
            "analysis_dataframe_imputed": False,
            "temporary_pca_imputation_values": pca.get("temporarily_imputed_values"),
            "log2_transformation_applied": scale.strip().lower() == "log2",
            "normalization_or_scaling_applied_by_this_run": False,
            "batch_covariate_recommendation": batch.get("recommendation"),
            "batch_covariate_rationale": batch.get("rationale"),
        },
        "distribution": distribution,
        "pca": pca,
        "batch": batch,
        "correlation": correlation,
        "clustering": clustering,
        "outliers": outliers,
        "warnings": warnings,
        "artifacts": artifacts,
    }


def _display(value: Any) -> str:
    if value is None:
        return "Not available"
    if isinstance(value, float):
        return f"{value:.6g}"
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False)
    return str(value)


def render_diagnosis_report(summary: dict[str, Any]) -> str:
    """Render a deterministic Markdown report from the structured summary."""
    inputs, prep = summary["inputs"], summary["preprocessing"]
    distribution, pca = summary["distribution"], summary["pca"]
    batch, correlation = summary["batch"], summary["correlation"]
    clustering, outliers = summary["clustering"], summary["outliers"]
    lines = [
        f"# Diagnosis Report: {summary['comparison']}", "",
        "## Run Metadata", "",
        f"- Run ID: {summary['run_id']}",
        f"- Generated at: {summary['generated_at']}",
        f"- Input data: `{inputs['data']}`", f"- Sample metadata: `{inputs['metadata']}`",
        f"- Input dimensions: {inputs['rows']} rows × {inputs['columns']} columns",
        f"- Samples: {inputs['samples']}", f"- Groups: {_display(inputs['groups'])}",
        f"- Batches: {_display(inputs['batches'])}", "",
        "## Missing Values and Imputation", "",
        f"- Missing abundance values: {summary['missing_values']['count']} "
        f"({summary['missing_values']['percent']:.3f}%)",
        f"- Analysis dataframe imputed: {_display(prep['analysis_dataframe_imputed'])}",
        "- Temporary values imputed for PCA only: "
        f"{_display(prep['temporary_pca_imputation_values'])}", "",
        "## Distribution and Preprocessing Assessment", "",
        f"- Abundance scale: {_display(prep['abundance_scale'])}",
        f"- Median sample skewness: {_display(distribution.get('median_skewness'))}",
        f"- Strongly skewed sample proportion: {_display(distribution.get('proportion_strongly_skewed'))}",
        f"- Distribution interpretation: {_display(distribution.get('interpretation'))}",
        f"- Log2 transformation already applied: {_display(prep['log2_transformation_applied'])}",
        "- Normalization/scaling applied by this diagnostic run: "
        f"{_display(prep['normalization_or_scaling_applied_by_this_run'])}", "",
        "## Batch Assessment", "",
        "- Observed QC finding: batch-associated variation assessed from PCA and metadata",
        f"- Statistical handling: {_display(prep['batch_covariate_recommendation'])}",
        f"- Rationale: {_display(prep['batch_covariate_rationale'])}",
        "- Downstream model: abundance ~ group + batch when multiple batches are present",
        f"- Unique batch variance in PC1/PC2: {_display(batch.get('weighted_unique_batch_r_squared_pc1_pc2_subspace'))}", "",
        "## PCA and QC Assessment", "",
        f"- PCA features initially available: {_display(pca.get('initial_features'))}",
        f"- PCA features used: {_display(pca.get('features_used'))}",
        f"- PC1 variance explained: {_display(pca.get('pc1_variance'))}",
        f"- PC2 variance explained: {_display(pca.get('pc2_variance'))}",
        f"- Median pairwise correlation: {_display(correlation.get('median_pairwise_correlation'))}",
        f"- Potential low-correlation samples: {_display(correlation.get('potential_low_correlation_samples'))}",
        f"- Potential outliers: {_display(outliers.get('flagged_samples'))}",
        f"- Outlier rule: {_display(outliers.get('flag_rule'))}",
        f"- Clustering settings: metric={_display(clustering.get('metric'))}, "
        f"linkage={_display(clustering.get('linkage_method'))}", "",
        "## Warnings", "",
    ]
    lines.extend(
        [f"- {warning}" for warning in summary["warnings"]]
        or ["- No additional deterministic warnings."]
    )
    lines.extend(["", "## Run Artifacts", ""])
    for label in ("figures", "tables"):
        lines.append(f"### {label.title()}")
        lines.append("")
        lines.extend([f"- [{path}]({path})" for path in summary["artifacts"][label]] or ["- None"])
        lines.append("")
    lines.extend([
        "## Diagnostic Scope", "",
        "This run is deterministic and diagnostic only. It does not transform, "
        "normalize, impute, subtract batch effects, exclude samples, or run downstream analyses.",
    ])
    return "\n".join(lines).rstrip() + "\n"


def print_terminal_summary(summary: dict[str, Any]) -> None:
    """Print a concise view of the same structured summary used for Markdown."""
    prep = summary["preprocessing"]
    print("\n=== Versioned diagnosis summary ===")
    print(f"Run ID: {summary['run_id']}")
    print(f"Comparison: {summary['comparison']}")
    print(f"Input dimensions: {summary['inputs']['rows']} x {summary['inputs']['columns']}")
    print(f"Groups: {_display(summary['inputs']['groups'])}")
    print(f"Missing abundance values: {summary['missing_values']['count']}")
    print(f"Abundance scale: {prep['abundance_scale']}")
    print(f"Batch model guidance: {_display(prep['batch_covariate_recommendation'])}")
    print(f"Potential outliers: {_display(summary['outliers'].get('flagged_samples'))}")


def run_diagnosis(
    *,
    comparison: str = "DM_vs_NDM",
    data_path: Path = normalization_diagnostics.PROCESSED_DATA_PATH,
    metadata_path: Path = normalization_diagnostics.METADATA_PATH,
    abundance_scale: str = "linear",
    diagnostics_root: Path = DIAGNOSTICS_ROOT,
    runner: Callable[[pd.DataFrame, pd.DataFrame, str], dict[str, Any]] = (
        normalization_diagnostics.run_normalization_diagnostics
    ),
    figures_source: Path = normalization_diagnostics.FIGURES_DIR,
    tables_source: Path = normalization_diagnostics.RESULTS_DIR,
) -> tuple[dict[str, Any], Path]:
    """Execute and persist one non-overwriting diagnosis run."""
    run_id, run_dir = reserve_next_run_directory(diagnostics_root)
    if data_path.suffix.lower() == ".csv":
        if not data_path.exists():
            raise FileNotFoundError(f"QC-approved dataframe not found: {data_path}")
        if not metadata_path.exists():
            raise FileNotFoundError(f"Sample metadata not found: {metadata_path}")
        data = pd.read_csv(data_path)
        metadata = pd.read_csv(metadata_path)
        data.attrs["abundance_scale"] = abundance_scale
        normalization_diagnostics.resolve_abundance_columns(data, metadata)
    else:
        data, metadata = normalization_diagnostics.load_qc_data(data_path, metadata_path)
    results = runner(data, metadata, run_id)
    artifacts = collect_run_artifacts(
        run_id, run_dir, figures_source=figures_source, tables_source=tables_source
    )
    generated_at = datetime.now(timezone.utc).isoformat()
    summary = build_diagnosis_summary(
        run_id=run_id, comparison=comparison, generated_at=generated_at,
        data_path=data_path, metadata_path=metadata_path, data=data,
        metadata=metadata, diagnostic_results=results, artifacts=artifacts,
    )
    report_path = run_dir / f"diagnosis_report_{comparison}.md"
    report_path.write_text(render_diagnosis_report(summary), encoding="utf-8")
    latest = {
        "run_id": run_id,
        "report": str(report_path.relative_to(diagnostics_root)),
        "generated_at": generated_at,
    }
    latest_path = diagnostics_root / "latest_diagnosis.json"
    temporary = latest_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(latest, indent=2) + "\n", encoding="utf-8")
    temporary.replace(latest_path)
    print_terminal_summary(summary)
    print(f"Diagnosis report: {report_path}")
    return summary, report_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", default="DM_vs_NDM")
    parser.add_argument("--data", type=Path, default=normalization_diagnostics.PROCESSED_DATA_PATH)
    parser.add_argument("--metadata", type=Path, default=normalization_diagnostics.METADATA_PATH)
    parser.add_argument(
        "--abundance-scale", choices=("linear", "linear_scaled", "log2", "unknown"),
        default="linear", help="Abundance scale assigned when --data is a CSV",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_diagnosis(
        comparison=args.comparison, data_path=args.data,
        metadata_path=args.metadata, abundance_scale=args.abundance_scale,
    )


if __name__ == "__main__":
    main()
