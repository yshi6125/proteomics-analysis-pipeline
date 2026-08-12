"""Differential protein-abundance analysis for the final approved analysis matrix.

Workflow: validate final approved analysis matrix -> per-protein linear models -> BH FDR ->
transparent candidate thresholds -> result tables -> volcano plot -> heatmap.

This module assumes all approved upstream preprocessing is complete. It performs
no QC, normalization, transformation, imputation for statistical testing, PCA,
normalization diagnostics, batch correction, or pathway enrichment.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist
from statsmodels.stats.multitest import multipletests


PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
FINAL_DATA_PATH = PROCESSED_DIR / "data_log2_batch_corrected.pkl"
METADATA_PATH = PROCESSED_DIR / "sample_metadata.csv"
RESULTS_DIR = PROJECT_ROOT / "results" / "differential_abundance"
FIGURES_DIR = PROJECT_ROOT / "figures" / "differential_abundance"

REQUIRED_METADATA_COLUMNS = ["sample_id", "group", "batch", "original_column"]
RESULT_COLUMNS = [
    "row_position",
    "Accession",
    "Gene.name",
    "Name",
    "PeptideCounts",
    "n_reference",
    "n_comparison",
    "mean_reference",
    "mean_comparison",
    "log2FC",
    "group_coefficient",
    "standard_error",
    "t_statistic",
    "p_value",
    "q_value",
]
UNTESTED_RESULT_COLUMNS = [
    "row_position",
    "Accession",
    "Gene.name",
    "Name",
    "PeptideCounts",
    "n_reference",
    "n_comparison",
    "reason_not_tested",
]


def _ensure_output_directories() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def load_final_analysis_data(
    data_path: Path = FINAL_DATA_PATH,
    metadata_path: Path = METADATA_PATH,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load the final approved dataframe and standardized sample metadata."""
    if not data_path.exists():
        raise FileNotFoundError(f"Final processed dataframe not found: {data_path}")
    if not metadata_path.exists():
        raise FileNotFoundError(f"Sample metadata not found: {metadata_path}")
    return pd.read_pickle(data_path), pd.read_csv(metadata_path)


def validate_final_matrix(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    reference_group: str,
    comparison_group: str,
) -> list[str]:
    """Validate the final approved analysis matrix and return ordered sample IDs."""
    missing_columns = [
        column for column in REQUIRED_METADATA_COLUMNS if column not in metadata.columns
    ]
    if missing_columns:
        raise ValueError(f"Metadata is missing required columns: {missing_columns}")
    if metadata.empty:
        raise ValueError("Metadata must contain at least one sample row.")
    required = metadata.loc[:, REQUIRED_METADATA_COLUMNS]
    if required.isna().any().any():
        columns = required.columns[required.isna().any()].tolist()
        raise ValueError(f"Required metadata columns contain missing values: {columns}")
    for column in REQUIRED_METADATA_COLUMNS:
        if metadata[column].map(
            lambda value: isinstance(value, str) and not value.strip()
        ).any():
            raise ValueError(f"Metadata column {column!r} contains empty strings.")
    sample_ids = metadata["sample_id"].tolist()
    if not all(isinstance(sample_id, str) and sample_id for sample_id in sample_ids):
        raise ValueError("Every metadata sample_id must be a non-empty string.")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Metadata sample_id values must be unique.")
    if metadata["original_column"].duplicated().any():
        raise ValueError("Metadata original_column values must be unique.")
    if df.empty:
        raise ValueError("Final analysis dataframe contains no protein rows.")
    if df.columns.duplicated().any():
        duplicates = df.columns[df.columns.duplicated(keep=False)].tolist()
        raise ValueError(f"Final dataframe contains duplicate columns: {duplicates}")
    missing_samples = [sample_id for sample_id in sample_ids if sample_id not in df.columns]
    if missing_samples:
        raise ValueError(f"Final dataframe is missing sample columns: {missing_samples}")
    abundance = df.loc[:, sample_ids]
    nonnumeric = [
        sample_id
        for sample_id in sample_ids
        if not pd.api.types.is_numeric_dtype(abundance[sample_id])
    ]
    if nonnumeric:
        raise ValueError(f"Final abundance columns must be numeric: {nonnumeric}")
    infinite_count = int(np.isinf(abundance.to_numpy(dtype=float)).sum())
    if infinite_count:
        raise ValueError(
            f"Final abundance matrix contains {infinite_count} infinite values."
        )
    if reference_group == comparison_group:
        raise ValueError("reference_group and comparison_group must be different.")
    observed_groups = set(metadata["group"].astype(str))
    missing_groups = [
        group
        for group in (reference_group, comparison_group)
        if group not in observed_groups
    ]
    if missing_groups:
        raise ValueError(f"Requested groups are absent from metadata: {missing_groups}")
    return sample_ids


def run_differential_abundance(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    reference_group: str = "NDM",
    comparison_group: str = "DM",
    min_observed_per_group: int = 4,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Fit ``abundance ~ group`` independently for every eligible protein."""
    if (
        isinstance(min_observed_per_group, bool)
        or not isinstance(min_observed_per_group, (int, np.integer))
        or min_observed_per_group < 1
    ):
        raise ValueError("min_observed_per_group must be a positive integer.")
    sample_ids = validate_final_matrix(
        df, metadata, reference_group, comparison_group
    )
    groups = metadata.set_index("sample_id").loc[sample_ids, "group"].astype(str)
    analysis_samples = groups.index[
        groups.isin([reference_group, comparison_group])
    ].tolist()
    analysis_groups = groups.loc[analysis_samples]
    reference_samples = analysis_groups.index[
        analysis_groups.eq(reference_group)
    ].tolist()
    comparison_samples = analysis_groups.index[
        analysis_groups.eq(comparison_group)
    ].tolist()

    tested_rows: list[dict[str, Any]] = []
    untested_rows: list[dict[str, Any]] = []
    annotation_columns = ["Accession", "Gene.name", "Name", "PeptideCounts"]

    for row_position, (_, row) in enumerate(df.iterrows()):
        annotations = {
            column: row[column] if column in df.columns else np.nan
            for column in annotation_columns
        }
        reference_values = pd.to_numeric(
            row.loc[reference_samples], errors="coerce"
        ).dropna()
        comparison_values = pd.to_numeric(
            row.loc[comparison_samples], errors="coerce"
        ).dropna()
        n_reference = int(len(reference_values))
        n_comparison = int(len(comparison_values))
        if (
            n_reference < min_observed_per_group
            or n_comparison < min_observed_per_group
        ):
            untested_rows.append(
                {
                    **annotations,
                    "row_position": row_position,
                    "n_reference": n_reference,
                    "n_comparison": n_comparison,
                    "reason_not_tested": (
                        "insufficient observations in one or both groups"
                    ),
                }
            )
            continue

        observed_samples = reference_values.index.tolist() + comparison_values.index.tolist()
        values = pd.concat([reference_values, comparison_values]).to_numpy(dtype=float)
        comparison_indicator = np.array(
            [0.0 if sample in reference_samples else 1.0 for sample in observed_samples]
        )
        design = sm.add_constant(comparison_indicator, has_constant="add")
        if len(values) <= design.shape[1] or np.linalg.matrix_rank(design) < design.shape[1]:
            untested_rows.append(
                {
                    **annotations,
                    "row_position": row_position,
                    "n_reference": n_reference,
                    "n_comparison": n_comparison,
                    "reason_not_tested": "group model was not estimable",
                }
            )
            continue

        try:
            model = sm.OLS(values, design).fit()
        except (ValueError, np.linalg.LinAlgError) as error:
            untested_rows.append(
                {
                    **annotations,
                    "row_position": row_position,
                    "n_reference": n_reference,
                    "n_comparison": n_comparison,
                    "reason_not_tested": f"group model failed: {error}",
                }
            )
            continue
        model_statistics = np.array(
            [model.params[1], model.bse[1], model.tvalues[1], model.pvalues[1]],
            dtype=float,
        )
        if not np.isfinite(model_statistics).all():
            untested_rows.append(
                {
                    **annotations,
                    "row_position": row_position,
                    "n_reference": n_reference,
                    "n_comparison": n_comparison,
                    "reason_not_tested": "group model returned non-finite statistics",
                }
            )
            continue
        mean_reference = float(reference_values.mean())
        mean_comparison = float(comparison_values.mean())
        log2fc = mean_comparison - mean_reference
        tested_rows.append(
            {
                "row_position": row_position,
                **annotations,
                "n_reference": n_reference,
                "n_comparison": n_comparison,
                "mean_reference": mean_reference,
                "mean_comparison": mean_comparison,
                "log2FC": log2fc,
                "group_coefficient": model_statistics[0],
                "standard_error": model_statistics[1],
                "t_statistic": model_statistics[2],
                "p_value": model_statistics[3],
            }
        )

    differential_results = pd.DataFrame(tested_rows)
    if differential_results.empty:
        differential_results = pd.DataFrame(columns=RESULT_COLUMNS)
    else:
        _, q_values, _, _ = multipletests(
            differential_results["p_value"].to_numpy(dtype=float),
            method="fdr_bh",
        )
        differential_results["q_value"] = q_values
        differential_results["_absolute_log2fc"] = differential_results[
            "log2FC"
        ].abs()
        differential_results = differential_results.sort_values(
            ["q_value", "_absolute_log2fc"],
            ascending=[True, False],
            kind="stable",
        ).drop(columns="_absolute_log2fc")
        differential_results = differential_results.loc[:, RESULT_COLUMNS]
    untested_results = pd.DataFrame(
        untested_rows,
        columns=UNTESTED_RESULT_COLUMNS,
    )
    return differential_results.reset_index(drop=True), untested_results.reset_index(drop=True)


def _threshold_definitions() -> list[dict[str, Any]]:
    return [
        {
            "level": 1,
            "label": "Level 1",
            "p_threshold": np.nan,
            "q_threshold": 0.05,
            "log2fc_threshold": 1.0,
            "interpretation": "FDR-significant and at least 2-fold change.",
        },
        {
            "level": 2,
            "label": "Level 2",
            "p_threshold": np.nan,
            "q_threshold": 0.05,
            "log2fc_threshold": 0.585,
            "interpretation": (
                "FDR-significant and at least approximately 1.5-fold change."
            ),
        },
        {
            "level": 3,
            "label": "Level 3",
            "p_threshold": np.nan,
            "q_threshold": 0.05,
            "log2fc_threshold": np.nan,
            "interpretation": "FDR-significant regardless of fold-change threshold.",
        },
        {
            "level": 4,
            "label": "Level 4 — Exploratory",
            "p_threshold": 0.05,
            "q_threshold": np.nan,
            "log2fc_threshold": 0.585,
            "interpretation": (
                "Nominally significant with at least approximately 1.5-fold change."
            ),
        },
        {
            "level": 5,
            "label": "Level 5 — Exploratory",
            "p_threshold": 0.05,
            "q_threshold": np.nan,
            "log2fc_threshold": np.nan,
            "interpretation": "Nominal significance only.",
        },
    ]


def _threshold_mask(results: pd.DataFrame, definition: dict[str, Any]) -> pd.Series:
    mask = pd.Series(True, index=results.index)
    if np.isfinite(definition["p_threshold"]):
        mask &= results["p_value"].lt(definition["p_threshold"])
    if np.isfinite(definition["q_threshold"]):
        mask &= results["q_value"].lt(definition["q_threshold"])
    if np.isfinite(definition["log2fc_threshold"]):
        mask &= results["log2FC"].abs().ge(definition["log2fc_threshold"])
    return mask


def report_candidate_thresholds(
    differential_results: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[int, pd.DataFrame]]:
    """Print and return all five predefined candidate-threshold levels."""
    definitions = _threshold_definitions()
    candidates_by_level: dict[int, pd.DataFrame] = {}
    summary_rows: list[dict[str, Any]] = []
    for definition in definitions:
        candidates = differential_results.loc[
            _threshold_mask(differential_results, definition)
        ].copy()
        level = int(definition["level"])
        candidates_by_level[level] = candidates
        summary_rows.append(
            {
                **definition,
                "number_of_candidates": len(candidates),
            }
        )
        p_text = (
            f"p < {definition['p_threshold']}"
            if np.isfinite(definition["p_threshold"])
            else "p threshold: not used"
        )
        q_text = (
            f"q < {definition['q_threshold']}"
            if np.isfinite(definition["q_threshold"])
            else "q threshold: not used"
        )
        fc_text = (
            f"|log2FC| >= {definition['log2fc_threshold']}"
            if np.isfinite(definition["log2fc_threshold"])
            else "|log2FC| threshold: not used"
        )
        print(f"\n{definition['label']}")
        print(f"Thresholds: {p_text}; {q_text}; {fc_text}")
        print(f"Interpretation: {definition['interpretation']}")
        print(f"Number of candidates: {len(candidates)}")
        if candidates.empty:
            print("Candidates: None")
        elif len(candidates) <= 20:
            names = candidates["Gene.name"].where(
                candidates["Gene.name"].notna(), candidates["Accession"]
            )
            print(f"Candidates: {names.astype(str).tolist()}")
        else:
            print("Candidates: too many to list; see the saved result table.")
    if differential_results.empty or not differential_results["q_value"].lt(0.05).any():
        print("\nNo proteins passed the FDR < 0.05 threshold.")
    return pd.DataFrame(summary_rows), candidates_by_level


def create_volcano_plot(
    differential_results: pd.DataFrame,
    reference_group: str,
    comparison_group: str,
    max_labels: int = 10,
) -> Path:
    """Create the primary q-value volcano plot with a limited set of labels."""
    if isinstance(max_labels, bool) or not isinstance(max_labels, (int, np.integer)):
        raise ValueError("max_labels must be a non-negative integer.")
    if max_labels < 0:
        raise ValueError("max_labels must be a non-negative integer.")
    _ensure_output_directories()
    output_path = FIGURES_DIR / f"volcano_{comparison_group}_vs_{reference_group}.png"
    figure, axis = plt.subplots(figsize=(10, 8))
    if differential_results.empty:
        axis.text(0.5, 0.5, "No proteins were successfully tested", ha="center")
    else:
        q_values = differential_results["q_value"].to_numpy(dtype=float)
        positive_q = q_values[q_values > 0]
        floor = float(positive_q.min() / 10) if len(positive_q) else np.finfo(float).tiny
        y_values = -np.log10(np.clip(q_values, floor, 1.0))
        absolute_fc = differential_results["log2FC"].abs()
        fdr = differential_results["q_value"].lt(0.05)
        categories = np.select(
            [fdr & absolute_fc.ge(0.585), fdr],
            ["FDR significant, |log2FC| ≥ 0.585", "FDR significant, smaller change"],
            default="Not FDR significant",
        )
        colors = {
            "FDR significant, |log2FC| ≥ 0.585": "#c62828",
            "FDR significant, smaller change": "#f9a825",
            "Not FDR significant": "#9e9e9e",
        }
        for category, color in colors.items():
            mask = categories == category
            axis.scatter(
                differential_results.loc[mask, "log2FC"],
                y_values[mask],
                s=22,
                alpha=0.75,
                color=color,
                label=category,
            )
        axis.axhline(-np.log10(0.05), color="black", linestyle="--", linewidth=1)
        for threshold, style in ((0.585, "--"), (1.0, ":")):
            axis.axvline(threshold, color="black", linestyle=style, linewidth=1)
            axis.axvline(-threshold, color="black", linestyle=style, linewidth=1)
        axis.legend(fontsize="small")
        label_candidates = differential_results.loc[fdr].copy()
        labels_are_exploratory = False
        if label_candidates.empty:
            print("No proteins passed FDR < 0.05.")
            label_candidates = differential_results.loc[
                differential_results["p_value"].lt(0.05)
            ].copy()
            labels_are_exploratory = not label_candidates.empty
        if max_labels and not label_candidates.empty:
            label_candidates["_absolute_log2fc"] = label_candidates["log2FC"].abs()
            ranking_column = "p_value" if labels_are_exploratory else "q_value"
            label_candidates = label_candidates.sort_values(
                [ranking_column, "_absolute_log2fc"],
                ascending=[True, False],
                kind="stable",
            ).head(max_labels)
            for _, result in label_candidates.iterrows():
                gene_name = result["Gene.name"]
                accession = result["Accession"]
                label = gene_name if pd.notna(gene_name) and str(gene_name).strip() else accession
                if pd.isna(label) or not str(label).strip():
                    label = f"row {int(result['row_position'])}"
                axis.annotate(
                    str(label),
                    (
                        result["log2FC"],
                        -np.log10(max(float(result["q_value"]), floor)),
                    ),
                    xytext=(4, 4),
                    textcoords="offset points",
                    fontsize=7,
                    alpha=0.9,
                )
            if labels_are_exploratory:
                axis.text(
                    0.01,
                    0.99,
                    "Protein labels are exploratory (nominal p < 0.05)",
                    transform=axis.transAxes,
                    ha="left",
                    va="top",
                    fontsize=8,
                    style="italic",
                )
    axis.set_title(f"Differential Abundance: {comparison_group} vs {reference_group}")
    axis.set_xlabel(f"log2FC ({comparison_group} - {reference_group})")
    axis.set_ylabel("-log10(q-value)")
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(figure)
    return output_path


def create_top_protein_heatmap(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    differential_results: pd.DataFrame,
    candidates_by_level: dict[int, pd.DataFrame],
    reference_group: str,
    comparison_group: str,
    top_n: int = 30,
) -> tuple[Path | None, str, list[int]]:
    """Create a visualization-only z-scored heatmap of prioritized proteins."""
    if isinstance(top_n, bool) or not isinstance(top_n, (int, np.integer)) or top_n < 1:
        raise ValueError("top_n must be a positive integer.")
    _ensure_output_directories()
    if differential_results.empty:
        return None, "No successfully tested proteins were available.", []

    selected_row_positions: list[int] = []

    def add_candidates(candidates: pd.DataFrame) -> None:
        for row_position in candidates["row_position"].tolist():
            row_position = int(row_position)
            if row_position not in selected_row_positions and 0 <= row_position < len(df):
                selected_row_positions.append(row_position)
            if len(selected_row_positions) >= top_n:
                break

    add_candidates(candidates_by_level[1])
    level_used = "Level 1: q < 0.05 and |log2FC| >= 1.0"
    if len(selected_row_positions) < top_n:
        add_candidates(candidates_by_level[2])
        level_used = "Levels 1–2: FDR-significant proteins with prioritized fold changes"
    if len(selected_row_positions) < top_n:
        add_candidates(candidates_by_level[3])
        level_used = "Levels 1–3: all FDR-significant proteins"
    if candidates_by_level[3].empty and len(selected_row_positions) < top_n:
        add_candidates(candidates_by_level[4])
        level_used = "Level 4 exploratory proteins because no FDR-significant proteins existed"
    if len(selected_row_positions) < top_n:
        add_candidates(differential_results)
        level_used += "; supplemented by strongest q-value and |log2FC| ranking"

    row_positions = selected_row_positions[:top_n]
    sample_ids = metadata["sample_id"].tolist()
    matrix = df.iloc[row_positions].loc[:, sample_ids].astype(float).reset_index(drop=True)
    row_means = matrix.mean(axis=1)
    row_standard_deviations = matrix.std(axis=1, ddof=0).replace(0, np.nan)
    z_scores = matrix.sub(row_means, axis=0).div(row_standard_deviations, axis=0)
    z_scores = z_scores.fillna(0.0)
    if len(z_scores) > 1:
        distances = pdist(z_scores.to_numpy(), metric="euclidean")
        if np.isfinite(distances).all():
            order = leaves_list(linkage(distances, method="average"))
            z_scores = z_scores.iloc[order]
            row_positions = [row_positions[int(position)] for position in order]

    protein_labels = []
    for row_position in row_positions:
        row = df.iloc[row_position]
        gene_name = row.get("Gene.name")
        accession = row.get("Accession")
        if pd.notna(gene_name) and str(gene_name).strip():
            label = str(gene_name)
        elif pd.notna(accession) and str(accession).strip():
            label = str(accession)
        else:
            label = f"row {row_position}"
        protein_labels.append(label)

    indexed_metadata = metadata.set_index("sample_id")
    sample_labels = [
        f"{sample_id} | {indexed_metadata.at[sample_id, 'group']} | "
        f"{indexed_metadata.at[sample_id, 'batch']}"
        for sample_id in sample_ids
    ]
    figure_height = max(7, 0.28 * len(z_scores) + 2)
    figure, axis = plt.subplots(figsize=(15, figure_height))
    image = axis.imshow(z_scores.to_numpy(), aspect="auto", cmap="coolwarm", vmin=-2.5, vmax=2.5)
    axis.set_xticks(range(len(sample_ids)), labels=sample_labels, rotation=90, fontsize=8)
    axis.set_yticks(range(len(protein_labels)), labels=protein_labels, fontsize=8)
    axis.set_title(
        f"Top Proteins: {comparison_group} vs {reference_group}\nSelection: {level_used}"
    )
    axis.set_xlabel("Sample | group | batch")
    axis.set_ylabel("Protein (row-wise z-score for visualization only)")
    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label("Protein-wise z-score")
    figure.tight_layout()
    output_path = FIGURES_DIR / f"top_protein_heatmap_{comparison_group}_vs_{reference_group}.png"
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(figure)
    return output_path, level_used, row_positions


def run_differential_abundance_analysis(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    reference_group: str = "NDM",
    comparison_group: str = "DM",
    min_observed_per_group: int = 4,
    top_n_heatmap: int = 30,
    max_volcano_labels: int = 10,
) -> tuple[
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    pd.DataFrame,
    dict[str, Any],
]:
    """Run the complete differential-abundance workflow."""
    _ensure_output_directories()
    sample_ids = validate_final_matrix(
        df, metadata, reference_group, comparison_group
    )
    differential_results, untested_results = run_differential_abundance(
        df,
        metadata,
        reference_group=reference_group,
        comparison_group=comparison_group,
        min_observed_per_group=min_observed_per_group,
    )
    threshold_summary, candidates_by_level = report_candidate_thresholds(
        differential_results
    )
    fdr_significant_results = differential_results.loc[
        differential_results["q_value"].lt(0.05)
    ].copy()
    nominal_results = differential_results.loc[
        differential_results["p_value"].lt(0.05)
    ].copy()
    exploratory_only_results = nominal_results.loc[
        nominal_results["q_value"].ge(0.05)
    ].copy()

    comparison_name = f"{comparison_group}_vs_{reference_group}"
    differential_results.to_csv(
        RESULTS_DIR / f"differential_abundance_{comparison_name}.csv", index=False
    )
    fdr_significant_results.to_csv(
        RESULTS_DIR / f"fdr_significant_proteins_{comparison_name}.csv", index=False
    )
    nominal_results.to_csv(
        RESULTS_DIR / f"nominal_p005_proteins_{comparison_name}.csv", index=False
    )
    exploratory_only_results.to_csv(
        RESULTS_DIR / f"exploratory_only_p005_proteins_{comparison_name}.csv",
        index=False,
    )
    untested_results.to_csv(
        RESULTS_DIR / f"untested_proteins_{comparison_name}.csv", index=False
    )
    threshold_summary.to_csv(
        RESULTS_DIR / f"candidate_threshold_summary_{comparison_name}.csv",
        index=False,
    )
    volcano_path = create_volcano_plot(
        differential_results,
        reference_group,
        comparison_group,
        max_labels=max_volcano_labels,
    )
    heatmap_path, heatmap_selection_rule, heatmap_row_positions = create_top_protein_heatmap(
        df,
        metadata,
        differential_results,
        candidates_by_level,
        reference_group,
        comparison_group,
        top_n=top_n_heatmap,
    )

    candidate_counts = {
        int(row.level): int(row.number_of_candidates)
        for row in threshold_summary.itertuples()
    }
    analysis_summary: dict[str, Any] = {
        "total_proteins_supplied": len(df),
        "proteins_successfully_tested": len(differential_results),
        "proteins_not_tested": len(untested_results),
        "number_of_samples": len(sample_ids),
        "reference_group": reference_group,
        "comparison_group": comparison_group,
        "model": "abundance ~ group",
        "fdr_method": "Benjamini-Hochberg (fdr_bh)",
        "candidate_counts": candidate_counts,
        "no_proteins_passed_fdr_005": fdr_significant_results.empty,
        "heatmap_selection_rule": heatmap_selection_rule,
        "nominal_p005_proteins": len(nominal_results),
        "exploratory_only_p005_proteins": len(exploratory_only_results),
        "heatmap_row_positions": heatmap_row_positions,
        "scaled": df.attrs.get("scaled") is True,
        "volcano_path": str(volcano_path),
        "heatmap_path": str(heatmap_path) if heatmap_path is not None else None,
    }

    print("\n=== Final differential-abundance report ===")
    print(f"Total proteins supplied: {len(df)}")
    print(f"Proteins successfully tested: {len(differential_results)}")
    print(f"Proteins not tested: {len(untested_results)}")
    print(f"Reference group: {reference_group}")
    print(f"Comparison group: {comparison_group}")
    print(f"Scaled: {df.attrs.get('scaled') is True}")
    print("FDR method: Benjamini-Hochberg (fdr_bh)")
    print(f"FDR significant (q < 0.05): {len(fdr_significant_results)}")
    print(f"Nominal p < 0.05 (includes FDR-significant proteins): {len(nominal_results)}")
    print(
        "Exploratory only (p < 0.05 and q >= 0.05): "
        f"{len(exploratory_only_results)}"
    )
    for row in threshold_summary.itertuples():
        print(
            f"{row.label}: {row.number_of_candidates} candidates — "
            f"{row.interpretation}"
        )
        if row.number_of_candidates == 0:
            print("Candidates: None")
    if fdr_significant_results.empty:
        print("No proteins passed the FDR < 0.05 threshold.")
    print(f"Heatmap selection rule: {heatmap_selection_rule}")

    return (
        differential_results,
        fdr_significant_results.reset_index(drop=True),
        nominal_results.reset_index(drop=True),
        exploratory_only_results.reset_index(drop=True),
        untested_results,
        threshold_summary,
        analysis_summary,
    )


def main() -> None:
    """Load the final approved matrix and run the default DM-versus-NDM analysis."""
    df, metadata = load_final_analysis_data()
    run_differential_abundance_analysis(df, metadata)


if __name__ == "__main__":
    main()
