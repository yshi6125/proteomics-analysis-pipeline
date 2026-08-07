"""Reusable normalization diagnostics for QC-approved abundance data.

The same diagnostics can be run before and after approved preprocessing. Sample
columns are always resolved from ``metadata["sample_id"]`` in metadata row order.

The module performs diagnosis only. It does not transform, normalize, scale,
impute the analysis dataframe, or batch-correct the abundance values.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any
import warnings

import numpy as np
import pandas as pd


PROJECT_ROOT: Path = Path(__file__).resolve().parent.parent
PROCESSED_DATA_PATH: Path = PROJECT_ROOT / "data" / "processed" / "data_qc.pkl"
METADATA_PATH: Path = PROJECT_ROOT / "data" / "processed" / "sample_metadata.csv"
FIGURES_DIR: Path = PROJECT_ROOT / "figures" / "normalization_diagnostics"
RESULTS_DIR: Path = PROJECT_ROOT / "results" / "normalization_diagnostics"
REQUIRED_METADATA_COLUMNS = ["sample_id", "group", "batch", "original_column"]

def _ensure_output_directories() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def _validate_output_prefix(output_prefix: str) -> str:
    """Validate a filename-safe, non-empty diagnostic-run prefix."""
    if not isinstance(output_prefix, str) or not output_prefix.strip():
        raise ValueError("output_prefix must be a non-empty string.")
    prefix = output_prefix.strip()
    if any(character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-" for character in prefix):
        raise ValueError(
            "output_prefix may contain only letters, numbers, '.', '_', and '-'."
        )
    return prefix


def _output_path(directory: Path, output_prefix: str, filename: str) -> Path:
    """Build a prefixed output path after validating the prefix."""
    return directory / f"{_validate_output_prefix(output_prefix)}_{filename}"


def load_qc_data(
    data_path: Path = PROCESSED_DATA_PATH,
    metadata_path: Path = METADATA_PATH,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load QC-approved data and standardized sample metadata."""
    if not data_path.exists():
        raise FileNotFoundError(f"QC-approved dataframe not found: {data_path}")
    if not metadata_path.exists():
        raise FileNotFoundError(f"Sample metadata not found: {metadata_path}")
    df = pd.read_pickle(data_path)
    metadata = pd.read_csv(metadata_path)
    resolve_abundance_columns(df, metadata)
    return df, metadata


def resolve_abundance_columns(
    df: pd.DataFrame, metadata: pd.DataFrame
) -> list[str]:
    """Resolve sample columns from authoritative metadata order."""
    missing_columns = [
        column for column in REQUIRED_METADATA_COLUMNS if column not in metadata.columns
    ]
    if missing_columns:
        raise ValueError(f"Metadata is missing required columns: {missing_columns}")
    duplicated_required_columns = [
        column
        for column in REQUIRED_METADATA_COLUMNS
        if list(metadata.columns).count(column) > 1
    ]
    if duplicated_required_columns:
        raise ValueError(
            "Required metadata column names must be unique. Duplicates: "
            f"{duplicated_required_columns}"
        )
    if metadata.empty:
        raise ValueError("Metadata must contain at least one sample row.")
    required_metadata = metadata.loc[:, REQUIRED_METADATA_COLUMNS]
    columns_with_missing_values = required_metadata.columns[
        required_metadata.isna().any()
    ].tolist()
    if columns_with_missing_values:
        raise ValueError(
            "Required metadata columns contain missing values: "
            f"{columns_with_missing_values}"
        )
    columns_with_empty_values = [
        column
        for column in REQUIRED_METADATA_COLUMNS
        if metadata[column].map(
            lambda value: isinstance(value, str) and not value.strip()
        ).any()
    ]
    if columns_with_empty_values:
        raise ValueError(
            "Required metadata columns contain empty string values: "
            f"{columns_with_empty_values}"
        )
    sample_ids = metadata["sample_id"].tolist()
    if not all(isinstance(value, str) and value.strip() for value in sample_ids):
        raise ValueError("metadata['sample_id'] must contain non-empty strings.")
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("metadata['sample_id'] values must be unique.")
    if metadata["original_column"].duplicated().any():
        raise ValueError("metadata['original_column'] values must be unique.")
    missing = [sample_id for sample_id in sample_ids if sample_id not in df.columns]
    if missing:
        raise ValueError(f"Dataframe is missing metadata sample columns: {missing}")
    return sample_ids


def copy_numeric_abundances(
    df: pd.DataFrame, metadata: pd.DataFrame
) -> pd.DataFrame:
    """Return a numeric abundance copy in authoritative metadata order."""
    abundance_cols = resolve_abundance_columns(df, metadata)
    raw_abundances = df.loc[:, abundance_cols].copy()
    abundance_df = raw_abundances.apply(pd.to_numeric, errors="coerce")
    malformed = raw_abundances.notna() & abundance_df.isna()
    if malformed.any().any():
        examples = {
            column: raw_abundances.loc[malformed[column], column]
            .astype(str)
            .unique()[:5]
            .tolist()
            for column in abundance_cols
            if malformed[column].any()
        }
        raise ValueError(
            "Abundance columns contain non-missing values that cannot be parsed as "
            f"numeric. Examples by sample: {examples}"
        )
    infinite_count = int(np.isinf(abundance_df.to_numpy(dtype=float)).sum())
    if infinite_count:
        raise ValueError(
            f"Abundance columns contain {infinite_count} infinite values."
        )
    return abundance_df


def assess_distribution(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    output_prefix: str,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Inspect and compare sample distributions on the current abundance scale."""
    import matplotlib.pyplot as plt

    _ensure_output_directories()
    prefix = _validate_output_prefix(output_prefix)
    abundance_scale = str(df.attrs.get("abundance_scale", "unknown"))
    abundance_df = copy_numeric_abundances(df, metadata)

    summary = pd.DataFrame(
        {
            "observed_count": abundance_df.count(),
            "median": abundance_df.median(),
            "mean": abundance_df.mean(),
            "standard_deviation": abundance_df.std(),
            "skewness": abundance_df.skew(),
            "missing_percent": abundance_df.isna().mean() * 100,
        }
    )
    summary.index.name = "sample"
    summary.to_csv(_output_path(RESULTS_DIR, prefix, "sample_distribution_summary.csv"))

    figure, axis = plt.subplots(figsize=(16, 7))
    observed_by_sample = [abundance_df[column].dropna() for column in abundance_df.columns]
    # Matplotlib renamed ``labels`` to ``tick_labels`` in newer releases.
    # Try the current argument first, then fall back for older installations.
    try:
        axis.boxplot(
            observed_by_sample,
            tick_labels=abundance_df.columns,
            showfliers=False,
        )
    except TypeError:
        axis.boxplot(  # type: ignore[call-arg]
            observed_by_sample,
            labels=abundance_df.columns,
            showfliers=False,
        )
    axis.set_title(f"Protein-Abundance Distributions ({abundance_scale} scale)")
    axis.set_xlabel("Sample")
    axis.set_ylabel(f"Abundance ({abundance_scale} scale)")
    axis.tick_params(axis="x", labelrotation=60)
    figure.tight_layout()
    figure.savefig(
        _output_path(FIGURES_DIR, prefix, "abundance_boxplots.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(13, 7))
    plotted = 0
    for column in abundance_df.columns:
        observed = abundance_df[column].dropna()
        if observed.nunique() < 2:
            continue
        axis.hist(
            observed,
            bins=60,
            density=True,
            histtype="step",
            linewidth=1.0,
            alpha=0.8,
            label=column,
        )
        plotted += 1
    axis.set_title(f"Protein-Abundance Density Profiles ({abundance_scale} scale)")
    axis.set_xlabel(f"Abundance ({abundance_scale} scale)")
    axis.set_ylabel("Density")
    if plotted:
        axis.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize="small")
    figure.tight_layout()
    figure.savefig(
        _output_path(FIGURES_DIR, prefix, "abundance_density.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(14, 6))
    summary["median"].plot(kind="bar", ax=axis)
    axis.set_title(f"Median Abundance by Sample ({abundance_scale} scale)")
    axis.set_xlabel("Sample")
    axis.set_ylabel("Median abundance")
    axis.tick_params(axis="x", labelrotation=60)
    figure.tight_layout()
    figure.savefig(
        _output_path(FIGURES_DIR, prefix, "sample_medians.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()
    plt.close(figure)

    valid_skewness = summary["skewness"].dropna()
    median_skewness = float(valid_skewness.median())
    proportion_strong_skew = float((valid_skewness > 1.0).mean())

    valid_medians = summary["median"].dropna()
    scale_key = abundance_scale.strip().lower()
    linear_scale = scale_key in {"linear", "linear_scaled"}
    log2_scale = scale_key == "log2"
    median_cv = float("nan")
    median_ratio = float("nan")
    median_log2_range = float("nan")
    median_log2_standard_deviation = float("nan")
    if linear_scale:
        median_cv = (
            float(valid_medians.std(ddof=1) / valid_medians.mean())
            if len(valid_medians) > 1 and valid_medians.mean() != 0
            else float("nan")
        )
        median_ratio = (
            float(valid_medians.max() / valid_medians.min())
            if not valid_medians.empty and valid_medians.min() > 0
            else float("nan")
        )
    elif log2_scale:
        median_log2_range = (
            float(valid_medians.max() - valid_medians.min())
            if not valid_medians.empty
            else float("nan")
        )
        median_log2_standard_deviation = (
            float(valid_medians.std(ddof=1))
            if len(valid_medians) > 1
            else float("nan")
        )

    strongly_skewed = proportion_strong_skew >= 0.60 and median_skewness > 1.0
    if strongly_skewed and scale_key in {"linear", "linear_scaled"}:
        skew_statement = (
            "Most samples are strongly right-skewed; a log2 transformation may be "
            "worth evaluating before PCA and correlation analysis. No transformation "
            "was applied by this diagnostic function."
        )
    elif strongly_skewed:
        skew_statement = (
            f"Most samples remain strongly right-skewed on the {abundance_scale} "
            "scale; review the distributions without assuming another transformation "
            "is appropriate."
        )
    else:
        skew_statement = (
            "The distributions are not consistently strongly right-skewed across "
            "samples."
        )

    if linear_scale:
        substantial_shift = (
            np.isfinite(median_cv)
            and np.isfinite(median_ratio)
            and median_cv >= 0.15
            and median_ratio >= 1.50
        )
        shift_statement = (
            "Sample medians differ substantially, which may indicate a global "
            "intensity shift. This plot alone cannot distinguish technical bias "
            "from a true global biological effect."
            if substantial_shift
            else "Sample medians are broadly comparable, with no clear large global "
            "shift by these diagnostics."
        )
    elif log2_scale:
        # A one-unit range on a log2 scale corresponds to a twofold span between
        # the smallest and largest sample medians.
        substantial_shift = (
            np.isfinite(median_log2_range) and median_log2_range >= 1.0
        )
        shift_statement = (
            "Sample medians span at least 1 log2 unit, indicating a twofold-or-larger "
            "range that should be reviewed for technical and biological structure."
            if substantial_shift
            else "Sample medians span less than 1 log2 unit, with no clear large "
            "global shift by this diagnostic."
        )
    else:
        shift_statement = (
            "The abundance scale is unknown, so sample-median differences are "
            "reported descriptively without applying linear- or log2-scale shift "
            "thresholds."
        )

    interpretation = f"{skew_statement} {shift_statement}"

    print("\n=== Distribution assessment ===")
    print(f"Abundance scale: {abundance_scale}")
    print(summary)
    print(f"Median sample skewness: {median_skewness:.3f}")
    print(f"Samples with skewness > 1: {proportion_strong_skew:.1%}")
    if linear_scale:
        print(f"Coefficient of variation of sample medians: {median_cv:.4f}")
        print(f"Maximum-to-minimum median ratio: {median_ratio:.4f}")
    elif log2_scale:
        print(f"Range of sample medians (log2 units): {median_log2_range:.4f}")
        print(
            "Standard deviation of sample medians (log2 units): "
            f"{median_log2_standard_deviation:.4f}"
        )
    else:
        print("Scale-dependent sample-median shift thresholds were not applied.")
    print(f"Conclusion: {interpretation}")

    diagnostics = {
        "abundance_scale": abundance_scale,
        "output_prefix": prefix,
        "median_skewness": median_skewness,
        "proportion_strongly_skewed": proportion_strong_skew,
        "median_cv": median_cv,
        "median_ratio": median_ratio,
        "median_log2_range": median_log2_range,
        "median_log2_standard_deviation": median_log2_standard_deviation,
        "interpretation": interpretation,
    }
    return abundance_df, summary, diagnostics


def _prepare_multivariate_matrix(
    abundance_df: pd.DataFrame,
    *,
    minimum_observed_fraction: float = 0.50,
) -> tuple[pd.DataFrame, dict[str, int]]:
    """Prepare a sample-by-protein matrix for visualization-only diagnostics.

    Proteins observed in fewer than the requested fraction of samples are removed.
    Remaining missing values are temporarily filled with the protein-wise median.
    The original dataframe is never modified.
    """
    if not 0 < minimum_observed_fraction <= 1:
        raise ValueError("minimum_observed_fraction must be between 0 and 1.")

    sample_by_protein = abundance_df.transpose().copy()
    initial_features = sample_by_protein.shape[1]

    minimum_count = int(np.ceil(minimum_observed_fraction * sample_by_protein.shape[0]))
    retained = sample_by_protein.count(axis=0) >= minimum_count
    sample_by_protein = sample_by_protein.loc[:, retained]
    removed_for_missingness = initial_features - sample_by_protein.shape[1]

    temporarily_imputed_values = int(sample_by_protein.isna().sum().sum())
    medians = sample_by_protein.median(axis=0)
    sample_by_protein = sample_by_protein.fillna(medians)

    variable = sample_by_protein.var(axis=0, ddof=0) > 0
    sample_by_protein = sample_by_protein.loc[:, variable]
    removed_zero_variance = int((~variable).sum())

    if sample_by_protein.shape[0] < 2 or sample_by_protein.shape[1] < 2:
        raise ValueError("At least two samples and two usable proteins are required.")

    info = {
        "initial_features": int(initial_features),
        "features_used": int(sample_by_protein.shape[1]),
        "temporarily_imputed_values": temporarily_imputed_values,
        "removed_for_missingness": int(removed_for_missingness),
        "removed_zero_variance": removed_zero_variance,
    }
    return sample_by_protein, info


def assess_pca(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    output_prefix: str,
    *,
    minimum_observed_fraction: float = 0.50,
) -> tuple[pd.DataFrame, Any, dict[str, Any]]:
    """Create PCA using temporary median filling for visualization only."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from sklearn.decomposition import PCA

    _ensure_output_directories()
    prefix = _validate_output_prefix(output_prefix)
    abundance_df = copy_numeric_abundances(df, metadata)
    sample_matrix, prep_info = _prepare_multivariate_matrix(
        abundance_df,
        minimum_observed_fraction=minimum_observed_fraction,
    )

    pca = PCA(n_components=2)
    coordinates_array = pca.fit_transform(sample_matrix)
    coordinates = pd.DataFrame(
        coordinates_array,
        index=sample_matrix.index,
        columns=["PC1", "PC2"],
    )
    coordinates.index.name = "sample"

    indexed_metadata = metadata.set_index("sample_id", drop=False).reindex(
        coordinates.index
    )
    coordinates = coordinates.join(
        indexed_metadata[[c for c in ["group", "batch"] if c in indexed_metadata]]
    )
    coordinates.to_csv(_output_path(RESULTS_DIR, prefix, "pca_coordinates.csv"))

    figure, axis = plt.subplots(figsize=(10, 8))
    group_values = coordinates["group"].dropna().astype(str).unique().tolist() if "group" in coordinates else []
    batch_values = coordinates["batch"].dropna().astype(str).unique().tolist() if "batch" in coordinates else []
    color_map_object = plt.get_cmap("tab10")
    colors = color_map_object(np.linspace(0, 1, max(len(group_values), 1)))
    color_map = dict(zip(group_values, colors, strict=False))
    marker_options = ["o", "s", "^", "D", "P", "X"]
    marker_map = {
        value: marker_options[index % len(marker_options)]
        for index, value in enumerate(batch_values)
    }

    for sample, row in coordinates.iterrows():
        group = str(row["group"]) if "group" in row and pd.notna(row["group"]) else "Unknown"
        batch = str(row["batch"]) if "batch" in row and pd.notna(row["batch"]) else "Unknown"
        axis.scatter(
            row["PC1"],
            row["PC2"],
            s=85,
            color=color_map.get(group, "gray"),
            marker=marker_map.get(batch, "o"),
        )
        axis.annotate(sample, (row["PC1"], row["PC2"]), xytext=(4, 4), textcoords="offset points", fontsize=8)

    legend_handles: list[Line2D] = []
    legend_handles.extend(
        Line2D([0], [0], marker="o", color="none", markerfacecolor=color_map[value], label=f"Group: {value}", markersize=8)
        for value in group_values
    )
    legend_handles.extend(
        Line2D([0], [0], marker=marker_map[value], color="black", linestyle="none", label=f"Batch: {value}", markersize=8)
        for value in batch_values
    )
    if legend_handles:
        axis.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1))

    pc1_variance = float(pca.explained_variance_ratio_[0])
    pc2_variance = float(pca.explained_variance_ratio_[1])
    axis.set_title("PCA of Proteomics Samples")
    axis.set_xlabel(f"PC1 ({pc1_variance:.1%} variance explained)")
    axis.set_ylabel(f"PC2 ({pc2_variance:.1%} variance explained)")
    figure.tight_layout()
    figure.savefig(
        _output_path(FIGURES_DIR, prefix, "pca_samples.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()
    plt.close(figure)

    interpretation = (
        "PCA displays the dominant sample-level variation. Review whether separation "
        "follows biological group, batch, or individual samples. PCA alone does not prove that "
        "batch correction is required."
    )

    print("\n=== PCA assessment ===")
    print(
        "Temporary protein-wise median filling was used only to calculate PCA; "
        "the original abundance dataframe was not changed."
    )
    print(f"Proteins initially available: {prep_info['initial_features']}")
    print(f"Proteins used for PCA: {prep_info['features_used']}")
    print(f"Values temporarily median-imputed: {prep_info['temporarily_imputed_values']}")
    print(f"Removed for excessive missingness: {prep_info['removed_for_missingness']}")
    print(f"Removed for zero variance: {prep_info['removed_zero_variance']}")
    print(f"PC1 variance explained: {pc1_variance:.2%}")
    print(f"PC2 variance explained: {pc2_variance:.2%}")
    print(f"Conclusion: {interpretation}")

    diagnostics = {
        **prep_info,
        "output_prefix": prefix,
        "pc1_variance": pc1_variance,
        "pc2_variance": pc2_variance,
        "interpretation": interpretation,
    }
    return coordinates, pca, diagnostics


def assess_correlation_heatmap(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    output_prefix: str,
    *,
    method: str = "pearson",
    minimum_pairwise_observations: int = 10,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Calculate pairwise sample correlations and plot a heatmap."""
    import matplotlib.pyplot as plt

    _ensure_output_directories()
    prefix = _validate_output_prefix(output_prefix)
    if method not in {"pearson", "spearman"}:
        raise ValueError("method must be 'pearson' or 'spearman'.")

    abundance_df = copy_numeric_abundances(df, metadata)
    correlation = abundance_df.corr(method=method, min_periods=minimum_pairwise_observations)
    correlation.to_csv(
        _output_path(RESULTS_DIR, prefix, f"sample_correlation_{method}.csv")
    )

    figure, axis = plt.subplots(figsize=(12, 10))
    image = axis.imshow(correlation.to_numpy(), vmin=-1, vmax=1, cmap="coolwarm")
    short_labels = correlation.columns.tolist()
    axis.set_xticks(range(len(short_labels)), labels=short_labels, rotation=90, fontsize=8)
    axis.set_yticks(range(len(short_labels)), labels=short_labels, fontsize=8)
    axis.set_title(f"Sample Correlation Heatmap ({method.title()})")
    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label("Correlation")
    figure.tight_layout()
    figure.savefig(
        _output_path(
            FIGURES_DIR, prefix, f"sample_correlation_heatmap_{method}.png"
        ),
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()
    plt.close(figure)

    off_diagonal = correlation.where(~np.eye(len(correlation), dtype=bool)).stack()
    median_correlation = float(off_diagonal.median()) if not off_diagonal.empty else float("nan")
    minimum_correlation = float(off_diagonal.min()) if not off_diagonal.empty else float("nan")

    sample_median_correlations = correlation.where(
        ~np.eye(len(correlation), dtype=bool)
    ).median(axis=1)
    sample_median_correlations.name = "median_correlation_to_other_samples"
    sample_median_correlations.to_csv(
        _output_path(RESULTS_DIR, prefix, "sample_median_correlations.csv")
    )

    # Calculate median absolute deviation manually because Series.mad() was
    # removed from modern pandas and represented mean absolute deviation rather
    # than the robust median absolute deviation needed here.
    correlation_center = float(sample_median_correlations.median())
    correlation_mad = float(
        (sample_median_correlations - correlation_center).abs().median()
    )

    modified_z_threshold = -3.5

    if np.isfinite(correlation_mad) and correlation_mad > 0:
        modified_z = (
            0.6745
            * (sample_median_correlations - correlation_center)
            / correlation_mad
        )

        low_samples = modified_z[
            modified_z < modified_z_threshold
        ].index.tolist()
    else:
        modified_z = pd.Series(
            np.nan,
            index=sample_median_correlations.index,
            dtype=float,
        )
        low_samples = []

    if median_correlation >= 0.90:
        interpretation = (
            "Most samples are highly correlated overall. Review the heatmap for any "
            "single sample with consistently lower correlation than the others."
        )
    elif median_correlation >= 0.75:
        interpretation = (
            "Overall sample correlation is moderate to high, but the heatmap should be "
            "reviewed for biological-group-, batch-, or sample-specific patterns."
        )
    else:
        interpretation = (
            "Overall correlations are relatively low. This may reflect strong biology, "
            "unequal missingness, scale-dependent skewness, or technical variation."
        )

    print("\n=== Correlation assessment ===")
    print(f"Correlation method: {method}")
    print(f"Median pairwise correlation: {median_correlation:.3f}")
    print(f"Minimum pairwise correlation: {minimum_correlation:.3f}")
    print(
        "Low-correlation criterion: modified z-score "
        f"< {modified_z_threshold}"
    )
    print(
        "Potential low-correlation samples: "
        + (str(low_samples) if low_samples else "None")
    )
    print(f"Conclusion: {interpretation}")

    diagnostics = {
        "output_prefix": prefix,
        "method": method,
        "median_pairwise_correlation": median_correlation,
        "minimum_pairwise_correlation": minimum_correlation,
        "correlation_center": correlation_center,
        "correlation_mad": correlation_mad,
        "low_correlation_modified_z_threshold": modified_z_threshold,
        "potential_low_correlation_samples": low_samples,
        "interpretation": interpretation,
    }

    return correlation, diagnostics


def assess_hierarchical_clustering(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    output_prefix: str,
    *,
    minimum_observed_fraction: float = 0.50,
    metric: str = "correlation",
    linkage_method: str = "average",
) -> tuple[np.ndarray, list[str], dict[str, Any]]:
    """Perform sample-level hierarchical clustering and create a dendrogram."""
    import matplotlib.pyplot as plt
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import pdist

    _ensure_output_directories()
    prefix = _validate_output_prefix(output_prefix)
    abundance_df = copy_numeric_abundances(df, metadata)
    sample_matrix, prep_info = _prepare_multivariate_matrix(
        abundance_df,
        minimum_observed_fraction=minimum_observed_fraction,
    )

    if metric == "correlation":
        distances = pdist(sample_matrix.to_numpy(), metric="correlation")
    else:
        distances = pdist(sample_matrix.to_numpy(), metric=metric)

    if not np.isfinite(distances).all():
        raise ValueError(
            "Hierarchical-clustering distances contain non-finite values. Review "
            "constant samples, missingness, or the selected distance metric."
        )

    linkage_matrix = linkage(distances, method=linkage_method)

    # Use metadata in the displayed leaf labels so biological-group and batch patterns can
    # be assessed directly from the dendrogram.
    matched_metadata = metadata.set_index("sample_id", drop=False).reindex(
        sample_matrix.index
    )
    missing_metadata_samples = matched_metadata.index[
        matched_metadata[[c for c in ["group", "batch"] if c in matched_metadata]].isna().all(axis=1)
    ].tolist() if not matched_metadata.empty else list(sample_matrix.index)
    if missing_metadata_samples:
        warnings.warn(
            "Metadata were unavailable for some clustering samples: "
            f"{missing_metadata_samples}. Their labels will show sample ID only.",
            stacklevel=2,
        )

    labels: list[str] = []
    for sample in sample_matrix.index:
        short_name = sample
        group = (
            str(matched_metadata.at[sample, "group"])
            if "group" in matched_metadata.columns
            and pd.notna(matched_metadata.at[sample, "group"])
            else "?"
        )
        batch = (
            str(matched_metadata.at[sample, "batch"])
            if "batch" in matched_metadata.columns
            and pd.notna(matched_metadata.at[sample, "batch"])
            else "?"
        )
        labels.append(f"{short_name} [{group}, {batch}]")

    figure, axis = plt.subplots(figsize=(16, 8))
    dendrogram(linkage_matrix, labels=labels, leaf_rotation=75, ax=axis)
    axis.set_title(
        f"Hierarchical Clustering of Samples ({linkage_method} linkage, {metric} distance)"
    )
    axis.set_xlabel("Sample")
    axis.set_ylabel("Distance")
    figure.tight_layout()
    figure.savefig(
        _output_path(FIGURES_DIR, prefix, "hierarchical_clustering_dendrogram.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()
    plt.close(figure)

    interpretation = (
        "Review whether samples cluster primarily by biological group, by batch, or whether one "
        "sample branches separately from all others. Clustering is descriptive and does "
        "not by itself establish a batch effect or justify sample removal."
    )

    print("\n=== Hierarchical-clustering assessment ===")
    print(f"Distance metric: {metric}")
    print(f"Linkage method: {linkage_method}")
    print(f"Proteins used: {prep_info['features_used']}")
    print(f"Values temporarily median-imputed: {prep_info['temporarily_imputed_values']}")
    print(f"Removed for excessive missingness: {prep_info['removed_for_missingness']}")
    print(f"Removed for zero variance: {prep_info['removed_zero_variance']}")
    print(
        "These preprocessing steps were temporary for clustering visualization; "
        "the original dataframe was not modified."
    )
    print("Dendrogram label format: sample [group, batch]")
    print(f"Conclusion: {interpretation}")

    diagnostics = {
        **prep_info,
        "output_prefix": prefix,
        "metric": metric,
        "linkage_method": linkage_method,
        "interpretation": interpretation,
    }
    return linkage_matrix, labels, diagnostics


def assess_outliers(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    output_prefix: str,
    pca_coordinates: pd.DataFrame | None = None,
    correlation_matrix: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Flag samples using multiple transparent diagnostic criteria.

    A sample is flagged as a potential outlier when at least two of the following
    criteria are met:
      1. unusually high missingness,
      2. unusually low median correlation to other samples,
      3. unusually large robust distance from the PCA center,
      4. unusually high or low sample median abundance.

    Flags require scientific review and do not justify automatic sample removal.
    """
    import matplotlib.pyplot as plt

    _ensure_output_directories()
    prefix = _validate_output_prefix(output_prefix)
    abundance_df = copy_numeric_abundances(df, metadata)

    if pca_coordinates is None:
        pca_coordinates, _, _ = assess_pca(df, metadata, prefix)
    if correlation_matrix is None:
        correlation_matrix, _ = assess_correlation_heatmap(df, metadata, prefix)

    sample_names = abundance_df.columns
    missing_percent = abundance_df.isna().mean(axis=0) * 100
    sample_median_abundance = abundance_df.median(axis=0)

    corr_no_diagonal = correlation_matrix.copy()
    np.fill_diagonal(corr_no_diagonal.values, np.nan)
    median_correlation = corr_no_diagonal.median(axis=1)

    pca_xy = pca_coordinates.loc[sample_names, ["PC1", "PC2"]].astype(float)
    pca_center = pca_xy.median(axis=0)
    pca_mad = (pca_xy - pca_center).abs().median(axis=0).replace(0, np.nan)
    robust_components = 0.6745 * (pca_xy - pca_center) / pca_mad
    robust_pca_distance = np.sqrt((robust_components**2).sum(axis=1))

    def robust_high_flag(series: pd.Series, threshold: float = 3.5) -> pd.Series:
        median = series.median()
        mad = (series - median).abs().median()
        if mad == 0 or pd.isna(mad):
            return pd.Series(False, index=series.index)
        robust_z = 0.6745 * (series - median) / mad
        return robust_z > threshold

    def robust_low_flag(series: pd.Series, threshold: float = -3.5) -> pd.Series:
        median = series.median()
        mad = (series - median).abs().median()
        if mad == 0 or pd.isna(mad):
            return pd.Series(False, index=series.index)
        robust_z = 0.6745 * (series - median) / mad
        return robust_z < threshold

    def robust_two_sided_flag(
        series: pd.Series, threshold: float = 3.5
    ) -> pd.Series:
        median = series.median()
        mad = (series - median).abs().median()
        if mad == 0 or pd.isna(mad):
            return pd.Series(False, index=series.index)
        robust_z = 0.6745 * (series - median) / mad
        return robust_z.abs() > threshold

    high_missingness = robust_high_flag(missing_percent)
    low_correlation = robust_low_flag(median_correlation)
    high_pca_distance = robust_pca_distance > 3.5
    abnormal_median_abundance = robust_two_sided_flag(sample_median_abundance)

    outlier_table = pd.DataFrame(
        {
            "missing_percent": missing_percent,
            "median_correlation_to_others": median_correlation,
            "robust_pca_distance": robust_pca_distance,
            "sample_median_abundance": sample_median_abundance,
            "high_missingness_flag": high_missingness,
            "low_correlation_flag": low_correlation,
            "high_pca_distance_flag": high_pca_distance,
            "abnormal_median_abundance_flag": abnormal_median_abundance,
        }
    )
    flag_columns = [
        "high_missingness_flag",
        "low_correlation_flag",
        "high_pca_distance_flag",
        "abnormal_median_abundance_flag",
    ]
    outlier_table["number_of_flags"] = outlier_table[flag_columns].sum(axis=1)
    outlier_table["potential_outlier"] = outlier_table["number_of_flags"] >= 2
    outlier_table.index.name = "sample"
    outlier_table.to_csv(
        _output_path(RESULTS_DIR, prefix, "sample_outlier_diagnostics.csv")
    )

    figure, axis = plt.subplots(figsize=(10, 8))
    for sample, row in pca_xy.iterrows():
        is_outlier = bool(outlier_table.at[sample, "potential_outlier"])
        axis.scatter(
            row["PC1"],
            row["PC2"],
            s=130 if is_outlier else 70,
            marker="X" if is_outlier else "o",
        )
        axis.annotate(sample, (row["PC1"], row["PC2"]), xytext=(4, 4), textcoords="offset points", fontsize=8)
    axis.set_title("PCA with Potential Outlier Flags")
    axis.set_xlabel("PC1")
    axis.set_ylabel("PC2")
    figure.tight_layout()
    figure.savefig(
        _output_path(FIGURES_DIR, prefix, "pca_outlier_flags.png"),
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()
    plt.close(figure)

    flagged_samples = outlier_table.index[outlier_table["potential_outlier"]].tolist()
    if flagged_samples:
        interpretation = (
            "One or more samples were flagged by at least two independent diagnostics. "
            "Review their abundance profiles, missingness, batch, and experimental "
            "notes before "
            "considering exclusion. Do not remove samples automatically."
        )
    else:
        interpretation = (
            "No sample was flagged by at least two of the selected diagnostics. This does "
            "not prove that every sample is problem-free, but no clear multimetric outlier "
            "was identified."
        )

    print("\n=== Outlier assessment ===")
    print(outlier_table.sort_values(["potential_outlier", "number_of_flags"], ascending=False))
    print(f"Potential outliers: {flagged_samples if flagged_samples else 'None'}")
    print(f"Conclusion: {interpretation}")

    diagnostics = {
        "output_prefix": prefix,
        "flagged_samples": flagged_samples,
        "number_flagged": len(flagged_samples),
        "flag_rule": "At least 2 of 4 diagnostic flags",
        "interpretation": interpretation,
    }
    return outlier_table, diagnostics


def run_normalization_diagnostics(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    output_prefix: str,
) -> dict[str, Any]:
    """Run the complete normalization-diagnostic workflow without modifying data."""
    prefix = _validate_output_prefix(output_prefix)
    abundance_scale = str(df.attrs.get("abundance_scale", "unknown"))
    print("\nStarting normalization diagnostic workflow.")
    print(f"Output prefix: {prefix}")
    print(f"Abundance scale: {abundance_scale}")
    print("No transformation, normalization, scaling, or batch correction will be applied.\n")

    resolve_abundance_columns(df, metadata)
    distribution_results = assess_distribution(df, metadata, prefix)
    pca_results = assess_pca(df, metadata=metadata, output_prefix=prefix)
    correlation_results = assess_correlation_heatmap(
        df, metadata=metadata, output_prefix=prefix, method="pearson"
    )
    clustering_results = assess_hierarchical_clustering(
        df, metadata=metadata, output_prefix=prefix
    )
    outlier_results = assess_outliers(
        df,
        metadata=metadata,
        output_prefix=prefix,
        pca_coordinates=pca_results[0],
        correlation_matrix=correlation_results[0],
    )

    print("\n=== Workflow complete ===")
    print(f"Figures saved to: {FIGURES_DIR}")
    print(f"Result tables saved to: {RESULTS_DIR}")
    print("No data modifications were performed.")

    return {
        "output_prefix": prefix,
        "abundance_scale": abundance_scale,
        "metadata": metadata,
        "distribution": distribution_results,
        "pca": pca_results,
        "correlation": correlation_results,
        "hierarchical_clustering": clustering_results,
        "outliers": outlier_results,
    }


def main() -> None:
    """Load QC-approved artifacts and run normalization diagnostics."""
    data, metadata = load_qc_data()
    run_normalization_diagnostics(
        data, metadata=metadata, output_prefix="qc_linear"
    )


if __name__ == "__main__":
    main()
