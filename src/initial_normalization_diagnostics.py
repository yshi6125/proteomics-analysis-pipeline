"""Proteomics data-quality diagnostics for the abundance table.

Workflow
--------
1. Check missing values
2. Inspect abundance distributions
3. PCA
4. Correlation heatmap
5. Hierarchical clustering
6. Outlier detection

The module performs diagnosis only. It does not transform, normalize, scale,
impute the analysis dataframe, or batch-correct the abundance values.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any
import re
import warnings

import numpy as np
import pandas as pd


SCRIPT_DIR: Path = Path(__file__).resolve().parent
# When stored at ``project/src/initial_normalization_diagnostics.py``, the project root is the
# parent of ``src``. The fallback also allows the module to work if it is placed
# directly in the project root.
PROJECT_ROOT: Path = SCRIPT_DIR.parent if SCRIPT_DIR.name == "src" else SCRIPT_DIR
RAW_DATA_PATH: Path = PROJECT_ROOT / "data" / "raw" / "abundantdata.xlsx"
FIGURES_DIR: Path = PROJECT_ROOT / "figures"
RESULTS_DIR: Path = PROJECT_ROOT / "results"

# Exact pattern for the 20 original abundance columns, e.g.
# Abundance.NDM_B1_1 or Abundance.DM_B2_20.
ABUNDANCE_PATTERN = re.compile(r"^Abundance\.(?:NDM|DM)_B\d+_\d+$")


def load_raw_data(file_path: Path = RAW_DATA_PATH) -> pd.DataFrame:
    """Load the raw Excel workbook without modifying its values."""
    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")
    return pd.read_excel(file_path)


def promote_first_row_to_header(data: pd.DataFrame) -> pd.DataFrame:
    """Return a copy using the first data row as column names."""
    if data.empty:
        raise ValueError("Cannot promote a header from an empty dataframe.")

    result = data.copy()
    result.columns = result.iloc[0].tolist()
    return result.iloc[1:].reset_index(drop=True)


def print_data_overview(data: pd.DataFrame) -> None:
    """Print the dataframe shape, column names, and first five rows."""
    print(f"DataFrame dimensions: {data.shape}")
    print(f"Column names: {data.columns.tolist()}")
    print("First five rows:")
    print(data.head(5))


def _ensure_output_directories() -> None:
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def resolve_abundance_columns(df: pd.DataFrame) -> list[str]:
    """Return only the original abundance columns using an exact name pattern."""
    abundance_cols = [
        column
        for column in df.columns
        if isinstance(column, str) and ABUNDANCE_PATTERN.fullmatch(column)
    ]
    if not abundance_cols:
        raise ValueError(
            "No original abundance columns were found. Expected names such as "
            "'Abundance.NDM_B1_1' or 'Abundance.DM_B2_20'."
        )
    return abundance_cols


def copy_numeric_abundances(df: pd.DataFrame) -> pd.DataFrame:
    """Return a numeric copy of the original abundance columns."""
    abundance_cols = resolve_abundance_columns(df)
    abundance_df = df.loc[:, abundance_cols].copy()
    return abundance_df.apply(pd.to_numeric, errors="coerce")


def create_sample_metadata(df: pd.DataFrame) -> pd.DataFrame:
    """Parse disease, batch, and sample number from abundance-column names."""
    rows: list[dict[str, Any]] = []

    for column in resolve_abundance_columns(df):
        sample_id = column.removeprefix("Abundance.")
        parts = sample_id.split("_")
        if len(parts) != 3:
            raise ValueError(f"Unexpected abundance-column format: {column}")

        disease, batch, sample_number = parts
        rows.append(
            {
                "sample": column,
                "sample_id": sample_id,
                "disease": disease,
                "batch": batch,
                "sample_number": int(sample_number),
            }
        )

    metadata = pd.DataFrame(rows).set_index("sample")
    return metadata


def check_missing_values(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Summarize missingness by sample and by protein and create QC plots."""
    import matplotlib.pyplot as plt

    _ensure_output_directories()
    abundance_df = copy_numeric_abundances(df)

    sample_missingness = pd.DataFrame(
        {
            "observed_count": abundance_df.count(axis=0),
            "missing_count": abundance_df.isna().sum(axis=0),
            "missing_percent": abundance_df.isna().mean(axis=0) * 100,
        }
    )
    sample_missingness.index.name = "sample"

    protein_missingness = pd.DataFrame(
        {
            "observed_count": abundance_df.count(axis=1),
            "missing_count": abundance_df.isna().sum(axis=1),
            "missing_percent": abundance_df.isna().mean(axis=1) * 100,
        },
        index=df.index,
    )
    if "Accession" in df.columns:
        protein_missingness.insert(0, "Accession", df["Accession"].values)

    sample_missingness.to_csv(RESULTS_DIR / "sample_missingness.csv")
    protein_missingness.to_csv(RESULTS_DIR / "protein_missingness.csv", index=False)

    figure, axis = plt.subplots(figsize=(14, 6))
    sample_missingness["missing_percent"].plot(kind="bar", ax=axis)
    axis.set_title("Missing Abundance Values by Sample")
    axis.set_xlabel("Sample")
    axis.set_ylabel("Missing values (%)")
    axis.tick_params(axis="x", labelrotation=60)
    figure.tight_layout()
    figure.savefig(FIGURES_DIR / "missingness_by_sample.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(10, 6))
    axis.hist(protein_missingness["missing_percent"], bins=21, edgecolor="black")
    axis.set_title("Distribution of Missingness Across Proteins")
    axis.set_xlabel("Missing values per protein (%)")
    axis.set_ylabel("Number of proteins")
    figure.tight_layout()
    figure.savefig(FIGURES_DIR / "missingness_by_protein.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(figure)

    max_sample_missing = float(sample_missingness["missing_percent"].max())
    median_sample_missing = float(sample_missingness["missing_percent"].median())
    sample_missing_spread = float(
        sample_missingness["missing_percent"].max()
        - sample_missingness["missing_percent"].min()
    )
    proteins_over_half_missing = int((protein_missingness["missing_percent"] > 50).sum())
    proportion_over_half_missing = float(
        (protein_missingness["missing_percent"] > 50).mean()
    )

    if max_sample_missing >= 40 or sample_missing_spread >= 25:
        interpretation = (
            "Missingness is high or strongly unequal across samples. Samples with much "
            "higher missingness may have lower detection depth or technical problems and "
            "should be reviewed before downstream analysis."
        )
    elif max_sample_missing >= 20 or sample_missing_spread >= 10:
        interpretation = (
            "Missingness is moderate or somewhat unequal across samples. Review the "
            "sample-level plot and determine whether missingness follows disease group, "
            "batch, or individual samples."
        )
    else:
        interpretation = (
            "Sample-level missingness appears relatively low and comparable. Missingness "
            "alone does not identify a clear problematic sample."
        )

    print("\n=== Missing-value assessment ===")
    print(sample_missingness.sort_values("missing_percent", ascending=False))
    print(f"Median sample missingness: {median_sample_missing:.2f}%")
    print(f"Maximum sample missingness: {max_sample_missing:.2f}%")
    print(f"Range of sample missingness: {sample_missing_spread:.2f} percentage points")
    print(
        "Proteins with >50% missing values: "
        f"{proteins_over_half_missing} ({proportion_over_half_missing:.1%})"
    )
    print(f"Conclusion: {interpretation}")

    diagnostics = {
        "median_sample_missing_percent": median_sample_missing,
        "max_sample_missing_percent": max_sample_missing,
        "sample_missingness_spread": sample_missing_spread,
        "proteins_over_half_missing": proteins_over_half_missing,
        "proportion_proteins_over_half_missing": proportion_over_half_missing,
        "interpretation": interpretation,
    }
    return sample_missingness, protein_missingness, diagnostics


def assess_distribution(
    df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Inspect raw abundance distributions without transforming or scaling data."""
    import matplotlib.pyplot as plt

    _ensure_output_directories()
    abundance_df = copy_numeric_abundances(df)

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
    summary.to_csv(RESULTS_DIR / "sample_distribution_summary.csv")

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
    axis.set_title("Raw Protein-Abundance Distributions")
    axis.set_xlabel("Sample")
    axis.set_ylabel("Raw abundance")
    axis.tick_params(axis="x", labelrotation=60)
    figure.tight_layout()
    figure.savefig(FIGURES_DIR / "raw_abundance_boxplots.png", dpi=300, bbox_inches="tight")
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
    axis.set_title("Raw Protein-Abundance Density Profiles")
    axis.set_xlabel("Raw abundance")
    axis.set_ylabel("Density")
    if plotted:
        axis.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize="small")
    figure.tight_layout()
    figure.savefig(FIGURES_DIR / "raw_abundance_density.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(14, 6))
    summary["median"].plot(kind="bar", ax=axis)
    axis.set_title("Median Raw Abundance by Sample")
    axis.set_xlabel("Sample")
    axis.set_ylabel("Median abundance")
    axis.tick_params(axis="x", labelrotation=60)
    figure.tight_layout()
    figure.savefig(FIGURES_DIR / "sample_medians_raw.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(figure)

    valid_skewness = summary["skewness"].dropna()
    median_skewness = float(valid_skewness.median())
    proportion_strong_skew = float((valid_skewness > 1.0).mean())

    valid_medians = summary["median"].dropna()
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

    skew_statement = (
        "Most samples are strongly right-skewed; a log2 transformation may be worth "
        "evaluating before PCA and correlation analysis. No transformation was applied."
        if proportion_strong_skew >= 0.60 and median_skewness > 1.0
        else "The raw distributions are not consistently strongly right-skewed across samples."
    )

    shift_statement = (
        "Sample medians differ substantially, which may indicate a global intensity shift. "
        "This plot alone cannot distinguish technical bias from a true global biological effect."
        if np.isfinite(median_cv)
        and np.isfinite(median_ratio)
        and median_cv >= 0.15
        and median_ratio >= 1.50
        else "Sample medians are broadly comparable, with no clear large global shift by these diagnostics."
    )

    interpretation = f"{skew_statement} {shift_statement}"

    print("\n=== Distribution assessment ===")
    print(summary)
    print(f"Median sample skewness: {median_skewness:.3f}")
    print(f"Samples with skewness > 1: {proportion_strong_skew:.1%}")
    print(f"Coefficient of variation of sample medians: {median_cv:.4f}")
    print(f"Maximum-to-minimum median ratio: {median_ratio:.4f}")
    print(f"Conclusion: {interpretation}")

    diagnostics = {
        "median_skewness": median_skewness,
        "proportion_strongly_skewed": proportion_strong_skew,
        "median_cv": median_cv,
        "median_ratio": median_ratio,
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
        "removed_for_missingness": int(removed_for_missingness),
        "removed_zero_variance": removed_zero_variance,
    }
    return sample_by_protein, info


def assess_pca(
    df: pd.DataFrame,
    metadata: pd.DataFrame | None = None,
    *,
    minimum_observed_fraction: float = 0.50,
) -> tuple[pd.DataFrame, Any, dict[str, Any]]:
    """Create PCA using temporary median filling for visualization only."""
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from sklearn.decomposition import PCA

    _ensure_output_directories()
    abundance_df = copy_numeric_abundances(df)
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

    if metadata is None:
        metadata = create_sample_metadata(df)
    metadata = metadata.reindex(coordinates.index)
    coordinates = coordinates.join(metadata[[c for c in ["disease", "batch"] if c in metadata]])
    coordinates.to_csv(RESULTS_DIR / "pca_coordinates.csv")

    figure, axis = plt.subplots(figsize=(10, 8))
    disease_values = coordinates["disease"].dropna().astype(str).unique().tolist() if "disease" in coordinates else []
    batch_values = coordinates["batch"].dropna().astype(str).unique().tolist() if "batch" in coordinates else []
    color_map_object = plt.get_cmap("tab10")
    colors = color_map_object(np.linspace(0, 1, max(len(disease_values), 1)))
    color_map = dict(zip(disease_values, colors, strict=False))
    marker_options = ["o", "s", "^", "D", "P", "X"]
    marker_map = {
        value: marker_options[index % len(marker_options)]
        for index, value in enumerate(batch_values)
    }

    for sample, row in coordinates.iterrows():
        disease = str(row["disease"]) if "disease" in row and pd.notna(row["disease"]) else "Unknown"
        batch = str(row["batch"]) if "batch" in row and pd.notna(row["batch"]) else "Unknown"
        axis.scatter(
            row["PC1"],
            row["PC2"],
            s=85,
            color=color_map.get(disease, "gray"),
            marker=marker_map.get(batch, "o"),
        )
        axis.annotate(sample.removeprefix("Abundance."), (row["PC1"], row["PC2"]), xytext=(4, 4), textcoords="offset points", fontsize=8)

    legend_handles: list[Line2D] = []
    legend_handles.extend(
        Line2D([0], [0], marker="o", color="none", markerfacecolor=color_map[value], label=f"Disease: {value}", markersize=8)
        for value in disease_values
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
    figure.savefig(FIGURES_DIR / "pca_samples.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(figure)

    interpretation = (
        "PCA displays the dominant sample-level variation. Review whether separation "
        "follows disease, batch, or individual samples. PCA alone does not prove that "
        "batch correction is required."
    )

    print("\n=== PCA assessment ===")
    print(
        "Temporary protein-wise median filling was used only to calculate PCA; "
        "the original abundance dataframe was not changed."
    )
    print(f"Proteins initially available: {prep_info['initial_features']}")
    print(f"Proteins used for PCA: {prep_info['features_used']}")
    print(f"Removed for excessive missingness: {prep_info['removed_for_missingness']}")
    print(f"Removed for zero variance: {prep_info['removed_zero_variance']}")
    print(f"PC1 variance explained: {pc1_variance:.2%}")
    print(f"PC2 variance explained: {pc2_variance:.2%}")
    print(f"Conclusion: {interpretation}")

    diagnostics = {
        **prep_info,
        "pc1_variance": pc1_variance,
        "pc2_variance": pc2_variance,
        "interpretation": interpretation,
    }
    return coordinates, pca, diagnostics


def assess_correlation_heatmap(
    df: pd.DataFrame,
    *,
    method: str = "pearson",
    minimum_pairwise_observations: int = 10,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Calculate pairwise sample correlations and plot a heatmap."""
    import matplotlib.pyplot as plt

    _ensure_output_directories()
    if method not in {"pearson", "spearman"}:
        raise ValueError("method must be 'pearson' or 'spearman'.")

    abundance_df = copy_numeric_abundances(df)
    correlation = abundance_df.corr(method=method, min_periods=minimum_pairwise_observations)
    correlation.to_csv(RESULTS_DIR / f"sample_correlation_{method}.csv")

    figure, axis = plt.subplots(figsize=(12, 10))
    image = axis.imshow(correlation.to_numpy(), vmin=-1, vmax=1, cmap="coolwarm")
    short_labels = [column.removeprefix("Abundance.") for column in correlation.columns]
    axis.set_xticks(range(len(short_labels)), labels=short_labels, rotation=90, fontsize=8)
    axis.set_yticks(range(len(short_labels)), labels=short_labels, fontsize=8)
    axis.set_title(f"Sample Correlation Heatmap ({method.title()})")
    colorbar = figure.colorbar(image, ax=axis)
    colorbar.set_label("Correlation")
    figure.tight_layout()
    figure.savefig(FIGURES_DIR / f"sample_correlation_heatmap_{method}.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(figure)

    off_diagonal = correlation.where(~np.eye(len(correlation), dtype=bool)).stack()
    median_correlation = float(off_diagonal.median()) if not off_diagonal.empty else float("nan")
    minimum_correlation = float(off_diagonal.min()) if not off_diagonal.empty else float("nan")

    sample_median_correlations = correlation.where(
        ~np.eye(len(correlation), dtype=bool)
    ).median(axis=1)
    sample_median_correlations.name = "median_correlation_to_other_samples"
    sample_median_correlations.to_csv(RESULTS_DIR / "sample_median_correlations.csv")

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
            "reviewed for disease-, batch-, or sample-specific patterns."
        )
    else:
        interpretation = (
            "Overall correlations are relatively low. This may reflect strong biology, "
            "unequal missingness, raw-scale skewness, or technical variation."
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
    metadata: pd.DataFrame | None = None,
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
    abundance_df = copy_numeric_abundances(df)
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

    # Use metadata in the displayed leaf labels so disease and batch patterns can
    # be assessed directly from the dendrogram.
    if metadata is None:
        metadata = create_sample_metadata(df)
    matched_metadata = metadata.reindex(sample_matrix.index)
    missing_metadata_samples = matched_metadata.index[
        matched_metadata[[c for c in ["disease", "batch"] if c in matched_metadata]].isna().all(axis=1)
    ].tolist() if not matched_metadata.empty else list(sample_matrix.index)
    if missing_metadata_samples:
        warnings.warn(
            "Metadata were unavailable for some clustering samples: "
            f"{missing_metadata_samples}. Their labels will show sample ID only.",
            stacklevel=2,
        )

    labels: list[str] = []
    for sample in sample_matrix.index:
        short_name = sample.removeprefix("Abundance.")
        disease = (
            str(matched_metadata.at[sample, "disease"])
            if "disease" in matched_metadata.columns
            and pd.notna(matched_metadata.at[sample, "disease"])
            else "?"
        )
        batch = (
            str(matched_metadata.at[sample, "batch"])
            if "batch" in matched_metadata.columns
            and pd.notna(matched_metadata.at[sample, "batch"])
            else "?"
        )
        labels.append(f"{short_name} [{disease}, {batch}]")

    figure, axis = plt.subplots(figsize=(16, 8))
    dendrogram(linkage_matrix, labels=labels, leaf_rotation=75, ax=axis)
    axis.set_title(
        f"Hierarchical Clustering of Samples ({linkage_method} linkage, {metric} distance)"
    )
    axis.set_xlabel("Sample")
    axis.set_ylabel("Distance")
    figure.tight_layout()
    figure.savefig(FIGURES_DIR / "hierarchical_clustering_dendrogram.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(figure)

    interpretation = (
        "Review whether samples cluster primarily by disease, by batch, or whether one "
        "sample branches separately from all others. Clustering is descriptive and does "
        "not by itself establish a batch effect or justify sample removal."
    )

    print("\n=== Hierarchical-clustering assessment ===")
    print(f"Distance metric: {metric}")
    print(f"Linkage method: {linkage_method}")
    print(f"Proteins used: {prep_info['features_used']}")
    print("Dendrogram label format: sample [disease, batch]")
    print(f"Conclusion: {interpretation}")

    diagnostics = {
        **prep_info,
        "metric": metric,
        "linkage_method": linkage_method,
        "interpretation": interpretation,
    }
    return linkage_matrix, labels, diagnostics


def assess_outliers(
    df: pd.DataFrame,
    pca_coordinates: pd.DataFrame | None = None,
    correlation_matrix: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Flag samples using multiple transparent diagnostic criteria.

    A sample is flagged as a potential outlier when at least two of the following
    criteria are met:
      1. unusually high missingness,
      2. unusually low median correlation to other samples,
      3. unusually large robust distance from the PCA center.

    Flags require scientific review and do not justify automatic sample removal.
    """
    import matplotlib.pyplot as plt

    _ensure_output_directories()
    abundance_df = copy_numeric_abundances(df)

    if pca_coordinates is None:
        pca_coordinates, _, _ = assess_pca(df)
    if correlation_matrix is None:
        correlation_matrix, _ = assess_correlation_heatmap(df)

    sample_names = abundance_df.columns
    missing_percent = abundance_df.isna().mean(axis=0) * 100

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

    high_missingness = robust_high_flag(missing_percent)
    low_correlation = robust_low_flag(median_correlation)
    high_pca_distance = robust_pca_distance > 3.5

    outlier_table = pd.DataFrame(
        {
            "missing_percent": missing_percent,
            "median_correlation_to_others": median_correlation,
            "robust_pca_distance": robust_pca_distance,
            "high_missingness_flag": high_missingness,
            "low_correlation_flag": low_correlation,
            "high_pca_distance_flag": high_pca_distance,
        }
    )
    flag_columns = [
        "high_missingness_flag",
        "low_correlation_flag",
        "high_pca_distance_flag",
    ]
    outlier_table["number_of_flags"] = outlier_table[flag_columns].sum(axis=1)
    outlier_table["potential_outlier"] = outlier_table["number_of_flags"] >= 2
    outlier_table.index.name = "sample"
    outlier_table.to_csv(RESULTS_DIR / "sample_outlier_diagnostics.csv")

    figure, axis = plt.subplots(figsize=(10, 8))
    for sample, row in pca_xy.iterrows():
        is_outlier = bool(outlier_table.at[sample, "potential_outlier"])
        axis.scatter(
            row["PC1"],
            row["PC2"],
            s=130 if is_outlier else 70,
            marker="X" if is_outlier else "o",
        )
        axis.annotate(sample.removeprefix("Abundance."), (row["PC1"], row["PC2"]), xytext=(4, 4), textcoords="offset points", fontsize=8)
    axis.set_title("PCA with Potential Outlier Flags")
    axis.set_xlabel("PC1")
    axis.set_ylabel("PC2")
    figure.tight_layout()
    figure.savefig(FIGURES_DIR / "pca_outlier_flags.png", dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(figure)

    flagged_samples = outlier_table.index[outlier_table["potential_outlier"]].tolist()
    if flagged_samples:
        interpretation = (
            "One or more samples were flagged by at least two independent diagnostics. "
            "Review their raw data, missingness, batch, and experimental notes before "
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
        "flagged_samples": flagged_samples,
        "number_flagged": len(flagged_samples),
        "flag_rule": "At least 2 of 3 diagnostic flags",
        "interpretation": interpretation,
    }
    return outlier_table, diagnostics


def run_diagnostic_workflow(
    df: pd.DataFrame,
    metadata: pd.DataFrame | None = None,
) -> dict[str, Any]:
    """Run all diagnosis-only steps in the requested order."""
    print("\nStarting proteomics diagnostic workflow.")
    print("No transformation, normalization, scaling, or batch correction will be applied.\n")

    if metadata is None:
        metadata = create_sample_metadata(df)

    missing_results = check_missing_values(df)
    distribution_results = assess_distribution(df)
    pca_results = assess_pca(df, metadata=metadata)
    correlation_results = assess_correlation_heatmap(df, method="pearson")
    clustering_results = assess_hierarchical_clustering(df, metadata=metadata)
    outlier_results = assess_outliers(
        df,
        pca_coordinates=pca_results[0],
        correlation_matrix=correlation_results[0],
    )

    print("\n=== Workflow complete ===")
    print(f"Figures saved to: {FIGURES_DIR}")
    print(f"Result tables saved to: {RESULTS_DIR}")
    print("No data modifications were performed.")

    return {
        "metadata": metadata,
        "missing_values": missing_results,
        "distribution": distribution_results,
        "pca": pca_results,
        "correlation": correlation_results,
        "hierarchical_clustering": clustering_results,
        "outliers": outlier_results,
    }


def run_preprocessing_diagnostics(
    df: pd.DataFrame,
    metadata: pd.DataFrame | None = None,
) -> dict[str, Any]:
    """Backward-compatible name for the diagnosis-only workflow.

    This function intentionally performs no log2 transformation, scaling, or
    batch correction. Those execution functions were removed to keep this file
    focused on quality-control diagnosis, as requested.
    """
    return run_diagnostic_workflow(df, metadata=metadata)


def main() -> None:
    """Load the Excel file, promote the header, and run all diagnostics."""
    data = load_raw_data(RAW_DATA_PATH)
    data = promote_first_row_to_header(data)
    print_data_overview(data)

    metadata = create_sample_metadata(data)
    print("\nParsed sample metadata:")
    print(metadata)

    run_diagnostic_workflow(data, metadata=metadata)


if __name__ == "__main__":
    main()
