"""Optional preprocessing operations for quantitative proteomics abundances.

This module contains transformations that modify abundance values in a returned
copy. It deliberately contains no quality-control workflow and executes nothing
automatically. Review the diagnostic results before calling these functions.

The intended abundance-scale workflow is::

    QC
      -> optional conservative scaling (linear scale only)
      -> log2 transformation
      -> downstream statistics with batch included in the model

``DataFrame.attrs["abundance_scale"]`` records the current abundance scale across
these in-memory operations. Batch effects are not removed here; preparation batch
is included as a covariate in the downstream differential-abundance model.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

import numpy as np
import pandas as pd


REQUIRED_METADATA_COLUMNS = {"sample_id", "group", "batch", "original_column"}


def _validate_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    """Validate and copy standardized sample metadata.

    Parameters
    ----------
    metadata
        Sample metadata following the standardized metadata contract.

    Returns
    -------
    pandas.DataFrame
        Validated metadata copy with its original row order preserved.

    Raises
    ------
    ValueError
        If required columns or values are missing, or identifiers are duplicated.
    """
    if metadata is None:
        raise ValueError("metadata is required to resolve abundance columns.")

    missing_metadata_columns = sorted(
        REQUIRED_METADATA_COLUMNS.difference(metadata.columns)
    )
    if missing_metadata_columns:
        raise ValueError(
            "Metadata is missing required columns: "
            f"{missing_metadata_columns}. Required columns are sample_id, group, "
            "batch, and original_column."
        )

    duplicated_required_columns = sorted(
        column
        for column in REQUIRED_METADATA_COLUMNS
        if list(metadata.columns).count(column) > 1
    )
    if duplicated_required_columns:
        raise ValueError(
            "Required metadata column names must be unique. Duplicates: "
            f"{duplicated_required_columns}."
        )

    # Normalize user-entered labels before testing for missing/empty values and
    # duplicates.  In particular, ``" S1 "`` and ``"S1"`` must identify the
    # same sample rather than two superficially distinct samples.
    metadata_copy = metadata.copy()
    for column in ("sample_id", "group", "batch", "original_column"):
        metadata_copy[column] = metadata_copy[column].map(
            lambda value: value.strip() if isinstance(value, str) else value
        )
    required_metadata = metadata_copy.loc[
        :, ["sample_id", "group", "batch", "original_column"]
    ]
    columns_with_missing_values = required_metadata.columns[
        required_metadata.isna().any()
    ].tolist()
    if columns_with_missing_values:
        raise ValueError(
            "Required metadata columns contain missing values: "
            f"{columns_with_missing_values}."
        )

    sample_ids = metadata_copy["sample_id"].tolist()
    if not sample_ids:
        raise ValueError("Metadata does not contain any sample rows.")
    if not all(isinstance(sample_id, str) and sample_id for sample_id in sample_ids):
        raise ValueError("Every metadata sample_id must be a non-empty string.")
    for column in ("group", "batch", "original_column"):
        empty_string_count = int(
            metadata_copy[column].map(
                lambda value: isinstance(value, str) and not value
            ).sum()
        )
        if empty_string_count:
            raise ValueError(
                f"Metadata {column} values must not be empty strings after "
                "whitespace is removed."
            )
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Metadata sample_id values must be unique.")
    if metadata_copy["original_column"].duplicated().any():
        raise ValueError("Metadata original_column values must be unique.")
    return metadata_copy


def _resolve_abundance_columns(
    df: pd.DataFrame, metadata: pd.DataFrame
) -> list[str]:
    """Resolve abundance columns from the standardized sample metadata.

    Parameters
    ----------
    df
        Proteomics dataframe containing standardized sample-ID columns.
    metadata
        Sample metadata containing ``sample_id``, ``group``, ``batch``, and
        ``original_column``. Its row order defines abundance-column order.

    Returns
    -------
    list of str
        Sample-ID column names in authoritative metadata order.

    Raises
    ------
    ValueError
        If the metadata contract is invalid or a declared sample column is absent.
    """
    validated_metadata = _validate_metadata(metadata)
    abundance_cols = validated_metadata["sample_id"].tolist()

    missing_abundance_cols = [
        sample_id for sample_id in abundance_cols if sample_id not in df.columns
    ]
    if missing_abundance_cols:
        raise ValueError(
            "Dataframe is missing abundance columns declared by metadata: "
            f"{missing_abundance_cols}."
        )
    duplicated_dataframe_columns = df.columns[df.columns.duplicated()].tolist()
    duplicated_abundance_cols = [
        sample_id
        for sample_id in abundance_cols
        if sample_id in duplicated_dataframe_columns
    ]
    if duplicated_abundance_cols:
        raise ValueError(
            "Dataframe abundance-column names must be unique. Duplicates: "
            f"{duplicated_abundance_cols}."
        )
    return abundance_cols


def _copy_numeric_abundances(
    df: pd.DataFrame, abundance_cols: Sequence[str]
) -> pd.DataFrame:
    """Create a numeric abundance copy without changing the input dataframe.

    Parameters
    ----------
    df
        Source proteomics dataframe.
    abundance_cols
        Abundance columns to copy and convert.

    Returns
    -------
    pandas.DataFrame
        Numeric abundance values with unparseable entries represented as missing.

    Raises
    ------
    ValueError
        If malformed non-missing text or infinity is present in an abundance column.
    """
    raw_abundances = df.loc[:, list(abundance_cols)].copy()
    abundance_df = raw_abundances.apply(pd.to_numeric, errors="coerce")

    # Distinguish genuinely missing values from non-missing text that failed parsing.
    malformed_mask = raw_abundances.notna() & abundance_df.isna()
    malformed_count = int(malformed_mask.sum().sum())
    if malformed_count:
        malformed_columns = malformed_mask.columns[
            malformed_mask.any(axis=0)
        ].tolist()
        malformed_examples = {
            column: raw_abundances.loc[malformed_mask[column], column]
            .astype(str)
            .unique()[:3]
            .tolist()
            for column in malformed_columns
        }
        raise ValueError(
            f"Abundance data contain {malformed_count} non-missing values that could "
            f"not be parsed as numeric. Examples by column: {malformed_examples}."
        )

    infinite_mask = np.isinf(abundance_df.to_numpy(dtype=float))
    infinite_count = int(infinite_mask.sum())
    if infinite_count:
        infinite_columns = abundance_df.columns[infinite_mask.any(axis=0)].tolist()
        raise ValueError(
            f"Abundance data contain {infinite_count} infinite values in columns "
            f"{infinite_columns}. Replace or investigate them before preprocessing."
        )
    return abundance_df


def apply_log2_transformation(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Apply log2 transformation to abundance columns in a returned copy.

    Parameters
    ----------
    df
        Proteomics dataframe. Non-abundance columns and row order are preserved.
    metadata
        Standardized sample metadata. ``sample_id`` defines the abundance columns
        and their authoritative order.
    Returns
    -------
    transformed_df : pandas.DataFrame
        New dataframe with transformed values replacing only abundance columns.
    summary : dict
        Counts and parameters describing the transformation.

    Raises
    ------
    ValueError
        If the data are already marked as log2-scaled, or if zero, negative, or
        infinite abundance values are present.

    Notes
    -----
    This operation assumes abundance values are on a linear, strictly positive
    scale. Missing values remain missing. Nonpositive values are never adjusted.
    """
    # Work from numeric copies so the caller's dataframe remains untouched.
    resolved_scale = df.attrs.get("abundance_scale", "linear")
    if isinstance(resolved_scale, str) and resolved_scale.strip().lower() == "log2":
        raise ValueError(
            "Log2 transformation cannot be applied because the dataframe is already "
            "marked as log2-scaled."
        )

    abundance_cols = _resolve_abundance_columns(df, metadata)
    abundance_df = _copy_numeric_abundances(df, abundance_cols)

    positive_count = int(abundance_df.gt(0).sum().sum())
    zero_count = int(abundance_df.eq(0).sum().sum())
    negative_count = int(abundance_df.lt(0).sum().sum())
    missing_count = int(abundance_df.isna().sum().sum())
    observed_value_count = int(abundance_df.notna().sum().sum())

    if observed_value_count == 0:
        raise ValueError(
            "No observed numeric abundance values are available for log2 "
            "transformation."
        )

    observed_values = abundance_df.stack()
    minimum_observed_abundance = float(observed_values.min())
    maximum_observed_abundance = float(observed_values.max())

    print(f"Observed abundance values: {observed_value_count}")
    print(f"Positive abundance values: {positive_count}")
    print(f"Zero abundance values: {zero_count}")
    print(f"Negative abundance values: {negative_count}")
    print(f"Missing abundance values: {missing_count}")
    print(f"Minimum observed abundance: {minimum_observed_abundance}")
    print(f"Maximum observed abundance: {maximum_observed_abundance}")

    if negative_count:
        raise ValueError(
            "Log2 transformation is not defined directly for negative abundance "
            "values. Their origin must be reviewed; this function will not shift them."
        )
    if zero_count:
        raise ValueError(
            "Log2 transformation is not defined for zero abundance values. Determine "
            "whether they are true zeros, missing values, or below-detection values "
            "before preprocessing; this function will not replace them."
        )

    # All observed values are finite and strictly positive at this point.
    transformed_abundances = np.log2(abundance_df)

    transformed_df = df.copy()
    transformed_df[abundance_cols] = transformed_abundances
    transformed_df.attrs["abundance_scale"] = "log2"
    transformed_count = int(transformed_abundances.notna().sum().sum())
    print(
        f"Log2 transformation completed for {transformed_count} values across "
        f"{len(abundance_cols)} abundance columns. The original dataframe was not "
        "modified."
    )
    print(
        "Assumption: all observed abundance values are finite, strictly positive "
        "measurements on a linear scale."
    )

    summary: dict[str, Any] = {
        "method": "log2",
        "number_of_samples": len(abundance_cols),
        "number_of_rows": len(df),
        "sample_ids": abundance_cols,
        "observed_value_count": observed_value_count,
        "transformed_value_count": transformed_count,
        "positive_count": positive_count,
        "zero_count": zero_count,
        "negative_count": negative_count,
        "missing_count": missing_count,
        "infinite_count": 0,
        "minimum_observed_abundance": minimum_observed_abundance,
        "maximum_observed_abundance": maximum_observed_abundance,
    }
    return transformed_df, summary


def apply_conservative_scaling(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    method: str = "median",
    reference_proteins: Sequence[Any] | None = None,
    protein_id_col: str = "Accession",
    target_value: float | None = None,
    abundance_scale: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Apply one multiplicative scaling factor to each abundance sample.

    Parameters
    ----------
    df
        Proteomics dataframe. Metadata columns and row order are preserved.
    metadata
        Standardized sample metadata. ``sample_id`` defines the abundance columns
        and their authoritative order.
    method
        ``"median"`` uses all observed proteins; ``"reference_proteins"`` uses
        only explicitly supplied stable proteins.
    reference_proteins
        Stable protein identifiers required for reference-protein scaling.
    protein_id_col
        Column containing protein identifiers used to locate references.
    target_value
        Optional positive target median. By default, the median of the
        sample-specific medians is used.
    abundance_scale
        Scale of the supplied abundance values. If omitted, the function uses the
        dataframe's ``abundance_scale`` attribute and otherwise assumes linear data.

    Returns
    -------
    scaled_df : pandas.DataFrame
        New dataframe with scaled abundance values and unchanged metadata.
    scaling_factors : pandas.DataFrame
        Sample medians, target, factors, and method.
    summary : dict
        Scaling settings and samples receiving diagnostic warnings.

    Raises
    ------
    ValueError
        If the method, references, target, sample medians, or abundance scale are
        invalid.

    Notes
    -----
    Scaling must precede log2 transformation. Exactly one factor is applied to every
    observed value within each sample.
    This does not force complete distributions to match. Median scaling assumes
    most proteins have comparable central abundance. Reference scaling assumes
    the chosen proteins are stable across biology and technical batches.
    """
    resolved_scale = (
        abundance_scale
        if abundance_scale is not None
        else df.attrs.get("abundance_scale", "linear")
    )
    if isinstance(resolved_scale, str) and resolved_scale.strip().lower() == "log2":
        raise ValueError(
            "Multiplicative sample scaling must be applied before log2 transformation."
        )

    if method not in {"median", "reference_proteins"}:
        raise ValueError("method must be 'median' or 'reference_proteins'.")

    abundance_cols = _resolve_abundance_columns(df, metadata)
    validated_metadata = _validate_metadata(metadata).reset_index(drop=True)
    abundance_df = _copy_numeric_abundances(df, abundance_cols)

    # Select the measurements that define the sample-specific reference medians.
    found_references: list[Any] = []
    missing_references: list[Any] = []
    if method == "median":
        reference_abundances = abundance_df
        print(
            "Median scaling selected. Sample medians were calculated from all "
            "observed abundance values."
        )
        print(
            "Assumption: the central abundance of most proteins should be comparable "
            "and observed global shifts are primarily technical rather than biological."
        )
    else:
        if reference_proteins is None or len(reference_proteins) < 2:
            raise ValueError(
                "At least two validated stable reference proteins are required for "
                "method='reference_proteins'."
            )
        if protein_id_col not in df.columns:
            raise ValueError(f"Protein identifier column not found: {protein_id_col}")

        requested_references = list(reference_proteins)
        if len(requested_references) != len(set(requested_references)):
            duplicate_requests = pd.Series(requested_references)[
                pd.Series(requested_references).duplicated(keep=False)
            ].unique().tolist()
            raise ValueError(
                "reference_proteins must contain unique identifiers. Duplicate "
                f"requests: {duplicate_requests}."
            )
        requested_set = set(requested_references)
        observed_ids = set(df[protein_id_col].dropna())
        found_references = [
            identifier for identifier in requested_references if identifier in observed_ids
        ]
        missing_references = [
            identifier for identifier in requested_references if identifier not in observed_ids
        ]
        print(f"Reference proteins found: {found_references}")
        print(f"Reference proteins missing: {missing_references}")
        if len(found_references) < 2:
            raise ValueError(
                "Fewer than two unique requested reference proteins were matched. "
                f"Found: {found_references}; missing: {missing_references}."
            )

        reference_mask = df[protein_id_col].isin(requested_set)
        matched_identifiers = df.loc[reference_mask, protein_id_col]
        duplicated_matches = matched_identifiers[
            matched_identifiers.duplicated(keep=False)
        ].unique().tolist()
        if duplicated_matches:
            raise ValueError(
                "Each reference protein must match exactly one dataframe row. "
                f"Duplicate matched identifiers: {duplicated_matches}."
            )
        reference_abundances = abundance_df.loc[reference_mask]
        if reference_abundances.empty:
            raise ValueError("No abundance rows matched the requested reference proteins.")
        print(
            "Reference-protein scaling selected. Sample medians were calculated only "
            "from the requested stable proteins."
        )
        print(
            "Assumption: the selected reference proteins remain stable across the "
            "biological conditions and technical batches being compared."
        )
    if method == "reference_proteins":
        reference_counts = reference_abundances.count(axis=0)

        insufficient_samples = reference_counts[
            reference_counts < 2
        ].index.tolist()

        if insufficient_samples:
            raise ValueError(
                "Fewer than two reference proteins have observed values in these "
                f"samples: {insufficient_samples}"
            )
    sample_medians = reference_abundances.median(axis=0)
    invalid_medians = sample_medians[
        sample_medians.isna() | ~np.isfinite(sample_medians) | sample_medians.le(0)
    ]
    if not invalid_medians.empty:
        raise ValueError(
            "Scaling requires a positive, finite observed median for every sample. "
            f"Invalid sample medians: {invalid_medians.to_dict()}"
        )

    # Resolve the common target without choosing any distribution-forcing method.
    if target_value is None:
        resolved_target = float(sample_medians.median())
    elif (
        isinstance(target_value, bool)
        or not isinstance(target_value, (int, float, np.number))
        or not np.isfinite(target_value)
        or target_value <= 0
    ):
        raise ValueError("target_value must be a positive finite number.")
    else:
        resolved_target = float(target_value)

    factors = resolved_target / sample_medians
    invalid_factors = factors[~np.isfinite(factors)]
    if not invalid_factors.empty:
        raise ValueError(
            "Calculated scaling factors must be finite. Invalid factors: "
            f"{invalid_factors.to_dict()}"
        )
    factor_metrics = pd.DataFrame(
        {
            "sample_id": abundance_cols,
            "sample_median": sample_medians.reindex(abundance_cols).to_numpy(),
            "target_median": resolved_target,
            "scaling_factor": factors.reindex(abundance_cols).to_numpy(),
            "method": method,
        }
    )
    calculated_columns = {
        "sample_median",
        "target_median",
        "scaling_factor",
        "method",
    }
    metadata_collisions = sorted(calculated_columns.intersection(metadata.columns))
    if metadata_collisions:
        raise ValueError(
            "Metadata columns conflict with calculated scaling-factor fields: "
            f"{metadata_collisions}."
        )
    scaling_factors = validated_metadata.merge(
        factor_metrics,
        on="sample_id",
        how="left",
        sort=False,
        validate="one_to_one",
    )

    # Apply exactly one factor per sample; pandas naturally preserves missing values.
    scaled_abundances = abundance_df.mul(factors, axis="columns")
    scaled_infinite_count = int(
        np.isinf(scaled_abundances.to_numpy(dtype=float)).sum()
    )
    if scaled_infinite_count:
        raise ValueError(
            "Scaling produced infinite abundance values. Review the input magnitude "
            "and calculated factors before proceeding."
        )
    scaled_df = df.copy()
    scaled_df[abundance_cols] = scaled_abundances
    scaled_df.attrs["abundance_scale"] = "linear_scaled"

    warning_mask = scaling_factors["scaling_factor"].lt(0.5) | scaling_factors[
        "scaling_factor"
    ].gt(2.0)
    warning_samples = scaling_factors.loc[warning_mask, "sample_id"].tolist()
    print("Scaling factors:")
    print(scaling_factors.to_string(index=False))
    if warning_samples:
        print(
            "Warning: scaling factors below 0.5 or above 2.0 were found for "
            f"{warning_samples}. This threshold is a diagnostic warning, not an "
            "automatic exclusion criterion."
        )
    print(
        "Conservative scaling completed. Exactly one multiplicative factor was "
        "applied per sample. The original dataframe was not modified."
    )

    summary: dict[str, Any] = {
        "method": method,
        "target_value": resolved_target,
        "number_of_samples": len(abundance_cols),
        "number_of_rows": len(df),
        "sample_ids": abundance_cols,
        "minimum_scaling_factor": float(factors.min()),
        "maximum_scaling_factor": float(factors.max()),
        "warning_samples": warning_samples,
        "reference_proteins_requested": (
            list(reference_proteins) if reference_proteins is not None else []
        ),
        "reference_proteins_found": found_references,
        "reference_proteins_missing": missing_references,
        "infinite_count": 0,
    }
    return scaled_df, scaling_factors, summary


# Manual workflow only; preprocessing is never executed automatically.
# Begin with the dataframe returned by QC.
#
# df_scaled, scaling_factors, scaling_info = apply_conservative_scaling(
#     df_qc,
#     metadata=metadata,
#     method="median",
# )
#
# df_log2, log2_info = apply_log2_transformation(
#     df_scaled,
#     metadata=metadata,
# )
#
# Save df_log2 as the production differential-abundance input. Batch is modeled
# simultaneously with group in differential_abundance.py; do not subtract it here.
