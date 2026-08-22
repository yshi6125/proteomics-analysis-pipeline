"""One-time loading, structural cleanup, and QC for proteomics input data.

This module prepares the raw workbook for normalization diagnostics and approved
preprocessing. It does not perform normalization, imputation, PCA, clustering,
correlation analysis, or downstream statistical-model specification.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import Any
import re

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
RAW_DATA_PATH = PROJECT_ROOT / "data" / "raw" / "abundantdata.xlsx"
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
QC_CSV_PATH = PROCESSED_DIR / "data_qc.csv"
QC_PICKLE_PATH = PROCESSED_DIR / "data_qc.pkl"
METADATA_PATH = PROCESSED_DIR / "sample_metadata.csv"
SAMPLE_QC_PATH = PROCESSED_DIR / "sample_qc.csv"
PROTEIN_QC_PATH = PROCESSED_DIR / "protein_qc.csv"

REQUIRED_METADATA_COLUMNS = ["sample_id", "group", "batch", "original_column"]
ABUNDANCE_PATTERN = re.compile(r"^Abundance\.(?:NDM|DM)_B\d+_\d+$")
ANALYSIS_NUMERIC_COLUMNS = [
    "PeptideCounts",
    "Number of samples imputed",
    "Fold-change (DM/NDM",
    "pValue",
    "qValue",
    "tStatistic",
]


def ensure_processed_directory() -> None:
    """Create the processed-data directory if it does not exist."""
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)


def load_raw_data(file_path: Path = RAW_DATA_PATH) -> pd.DataFrame:
    """Load the raw Excel workbook without modifying values."""
    if not file_path.exists():
        raise FileNotFoundError(f"Input workbook not found: {file_path}")
    return pd.read_excel(file_path)


def promote_first_row_to_header(data: pd.DataFrame) -> pd.DataFrame:
    """Return a copy whose first data row supplies the column names."""
    if data.empty:
        raise ValueError("Cannot promote a header from an empty dataframe.")
    candidate_header = data.iloc[0]
    if candidate_header.isna().any():
        positions = np.flatnonzero(candidate_header.isna().to_numpy()).tolist()
        raise ValueError(
            "First-row header contains missing values at column positions: "
            f"{positions}"
        )
    if not candidate_header.map(lambda value: isinstance(value, str)).all():
        positions = np.flatnonzero(
            ~candidate_header.map(lambda value: isinstance(value, str)).to_numpy()
        ).tolist()
        raise ValueError(
            "Every promoted header value must be a string. Non-string values occur "
            f"at column positions: {positions}"
        )
    header = candidate_header.astype(str).str.strip().tolist()
    blank_positions = [position for position, value in enumerate(header) if not value]
    if blank_positions:
        raise ValueError(
            f"First-row header contains blank names at positions: {blank_positions}"
        )
    header_index = pd.Index(header)
    if header_index.has_duplicates:
        duplicates = header_index[header_index.duplicated(keep=False)].unique().tolist()
        raise ValueError(f"Promoted header contains duplicate names: {duplicates}")
    if "Accession" not in header or not any(
        value.startswith("Abundance.") for value in header
    ):
        raise ValueError(
            "First row does not look like the expected proteomics header: it must "
            "contain 'Accession' and at least one 'Abundance.*' column."
        )
    promoted = data.copy()
    promoted.columns = header
    return promoted.iloc[1:].reset_index(drop=True)


def identify_raw_abundance_columns(df: pd.DataFrame) -> list[str]:
    """Identify raw abundance columns in their workbook order."""
    abundance_like_columns = [
        column
        for column in df.columns
        if isinstance(column, str) and column.startswith("Abundance.")
    ]
    malformed_columns = [
        column
        for column in abundance_like_columns
        if ABUNDANCE_PATTERN.fullmatch(column) is None
    ]
    if malformed_columns:
        raise ValueError(
            "Malformed Abundance.* column names were found: "
            f"{malformed_columns}. Expected 'Abundance.<group>_B<number>_<number>' "
            "for the supported NDM and DM groups."
        )
    columns = [
        column for column in abundance_like_columns if ABUNDANCE_PATTERN.fullmatch(column)
    ]
    if not columns:
        raise ValueError(
            "No raw abundance columns matched names such as "
            "'Abundance.NDM_B1_1' or 'Abundance.DM_B2_20'."
        )
    return columns


def create_sample_metadata(raw_abundance_columns: Sequence[str]) -> pd.DataFrame:
    """Create metadata matching the contract required by normalization.py."""
    rows: list[dict[str, Any]] = []
    for original_column in raw_abundance_columns:
        raw_sample_id = original_column.removeprefix("Abundance.")
        parts = raw_sample_id.split("_")
        if len(parts) != 3:
            raise ValueError(f"Unexpected abundance-column format: {original_column}")
        group, batch, sample_number_text = parts
        sample_number = int(sample_number_text)
        sample_id = f"{group}_{batch}_{sample_number:02d}"
        rows.append(
            {
                "sample_id": sample_id,
                "group": group,
                "batch": batch,
                "original_column": original_column,
                "sample_number": sample_number,
            }
        )

    metadata = pd.DataFrame(rows)
    if metadata.empty:
        raise ValueError("Cannot create metadata without abundance columns.")
    if metadata["sample_id"].duplicated().any():
        raise ValueError("Standardized sample_id values must be unique.")
    if metadata["original_column"].duplicated().any():
        raise ValueError("Metadata original_column values must be unique.")
    return metadata


def standardize_sample_columns(
    df: pd.DataFrame, metadata: pd.DataFrame
) -> pd.DataFrame:
    """Rename raw abundance columns to standardized sample IDs."""
    missing = [
        column
        for column in metadata["original_column"]
        if column not in df.columns
    ]
    if missing:
        raise ValueError(f"Raw dataframe is missing declared abundance columns: {missing}")
    rename_map = dict(
        zip(metadata["original_column"], metadata["sample_id"], strict=True)
    )
    standardized = df.rename(columns=rename_map).copy()
    if standardized.columns.duplicated().any():
        duplicates = standardized.columns[
            standardized.columns.duplicated(keep=False)
        ].tolist()
        raise ValueError(f"Standardized dataframe has duplicate columns: {duplicates}")
    return standardized


def _convert_numeric_columns(
    df: pd.DataFrame, metadata: pd.DataFrame
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Convert abundance and available analysis columns, rejecting malformed values."""
    sample_ids = metadata["sample_id"].tolist()
    numeric_columns = sample_ids + [
        column for column in ANALYSIS_NUMERIC_COLUMNS if column in df.columns
    ]
    converted = df.copy()
    malformed_by_column: dict[str, list[str]] = {}

    for column in numeric_columns:
        raw = converted[column]
        numeric = pd.to_numeric(raw, errors="coerce")
        malformed = raw.notna() & numeric.isna()
        if malformed.any():
            malformed_by_column[column] = (
                raw.loc[malformed].astype(str).unique()[:5].tolist()
            )
        converted[column] = numeric

    if malformed_by_column:
        raise ValueError(
            "Non-missing values could not be parsed as numeric. Examples by "
            f"column: {malformed_by_column}"
        )

    infinite_by_column = {
        column: int(np.isinf(converted[column].to_numpy(dtype=float)).sum())
        for column in numeric_columns
        if np.isinf(converted[column].to_numpy(dtype=float)).any()
    }
    if infinite_by_column:
        raise ValueError(
            "Infinite numeric values must be resolved before approval: "
            f"{infinite_by_column}"
        )

    return converted, {
        "numeric_columns": numeric_columns,
        "malformed_value_count": 0,
        "infinite_value_count": 0,
    }


def assess_qc_issues(
    df: pd.DataFrame, metadata: pd.DataFrame
) -> dict[str, Any]:
    """Assess data-validity and protein/sample QC issues without exclusions."""
    sample_ids = metadata["sample_id"].tolist()
    abundance = df.loc[:, sample_ids]
    duplicate_rows = df.duplicated(keep=False)
    duplicate_accessions = (
        df["Accession"].duplicated(keep=False)
        if "Accession" in df.columns
        else pd.Series(False, index=df.index)
    )

    peptide_summary: dict[str, Any] = {
        "column_present": "PeptideCounts" in df.columns,
        "missing_count": 0,
        "nonpositive_count": 0,
        "minimum": None,
    }
    if "PeptideCounts" in df.columns:
        peptides = df["PeptideCounts"]
        observed_peptides = peptides.dropna()
        peptide_summary.update(
            {
                "missing_count": int(peptides.isna().sum()),
                "nonpositive_count": int(peptides.le(0).sum()),
                "minimum": (
                    float(observed_peptides.min())
                    if not observed_peptides.empty
                    else None
                ),
            }
        )

    return {
        "number_of_rows": len(df),
        "number_of_samples": len(sample_ids),
        "sample_ids": sample_ids,
        "missing_abundance_count": int(abundance.isna().sum().sum()),
        "missing_abundance_by_sample": abundance.isna().sum().to_dict(),
        "zero_abundance_count": int(abundance.eq(0).sum().sum()),
        "negative_abundance_count": int(abundance.lt(0).sum().sum()),
        "malformed_abundance_count": 0,
        "infinite_abundance_count": 0,
        "duplicate_row_count": int(duplicate_rows.sum()),
        "duplicate_accession_row_count": int(duplicate_accessions.sum()),
        "duplicate_accessions": (
            df.loc[duplicate_accessions, "Accession"].dropna().unique().tolist()
            if "Accession" in df.columns
            else []
        ),
        "peptide_count_qc": peptide_summary,
    }


def apply_nonpositive_value_policies(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    zero_policy: str,
    negative_policy: str,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Apply explicit policies to zero and negative abundance values.

    ``zero_policy`` accepts ``"error"``, ``"set_missing"``, or ``"allow"``.
    ``negative_policy`` accepts ``"error"``, ``"exclude_rows"``, or ``"allow"``.
    No pseudocount or value shift is performed.
    """
    allowed_zero_policies = {"error", "set_missing", "allow"}
    allowed_negative_policies = {"error", "exclude_rows", "allow"}
    if zero_policy not in allowed_zero_policies:
        raise ValueError(
            "zero_policy must be one of "
            f"{sorted(allowed_zero_policies)}; received {zero_policy!r}."
        )
    if negative_policy not in allowed_negative_policies:
        raise ValueError(
            "negative_policy must be one of "
            f"{sorted(allowed_negative_policies)}; received {negative_policy!r}."
        )

    resolved = df.copy()
    sample_ids = metadata["sample_id"].tolist()
    abundance = resolved.loc[:, sample_ids]
    zero_mask = abundance.eq(0)
    negative_mask = abundance.lt(0)
    zero_count = int(zero_mask.sum().sum())
    negative_count = int(negative_mask.sum().sum())
    negative_row_mask = negative_mask.any(axis=1)
    negative_row_indices = resolved.index[negative_row_mask].tolist()

    if zero_count and zero_policy == "error":
        raise ValueError(
            f"Found {zero_count} zero abundance values. Choose zero_policy='allow' "
            "or zero_policy='set_missing' only after scientific review."
        )
    if negative_count and negative_policy == "error":
        raise ValueError(
            f"Found {negative_count} negative abundance values. Choose "
            "negative_policy='allow' or negative_policy='exclude_rows' only after "
            "scientific review."
        )
    if zero_policy == "set_missing":
        abundance = abundance.mask(zero_mask)
    resolved.loc[:, sample_ids] = abundance
    if negative_policy == "exclude_rows":
        resolved = resolved.loc[~negative_row_mask].copy()

    return resolved, {
        "zero_policy": zero_policy,
        "negative_policy": negative_policy,
        "zero_values_found": zero_count,
        "negative_values_found": negative_count,
        "zero_values_set_missing": zero_count if zero_policy == "set_missing" else 0,
        "negative_rows_excluded": (
            len(negative_row_indices) if negative_policy == "exclude_rows" else 0
        ),
        "negative_row_indices": (
            negative_row_indices if negative_policy == "exclude_rows" else []
        ),
    }


def create_qc_tables(
    df: pd.DataFrame, metadata: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create ordered sample-level and protein-level QC tables."""
    sample_ids = metadata["sample_id"].tolist()
    abundance = df.loc[:, sample_ids]

    sample_qc = metadata.loc[:, REQUIRED_METADATA_COLUMNS].copy()
    sample_qc["observed_count"] = abundance.notna().sum(axis=0).to_numpy()
    sample_qc["missing_count"] = abundance.isna().sum(axis=0).to_numpy()
    sample_qc["missing_percent"] = (
        abundance.isna().mean(axis=0).mul(100).to_numpy()
    )
    sample_qc["zero_count"] = abundance.eq(0).sum(axis=0).to_numpy()
    sample_qc["negative_count"] = abundance.lt(0).sum(axis=0).to_numpy()
    sample_qc["minimum_observed"] = abundance.min(axis=0).to_numpy()
    sample_qc["median_observed"] = abundance.median(axis=0).to_numpy()
    sample_qc["maximum_observed"] = abundance.max(axis=0).to_numpy()

    protein_qc = pd.DataFrame({"row_index": df.index})
    for column in ("Accession", "Gene.name", "PeptideCounts"):
        if column in df.columns:
            protein_qc[column] = df[column].to_numpy()
    protein_qc["observed_sample_count"] = abundance.notna().sum(axis=1).to_numpy()
    protein_qc["missing_sample_count"] = abundance.isna().sum(axis=1).to_numpy()
    protein_qc["missing_percent"] = abundance.isna().mean(axis=1).mul(100).to_numpy()
    protein_qc["zero_count"] = abundance.eq(0).sum(axis=1).to_numpy()
    protein_qc["negative_count"] = abundance.lt(0).sum(axis=1).to_numpy()
    for group in metadata["group"].drop_duplicates().tolist():
        group_sample_ids = metadata.loc[
            metadata["group"].eq(group), "sample_id"
        ].tolist()
        protein_qc[f"observed_{group}"] = (
            abundance.loc[:, group_sample_ids].notna().sum(axis=1).to_numpy()
        )
    protein_qc["duplicate_full_row"] = df.duplicated(keep=False).to_numpy()
    protein_qc["duplicate_accession"] = (
        df["Accession"].duplicated(keep=False).to_numpy()
        if "Accession" in df.columns
        else False
    )
    return sample_qc, protein_qc


def validate_final_contract(df: pd.DataFrame, metadata: pd.DataFrame) -> list[str]:
    """Validate the final dataframe–metadata contract before persistence."""
    if df.empty:
        raise ValueError(
            "Final QC dataframe must contain at least one protein row. Review the "
            "approved exclusion criteria before saving."
        )
    missing_metadata_columns = [
        column for column in REQUIRED_METADATA_COLUMNS if column not in metadata.columns
    ]
    if missing_metadata_columns:
        raise ValueError(
            f"Final metadata is missing required columns: {missing_metadata_columns}"
        )
    if metadata.empty:
        raise ValueError("Final metadata must contain at least one sample.")
    required = metadata.loc[:, REQUIRED_METADATA_COLUMNS]
    if required.isna().any().any():
        columns = required.columns[required.isna().any()].tolist()
        raise ValueError(f"Final metadata contains missing required values: {columns}")
    for column in REQUIRED_METADATA_COLUMNS:
        invalid = ~metadata[column].map(
            lambda value: isinstance(value, str) and bool(value.strip())
        )
        if invalid.any():
            raise ValueError(
                f"Final metadata {column} values must be non-empty strings."
            )
    sample_ids = metadata["sample_id"].tolist()
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Final metadata sample_id values must be unique.")
    if metadata["original_column"].duplicated().any():
        raise ValueError("Final metadata original_column values must be unique.")
    if df.columns.duplicated().any():
        duplicates = df.columns[df.columns.duplicated(keep=False)].tolist()
        raise ValueError(f"Final dataframe contains duplicate columns: {duplicates}")
    missing_samples = [sample_id for sample_id in sample_ids if sample_id not in df.columns]
    if missing_samples:
        raise ValueError(
            f"Final dataframe is missing metadata sample columns: {missing_samples}"
        )
    nonnumeric_samples = [
        sample_id
        for sample_id in sample_ids
        if not pd.api.types.is_numeric_dtype(df[sample_id])
    ]
    if nonnumeric_samples:
        raise ValueError(
            f"Final abundance columns must be numeric: {nonnumeric_samples}"
        )
    infinite_count = int(
        np.isinf(df.loc[:, sample_ids].to_numpy(dtype=float)).sum()
    )
    if infinite_count:
        raise ValueError(
            f"Final abundance data contain {infinite_count} infinite values."
        )
    if df.attrs.get("abundance_scale") != "linear":
        raise ValueError(
            "Final QC dataframe must have abundance_scale='linear' before saving."
        )
    return sample_ids


def apply_approved_exclusions(
    df: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    exclude_sample_ids: Sequence[str] | None = None,
    exclude_row_indices: Sequence[Any] | None = None,
    exclude_accessions: Sequence[Any] | None = None,
    minimum_peptide_count: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Apply only exclusions explicitly supplied by the user.

    When ``minimum_peptide_count`` is requested, proteins with missing peptide
    counts are excluded because they cannot be shown to meet the cutoff.
    """
    approved_df = df.copy()
    approved_metadata = metadata.copy()

    requested_samples = list(exclude_sample_ids or [])
    unknown_samples = sorted(
        set(requested_samples).difference(approved_metadata["sample_id"])
    )
    if unknown_samples:
        raise ValueError(f"Requested sample exclusions were not found: {unknown_samples}")
    if requested_samples:
        approved_df = approved_df.drop(columns=requested_samples)
        approved_metadata = approved_metadata.loc[
            ~approved_metadata["sample_id"].isin(requested_samples)
        ].reset_index(drop=True)

    row_mask = pd.Series(False, index=approved_df.index)
    requested_indices = list(exclude_row_indices or [])
    missing_indices = [index for index in requested_indices if index not in approved_df.index]
    if missing_indices:
        raise ValueError(f"Requested row indices were not found: {missing_indices}")
    row_mask |= approved_df.index.isin(requested_indices)

    requested_accessions = list(exclude_accessions or [])
    if requested_accessions:
        if "Accession" not in approved_df.columns:
            raise ValueError("Cannot exclude accessions: Accession column is absent.")
        missing_accessions = sorted(
            set(requested_accessions).difference(approved_df["Accession"].dropna())
        )
        if missing_accessions:
            raise ValueError(
                f"Requested accession exclusions were not found: {missing_accessions}"
            )
        row_mask |= approved_df["Accession"].isin(requested_accessions)

    missing_peptide_count_rows = 0
    below_minimum_peptide_count_rows = 0
    if minimum_peptide_count is not None:
        if "PeptideCounts" not in approved_df.columns:
            raise ValueError("Cannot use minimum_peptide_count without PeptideCounts.")
        if (
            isinstance(minimum_peptide_count, bool)
            or not isinstance(minimum_peptide_count, (int, np.integer))
            or minimum_peptide_count < 1
        ):
            raise ValueError("minimum_peptide_count must be a positive integer.")
        missing_peptide_mask = approved_df["PeptideCounts"].isna()
        below_minimum_mask = approved_df["PeptideCounts"].lt(
            minimum_peptide_count
        )
        missing_peptide_count_rows = int(missing_peptide_mask.sum())
        below_minimum_peptide_count_rows = int(below_minimum_mask.sum())
        row_mask |= missing_peptide_mask | below_minimum_mask

    excluded_rows = int(row_mask.sum())
    approved_df = approved_df.loc[~row_mask].reset_index(drop=True)
    if approved_df.empty:
        raise ValueError(
            "Approved exclusions removed every protein row. Revise the row, "
            "accession, peptide-count, or negative-value exclusion criteria."
        )
    if approved_metadata.empty:
        raise ValueError("Approved exclusions removed every sample.")
    return approved_df, approved_metadata, {
        "excluded_sample_ids": requested_samples,
        "excluded_row_count": excluded_rows,
        "excluded_row_indices": requested_indices,
        "excluded_accessions": requested_accessions,
        "minimum_peptide_count": minimum_peptide_count,
        "missing_peptide_count_rows_excluded": missing_peptide_count_rows,
        "below_minimum_peptide_count_rows_excluded": (
            below_minimum_peptide_count_rows
        ),
    }


def save_qc_outputs(
    qc_df: pd.DataFrame,
    metadata: pd.DataFrame,
    sample_qc: pd.DataFrame,
    protein_qc: pd.DataFrame,
) -> None:
    """Save the QC-approved dataframe and standardized metadata."""
    validate_final_contract(qc_df, metadata)
    ensure_processed_directory()
    qc_df.to_csv(QC_CSV_PATH, index=False)
    qc_df.to_pickle(QC_PICKLE_PATH)
    metadata.to_csv(METADATA_PATH, index=False)
    sample_qc.to_csv(SAMPLE_QC_PATH, index=False)
    protein_qc.to_csv(PROTEIN_QC_PATH, index=False)


def run_qc_workflow(
    file_path: Path = RAW_DATA_PATH,
    *,
    promote_header: bool = True,
    exclude_sample_ids: Sequence[str] | None = None,
    exclude_row_indices: Sequence[Any] | None = None,
    exclude_accessions: Sequence[Any] | None = None,
    minimum_peptide_count: int | None = None,
    zero_policy: str = "error",
    negative_policy: str = "error",
    save_outputs: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Run the complete one-time QC preparation workflow."""
    raw_df = load_raw_data(file_path)
    structured_df = (
        promote_first_row_to_header(raw_df) if promote_header else raw_df.copy()
    )
    raw_abundance_columns = identify_raw_abundance_columns(structured_df)
    metadata = create_sample_metadata(raw_abundance_columns)
    standardized_df = standardize_sample_columns(structured_df, metadata)
    numeric_df, conversion_summary = _convert_numeric_columns(
        standardized_df, metadata
    )
    issues_before_exclusions = assess_qc_issues(numeric_df, metadata)
    policy_df, nonpositive_policy_summary = apply_nonpositive_value_policies(
        numeric_df,
        metadata,
        zero_policy=zero_policy,
        negative_policy=negative_policy,
    )
    qc_df, approved_metadata, exclusion_summary = apply_approved_exclusions(
        policy_df,
        metadata,
        exclude_sample_ids=exclude_sample_ids,
        exclude_row_indices=exclude_row_indices,
        exclude_accessions=exclude_accessions,
        minimum_peptide_count=minimum_peptide_count,
    )
    issues_after_exclusions = assess_qc_issues(qc_df, approved_metadata)
    qc_df.attrs["abundance_scale"] = "linear"
    validate_final_contract(qc_df, approved_metadata)
    sample_qc, protein_qc = create_qc_tables(qc_df, approved_metadata)

    if save_outputs:
        save_qc_outputs(qc_df, approved_metadata, sample_qc, protein_qc)

    summary = {
        "input_path": str(file_path),
        "header_promoted": promote_header,
        "conversion": conversion_summary,
        "issues_before_exclusions": issues_before_exclusions,
        "nonpositive_value_policies": nonpositive_policy_summary,
        "approved_exclusions": exclusion_summary,
        "issues_after_exclusions": issues_after_exclusions,
        "sample_qc": sample_qc,
        "protein_qc": protein_qc,
        "final_contract_validated": True,
        "outputs_saved": save_outputs,
    }
    return qc_df, approved_metadata, summary


def main() -> None:
    """Prepare the raw workbook without applying any automatic exclusions."""
    qc_df, metadata, summary = run_qc_workflow()
    print(f"QC-approved dataframe shape: {qc_df.shape}")
    print(f"Number of approved samples: {len(metadata)}")
    print(f"QC issue summary: {summary['issues_after_exclusions']}")
    print(f"Saved dataframe: {QC_CSV_PATH}")
    print(f"Saved dataframe with attributes: {QC_PICKLE_PATH}")
    print(f"Saved metadata: {METADATA_PATH}")
    print(f"Saved sample-level QC table: {SAMPLE_QC_PATH}")
    print(f"Saved protein-level QC table: {PROTEIN_QC_PATH}")


if __name__ == "__main__":
    main()
