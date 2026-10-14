"""Protein evidence and memberships for the single Levels 1–5 production Route A."""
from __future__ import annotations
from collections.abc import Sequence
from pathlib import Path
from typing import Any
import numpy as np
import pandas as pd
from src.pathway_review_io import ROOT
DIFFERENTIAL_RESULTS_DIR = ROOT / 'results' / 'differential_abundance'

REQUIRED_DIFFERENTIAL_COLUMNS = [
    "row_position",
    "Accession",
    "Gene.name",
    "log2FC",
    "t_statistic",
    "p_value",
    "q_value",
]


CANDIDATE_ORIGINS = ["Level 1", "Level 2", "Level 3", "Level 4", "Level 5"]


ORIGIN_PRIORITY = {origin: priority for priority, origin in enumerate(CANDIDATE_ORIGINS, 1)}


ENRICHMENT_COLUMNS = [
    "analysis_type",
    "analysis_label",
    "direction",
    "pathway_library",
    "pathway",
    "enrichment_p_value",
    "enrichment_adjusted_p_value",
    "overlap_count",
    "overlap_size",
    "overlap",
    "input_candidate_count",
    "gene_ratio",
    "candidate_genes",
]


MEMBER_COLUMNS = [
    "analysis_type",
    "analysis_label",
    "direction",
    "pathway_library",
    "pathway",
    "candidate_gene",
    "row_position",
    "Gene.name",
    "Accession",
    "candidate_origin",
    "log2FC",
    "t_statistic",
    "p_value",
    "q_value",
]


def load_differential_abundance_results(
    comparison: str = "DM_vs_NDM",
    results_dir: Path = DIFFERENTIAL_RESULTS_DIR,
) -> pd.DataFrame:
    """Load the complete differential-abundance table for one comparison."""
    path = results_dir / f"differential_abundance_{comparison}.csv"
    if not path.exists():
        raise FileNotFoundError(f"Differential-abundance result table not found: {path}")
    return pd.read_csv(path)


def validate_differential_results(results: pd.DataFrame) -> None:
    """Validate fields required for candidate-based pathway analysis."""
    missing = [column for column in REQUIRED_DIFFERENTIAL_COLUMNS if column not in results]
    if missing:
        raise ValueError(f"Differential-abundance results are missing columns: {missing}")
    if results.empty:
        raise ValueError("Differential-abundance results contain no tested proteins.")
    if results["row_position"].isna().any() or results["row_position"].duplicated().any():
        raise ValueError("row_position must be non-missing and unique.")
    numeric_columns = ["log2FC", "t_statistic", "p_value", "q_value"]
    numeric = results.loc[:, numeric_columns].apply(pd.to_numeric, errors="coerce")
    malformed = numeric.isna() & results.loc[:, numeric_columns].notna()
    if malformed.any().any():
        columns = malformed.columns[malformed.any()].tolist()
        raise ValueError(f"Differential statistics contain malformed values: {columns}")
    if numeric.isna().any().any() or not np.isfinite(numeric.to_numpy()).all():
        raise ValueError("Differential statistics must be finite and non-missing.")
    if ((numeric["p_value"] < 0) | (numeric["p_value"] > 1)).any():
        raise ValueError("p_value values must be between 0 and 1.")
    if ((numeric["q_value"] < 0) | (numeric["q_value"] > 1)).any():
        raise ValueError("q_value values must be between 0 and 1.")


def assign_candidate_origin(results: pd.DataFrame) -> pd.DataFrame:
    """Assign each protein to exactly one strictest candidate-origin level."""
    validate_differential_results(results)
    candidates = results.copy()
    absolute_fc = candidates["log2FC"].abs()
    fdr = candidates["q_value"].lt(0.05)
    nominal_only = candidates["p_value"].lt(0.05) & candidates["q_value"].ge(0.05)
    conditions = [
        fdr & absolute_fc.ge(1.0),
        fdr & absolute_fc.ge(0.585),
        fdr,
        nominal_only & absolute_fc.ge(0.585),
        nominal_only,
    ]
    candidates["candidate_origin"] = np.select(
        conditions,
        CANDIDATE_ORIGINS,
        default="Not candidate",
    )
    gene_names = candidates["Gene.name"].astype("string").str.strip().str.upper()
    candidates["analysis_gene"] = gene_names.mask(gene_names.eq(""), pd.NA)
    return candidates


def _split_genes(value: Any) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).replace(";", ",")
    return [gene.strip().upper() for gene in text.split(",") if gene.strip()]


def _representative_candidate_genes(
    candidates: pd.DataFrame,
    candidate_origins: Sequence[str],
    direction: str,
) -> pd.DataFrame:
    """Select genes globally, then retain representatives in one direction."""
    if direction not in {"increased", "decreased"}:
        raise ValueError(f"Unknown direction: {direction!r}")
    representatives = _representative_candidate_genes_all_directions(
        candidates, candidate_origins
    )
    direction_mask = (
        representatives["log2FC"].gt(0)
        if direction == "increased"
        else representatives["log2FC"].lt(0)
    )
    return representatives.loc[direction_mask].copy()


def _representative_candidate_genes_all_directions(
    candidates: pd.DataFrame,
    candidate_origins: Sequence[str],
) -> pd.DataFrame:
    """Select the strongest candidate row for each usable gene symbol."""
    representatives = candidates.loc[
        candidates["candidate_origin"].isin(candidate_origins)
        & candidates["analysis_gene"].notna()
    ].copy()
    if representatives.empty:
        return representatives
    representatives["_origin_priority"] = representatives["candidate_origin"].map(
        ORIGIN_PRIORITY
    )
    representatives = representatives.sort_values(
        ["analysis_gene", "_origin_priority", "q_value", "p_value", "row_position"],
        ascending=True,
        kind="stable",
    ).drop_duplicates("analysis_gene", keep="first")
    return representatives.drop(columns="_origin_priority")


def build_pathway_candidate_members(
    enrichment_results: pd.DataFrame,
    candidates: pd.DataFrame,
) -> pd.DataFrame:
    """Build long-format pathway-to-candidate evidence records."""
    rows: list[dict[str, Any]] = []
    for pathway in enrichment_results.itertuples(index=False):
        allowed_origins = CANDIDATE_ORIGINS
        candidate_subset = _representative_candidate_genes(
            candidates, allowed_origins, pathway.direction
        )
        member_genes = set(_split_genes(pathway.candidate_genes))
        members = candidate_subset.loc[
            candidate_subset["analysis_gene"].isin(member_genes)
        ]
        for _, member in members.iterrows():
            rows.append(
                {
                    "analysis_type": pathway.analysis_type,
                    "analysis_label": pathway.analysis_label,
                    "direction": pathway.direction,
                    "pathway_library": pathway.pathway_library,
                    "pathway": pathway.pathway,
                    "candidate_gene": member["analysis_gene"],
                    "row_position": member["row_position"],
                    "Gene.name": member["Gene.name"],
                    "Accession": member["Accession"],
                    "candidate_origin": member["candidate_origin"],
                    "log2FC": member["log2FC"],
                    "t_statistic": member["t_statistic"],
                    "p_value": member["p_value"],
                    "q_value": member["q_value"],
                }
            )
    return pd.DataFrame(rows, columns=MEMBER_COLUMNS)

