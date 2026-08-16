"""Candidate-based pathway enrichment for differential-abundance results.

Workflow: complete differential-abundance result table -> mutually exclusive
candidate origins -> primary FDR-supported enrichment -> exploratory enrichment
-> candidate-member evidence -> pathway priorities -> separate summary figures.

This module does not perform QC, preprocessing, differential-abundance modeling,
or multiple-testing correction of protein-level results. A future preranked GSEA
workflow can use ``t_statistic`` from every successfully tested protein without
requiring protein-level significance.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DIFFERENTIAL_RESULTS_DIR = PROJECT_ROOT / "results" / "differential_abundance"
RESULTS_DIR = PROJECT_ROOT / "results" / "pathway_analysis"
FIGURES_DIR = PROJECT_ROOT / "figures" / "pathway_analysis"

DEFAULT_PATHWAY_LIBRARIES = (
    "GO_Biological_Process_2023",
    "Reactome_2022",
    "KEGG_2021_Human",
)
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
RANKING_METRIC = "t_statistic"
ANALYSIS_LABELS = {
    "primary_fdr": "PRIMARY — FDR-supported pathway analysis",
    "exploratory_all_candidates": (
        "EXPLORATORY — includes nominal p < 0.05 candidates"
    ),
}

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
PRIORITY_COLUMNS = ENRICHMENT_COLUMNS + [
    "Level1_count",
    "Level2_count",
    "Level3_count",
    "Level4_count",
    "Level5_count",
    "strongest_candidate_origin",
    "pathway_priority",
]


def _ensure_output_directories() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)


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


def _parse_overlap(value: Any) -> tuple[int, int | float]:
    text = str(value)
    if "/" in text:
        numerator, denominator = text.split("/", maxsplit=1)
        try:
            return int(float(numerator)), int(float(denominator))
        except ValueError:
            pass
    try:
        return int(float(value)), np.nan
    except (TypeError, ValueError):
        return 0, np.nan


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


def _default_enrichr_runner(**kwargs: Any) -> Any:
    try:
        import gseapy as gp
    except ImportError as error:
        raise ImportError(
            "Pathway enrichment requires gseapy. Install project dependencies "
            "with `pip install -r requirements.txt`."
        ) from error
    return gp.enrichr(**kwargs)


def _normalize_enrichr_results(
    raw_results: pd.DataFrame,
    analysis_type: str,
    direction: str,
    input_candidate_count: int,
) -> pd.DataFrame:
    """Normalize Enrichr output and add candidate-list overlap context.

    ``overlap_count`` is the number of submitted candidate genes found in the
    pathway; ``overlap_size`` is the total number of genes annotated to that
    pathway. ``gene_ratio`` uses the submitted candidate count as its denominator.
    """
    if raw_results is None or raw_results.empty:
        return pd.DataFrame(columns=ENRICHMENT_COLUMNS)
    aliases = {
        "Gene_set": "pathway_library",
        "Term": "pathway",
        "P-value": "enrichment_p_value",
        "Adjusted P-value": "enrichment_adjusted_p_value",
        "Genes": "candidate_genes",
    }
    missing = [column for column in aliases if column not in raw_results]
    if missing:
        raise ValueError(f"GSEApy Enrichr output is missing columns: {missing}")
    normalized = raw_results.rename(columns=aliases).copy()
    for column in ("enrichment_p_value", "enrichment_adjusted_p_value"):
        normalized[column] = pd.to_numeric(normalized[column], errors="coerce")
        if normalized[column].isna().any() or not np.isfinite(normalized[column]).all():
            raise ValueError(f"GSEApy Enrichr output contains invalid {column} values.")
    overlaps = normalized.get("Overlap", pd.Series(0, index=normalized.index)).map(
        _parse_overlap
    )
    normalized["overlap_count"] = overlaps.map(lambda value: value[0]).astype(int)
    normalized["overlap_size"] = pd.array(
        overlaps.map(lambda value: value[1]), dtype="Int64"
    )
    if normalized["overlap_count"].gt(input_candidate_count).any():
        raise ValueError(
            "GSEApy Enrichr overlap_count exceeds the submitted candidate count."
        )
    normalized["overlap"] = overlaps.map(
        lambda value: (
            f"{value[0]}/{int(value[1])}"
            if np.isfinite(value[1])
            else str(value[0])
        )
    )
    normalized["input_candidate_count"] = input_candidate_count
    normalized["gene_ratio"] = (
        normalized["overlap_count"] / input_candidate_count
    )
    normalized["analysis_type"] = analysis_type
    normalized["analysis_label"] = ANALYSIS_LABELS[analysis_type]
    normalized["direction"] = direction
    normalized["candidate_genes"] = normalized["candidate_genes"].map(
        lambda value: ";".join(_split_genes(value))
    )
    return normalized.loc[:, ENRICHMENT_COLUMNS]


def run_candidate_enrichment(
    candidates: pd.DataFrame,
    candidate_origins: Sequence[str],
    analysis_type: str,
    pathway_libraries: Sequence[str] = DEFAULT_PATHWAY_LIBRARIES,
    organism: str = "Human",
    enrichr_runner: Callable[..., Any] | None = None,
) -> pd.DataFrame:
    """Run separate increased/decreased Enrichr analyses for selected origins."""
    invalid_origins = sorted(set(candidate_origins).difference(CANDIDATE_ORIGINS))
    if invalid_origins:
        raise ValueError(f"Unknown candidate origins: {invalid_origins}")
    if analysis_type not in ANALYSIS_LABELS:
        raise ValueError(f"Unknown analysis_type: {analysis_type!r}")
    if not pathway_libraries:
        raise ValueError("At least one pathway library is required.")
    runner = enrichr_runner or _default_enrichr_runner
    outputs: list[pd.DataFrame] = []
    for direction in ("increased", "decreased"):
        representatives = _representative_candidate_genes(
            candidates, candidate_origins, direction
        )
        genes = representatives["analysis_gene"].astype(str).tolist()
        if not genes:
            continue
        enrichment = runner(
            gene_list=genes,
            gene_sets=list(pathway_libraries),
            organism=organism.strip().lower(),
            outdir=None,
            no_plot=True,
        )
        raw_results = enrichment.results if hasattr(enrichment, "results") else enrichment
        outputs.append(
            _normalize_enrichr_results(
                raw_results,
                analysis_type,
                direction,
                input_candidate_count=len(genes),
            )
        )
    if not outputs:
        return pd.DataFrame(columns=ENRICHMENT_COLUMNS)
    return pd.concat(outputs, ignore_index=True)


def build_pathway_candidate_members(
    enrichment_results: pd.DataFrame,
    candidates: pd.DataFrame,
) -> pd.DataFrame:
    """Build long-format pathway-to-candidate evidence records."""
    rows: list[dict[str, Any]] = []
    for pathway in enrichment_results.itertuples(index=False):
        allowed_origins = (
            CANDIDATE_ORIGINS[:3]
            if pathway.analysis_type == "primary_fdr"
            else CANDIDATE_ORIGINS
        )
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


def add_pathway_priorities(
    enrichment_results: pd.DataFrame,
    pathway_members: pd.DataFrame,
) -> pd.DataFrame:
    """Add evidence counts and priority without modifying enrichment statistics."""
    if enrichment_results.empty:
        return pd.DataFrame(columns=PRIORITY_COLUMNS)
    keys = [
        "analysis_type",
        "analysis_label",
        "direction",
        "pathway_library",
        "pathway",
    ]
    if pathway_members.empty:
        counts = pd.DataFrame(columns=keys + CANDIDATE_ORIGINS)
    else:
        counts = (
            pathway_members.groupby(keys + ["candidate_origin"], dropna=False)
            .size()
            .unstack(fill_value=0)
            .reset_index()
        )
    summary = enrichment_results.merge(counts, on=keys, how="left")
    for level in CANDIDATE_ORIGINS:
        if level not in summary.columns:
            summary[level] = 0
        else:
            summary[level] = summary[level].fillna(0).astype(int)
    present = summary.loc[:, CANDIDATE_ORIGINS].gt(0)
    summary["strongest_candidate_origin"] = present.idxmax(axis=1).where(
        present.any(axis=1), "No mapped candidate"
    )
    summary["pathway_priority"] = summary["strongest_candidate_origin"].map(
        lambda origin: (
            f"Priority {ORIGIN_PRIORITY[origin]}"
            + (" — Exploratory" if ORIGIN_PRIORITY[origin] >= 4 else "")
            if origin in ORIGIN_PRIORITY
            else "Unassigned"
        )
    )
    summary = summary.rename(
        columns={origin: f"Level{number}_count" for number, origin in enumerate(CANDIDATE_ORIGINS, 1)}
    )
    summary["_priority_number"] = summary["strongest_candidate_origin"].map(
        ORIGIN_PRIORITY
    ).fillna(99)
    summary = summary.sort_values(
        ["analysis_type", "_priority_number", "enrichment_adjusted_p_value", "overlap_count"],
        ascending=[True, True, True, False],
        kind="stable",
    ).drop(columns="_priority_number")
    return summary.loc[:, PRIORITY_COLUMNS].reset_index(drop=True)


def create_pathway_summary_plot(
    pathway_summary: pd.DataFrame,
    analysis_type: str,
    comparison: str,
    top_n: int = 20,
    adjusted_p_threshold: float = 0.05,
) -> Path:
    """Plot top pathways for one clearly identified analysis type."""
    if isinstance(top_n, bool) or not isinstance(top_n, (int, np.integer)) or top_n < 1:
        raise ValueError("top_n must be a positive integer.")
    if not 0 < adjusted_p_threshold <= 1:
        raise ValueError("adjusted_p_threshold must be greater than 0 and at most 1.")
    _ensure_output_directories()
    selected = pathway_summary.loc[
        pathway_summary["analysis_type"].eq(analysis_type)
        & pathway_summary["enrichment_adjusted_p_value"].lt(adjusted_p_threshold)
    ].head(top_n).copy()
    figure, axis = plt.subplots(figsize=(12, max(6, 0.42 * max(len(selected), 1) + 2)))
    if selected.empty:
        axis.text(
            0.5,
            0.5,
            f"No pathways passed adjusted p < {adjusted_p_threshold:g}",
            ha="center",
            va="center",
        )
        axis.set_axis_off()
    else:
        adjusted = selected["enrichment_adjusted_p_value"].astype(float)
        positive = adjusted[adjusted.gt(0)]
        floor = positive.min() / 10 if not positive.empty else np.finfo(float).tiny
        strength = -np.log10(adjusted.clip(lower=floor, upper=1.0))
        priority_numbers = selected["strongest_candidate_origin"].map(ORIGIN_PRIORITY).fillna(6)
        colors = plt.cm.viridis_r((priority_numbers - 1) / 5)
        labels = (
            selected["pathway"].astype(str)
            + " [" + selected["direction"].astype(str) + "]"
            + " — " + selected["pathway_priority"].astype(str)
        )
        positions = np.arange(len(selected))
        axis.barh(positions, strength, color=colors)
        axis.set_yticks(positions, labels=labels, fontsize=8)
        axis.invert_yaxis()
        axis.set_xlabel("-log10(adjusted pathway p-value)")
    display_name = ANALYSIS_LABELS[analysis_type]
    axis.set_title(f"{display_name}\n{comparison}")
    figure.tight_layout()
    filename = (
        f"primary_fdr_pathways_{comparison}.png"
        if analysis_type == "primary_fdr"
        else f"exploratory_all_candidates_pathways_{comparison}.png"
    )
    output_path = FIGURES_DIR / filename
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(figure)
    return output_path


def _print_analysis_report(
    label: str,
    origins: Sequence[str],
    candidates: pd.DataFrame,
    summary: pd.DataFrame,
    adjusted_p_threshold: float,
    top_n: int = 10,
) -> None:
    selected = _representative_candidate_genes_all_directions(candidates, origins)
    significant = summary.loc[summary["enrichment_adjusted_p_value"].lt(adjusted_p_threshold)]
    print(f"\n{label}:")
    print(f"Candidates used: {', '.join(origins)}")
    print(f"Total unique candidates: {len(selected)}")
    for direction in ("increased", "decreased"):
        submitted_count = len(
            _representative_candidate_genes(candidates, origins, direction)
        )
        print(f"Submitted unique genes ({direction}): {submitted_count}")
    print(f"Significant pathways (adjusted p < {adjusted_p_threshold}): {len(significant)}")
    for priority in range(1, max(ORIGIN_PRIORITY[origin] for origin in origins) + 1):
        count = significant["pathway_priority"].str.startswith(
            f"Priority {priority}", na=False
        ).sum()
        print(f"Priority {priority} pathways: {count}")
    for pathway in significant.head(top_n).itertuples(index=False):
        counts = ", ".join(
            f"L{level}={getattr(pathway, f'Level{level}_count')}" for level in range(1, 6)
        )
        print(
            f"- {pathway.pathway} | adjusted p={pathway.enrichment_adjusted_p_value:.4g} "
            f"| overlap={pathway.overlap} | gene ratio={pathway.gene_ratio:.3f} "
            f"| {pathway.pathway_priority} | {counts} | {pathway.candidate_genes}"
        )


def run_pathway_analysis(
    comparison: str = "DM_vs_NDM",
    pathway_libraries: Sequence[str] = DEFAULT_PATHWAY_LIBRARIES,
    organism: str = "Human",
    adjusted_p_threshold: float = 0.05,
    top_n_plot: int = 20,
    differential_results: pd.DataFrame | None = None,
    enrichr_runner: Callable[..., Any] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Run primary and exploratory candidate-based pathway enrichment."""
    if not 0 < adjusted_p_threshold <= 1:
        raise ValueError("adjusted_p_threshold must be greater than 0 and at most 1.")
    _ensure_output_directories()
    differential_results = (
        load_differential_abundance_results(comparison)
        if differential_results is None
        else differential_results
    )
    candidates = assign_candidate_origin(differential_results)
    primary = run_candidate_enrichment(
        candidates,
        CANDIDATE_ORIGINS[:3],
        "primary_fdr",
        pathway_libraries,
        organism,
        enrichr_runner,
    )
    exploratory = run_candidate_enrichment(
        candidates,
        CANDIDATE_ORIGINS,
        "exploratory_all_candidates",
        pathway_libraries,
        organism,
        enrichr_runner,
    )
    combined = pd.concat([primary, exploratory], ignore_index=True)
    members = build_pathway_candidate_members(combined, candidates)
    priority_summary = add_pathway_priorities(combined, members)

    primary.to_csv(RESULTS_DIR / f"primary_fdr_pathways_{comparison}.csv", index=False)
    exploratory.to_csv(
        RESULTS_DIR / f"exploratory_all_candidates_pathways_{comparison}.csv", index=False
    )
    members.to_csv(RESULTS_DIR / f"pathway_candidate_members_{comparison}.csv", index=False)
    priority_summary.to_csv(
        RESULTS_DIR / f"pathway_priority_summary_{comparison}.csv", index=False
    )
    primary_plot = create_pathway_summary_plot(
        priority_summary,
        "primary_fdr",
        comparison,
        top_n_plot,
        adjusted_p_threshold,
    )
    exploratory_plot = create_pathway_summary_plot(
        priority_summary,
        "exploratory_all_candidates",
        comparison,
        top_n_plot,
        adjusted_p_threshold,
    )

    gene_candidates = _representative_candidate_genes_all_directions(
        candidates, CANDIDATE_ORIGINS
    )
    counts = gene_candidates["candidate_origin"].value_counts()
    candidate_mask = candidates["candidate_origin"].isin(CANDIDATE_ORIGINS)
    missing_gene_symbol_count = int(
        (candidate_mask & candidates["analysis_gene"].isna()).sum()
    )
    print("\nCandidate origins:")
    for level in CANDIDATE_ORIGINS:
        suffix = " only" if level != "Level 1" else ""
        print(f"{level}{suffix}: {int(counts.get(level, 0))}")
    print(
        "Candidate proteins excluded from Enrichr because no usable gene symbol "
        f"was available: {missing_gene_symbol_count}"
    )
    _print_analysis_report(
        "PRIMARY ANALYSIS",
        CANDIDATE_ORIGINS[:3],
        candidates,
        priority_summary.loc[priority_summary["analysis_type"].eq("primary_fdr")],
        adjusted_p_threshold,
    )
    _print_analysis_report(
        "EXPLORATORY ANALYSIS",
        CANDIDATE_ORIGINS,
        candidates,
        priority_summary.loc[
            priority_summary["analysis_type"].eq("exploratory_all_candidates")
        ],
        adjusted_p_threshold,
    )
    analysis_summary = {
        "comparison": comparison,
        "organism": organism,
        "pathway_libraries": list(pathway_libraries),
        "candidate_origin_counts": {
            level: int(counts.get(level, 0)) for level in CANDIDATE_ORIGINS
        },
        "candidate_proteins_excluded_missing_gene_symbol": missing_gene_symbol_count,
        "submitted_gene_counts": {
            analysis_type: {
                direction: len(
                    _representative_candidate_genes(candidates, origins, direction)
                )
                for direction in ("increased", "decreased")
            }
            for analysis_type, origins in (
                ("primary_fdr", CANDIDATE_ORIGINS[:3]),
                ("exploratory_all_candidates", CANDIDATE_ORIGINS),
            )
        },
        "primary_candidates": len(
            _representative_candidate_genes_all_directions(
                candidates, CANDIDATE_ORIGINS[:3]
            )
        ),
        "exploratory_all_candidates": len(gene_candidates),
        "primary_significant_pathways": int(primary["enrichment_adjusted_p_value"].lt(adjusted_p_threshold).sum()),
        "exploratory_significant_pathways": int(exploratory["enrichment_adjusted_p_value"].lt(adjusted_p_threshold).sum()),
        "primary_plot": str(primary_plot),
        "exploratory_plot": str(exploratory_plot),
        "future_preranked_gsea_metric": RANKING_METRIC,
    }
    return candidates, primary, exploratory, members, priority_summary, analysis_summary


def main() -> None:
    """Run the default pathway analysis from the saved DM-versus-NDM results."""
    run_pathway_analysis()


if __name__ == "__main__":
    main()
