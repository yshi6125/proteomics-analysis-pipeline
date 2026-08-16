"""Build separate structured evidence JSON files for pathway-analysis branches.

This module combines branch-specific consolidation outputs with pathway candidate
members. It performs no filtering, scoring, semantic interpretation, or LLM calls.
"""

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
CONSOLIDATION_RESULTS_DIR = PROJECT_ROOT / "results" / "pathway_consolidation"
PATHWAY_RESULTS_DIR = PROJECT_ROOT / "results" / "pathway_analysis"

PATHWAY_KEY_COLUMNS = ["analysis_type", "direction", "pathway_library", "pathway"]
REQUIRED_CLUSTER_COLUMNS = PATHWAY_KEY_COLUMNS + [
    "cluster_id",
    "cluster_size",
    "is_singleton_cluster",
    "enrichment_adjusted_p_value",
]
REQUIRED_CLUSTER_COLUMNS += ["is_representative", "gene_ratio"]
REQUIRED_MEMBER_COLUMNS = PATHWAY_KEY_COLUMNS + ["candidate_gene", "candidate_origin"]
ALL_CANDIDATE_LEVELS = ["Level 1", "Level 2", "Level 3", "Level 4", "Level 5"]
PIPELINE_VERSION = "v1"
LEVEL_PRIORITY = {level: position for position, level in enumerate(ALL_CANDIDATE_LEVELS, 1)}
BRANCH_CONFIGURATIONS = {
    "primary": {
        "analysis_type": "primary_fdr",
        "allowed_levels": ALL_CANDIDATE_LEVELS[:3],
    },
    "exploratory": {
        "analysis_type": "exploratory_all_candidates",
        "allowed_levels": ALL_CANDIDATE_LEVELS,
    },
}


def _load_branch_inputs(
    comparison: str,
    branch_dir: Path,
    pathway_results_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    cluster_path = branch_dir / f"pathway_clusters_{comparison}.csv"
    member_path = pathway_results_dir / f"pathway_candidate_members_{comparison}.csv"
    metadata_path = branch_dir / f"consolidation_metadata_{comparison}.json"
    for label, path in (
        ("Pathway clusters", cluster_path),
        ("Pathway candidate members", member_path),
        ("Consolidation metadata", metadata_path),
    ):
        if not path.exists():
            raise FileNotFoundError(f"{label} file not found: {path}")
    return (
        pd.read_csv(cluster_path),
        pd.read_csv(member_path),
        json.loads(metadata_path.read_text(encoding="utf-8")),
    )


def _validate_consolidation_metadata(
    metadata: dict[str, Any],
    comparison: str,
    analysis_type: str,
    number_of_clusters: int,
) -> None:
    required = [
        "comparison",
        "analysis_type",
        "jaccard_threshold",
        "minimum_candidate_genes",
        "adjusted_p_threshold",
        "louvain_seed",
        "number_of_clusters",
    ]
    missing = [field for field in required if field not in metadata]
    if missing:
        raise ValueError(f"Consolidation metadata is missing fields: {missing}")
    if metadata["comparison"] != comparison:
        raise ValueError("Consolidation metadata comparison does not match the request.")
    if metadata["analysis_type"] != analysis_type:
        raise ValueError("Consolidation metadata analysis_type does not match the branch.")
    if int(metadata["number_of_clusters"]) != number_of_clusters:
        raise ValueError("Consolidation metadata cluster count does not match the audit CSV.")


def _validate_columns(frame: pd.DataFrame, required: list[str], label: str) -> None:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{label} is missing required columns: {missing}")


def _validate_branch_inputs(
    clusters: pd.DataFrame,
    members: pd.DataFrame,
    analysis_type: str,
    allowed_levels: list[str],
) -> None:
    _validate_columns(clusters, REQUIRED_CLUSTER_COLUMNS, "Pathway clusters")
    _validate_columns(members, REQUIRED_MEMBER_COLUMNS, "Pathway candidate members")
    if clusters.empty:
        raise ValueError("Branch pathway clusters must not be empty.")
    observed = set(clusters["analysis_type"].astype(str))
    if observed != {analysis_type}:
        raise ValueError(
            f"clusters must contain only analysis_type={analysis_type!r}; "
            f"found {sorted(observed)}."
        )
    if clusters["cluster_id"].isna().any():
        raise ValueError("cluster_id must be non-missing.")
    for cluster_id, group in clusters.groupby("cluster_id", sort=False):
        if group["analysis_type"].nunique() != 1 or group["direction"].nunique() != 1:
            raise ValueError(
                f"Cluster {cluster_id!r} contains multiple analysis types or directions."
            )
        if int(group["cluster_size"].iloc[0]) != len(group):
            raise ValueError(f"Cluster size is inconsistent for {cluster_id!r}.")
        if group["is_representative"].sum() != 1:
            raise ValueError(f"Cluster {cluster_id!r} must have one representative.")
    for column in ("enrichment_adjusted_p_value", "gene_ratio"):
        clusters[column] = pd.to_numeric(clusters[column], errors="coerce")
        if clusters[column].isna().any():
            raise ValueError(f"Evidence column {column!r} must be numeric.")


def _pathway_keys(frame: pd.DataFrame) -> pd.MultiIndex:
    return pd.MultiIndex.from_frame(frame.loc[:, PATHWAY_KEY_COLUMNS].astype(str))


def _strongest_gene_levels(
    member_rows: pd.DataFrame,
    allowed_levels: list[str],
) -> dict[str, str]:
    """Assign each unique gene its strongest allowed candidate level."""
    allowed = set(allowed_levels)
    levels: dict[str, str] = {}
    for row in member_rows.itertuples(index=False):
        gene = str(row.candidate_gene).strip().upper()
        origin = str(row.candidate_origin).strip()
        if not gene or origin not in allowed:
            continue
        if gene not in levels or LEVEL_PRIORITY[origin] < LEVEL_PRIORITY[levels[gene]]:
            levels[gene] = origin
    return levels


def build_branch_evidence(
    comparison: str,
    branch: str,
    analysis_type: str,
    allowed_levels: list[str],
    consolidation_results_dir: Path = CONSOLIDATION_RESULTS_DIR,
    pathway_results_dir: Path = PATHWAY_RESULTS_DIR,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build and save evidence for one independent consolidation branch."""
    branch_dir = consolidation_results_dir / branch
    clusters, members, consolidation_metadata = _load_branch_inputs(
        comparison, branch_dir, pathway_results_dir
    )
    members = members.loc[
        members["analysis_type"].astype(str).eq(analysis_type)
    ].copy()
    _validate_branch_inputs(
        clusters, members, analysis_type, allowed_levels
    )
    _validate_consolidation_metadata(
        consolidation_metadata,
        comparison,
        analysis_type,
        clusters["cluster_id"].nunique(),
    )
    cluster_pathway_keys = set(_pathway_keys(clusters).tolist())
    relevant_members = members.loc[
        _pathway_keys(members).isin(cluster_pathway_keys)
    ].copy()
    invalid_origins = sorted(
        set(relevant_members["candidate_origin"].dropna().astype(str)).difference(
            allowed_levels
        )
    )
    if invalid_origins:
        raise ValueError(
            f"Retained pathway members for {analysis_type!r} contain disallowed "
            f"candidate levels: {invalid_origins}"
        )
    cluster_evidence: list[dict[str, Any]] = []
    all_genes: set[str] = set()
    for cluster_id, cluster in clusters.groupby("cluster_id", sort=True):
        representative = cluster.loc[cluster["is_representative"]].iloc[0]
        cluster_keys = set(_pathway_keys(cluster).tolist())
        cluster_members = relevant_members.loc[
            _pathway_keys(relevant_members).isin(cluster_keys)
        ]
        gene_levels = _strongest_gene_levels(cluster_members, allowed_levels)
        level_genes = {
            level: sorted(
                gene for gene, origin in gene_levels.items() if origin == level
            )
            for level in allowed_levels
        }
        candidate_genes = [
            gene for level in allowed_levels for gene in level_genes[level]
        ]
        if not candidate_genes:
            raise ValueError(
                f"No allowed candidate genes were reconstructed for cluster {cluster_id!r}."
            )
        all_genes.update(candidate_genes)
        member_pathways = [
            {
                "library": str(row.pathway_library),
                "name": str(row.pathway),
            }
            for row in cluster.sort_values(
                ["pathway_library", "pathway"], kind="stable"
            ).itertuples(index=False)
        ]
        entry: dict[str, Any] = {
            "cluster_id": str(cluster_id),
            "analysis_type": analysis_type,
            "direction": str(representative["direction"]),
            "representative_pathway": {
                "library": str(representative["pathway_library"]),
                "name": str(representative["pathway"]),
            },
            "member_pathways": member_pathways,
            "cluster_size": len(member_pathways),
            "genes": [
                {"gene": gene, "level": LEVEL_PRIORITY[gene_levels[gene]]}
                for gene in candidate_genes
            ],
            "n_genes": len(candidate_genes),
            "level_counts": {
                str(level_number): len(level_genes.get(f"Level {level_number}", []))
                for level_number in range(1, 6)
            },
        }
        entry.update(
            {
                "best_adjusted_p": float(
                    cluster["enrichment_adjusted_p_value"].min()
                ),
                "representative_gene_ratio": float(representative["gene_ratio"]),
                "is_singleton": bool(representative["is_singleton_cluster"]),
            }
        )
        cluster_evidence.append(entry)

    evidence = {
        "metadata": {
            "comparison": comparison,
            "analysis_type": analysis_type,
            "pipeline_version": PIPELINE_VERSION,
            "jaccard_threshold": consolidation_metadata["jaccard_threshold"],
            "minimum_candidate_genes": consolidation_metadata[
                "minimum_candidate_genes"
            ],
            "adjusted_p_threshold": consolidation_metadata[
                "adjusted_p_threshold"
            ],
            "louvain_seed": consolidation_metadata["louvain_seed"],
            "generated_at": datetime.now(timezone.utc).isoformat(),
        },
        "number_of_clusters": len(cluster_evidence),
        "clusters": cluster_evidence,
    }
    output_path = branch_dir / f"pathway_evidence_{comparison}.json"
    output_path.write_text(
        json.dumps(evidence, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    gene_count_distribution = dict(
        sorted(
            Counter(
                cluster["n_genes"]
                for cluster in cluster_evidence
            ).items()
        )
    )
    cluster_candidate_level_occurrence_distribution = {
        level: sum(
            cluster["level_counts"][str(level_number)]
            for cluster in cluster_evidence
        )
        for level_number, level in enumerate(allowed_levels, start=1)
    }
    branch_gene_levels = _strongest_gene_levels(relevant_members, allowed_levels)
    unique_candidate_level_distribution = {
        level: sum(origin == level for origin in branch_gene_levels.values())
        for level in allowed_levels
    }
    summary = {
        "analysis_type": analysis_type,
        "total_clusters": len(cluster_evidence),
        "total_unique_candidate_genes": len(all_genes),
        "singleton_clusters": sum(
            cluster["is_singleton"] for cluster in cluster_evidence
        ),
        "unique_genes_per_cluster_distribution": gene_count_distribution,
        "cluster_candidate_level_occurrence_distribution": (
            cluster_candidate_level_occurrence_distribution
        ),
        "unique_candidate_level_distribution": unique_candidate_level_distribution,
    }
    print(f"\nAnalysis type: {analysis_type}")
    print(f"Total clusters: {summary['total_clusters']}")
    print(f"Total unique candidate genes: {summary['total_unique_candidate_genes']}")
    print(f"Singleton clusters: {summary['singleton_clusters']}")
    print(
        "Distribution of unique genes per cluster (genes: clusters): "
        f"{gene_count_distribution}"
    )
    print(
        "Cluster-level candidate-level occurrences "
        f"(genes may repeat across clusters): {cluster_candidate_level_occurrence_distribution}"
    )
    print(
        "Unique candidate-level distribution across the branch: "
        f"{unique_candidate_level_distribution}"
    )
    return evidence, summary


def build_pathway_evidence(
    comparison: str = "DM_vs_NDM",
    consolidation_results_dir: Path = CONSOLIDATION_RESULTS_DIR,
    pathway_results_dir: Path = PATHWAY_RESULTS_DIR,
) -> dict[str, tuple[dict[str, Any], dict[str, Any]]]:
    """Build independent primary and exploratory pathway evidence outputs."""
    outputs = {
        branch: build_branch_evidence(
            comparison=comparison,
            branch=branch,
            analysis_type=configuration["analysis_type"],
            allowed_levels=list(configuration["allowed_levels"]),
            consolidation_results_dir=consolidation_results_dir,
            pathway_results_dir=pathway_results_dir,
        )
        for branch, configuration in BRANCH_CONFIGURATIONS.items()
    }
    primary = outputs["primary"][1]
    exploratory = outputs["exploratory"][1]
    print("\nFinal evidence comparison:")
    print(
        f"Primary: {primary['total_clusters']} clusters, "
        f"{primary['total_unique_candidate_genes']} unique genes"
    )
    print(
        f"Exploratory: {exploratory['total_clusters']} clusters, "
        f"{exploratory['total_unique_candidate_genes']} unique genes"
    )
    return outputs


def main() -> None:
    """Build default primary and exploratory DM-versus-NDM evidence JSON files."""
    build_pathway_evidence()


if __name__ == "__main__":
    main()
