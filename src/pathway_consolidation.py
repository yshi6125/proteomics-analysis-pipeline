"""Consolidate redundant primary and exploratory pathways independently.

Workflow: significant primary pathways -> candidate-gene sets -> pairwise Jaccard
network -> Louvain communities -> deterministic representative pathways.

All member pathways are retained. This module performs no semantic or LLM-based
interpretation.
"""

from __future__ import annotations

from collections import Counter
from itertools import combinations
import json
from pathlib import Path
from typing import Any

import networkx as nx
import pandas as pd
from networkx.algorithms.community import louvain_communities


PROJECT_ROOT = Path(__file__).resolve().parent.parent
PATHWAY_RESULTS_DIR = PROJECT_ROOT / "results" / "pathway_analysis"
RESULTS_DIR = PROJECT_ROOT / "results" / "pathway_consolidation"

JACCARD_THRESHOLD = 0.25
MIN_CANDIDATE_GENES = 2
LOUVAIN_SEED = 42
PATHWAY_KEY_COLUMNS = ["analysis_type", "direction", "pathway_library", "pathway"]
REQUIRED_PATHWAY_COLUMNS = PATHWAY_KEY_COLUMNS + [
    "enrichment_adjusted_p_value",
    "overlap_size",
    "gene_ratio",
]
REQUIRED_MEMBER_COLUMNS = PATHWAY_KEY_COLUMNS + ["candidate_gene"]


def _load_analysis_inputs(
    pathway_input_path: Path,
    member_input_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not pathway_input_path.exists():
        raise FileNotFoundError(f"Pathway results not found: {pathway_input_path}")
    member_path = member_input_path
    if not member_path.exists():
        raise FileNotFoundError(f"Pathway candidate members not found: {member_path}")
    return pd.read_csv(pathway_input_path), pd.read_csv(member_path)


def _validate_inputs(pathways: pd.DataFrame, members: pd.DataFrame) -> None:
    missing_pathway = [column for column in REQUIRED_PATHWAY_COLUMNS if column not in pathways]
    missing_member = [column for column in REQUIRED_MEMBER_COLUMNS if column not in members]
    if missing_pathway:
        raise ValueError(f"Primary pathway results are missing columns: {missing_pathway}")
    if missing_member:
        raise ValueError(f"Pathway member results are missing columns: {missing_member}")
    if pathways.duplicated(PATHWAY_KEY_COLUMNS).any():
        raise ValueError("Primary pathway identities must be unique.")
    for column in ("enrichment_adjusted_p_value", "overlap_size", "gene_ratio"):
        converted = pd.to_numeric(pathways[column], errors="coerce")
        if converted.isna().any():
            raise ValueError(f"Primary pathway column {column!r} must be numeric and non-missing.")
        pathways[column] = converted


def _pathway_identity(row: pd.Series) -> tuple[str, str, str, str]:
    return tuple(str(row[column]) for column in PATHWAY_KEY_COLUMNS)


def reconstruct_pathway_gene_sets(
    pathways: pd.DataFrame,
    members: pd.DataFrame,
) -> dict[tuple[str, str, str, str], frozenset[str]]:
    """Reconstruct unique candidate-gene sets for all retained pathways."""
    retained_keys = {
        _pathway_identity(row) for _, row in pathways.loc[:, PATHWAY_KEY_COLUMNS].iterrows()
    }
    gene_sets: dict[tuple[str, str, str, str], set[str]] = {
        key: set() for key in retained_keys
    }
    for _, member in members.iterrows():
        key = _pathway_identity(member)
        if key not in gene_sets or pd.isna(member["candidate_gene"]):
            continue
        gene = str(member["candidate_gene"]).strip().upper()
        if gene:
            gene_sets[key].add(gene)
    missing_genes = [key for key, genes in gene_sets.items() if not genes]
    if missing_genes:
        raise ValueError(
            f"No candidate genes were reconstructed for {len(missing_genes)} retained pathways."
        )
    return {key: frozenset(genes) for key, genes in gene_sets.items()}


def calculate_jaccard(first: frozenset[str], second: frozenset[str]) -> float:
    """Calculate Jaccard similarity between two non-empty gene sets."""
    union = first | second
    return len(first & second) / len(union) if union else 0.0


def build_similarity_network(
    gene_sets: dict[tuple[str, str, str, str], frozenset[str]],
    jaccard_threshold: float = JACCARD_THRESHOLD,
) -> nx.Graph:
    """Build an undirected pathway network using thresholded Jaccard similarity."""
    if not 0 <= jaccard_threshold <= 1:
        raise ValueError("jaccard_threshold must be between 0 and 1.")
    graph = nx.Graph()
    nodes = sorted(gene_sets)
    graph.add_nodes_from(nodes)
    for first, second in combinations(nodes, 2):
        if first[0] != second[0] or first[1] != second[1]:
            continue
        similarity = calculate_jaccard(gene_sets[first], gene_sets[second])
        if similarity >= jaccard_threshold:
            graph.add_edge(first, second, weight=similarity, jaccard_similarity=similarity)
    return graph


def cluster_similarity_network(
    graph: nx.Graph,
    seed: int = LOUVAIN_SEED,
) -> list[set[tuple[str, str, str, str]]]:
    """Run Louvain on connected pathways and retain isolates as singletons."""
    isolates = sorted(nx.isolates(graph))
    connected_nodes = sorted(set(graph.nodes).difference(isolates))
    communities: list[set[tuple[str, str, str, str]]] = []
    if connected_nodes:
        communities.extend(
            set(community)
            for community in louvain_communities(
                graph.subgraph(connected_nodes).copy(),
                weight="weight",
                seed=seed,
            )
        )
    communities.extend({node} for node in isolates)
    return communities


def _rank_cluster_members(cluster: pd.DataFrame) -> pd.DataFrame:
    """Apply the deterministic representative-pathway ranking."""
    return cluster.sort_values(
        [
            "enrichment_adjusted_p_value",
            "overlap_size",
            "gene_ratio",
            "candidate_gene_count",
            "pathway_library",
            "pathway",
            "direction",
            "analysis_type",
        ],
        ascending=[True, False, False, False, True, True, True, True],
        kind="stable",
    )


def consolidate_analysis(
    comparison: str,
    analysis_type: str,
    pathway_input_path: Path,
    member_input_path: Path,
    output_dir: Path,
    jaccard_threshold: float = JACCARD_THRESHOLD,
    min_candidate_genes: int = MIN_CANDIDATE_GENES,
    adjusted_p_threshold: float = 0.05,
    louvain_seed: int = LOUVAIN_SEED,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Consolidate one pathway-analysis branch without crossing analysis types."""
    if not 0 < adjusted_p_threshold <= 1:
        raise ValueError("adjusted_p_threshold must be greater than 0 and at most 1.")
    if (
        isinstance(min_candidate_genes, bool)
        or not isinstance(min_candidate_genes, int)
        or min_candidate_genes < 1
    ):
        raise ValueError("min_candidate_genes must be a positive integer.")
    pathways, members = _load_analysis_inputs(pathway_input_path, member_input_path)
    _validate_inputs(pathways, members)
    observed_types = set(pathways["analysis_type"].astype(str))
    if observed_types != {analysis_type}:
        raise ValueError(
            f"Pathway input must contain only analysis_type={analysis_type!r}; "
            f"found {sorted(observed_types)}."
        )
    members = members.loc[members["analysis_type"].astype(str).eq(analysis_type)].copy()
    significant = pathways.loc[
        pathways["enrichment_adjusted_p_value"].lt(adjusted_p_threshold)
    ].copy()
    if significant.empty:
        raise ValueError(
            f"No {analysis_type} pathways passed adjusted p < {adjusted_p_threshold:g}."
        )
    gene_sets = reconstruct_pathway_gene_sets(significant, members)
    significant["_pathway_identity"] = significant.apply(_pathway_identity, axis=1)
    significant["candidate_gene_count"] = significant["_pathway_identity"].map(
        lambda identity: len(gene_sets[identity])
    )
    significant["candidate_genes_reconstructed"] = significant[
        "_pathway_identity"
    ].map(lambda identity: ";".join(sorted(gene_sets[identity])))
    excluded = significant.loc[
        significant["candidate_gene_count"].lt(min_candidate_genes)
    ].copy()
    retained = significant.loc[
        significant["candidate_gene_count"].ge(min_candidate_genes)
    ].copy()
    output_prefix = output_dir
    output_prefix.mkdir(parents=True, exist_ok=True)
    excluded_columns = [
        "analysis_type",
        "direction",
        "pathway_library",
        "pathway",
        "enrichment_adjusted_p_value",
        "candidate_gene_count",
        "candidate_genes_reconstructed",
    ]
    excluded_path = output_prefix / f"excluded_low_support_pathways_{comparison}.csv"
    if excluded.empty:
        excluded_path.unlink(missing_ok=True)
    else:
        excluded.loc[:, excluded_columns].to_csv(excluded_path, index=False)
    if retained.empty:
        raise ValueError(
            f"No {analysis_type} pathways had at least {min_candidate_genes} "
            "reconstructed candidate genes."
        )
    retained_gene_sets = {
        identity: gene_sets[identity] for identity in retained["_pathway_identity"]
    }
    graph = build_similarity_network(retained_gene_sets, jaccard_threshold)
    communities = cluster_similarity_network(graph, seed=louvain_seed)

    ranked_communities: list[tuple[tuple[str, str, str, str], set[tuple[str, str, str, str]]]] = []
    for community in communities:
        cluster = retained.loc[retained["_pathway_identity"].isin(community)]
        representative = _rank_cluster_members(cluster).iloc[0]["_pathway_identity"]
        ranked_communities.append((representative, community))
    ranked_communities.sort(key=lambda item: item[0])

    cluster_rows: list[pd.DataFrame] = []
    representative_rows: list[pd.Series] = []
    for cluster_number, (representative, community) in enumerate(ranked_communities, start=1):
        cluster = retained.loc[retained["_pathway_identity"].isin(community)].copy()
        cluster = _rank_cluster_members(cluster)
        cluster["cluster_id"] = f"C{cluster_number:03d}"
        cluster["cluster_size"] = len(cluster)
        cluster["is_singleton_cluster"] = len(cluster) == 1
        cluster["is_representative"] = cluster["_pathway_identity"].map(
            lambda identity: identity == representative
        )
        cluster["candidate_genes_reconstructed"] = cluster["_pathway_identity"].map(
            lambda identity: ";".join(sorted(retained_gene_sets[identity]))
        )
        representative_row = cluster.loc[cluster["is_representative"]].iloc[0].copy()
        representative_row["member_pathways"] = ";".join(
            cluster["pathway_library"].astype(str)
            + " | "
            + cluster["pathway"].astype(str)
        )
        representative_rows.append(representative_row)
        cluster_rows.append(cluster)

    pathway_clusters = pd.concat(cluster_rows, ignore_index=True).drop(
        columns="_pathway_identity"
    )
    consolidated_pathways = pd.DataFrame(representative_rows).drop(
        columns="_pathway_identity"
    ).reset_index(drop=True)
    pathway_clusters.to_csv(
        output_prefix / f"pathway_clusters_{comparison}.csv", index=False
    )
    for obsolete_name in (
        f"consolidated_pathways_{comparison}.csv",
        f"cluster_summary_{comparison}.csv",
    ):
        (output_prefix / obsolete_name).unlink(missing_ok=True)

    cluster_sizes = [len(community) for _, community in ranked_communities]
    size_distribution = dict(sorted(Counter(cluster_sizes).items()))
    singleton_count = sum(size == 1 for size in cluster_sizes)
    summary = {
        "comparison": comparison,
        "analysis_type": analysis_type,
        "jaccard_threshold": jaccard_threshold,
        "adjusted_p_threshold": adjusted_p_threshold,
        "min_candidate_genes": min_candidate_genes,
        "significant_pathways_before_support_filter": len(significant),
        "excluded_low_support_pathways": len(excluded),
        "pathways_retained_for_clustering": len(retained),
        "significant_pathways": len(significant),
        "clusters": len(ranked_communities),
        "singleton_clusters": singleton_count,
        "cluster_size_distribution": size_distribution,
        "network_edges": graph.number_of_edges(),
    }
    consolidation_metadata = {
        "comparison": comparison,
        "analysis_type": analysis_type,
        "jaccard_threshold": jaccard_threshold,
        "minimum_candidate_genes": min_candidate_genes,
        "adjusted_p_threshold": adjusted_p_threshold,
        "louvain_seed": louvain_seed,
        "significant_pathways_before_support_filter": len(significant),
        "excluded_low_support_pathways": len(excluded),
        "retained_pathways": len(retained),
        "number_of_clusters": len(ranked_communities),
    }
    metadata_path = output_prefix / f"consolidation_metadata_{comparison}.json"
    metadata_path.write_text(
        json.dumps(consolidation_metadata, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"\nAnalysis type: {analysis_type}")
    print(f"Significant pathways before candidate-gene filtering: {len(significant)}")
    print(
        f"Pathways removed because candidate_gene_count < {min_candidate_genes}: "
        f"{len(excluded)}"
    )
    print(f"Pathways retained for redundancy analysis: {len(retained)}")
    print(f"Number of clusters: {len(ranked_communities)}")
    print(f"Number of singleton clusters: {singleton_count}")
    print(f"Cluster-size distribution (size: count): {size_distribution}")
    print(f"Network edge count: {graph.number_of_edges()}")
    return pathway_clusters, consolidated_pathways, summary


def consolidate_pathways(
    comparison: str = "DM_vs_NDM",
    jaccard_threshold: float = JACCARD_THRESHOLD,
    min_candidate_genes: int = MIN_CANDIDATE_GENES,
    adjusted_p_threshold: float = 0.05,
    pathway_results_dir: Path = PATHWAY_RESULTS_DIR,
    results_dir: Path = RESULTS_DIR,
    louvain_seed: int = LOUVAIN_SEED,
) -> dict[str, tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]]:
    """Run independent primary and exploratory pathway consolidation."""
    member_input_path = (
        pathway_results_dir / f"pathway_candidate_members_{comparison}.csv"
    )
    configurations = {
        "primary": (
            "primary_fdr",
            pathway_results_dir / f"primary_fdr_pathways_{comparison}.csv",
        ),
        "exploratory": (
            "exploratory_all_candidates",
            pathway_results_dir
            / f"exploratory_all_candidates_pathways_{comparison}.csv",
        ),
    }
    outputs: dict[str, tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]] = {}
    for branch, (analysis_type, pathway_input_path) in configurations.items():
        outputs[branch] = consolidate_analysis(
            comparison=comparison,
            analysis_type=analysis_type,
            pathway_input_path=pathway_input_path,
            member_input_path=member_input_path,
            output_dir=results_dir / branch,
            jaccard_threshold=jaccard_threshold,
            min_candidate_genes=min_candidate_genes,
            adjusted_p_threshold=adjusted_p_threshold,
            louvain_seed=louvain_seed,
        )
    primary_summary = outputs["primary"][2]
    exploratory_summary = outputs["exploratory"][2]
    print("\nFinal consolidation comparison:")
    print(
        "Primary: "
        f"{primary_summary['pathways_retained_for_clustering']} retained pathways "
        f"-> {primary_summary['clusters']} clusters"
    )
    print(
        "Exploratory: "
        f"{exploratory_summary['pathways_retained_for_clustering']} retained pathways "
        f"-> {exploratory_summary['clusters']} clusters"
    )
    return outputs


def main() -> None:
    """Run default consolidation for the saved DM-versus-NDM pathway results."""
    consolidate_pathways()


if __name__ == "__main__":
    main()
