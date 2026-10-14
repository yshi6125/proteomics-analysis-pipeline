"""Shared Jaccard network, seeded Louvain and deterministic representative ranking."""
from __future__ import annotations
from itertools import combinations
import networkx as nx
import pandas as pd
from networkx.algorithms.community import louvain_communities

JACCARD_THRESHOLD = 0.25


LOUVAIN_SEED = 42


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

