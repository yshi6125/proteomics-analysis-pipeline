"""Levels 1–5 local ORA with a per-library measured-proteome universe."""
from __future__ import annotations
import pandas as pd
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from src import protein_evidence as evidence
from src.pathway_redundancy import JACCARD_THRESHOLD, LOUVAIN_SEED, _rank_cluster_members, build_similarity_network, cluster_similarity_network

def run_local_ora(proteins, gene_sets):
    """Test all nonempty measured terms, including zero-overlap terms, in both directions."""
    candidates = evidence.assign_candidate_origin(proteins)
    measured = set(candidates.analysis_gene.dropna())
    analysis_type, origins = 'protein_centric', evidence.CANDIDATE_ORIGINS
    rows, audits = [], []
    for library, terms in sorted(gene_sets.items()):
        annotations = {term: {str(g).strip().upper() for g in genes} for term, genes in terms.items()}
        library_genes = set().union(*annotations.values()) if annotations else set()
        universe = measured & library_genes
        for gene in sorted(measured):
            audits.append(dict(pathway_library=library, analysis_gene=gene, in_universe=gene in universe))
        for direction in ('increased', 'decreased'):
            reps = evidence._representative_candidate_genes(candidates, origins, direction)
            submitted = set(reps.analysis_gene) & universe
            for term, genes in sorted(annotations.items()):
                members = genes & universe
                if not members:
                    continue
                overlap = sorted(submitted & members)
                n, k, draw, hit = len(universe), len(members), len(submitted), len(overlap)
                p = float(hypergeom.sf(hit - 1, n, k, draw)) if draw else 1.
                rows.append(dict(analysis_type=analysis_type,
                    analysis_label='Route A — Protein-Centric Analysis', direction=direction,
                    pathway_library=library, pathway=term, enrichment_p_value=p,
                    enrichment_adjusted_p_value=1., overlap_count=hit, overlap_size=len(genes),
                    overlap=f'{hit}/{len(genes)}', input_candidate_count=draw,
                    gene_ratio=hit / draw if draw else 0., candidate_genes=';'.join(overlap),
                    measured_pathway_size=k, universe_size=n,
                    table_a=hit, table_b=draw-hit, table_c=k-hit, table_d=n-draw-k+hit))
    columns = evidence.ENRICHMENT_COLUMNS + ['measured_pathway_size', 'universe_size',
                                          'table_a', 'table_b', 'table_c', 'table_d']
    results = pd.DataFrame(rows, columns=columns)
    if not results.empty:
        results['enrichment_adjusted_p_value'] = multipletests(results.enrichment_p_value, method='fdr_bh')[1]
    audit = pd.DataFrame(audits, columns=['pathway_library', 'analysis_gene', 'in_universe'])
    return results, audit


def consolidate_protein_pathways(results, adjusted_p_threshold=.05, min_candidate_genes=2,
                         jaccard_threshold=JACCARD_THRESHOLD, louvain_seed=LOUVAIN_SEED):
    """Preserve A's support, representative rank and graph logic, allowing empty outputs."""
    if not 0 < adjusted_p_threshold <= 1:
        raise ValueError('adjusted_p_threshold must be in (0,1]')
    if isinstance(min_candidate_genes, bool) or not isinstance(min_candidate_genes, int) or min_candidate_genes < 1:
        raise ValueError('min_candidate_genes must be a positive integer')
    all_results = results.copy()
    all_results['candidate_gene_count'] = all_results.candidate_genes.map(lambda s: len(evidence._split_genes(s)))
    all_results['candidate_genes_reconstructed'] = all_results.candidate_genes
    selected = all_results.loc[all_results.enrichment_adjusted_p_value.lt(adjusted_p_threshold)
        & all_results.candidate_gene_count.ge(min_candidate_genes)].copy()
    sets = {(r.analysis_type, r.direction, r.pathway_library, r.pathway):
            frozenset(evidence._split_genes(r.candidate_genes)) for r in selected.itertuples()}
    graph = build_similarity_network(sets, jaccard_threshold)
    # Zero similarity cannot supply a meaningful weighted connection at threshold zero.
    graph.remove_edges_from([(u, v) for u, v, d in graph.edges(data=True) if d['weight'] == 0])
    communities = cluster_similarity_network(graph, louvain_seed)
    ranked = []
    for community in communities:
        subset = selected.loc[[tuple(row[c] for c in ['analysis_type', 'direction', 'pathway_library', 'pathway'])
                               in community for _, row in selected.iterrows()]]
        representative = _rank_cluster_members(subset).iloc[0]
        identity = tuple(representative[c] for c in ['analysis_type', 'direction', 'pathway_library', 'pathway'])
        ranked.append((identity, community))
    assignments = {}
    for number, (representative, community) in enumerate(sorted(ranked), 1):
        for key in community:
            assignments[key] = (f'C{number:03d}', len(community), key == representative)
    selected['cluster_id'] = [assignments[(r.analysis_type, r.direction, r.pathway_library, r.pathway)][0]
                              for r in selected.itertuples()]
    selected['cluster_size'] = [assignments[(r.analysis_type, r.direction, r.pathway_library, r.pathway)][1]
                                for r in selected.itertuples()]
    selected['is_representative'] = [assignments[(r.analysis_type, r.direction, r.pathway_library, r.pathway)][2]
                                     for r in selected.itertuples()]
    selected['is_singleton_cluster'] = selected.cluster_size.eq(1)
    return selected.reset_index(drop=True)

