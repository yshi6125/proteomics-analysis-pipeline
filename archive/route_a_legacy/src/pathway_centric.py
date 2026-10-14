"""Threshold-independent competitive mean-t gene-set testing (Route B)."""
from __future__ import annotations

import hashlib
import json
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from src.pathway_analysis import load_differential_abundance_results
from src.pathway_consolidation import (
    JACCARD_THRESHOLD, LOUVAIN_SEED, build_similarity_network, cluster_similarity_network,
)

ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = ROOT / 'results' / 'pathway_centric'
DEFAULT_LIBRARIES = ('GO_Biological_Process_2023', 'GO_Cellular_Component_2023')
RESULT_COLUMNS = ['analysis_route', 'pathway_library', 'pathway', 'pathway_size',
                  'annotated_size', 'mean_t', 'direction', 'mixed_signs', 'p_value', 'q_value']


def prepare_proteins(results):
    """Audit invalid rows; collapse repeated gene symbols by mean t, never significance."""
    required = {'Gene.name', 't_statistic'}
    if not required.issubset(results.columns):
        raise ValueError(f'Missing columns: {sorted(required - set(results.columns))}')
    audit = results.copy().reset_index(drop=True)
    audit['input_row'] = np.arange(len(audit))
    audit['analysis_gene'] = audit['Gene.name'].astype('string').str.strip().str.upper()
    values = pd.to_numeric(audit['t_statistic'], errors='coerce')
    invalid_gene = audit['analysis_gene'].isna() | audit['analysis_gene'].eq('')
    invalid_stat = ~np.isfinite(values)
    audit['exclusion_reason'] = np.select(
        [invalid_gene.fillna(True).to_numpy(dtype=bool), invalid_stat],
        ['missing_gene_symbol', 'invalid_t_statistic'], default='')
    audit['eligible'] = audit['exclusion_reason'].eq('')
    audit['t_statistic'] = values
    eligible = audit.loc[audit['eligible']]
    stats = eligible.groupby('analysis_gene', sort=True)['t_statistic'].mean()
    if stats.empty:
        raise ValueError('No eligible proteins with gene symbols and finite t-statistics.')
    return stats, audit


def load_gene_sets(libraries=DEFAULT_LIBRARIES, organism='Human'):
    """Use the same Enrichr/GSEApy infrastructure as Route A."""
    import gseapy as gp
    return {library: gp.get_library(name=library, organism=organism) for library in libraries}


def test_pathways(stats, gene_sets, min_size=5, max_size=500, permutations=10000, seed=42):
    """Two-sided equal-tail Monte Carlo test of mean signed t versus random gene sets.

    Sampling is without replacement within each permutation. Nulls are shared by
    measured set size; plus-one tail estimates prevent zero p-values. BH covers
    all tested terms across libraries. No protein p/q/FC thresholds are used.
    """
    for name, value in [('min_size', min_size), ('max_size', max_size), ('permutations', permutations)]:
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f'{name} must be a positive integer.')
    if max_size < min_size:
        raise ValueError('max_size must be at least min_size.')
    universe = set(stats.index)
    values = stats.to_numpy(dtype=float)
    nulls = {}
    rows, memberships, sets_audit = [], [], []
    for library, terms in sorted(gene_sets.items()):
        for term, genes in sorted(terms.items()):
            annotated = {str(g).strip().upper() for g in genes if str(g).strip()}
            measured = sorted(annotated & universe)
            n = len(measured)
            reason = 'below_min_size' if n < min_size else 'above_max_size' if n > max_size else ''
            sets_audit.append(dict(pathway_library=library, pathway=term, annotated_size=len(annotated),
                                   measured_size=n, exclusion_reason=reason))
            if reason:
                continue
            if n not in nulls:
                rng = np.random.default_rng(np.random.SeedSequence([seed, n]))
                nulls[n] = np.array([values[rng.choice(len(values), n, replace=False)].mean()
                                    for _ in range(permutations)])
            member_stats = stats.loc[measured]
            mean = float(member_stats.mean())
            null = nulls[n]
            lower = (1 + np.count_nonzero(null <= mean)) / (permutations + 1)
            upper = (1 + np.count_nonzero(null >= mean)) / (permutations + 1)
            p = min(1., 2 * min(lower, upper))
            direction = 'positive' if mean > 0 else 'negative' if mean < 0 else 'non_directional'
            rows.append(['pathway_centric', library, term, n, len(annotated), mean, direction,
                         bool((member_stats > 0).any() and (member_stats < 0).any()), p, np.nan])
            memberships.extend(dict(pathway_library=library, pathway=term, analysis_gene=g,
                                    gene_mean_t=float(stats[g])) for g in measured)
    results = pd.DataFrame(rows, columns=RESULT_COLUMNS)
    if not results.empty:
        results['q_value'] = multipletests(results['p_value'], method='fdr_bh')[1]
    members = pd.DataFrame(memberships, columns=['pathway_library', 'pathway', 'analysis_gene', 'gene_mean_t'])
    return results, members, pd.DataFrame(sets_audit, columns=['pathway_library', 'pathway', 'annotated_size', 'measured_size', 'exclusion_reason'])


def select_candidates(results, q_cutoff=.05, p_cutoff=1.):
    for value in (q_cutoff, p_cutoff):
        if not 0 < value <= 1:
            raise ValueError('Pathway cutoffs must be greater than zero and at most one.')
    return results.loc[results.q_value.le(q_cutoff) & results.p_value.le(p_cutoff)].copy()


def consolidate_candidates(candidates, members, jaccard_threshold=JACCARD_THRESHOLD, seed=LOUVAIN_SEED):
    sets = {}
    for row in candidates.itertuples():
        key = ('pathway_centric', row.direction, row.pathway_library, row.pathway)
        subset = members.loc[members.pathway_library.eq(row.pathway_library) & members.pathway.eq(row.pathway)]
        sets[key] = frozenset(subset.analysis_gene)
        if not sets[key]:
            raise ValueError(f'Missing measured members for {key}')
    graph = build_similarity_network(sets, jaccard_threshold)
    graph.remove_edges_from([(u, v) for u, v, data in graph.edges(data=True) if data['weight'] == 0])
    communities = sorted(cluster_similarity_network(graph, seed), key=lambda c: min(c))
    assignments = {}
    for number, community in enumerate(communities, 1):
        for key in community:
            assignments[(key[2], key[3])] = f'C{number:03d}'
    clusters = candidates.copy()
    clusters['cluster_id'] = [assignments[(r.pathway_library, r.pathway)] for r in clusters.itertuples()]
    clusters = clusters.sort_values(['q_value', 'p_value', 'pathway_library', 'pathway'], kind='stable')
    clusters['is_representative'] = ~clusters.cluster_id.duplicated()
    representatives = clusters.loc[clusters.is_representative].copy()
    return clusters, representatives


def run_pathway_centric(comparison='DM_vs_NDM', *, differential_results=None, gene_sets=None,
                         output_dir=OUTPUT_DIR, min_size=5, max_size=500, permutations=10000,
                         seed=42, q_cutoff=.05, p_cutoff=1., jaccard_threshold=JACCARD_THRESHOLD,
                         louvain_seed=LOUVAIN_SEED, libraries=DEFAULT_LIBRARIES, organism='Human'):
    if Path(comparison).name != comparison or comparison in {'.', '..'}:
        raise ValueError('comparison must be a filename-safe label.')
    results = load_differential_abundance_results(comparison) if differential_results is None else differential_results
    stats, audit = prepare_proteins(results)
    gene_sets = load_gene_sets(libraries, organism) if gene_sets is None else gene_sets
    all_results, members, set_audit = test_pathways(stats, gene_sets, min_size, max_size, permutations, seed)
    candidates = select_candidates(all_results, q_cutoff, p_cutoff)
    clusters, representatives = consolidate_candidates(candidates, members, jaccard_threshold, louvain_seed)
    # Keep every eligible original protein row, including repeated gene mappings.
    protein_members = members.merge(audit.loc[audit.eligible], on='analysis_gene', how='left')
    metadata = dict(comparison=comparison, analysis_route='pathway_centric', gene_set_database=sorted(gene_sets),
                    organism=organism, protein_statistic='t_statistic (existing fitted group effect)',
                    number_of_eligible_proteins=int(audit.eligible.sum()), number_of_eligible_genes=len(stats),
                    number_of_excluded_proteins=int((~audit.eligible).sum()),
                    duplicate_gene_rule='arithmetic mean of all eligible protein t-statistics',
                    pathway_testing_method='competitive mean-t gene sampling, two-sided equal-tail Monte Carlo',
                    minimum_pathway_size=min_size, maximum_pathway_size=max_size,
                    pathway_q_cutoff=q_cutoff, pathway_p_cutoff=p_cutoff, multiple_testing_method='BH across all libraries',
                    permutations=permutations, permutation_seed=seed, jaccard_threshold=jaccard_threshold,
                    louvain_seed=louvain_seed, number_of_tested_pathways=len(all_results),
                    number_of_candidate_pathways=len(candidates), number_of_redundancy_clusters=len(representatives),
                    direction_definition='sign of mean t; mixed_signs flags both signs; not pathway activity',
                    clustering_policy='separate directions; all measured pathway genes; q then p then library/name representative',
                    limitations='Ordinary, not moderated t. Competitive gene sampling assumes exchangeable genes; does not preserve inter-protein correlation. Signed mean may cancel mixed regulation.')
    output_dir = Path(output_dir) / comparison
    output_dir.mkdir(parents=True, exist_ok=True)
    tables = dict(all_pathways=all_results, pathway_candidates=candidates, pathway_members=protein_members,
                  protein_eligibility=audit, gene_set_audit=set_audit, pathway_clusters=clusters,
                  representative_pathways=representatives)
    for name, table in tables.items():
        table.to_csv(output_dir / f'{name}.csv', index=False)
    snapshot = {lib: {term: sorted({str(g).strip().upper() for g in genes}) for term, genes in sorted(terms.items())}
                for lib, terms in sorted(gene_sets.items())}
    snapshot_text = json.dumps(snapshot, sort_keys=True, indent=2) + '\n'
    (output_dir / 'gene_sets.json').write_text(snapshot_text)
    metadata['package_versions'] = {name: version(name) for name in ['numpy', 'pandas', 'statsmodels', 'networkx', 'gseapy']}
    metadata['gene_sets_sha256'] = hashlib.sha256(snapshot_text.encode()).hexdigest()
    metadata['differential_input_sha256'] = hashlib.sha256(results.to_csv(index=False).encode()).hexdigest()
    (output_dir / 'analysis_metadata.json').write_text(json.dumps(metadata, indent=2, sort_keys=True) + '\n')
    return tables, metadata
