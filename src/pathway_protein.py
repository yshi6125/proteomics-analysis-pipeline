"""Canonical Route A — Protein-Centric Analysis: Levels 1–5, local measured-universe ORA."""
from __future__ import annotations

import argparse
import json
from importlib.metadata import version
from pathlib import Path

import pandas as pd

from src import protein_evidence as annotations
from src.pathway_ora import run_local_ora, consolidate_protein_pathways
from src.pathway_review_io import (
    ROOT, dataframe_hash, guard_output, load_definitions, provenance, safe_label, sha256_file,
)

OUTPUT_DIR = ROOT/'results/pathway_protein'


def default_snapshot(comparison):
    """Use the established local B snapshot; never silently download or change the library."""
    path = ROOT/'results/pathway_centric'/safe_label(comparison)/'gene_sets.json'
    if not path.exists():
        raise FileNotFoundError('Route A requires local definitions. Supply --gene-sets PATH to a pinned gene-set JSON.')
    return path


def run_protein_centric(comparison='DM_vs_NDM', *, gene_sets_path=None, differential_path=None,
                        output_root=OUTPUT_DIR, adjusted_p_threshold=.05, min_candidate_genes=2,
                        jaccard_threshold=.25, louvain_seed=42):
    """Generate Route A independently from the complete DA table and local definitions."""
    safe_label(comparison)
    gene_sets_path = Path(gene_sets_path) if gene_sets_path is not None else default_snapshot(comparison)
    differential_path = Path(differential_path or annotations.DIFFERENTIAL_RESULTS_DIR/f'differential_abundance_{comparison}.csv')
    sources = [differential_path, gene_sets_path]
    target = guard_output(Path(output_root)/comparison,
        [ROOT/'archive', ROOT/'results/pathway_centric', ROOT/'results/pathway_comparison'], sources)
    definitions = load_definitions(gene_sets_path)
    proteins = pd.read_csv(differential_path)
    assigned = annotations.assign_candidate_origin(proteins)
    results, universe = run_local_ora(proteins, definitions)
    clusters = consolidate_protein_pathways(results, adjusted_p_threshold, min_candidate_genes,
                                   jaccard_threshold, louvain_seed)
    members = annotations.build_pathway_candidate_members(results, assigned)
    # Labels change after inference and clustering; no statistic or gene set changes.
    for frame in [results, clusters, members]:
        frame['analysis_type'] = 'protein_centric'
        frame['analysis_label'] = 'Route A — Protein-Centric Analysis'
    target.mkdir(parents=True, exist_ok=True)
    tables = dict(all_pathways=results, pathway_clusters=clusters,
        pathway_candidates=clusters.drop(columns=['cluster_id','cluster_size','is_representative','is_singleton_cluster']),
        representative_pathways=clusters.loc[clusters.is_representative], universe_genes=universe,
        candidate_members=members, protein_levels=assigned)
    for name, table in tables.items():
        table.to_csv(target/f'{name}.csv', index=False)
    (target/'gene_sets.json').write_text(json.dumps(definitions,indent=2,sort_keys=True)+'\n')
    metadata = dict(comparison=comparison, analysis_route='protein_centric', display_name='Route A',
        protein_evidence_scope='Levels 1–5', protein_candidate_levels=[1,2,3,4,5],
        enrichment_method='local ORA', enrichment_test='one-sided upper-tail hypergeometric (Fisher greater equivalent)',
        enrichment_background='measured proteome / eligible measured genes',
        universe_definition='unique valid symbols in complete successfully tested DA table intersect library annotation union',
        duplicate_rule='universe unique genes; candidates strongest origin, q, p, row position before direction',
        absent_library_gene_rule='excluded from that library universe and submitted list; retained in universe audit',
        multiple_testing='BH over all nonempty measured terms, libraries and both directions including zero hits',
        adjusted_p_threshold=adjusted_p_threshold,cutoff_operator='<',minimum_candidate_genes=min_candidate_genes,
        jaccard_threshold=jaccard_threshold,louvain_seed=louvain_seed,libraries=sorted(definitions),
        sources=provenance(sources),differential_input_sha256=dataframe_hash(proteins),
        gene_sets_sha256=sha256_file(target/'gene_sets.json'),
        implementation_hashes=provenance([Path(__file__), ROOT/'src/pathway_ora.py',ROOT/'src/protein_evidence.py',ROOT/'src/pathway_redundancy.py']),
        source_policy='local JSON only; supplied versioned library labels, not independently verified release dates',
        package_versions={p:version(p) for p in ['scipy','statsmodels','pandas','networkx']},
        tested_directional_terms=len(results),candidate_directional_terms=len(clusters),
        candidate_pathways=len(clusters[['pathway_library','pathway']].drop_duplicates()),
        clusters=int(clusters.cluster_id.nunique()),
        method_implementation='src.pathway_ora.run_local_ora; Levels 1–5 inference')
    (target/'analysis_metadata.json').write_text(json.dumps(metadata,indent=2,sort_keys=True)+'\n')
    print(f'Route A — Protein-Centric Analysis: {target}')
    return target,metadata


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison',default='DM_vs_NDM')
    parser.add_argument('--gene-sets',type=Path,help='Pinned local JSON; defaults to the existing Route B snapshot')
    parser.add_argument('--differential-results',type=Path)
    parser.add_argument('--adjusted-p-threshold',type=float,default=.05)
    parser.add_argument('--min-candidate-genes',type=int,default=2)
    parser.add_argument('--jaccard-threshold',type=float,default=.25)
    parser.add_argument('--louvain-seed',type=int,default=42)
    args=parser.parse_args()
    run_protein_centric(args.comparison,gene_sets_path=args.gene_sets,differential_path=args.differential_results,
        adjusted_p_threshold=args.adjusted_p_threshold,min_candidate_genes=args.min_candidate_genes,
        jaccard_threshold=args.jaccard_threshold,louvain_seed=args.louvain_seed)


if __name__=='__main__':
    main()
