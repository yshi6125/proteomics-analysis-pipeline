"""Frozen-candidate Route B graph diagnostics and threshold sensitivity; no inference selection changes."""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from importlib.metadata import version
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from src.pathway_centric import consolidate_candidates
from src.pathway_redundancy import build_similarity_network, calculate_jaccard
from src.pathway_direction_diagnostics import describe_pathways, attach_fields, enrich_comparison
from src.pathway_review_io import ROOT, dataframe_hash, load_definitions, provenance, safe_label, sha256_file

THRESHOLDS = (.15,.20,.25,.30,.35,.40)
# Reporting-only membership lists, identified from the user's current cluster examples.
# They are never consulted by graph construction, Louvain or representative selection.
GROUP_CLUSTERS = {'respiration_OXPHOS_electron_transport': ['C001','C002','C003'],
                  'fatty_acid_metabolism': ['C004','C005'],
                  'coagulation_fibrinolysis_platelet_secretion': ['C008','C009','C011']}


def gene_memberships(candidates, members):
    sets = {}
    for row in candidates.itertuples():
        subset = members.loc[members.pathway_library.eq(row.pathway_library) & members.pathway.eq(row.pathway)]
        genes = frozenset(subset.analysis_gene.dropna())
        if len(genes) != row.pathway_size:
            raise ValueError('Frozen member table does not reconstruct the saved candidate size')
        sets[('pathway_centric',row.direction,row.pathway_library,row.pathway)] = genes
    return sets


def sensitivity_tables(candidates, members, canonical_clusters, thresholds=THRESHOLDS, seed=42, canonical_threshold=.25, report_group_clusters=None):
    """Use existing B consolidation unmodified; candidates and memberships are read-only."""
    if candidates.duplicated(['pathway_library','pathway']).any():
        raise ValueError('Candidate pathway identities must be unique')
    sets = gene_memberships(candidates, members)
    keys = sorted(sets)
    pairs = []
    for first,second in combinations(keys,2):
        overlap = sets[first] & sets[second]
        pairs.append(dict(first_library=first[2],first_pathway=first[3],second_library=second[2],second_pathway=second[3],
            first_direction=first[1],second_direction=second[1],same_direction=first[1]==second[1],
            intersection_count=len(overlap),union_count=len(sets[first]|sets[second]),
            measured_gene_jaccard=calculate_jaccard(sets[first],sets[second]),overlap_genes=';'.join(sorted(overlap))))
    pairwise = pd.DataFrame(pairs)
    summaries, assignments, representative_rows = [], [], []
    candidate_fingerprint=dataframe_hash(candidates)
    membership_fingerprint=dataframe_hash(members)
    graphs, runs = {}, {}
    for threshold in thresholds:
        graph = build_similarity_network(sets,threshold)
        graph.remove_edges_from([(u,v) for u,v,d in graph.edges(data=True) if d['weight']==0])
        clusters,reps = consolidate_candidates(candidates,members,jaccard_threshold=threshold,seed=seed)
        # Prove consolidation's returned rows still have the same identities and core evidence.
        before=candidates.sort_values(['pathway_library','pathway']).reset_index(drop=True)
        after=clusters.sort_values(['pathway_library','pathway']).reset_index(drop=True)
        pd.testing.assert_frame_equal(before,after[candidates.columns],check_exact=True)
        components=sorted((set(group) for group in nx.connected_components(graph)),key=lambda group:min(group))
        component_map={key:f'CC{i:03d}' for i,group in enumerate(components,1) for key in group}
        sizes=clusters.cluster_id.value_counts()
        if dataframe_hash(candidates)!=candidate_fingerprint or dataframe_hash(members)!=membership_fingerprint:
            raise RuntimeError('Frozen discovery inputs changed during a sensitivity run')
        summaries.append(dict(threshold=threshold,candidate_pathways=len(candidates),
            candidate_evidence_sha256=candidate_fingerprint,pathway_membership_sha256=membership_fingerprint,edges=graph.number_of_edges(),
            connected_components=len(components),clusters=len(sizes),singleton_clusters=int((sizes==1).sum()),
            largest_cluster=int(sizes.max()) if len(sizes) else 0,
            cluster_size_distribution=json.dumps(dict(sorted(Counter(map(int,sizes)).items()))),
            representative_pathways=json.dumps(reps[['pathway_library','pathway']].to_dict('records'),sort_keys=True)))
        for row in clusters.itertuples():
            key=('pathway_centric',row.direction,row.pathway_library,row.pathway)
            assignments.append(dict(threshold=threshold,pathway_library=row.pathway_library,pathway=row.pathway,
                measured_pathway_size=row.pathway_size,cluster_id=row.cluster_id,
                connected_component=component_map[key],is_representative=row.is_representative,
                mean_t=row.mean_t,p_value=row.p_value,q_value=row.q_value))
        copy=reps.copy();copy['threshold']=threshold;representative_rows.append(copy)
        graphs[threshold]=graph;runs[threshold]=clusters
    if canonical_threshold not in graphs:
        raise ValueError('Sensitivity thresholds must include the canonical threshold')
    canonical=runs[canonical_threshold]
    cols=list(canonical_clusters.columns)
    pd.testing.assert_frame_equal(canonical.sort_values(['pathway_library','pathway']).reset_index(drop=True)[cols],
        canonical_clusters.sort_values(['pathway_library','pathway']).reset_index(drop=True),check_exact=True)
    graph=graphs[canonical_threshold]
    diag=[]
    saved=canonical_clusters.set_index(['pathway_library','pathway'])
    for key in keys:
        peers=[other for other in keys if other!=key]
        # Maximum and strongest report raw overlaps, including opposite directions (not graph-eligible).
        ranked=sorted(peers,key=lambda other:(-calculate_jaccard(sets[key],sets[other]),other))
        same=sorted([other for other in peers if other[1]==key[1]],key=lambda other:(-calculate_jaccard(sets[key],sets[other]),other))
        strongest=ranked[0] if ranked else None
        neighbors=list(graph.neighbors(key));weights=[graph[key][other]['weight'] for other in neighbors]
        original=saved.loc[(key[2],key[3])]
        diag.append(dict(pathway_library=key[2],pathway=key[3],observed_direction=key[1],
            measured_pathway_size=len(sets[key]),cluster_id=original.cluster_id,is_representative=original.is_representative,
            graph_neighbors=len(neighbors),maximum_jaccard_similarity=calculate_jaccard(sets[key],sets[strongest]) if strongest else 0.,
            mean_jaccard_to_connected_neighbors=float(np.mean(weights)) if weights else np.nan,
            strongest_overlap_library=strongest[2] if strongest else '',strongest_overlap_pathway=strongest[3] if strongest else '',
            strongest_overlap_jaccard=calculate_jaccard(sets[key],sets[strongest]) if strongest else 0.,
            maximum_same_direction_jaccard=calculate_jaccard(sets[key],sets[same[0]]) if same else 0.))
    assignment=pd.DataFrame(assignments)
    current=assignment.loc[assignment.threshold.eq(canonical_threshold)].copy()
    current['component_size']=current.groupby('connected_component').pathway.transform('size')
    current['louvain_clusters_in_component']=current.groupby('connected_component').cluster_id.transform('nunique')
    stability=[]
    for first,second in combinations(thresholds,2):
        labels1=runs[first].set_index(['pathway_library','pathway']).cluster_id
        labels2=runs[second].set_index(['pathway_library','pathway']).cluster_id
        aligned=labels1.sort_index();other=labels2.reindex(aligned.index)
        stability.append(dict(first_threshold=first,second_threshold=second,
            adjusted_rand_index=float(adjusted_rand_score(aligned,other))))
    transitions=[]
    reference_sets={label:set(map(tuple,g[['pathway_library','pathway']].to_numpy())) for label,g in canonical_clusters.groupby('cluster_id')}
    for threshold,run in runs.items():
        targets={label:set(map(tuple,g[['pathway_library','pathway']].to_numpy())) for label,g in run.groupby('cluster_id')}
        for label,original in sorted(reference_sets.items()):
            destinations={name:sorted(original & group) for name,group in targets.items() if original & group}
            external=any(targets[name]-original for name in destinations)
            status=('unchanged' if not external else 'merged') if len(destinations)==1 else ('split_and_recombined' if external else 'split')
            transitions.append(dict(threshold=threshold,canonical_cluster=label,status=status,
                original_pathway_count=len(original),destination_count=len(destinations),
                destination_memberships=json.dumps(destinations,sort_keys=True)))
    grouped=[]
    for group,ids in (report_group_clusters or {}).items():
        selected={key for key in keys if saved.loc[(key[2],key[3])].cluster_id in ids}
        if pairwise.empty:continue
        subset=pairwise.loc[[('pathway_centric',r.first_direction,r.first_library,r.first_pathway) in selected and
                              ('pathway_centric',r.second_direction,r.second_library,r.second_pathway) in selected
                              for r in pairwise.itertuples()]].copy()
        subset['diagnostic_group']=group;grouped.append(subset)
    return dict(current_graph_diagnostics=pd.DataFrame(diag),jaccard_threshold_summary=pd.DataFrame(summaries),
        pathway_cluster_sensitivity=assignment,pathway_cluster_assignments_wide=assignment.pivot(
            index=['pathway_library','pathway'],columns='threshold',values='cluster_id').reset_index(),
        pathway_pairwise_jaccard=pairwise,connected_components_vs_louvain=current,
        representative_pathways_by_threshold=pd.concat(representative_rows,ignore_index=True),
        cluster_partition_stability=pd.DataFrame(stability),canonical_cluster_transitions=pd.DataFrame(transitions),
        diagnostic_group_pairwise_jaccard=pd.concat(grouped,ignore_index=True) if grouped else pd.DataFrame())


def run_diagnostics(comparison='DM_vs_NDM'):
    safe_label(comparison)
    directory=ROOT/'results/pathway_centric'/comparison
    differential=ROOT/'results/differential_abundance'/f'differential_abundance_{comparison}.csv'
    source_files=list(directory.glob('*.csv'))+[directory/'analysis_metadata.json',directory/'gene_sets.json',differential]
    a_dir=ROOT/'results/pathway_protein'/comparison
    comparison_dir=ROOT/'results/pathway_comparison'/comparison
    source_files += [p for d in [a_dir,comparison_dir] for p in d.glob('*') if p.is_file()]
    before=provenance(source_files)
    metadata=json.loads((directory/'analysis_metadata.json').read_text())
    proteins=pd.read_csv(differential)
    if dataframe_hash(proteins)!=metadata['differential_input_sha256'] or sha256_file(directory/'gene_sets.json')!=metadata['gene_sets_sha256']:
        raise ValueError('Frozen Route B input/snapshot hashes do not match')
    all_results=pd.read_csv(directory/'all_pathways.csv')
    candidates=pd.read_csv(directory/'pathway_candidates.csv')
    members=pd.read_csv(directory/'pathway_members.csv')
    canonical=pd.read_csv(directory/'pathway_clusters.csv')
    representatives=pd.read_csv(directory/'representative_pathways.csv')
    # Discovery/counts and membership remain frozen, not regenerated from a new cutoff.
    if len(candidates)!=metadata['number_of_candidate_pathways']:
        raise ValueError('Frozen candidate count differs from metadata')
    definitions=load_definitions(directory/'gene_sets.json')
    enriched,annotations,stats,audit=describe_pathways(all_results,proteins,definitions,
        permutations=metadata['permutations'],seed=metadata['permutation_seed'])
    if len(stats)!=metadata['number_of_eligible_genes'] or len(all_results)!=metadata['number_of_tested_pathways']:
        raise ValueError('Frozen eligible/tested counts disagree')
    # Verify every tested membership, not just the selected pathways.
    for row in all_results.itertuples():
        expected=set(definitions[row.pathway_library][row.pathway]) & set(stats.index)
        actual=set(members.loc[members.pathway_library.eq(row.pathway_library)&members.pathway.eq(row.pathway)].analysis_gene)
        if expected!=actual:raise ValueError('Frozen tested pathway membership mismatch')
    tables=sensitivity_tables(candidates,members,canonical,seed=metadata['louvain_seed'],canonical_threshold=metadata['jaccard_threshold'],
        report_group_clusters=GROUP_CLUSTERS if comparison=='DM_vs_NDM' else None)
    tables.update(all_pathways_with_direction=enriched,pathway_candidates_with_direction=attach_fields(candidates,annotations),
                  representative_pathways_with_direction=attach_fields(representatives,annotations))
    if (comparison_dir/'pathway_comparison.csv').exists():
        tables['pathway_comparison_with_direction']=enrich_comparison(pd.read_csv(comparison_dir/'pathway_comparison.csv'),annotations)
    target=directory/'redundancy_sensitivity';target.mkdir(exist_ok=True)
    for name,table in tables.items():table.to_csv(target/f'{name}.csv',index=False)
    info=dict(comparison=comparison,inputs=before,canonical_threshold=metadata['jaccard_threshold'],
        thresholds=list(THRESHOLDS),louvain_seed=metadata['louvain_seed'],candidate_pathways=len(candidates),
        tested_pathways=len(all_results),eligible_genes=len(stats),
        frozen_inference='mean t, p, q, candidate identities and measured memberships unchanged; no pathway reselection',
        direction='observed_direction = sign of mean_t; existing direction retained including non_directional for zero',
        null_tail='smaller inclusive empirical tail; equal tails balanced_or_tied; descriptive, not significance or activity',
        coherence='counts/fractions from unique eligible gene means; no effect on testing/clustering/representatives',
        tail_recovery='exact recorded size-seeded null replay, with frozen p-value agreement check; no p/q replacement',
        group_reporting='fixed original cluster lists from user examples; used only to subset pairwise reporting, never graph construction',
        graph_policy='measured-gene Jaccard; separate observed directions; original seeded Louvain and representative ranking',
        software_versions={p:version(p) for p in ['numpy','pandas','networkx','scikit-learn']},
        implementation_hashes=provenance([Path(__file__),ROOT/'src/pathway_direction_diagnostics.py']))
    if before!=provenance(source_files):raise RuntimeError('A frozen source changed during diagnostics')
    info['all_frozen_source_hashes_unchanged']=True
    (target/'redundancy_sensitivity_metadata.json').write_text(json.dumps(info,indent=2,sort_keys=True)+'\n')
    print(tables['jaccard_threshold_summary'][['threshold','edges','connected_components','clusters','singleton_clusters','largest_cluster']].to_string(index=False))
    print(f'Diagnostic views only: {target}; canonical files unchanged')
    return target,tables,info


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison',default='DM_vs_NDM')
    run_diagnostics(parser.parse_args().comparison)


if __name__=='__main__':main()
