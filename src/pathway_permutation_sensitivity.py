"""Route B permutation-count sensitivity using unchanged inference; no clustering runs."""
from __future__ import annotations

import argparse
from itertools import combinations
import json
from importlib.metadata import version
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from src import pathway_centric as b
from src.pathway_review_io import ROOT, dataframe_hash, provenance, safe_label, sha256_file

COUNTS = (10000, 50000, 100000)
KEYS = ['pathway_library','pathway']
CORE_FIXED = KEYS+['analysis_route','pathway_size','annotated_size','mean_t','direction','mixed_signs']


def candidate_keys(frame, cutoff=.05):
    """Strict requested criterion; no alteration to the production selector."""
    return set(map(tuple,frame.loc[frame.q_value.lt(cutoff),KEYS].to_numpy()))


def compare_runs(runs, cutoff=.05):
    """Compare identical testing families by exact library/term, using tie-aware ranks."""
    counts=sorted(runs)
    reference=runs[counts[0]].sort_values(KEYS).reset_index(drop=True)
    wide=reference[CORE_FIXED].copy()
    baseline=candidate_keys(reference,cutoff)
    summary=[];changes=[];correlations=[]
    aligned={}
    for count in counts:
        frame=runs[count].sort_values(KEYS).reset_index(drop=True)
        pd.testing.assert_frame_equal(frame[CORE_FIXED],reference[CORE_FIXED],check_exact=True)
        aligned[count]=frame
        candidates=candidate_keys(frame,cutoff)
        wide[f'p_{count}']=frame.p_value
        wide[f'q_{count}']=frame.q_value
        wide[f'candidate_{count}']=frame.q_value.lt(cutoff)
        wide[f'p_rank_{count}']=frame.p_value.rank(method='average')
        wide[f'q_rank_{count}']=frame.q_value.rank(method='average')
        intersection=len(candidates & baseline);union=len(candidates | baseline)
        pdelta=frame.p_value-reference.p_value;qdelta=frame.q_value-reference.q_value
        summary.append(dict(permutations=count,tested_pathways=len(frame),q_lt_005_candidates=len(candidates),
            overlap_with_10k=intersection,candidate_jaccard_with_10k=intersection/union if union else 1.,
            retained_fraction_of_10k=intersection/len(baseline) if baseline else np.nan,
            entering_vs_10k=len(candidates-baseline),leaving_vs_10k=len(baseline-candidates),
            exact_q_boundary_cases=int(frame.q_value.eq(cutoff).sum()),minimum_attainable_two_sided_p=2/(count+1),
            median_absolute_p_change_vs_10k=float(pdelta.abs().median()),max_absolute_p_change_vs_10k=float(pdelta.abs().max()),
            median_absolute_q_change_vs_10k=float(qdelta.abs().median()),max_absolute_q_change_vs_10k=float(qdelta.abs().max())))
        wide[f'p_change_{count}_vs_10k']=pdelta
        wide[f'q_change_{count}_vs_10k']=qdelta
        for key in sorted(candidates-baseline):changes.append(dict(permutations=count,transition='entering_vs_10k',pathway_library=key[0],pathway=key[1]))
        for key in sorted(baseline-candidates):changes.append(dict(permutations=count,transition='leaving_vs_10k',pathway_library=key[0],pathway=key[1]))
    for first,second in combinations(counts,2):
        one=aligned[first];two=aligned[second]
        ca=candidate_keys(one,cutoff);cb=candidate_keys(two,cutoff)
        selected=one.q_value.lt(cutoff) | two.q_value.lt(cutoff)
        row=dict(first_permutations=first,second_permutations=second,
            spearman_p_all=float(spearmanr(one.p_value,two.p_value).statistic),
            spearman_q_all=float(spearmanr(one.q_value,two.q_value).statistic),
            candidate_intersection=len(ca&cb),candidate_union=len(ca|cb),
            candidate_jaccard=len(ca&cb)/len(ca|cb) if ca|cb else 1.,
            only_in_first=len(ca-cb),only_in_second=len(cb-ca))
        for metric in ['p_value','q_value']:
            row['spearman_'+metric+'_candidate_union']=float(spearmanr(one.loc[selected,metric],two.loc[selected,metric]).statistic) if selected.sum()>1 else np.nan
        correlations.append(row)
    wide['candidate_in_any_run']=wide[[f'candidate_{count}' for count in counts]].any(axis=1)
    wide['candidate_in_all_runs']=wide[[f'candidate_{count}' for count in counts]].all(axis=1)
    wide['candidate_status_changes']=wide[[f'candidate_{count}' for count in counts]].nunique(axis=1).gt(1)
    wide['minimum_distance_to_q_cutoff']=wide[[f'q_{count}' for count in counts]].sub(cutoff).abs().min(axis=1)
    # Prespecified descriptive boundary window plus every changing identity.
    near=wide.loc[(wide[[f'q_{count}' for count in counts]].ge(.025)&wide[[f'q_{count}' for count in counts]].le(.075)).any(axis=1)
                  |wide.candidate_status_changes].sort_values(['minimum_distance_to_q_cutoff']+KEYS,kind='stable')
    nearest=wide.sort_values(['minimum_distance_to_q_cutoff']+KEYS,kind='stable').head(30)
    transitions=pd.DataFrame(changes,columns=['permutations','transition']+KEYS)
    return dict(pathway_permutation_comparison=wide,permutation_count_summary=pd.DataFrame(summary),
        candidate_transitions=transitions,pathway_rank_correlations=pd.DataFrame(correlations),
        near_q_cutoff=near,nearest_30_to_q_cutoff=nearest)


def run_sensitivity(comparison='DM_vs_NDM', counts=COUNTS):
    safe_label(comparison)
    if tuple(counts)!=COUNTS:raise ValueError('Production sensitivity uses exactly 10000, 50000, 100000')
    folder=ROOT/'results/pathway_centric'/comparison
    differential=ROOT/'results/differential_abundance'/f'differential_abundance_{comparison}.csv'
    sources=[p for d in [folder,ROOT/'results/pathway_protein'/comparison,ROOT/'results/pathway_comparison'/comparison]
             for p in d.glob('*') if p.is_file()]+[differential]
    before=provenance(sources)
    metadata=json.loads((folder/'analysis_metadata.json').read_text())
    if metadata['permutations']!=10000 or metadata['pathway_q_cutoff']!=.05 or metadata['pathway_p_cutoff']!=1.:
        raise ValueError('Unexpected canonical baseline settings; review before proceeding')
    proteins=pd.read_csv(differential)
    if dataframe_hash(proteins)!=metadata['differential_input_sha256'] or sha256_file(folder/'gene_sets.json')!=metadata['gene_sets_sha256']:
        raise ValueError('Frozen input/snapshot hash mismatch')
    # Read the literal local resource, exactly as the canonical offline call does.
    definitions=json.loads((folder/'gene_sets.json').read_text())
    stats,audit=b.prepare_proteins(proteins)
    if len(stats)!=metadata['number_of_eligible_genes'] or int(audit.eligible.sum())!=metadata['number_of_eligible_proteins']:
        raise ValueError('Eligibility changed from canonical')
    saved=pd.read_csv(folder/'all_pathways.csv').sort_values(KEYS).reset_index(drop=True)
    frozen_candidates=pd.read_csv(folder/'pathway_candidates.csv')
    saved_members=pd.read_csv(folder/'pathway_members.csv')[KEYS+['analysis_gene','gene_mean_t']].drop_duplicates()
    saved_members=saved_members.sort_values(KEYS+['analysis_gene']).reset_index(drop=True)
    target=folder/'permutation_sensitivity';target.mkdir(exist_ok=True)
    runs={};timings={};membership_reference=None
    for count in counts:
        print(f'Running unchanged Route B inference: {count:,} permutations; no clustering',flush=True)
        start=time.monotonic()
        result,members,set_audit=b.test_pathways(stats,definitions,
            min_size=metadata['minimum_pathway_size'],max_size=metadata['maximum_pathway_size'],
            permutations=count,seed=metadata['permutation_seed'])
        timings[count]=time.monotonic()-start
        ordered=result.sort_values(KEYS).reset_index(drop=True)
        pd.testing.assert_frame_equal(ordered[CORE_FIXED],saved[CORE_FIXED],check_exact=False,rtol=1e-12,atol=1e-14)
        member_order=members.sort_values(KEYS+['analysis_gene']).reset_index(drop=True)
        pd.testing.assert_frame_equal(member_order,saved_members,check_exact=False,rtol=1e-12,atol=1e-14)
        if membership_reference is None:membership_reference=member_order
        else:pd.testing.assert_frame_equal(member_order,membership_reference,check_exact=True)
        if count==10000:
            pd.testing.assert_frame_equal(ordered,saved,check_exact=False,rtol=1e-12,atol=1e-14)
            if candidate_keys(result)!=set(map(tuple,frozen_candidates[KEYS].to_numpy())):
                raise ValueError('Requested strict 10k candidates differ from frozen canonical candidates')
        runs[count]=result
        result.to_csv(target/f'all_pathways_{count}.csv',index=False)
        result.loc[result.q_value.lt(.05)].to_csv(target/f'pathway_candidates_{count}.csv',index=False)
        print(f'Completed {count:,}: {len(result)} tests; {int(result.q_value.lt(.05).sum())} strict-q candidates; {timings[count]:.1f}s',flush=True)
    tables=compare_runs(runs)
    for name,frame in tables.items():frame.to_csv(target/f'{name}.csv',index=False)
    stats.rename('gene_mean_t').reset_index().to_csv(target/'eligible_gene_statistics.csv',index=False)
    membership_reference.to_csv(target/'frozen_gene_memberships.csv',index=False)
    if before!=provenance(sources):raise RuntimeError('A canonical source changed during the analysis')
    info=dict(comparison=comparison,permutation_counts=list(counts),source_hashes=before,all_canonical_source_hashes_unchanged=True,
        eligible_genes=len(stats),eligible_proteins=int(audit.eligible.sum()),tested_pathways=len(saved),
        gene_statistics_sha256=dataframe_hash(stats.rename('gene_mean_t').reset_index()),
        gene_memberships_sha256=dataframe_hash(membership_reference),
        unchanged_canonical_settings=metadata,changed_parameter='permutations only',
        candidate_rule='q < 0.05 as requested; canonical uses <=, exact-boundary cases separately reported; 10k identities verified equal',
        seed_strategy='SeedSequence([42, measured_size]); same eligible gene ordering; longer runs extend each size-specific RNG stream',
        rank_comparison='Spearman on p/q with average ranks for ties; all terms and pairwise candidate unions',
        boundary_window='any run q in [0.025,0.075], plus all entering/leaving terms; nearest 30 also retained',
        uncertainty_scope='one fixed seed and nested streams; not an independent multi-seed replication study',
        redundancy='not invoked or evaluated',runtime_seconds=timings,
        implementation_hashes=provenance([Path(__file__),ROOT/'src/pathway_centric.py']),
        software_versions={name:version(name) for name in ['numpy','pandas','scipy','statsmodels']})
    (target/'permutation_sensitivity_metadata.json').write_text(json.dumps(info,indent=2,sort_keys=True)+'\n')
    print(tables['permutation_count_summary'].to_string(index=False),flush=True)
    return target,tables,info


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comparison',default='DM_vs_NDM')
    run_sensitivity(parser.parse_args().comparison)


if __name__=='__main__':main()
