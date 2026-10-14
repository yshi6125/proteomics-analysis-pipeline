import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src import pathway_centric as b
from src.pathway_direction_diagnostics import describe_pathways, enrich_comparison, FIELDS
from src.pathway_redundancy_sensitivity import sensitivity_tables, gene_memberships, THRESHOLDS
from src.pathway_review_io import ROOT, provenance


def fixture():
    p=pd.DataFrame({'Gene.name':list('ABCDEF')+['A'], 't_statistic':[4.,3.,-4.,-3.,0.,.1,2.]})
    defs={'GO':{'up':list('AB'),'up2':list('ABF'),'down':list('CD'),'zero':['E'],'mixed':list('AC')}}
    stats,_=b.prepare_proteins(p)
    all_results,members,_=b.test_pathways(stats,defs,min_size=1,permutations=99)
    candidates=b.select_candidates(all_results,q_cutoff=1.)
    clusters,reps=b.consolidate_candidates(candidates,members)
    return p,defs,all_results,members,candidates,clusters,reps


def test_direction_coherence_and_core_statistics_unchanged():
    p,defs,results,members,candidates,clusters,reps=fixture()
    before=results.copy(deep=True)
    enriched,ann,stats,audit=describe_pathways(results,p,defs,permutations=99)
    pd.testing.assert_frame_equal(enriched[results.columns],before,check_exact=True)
    pd.testing.assert_frame_equal(results,before,check_exact=True)
    by=enriched.set_index('pathway')
    assert by.loc['up','n_positive_members']==2 and by.loc['up','dominant_direction_fraction']==1
    assert by.loc['down','n_negative_members']==2
    assert by.loc['zero','observed_direction']=='zero' and by.loc['zero','direction']=='non_directional'
    assert by.loc['zero','n_zero_members']==1 and by.loc['zero','dominant_direction_fraction']==0
    assert by.loc['mixed','n_negative_members']==1 and by.loc['mixed','n_positive_members']==1
    assert by.loc['mixed','dominant_direction_fraction']==.5
    assert ((enriched.n_positive_members+enriched.n_negative_members+enriched.n_zero_members)==enriched.pathway_size).all()
    reproduced=np.minimum(1,2*enriched[['empirical_lower_tail_probability','empirical_upper_tail_probability']].min(axis=1))
    np.testing.assert_allclose(reproduced,results.p_value,rtol=1e-12)
    again,_,_,_=describe_pathways(results,p,defs,permutations=99)
    pd.testing.assert_frame_equal(enriched,again,check_exact=True)


def test_observed_sign_can_differ_from_null_tail():
    p=pd.DataFrame({'Gene.name':list('ABCDEF'),'t_statistic':[.1,.2,3.,4.,5.,6.]})
    defs={'GO':{'low':['A','B']}}
    stats,_=b.prepare_proteins(p)
    results,_,_=b.test_pathways(stats,defs,min_size=2,permutations=499)
    enriched,_,_,_=describe_pathways(results,p,defs,permutations=499)
    assert enriched.iloc[0].observed_direction=='positive'
    assert enriched.iloc[0].null_tail=='lower_than_random'


def test_balanced_tail_is_explicit():
    p=pd.DataFrame({'Gene.name':list('ABCD'),'t_statistic':[1.]*4})
    defs={'GO':{'constant':['A','B']}}
    stats,_=b.prepare_proteins(p)
    results,_,_=b.test_pathways(stats,defs,min_size=2,permutations=99)
    enriched,_,_,_=describe_pathways(results,p,defs,permutations=99)
    assert enriched.iloc[0].null_tail=='balanced_or_tied'
    assert enriched.iloc[0].p_value==1


def test_frozen_statistics_mismatch_rejected():
    p,defs,results,*_=fixture()
    altered=results.copy();altered.loc[0,'p_value']=.123456
    with pytest.raises(ValueError,match='reproduce'):describe_pathways(altered,p,defs,permutations=99)


def test_sweep_reproducible_and_every_core_input_identical():
    p,defs,results,members,candidates,clusters,reps=fixture()
    before=candidates.copy(deep=True);member_before=members.copy(deep=True)
    first=sensitivity_tables(candidates,members,clusters)
    second=sensitivity_tables(candidates,members,clusters)
    for name in first:pd.testing.assert_frame_equal(first[name],second[name],check_exact=True)
    pd.testing.assert_frame_equal(candidates,before,check_exact=True)
    pd.testing.assert_frame_equal(members,member_before,check_exact=True)
    expected=candidates[['pathway_library','pathway','mean_t','p_value','q_value']].sort_values(['pathway_library','pathway']).reset_index(drop=True)
    for threshold in THRESHOLDS:
        actual=first['pathway_cluster_sensitivity'].query('threshold == @threshold')[expected.columns].sort_values(['pathway_library','pathway']).reset_index(drop=True)
        pd.testing.assert_frame_equal(actual,expected,check_exact=True)
    assert set(first['jaccard_threshold_summary'].candidate_pathways)=={len(candidates)}
    assert first['jaccard_threshold_summary'].candidate_evidence_sha256.nunique()==1
    assert first['jaccard_threshold_summary'].pathway_membership_sha256.nunique()==1


def test_descriptive_fields_do_not_change_clusters_or_representatives():
    from src.pathway_direction_diagnostics import attach_fields
    p,defs,results,members,candidates,clusters,reps=fixture()
    _,ann,_,_=describe_pathways(results,p,defs,permutations=99)
    diagnostic_candidates=attach_fields(candidates,ann)
    other,_=b.consolidate_candidates(diagnostic_candidates,members)
    pd.testing.assert_frame_equal(other[clusters.columns],clusters,check_exact=True)
    comparison=pd.DataFrame({'pathway_library':['GO'],'pathway_identifier':['up'],'category':['B_ONLY']})
    enriched=enrich_comparison(comparison,ann)
    pd.testing.assert_frame_equal(enriched[comparison.columns],comparison)
    assert enriched.iloc[0].route_b_n_positive_members==2


def test_actual_frozen_graph_and_candidate_sweep():
    folder=ROOT/'archive/route_b_10k_canonical_20261014/results/pathway_centric/DM_vs_NDM'
    candidates=pd.read_csv(folder/'pathway_candidates.csv')
    members=pd.read_csv(folder/'pathway_members.csv')
    clusters=pd.read_csv(folder/'pathway_clusters.csv')
    before=provenance([p for p in folder.glob('*') if p.is_file()])
    tables=sensitivity_tables(candidates,members,clusters)
    current=tables['jaccard_threshold_summary'].query('threshold == 0.25').iloc[0]
    assert current.candidate_pathways==28 and current.clusters==11
    assert before==provenance([Path(p) for p in before])
    assignments=tables['pathway_cluster_sensitivity']
    for threshold in THRESHOLDS:
        assert assignments.query('threshold == @threshold').shape[0]==28
    assert gene_memberships(candidates,members)==gene_memberships(candidates,members.iloc[::-1])


def test_saved_diagnostic_views_keep_all_frozen_fields():
    folder=ROOT/'archive/route_b_10k_canonical_20261014/results/pathway_centric/DM_vs_NDM'
    diag=ROOT/'results/pathway_centric/DM_vs_NDM/redundancy_sensitivity'
    for original,enriched in [('all_pathways','all_pathways_with_direction'),('pathway_candidates','pathway_candidates_with_direction'),
                              ('representative_pathways','representative_pathways_with_direction')]:
        core=pd.read_csv(folder/(original+'.csv'))
        view=pd.read_csv(diag/(enriched+'.csv'))
        pd.testing.assert_frame_equal(view[core.columns],core,check_exact=True)
        assert set(FIELDS)<=set(view.columns)
    all_results=pd.read_csv(folder/'all_pathways.csv')
    metadata=json.loads((folder/'analysis_metadata.json').read_text())
    assert len(all_results)==metadata['number_of_tested_pathways']==1617
    assert len(pd.read_csv(folder/'pathway_candidates.csv'))==28
    manifest=json.loads((diag/'redundancy_sensitivity_metadata.json').read_text())
    assert manifest['eligible_genes']==1946 and manifest['all_frozen_source_hashes_unchanged']
    from provenance_helpers import assert_recorded_sources
    assert_recorded_sources(manifest['inputs'])
