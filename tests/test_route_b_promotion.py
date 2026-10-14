"""Verify promoted production against independent saved sensitivity and byte guards."""
import hashlib
import inspect
import json
import numpy as np
import pandas as pd
from src import pathway_centric as b
from src.pathway_review_io import ROOT
from provenance_helpers import PREVIOUS

CURRENT = ROOT/'results/pathway_centric/DM_vs_NDM'
OLD = PREVIOUS/'results/pathway_centric/DM_vs_NDM'
COMPARISON = ROOT/'results/pathway_comparison/DM_vs_NDM'


def test_production_defaults_and_strict_boundary():
    for function in [b.test_pathways,b.run_pathway_centric]:
        assert inspect.signature(function).parameters['permutations'].default == 100000
    result=pd.read_csv(CURRENT/'all_pathways.csv').head(3).copy()
    result['q_value']=[.049,.05,.051]
    assert len(b.select_candidates(result)) == 1


def test_production_matches_saved_100k_inference():
    actual=pd.read_csv(CURRENT/'all_pathways.csv')
    expected=pd.read_csv(CURRENT/'permutation_sensitivity/all_pathways_100000.csv')
    pd.testing.assert_frame_equal(actual[expected.columns],expected,check_exact=True)
    candidates=pd.read_csv(CURRENT/'pathway_candidates.csv')
    expected_candidates=pd.read_csv(CURRENT/'permutation_sensitivity/pathway_candidates_100000.csv')
    assert set(zip(candidates.pathway_library,candidates.pathway)) == set(zip(expected_candidates.pathway_library,expected_candidates.pathway))
    assert len(actual)==1617 and len(candidates)==28


def test_fixed_inputs_settings_and_protected_records():
    old=json.loads((OLD/'analysis_metadata.json').read_text())
    new=json.loads((CURRENT/'analysis_metadata.json').read_text())
    keys=['differential_input_sha256','gene_sets_sha256','duplicate_gene_rule','minimum_pathway_size','maximum_pathway_size','permutation_seed','jaccard_threshold','louvain_seed','multiple_testing_method','pathway_testing_method','protein_statistic','pathway_q_cutoff','pathway_p_cutoff','number_of_eligible_genes','number_of_eligible_proteins']
    assert {k:new[k] for k in keys}=={k:old[k] for k in keys}
    assert old['permutations']==10000 and new['permutations']==100000
    for filename in ['gene_sets.json','protein_eligibility.csv','gene_set_audit.csv','pathway_members.csv']:
        assert (OLD/filename).read_bytes()==(CURRENT/filename).read_bytes()
    manifest=json.loads((PREVIOUS/'MANIFEST.json').read_text())
    for path,digest in manifest['protected_route_a_and_validation_hashes'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest


def test_production_diagnostics_are_consistent_and_propagated():
    for filename in ['all_pathways.csv','pathway_candidates.csv','pathway_clusters.csv','representative_pathways.csv']:
        frame=pd.read_csv(CURRENT/filename)
        assert set(b.DIRECTION_COLUMNS).issubset(frame.columns)
        np.testing.assert_allclose(frame.p_value,np.minimum(1,2*np.minimum(frame.empirical_lower_tail_probability,frame.empirical_upper_tail_probability)),rtol=1e-12,atol=1e-15)
        assert (frame.n_positive_members+frame.n_negative_members+frame.n_zero_members).equals(frame.pathway_size)
    for filename in ['pathway_comparison.csv','pathway_protein_evidence.csv']:
        frame=pd.read_csv(COMPARISON/filename)
        assert {'route_b_'+field for field in b.DIRECTION_COLUMNS}.issubset(frame.columns)


def test_candidate_transitions_and_deterministic_outputs():
    previous=pd.read_csv(OLD/'pathway_candidates.csv');current=pd.read_csv(CURRENT/'pathway_candidates.csv')
    assert set(current.pathway)-set(previous.pathway)=={'Fatty Acid Beta-Oxidation Using acyl-CoA Dehydrogenase (GO:0033539)'}
    assert set(previous.pathway)-set(current.pathway)=={'Protein Secretion (GO:0009306)'}
    clusters=pd.read_csv(CURRENT/'pathway_clusters.csv');representatives=pd.read_csv(CURRENT/'representative_pathways.csv')
    assert clusters.cluster_id.nunique()==9 and (clusters.groupby('cluster_id').size()==1).sum()==4
    assert len(representatives)==9
    for _,group in clusters.groupby('cluster_id'):
        expected=group.sort_values(['q_value','p_value','pathway_library','pathway']).iloc[0]
        assert group.loc[group.is_representative,'pathway'].tolist()==[expected.pathway]
    comparison=pd.read_csv(COMPARISON/'pathway_comparison.csv')
    assert comparison['category'].value_counts().to_dict()=={'B_ONLY':23,'A_AND_B':5}
