"""Moderation validation: independent formula checks and production immutability."""
import json
import subprocess
import numpy as np
import pandas as pd
import pytest
from scipy.stats import t as t_distribution
from statsmodels.stats.multitest import multipletests
from src import pathway_centric as b
from src import pathway_moderated_t_sensitivity as m
from src.pathway_review_io import ROOT, provenance

FOLDER=ROOT/'results/pathway_centric/DM_vs_NDM/moderated_t_sensitivity'


def test_statistic_thresholds_and_tie_aware_ranking():
    f,summary=m.compare_statistics(pd.DataFrame({'ordinary_t':[4.,-2.,1.,0.], 'moderated_t':[2.,-3.,-1.,0.]}))
    assert summary['sign_changes']==1
    assert summary['large_statistic_changes']==3
    assert f.ordinary_absolute_t_rank.tolist()==[1.,2.,3.,4.]
    assert summary['pearson_t']<1


def test_pathway_alignment_boundary_and_direction():
    ordinary=pd.DataFrame({'pathway_library':['GO']*3,'pathway':['a','b','c'],'pathway_size':[5]*3,'mean_t':[1.,-2.,3.],
        'p_value':[.01,.02,.04],'q_value':[.049,.05,.051],'observed_direction':['positive','negative','positive']})
    for field in b.DIRECTION_COLUMNS[1:]:ordinary[field]=0.
    moderated=ordinary.iloc[::-1].copy()
    moderated['q_value']=[.04,.049,.051]
    moderated.loc[moderated.pathway.eq('a'),'mean_t']=-1.
    moderated.loc[moderated.pathway.eq('a'),'observed_direction']='negative'
    wide,summary=m.compare_pathways(ordinary,moderated)
    assert summary['ordinary_candidates']==1 and summary['moderated_candidates']==2
    assert summary['shared']==0 and summary['candidate_jaccard']==0
    assert summary['direction_agreement_all']==2
    assert wide.set_index('pathway').loc['b','transition']=='moderated_only'
    with pytest.raises(ValueError,match='Testing families differ'):
        m.compare_pathways(ordinary,moderated.iloc[:2])


def test_limma_posterior_variance_t_and_p_formulas():
    f=pd.read_csv(FOLDER/'limma_statistics.csv');h=pd.read_csv(FOLDER/'limma_hyperparameters.csv').iloc[0]
    expected_variance=(f.df_residual*f.residual_variance+h.df_prior*h.s2_prior)/(f.df_residual+h.df_prior)
    np.testing.assert_allclose(f.posterior_variance,expected_variance,rtol=1e-12)
    np.testing.assert_allclose(f.moderated_t,f.ordinary_t_limma*np.sqrt(f.residual_variance/f.posterior_variance),rtol=1e-12)
    np.testing.assert_allclose(f.moderated_p,2*t_distribution.sf(f.moderated_t.abs(),f.df_total),rtol=1e-12)
    np.testing.assert_allclose(f.moderated_q,multipletests(f.moderated_p,method='fdr_bh')[1],rtol=1e-12)
    assert len(f)==1965 and set(f.df_residual)=={17}
    assert np.isfinite(f.moderated_t).all()


def test_saved_design_and_protein_identifiers_match_canonical():
    design=pd.read_csv(FOLDER/'design_matrix.csv').set_index('sample_id')
    samples=pd.read_csv(FOLDER/'analysis_samples.csv').set_index('sample_id')
    assert list(design.columns)==['intercept','group__DM','batch__B2']
    np.testing.assert_array_equal(design.group__DM,samples.group.eq('DM').astype(float))
    np.testing.assert_array_equal(design.batch__B2,samples.batch.eq('B2').astype(float))
    assert np.linalg.matrix_rank(design)==3 and len(design)==20
    old=pd.read_csv(ROOT/'results/differential_abundance/differential_abundance_DM_vs_NDM.csv').sort_values('row_position').reset_index(drop=True)
    new=pd.read_csv(FOLDER/'moderated_differential_statistics.csv')
    pd.testing.assert_frame_equal(old[['row_position','Accession','Gene.name']],new[['row_position','Accession','Gene.name']])
    fit=pd.read_csv(FOLDER/'limma_statistics.csv')
    np.testing.assert_allclose(old.t_statistic,fit.ordinary_t_limma,rtol=1e-10,atol=1e-12)
    np.testing.assert_allclose(old.log2FC,fit.coefficient,rtol=1e-10,atol=1e-12)


def test_saved_pathway_statistics_use_exact_gene_means_and_bh():
    genes=pd.read_csv(FOLDER/'gene_statistic_comparison.csv').set_index('analysis_gene')
    members=pd.read_csv(FOLDER/'moderated_pathway_members.csv')
    results=pd.read_csv(FOLDER/'moderated_all_pathways.csv').set_index(m.KEYS)
    np.testing.assert_allclose(members.gene_mean_t,genes.loc[members.analysis_gene,'moderated_t'],rtol=1e-12)
    means=members.groupby(m.KEYS).gene_mean_t.mean()
    np.testing.assert_allclose(results.loc[means.index,'mean_t'],means,rtol=1e-12,atol=1e-14)
    np.testing.assert_allclose(results.p_value,np.minimum(1,2*np.minimum(results.empirical_lower_tail_probability,results.empirical_upper_tail_probability)),rtol=1e-12,atol=1e-15)
    np.testing.assert_allclose(results.q_value,multipletests(results.p_value,method='fdr_bh')[1],rtol=1e-12)
    actual=pd.read_csv(FOLDER/'moderated_pathway_candidates.csv')
    assert set(actual.pathway)==set(results.reset_index().loc[results.q_value.to_numpy()<.05,'pathway'])
    original=pd.read_csv(ROOT/'results/pathway_centric/DM_vs_NDM/pathway_members.csv')
    pd.testing.assert_frame_equal(members[m.KEYS+['analysis_gene']],original[m.KEYS+['analysis_gene']].drop_duplicates().reset_index(drop=True))


def test_production_is_byte_unchanged_and_no_clustering():
    meta=json.loads((FOLDER/'sensitivity_metadata.json').read_text())
    assert meta['protected_source_hashes']==provenance(meta['protected_source_hashes'])
    assert meta['protected_files_unchanged'] and not meta['clustering_run']
    assert meta['canonical_settings']['permutations']==100000
    assert meta['proteins_successfully_tested']==1965 and meta['proteins_excluded']==0
    assert meta['eligible_protein_rows']==1947 and meta['eligible_genes']==1946
    assert not list(FOLDER.glob('*cluster*'))


def test_validation_refuses_to_overwrite_existing_experiment():
    with pytest.raises(ValueError,match='Refusing to overwrite'):
        m.run_sensitivity()


def test_saved_null_samples_reproduce_disagreement_empirical_p_values():
    samples=np.load(FOLDER/'changing_pathway_nulls.npz')
    diagnostics=pd.read_csv(FOLDER/'changing_pathway_null_diagnostics.csv')
    wide=pd.read_csv(FOLDER/'pathway_comparison.csv').set_index(m.KEYS)
    for _,row in diagnostics.iterrows():
        pathway=wide.loc[(row.pathway_library,row.pathway)]
        size=int(pathway.ordinary_pathway_size)
        null=samples[f'{row.method}_size_{size}']
        assert len(null)==100000
        lower=(1+np.count_nonzero(null<=row.mean_t))/100001
        upper=(1+np.count_nonzero(null>=row.mean_t))/100001
        np.testing.assert_allclose(row.empirical_p,min(1,2*min(lower,upper)),rtol=1e-12,atol=1e-15)
        np.testing.assert_allclose(row.empirical_p,pathway[row.method+'_p'],rtol=1e-12,atol=1e-15)


def test_theme_anchors_and_status_disagreements_are_borderline():
    from src.pathway_moderated_t_review import summarize_themes
    wide=pd.read_csv(FOLDER/'pathway_comparison.csv')
    themes,_=summarize_themes(wide)
    assert len(themes)==6
    assert (themes.ordinary_candidates==themes.anchor_terms).all()
    assert (themes.moderated_candidates==themes.anchor_terms).all()
    assert (themes.direction_agreements==themes.anchor_terms).all()
    changes=pd.read_csv(FOLDER/'candidate_disagreements.csv')
    assert len(changes)==4 and set(changes.transition)=={'moderated_only'}
    assert changes.ordinary_near_cutoff.all() and changes.moderated_near_cutoff.all()
    assert (changes.members_with_large_t_change==0).all()
    assert wide.loc[wide.transition.ne('neither'),'direction_agreement'].all()


def test_disagreement_concentration_and_duplicate_gene_mean():
    genes=pd.read_csv(FOLDER/'gene_statistic_comparison.csv').set_index('analysis_gene')
    proteins=pd.read_csv(FOLDER/'protein_statistic_comparison.csv')
    duplicate=proteins.loc[proteins.analysis_gene.eq('HLA-A')]
    assert len(duplicate)==2
    np.testing.assert_allclose(genes.loc['HLA-A','moderated_t'],duplicate.moderated_t.mean(),rtol=1e-12)
    detail=pd.read_csv(FOLDER/'disagreement_gene_evidence.csv')
    summary=pd.read_csv(FOLDER/'candidate_disagreements.csv')
    for _,row in summary.iterrows():
        member=detail.loc[detail.pathway.eq(row.pathway)&detail.pathway_library.eq(row.pathway_library)]
        top=member.nlargest(3,'absolute_t_change')
        assert len(member)==row.pathway_size
        np.testing.assert_allclose(row.top3_share_absolute_member_changes,top.absolute_t_change.sum()/member.absolute_t_change.sum(),rtol=1e-12)
        np.testing.assert_allclose(row.mean_t_change,member.t_change.mean(),rtol=1e-12,atol=1e-14)
