import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from src import pathway_analysis as a
from src import pathway_centric as b
from src import pathway_comparison as comparison
from src import pathway_measured_ora as measured
from src.pathway_review_io import dataframe_hash, guard_output, load_definitions, provenance


def proteins():
    return pd.DataFrame(dict(row_position=range(9), Accession=[f'P{i}' for i in range(9)],
        **{'Gene.name': ['A','B','C','D','E','F','G','A','UNANNOTATED']},
        log2FC=[2.,.8,.2,-2.,-.7,-.2,0.,.1,3.], t_statistic=[4.,3.,2.,-4.,-3.,-2.,.1,.2,5.],
        p_value=[.001,.002,.003,.01,.02,.03,.9,.8,.001],
        q_value=[.01,.02,.03,.1,.2,.3,.99,.9,.01]))


SETS = {'GO_2023': {'both (GO:0000001)': ['A','B','X'], 'a_only (GO:0000002)': ['B','C'],
                   'b_only (GO:0000003)': ['D','E','F'], 'neutral (GO:0000004)': ['G']}}


def saved_tables():
    p = proteins()
    all_a, _ = measured.measured_enrichment(p, SETS, 'primary')
    all_a['enrichment_adjusted_p_value'] = .5
    all_a['enrichment_p_value'] = .2
    select = all_a.pathway.str.contains('both|a_only') & all_a.direction.eq('increased')
    all_a.loc[select,'enrichment_adjusted_p_value'] = .01
    all_a.loc[select,'enrichment_p_value'] = .001
    ac = measured.consolidate_measured(all_a)
    am = a.build_pathway_candidate_members(all_a, a.assign_candidate_origin(p))
    stats, _ = b.prepare_proteins(p)
    all_b, _, _ = b.test_pathways(stats, SETS, min_size=1, permutations=99)
    all_b['q_value'] = .5
    all_b.loc[all_b.pathway.str.contains('both|b_only'),'q_value'] = .01
    bcandidate = b.select_candidates(all_b)
    membership = pd.DataFrame([dict(pathway_library=lib,pathway=term,analysis_gene=gene)
        for lib,terms in SETS.items() for term,genes in terms.items() for gene in genes if gene in stats.index])
    bc, _ = b.consolidate_candidates(bcandidate,membership)
    return p, all_a, ac, am, all_b, bcandidate, bc


def build(tables=None, **kwargs):
    return comparison.build_comparison(*(tables or saved_tables()), SETS, ['GO_2023'], **kwargs)


def write_inputs(tmp_path):
    p = proteins()
    differential = tmp_path/'proteins.csv'
    p.to_csv(differential,index=False)
    snapshot = tmp_path/'sets.json'
    snapshot.write_text(json.dumps(SETS))
    a_dir,_=measured.run_measured_ora('example',gene_sets_path=snapshot,differential_path=differential,
        output_root=tmp_path/'measured',adjusted_p_threshold=1.)
    b.run_pathway_centric('example',differential_results=pd.read_csv(differential),gene_sets=SETS,
        output_dir=tmp_path/'b',min_size=1,permutations=99,q_cutoff=1.)
    return differential,a_dir,tmp_path/'b/example',snapshot


def test_measured_contingency_uses_exact_experimental_universe():
    all_a, audit = measured.measured_enrichment(proteins(),SETS,'primary')
    term = all_a.query("direction == 'increased'").iloc[0]
    assert term.universe_size == 7  # A..G; duplicate A once, X absent; UNANNOTATED excluded
    assert term.input_candidate_count == 3
    assert term.measured_pathway_size == 2
    assert term.overlap_count == 2
    assert not audit.set_index('analysis_gene').loc['UNANNOTATED','in_universe']
    expected = fisher_exact([[2,1],[0,4]],alternative='greater').pvalue
    assert term.enrichment_p_value == pytest.approx(expected)
    np.testing.assert_allclose(all_a.enrichment_adjusted_p_value,
                              multipletests(all_a.enrichment_p_value,method='fdr_bh')[1])
    assert len(all_a) == 8  # four measured terms x both directions, including zero hits


def test_universe_is_library_specific_and_unmeasured_annotations_do_not_affect_p():
    all_a,_=measured.measured_enrichment(proteins(),SETS,'primary')
    alternative = json.loads(json.dumps(SETS))
    alternative['GO_2023']['both (GO:0000001)'] += [f'UNMEASURED{i}' for i in range(50)]
    alternative['other']={'one': ['A','D','UNANNOTATED']}
    changed,_=measured.measured_enrichment(proteins(),alternative,'primary')
    pd.testing.assert_series_equal(all_a.enrichment_p_value,
        changed.loc[changed.pathway_library.eq('GO_2023'),'enrichment_p_value'].reset_index(drop=True))
    assert set(changed.loc[changed.pathway_library.eq('other'),'universe_size']) == {3}


def test_categories_counts_and_protein_traceability():
    table,evidence,summary=build()
    assert table.set_index('pathway_identifier').category.to_dict() == {
        'GO:0000001':'A_AND_B','GO:0000002':'A_ONLY','GO:0000003':'B_ONLY'}
    row=table.set_index('pathway_identifier').loc['GO:0000001']
    assert row.total_measured_protein_count == 3
    assert row.total_measured_gene_count == 2
    assert row.level_1_protein_count == 1 and row.level_1_gene_count == 1
    assert row.non_candidate_protein_count == 1 and row.non_candidate_gene_count == 0
    assert evidence.query("pathway_identifier == 'GO:0000001'").Accession.tolist() == ['P0','P1','P7']
    assert set(evidence.query("pathway_identifier == 'GO:0000001' and route_a_contributing").row_position) == {0,1}
    assert summary['categories'] == {'A_AND_B':1,'A_ONLY':1,'B_ONLY':1}
    assert summary['shared_libraries'] == ['GO_2023']


def test_matching_exact_versions_and_ids_and_deterministic_shuffle():
    assert comparison.pathway_key('GO_2023','New label (GO:0000001)') == ('GO_2023','GO:0000001')
    assert comparison.pathway_key('GO_2025','New label (GO:0000001)') != ('GO_2023','GO:0000001')
    assert comparison.pathway_key('custom','Mixed Case') != comparison.pathway_key('custom','mixed case')
    first,evidence,_=build()
    shuffled=[frame.sample(frac=1,random_state=3) for frame in saved_tables()]
    second,other,_=build(shuffled)
    pd.testing.assert_frame_equal(first,second)
    pd.testing.assert_frame_equal(evidence.sort_values(['pathway_identifier','row_position']).reset_index(drop=True),
                                  other.sort_values(['pathway_identifier','row_position']).reset_index(drop=True))


def test_unshared_libraries_are_flagged_not_forced_matches():
    tables=list(saved_tables())
    for i in [4,5,6]:
        tables[i]=tables[i].assign(pathway_library='GO_2025')
    table,_,summary=comparison.build_comparison(*tables, {'GO_2025':SETS['GO_2023']},['GO_2023'])
    assert not table.directly_comparable.any()
    assert summary['categories']['A_AND_B'] == 0
    assert summary['shared_libraries'] == []


def test_different_local_definitions_are_flagged():
    changed=json.loads(json.dumps(SETS))
    changed['GO_2023']['both (GO:0000001)']=['A','C']
    table,detail,_=build(a_definitions=changed)
    row=table.set_index('pathway_identifier').loc['GO:0000001']
    assert row.definition_agreement == 'different'
    assert not row.directly_comparable
    members=detail.query("pathway_identifier == 'GO:0000001'").set_index('analysis_gene')
    assert members.loc['B','in_b_definition']
    assert not members.loc['B','in_a_definition']


def test_both_a_directions_retained_without_false_direction_agreement():
    tables=list(saved_tables())
    all_a=tables[1].copy()
    both=all_a.pathway.str.contains('both')
    all_a.loc[both,'enrichment_adjusted_p_value']=.001
    all_a.loc[both & all_a.direction.eq('decreased'),'candidate_genes']='D;E'
    tables[1]=all_a
    tables[2]=measured.consolidate_measured(all_a)
    table,_,_=build(tables)
    row=table.set_index('pathway_identifier').loc['GO:0000001']
    assert row.route_a_direction == 'decreased;increased'
    assert row.direction_agreement == 'not_meaningful'
    assert len(json.loads(row.route_a_cluster_records)) == 2


def test_levels_do_not_modify_b_statistics_or_selection(tmp_path):
    p=proteins()
    first,_=b.run_pathway_centric('x',differential_results=p,gene_sets=SETS,output_dir=tmp_path/'first',
                                min_size=1,permutations=99)
    changed=p.assign(p_value=.8,q_value=.9,log2FC=.1)
    second,_=b.run_pathway_centric('x',differential_results=changed,gene_sets=SETS,output_dir=tmp_path/'second',
                                 min_size=1,permutations=99)
    for name in ['all_pathways','pathway_candidates','pathway_clusters','representative_pathways']:
        pd.testing.assert_frame_equal(first[name],second[name])


def test_frozen_pre_update_route_b_results():
    """Frozen numerical/cluster oracle generated from the pre-update B source."""
    p=proteins()
    fixture=json.loads((Path(__file__).parent/'fixtures/pathway_b_regression.json').read_text())
    stats,_=b.prepare_proteins(p)
    all_b,members,_=b.test_pathways(stats,SETS,min_size=1,permutations=99)
    candidates=b.select_candidates(all_b,q_cutoff=1.)
    clusters,reps=b.consolidate_candidates(candidates,members)
    for name,frame in [('all_pathways',all_b),('candidates',candidates),('clusters',clusters),('representatives',reps)]:
        assert json.loads(frame.to_json(orient='records',double_precision=15)) == fixture[name]


def test_saved_comparison_regenerates_offline_without_touching_sources(tmp_path,monkeypatch):
    differential,a_dir,b_dir,snapshot=write_inputs(tmp_path)
    protected=list(a_dir.rglob('*'))+list(b_dir.rglob('*'))+[differential,snapshot]
    before=provenance([path for path in protected if path.is_file()])
    monkeypatch.setattr(b,'test_pathways',lambda *args,**kwargs: pytest.fail('must not retest B'))
    monkeypatch.setattr(measured,'measured_enrichment',lambda *args,**kwargs: pytest.fail('must not retest A'))
    import urllib.request
    monkeypatch.setattr(urllib.request,'urlopen',lambda *args,**kwargs: pytest.fail('no network allowed'))
    kwargs=dict(a_mode='measured',a_dir=a_dir,b_dir=b_dir,differential_path=differential)
    one,summary=comparison.run_comparison('example',output_root=tmp_path/'review1',**kwargs)
    two,_=comparison.run_comparison('example',output_root=tmp_path/'review2',**kwargs)
    for filename in ['pathway_comparison.csv','pathway_protein_evidence.csv','comparison_summary.json','comparison_metadata.json']:
        assert (one/filename).read_bytes() == (two/filename).read_bytes()
    assert before == provenance([Path(path) for path in before])
    metadata=json.loads((one/'comparison_metadata.json').read_text())
    assert metadata['common_input_verified_for_A']


def test_output_guards_and_hash_mismatch(tmp_path):
    differential,a_dir,b_dir,snapshot=write_inputs(tmp_path)
    with pytest.raises(ValueError,match='overlap'):
        comparison.run_comparison('example',a_mode='measured',a_dir=a_dir,b_dir=b_dir,
            differential_path=differential,output_root=b_dir)
    p=pd.read_csv(differential)
    p.loc[0,'t_statistic']=7.
    p.to_csv(differential,index=False)
    with pytest.raises(ValueError,match='hash'):
        comparison.run_comparison('example',a_mode='measured',a_dir=a_dir,b_dir=b_dir,
            differential_path=differential,output_root=tmp_path/'review')
    with pytest.raises(ValueError,match='overlap'):
        guard_output(a_dir,[a_dir])


def test_empty_candidates_and_bad_definitions(tmp_path):
    tables=list(saved_tables())
    tables[2]=tables[2].iloc[:0]
    tables[5]=tables[5].iloc[:0]
    tables[6]=tables[6].iloc[:0]
    table,evidence,summary=build(tables)
    assert table.empty and evidence.empty
    assert 'route_b_q_value' in table and 'row_position' in evidence
    assert summary['categories'] == {'A_AND_B':0,'A_ONLY':0,'B_ONLY':0}
    path=tmp_path/'bad.json'
    path.write_text('{"GO":{"term":"ABC"}}')
    with pytest.raises(ValueError):
        load_definitions(path)


def test_ambiguous_identity_rejected():
    with pytest.raises(ValueError,match='Ambiguous'):
        comparison.pathway_key('GO','GO:0000001 and GO:0000002')
    tables=list(saved_tables())
    tables[4]=pd.concat([tables[4],tables[4].iloc[:1]],ignore_index=True)
    with pytest.raises(ValueError,match='Duplicate'):
        build(tables)


def test_stale_legacy_members_are_rejected():
    p,all_a,ac,am,*_=saved_tables()
    comparison.validate_a_members(p,am)
    changed=am.copy()
    changed.loc[changed.index[0],'t_statistic']=100.
    with pytest.raises(ValueError,match='statistics'):
        comparison.validate_a_members(p,changed)


def test_measured_output_cannot_replace_legacy_or_b(tmp_path):
    from src.pathway_review_io import ROOT
    differential=tmp_path/'p.csv'
    proteins().to_csv(differential,index=False)
    snapshot=tmp_path/'sets.json'
    snapshot.write_text(json.dumps(SETS))
    for root in [ROOT/'results/pathway_analysis', ROOT/'results/pathway_centric']:
        with pytest.raises(ValueError,match='overlap'):
            measured.run_measured_ora('example',gene_sets_path=snapshot,differential_path=differential,output_root=root)
