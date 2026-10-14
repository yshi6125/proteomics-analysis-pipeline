import ast
import hashlib
import json
import os
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
import pytest
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from src import protein_evidence as evidence
from src import pathway_ora as ora
from src import pathway_protein as a
from src import pathway_centric as b
from src import pathway_comparison as c
from src import pathway_routes as routes
from src.pathway_review_io import ROOT, provenance


def proteins():
    return pd.DataFrame(dict(row_position=range(9), Accession=[f'P{i}' for i in range(9)],
        **{'Gene.name':['A','B','C','D','E','F','G','A','UNANNOTATED']},
        log2FC=[2.,.8,.2,-2.,-.7,-.2,0.,.1,3.],t_statistic=[4.,3.,2.,-4.,-3.,-2.,.1,.2,5.],
        p_value=[.001,.002,.003,.01,.02,.03,.9,.8,.001],q_value=[.01,.02,.03,.1,.2,.3,.99,.9,.01]))


SETS={'GO_2023':{'both (GO:0000001)':['A','B','X'],'a_only (GO:0000002)':['B','C'],
                 'b_only (GO:0000003)':['D','E','F'],'neutral (GO:0000004)':['G']}}


def saved_tables():
    p=proteins()
    all_a,_=ora.run_local_ora(p,SETS)
    all_a['enrichment_adjusted_p_value']=.5
    all_a['enrichment_p_value']=.2
    select=all_a.pathway.str.contains('both|a_only') & all_a.direction.eq('increased')
    all_a.loc[select,'enrichment_adjusted_p_value']=.01
    all_a.loc[select,'enrichment_p_value']=.001
    ac=ora.consolidate_protein_pathways(all_a)
    am=evidence.build_pathway_candidate_members(all_a,evidence.assign_candidate_origin(p))
    stats,_=b.prepare_proteins(p)
    all_b,mem,_=b.test_pathways(stats,SETS,min_size=1,permutations=99)
    all_b['q_value']=.5
    all_b.loc[all_b.pathway.str.contains('both|b_only'),'q_value']=.01
    candidates=b.select_candidates(all_b)
    bc,_=b.consolidate_candidates(candidates,mem)
    return p,all_a,ac,am,all_b,candidates,bc


def test_experimental_universe_and_fisher_bh():
    results,audit=ora.run_local_ora(proteins(),SETS)
    term=results.query("direction == 'increased'").iloc[0]
    assert term.universe_size==7 and term.input_candidate_count==3
    assert not audit.set_index('analysis_gene').loc['UNANNOTATED','in_universe']
    assert term.enrichment_p_value==pytest.approx(fisher_exact([[2,1],[0,4]],alternative='greater').pvalue)
    np.testing.assert_allclose(results.enrichment_adjusted_p_value,multipletests(results.enrichment_p_value,method='fdr_bh')[1])
    assert len(results)==8


def test_comparison_categories_traceability_and_duplicate_counts():
    table,detail,summary=c.build_comparison(*saved_tables(),SETS,['GO_2023'])
    assert table.set_index('pathway_identifier').category.to_dict()=={'GO:0000001':'A_AND_B','GO:0000002':'A_ONLY','GO:0000003':'B_ONLY'}
    row=table.set_index('pathway_identifier').loc['GO:0000001']
    assert row.total_measured_protein_count==3 and row.total_measured_gene_count==2
    assert row.non_candidate_protein_count==1 and row.non_candidate_gene_count==0
    assert detail.query("pathway_identifier == 'GO:0000001'").Accession.tolist()==['P0','P1','P7']
    assert summary['categories']=={'A_AND_B':1,'A_ONLY':1,'B_ONLY':1}


def test_exact_matching_and_shuffle_determinism():
    assert c.pathway_key('GO_2023','Renamed (GO:0000001)')==('GO_2023','GO:0000001')
    assert c.pathway_key('GO_2025','Renamed (GO:0000001)')!=('GO_2023','GO:0000001')
    one,detail,_=c.build_comparison(*saved_tables(),SETS,['GO_2023'])
    shuffled=[frame.sample(frac=1,random_state=3) for frame in saved_tables()]
    two,other,_=c.build_comparison(*shuffled,SETS,['GO_2023'])
    pd.testing.assert_frame_equal(one,two)
    pd.testing.assert_frame_equal(detail,other)


def test_versions_and_definition_differences_not_forced_matches():
    tables=list(saved_tables())
    for i in [4,5,6]:tables[i]=tables[i].assign(pathway_library='GO_2025')
    table,_,summary=c.build_comparison(*tables,{'GO_2025':SETS['GO_2023']},['GO_2023'])
    assert not table.directly_comparable.any() and summary['categories']['A_AND_B']==0
    changed=json.loads(json.dumps(SETS));changed['GO_2023']['both (GO:0000001)']=['A','C']
    table,_,_=c.build_comparison(*saved_tables(),SETS,['GO_2023'],a_definitions=changed)
    assert not table.set_index('pathway_identifier').loc['GO:0000001'].directly_comparable


def test_frozen_b_regression_and_levels_do_not_affect_selection():
    fixture=json.loads((ROOT/'tests/fixtures/pathway_b_regression.json').read_text())
    stats,_=b.prepare_proteins(proteins())
    results,mem,_=b.test_pathways(stats,SETS,min_size=1,permutations=99)
    candidates=results.copy()
    clusters,reps=b.consolidate_candidates(candidates,mem)
    for name,frame in [('all_pathways',results),('candidates',candidates),('clusters',clusters),('representatives',reps)]:
        assert json.loads(frame[list(fixture[name][0])].to_json(orient='records',double_precision=15))==fixture[name]
    changed=proteins().assign(p_value=.9,q_value=.99,log2FC=.1)
    stats2,_=b.prepare_proteins(changed)
    results2,_,_=b.test_pathways(stats2,SETS,min_size=1,permutations=99)
    pd.testing.assert_frame_equal(results,results2)


def test_invalid_b_statistics_are_audited():
    p=proteins();p.loc[0,'t_statistic']=np.inf;p.loc[1,'Gene.name']=None
    stats,audit=b.prepare_proteins(p)
    assert set(audit.exclusion_reason)=={'','missing_gene_symbol','invalid_t_statistic'}
    assert not audit.iloc[:2].eligible.any()


def test_representative_ties_and_seed():
    stats,_=b.prepare_proteins(proteins())
    result,mem,_=b.test_pathways(stats,{'GO':{'z':['A','B'],'a':['A','B']}},min_size=2,permutations=99)
    clusters,reps=b.consolidate_candidates(result,mem)
    assert list(reps.pathway)==['a']
    other,_=b.consolidate_candidates(result.iloc[::-1],mem.iloc[::-1])
    pd.testing.assert_frame_equal(clusters.reset_index(drop=True),other.reset_index(drop=True))


def test_empty_comparison_has_stable_schema():
    tables=list(saved_tables())
    for i in [2,5,6]:tables[i]=tables[i].iloc[:0]
    table,detail,summary=c.build_comparison(*tables,SETS,['GO_2023'])
    assert table.empty and detail.empty and 'route_b_q_value' in table and 'row_position' in detail
    assert summary['categories']=={'A_AND_B':0,'A_ONLY':0,'B_ONLY':0}


def test_stale_members_rejected():
    p,_,_,am,*_=saved_tables();c.validate_a_members(p,am)
    am.loc[am.index[0],'t_statistic']=100.
    with pytest.raises(ValueError,match='statistics'):c.validate_a_members(p,am)


def test_active_b_numeric_function_ast_unchanged():
    previous=ast.parse((ROOT/'archive/route_a_legacy/src/pathway_centric.py').read_text())
    current=ast.parse((ROOT/'src/pathway_centric.py').read_text())
    old={n.name:ast.dump(n,include_attributes=False) for n in previous.body if isinstance(n,ast.FunctionDef)}
    new={n.name:ast.dump(n,include_attributes=False) for n in current.body if isinstance(n,ast.FunctionDef)}
    for name in ['prepare_proteins','load_gene_sets','consolidate_candidates']:
        assert old[name]==new[name]
    # New fields/defaults may change packaging; all inference arithmetic assignments stay identical.
    def arithmetic(tree):
        fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='test_pathways')
        names={'universe','values','rng','member_stats','mean','null','lower','upper','p','direction'}
        return [ast.dump(n,include_attributes=False) for n in ast.walk(fn) if isinstance(n,ast.Assign)
                and any(isinstance(t,ast.Name) and t.id in names for t in n.targets)]
    assert arithmetic(previous)==arithmetic(current)


def test_archive_manifest_preserves_every_byte():
    records=json.loads((ROOT/'archive/route_a_legacy/ARCHIVE_MANIFEST.json').read_text())
    for item in records:
        assert hashlib.sha256((ROOT/item['archived_path']).read_bytes()).hexdigest()==item['sha256']


def test_canonical_outputs_and_b_not_modified():
    records=json.loads((ROOT/'archive/route_a_legacy/ARCHIVE_MANIFEST.json').read_text())
    baseline=[r for r in records if '/baseline/results/' in r['archived_path']]
    assert baseline
    for r in baseline:
        assert hashlib.sha256((ROOT/r['archived_path']).read_bytes()).hexdigest()==r['sha256']
        if r['original_path'].startswith('results/pathway_protein/'):
            assert hashlib.sha256((ROOT/r['original_path']).read_bytes()).hexdigest()==r['sha256']


def test_production_generation_without_archive_or_obsolete_results(tmp_path):
    """An isolated production-only installation proves there is no archive dependency."""
    import shutil
    isolated=tmp_path/'project';(isolated/'src').mkdir(parents=True)
    names=['__init__.py','pathway_review_io.py','protein_evidence.py','pathway_redundancy.py',
           'pathway_ora.py','pathway_protein.py','pathway_centric.py','pathway_comparison.py','pathway_routes.py','pathway_direction_diagnostics.py']
    for name in names:shutil.copy2(ROOT/'src'/name,isolated/'src'/name)
    proteins().to_csv(isolated/'p.csv',index=False)
    (isolated/'sets.json').write_text(json.dumps(SETS))
    code="""
import pandas as pd,json
from pathlib import Path
from src.pathway_protein import run_protein_centric
from src.pathway_centric import run_pathway_centric
from src.pathway_comparison import run_comparison
run_protein_centric('x',gene_sets_path='sets.json',differential_path='p.csv',adjusted_p_threshold=1.)
run_pathway_centric('x',differential_results=pd.read_csv('p.csv'),gene_sets=json.loads(Path('sets.json').read_text()),min_size=1,permutations=99,q_cutoff=1.)
run_comparison('x',differential_path='p.csv')
"""
    subprocess.run(['python','-c',code],cwd=isolated,env={**os.environ,'MPLCONFIGDIR':'/private/tmp/proteomics-mpl'},check=True,capture_output=True,text=True)
    for name in ['pathway_analysis','pathway_consolidation','pathway_analysis_measured']:
        assert not (isolated/'results'/name).exists()
    assert not (isolated/'archive').exists()
    assert (isolated/'results/pathway_comparison/x/pathway_comparison.csv').exists()
    original=tmp_path/'reference_installation'
    shutil.copytree(ROOT/'archive/route_a_legacy/src',original/'src')
    shutil.copy2(isolated/'p.csv',original/'p.csv')
    shutil.copy2(isolated/'sets.json',original/'sets.json')
    script="from src.pathway_protein import run_protein_centric; run_protein_centric('x',gene_sets_path='sets.json',differential_path='p.csv',adjusted_p_threshold=1.)"
    subprocess.run(['python','-c',script],cwd=original,env={**os.environ,'MPLCONFIGDIR':'/private/tmp/proteomics-mpl'},check=True,capture_output=True,text=True)
    for filename in ['all_pathways.csv','pathway_candidates.csv','pathway_clusters.csv','representative_pathways.csv','candidate_members.csv','universe_genes.csv','protein_levels.csv']:
        previous=pd.read_csv(original/'results/pathway_protein/x'/filename)
        current=pd.read_csv(isolated/'results/pathway_protein/x'/filename)
        pd.testing.assert_frame_equal(previous,current,check_exact=True)



def test_comparison_regeneration_preserves_inputs_and_review(tmp_path,monkeypatch):
    from shutil import copytree
    adir=tmp_path/'a';bdir=tmp_path/'b'
    copytree(ROOT/'results/pathway_protein/DM_vs_NDM',adir)
    copytree(ROOT/'results/pathway_centric/DM_vs_NDM',bdir)
    before=provenance([p for d in [adir,bdir] for p in d.rglob('*') if p.is_file()])
    monkeypatch.setattr(b,'test_pathways',lambda *args,**kwargs:pytest.fail('comparison must not test pathways'))
    target,_=c.run_comparison(a_dir=adir,b_dir=bdir,output_root=tmp_path/'review')
    current=pd.read_csv(target/'pathway_comparison.csv');old=pd.read_csv(ROOT/'results/pathway_comparison/DM_vs_NDM/pathway_comparison.csv')
    pd.testing.assert_frame_equal(current,old,check_exact=True)
    selection=target/'pathway_review_selection.csv';selection.write_text(selection.read_text()+'\n')
    selected=selection.read_bytes()
    c.run_comparison(a_dir=adir,b_dir=bdir,output_root=tmp_path/'review')
    assert selection.read_bytes()==selected and before==provenance([Path(p) for p in before])
    with pytest.raises(ValueError,match='overlap'):c.run_comparison(a_dir=adir,b_dir=bdir,output_root=adir)


def test_cli_and_dispatch_only_production_options(monkeypatch,tmp_path):
    import inspect
    assert not {'a_mode','branch'} & set(inspect.signature(c.run_comparison).parameters)
    env={**os.environ,'MPLCONFIGDIR':'/private/tmp/proteomics-mpl'}
    for name in ['pathway_protein','pathway_routes','pathway_comparison']:
        text=subprocess.check_output(['python','-m','src.'+name,'--help'],env=env,text=True)
        assert '--a-mode' not in text and '--a-branch' not in text
    snapshot=tmp_path/'sets.json';snapshot.write_text(json.dumps(SETS));calls=[]
    monkeypatch.setattr(routes,'run_protein_centric',lambda **kw:calls.append('A'))
    monkeypatch.setattr(routes,'run_pathway_centric',lambda **kw:calls.append('B'))
    routes.run_routes('both',gene_sets_path=snapshot)
    assert calls==['A','B']
