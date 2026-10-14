import json
import os
from pathlib import Path
import subprocess

import pandas as pd
import pytest

from src import pathway_measured_ora as reference
from src import pathway_protein as production
from src import pathway_centric as b
from src import pathway_comparison as comparison
from src.pathway_review_io import provenance


def inputs(tmp_path):
    proteins=pd.DataFrame({'row_position':range(6),'Accession':list('abcdef'),'Gene.name':list('ABCDEF'),
        'log2FC':[2.,.8,.2,-1.,-.1,0.], 't_statistic':[4.,3.,2.,-3.,-.2,.1],
        'p_value':[.001,.002,.003,.01,.02,.9],'q_value':[.01,.02,.03,.1,.2,.99]})
    differential=tmp_path/'proteins.csv'
    proteins.to_csv(differential,index=False)
    sets={'GO_2023':{'up (GO:0000001)':list('ABC'),'down (GO:0000002)':list('DE'),'all':list('ABCDEF')}}
    snapshot=tmp_path/'sets.json'
    snapshot.write_text(json.dumps(sets))
    return differential,snapshot,sets


def generate(tmp_path):
    differential,snapshot,sets=inputs(tmp_path)
    old,_=reference.run_measured_ora('example',gene_sets_path=snapshot,differential_path=differential,
                                   output_root=tmp_path/'reference',adjusted_p_threshold=1.)
    before=provenance([p for p in old.rglob('*') if p.is_file()])
    new,metadata=production.run_protein_centric('example',gene_sets_path=snapshot,differential_path=differential,
                                              output_root=tmp_path/'production',adjusted_p_threshold=1.)
    return differential,sets,old,new,metadata,before


def test_canonical_a_matches_selected_reference_inference_and_memberships(tmp_path):
    differential,sets,old,new,metadata,before=generate(tmp_path)
    for filename in ['all_pathways.csv','pathway_candidates.csv','pathway_clusters.csv',
                     'representative_pathways.csv','universe_genes.csv','candidate_members.csv']:
        original=pd.read_csv(old/'exploratory'/filename)
        current=pd.read_csv(new/filename)
        columns=[c for c in original.columns if c not in {'analysis_type','analysis_label'}]
        pd.testing.assert_frame_equal(original[columns],current[columns],check_exact=True)
    assert set(pd.read_csv(new/'all_pathways.csv').analysis_type)=={'protein_centric'}
    assert metadata['protein_candidate_levels']==[1,2,3,4,5]
    assert metadata['analysis_route']=='protein_centric'
    assert metadata['enrichment_method']=='local ORA'
    assert metadata['enrichment_background']=='measured proteome / eligible measured genes'
    assert not (new/'primary').exists() and not (new/'exploratory').exists()
    assert before==provenance([Path(p) for p in before])


def test_canonical_comparison_matches_selected_reference_and_keeps_review(tmp_path):
    differential,sets,old,new,metadata,before=generate(tmp_path)
    b_dir=tmp_path/'b/example'
    b.run_pathway_centric('example',gene_sets=sets,differential_results=pd.read_csv(differential),
                         output_dir=tmp_path/'b',min_size=1,permutations=99,q_cutoff=1.)
    historical,hs=comparison.run_comparison('example',a_mode='measured',branch='exploratory',a_dir=old,
        b_dir=b_dir,differential_path=differential,output_root=tmp_path/'reviews')
    protected=provenance([p for directory in [historical,b_dir,new] for p in directory.rglob('*') if p.is_file()])
    canonical,cs=comparison.run_comparison('example',a_dir=new,b_dir=b_dir,differential_path=differential,
                                         output_root=tmp_path/'reviews')
    h=pd.read_csv(historical/'pathway_comparison.csv').drop(columns='route_a_branch')
    c=pd.read_csv(canonical/'pathway_comparison.csv').drop(columns=['route_a_name','route_b_name'])
    pd.testing.assert_frame_equal(h,c,check_exact=True)
    assert cs['categories']==hs['categories']
    assert cs['route_a_name']=='Route A' and cs['route_b_name']=='Route B'
    assert protected==provenance([Path(p) for p in protected])
    review=canonical/'pathway_review_selection.csv'
    text=review.read_text().replace(',,\n', ',both,checked\n',1)
    review.write_text(text)
    comparison.run_comparison('example',a_dir=new,b_dir=b_dir,differential_path=differential,
                              output_root=tmp_path/'reviews')
    assert review.read_text()==text
    meta=json.loads((canonical/'comparison_metadata.json').read_text())
    assert 'route_a_mode' not in meta and 'route_a_branch' not in meta


def test_production_cli_has_no_historical_mode_options():
    env={**os.environ,'MPLCONFIGDIR':'/private/tmp/proteomics-mpl'}
    for module in ['src.pathway_protein','src.pathway_comparison','src.pathway_routes']:
        output=subprocess.check_output(['python','-m',module,'--help'],env=env,text=True)
        assert '--a-mode' not in output and '--a-branch' not in output
        assert 'measured_exploratory' not in output
    output=subprocess.check_output(['python','-m','src.pathway_reference','--help'],env=env,text=True)
    assert '--a-mode' in output and 'Historical/reference' in output


def test_reference_primary_is_not_a_production_scope(tmp_path):
    with pytest.raises(ValueError,match='only the Levels'):
        comparison.run_comparison('example',branch='primary')


def test_missing_local_resource_has_no_silent_network_fallback(monkeypatch,tmp_path):
    monkeypatch.setattr(production,'ROOT',tmp_path)
    with pytest.raises(FileNotFoundError,match='requires local definitions'):
        production.default_snapshot('example')


def test_pre_change_source_manifest_is_intact():
    import hashlib
    root=Path(__file__).resolve().parent.parent
    manifest=json.loads((root/'docs/reference/source/SOURCE_MANIFEST.json').read_text())
    for original,record in manifest.items():
        assert original.startswith('src/')
        assert hashlib.sha256((root/record['snapshot']).read_bytes()).hexdigest()==record['sha256']
