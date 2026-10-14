import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from src import pathway_centric as b
from src import pathway_permutation_sensitivity as sensitivity
from src.pathway_review_io import ROOT, provenance


def example():
    p=pd.DataFrame({'Gene.name':list('ABCDEFGH'),'t_statistic':[4.,3.,2.,1.,-4.,-3.,-2.,0.]})
    stats,_=b.prepare_proteins(p)
    definitions={'GO':{'up':list('ABC'),'down':list('EFG'),'mixed':list('ADE'),'all':list('ABCDEFGH')}}
    return stats,definitions


def test_inference_only_and_memberships_fixed(monkeypatch):
    monkeypatch.setattr(b,'consolidate_candidates',lambda *a,**kw:pytest.fail('clustering forbidden'))
    stats,defs=example();runs={};reference=None
    for count in [100,200,300]:
        results,members,_=b.test_pathways(stats,defs,min_size=1,permutations=count,seed=42)
        if reference is None:reference=members
        else:pd.testing.assert_frame_equal(reference,members,check_exact=True)
        runs[count]=results
    tables=sensitivity.compare_runs(runs)
    assert set(tables['permutation_count_summary'].tested_pathways)=={4}
    assert set(tables['pathway_permutation_comparison'].pathway)==set(defs['GO'])


def test_strict_cutoff_and_entry_exit_labels():
    stats,defs=example()
    result,_,_=b.test_pathways(stats,defs,min_size=1,permutations=100)
    one=result.copy();two=result.copy();three=result.copy()
    one['q_value']=[.049,.05,.2,.3]
    two['q_value']=[.051,.049,.2,.3]
    three['q_value']=[.04,.03,.2,.3]
    tables=sensitivity.compare_runs({10000:one,50000:two,100000:three})
    summary=tables['permutation_count_summary'].set_index('permutations')
    assert summary.loc[10000,'q_lt_005_candidates']==1
    assert summary.loc[10000,'exact_q_boundary_cases']==1
    assert summary.loc[50000,'entering_vs_10k']==1
    assert summary.loc[50000,'leaving_vs_10k']==1
    assert summary.loc[100000,'entering_vs_10k']==1
    assert tables['near_q_cutoff'].shape[0]>=2


def test_changed_testing_family_or_statistic_rejected():
    stats,defs=example()
    result,_,_=b.test_pathways(stats,defs,min_size=1,permutations=100)
    changed=result.copy();changed.loc[0,'mean_t']+=.1
    with pytest.raises(AssertionError):sensitivity.compare_runs({10000:result,50000:changed})


def test_result_alignment_is_deterministic():
    stats,defs=example();runs={}
    for count in [100,200,300]:runs[count]=b.test_pathways(stats,defs,min_size=1,permutations=count)[0]
    first=sensitivity.compare_runs(runs)
    second=sensitivity.compare_runs({n:r.iloc[::-1] for n,r in runs.items()})
    for name in first:pd.testing.assert_frame_equal(first[name],second[name],check_exact=True)


def test_saved_runs_and_all_canonical_hashes_unchanged():
    folder=ROOT/'results/pathway_centric/DM_vs_NDM/permutation_sensitivity'
    info=json.loads((folder/'permutation_sensitivity_metadata.json').read_text())
    from provenance_helpers import assert_recorded_sources
    assert_recorded_sources(info['source_hashes'])
    assert info['all_canonical_source_hashes_unchanged'] and info['redundancy']=='not invoked or evaluated'
    assert info['eligible_genes']==1946 and info['tested_pathways']==1617
    runs={n:pd.read_csv(folder/f'all_pathways_{n}.csv') for n in sensitivity.COUNTS}
    tables=sensitivity.compare_runs(runs)
    for run in runs.values():
        pd.testing.assert_frame_equal(run[sensitivity.CORE_FIXED],runs[10000][sensitivity.CORE_FIXED],check_exact=True)
    reference=pd.read_csv(ROOT/'archive/route_b_10k_canonical_20261014/results/pathway_centric/DM_vs_NDM/all_pathways.csv')
    pd.testing.assert_frame_equal(runs[10000],reference,check_exact=False,rtol=1e-12,atol=1e-14)
