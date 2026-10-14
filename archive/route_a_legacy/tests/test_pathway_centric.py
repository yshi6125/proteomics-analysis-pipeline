import importlib.util
import subprocess

import numpy as np
import pandas as pd
import pytest
from statsmodels.stats.multitest import multipletests

from src import pathway_centric as b
from src import pathway_routes as routes
from src import pathway_analysis as a


def proteins():
    return pd.DataFrame({'Gene.name': ['A', 'B', 'C', 'D', 'E', 'F', None, 'BAD', 'A'],
                         't_statistic': [4., 3., -4., -3., 0., .1, 2., np.inf, 2.],
                         'p_value': [.9] * 9, 'q_value': [.99] * 9})


SETS = {'GO': {'up': ['A', 'B'], 'up2': ['A', 'B', 'F'], 'down': ['C', 'D'], 'tiny': ['E']}}


def test_eligibility_membership_bh_and_reproducibility(tmp_path):
    tables, meta = b.run_pathway_centric(differential_results=proteins(), gene_sets=SETS,
        min_size=2, permutations=499, q_cutoff=1, output_dir=tmp_path)
    assert meta['number_of_eligible_proteins'] == 7
    assert meta['number_of_eligible_genes'] == 6
    assert set(tables['protein_eligibility'].exclusion_reason) == {'', 'missing_gene_symbol', 'invalid_t_statistic'}
    assert tables['pathway_members'].query("pathway == 'up'").shape[0] == 3
    results = tables['all_pathways']
    np.testing.assert_allclose(results.q_value, multipletests(results.p_value, method='fdr_bh')[1])
    assert results.set_index('pathway').loc['up', 'mean_t'] == 3
    tables2, _ = b.run_pathway_centric(differential_results=proteins().assign(p_value=0., q_value=0.),
        gene_sets=SETS, min_size=2, permutations=499, q_cutoff=1, output_dir=tmp_path / 'repeat')
    pd.testing.assert_frame_equal(results, tables2['all_pathways'])
    pd.testing.assert_frame_equal(tables['pathway_clusters'], tables2['pathway_clusters'])
    assert (tmp_path / 'DM_vs_NDM' / 'gene_sets.json').exists()


def test_null_and_empty_candidates(tmp_path):
    stats = pd.Series([1.] * 5, index=list('ABCDE'))
    results, _, _ = b.test_pathways(stats, {'GO': {'constant': ['A', 'B']}}, min_size=2, permutations=50)
    assert results.iloc[0].p_value == 1
    tables, meta = b.run_pathway_centric(differential_results=proteins(), gene_sets={'GO': {'tiny': ['A']}},
                                        output_dir=tmp_path, min_size=2, permutations=20)
    assert meta['number_of_tested_pathways'] == 0
    assert meta['number_of_redundancy_clusters'] == 0
    assert pd.read_csv(tmp_path / 'DM_vs_NDM' / 'pathway_clusters.csv').empty


def test_representative_tie_and_shuffle():
    stats, _ = b.prepare_proteins(proteins())
    results, members, _ = b.test_pathways(stats, {'GO': {'z': ['A', 'B'], 'a': ['A', 'B']}}, min_size=2, permutations=50)
    clusters, representatives = b.consolidate_candidates(results, members)
    assert list(representatives.pathway) == ['a']
    other, _ = b.consolidate_candidates(results.iloc[::-1], members.iloc[::-1])
    pd.testing.assert_frame_equal(clusters.reset_index(drop=True), other.reset_index(drop=True))
    assert len(clusters) == 2


def test_route_dispatch_isolation(monkeypatch, tmp_path):
    calls = []
    snapshot = tmp_path/'sets.json'
    snapshot.write_text('{"GO": {"term": ["A"]}}')
    monkeypatch.setattr(routes, 'run_protein_centric', lambda **kw: calls.append(('a', kw)))
    monkeypatch.setattr(routes, 'run_pathway_centric', lambda **kw: calls.append(('b', kw)))
    routes.run_routes('both', 'example', gene_sets_path=snapshot, q_cutoff=.2)
    assert calls == [('a', {'comparison': 'example', 'gene_sets_path': snapshot}),
                     ('b', {'comparison': 'example', 'q_cutoff': .2, 'gene_sets': {'GO': {'term': ['A']}}})]
    assert b.OUTPUT_DIR != a.RESULTS_DIR
    calls.clear()
    routes.run_routes('pathway')
    assert [call[0] for call in calls] == ['b']


def test_route_a_matches_git_baseline(tmp_path, monkeypatch):
    """Run the pre-change implementation and current Route A on identical fixtures."""
    baseline_path = tmp_path / 'baseline.py'
    baseline_path.write_bytes(subprocess.check_output(['git', 'show', 'HEAD:src/pathway_analysis.py']))
    spec = importlib.util.spec_from_file_location('baseline_pathways', baseline_path)
    baseline = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(baseline)
    frame = pd.DataFrame({'row_position': range(6), 'Accession': list('abcdef'),
                          'Gene.name': list('ABCDEF'), 'log2FC': [2., 1., .2, -1., -.2, .1],
                          't_statistic': [4., 3., 2., -3., -2., .1],
                          'p_value': [.001, .002, .003, .01, .02, .9],
                          'q_value': [.01, .02, .03, .1, .2, .99]})
    def runner(**kw):
        class Response:
            results = pd.DataFrame({'Gene_set': ['GO'], 'Term': ['term'], 'P-value': [.01],
                'Adjusted P-value': [.02], 'Genes': [';'.join(kw['gene_list'])], 'Overlap': ['2/10']})
        return Response()
    outputs = []
    for module, name in [(baseline, 'before'), (a, 'after')]:
        monkeypatch.setattr(module, 'RESULTS_DIR', tmp_path / name)
        monkeypatch.setattr(module, 'FIGURES_DIR', tmp_path / name / 'figures')
        # Avoid rendering figures in a numerical/output regression check.
        monkeypatch.setattr(module, 'create_pathway_summary_plot', lambda *args, **kw: tmp_path / 'plot.png')
        outputs.append(module.run_pathway_analysis(differential_results=frame, pathway_libraries=['GO'], enrichr_runner=runner))
    for before, after in zip(outputs[0][:-1], outputs[1][:-1]):
        pd.testing.assert_frame_equal(before, after)
    assert outputs[0][-1] == outputs[1][-1]
    for before in (tmp_path / 'before').glob('*.csv'):
        assert before.read_bytes() == (tmp_path / 'after' / before.name).read_bytes()
