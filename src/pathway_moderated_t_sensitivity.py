"""Isolated limma ordinary/moderated-t validation; never clusters or writes production."""
from __future__ import annotations
import json
import subprocess
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from src import differential_abundance as da, pathway_centric as b
from src.pathway_review_io import ROOT, provenance, dataframe_hash, sha256_file

KEYS=['pathway_library','pathway']
# Descriptive thresholds, chosen before seeing this experiment; neither selects pathways.
LARGE_T_CHANGE=1.0
LARGE_RANK_FRACTION=.10
BORDERLINE_WINDOW=(.025,.075)


def compare_statistics(frame, ordinary='ordinary_t', moderated='moderated_t'):
    frame=frame.copy()
    a=frame[ordinary];z=frame[moderated]
    frame['t_change']=z-a
    frame['absolute_t_change']=frame.t_change.abs()
    frame['ordinary_absolute_t_rank']=a.abs().rank(ascending=False,method='average')
    frame['moderated_absolute_t_rank']=z.abs().rank(ascending=False,method='average')
    frame['absolute_t_rank_change']=frame.moderated_absolute_t_rank-frame.ordinary_absolute_t_rank
    frame['sign_changed']=np.sign(a)!=np.sign(z)
    frame['large_statistic_change']=frame.absolute_t_change.ge(LARGE_T_CHANGE)
    frame['large_rank_change']=frame.absolute_t_rank_change.abs().ge(LARGE_RANK_FRACTION*len(frame))
    summary={'n':len(frame),'pearson_t':float(pearsonr(a,z).statistic),
        'spearman_signed_t':float(spearmanr(a,z).statistic),
        'spearman_absolute_t':float(spearmanr(a.abs(),z.abs()).statistic),
        'sign_changes':int(frame.sign_changed.sum()),'large_statistic_changes':int(frame.large_statistic_change.sum()),
        'large_rank_changes':int(frame.large_rank_change.sum()),'median_absolute_t_change':float(frame.absolute_t_change.median()),
        'max_absolute_t_change':float(frame.absolute_t_change.max())}
    return frame,summary


def compare_pathways(ordinary,moderated):
    fields=['pathway_size','mean_t','p_value','q_value','observed_direction']+b.DIRECTION_COLUMNS[1:]
    a=ordinary[KEYS+fields].rename(columns={x:'ordinary_'+x for x in fields})
    z=moderated[KEYS+fields].rename(columns={x:'moderated_'+x for x in fields})
    wide=a.merge(z,on=KEYS,validate='one_to_one')
    if len(wide)!=len(a) or len(wide)!=len(z):raise ValueError('Testing families differ')
    assert wide.ordinary_pathway_size.equals(wide.moderated_pathway_size)
    wide=wide.rename(columns={'ordinary_p_value':'ordinary_p','moderated_p_value':'moderated_p','ordinary_q_value':'ordinary_q','moderated_q_value':'moderated_q'})
    for method in ['ordinary','moderated']:
        wide[method+'_candidate']=wide[method+'_q'].lt(.05)
        wide[method+'_p_rank']=wide[method+'_p'].rank(method='average')
        wide[method+'_q_rank']=wide[method+'_q'].rank(method='average')
    wide['direction_agreement']=wide.ordinary_observed_direction.eq(wide.moderated_observed_direction)
    wide['candidate_status_changed']=wide.ordinary_candidate.ne(wide.moderated_candidate)
    wide['transition']=np.select([wide.ordinary_candidate & wide.moderated_candidate,wide.ordinary_candidate,wide.moderated_candidate],['shared','ordinary_only','moderated_only'],default='neither')
    summary={name:int(wide.transition.eq(name).sum()) for name in ['shared','ordinary_only','moderated_only']}
    summary.update(tested_pathways=len(wide),ordinary_candidates=int(wide.ordinary_candidate.sum()),moderated_candidates=int(wide.moderated_candidate.sum()),
        candidate_jaccard=summary['shared']/(summary['shared']+summary['ordinary_only']+summary['moderated_only']),
        direction_agreement_all=int(wide.direction_agreement.sum()),direction_agreement_candidate_union=int(wide.loc[wide.transition.ne('neither'),'direction_agreement'].sum()))
    for field in ['mean_t','p','q']:
        summary['pearson_'+field]=float(pearsonr(wide['ordinary_'+field],wide['moderated_'+field]).statistic)
        summary['spearman_'+field]=float(spearmanr(wide['ordinary_'+field],wide['moderated_'+field]).statistic)
    union=wide.transition.ne('neither')
    for field in ['p','q']:
        summary['spearman_'+field+'_candidate_union']=float(spearmanr(wide.loc[union,'ordinary_'+field],wide.loc[union,'moderated_'+field]).statistic)
    return wide,summary


def explain_disagreements(wide,members,genes,proteins):
    gene_values=genes.set_index('analysis_gene')
    details=[]; summaries=[]
    for _,row in wide.loc[wide.candidate_status_changed].iterrows():
        selected=members.loc[members.pathway_library.eq(row.pathway_library)&members.pathway.eq(row.pathway),'analysis_gene']
        frame=gene_values.loc[selected].reset_index()
        frame['pathway_library']=row.pathway_library;frame['pathway']=row.pathway
        n=len(frame); total=frame.t_change.sum();absolute=frame.t_change.abs().sum()
        frame['contribution_to_mean_t_change']=frame.t_change/n
        frame['ordinary_contribution_to_mean_t']=frame.ordinary_t/n
        frame['moderated_contribution_to_mean_t']=frame.moderated_t/n
        order=frame.t_change.abs().sort_values(ascending=False).index
        top=frame.loc[order].head(3)
        # Descriptive concentration of changes, not re-testing after removal.
        summaries.append({**{k:row[k] for k in KEYS},'transition':row.transition,'ordinary_q':row.ordinary_q,'moderated_q':row.moderated_q,
            'pathway_size':n,'mean_t_change':row.moderated_mean_t-row.ordinary_mean_t,
            'ordinary_near_cutoff':BORDERLINE_WINDOW[0]<=row.ordinary_q<=BORDERLINE_WINDOW[1],
            'moderated_near_cutoff':BORDERLINE_WINDOW[0]<=row.moderated_q<=BORDERLINE_WINDOW[1],
            'members_with_large_t_change':int(frame.large_statistic_change.sum()),
            'fraction_members_with_large_t_change':float(frame.large_statistic_change.mean()),
            'median_absolute_member_t_change':float(frame.absolute_t_change.median()),
            'max_absolute_member_t_change':float(frame.absolute_t_change.max()),
            'top1_share_absolute_member_changes':float(frame.absolute_t_change.max()/absolute) if absolute else 0.,
            'top3_share_absolute_member_changes':float(top.absolute_t_change.sum()/absolute) if absolute else 0.,
            'top3_signed_mean_change_contribution':float(top.t_change.sum()/n),
            'mean_change_excluding_top3_descriptive':float((total-top.t_change.sum())/(n-len(top))) if n>len(top) else np.nan,
            'top3_changed_genes':';'.join(top.analysis_gene)})
        details.append(frame)
    summary=pd.DataFrame(summaries)
    detail=pd.concat(details,ignore_index=True) if details else pd.DataFrame()
    if not detail.empty:
        protein_detail=detail[KEYS+['analysis_gene']].merge(proteins,on='analysis_gene',validate='many_to_many')
    else:protein_detail=pd.DataFrame()
    return summary,detail,protein_detail


def run_sensitivity():
    folder=ROOT/'results/pathway_centric/DM_vs_NDM';target=folder/'moderated_t_sensitivity'
    if target.exists() and any(target.iterdir()):raise ValueError('Refusing to overwrite an existing validation experiment')
    dirs=[ROOT/'results/pathway_protein',ROOT/'results/pathway_comparison',ROOT/'results/differential_abundance',ROOT/'data/processed']
    protected=[p for d in dirs for p in d.rglob('*') if p.is_file()]+[p for p in folder.rglob('*') if p.is_file() and not p.is_relative_to(target)]
    before=provenance(protected)
    meta=json.loads((folder/'analysis_metadata.json').read_text())
    assert meta['permutations']==100000 and meta['pathway_q_cutoff']==.05 and meta['pathway_p_cutoff']==1.
    ordinary=pd.read_csv(ROOT/'results/differential_abundance/differential_abundance_DM_vs_NDM.csv')
    assert dataframe_hash(ordinary)==meta['differential_input_sha256']
    assert sha256_file(folder/'gene_sets.json')==meta['gene_sets_sha256']
    data,samples=da.load_final_analysis_data()
    ref=json.loads((ROOT/'results/differential_abundance/analysis_metadata_DM_vs_NDM.json').read_text())
    assert ref['reference_group']=='NDM' and ref['comparison_group']=='DM' and ref['model']=='abundance ~ group + batch'
    reproduced,excluded=da.run_differential_abundance(data,samples,reference_group='NDM',comparison_group='DM',min_observed_per_group=4)
    a=ordinary.sort_values('row_position').reset_index(drop=True);r=reproduced.sort_values('row_position').reset_index(drop=True)
    assert a.row_position.equals(r.row_position)
    for field in ['log2FC','t_statistic','p_value']:
        np.testing.assert_allclose(a[field],r[field],rtol=1e-10,atol=1e-12)
    ids=da.validate_final_matrix(data,samples,'NDM','DM')
    sm=samples.set_index('sample_id').loc[ids];sm=sm.loc[sm.group.isin(['NDM','DM'])]
    design=pd.concat([pd.DataFrame({'intercept':1.,'group__DM':sm.group.eq('DM').astype(float)},index=sm.index),pd.get_dummies(sm.batch.astype(str),prefix='batch',prefix_sep='__',drop_first=True,dtype=float)],axis=1)
    target.mkdir(parents=True,exist_ok=True)
    (target/'protected_files_before.json').write_text(json.dumps(before,indent=2)+'\n')
    # Use exactly successfully tested proteins; no new imputation and no changed eligibility.
    abundance=data.iloc[a.row_position.astype(int)][list(sm.index)].apply(pd.to_numeric,errors='coerce')
    assert np.isfinite(abundance.to_numpy()).all() # actual approved matrix is complete
    abundance.insert(0,'row_position',a.row_position.to_numpy())
    abundance.to_csv(target/'abundance_matrix.csv',index=False)
    design.rename_axis('sample_id').reset_index().to_csv(target/'design_matrix.csv',index=False)
    sm.reset_index().to_csv(target/'analysis_samples.csv',index=False)
    excluded.to_csv(target/'excluded_proteins.csv',index=False)
    result=subprocess.run(['Rscript',str(ROOT/'scripts/moderated_t_limma.R'),str(target)],capture_output=True,text=True,check=True)
    (target/'limma_execution.log').write_text(result.stdout+result.stderr)
    fit=pd.read_csv(target/'limma_statistics.csv')
    joined=a.merge(fit,on='row_position',validate='one_to_one')
    np.testing.assert_allclose(joined.log2FC,joined.coefficient,rtol=1e-10,atol=1e-12)
    np.testing.assert_allclose(joined.t_statistic,joined.ordinary_t_limma,rtol=1e-10,atol=1e-12)
    assert len(joined)==len(a) and np.isfinite(joined.moderated_t).all()
    annotated=joined.rename(columns={'t_statistic':'ordinary_t'})
    annotated['analysis_gene']=annotated['Gene.name'].astype('string').str.strip().str.upper()
    protein_table,protein_summary=compare_statistics(annotated)
    protein_table.to_csv(target/'protein_statistic_comparison.csv',index=False)
    mod=a.copy();mod['t_statistic']=joined.moderated_t.to_numpy();mod['p_value']=joined.moderated_p.to_numpy();mod['q_value']=joined.moderated_q.to_numpy()
    mod.to_csv(target/'moderated_differential_statistics.csv',index=False)
    ordinary_genes,ordinary_audit=b.prepare_proteins(a);moderated_genes,moderated_audit=b.prepare_proteins(mod)
    pd.testing.assert_series_equal(ordinary_audit.eligible,moderated_audit.eligible)
    assert ordinary_genes.index.equals(moderated_genes.index)
    gene_table,gene_summary=compare_statistics(pd.DataFrame({'analysis_gene':ordinary_genes.index,'ordinary_t':ordinary_genes.values,'moderated_t':moderated_genes.values}))
    gene_table.to_csv(target/'gene_statistic_comparison.csv',index=False)
    print('OLS reproduced; limma ordinary statistics matched. Running 100k moderated Route B inference only.',flush=True)
    start=time.monotonic()
    results,members,audit=b.test_pathways(moderated_genes,json.loads((folder/'gene_sets.json').read_text()),min_size=meta['minimum_pathway_size'],max_size=meta['maximum_pathway_size'],permutations=meta['permutations'],seed=meta['permutation_seed'])
    runtime=time.monotonic()-start
    results.to_csv(target/'moderated_all_pathways.csv',index=False)
    b.select_candidates(results).to_csv(target/'moderated_pathway_candidates.csv',index=False)
    members.to_csv(target/'moderated_pathway_members.csv',index=False)
    audit.to_csv(target/'gene_set_audit.csv',index=False)
    saved_members=pd.read_csv(folder/'pathway_members.csv')
    pd.testing.assert_frame_equal(members[KEYS+['analysis_gene']].sort_values(KEYS+['analysis_gene']).reset_index(drop=True),saved_members[KEYS+['analysis_gene']].drop_duplicates().sort_values(KEYS+['analysis_gene']).reset_index(drop=True))
    wide,pathway_summary=compare_pathways(pd.read_csv(folder/'all_pathways.csv'),results)
    wide.to_csv(target/'pathway_comparison.csv',index=False)
    wide.loc[wide.transition.ne('neither')].to_csv(target/'candidate_union.csv',index=False)
    changes,member_details,protein_details=explain_disagreements(wide,members,gene_table,protein_table)
    changes.to_csv(target/'candidate_disagreements.csv',index=False)
    member_details.to_csv(target/'disagreement_gene_evidence.csv',index=False)
    protein_details.to_csv(target/'disagreement_protein_evidence.csv',index=False)
    assert before==provenance(protected),'Protected production/input files changed'
    summary={'proteins':protein_summary,'genes':gene_summary,'pathways':pathway_summary}
    (target/'comparison_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    metadata={'experiment':'ordinary versus standard limma moderated t; validation only','software':'R 4.2.1; limma 3.54.0',
        'model':'abundance ~ group + batch','design_columns':list(design.columns),'design_rank':int(np.linalg.matrix_rank(design)),
        'contrast':'group__DM coefficient; DM minus NDM adjusted for batch B2 versus B1','samples':len(sm),
        'missing_value_handling':'Same approved normalized/log2 matrix, before batch subtraction; no new imputation. Same >=4 observed values per group, full-rank observed-row design and positive residual degrees of freedom gate. Actual selected matrix has no missing/nonfinite entries; lmFit uses all 20 observations per protein.',
        'proteins_successfully_tested':len(a),'proteins_excluded':len(excluded),'eligible_protein_rows':int(moderated_audit.eligible.sum()),'eligible_genes':len(moderated_genes),
        'moderation':'lmFit(method=ls); eBayes(trend=FALSE, robust=FALSE, proportion=0.01, stdev.coef.lim=c(0.1,4)); no weights, blocks or alternate contrasts. Moderation across all successfully tested proteins, including rows without gene symbols.',
        'canonical_settings':meta,'descriptive_thresholds':{'large_absolute_t_change':LARGE_T_CHANGE,'large_absolute_t_rank_change_fraction':LARGE_RANK_FRACTION,'borderline_q_window':BORDERLINE_WINDOW},
        'ols_coefficient_and_t_reproduction_passed':True,'eligible_gene_and_pathway_memberships_identical':True,'clustering_run':False,
        'pathway_runtime_seconds':runtime,'protected_source_hashes':before,'protected_files_unchanged':True,
        'implementation_hashes':provenance([Path(__file__),ROOT/'scripts/moderated_t_limma.R',ROOT/'src/pathway_centric.py',ROOT/'src/differential_abundance.py']),
        'limma_documentation':'https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf'}
    (target/'sensitivity_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)
    return target

if __name__=='__main__':run_sensitivity()
