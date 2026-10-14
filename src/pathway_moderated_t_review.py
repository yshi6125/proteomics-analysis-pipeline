"""Read completed validation tables and write descriptive review into their directory."""
from __future__ import annotations
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
from src.pathway_moderated_t_sensitivity import KEYS, compare_statistics
from src.pathway_review_io import ROOT, provenance

THEME_ANCHORS={
 'Mitochondrial respiration / OXPHOS':['GO:0019646','GO:0009060','GO:0045333','GO:0042775','GO:0006119','GO:0015980'],
 'Electron transport':['GO:0019646','GO:0022900','GO:0022904','GO:0042775','GO:0006122','GO:0005750','GO:0045275'],
 'Fatty-acid oxidation / metabolism':['GO:0006635','GO:0033539','GO:0009062','GO:0006631','GO:0019395'],
 'Coagulation / fibrinolysis':['GO:0051917','GO:0030195','GO:0051918','GO:0030194'],
 'Platelet granules':['GO:0031091','GO:0031093'],
 'Extracellular matrix':['GO:0062023'],
}


def summarize_themes(wide):
    """Fixed canonical anchors; overlapping families, no theme-level statistical test."""
    rows=[];details=[]
    for theme,anchors in THEME_ANCHORS.items():
        selected=wide.loc[wide.pathway.str.extract(r'(GO:\d+)')[0].isin(anchors)].copy()
        assert len(selected)==len(anchors),f'Missing canonical theme anchor for {theme}'
        selected['theme']=theme;details.append(selected)
        rows.append({'theme':theme,'anchor_terms':len(selected),'ordinary_candidates':int(selected.ordinary_candidate.sum()),
            'moderated_candidates':int(selected.moderated_candidate.sum()),'retained_candidates':int((selected.ordinary_candidate & selected.moderated_candidate).sum()),
            'direction_agreements':int(selected.direction_agreement.sum()),'max_ordinary_q':float(selected.ordinary_q.max()),'max_moderated_q':float(selected.moderated_q.max())})
    return pd.DataFrame(rows),pd.concat(details,ignore_index=True)


def changing_null_diagnostics(wide,genes,permutations=100000,seed=42):
    """Reproduce only changing-size nulls; paired identical sets; no inference changes."""
    selected=wide.loc[wide.candidate_status_changed]
    genes=genes.sort_values('analysis_gene')
    ordinary=genes.ordinary_t.to_numpy();moderated=genes.moderated_t.to_numpy()
    nulls={};summary=[];terms=[]
    for size in sorted(selected.ordinary_pathway_size.unique()):
        size=int(size);rng=np.random.default_rng(np.random.SeedSequence([seed,size]))
        pair=np.empty((2,permutations))
        for j in range(permutations):
            indices=rng.choice(len(genes),size,replace=False)
            pair[:,j]=[ordinary[indices].mean(),moderated[indices].mean()]
        for i,method in enumerate(['ordinary','moderated']):
            nulls[f'{method}_size_{size}']=pair[i]
            summary.append({'method':method,'pathway_size':size,'null_mean':float(pair[i].mean()),'null_sd':float(pair[i].std(ddof=1)),'null_q001':float(np.quantile(pair[i],.001)),'null_q01':float(np.quantile(pair[i],.01))})
        for _,row in selected.loc[selected.ordinary_pathway_size.eq(size)].iterrows():
            for i,method in enumerate(['ordinary','moderated']):
                observed=row[method+'_mean_t'];lower=(1+np.count_nonzero(pair[i]<=observed))/(permutations+1);upper=(1+np.count_nonzero(pair[i]>=observed))/(permutations+1)
                p=min(1.,2*min(lower,upper))
                np.testing.assert_allclose(p,row[method+'_p'],rtol=1e-12,atol=1e-15)
                terms.append({**{k:row[k] for k in KEYS},'method':method,'mean_t':observed,'null_mean':float(pair[i].mean()),'null_sd':float(pair[i].std(ddof=1)),
                    'descriptive_standardized_departure':float((observed-pair[i].mean())/pair[i].std(ddof=1)),
                    'lower_tail':lower,'upper_tail':upper,'empirical_p':p})
    return pd.DataFrame(summary),pd.DataFrame(terms),nulls


def render_review(target):
    # Exact validation directory only, never production output roots.
    if target.resolve()!= (ROOT/'results/pathway_centric/DM_vs_NDM/moderated_t_sensitivity').resolve():raise ValueError('Unexpected sensitivity output directory')
    meta=json.loads((target/'sensitivity_metadata.json').read_text());before=meta['protected_source_hashes']
    assert before==provenance(before)
    wide=pd.read_csv(target/'pathway_comparison.csv');genes=pd.read_csv(target/'gene_statistic_comparison.csv');proteins=pd.read_csv(target/'protein_statistic_comparison.csv')
    themes,anchors=summarize_themes(wide);themes.to_csv(target/'theme_stability.csv',index=False);anchors.to_csv(target/'theme_anchor_pathways.csv',index=False)
    null_summary,null_terms,nulls=changing_null_diagnostics(wide,genes,meta['canonical_settings']['permutations'],meta['canonical_settings']['permutation_seed'])
    null_summary.to_csv(target/'changing_size_null_summary.csv',index=False);null_terms.to_csv(target/'changing_pathway_null_diagnostics.csv',index=False)
    np.savez_compressed(target/'changing_pathway_nulls.npz',**nulls)
    eligible=proteins.loc[proteins.analysis_gene.notna() & proteins.analysis_gene.ne('')]
    eligible_comparison,eligible_summary=compare_statistics(eligible)
    eligible_comparison.to_csv(target/'eligible_protein_statistic_comparison.csv',index=False)
    summary=json.loads((target/'comparison_summary.json').read_text());summary['eligible_proteins']=eligible_summary
    summary['gene_background']={method:{'mean':float(genes[method+'_t'].mean()),'sd':float(genes[method+'_t'].std(ddof=1))} for method in ['ordinary','moderated']}
    flips=wide.loc[~wide.direction_agreement]
    summary['nonsignificant_pathway_direction_changes']={'count':len(flips),'candidate_in_either_method':int((flips.ordinary_candidate|flips.moderated_candidate).sum()),'maximum_absolute_mean_t':float(flips[['ordinary_mean_t','moderated_mean_t']].abs().to_numpy().max()),'minimum_q_either_method':float(flips[['ordinary_q','moderated_q']].to_numpy().min())}
    summary['classification']='GENERALLY ROBUST WITH BORDERLINE DIFFERENCES'
    (target/'comparison_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    changes=pd.read_csv(target/'candidate_disagreements.csv')
    metrics=[]
    for unit in ['proteins','eligible_proteins','genes']:
        s=summary[unit];metrics.append(f"| {unit.replace('_',' ')} | {s['n']} | {s['pearson_t']:.6f} | {s['spearman_signed_t']:.6f} | {s['spearman_absolute_t']:.6f} | {s['sign_changes']} | {s['large_statistic_changes']} | {s['large_rank_changes']} |")
    p=summary['pathways'];hyper=pd.read_csv(target/'limma_hyperparameters.csv').iloc[0]
    lines=['# Ordinary versus limma moderated t: Route B validation','','**Classification: GENERALLY ROBUST WITH BORDERLINE DIFFERENCES.** All 28 canonical candidates are retained, all six requested themes retain every canonical anchor, and all candidate-union directions agree. Moderation adds four near-threshold metabolic/ATP-synthesis terms; it reveals no material problem requiring a production change. Keep canonical ordinary-t Route B at 100,000 permutations. This conclusion concerns this specific standard limma sensitivity experiment, not every possible moderation model.','','## Design, software and missing values','','R 4.2.1; limma 3.54.0. `lmFit(method="ls")`, followed by `eBayes(trend=FALSE, robust=FALSE, proportion=0.01, stdev.coef.lim=c(0.1,4))`. No weights, blocking, peptide-count trend, alternate contrast or tuning to improve agreement. The established standard empirical-Bayes method moderates residual variances, leaving fitted group coefficients unchanged. [Limma reference manual](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf). The installed-version session record and RDS fitted object are saved.','','Input: the same approved `data/processed/data_log2.csv` before batch subtraction, using only its 20 sample-abundance columns (not embedded historical statistic columns), and `sample_metadata.csv`. Design columns: intercept, group__DM, batch__B2; full rank 3; reference NDM and batch B1. Contrast is the adjusted DM-minus-NDM coefficient. There are 5 samples per group/batch cell.','','The existing gate (at least 4 observed abundances per group, an estimable full-rank observed design and positive residual degrees of freedom) is preserved. No new imputation or filtering: the actual matrix is complete. All 1,965 proteins were successfully tested; zero excluded; residual df=17 for all. Same Route B mapping: 1,947 eligible rows / 1,946 genes; 18 rows lack valid gene symbols. Moderation uses all successfully tested proteins before mapping. HLA-A remains the sole duplicate and is still averaged arithmetically.','','OLS recomputation and limma unmoderated coefficients/t values reproduce production within rtol=1e-10, atol=1e-12. Estimated prior df='+format(hyper.df_prior,'.8g')+', prior variance='+format(hyper.s2_prior,'.8g')+'. Posterior variance, moderated t, p and BH protein q formulas were checked independently.','','Route B uses the same local gene_sets.json, measured sizes 5–500, signed mean t, equal-size uniform random sets without replacement, size-specific SeedSequence([42,size]), 100,000 permutations, inclusive plus-one two-sided p, global BH over 1,617 pathways, strict q<0.05 and p<=1. Neither Jaccard nor Louvain was run or changed.','','## Protein and gene statistics','','Large statistic change means |moderated t − ordinary t| ≥1. Large rank change means a shift of at least 10% of the compared population in the average-tie rank of |t|. Signed and absolute-t rank correlations are both reported. These are descriptive thresholds fixed before inspecting results, not selection criteria.','','| Population | N | Pearson t | Spearman signed t | Spearman absolute t | Sign changes | Large t changes | Large rank changes |','|---|---:|---:|---:|---:|---:|---:|---:|']+metrics+['','Only four proteins have |Δt|≥1: HSPD1, ATP5A1, PRDX1 and COMMD4. All preserve sign. HSPD1 has the largest change (−5.296573 → −3.095780), reflecting an increased residual variance estimate after moderation. Per-protein coefficients, ordinary/moderated t, variances, df, rank and sign changes are saved in protein_statistic_comparison.csv; gene statistics are saved separately.','','## Exact pathway stability','',f"Ordinary candidates={p['ordinary_candidates']}; moderated={p['moderated_candidates']}; intersection={p['shared']}; ordinary-only={p['ordinary_only']}; moderated-only={p['moderated_only']}; candidate-set Jaccard={p['candidate_jaccard']:.3f}. No canonical candidate is lost. Candidate identities/statuses are recorded in candidate_union.csv and the complete side-by-side pathway_comparison.csv.",'','| Measure | Pearson | Spearman / rank |','|---|---:|---:|']
    for field in ['mean_t','p','q']:lines.append(f"| Pathway {field} | {p['pearson_'+field]:.6f} | {p['spearman_'+field]:.6f} |")
    lines += ['',f"Candidate-union p-rank correlation={p['spearman_p_candidate_union']:.6f}; q-rank={p['spearman_q_candidate_union']:.6f}. Absolute observed-mean directions agree for {p['direction_agreement_all']}/1,617 pathways ({p['direction_agreement_all']/1617:.2%}) and all 32/32 candidate-union pathways. The 52 mean-sign changes affect no selected pathway; gene/protein signs do not change. Changes in weighting magnitude can shift a mixed-sign pathway mean across zero.",'','## Theme/family stability','','These fixed canonical GO anchors group related terms for descriptive review; themes overlap and are not separately tested or given combined p-values. Mitochondrial membrane/matrix support is additionally retained in the complete catalog.','','| Theme | Ordinary anchors selected | Moderated anchors selected | Direction agreements |','|---|---:|---:|---:|']
    for _,r in themes.iterrows():lines.append(f"| {r.theme} | {r.ordinary_candidates}/{r.anchor_terms} | {r.moderated_candidates}/{r.anchor_terms} | {r.direction_agreements}/{r.anchor_terms} |")
    lines += ['','All six themes remain supported. Two added proton-motive ATP-synthesis terms describe the same measured gene set, so they are not independent new biological programs. Added dicarboxylic-acid metabolism and acetyl-CoA biosynthesis terms extend related metabolic support. Statistical selection and negative/positive observed means do not establish pathway activation/inhibition or causation.','','## Candidate disagreements','','All four additions have ordinary q=0.052552 and moderated q between 0.044606 and 0.048510. They fall within the prespecified q=0.025–0.075 boundary window under both methods. No member gene changes by |Δt|≥1.','','| Added term | Ordinary q | Moderated q | Members | Largest absolute Δt | Top 3 share of absolute member changes | Top changed genes |','|---|---:|---:|---:|---:|---:|---|']
    for _,r in changes.iterrows():lines.append(f"| {r.pathway} | {r.ordinary_q:.6f} | {r.moderated_q:.6f} | {r.pathway_size} | {r.max_absolute_member_t_change:.6f} | {r.top3_share_absolute_member_changes:.1%} | {r.top3_changed_genes} |")
    lines += ['','The two ATP-synthesis terms share 18 genes. NDUFA10 accounts for 55.5% of absolute member changes (ordinary −4.776388 → moderated −3.970548); the top three account for 70.4%, with other genes shifting in the opposite direction. Acetyl-CoA biosynthesis has 7 genes; DLD accounts for 64.6%, and DLD/ACSS2/ACLY for 88.2%. Dicarboxylic-acid metabolism has 20 genes; its largest change accounts for 22.8%, and MDH2/IDH1/ASS1 together 51.1%, indicating a more distributed change. Thus the first two families have concentrated modest changes; the dicarboxylic-acid term has broader modest changes. None reflects a large/sign-changing protein statistic.','','These concentration measures describe changes to the observed mean, not causal attribution of crossing q=0.05. All four pathway means actually move toward zero while p-values improve. The competitive null also changes: eligible-gene mean shifts from −0.119608 to −0.107959 and sd from 1.395673 to 1.371834. Paired size-specific 100k nulls for sizes 7, 18 and 20 reproduce both methods’ empirical p-values and are saved in changing_pathway_nulls.npz, changing_size_null_summary.csv and changing_pathway_null_diagnostics.csv. Changes in null tails and global BH ranking contribute to selection, so attributing crossing solely to one protein would be unjustified. No leave-one-out inference or method adjustment was performed.','','Individual protein rows and gene-level means for every changing pathway are in disagreement_protein_evidence.csv and disagreement_gene_evidence.csv. Mean-change concentration and descriptive summaries are in candidate_disagreements.csv.','','## Production protection','','All protected production/input files match the pre-experiment SHA-256 manifest. This covers canonical differential analysis, Route A, canonical 100k Route B (all 28 candidates and Jaccard communities/representatives), canonical A/B comparison, processed inputs and existing validation records. No production file was changed. Hashes are saved in protected_files_before.json and sensitivity_metadata.json. The first incomplete run was retained in a temporary directory after a validation join guard was corrected to compare unique gene memberships; the completed deterministic rerun is the delivered experiment.']
    (target/'moderated_t_sensitivity_report.md').write_text('\n'.join(lines)+'\n')
    assert before==provenance(before)
    review_meta={'theme_anchor_definitions':THEME_ANCHORS,'theme_policy':'Fixed canonical GO anchors, overlapping descriptive families; no combined significance','null_diagnostics_policy':'Paired same-seed subsets for changing sizes only, reproduces saved p values; no alternate inference','classification':summary['classification'],'protected_files_unchanged':True,'implementation_hashes':provenance([Path(__file__)])}
    (target/'review_metadata.json').write_text(json.dumps(review_meta,indent=2)+'\n')
    return summary

if __name__=='__main__':
    result=render_review(ROOT/'results/pathway_centric/DM_vs_NDM/moderated_t_sensitivity')
    print(json.dumps(result,indent=2))
