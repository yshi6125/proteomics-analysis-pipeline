# Ordinary versus limma moderated t: Route B validation

**Classification: GENERALLY ROBUST WITH BORDERLINE DIFFERENCES.** All 28 canonical candidates are retained, all six requested themes retain every canonical anchor, and all candidate-union directions agree. Moderation adds four near-threshold metabolic/ATP-synthesis terms; it reveals no material problem requiring a production change. Keep canonical ordinary-t Route B at 100,000 permutations. This conclusion concerns this specific standard limma sensitivity experiment, not every possible moderation model.

## Design, software and missing values

R 4.2.1; limma 3.54.0. `lmFit(method="ls")`, followed by `eBayes(trend=FALSE, robust=FALSE, proportion=0.01, stdev.coef.lim=c(0.1,4))`. No weights, blocking, peptide-count trend, alternate contrast or tuning to improve agreement. The established standard empirical-Bayes method moderates residual variances, leaving fitted group coefficients unchanged. [Limma reference manual](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf). The installed-version session record and RDS fitted object are saved.

Input: the same approved `data/processed/data_log2.csv` before batch subtraction, using only its 20 sample-abundance columns (not embedded historical statistic columns), and `sample_metadata.csv`. Design columns: intercept, group__DM, batch__B2; full rank 3; reference NDM and batch B1. Contrast is the adjusted DM-minus-NDM coefficient. There are 5 samples per group/batch cell.

The existing gate (at least 4 observed abundances per group, an estimable full-rank observed design and positive residual degrees of freedom) is preserved. No new imputation or filtering: the actual matrix is complete. All 1,965 proteins were successfully tested; zero excluded; residual df=17 for all. Same Route B mapping: 1,947 eligible rows / 1,946 genes; 18 rows lack valid gene symbols. Moderation uses all successfully tested proteins before mapping. HLA-A remains the sole duplicate and is still averaged arithmetically.

OLS recomputation and limma unmoderated coefficients/t values reproduce production within rtol=1e-10, atol=1e-12. Estimated prior df=3.1570374, prior variance=0.38849992. Posterior variance, moderated t, p and BH protein q formulas were checked independently.

Route B uses the same local gene_sets.json, measured sizes 5–500, signed mean t, equal-size uniform random sets without replacement, size-specific SeedSequence([42,size]), 100,000 permutations, inclusive plus-one two-sided p, global BH over 1,617 pathways, strict q<0.05 and p<=1. Neither Jaccard nor Louvain was run or changed.

## Protein and gene statistics

Large statistic change means |moderated t − ordinary t| ≥1. Large rank change means a shift of at least 10% of the compared population in the average-tie rank of |t|. Signed and absolute-t rank correlations are both reported. These are descriptive thresholds fixed before inspecting results, not selection criteria.

| Population | N | Pearson t | Spearman signed t | Spearman absolute t | Sign changes | Large t changes | Large rank changes |
|---|---:|---:|---:|---:|---:|---:|---:|
| proteins | 1965 | 0.994838 | 0.998328 | 0.993549 | 0 | 4 | 44 |
| eligible proteins | 1947 | 0.994813 | 0.998308 | 0.993500 | 0 | 4 | 44 |
| genes | 1946 | 0.994810 | 0.998308 | 0.993499 | 0 | 4 | 44 |

Only four proteins have |Δt|≥1: HSPD1, ATP5A1, PRDX1 and COMMD4. All preserve sign. HSPD1 has the largest change (−5.296573 → −3.095780), reflecting an increased residual variance estimate after moderation. Per-protein coefficients, ordinary/moderated t, variances, df, rank and sign changes are saved in protein_statistic_comparison.csv; gene statistics are saved separately.

## Exact pathway stability

Ordinary candidates=28; moderated=32; intersection=28; ordinary-only=0; moderated-only=4; candidate-set Jaccard=0.875. No canonical candidate is lost. Candidate identities/statuses are recorded in candidate_union.csv and the complete side-by-side pathway_comparison.csv.

| Measure | Pearson | Spearman / rank |
|---|---:|---:|
| Pathway mean_t | 0.995461 | 0.994962 |
| Pathway p | 0.979186 | 0.981001 |
| Pathway q | 0.989948 | 0.980562 |

Candidate-union p-rank correlation=0.952293; q-rank=0.952528. Absolute observed-mean directions agree for 1565/1,617 pathways (96.78%) and all 32/32 candidate-union pathways. The 52 mean-sign changes affect no selected pathway; gene/protein signs do not change. Changes in weighting magnitude can shift a mixed-sign pathway mean across zero.

## Theme/family stability

These fixed canonical GO anchors group related terms for descriptive review; themes overlap and are not separately tested or given combined p-values. Mitochondrial membrane/matrix support is additionally retained in the complete catalog.

| Theme | Ordinary anchors selected | Moderated anchors selected | Direction agreements |
|---|---:|---:|---:|
| Mitochondrial respiration / OXPHOS | 6/6 | 6/6 | 6/6 |
| Electron transport | 7/7 | 7/7 | 7/7 |
| Fatty-acid oxidation / metabolism | 5/5 | 5/5 | 5/5 |
| Coagulation / fibrinolysis | 4/4 | 4/4 | 4/4 |
| Platelet granules | 2/2 | 2/2 | 2/2 |
| Extracellular matrix | 1/1 | 1/1 | 1/1 |

All six themes remain supported. Two added proton-motive ATP-synthesis terms describe the same measured gene set, so they are not independent new biological programs. Added dicarboxylic-acid metabolism and acetyl-CoA biosynthesis terms extend related metabolic support. Statistical selection and negative/positive observed means do not establish pathway activation/inhibition or causation.

## Candidate disagreements

All four additions have ordinary q=0.052552 and moderated q between 0.044606 and 0.048510. They fall within the prespecified q=0.025–0.075 boundary window under both methods. No member gene changes by |Δt|≥1.

| Added term | Ordinary q | Moderated q | Members | Largest absolute Δt | Top 3 share of absolute member changes | Top changed genes |
|---|---:|---:|---:|---:|---:|---|
| Dicarboxylic Acid Metabolic Process (GO:0043648) | 0.052552 | 0.048510 | 20 | 0.375957 | 51.1% | MDH2;IDH1;ASS1 |
| Proton Motive Force-Driven ATP Synthesis (GO:0015986) | 0.052552 | 0.045901 | 18 | 0.805839 | 70.4% | NDUFA10;SDHA;NDUFS8 |
| Proton Motive Force-Driven Mitochondrial ATP Synthesis (GO:0042776) | 0.052552 | 0.045901 | 18 | 0.805839 | 70.4% | NDUFA10;SDHA;NDUFS8 |
| acetyl-CoA Biosynthetic Process (GO:0006085) | 0.052552 | 0.044606 | 7 | 0.400403 | 88.2% | DLD;ACSS2;ACLY |

The two ATP-synthesis terms share 18 genes. NDUFA10 accounts for 55.5% of absolute member changes (ordinary −4.776388 → moderated −3.970548); the top three account for 70.4%, with other genes shifting in the opposite direction. Acetyl-CoA biosynthesis has 7 genes; DLD accounts for 64.6%, and DLD/ACSS2/ACLY for 88.2%. Dicarboxylic-acid metabolism has 20 genes; its largest change accounts for 22.8%, and MDH2/IDH1/ASS1 together 51.1%, indicating a more distributed change. Thus the first two families have concentrated modest changes; the dicarboxylic-acid term has broader modest changes. None reflects a large/sign-changing protein statistic.

These concentration measures describe changes to the observed mean, not causal attribution of crossing q=0.05. All four pathway means actually move toward zero while p-values improve. The competitive null also changes: eligible-gene mean shifts from −0.119608 to −0.107959 and sd from 1.395673 to 1.371834. Paired size-specific 100k nulls for sizes 7, 18 and 20 reproduce both methods’ empirical p-values and are saved in changing_pathway_nulls.npz, changing_size_null_summary.csv and changing_pathway_null_diagnostics.csv. Changes in null tails and global BH ranking contribute to selection, so attributing crossing solely to one protein would be unjustified. No leave-one-out inference or method adjustment was performed.

Individual protein rows and gene-level means for every changing pathway are in disagreement_protein_evidence.csv and disagreement_gene_evidence.csv. Mean-change concentration and descriptive summaries are in candidate_disagreements.csv.

## Production protection

All protected production/input files match the pre-experiment SHA-256 manifest. This covers canonical differential analysis, Route A, canonical 100k Route B (all 28 candidates and Jaccard communities/representatives), canonical A/B comparison, processed inputs and existing validation records. No production file was changed. Hashes are saved in protected_files_before.json and sensitivity_metadata.json. The first incomplete run was retained in a temporary directory after a validation join guard was corrected to compare unique gene memberships; the completed deterministic rerun is the delivered experiment.
## Tests and final verification

`MPLCONFIGDIR=/private/tmp/proteomics-mpl python -m pytest -q`: **106 passed in 13.91s; zero failures**, including 10 new moderation-validation tests and all 96 existing tests. Checks cover independent variance/t/p/BH formulas, OLS/design/contrast equivalence, mapping and duplicate means, unchanged pathway family/memberships, paired-null empirical tails, strict cutoff/ranks, theme anchors, change concentration, overwrite/no-clustering guards and protected hashes. The repeated moderated inference is exactly identical to the initial computed pathway results. Final hashes verify all 68 protected files unchanged and saved implementation hashes accurate. Machine-readable test results are in test_results.json.
