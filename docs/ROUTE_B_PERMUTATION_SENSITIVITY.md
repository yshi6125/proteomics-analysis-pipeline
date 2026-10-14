# Route B permutation-count sensitivity

DM_vs_NDM, October 14, 2026. This analysis changes only the number of permutations
in the existing Route B inference function. No redundancy or clustering function
was invoked or evaluated. All canonical Route A, Route B and comparison files
remain unchanged, verified by their SHA-256 hashes.

## Fixed inputs and procedure

Same local `results/pathway_centric/DM_vs_NDM/gene_sets.json`; same successfully
tested protein CSV; 1,947 eligible protein rows and 1,946 unique gene means;
same duplicate averaging, eligible gene ordering and size limits 5–500;
same mean-t statistic, size-specific `SeedSequence([42, size])`, random subsets
without replacement, inclusive plus-one two-sided tail calculation and global BH
family across BP/CC terms. Every run tests the same 1,617 pathway identities and
measured memberships. Increasing counts extends the same RNG streams.

The requested selection is strictly **q < .05**. Canonical code uses q <= .05;
there are no exact q=.05 cases in any of these three runs, and replayed 10k
candidate identities match the frozen 28 exactly. The canonical p<=1 condition
is automatically met by all empirical p-values.

The 10k replay agrees with frozen mean t/p/q within CSV round-trip precision.
Across the three fresh runs, fixed statistic columns and gene memberships are
checked for exact equality. Core inference source was not changed.

## Candidate stability

| Permutations | Tested | Candidates | Shared with 10k | Entering vs 10k | Leaving vs 10k | Candidate-set Jaccard with 10k |
|---:|---:|---:|---:|---:|---:|---:|
| 10,000 | 1,617 | 28 | 28 | 0 | 0 | 1.000000 |
| 50,000 | 1,617 | 27 | 27 | 0 | 1 | .964286 |
| 100,000 | 1,617 | 28 | 27 | 1 | 1 | .931034 |

All 20 original candidates with q<.025 remain candidates at every count.
Twenty-seven pathways are selected in all three runs. Fifty and 100k share all
27 of the 50k candidates; 100k adds one marginal term. Identical counts at 10k
and 100k therefore do **not** imply identical candidate lists.

### Changing identities

| Pathway | Mean t (fixed) | 10k p | 10k q | 50k p | 50k q | 100k p | 100k q |
|---|---:|---:|---:|---:|---:|---:|---:|
| Protein Secretion (GO:0009306) | 1.333261 | .00079992 | .04619538 | .00115998 | .05879882 | .00115999 | .05683943 |
| Fatty Acid Beta-Oxidation Using acyl-CoA Dehydrogenase (GO:0033539) | -2.066804 | .00099990 | .05389461 | .00095998 | .05543889 | .00077999 | .04504455 |

Protein secretion leaves at both larger counts. The acyl-CoA dehydrogenase term
enters only at 100k. Direction and measured memberships of both remain fixed.
Additional fatty-acid evidence reinforces an already supported annotation family;
losing secretion means that specific marginal pathway should not be described as
robust to permutation count. No biological activation/inhibition or causal
conclusion follows from these changes.

## Rank and p/q changes

Tie-aware Spearman correlation across all 1,617 tested terms:

| Pair | p-rank correlation | q-rank correlation | q-rank correlation within candidate union |
|---|---:|---:|---:|
| 10k / 50k | .999636 | .999489 | .866680 |
| 10k / 100k | .999608 | .999445 | .825083 |
| 50k / 100k | .999969 | .999866 | .893411 |

All-term ranking is very stable. Ranking among discoveries is less stable,
partly because larger runs resolve the many 10k minimum-tail ties; the original
strong-candidate membership remains stable despite reordered significance ranks.
Do not mistake resolution of tied p/q ranks for a change in observed mean t.

| Run versus 10k | Median absolute p change | Maximum absolute p change | Median absolute q change | Maximum absolute q change |
|---|---:|---:|---:|---:|
| 50k | .004038 | .028097 | .003103 | .058486 |
| 100k | .004597 | .030824 | .003733 | .061843 |

These summary changes cover all terms, not just discoveries. Full signed
per-pathway differences are in the side-by-side table. The minimum attainable
two-sided p decreases from approximately .000200 to .000040 to .000020.

## Particular attention near q=.05

At 100k, the following remain just above the cutoff:

- Dicarboxylic Acid Metabolic Process: q=.052552.
- Proton Motive Force-Driven ATP Synthesis: q=.052552.
- Proton Motive Force-Driven Mitochondrial ATP Synthesis: q=.052552.
- acetyl-CoA Biosynthetic Process: q=.052552.

Protein secretion remains borderline but excluded (q=.056839); the newly included
acyl-CoA dehydrogenase term is also marginal (q=.045045). Mitochondrial envelope
remains selected (q=.037131), while electron transport chain moves further below
the cutoff (q=.025872). mRNA splicing falls from q=.121263 to .059924 but remains
excluded. Thus the nearest-cutoff review must include terms outside the original
10k near-boundary band as well.

`near_q_cutoff.csv` includes every term with q in [.025,.075] in any run plus all
changing identities. `nearest_30_to_q_cutoff.csv` provides a second cutoff-focused
view. These diagnostic windows are not candidate-selection rules.

## Recommendation before production freeze

**Use 100,000 permutations for the final production inference run**, subject to
selection of that setting. This analysis has not changed the default or canonical
outputs. Runtime here was approximately 6.8 seconds at 10k, 32.2 seconds at 50k and
64.7 seconds at 100k, so improved resolution is practical for this dataset.

Ten thousand is sufficiently stable for the broad association patterns and the
stronger candidates, but **not sufficiently stable to freeze the exact candidate
catalog or treat q≈.05 calls as robust**. The larger counts do not overturn the
main existing signal, yet they materially change two individual threshold decisions.
Fifty thousand is also not fully identical to 100k, so retaining 10k merely because
both 10k and 100k produce 28 candidates would hide the difference.

These are nested runs from one fixed seed, not independent seed replications or
proof that every 100k borderline call has converged. Even at 100k, q near .05
should be reported as marginal. No extra seed, new test, threshold, correlation
adjustment, moderation or enrichment method was introduced.

## Output files and reproduction

All outputs are under:
`results/pathway_centric/DM_vs_NDM/permutation_sensitivity/`.

- `all_pathways_10000.csv`, `all_pathways_50000.csv`, `all_pathways_100000.csv`:
  full unchanged-method test results.
- `pathway_candidates_<count>.csv`: strict q<.05 candidate identities and evidence.
- `pathway_permutation_comparison.csv`: all 1,617 terms with p/q/ranks and signed
  p/q changes side by side.
- `permutation_count_summary.csv`: counts, overlap, boundary cases and changes.
- `candidate_transitions.csv`: entering/leaving identities relative to 10k.
- `pathway_rank_correlations.csv`: all-term and candidate-union rank comparisons.
- `near_q_cutoff.csv`, `nearest_30_to_q_cutoff.csv`: boundary review tables.
- `eligible_gene_statistics.csv`, `frozen_gene_memberships.csv`: fixed inputs for
  inspecting inference provenance, with hashes recorded in metadata.
- `permutation_sensitivity_metadata.json`: canonical settings, input hashes,
  runtime, strict-cutoff distinction, seed strategy, code/environment provenance
  and confirmation that clustering was not invoked.
- `permutation_sensitivity_report.md`: this report alongside the results.

Offline reproduction:

```bash
python -m src.pathway_permutation_sensitivity --comparison DM_vs_NDM
```

Five focused tests validate strict-cutoff transitions, deterministic matching,
identical testing families/memberships, fixed statistics, 10k replay, canonical
source hashes and absence of consolidation calls. Only those tests were run for
this task; no redundancy tests or clustering evaluation were run.
