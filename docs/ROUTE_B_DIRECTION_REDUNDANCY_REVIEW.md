# Route B direction and redundancy review

October 14, 2026. Dataset: DM_vs_NDM. Frozen discovery: 1,946 eligible genes,
1,617 tested pathways, 28 candidates. Route A and all canonical A/B result files
are unchanged. No new threshold has been selected.

## A. Direction diagnostics

The enriched diagnostic views expose:

| Field | Definition |
|---|---|
| `observed_direction` | Sign of observed pathway mean t: positive, negative, zero |
| `direction` | Existing field retained; zero was named non_directional |
| `null_tail` | Smaller empirical tail: higher_than_random or lower_than_random; exact ties balanced_or_tied |
| `empirical_lower_tail_probability` | `(1 + count(null_mean <= observed_mean))/(permutations+1)` |
| `empirical_upper_tail_probability` | `(1 + count(null_mean >= observed_mean))/(permutations+1)` |
| `n_positive_members`, `n_negative_members`, `n_zero_members` | Counts of unique eligible gene means >0, <0, =0 |
| `positive_fraction`, `negative_fraction` | Corresponding counts / measured pathway size |
| `dominant_direction_fraction` | max(positive count, negative count) / measured pathway size |

Duplicate-gene means are identical to the frozen calculation; protein rows are
not counted multiple times. These fields never enter testing, correction,
candidate selection, graph construction or representative ranking. A smaller tail
label alone does not mean statistical significance. Neither observed direction
nor tail labels imply activation/inhibition.

The old tables did not save their lower/upper probabilities. The diagnostic command
replays precisely the recorded size-seeded random-set nulls (10,000 draws, seed 42),
then verifies the resulting two-sided probabilities agree with **every** stored
p-value. It does not replace p/q, change sample count/procedure or reselect pathways.
No differential analysis is refit. All frozen source hashes are checked afterward.

Among tested terms, 183 have **negative observed mean but higher-than-random tail**;
659 have negative mean/lower tail and 775 positive mean/higher tail. This shows why
absolute observed sign and competitive-null direction must be separate.
All 28 candidates have consistent signs/tails in this dataset (20 negative/lower,
8 positive/higher). Candidate dominant-direction fractions range from 0.6328125
to 1.0. No coherence cutoff is applied.

## B. Current canonical graph

At Jaccard 0.25, seeded Louvain 42:

- 28 nodes and 46 edges.
- 10 connected components and 11 communities.
- Five singleton communities; largest community has eight pathways.
- Community-size distribution: five of size 1, two of size 2, one of size 3,
  two of size 4, one of size 8.

The graph is weighted by measured unique-gene Jaccard. All nodes are retained;
edges require the same existing direction and thresholded overlap. Cross-direction
raw overlaps are retained in pairwise diagnostics but cannot become graph edges.

`current_graph_diagnostics.csv` reports every candidate's size, cluster,
representative flag, degree, maximum raw overlap, strongest overlap partner,
mean connected-neighbor overlap and maximum same-direction overlap. A zero-degree
node's neighbor mean is blank, not falsely zero. Strongest-overlap ties use stable
identities, not a biological/AI judgement.

## C. Threshold sensitivity

Only redundancy was rerun; every row uses the same candidates and memberships.

| Jaccard | Edges | Components | Clusters | Singleton clusters | Largest cluster |
|---:|---:|---:|---:|---:|---:|
| .15 | 73 | 6 | 9 | 3 | 9 |
| .20 | 56 | 9 | 10 | 4 | 9 |
| .25 | 46 | 10 | 11 | 5 | 8 |
| .30 | 42 | 10 | 11 | 5 | 8 |
| .35 | 38 | 10 | 12 | 5 | 5 |
| .40 | 35 | 10 | 12 | 5 | 5 |

Seed is always 42. Summary CSV includes distributions, representatives and identical
candidate-evidence/membership hashes for each threshold. Long/wide assignment
CSV files show all pathways. Component and representative files remain separate
from canonical outputs. Numeric Cxxx labels can change merely because communities
are numbered again; use memberships and adjusted Rand index (ARI), not label
text, to assess stability.

## D. Membership-based cluster stability

The .25 and .30 partitions are **identical**, ARI=1, as are .35 and .40.
Relative to .25, ARI is .906302 at .20, .874525 at .15, and .778947 at .35/.40.
These are descriptive comparisons over the same 28 pathways, not significance tests
or proof of a universal optimal threshold. Only one requested seed is evaluated;
no stochastic seed sweep or algorithm replacement was performed.

Using original canonical cluster identities:

- C002 (the two small electron-transport terms) stays together and separate at
  every tested threshold.
- C003 (energy derivation plus three Complex III terms) stays together and separate
  at every threshold.
- The C004 fatty-acid oxidation trio stays together throughout; C005 joins it only
  at .15.
- C006 mitochondrial envelope joins C001 at .15/.20 and is a singleton at .25–.40.
- C001's five respiration/BP terms stay together. Its three broad membrane/CC
  terms split into a separate community at .35/.40.
- Coagulation/fibrinolysis, platelet/granule, protein secretion, extracellular
  matrix and mitochondrial matrix retain their own respective communities
  throughout. Some receive different numeric cluster labels.

`canonical_cluster_transitions.csv` identifies unchanged, merged, split and
split/recombined memberships relative to .25 without being confused by numbering.

## E. Specific measured-overlap findings

### Respiration, OXPHOS and electron transport

C001 and C003 have real graph connections. Aerobic electron transport versus
energy derivation shares **14/30 union genes**, J=.466667; ATP-coupled transport
versus energy derivation is .451613. The two six-gene Complex III cellular-component annotations
have identical measured membership (J=1); the six-gene ubiquinol-to-cytochrome-c
process shares five genes with each (5/7 union, J=.714286). Their tight internal overlap and the
energy-term bridges are consistent with Louvain retaining a separate community,
even though all twelve C001/C003 pathways are connected.

By contrast, C002's maximum overlap with C001 is **6/36=.166667**, for electron
transport chain versus cellular respiration. C002 versus C003 has **zero** shared
measured genes. Lowering to .15 creates a connection to the broader respiration
component but does not merge C002's Louvain community. These are not simply
several edges narrowly missing .25.

### Fatty-acid metabolism

The broad fatty-acid metabolism term differs substantially from the oxidation
trio: overlap with catabolic process is **6/39=.153846**, with oxidation **5/40=.125**,
and with beta-oxidation **4/43=.093023**. Only .15 connects/merges the broad term
through its strongest catabolic overlap. Biological similarity of names does not
establish sufficient measured membership redundancy.

### Coagulation, fibrinolysis, platelet and secretion

The strongest coagulation/platelet bridge is negative regulation of coagulation
versus platelet alpha-granule lumen, **6/38=.157895**. It creates a shared component
at .15, but Louvain still preserves separate communities. Coagulation versus
protein secretion peaks at **1/22=.045455**; secretion versus platelet lumen is
**1/40=.025**. Shared clinical/biological discussion does not imply statistical
redundancy among these measured sets.

The diagnostic group file uses the user-specified original cluster sets solely
as report subsets. Names, groups, embeddings, literature and GO hierarchy do not
enter any graph or community calculation.

## F. Components versus Louvain

At .25, nine of ten components each correspond to a single community. The only
split is the twelve-node C001/C003 component into eight- and four-node communities.
C002, C004/C005, platelet/coagulation/secretion and broad mitochondrial singletons
are mostly separated by **missing eligible Jaccard edges**. C001 versus C003 is
specifically a **Louvain partition within an already connected component**.

At .15, components also combine the small electron-transport pair with the main
respiration component, and platelet with coagulation; Louvain still separates
them. Lower edge threshold alone therefore does not implement the desired broad
biological theme organization—and that is not the job of statistical redundancy.

## G. Recommendation: KEEP 0.25

Retain the canonical threshold and Louvain method/seed. The .25–.30 stable plateau
is evidence for a modestly robust membership partition. Lower thresholds merge
some low-overlap sets but do not resolve the observed respiration subdivisions;
higher thresholds split the membrane community without reducing component count.
Fewer communities is not itself an improvement. This dataset provides no compelling
reason to alter the selected threshold or representative rule.

Broader mitochondrial/metabolic themes may be organized in a future interpretation
layer, separately from statistical redundancy. No such grouping or AI clustering
is implemented here.

## H. Outputs and reproducibility

Run offline from frozen local inputs:

```bash
python -m src.pathway_redundancy_sensitivity --comparison DM_vs_NDM
```

Output directory:
`results/pathway_centric/DM_vs_NDM/redundancy_sensitivity/`.
It contains the six requested files, plus:

- `pathway_cluster_assignments_wide.csv`
- `representative_pathways_by_threshold.csv`
- `cluster_partition_stability.csv`
- `canonical_cluster_transitions.csv`
- `diagnostic_group_pairwise_jaccard.csv`
- `all_pathways_with_direction.csv`
- `pathway_candidates_with_direction.csv`
- `representative_pathways_with_direction.csv`
- `pathway_comparison_with_direction.csv`

The enriched views join new fields to unchanged saved rows. Comparison fields are
prefixed `route_b_`; A fields/categories remain unchanged. This sidecar design
protects the frozen canonical files, particularly `pathway_clusters.csv`, pending
review. Running this command repeats diagnostic/null replay and consolidation only;
it never overwrites canonical results. File hashes, seeds, policies, software and
implementation hashes are recorded in `redundancy_sensitivity_metadata.json`.

## I. Validation

Tests cover deterministic direction/tails/coherence (including zero and tied tails),
observed-sign versus null-tail disagreement, frozen p recovery, descriptive fields
having no cluster/representative effect, every threshold's identical discovery
inputs, fixed-seed reproducibility, actual 28-candidate/11-cluster reconstruction,
sidecar propagation and all frozen source hashes. Existing A/B regression tests
remain active. The complete suite passed: **86 tests**, including eight new diagnostic tests. No core statistical source module was modified for this task.
