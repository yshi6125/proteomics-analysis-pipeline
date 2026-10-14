# Route B 100k production promotion — 2026-10-14

Canonical/default Route B now uses 100,000 permutations. All 1,617 tested terms use the same differential input, 1,947 eligible protein rows / 1,946 eligible genes, arithmetic-mean duplicate aggregation, local gene sets, measured sizes 5–500, size-specific SeedSequence([42, size]) random streams, mean-t statistic, inclusive plus-one two-sided empirical test and global BH. Jaccard remains 0.25; Louvain seed remains 42.

The only effective inferential change is 10k → 100k sampling. The requested strict q<0.05 criterion is now explicit (the former selector used <=); there are no exact q=0.05 cases in either run, so this has no additional selection effect. Diagnostics are descriptive only.

## Results

28 candidates: 27 retained, one added and one removed. Added: Fatty Acid Beta-Oxidation Using acyl-CoA Dehydrogenase (GO:0033539). Removed: Protein Secretion (GO:0009306).

The removed protein-secretion term changes q from 0.0461953804619538 to 0.0568394316056839; the added acyl-CoA dehydrogenase term changes q from 0.0538946105389461 to 0.0450445495545044.

A_AND_B = 5; A_ONLY = 0; B_ONLY = 23. There are 9 Route B communities, including 4 singleton communities (previously 11 communities). New candidate membership links the fatty-acid oxidation and electron-transport communities; secretion loses its singleton. Representatives are selected by the unchanged q, p, library/name ordering.

## Complete candidate catalog and representatives

An asterisk marks the representative. BP = GO_Biological_Process_2023; CC = GO_Cellular_Component_2023.

### C001

- * Aerobic Electron Transport Chain (GO:0019646) [BP]; q=0.00323396766
- Cellular Respiration (GO:0045333) [BP]; q=0.00323396766
- Mitochondrial ATP Synthesis Coupled Electron Transport (GO:0042775) [BP]; q=0.00323396766
- Mitochondrial Inner Membrane (GO:0005743) [CC]; q=0.00323396766
- Mitochondrial Membrane (GO:0031966) [CC]; q=0.00323396766
- Organelle Inner Membrane (GO:0019866) [CC]; q=0.00323396766
- Oxidative Phosphorylation (GO:0006119) [BP]; q=0.00431195688
- Aerobic Respiration (GO:0009060) [BP]; q=0.005389946101

### C002

- * Fatty Acid Beta-Oxidation (GO:0006635) [BP]; q=0.00323396766
- Fatty Acid Catabolic Process (GO:0009062) [BP]; q=0.00431195688
- Fatty Acid Oxidation (GO:0019395) [BP]; q=0.00431195688
- Respiratory Electron Transport Chain (GO:0022904) [BP]; q=0.02560224398
- Electron Transport Chain (GO:0022900) [BP]; q=0.02587174128
- Fatty Acid Beta-Oxidation Using acyl-CoA Dehydrogenase (GO:0033539) [BP]; q=0.04504454955

### C003

- * Energy Derivation By Oxidation Of Organic Compounds (GO:0015980) [BP]; q=0.00431195688
- Mitochondrial Respiratory Chain Complex III (GO:0005750) [CC]; q=0.007699923001
- Respiratory Chain Complex III (GO:0045275) [CC]; q=0.007699923001
- Mitochondrial Electron Transport, Ubiquinol To Cytochrome C (GO:0006122) [BP]; q=0.02736434174

### C004

- * Fatty Acid Metabolic Process (GO:0006631) [BP]; q=0.00323396766

### C005

- * Mitochondrial Envelope (GO:0005740) [CC]; q=0.0371307398

### C006

- * Mitochondrial Matrix (GO:0005759) [CC]; q=0.00323396766

### C007

- * Regulation Of Fibrinolysis (GO:0051917) [BP]; q=0.00431195688
- Negative Regulation Of Blood Coagulation (GO:0030195) [BP]; q=0.005389946101
- Negative Regulation Of Fibrinolysis (GO:0051918) [BP]; q=0.005389946101
- Positive Regulation Of Blood Coagulation (GO:0030194) [BP]; q=0.02530931212

### C008

- * Collagen-Containing Extracellular Matrix (GO:0062023) [CC]; q=0.00323396766

### C009

- * Platelet Alpha Granule Lumen (GO:0031093) [CC]; q=0.007699923001
- Platelet Alpha Granule (GO:0031091) [CC]; q=0.008819911801

## Interpretation and provenance

Major themes remain mitochondrial respiration/oxidative phosphorylation, fatty-acid metabolism, coagulation/fibrinolysis, platelet granules and extracellular matrix. The exact marginal catalog changes: protein secretion is no longer significant and acyl-CoA dehydrogenase beta-oxidation enters. These are statistical support themes, not activity or causal conclusions.

All 9 Route A files and all 31 saved permutation/redundancy sensitivity files match their pre-promotion SHA-256 hashes byte for byte. Route A was not rerun or edited. The old 10k production and comparison files are retained in `archive/route_b_10k_canonical_20261014/` with their manifest. All 10k/50k/100k sensitivity outputs remain unchanged.

The ten diagnostic fields are present in production all-pathway, candidate, cluster and representative views and prefixed `route_b_` in comparison/protein-evidence views. Retained review decisions are carried forward; removed selections are saved separately.

Machine-readable promotion receipt: `results/pathway_centric/DM_vs_NDM/production_promotion.json`.

## Validation

`MPLCONFIGDIR=/private/tmp/proteomics-mpl python -m pytest -q`: **96 passed in 14.42s**, no failures. This includes exact equality against saved 100k inference, protected-byte/input guards, unchanged inference arithmetic, strict cutoff boundary behavior, diagnostic consistency/propagation, deterministic representatives, candidate transitions, review preservation, and historical regression tests. Canonical Route B and comparison regeneration both completed successfully.
