# Proteomics Analysis Pipeline

## Overview

This repository contains a reproducible V1 workflow for quantitative proteomics,
from input quality control through disease-specific, evidence-grounded pathway
reports. Deterministic analysis produces structured pathway evidence before any
LLM is used.

## End-to-End Workflow

```text
Proteomics input
  -> QC, missingness, and transformation diagnostics
  -> diagnosis-driven normalization and log2 transformation
  -> PCA plus group/batch diagnosis
  -> differential protein analysis (abundance ~ group + batch)
  -> candidate prioritization (Levels 1-5)
  -> pathway enrichment
  -> pathway support filtering
  -> Jaccard network construction and Louvain redundancy reduction
  -> evidence aggregation
  -> separate primary and exploratory branches
  -> disease-specific PubMed retrieval (previous five years)
  -> batched LLM literature assessment and interpretation
  -> deterministic dataset-evidence and disease-relevance classification
  -> scientific reviewer and bounded correction pass
  -> final pathway report
```

The analysis writes processed data and audit tables under `data/processed/` and
`results/`, with diagnostic figures under `figures/`. Primary evidence uses
Levels 1-3; the exploratory branch additionally includes Levels 4-5. Pathway
redundancy is reduced independently within each branch using candidate-gene
Jaccard similarity, a pathway network, and seeded Louvain clustering. The
resulting evidence JSON is the fixed input to the reporting layer.

## Design Principles

- Preprocessing is diagnosis-driven: diagnostics inform explicit normalization,
  transformation, and statistical-design decisions rather than applying automatic
  transformations.
- Batch-associated variation is diagnosed before analysis and batch is modeled
  directly as a nuisance covariate in differential abundance. Batch effects are
  not subtracted from the production analysis matrix before inference.
- QC, differential analysis, enrichment, clustering, and evidence aggregation
  remain deterministic and precede LLM interpretation.
- Primary and exploratory evidence are never combined or used to upgrade one
  another.
- Dataset Evidence Strength is independent of literature-based Disease Relevance;
  both classifications are computed in Python using fixed rules.
- Enrichment and protein abundance do not imply pathway activation, inhibition,
  metabolic flux, or causality.
- AI interpretation is restricted to supplied evidence and checked by a bounded
  Scientific Reviewer Agent before the final report is saved.

## Installation

```bash
python -m pip install -r requirements.txt
```

## Usage

Run the deterministic stages in order. `normalization.py` provides explicit
preprocessing functions; apply the approved normalization and log2 transformation
after reviewing the diagnostics. Differential abundance reads
`data/processed/data_log2.csv` and fits the biological contrast and preparation
batch simultaneously. Batch is included as a nuisance covariate; no
batch-corrected abundance matrix is created or used for inference.

```bash
python src/qc.py
python src/diagnosis.py --comparison DM_vs_NDM
python src/differential_abundance.py
python src/pathway_analysis.py
python src/pathway_consolidation.py
python src/evidence_builder.py
```

Each diagnosis execution is saved without overwriting prior runs under
`results/diagnostics/diagnosis_###/`, including a deterministic Markdown report,
the figures and tables produced by that run, and a `latest_diagnosis.json`
pointer. Diagnosis history is separate from downstream production outputs.

Diagnosis/QC is intentionally versioned because preprocessing diagnostics may be
run repeatedly while evaluating transformations, scaling, batch effects, and
potential outliers before committing to the final analysis strategy:

```text
results/diagnostics/
├── diagnosis_001/
│   ├── diagnosis_report_DM_vs_NDM.md
│   ├── figures/
│   └── tables/
├── diagnosis_002/
│   ├── diagnosis_report_DM_vs_NDM.md
│   ├── figures/
│   └── tables/
├── diagnosis_003/
│   ├── diagnosis_report_DM_vs_NDM.md
│   ├── figures/
│   └── tables/
├── diagnosis_004/
│   ├── diagnosis_report_DM_vs_NDM.md
│   ├── figures/
│   └── tables/
└── latest_diagnosis.json
```

These directories are exploratory QC and preprocessing records. Their numbering
does not create new versions of the final differential-abundance, pathway,
evidence, or AI-report outputs; production analysis remains separate and begins
only after the preprocessing strategy has been approved.

For datasets with two or more valid preparation batches, differential abundance
fits `abundance ~ group + batch`; with one batch it records and uses the
group-only fallback. Perfect group/batch confounding is rejected. Future studies
may require additional covariates selected for their specific clinical design;
the pipeline does not add age, sex, BMI, or other covariates automatically.

The final two commands create independent primary and exploratory cluster tables
and `pathway_evidence_<comparison>.json` files. Their default comparison is
`DM_vs_NDM`.

Preview the disease-specific PubMed queries without searching PubMed, generating
a report, or invoking the reviewer:

```bash
python src/generate_pathway_report.py \
  --comparison DM_vs_NDM \
  --disease "Type 2 diabetes" \
  --branch exploratory \
  --provider gemini \
  --model <supported-gemini-model> \
  --preview-queries
```

Generate a reviewed report with Gemini and four clusters per literature batch:

```bash
python src/generate_pathway_report.py \
  --comparison DM_vs_NDM \
  --disease "Type 2 diabetes" \
  --branch exploratory \
  --provider gemini \
  --model <supported-gemini-model> \
  --literature-batch-size 4
```

Use `--provider openai --model <openai-model>` for OpenAI. `--refresh-cache`
bypasses cached themes, literature assessments, and reviews. `--skip-review`
disables scientific review for development or debugging. Reports are written to
`results/pathway_reports/<branch>/`. Final reports may be committed as example
outputs; caches, temporary drafts, and reviewer working records are not tracked.

## Environment Variables

- `GEMINI_API_KEY`: required when `--provider gemini` is selected.
- `OPENAI_API_KEY`: required when `--provider openai` is selected.
- `NCBI_EMAIL`: optional contact address included in NCBI E-utilities requests.
- `NCBI_API_KEY`: optional; enables the higher NCBI request rate. Without it,
  requests are conservatively throttled to comply with the lower public rate.

Do not commit environment-variable values or local credential files.

## Key Outputs

- `results/diagnostics/diagnosis_###/`: immutable diagnosis reports with the
  figures and tables generated during each QC/preprocessing assessment.
- `results/differential_abundance/`: protein statistics and Level 1-5 candidate
  summaries.
- `results/pathway_analysis/`: primary and exploratory enrichment results plus
  pathway-member evidence.
- `results/pathway_consolidation/<branch>/`: support-filtered clusters,
  consolidation metadata, and pathway evidence JSON.
- `results/pathway_reports/<branch>/`: generated reports and reviewer provenance
  at runtime.

Differential-abundance figures are organized into two evidence views:

```text
Differential Abundance
├── Primary candidates (Levels 1–3)
│   ├── figures/differential_abundance/volcano_<comparison>.png
│   └── figures/differential_abundance/top_protein_heatmap_<comparison>.png
└── Exploratory candidates (exclusive Levels 4–5)
    ├── figures/differential_abundance/exploratory/volcano_<comparison>.png
    ├── figures/differential_abundance/exploratory/heatmap_all_<comparison>.png
    └── figures/differential_abundance/exploratory/heatmap_top_<N>_<comparison>.png
```

Primary figures emphasize higher-confidence candidates. Exploratory figures
show nominal Level 4–5 candidates retained for hypothesis generation and
exploratory pathway analysis. When a top-N exploratory heatmap is produced, it
ranks Level 4 before Level 5, then smaller p-value, larger absolute log2FC, and a
stable gene/accession tie-break; the complete heatmap still contains every
exclusive Level 4–5 candidate.

Example outputs: [primary Type 2 diabetes pathway report](results/pathway_reports/primary/dataset_report_DM_vs_NDM.md)
and [exploratory Type 2 diabetes pathway report](results/pathway_reports/exploratory/dataset_report_DM_vs_NDM.md).

## Limitations

- Many exploratory pathways may be supported mainly by lower-tier candidates.
- Literature retrieval depends on the biological theme and disease terminology
  used in PubMed queries.
- Protein abundance and enrichment do not establish pathway activity, metabolic
  flux, or causality.
- AI-generated interpretations are evidence-grounded and reviewer-checked, but
  still require expert scientific review.

## Data Availability

The repository uses an anonymized, simulated, or publicly available example
dataset. Proprietary or unpublished research data are not included.

## Route B: independent pathway-centric analysis

The existing commands, Route A candidate thresholds, primary/exploratory outputs,
consolidation, and downstream reports remain unchanged. A new dispatcher runs
one or both routes after differential analysis:

```bash
python -m src.pathway_routes --route pathway --comparison DM_vs_NDM
python -m src.pathway_routes --route both --comparison DM_vs_NDM
python -m src.pathway_routes --route protein --comparison DM_vs_NDM
```

Route B defaults to Enrichr GO Biological Process and Cellular Component 2023
libraries, obtained with the existing GSEApy dependency. Use `--libraries` and
`--organism` for another collection; unavailable libraries fail explicitly.
For offline execution or a pinned database, supply `--gene-sets sets.json`, with
JSON shaped as `{"library": {"pathway": ["GENE1", "GENE2"]}}`.

Route B reads the complete differential table and includes every row with a
usable gene symbol and finite group-effect `t_statistic`, regardless of protein
p-value, q-value, fold change, or candidate level. Missing gene symbols and invalid
t-statistics are recorded in an eligibility audit. Gene symbols are trimmed and
uppercased, matching Route A's mapping convention. Repeated gene symbols contribute
the arithmetic mean of their eligible protein t-statistics, avoiding selection by
significance or counting a gene more than once. All contributing original protein
rows are retained in memberships. Accessions alone cannot map to these gene-symbol
libraries and are audited as missing gene symbols.

The test statistic is mean signed t. Its competitive null comprises uniformly
sampled subsets of the eligible measured gene universe, without replacement,
matched to each pathway's measured size. Default: 10,000 permutations with seed
42. Two-sided p-values are twice the smaller inclusive tail probability, capped
at one, with plus-one Monte Carlo correction. BH adjustment covers every tested
term across all requested libraries in a single family. Change `--permutations`
for better tail resolution (minimum possible two-sided p is 2/(permutations+1)).
Size-matched nulls are shared and seeded by size, so library ordering does not
change results. Size limits apply to measured genes, not total database membership.

This is inspired by PIANO mean-statistic/gene-sampling analysis, as described in
[the PIANO vignette](https://bioconductor.org/packages/release/bioc/vignettes/piano/inst/doc/piano-vignette.pdf),
but is not an exact reproduction of the reference study. Existing model t-values
are ordinary, not moderated. Gene sampling assumes gene exchangeability and does
not preserve correlations between proteins; pathway p-values are competitive
relative to the measured background, rather than a correlation-adjusted test of
no group effect. Opposite-sign changes can cancel in the signed mean. Positive
and negative labels reflect the mean t sign, while `mixed_signs` explicitly marks
pathways containing both signs. These labels do not imply activation/inhibition.

Options: `--min-size 5`, `--max-size 500`, `--q-cutoff 0.05`, `--p-cutoff 1`,
`--permutations 10000`, `--seed 42`, `--jaccard-threshold 0.25`, and
`--louvain-seed 42`. Pathway candidates meet both inclusive p/q cutoffs.
The dispatcher options apply to Route B; Route A uses its unchanged defaults.

Route B writes only under `results/pathway_centric/<comparison>/`:

- `all_pathways.csv`: every tested term and its mean t, p, q, and direction.
- `pathway_candidates.csv`: selected terms.
- `pathway_members.csv`: all measured contributing protein rows for tested terms.
- `protein_eligibility.csv` and `gene_set_audit.csv`: input and size exclusions.
- `pathway_clusters.csv`: every candidate with its cluster and representative flag.
- `representative_pathways.csv`: one representative per cluster.
- `gene_sets.json`: full normalized database snapshot, including unmeasured genes.
- `analysis_metadata.json`: parameters, counts, limitations, and input/snapshot hashes.

Clustering reuses Route A's Jaccard network and seeded Louvain functions, using
all measured genes of selected pathways. Directions are clustered separately.
Isolates are retained. Representatives rank by lowest q, then p, then library and
term name; all redundant terms remain in the cluster table. No LLM biological
theme naming or cross-route evidence integration is added. Existing report builders
continue to consume Route A evidence only.

Cutoff changes can reuse saved `all_pathways.csv` and `pathway_members.csv` with
`select_candidates` and `consolidate_candidates` in `src.pathway_centric`, without
repeating differential analysis or pathway testing. Empty candidate sets produce
valid empty tables. Rerunning the same comparison replaces its Route B files;
use a separate comparison label to retain alternate parameter runs.

Validate with `python -m pytest -q`. Route A regression checks execute the Git
baseline and current implementations on the same mock enrichment fixture and
compare numerical tables, metadata, and CSV bytes; Route B checks cover eligibility,
BH, memberships, fixed-seed reproducibility, representative ties, and route dispatch.

## Controlled Route A background and A/B review checkpoint

Legacy Route A remains the default; Route B statistics and AI reporting are
unchanged. An optional local measured-background ORA variant writes separately:

```bash
python -m src.pathway_measured_ora --comparison DM_vs_NDM \
  --gene-sets results/pathway_centric/DM_vs_NDM/gene_sets.json
python -m src.pathway_comparison --comparison DM_vs_NDM \
  --a-mode measured --a-branch exploratory
```

Both commands are offline with the saved snapshot. Revised A results go to
`results/pathway_analysis_measured/DM_vs_NDM/`; review tables go to
`results/pathway_comparison/DM_vs_NDM/measured_exploratory/`. Use `--a-mode legacy`
to inspect the saved original ORA outputs. Comparisons preserve route statistics,
flag incompatible coverage, and stop at human review; no automatic reporting or
evidence integration is added. See [the comparison guide](docs/PATHWAY_COMPARISON_GUIDE.md)
for universe definition, matching/counting rules, provenance, limitations and the
ordered commands.
