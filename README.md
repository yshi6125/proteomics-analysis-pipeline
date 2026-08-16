# Proteomics Analysis Pipeline

## Overview

This repository contains a reproducible V1 workflow for quantitative proteomics,
from input quality control through disease-specific, evidence-grounded pathway
reports. Deterministic analysis produces structured pathway evidence before any
LLM is used.

## End-to-End Workflow

```text
Proteomics input
  -> QC and normalization diagnostics
  -> diagnosis-driven normalization and batch correction
  -> differential protein analysis
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

- Preprocessing is diagnosis-driven: diagnostics inform explicit normalization
  and batch-correction decisions rather than applying automatic transformations.
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
preprocessing functions; apply the approved transformations after reviewing the
diagnostics, then save the final analysis matrix used downstream.

```bash
python src/qc.py
python src/normalization_diagnostics.py
python src/differential_abundance.py
python src/pathway_analysis.py
python src/pathway_consolidation.py
python src/evidence_builder.py
```

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
`results/pathway_reports/<branch>/`; generated reports, review records, and caches
are intentionally not tracked.

## Environment Variables

- `GEMINI_API_KEY`: required when `--provider gemini` is selected.
- `OPENAI_API_KEY`: required when `--provider openai` is selected.
- `NCBI_EMAIL`: optional contact address included in NCBI E-utilities requests.
- `NCBI_API_KEY`: optional; enables the higher NCBI request rate. Without it,
  requests are conservatively throttled to comply with the lower public rate.

Do not commit environment-variable values or local credential files.

## Key Outputs

- `results/differential_abundance/`: protein statistics and Level 1-5 candidate
  summaries.
- `results/pathway_analysis/`: primary and exploratory enrichment results plus
  pathway-member evidence.
- `results/pathway_consolidation/<branch>/`: support-filtered clusters,
  consolidation metadata, and pathway evidence JSON.
- `results/pathway_reports/<branch>/`: generated reports and reviewer provenance
  at runtime.

Example output: [exploratory Type 2 diabetes pathway report](results/pathway_reports/exploratory/dataset_report_DM_vs_NDM.md).

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
