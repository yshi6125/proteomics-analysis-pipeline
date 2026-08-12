# Proteomics Analysis Pipeline

## Overview

This project demonstrates a reproducible Python workflow for quality control,
normalization assessment, differential protein-abundance analysis, and pathway
enrichment of quantitative proteomics data.

The pipeline is designed to convert a protein-level data table into analysis-ready data and interpretable quality-control outputs.

## Objectives

The workflow addresses the following questions:

1. What is the structure and quality of the input dataset?
2. Are there duplicated proteins, missing values, or invalid measurements?
3. How are peptide counts and protein intensities distributed?
4. Is normalization needed across samples?
5. Which proteins differ between biological groups?
6. Which biological pathways are represented by those proteins?

## Workflow

```text
raw workbook
    -> QC and standardized metadata
    -> normalization diagnostics and approved preprocessing
    -> differential protein-abundance analysis
    -> candidate-based pathway enrichment
```

## Repository Structure

```text
data/        Example or public input data
src/         Python analysis scripts
notebooks/   Step-by-step exploratory analysis
figures/     Quality-control and result figures
results/     Processed tables and analysis outputs
```

## Tools

- Python
- pandas
- NumPy
- matplotlib
- SciPy
- scikit-learn
- statsmodels
- GSEApy

## Current Analyses

- Data structure and abundance-value QC
- Missing-value, duplicate, and peptide-count assessment
- Normalization diagnostics before and after preprocessing
- Per-protein linear models with Benjamini-Hochberg FDR correction
- Volcano plots and prioritized protein heatmaps
- Primary FDR-supported and secondary exploratory pathway enrichment

## Run QC and Normalization Diagnostics

```bash
/opt/anaconda3/bin/python src/qc.py
/opt/anaconda3/bin/python src/normalization_diagnostics.py
```

## Run Differential Abundance

```bash
/opt/anaconda3/bin/python src/differential_abundance.py
```

Differential-abundance tables are saved in `results/differential_abundance/`,
and the volcano plot and protein heatmap are saved in
`figures/differential_abundance/`.

## Run Pathway Analysis

```bash
/opt/anaconda3/bin/python src/pathway_analysis.py
```

The pathway workflow uses unique, normalized gene symbols and analyzes increased
and decreased genes separately. It supports configurable GO Biological Process,
Reactome, and KEGG libraries. Pathway significance is based on
`enrichment_adjusted_p_value < 0.05`.

Outputs in `results/pathway_analysis/` include:

- `primary_fdr_pathways_<comparison>.csv`: pathways from protein candidates that
  passed protein-level FDR (Levels 1–3). Filter the pathway adjusted p-value at
  0.05 for the primary significant-pathway list.
- `exploratory_all_candidates_pathways_<comparison>.csv`: pathways from all
  Levels 1–5, including nominal protein candidates. Significant pathways in this
  table remain exploratory.
- `pathway_candidate_members_<comparison>.csv`: long-format gene-level evidence
  linking pathways to contributing proteins, candidate origins, fold changes,
  statistics, accessions, and source row positions.
- `pathway_priority_summary_<comparison>.csv`: combined primary and exploratory
  pathway summary with overlap, Gene Ratio, Level 1–5 evidence counts, strongest
  candidate origin, and pathway priority.

Separate primary and exploratory pathway figures are saved in
`figures/pathway_analysis/`. Pathway priority describes the strongest underlying
protein evidence and does not modify enrichment p-values.

The current committed pathway outputs were generated with GSEApy's local
over-representation engine using downloaded Enrichr GO, Reactome, and KEGG
libraries because the live Enrichr submission endpoint returned HTTP 429.

## Planned Analysis

- Preranked GSEA using all tested proteins ranked by differential-abundance
  `t_statistic`

## Data Availability

The repository uses an anonymized, simulated, or publicly available example dataset. Proprietary or unpublished research data are not included.
