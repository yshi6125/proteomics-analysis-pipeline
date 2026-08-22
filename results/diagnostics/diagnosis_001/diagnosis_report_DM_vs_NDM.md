# Diagnosis Report: DM_vs_NDM

## Run Metadata

- Run ID: diagnosis_001
- Generated at: 2026-08-18T21:31:53.429542+00:00
- Input data: `data/processed/data_qc.csv`
- Sample metadata: `data/processed/sample_metadata.csv`
- Input dimensions: 1965 rows × 30 columns
- Samples: 20
- Groups: {"DM": 10, "NDM": 10}
- Batches: {"B1": 10, "B2": 10}

## Missing Values and Imputation

- Missing abundance values: 0 (0.000%)
- Analysis dataframe imputed: False
- Temporary values imputed for PCA only: 0

## Distribution and Preprocessing Assessment

- Abundance scale: linear
- Median sample skewness: 26.0729
- Strongly skewed sample proportion: 1
- Distribution interpretation: Most samples are strongly right-skewed; a log2 transformation may be worth evaluating before PCA and correlation analysis. No transformation was applied by this diagnostic function. Sample medians differ substantially, which may indicate a global intensity shift. This plot alone cannot distinguish technical bias from a true global biological effect.
- Log2 transformation already applied: False
- Normalization/scaling applied by this diagnostic run: False

## Batch Assessment

- Recommendation: Consider batch correction
- Rationale: Batch shows statistical association with a leading component or uniquely explains at least 10% of the weighted variation within the PC1/PC2 subspace.
- Batch correction applied by this run: False
- Unique batch variance in PC1/PC2: 0.24966

## PCA and QC Assessment

- PCA features initially available: 1965
- PCA features used: 1965
- PC1 variance explained: 0.707042
- PC2 variance explained: 0.216088
- Median pairwise correlation: 0.969304
- Potential low-correlation samples: ["NDM_B1_01", "NDM_B1_04", "NDM_B2_07", "DM_B1_13", "DM_B1_15"]
- Potential outliers: []
- Outlier rule: At least 2 of 4 diagnostic flags
- Clustering settings: metric=correlation, linkage=average

## Warnings

- Batch assessment on linear-scale data should be repeated after an approved log2 transformation before a final correction decision.

## Run Artifacts

### Figures

- [figures/abundance_boxplots.png](figures/abundance_boxplots.png)
- [figures/abundance_density.png](figures/abundance_density.png)
- [figures/batch_assessment.png](figures/batch_assessment.png)
- [figures/hierarchical_clustering_dendrogram.png](figures/hierarchical_clustering_dendrogram.png)
- [figures/pca_outlier_flags.png](figures/pca_outlier_flags.png)
- [figures/pca_samples.png](figures/pca_samples.png)
- [figures/sample_correlation_heatmap_pearson.png](figures/sample_correlation_heatmap_pearson.png)
- [figures/sample_medians.png](figures/sample_medians.png)

### Tables

- [tables/batch_assessment.csv](tables/batch_assessment.csv)
- [tables/pca_coordinates.csv](tables/pca_coordinates.csv)
- [tables/sample_correlation_pearson.csv](tables/sample_correlation_pearson.csv)
- [tables/sample_distribution_summary.csv](tables/sample_distribution_summary.csv)
- [tables/sample_median_correlations.csv](tables/sample_median_correlations.csv)
- [tables/sample_outlier_diagnostics.csv](tables/sample_outlier_diagnostics.csv)

## Diagnostic Scope

This run is deterministic and diagnostic only. It does not transform, normalize, impute, batch-correct, exclude samples, or run downstream analyses.
