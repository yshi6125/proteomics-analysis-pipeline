# Diagnosis Report: DM_vs_NDM

## Run Metadata

- Run ID: diagnosis_004
- Generated at: 2026-08-22T20:11:54.921779+00:00
- Input data: `data/processed/data_log2.csv`
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

- Abundance scale: log2
- Median sample skewness: 0.745656
- Strongly skewed sample proportion: 0
- Distribution interpretation: The distributions are not consistently strongly right-skewed across samples. Sample medians span less than 1 log2 unit, with no clear large global shift by this diagnostic.
- Log2 transformation already applied: True
- Normalization/scaling applied by this diagnostic run: False

## Batch Assessment

- Observed QC finding: batch-associated variation assessed from PCA and metadata
- Statistical handling: Include batch as a model covariate
- Rationale: Batch is strongly associated with at least one leading component and uniquely explains at least 20% of the weighted variation within the PC1/PC2 subspace.
- Batch subtraction applied by this run: False
- Unique batch variance in PC1/PC2: 0.481755

## PCA and QC Assessment

- PCA features initially available: 1965
- PCA features used: 1965
- PC1 variance explained: 0.162191
- PC2 variance explained: 0.0884615
- Median pairwise correlation: 0.877287
- Potential low-correlation samples: []
- Potential outliers: []
- Outlier rule: At least 2 of 4 diagnostic flags
- Clustering settings: metric=correlation, linkage=average

## Warnings

- No additional deterministic warnings.

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

This run is deterministic and diagnostic only. It does not transform, normalize, impute, subtract batch effects, exclude samples, or run downstream analyses.
