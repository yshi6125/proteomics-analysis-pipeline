# Separated historical workflow archive

The active pipeline contains one Route A (Levels 1–5, measured-background local
ORA), one unchanged Route B and their comparison. The directories and modules in
this archive are historical/reference material, never imported/read by production.

## Contents and original identities

`ARCHIVE_MANIFEST.json` records every preserved file's original repository-relative
path, current archived path, SHA-256 and move/copy operation. Original result
metadata and terminology have not been rewritten. The original absolute paths in
historical provenance resolve via this manifest, relative to the former project
root. Previously generated canonical metadata may also refer to the archived
pre-cleanup source files; it remains truthful generation-time provenance.

- `src/`: complete pre-cleanup source package, including all former A writers,
  four-mode comparison/reference handling, and unchanged reporting scripts/agents.
- `tests/`, `docs/`, `prompts/`, `README_PRE_CLEANUP.md`, `requirements.txt`:
  pre-cleanup tests, documentation, prompts and dependencies.
- `data/processed/`, `results/differential_abundance/`: copies of shared input and
  protein statistics available at cleanup time.
- `results/pathway_analysis/`, `results/pathway_consolidation/`,
  `results/pathway_analysis_measured/`, `results/pathway_reports/`: moved former
  branch analysis, evidence and report artifacts.
- `results/pathway_comparison/<comparison>/<former_variant>/`: moved original
  four-variant comparisons.
- `figures/pathway_analysis/`: moved former ORA figures.
- `baseline/results/`: byte-for-byte snapshots of production A, B and standard
  comparison at cleanup time, used for regression/provenance checks.

## Inspect without running

Use the recorded original identity and metadata when interpreting a result. Check
its SHA-256 against the manifest before copying or using it for a regression.
Do not point a production command at an old branch folder. There is no active
historical dispatcher or historical mode in the production comparison API.

## Reconstruct in an isolated workspace

1. Create a separate empty working directory outside the production result trees.
2. Copy this archive's `src/`, `tests/`, `docs/`, `prompts/`, `data/`,
   `requirements.txt` and `results/` into it at those original relative paths.
3. Copy `baseline/results/pathway_protein/` and `baseline/results/pathway_centric/`
   to that workspace's `results/`, and copy the flat standard comparison files
   from `baseline/results/pathway_comparison/` without replacing variant folders.
4. Use the versions recorded in the applicable result metadata. Dependencies in
   requirements are ranges, not an environment lock. The preserved raw input
   results are sufficient for downstream replay; raw Excel regeneration additionally
   requires the original workbook (not duplicated here).
5. Inspect `README_PRE_CLEANUP.md` and preserved source CLI help. Explicit historical
   commands such as `python -m src.pathway_reference --stage comparison --a-mode
   measured --a-branch exploratory` run only in this reconstructed old layout.

Historical remote Enrichr reruns need internet and may differ if its service or
resources change. Local snapshot-based ORA and saved-result comparisons can run
offline. Writers overwrite reconstructed output files; keep the immutable archive
as the reference. Some older tests use `git show HEAD:src/pathway_analysis.py`;
strict test replay therefore also needs the original Git baseline, or an explicit
baseline supplied from the preserved source snapshots. Production tests no longer
rely on those old source/result locations.
