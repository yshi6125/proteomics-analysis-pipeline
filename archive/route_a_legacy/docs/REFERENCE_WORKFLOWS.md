# Historical/reference workflows

These alternatives are retained for inspection and reproducibility, outside the
canonical production workflow. Nothing was deleted or moved. Production Route A
is Levels 1–5, measured-background local ORA and is named simply Route A.

## Preserved paths

| Historical material | Original path |
|---|---|
| Remote A primary enrichment | `results/pathway_analysis/primary_fdr_pathways_<comparison>.csv` |
| Remote A exploratory enrichment | `results/pathway_analysis/exploratory_all_candidates_pathways_<comparison>.csv` |
| Remote A shared membership/priority | `results/pathway_analysis/` |
| Remote A clusters/evidence | `results/pathway_consolidation/primary/` and `/exploratory/` |
| Measured primary | `results/pathway_analysis_measured/<comparison>/primary/` |
| Selected pre-canonical method | `results/pathway_analysis_measured/<comparison>/exploratory/` |
| Shared measured metadata/resource/levels | `results/pathway_analysis_measured/<comparison>/` |
| Original four comparison outputs | `results/pathway_comparison/<comparison>/legacy_primary/`, `legacy_exploratory/`, `measured_primary/`, `measured_exploratory/` |
| Existing AI reports/caches/reviews | `results/pathway_reports/primary/` and `/exploratory/` |

These nested comparison folders remain beside the new standard comparison CSVs
at the parent directory; their contents are reference material. Their old naming
is preserved to avoid breaking recorded paths/hashes. Historical metadata may use
branch/mode names that are no longer production choices.

## Reproduction, explicitly reference-only

Inspect files first; these commands rerun historical writers and can replace their
original comparison-specific reference outputs. Remote results depend on Enrichr
availability/versioning and are not guaranteed identical across service changes.

```bash
# Reference remote ORA and both historical branches; internet required
python -m src.pathway_reference --stage legacy-analysis --comparison DM_vs_NDM

# Reference measured-background branches; offline
python -m src.pathway_reference --stage measured-analysis --comparison DM_vs_NDM \
  --gene-sets results/pathway_centric/DM_vs_NDM/gene_sets.json

# Reference comparison; offline, reads existing reference results
python -m src.pathway_reference --stage comparison --comparison DM_vs_NDM \
  --a-mode legacy --a-branch primary
```

For reference comparisons only, choose legacy/measured and primary/exploratory.
These switches do not exist in the production comparison CLI. The corresponding
Python reference functions remain available for tests and provenance. Old direct
script paths are retained but labeled reference-only in their module descriptions.

## Retained implementation and documentation

- `src/pathway_analysis.py`: original remote ORA writer and reusable candidate/member helpers.
- `src/pathway_consolidation.py`: reference branch writer and reusable graph/ranking functions.
- `src/pathway_measured_ora.py`: original two-branch writer and unchanged selected numerical functions.
- `src/pathway_reference.py`: explicit historical reproduction entry point.
- `src/evidence_builder.py`, report generation and agents: unchanged historical reporting schema/consumers.

Pre-canonical document snapshots (their relative links reflect original locations):

- `docs/reference/README_PRE_CANONICAL.md`
- `docs/reference/PIPELINE_ARCHITECTURE_PRE_CANONICAL.md`
- `docs/reference/PIPELINE_ONE_PAGE_PRE_CANONICAL.md`
- `docs/reference/PATHWAY_COMPARISON_GUIDE_PRE_CANONICAL.md`

`docs/PATHWAY_UPDATE_FILE_INVENTORY.md` records the earlier controlled-update files
and remains historical documentation. No old analysis or comparison result is
silently promoted/moved: canonical A is regenerated in its new tree with the same
selected numerical functions, new labels and explicit methodology metadata.

Exact pre-change source snapshots are retained as non-importable text under
`docs/reference/source/`, with original paths and SHA-256 values in
`SOURCE_MANIFEST.json`. This preserves inspection of the old dispatcher/comparison
and source fingerprints referenced by historical provenance; production uses the
current modules. For strict historical replay, use the snapshots in an isolated
checkout with the recorded inputs/dependencies rather than overwriting production
sources.
