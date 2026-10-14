# Production pathway cleanup report

## A. Components previously mixed with production

The six retired entry-point/writer modules were `src/pathway_analysis.py`,
`src/pathway_consolidation.py`, `src/pathway_measured_ora.py`,
`src/pathway_reference.py`, `src/evidence_builder.py` and
`src/generate_pathway_report.py`. They exposed remote ORA, two measured scopes,
four comparison alternatives or branch-specific reporting. Canonical A imported
numerical/annotation helpers from those modules but did **not** read their output
files. The comparison API still accepted historical mode/branch arguments.

Old result trees, four nested comparison folders, report branches, ORA figures,
reference documents and three reference test modules remained beside production.
These are now separated under the archive. The current production API has only
Route A and Route B and contains no historical selector.

## B. Exact relocation roots

Every individual file's old/new path and SHA-256 is in
`archive/route_a_legacy/ARCHIVE_MANIFEST.json`. Directory moves preserve all child
relative paths and original metadata content:

| Old repository-relative location | New repository-relative location |
|---|---|
| `results/pathway_analysis/` | `archive/route_a_legacy/results/pathway_analysis/` |
| `results/pathway_consolidation/` | `archive/route_a_legacy/results/pathway_consolidation/` |
| `results/pathway_analysis_measured/` | `archive/route_a_legacy/results/pathway_analysis_measured/` |
| `results/pathway_reports/` | `archive/route_a_legacy/results/pathway_reports/` |
| `figures/pathway_analysis/` | `archive/route_a_legacy/figures/pathway_analysis/` |
| `results/pathway_comparison/DM_vs_NDM/legacy_primary/` | `archive/route_a_legacy/results/pathway_comparison/DM_vs_NDM/legacy_primary/` |
| `results/pathway_comparison/DM_vs_NDM/legacy_exploratory/` | `archive/route_a_legacy/results/pathway_comparison/DM_vs_NDM/legacy_exploratory/` |
| `results/pathway_comparison/DM_vs_NDM/measured_primary/` | `archive/route_a_legacy/results/pathway_comparison/DM_vs_NDM/measured_primary/` |
| `results/pathway_comparison/DM_vs_NDM/measured_exploratory/` | `archive/route_a_legacy/results/pathway_comparison/DM_vs_NDM/measured_exploratory/` |
| Six retired `src/` modules named above | `archive/route_a_legacy/src/<same filename>` |
| `tests/test_pathway_review.py` | `archive/route_a_legacy/tests/test_pathway_review.py` |
| `tests/test_canonical_route_a.py` | `archive/route_a_legacy/tests/test_canonical_route_a.py` |
| `tests/test_pathway_centric.py` | `archive/route_a_legacy/tests/test_pathway_centric.py` |
| `docs/reference/` | `archive/route_a_legacy/docs/reference/` |
| `docs/REFERENCE_WORKFLOWS.md` | `archive/route_a_legacy/docs/REFERENCE_WORKFLOWS.md` |
| `docs/PATHWAY_UPDATE_FILE_INVENTORY.md` | `archive/route_a_legacy/docs/PATHWAY_UPDATE_FILE_INVENTORY.md` |
| `docs/CANONICAL_UPDATE_FILE_INVENTORY.md` | `archive/route_a_legacy/docs/CANONICAL_UPDATE_FILE_INVENTORY.md` |

Before moving/removing active copies, all pre-cleanup source, test/document packages,
processed input data and protein statistics were copied and hash-verified in the
archive. Production result snapshots are in its `baseline/results/`. The manifest
contains 148 preserved file records. No historical metadata identity was renamed.

## C. Active source and output locations

| Component | Active sources | Active outputs |
|---|---|---|
| Route A | `src/pathway_protein.py`, `src/pathway_ora.py`, `src/protein_evidence.py` | `results/pathway_protein/DM_vs_NDM/` |
| Route B | `src/pathway_centric.py` | `results/pathway_centric/DM_vs_NDM/` |
| Comparison | `src/pathway_comparison.py` | Flat files in `results/pathway_comparison/DM_vs_NDM/` |
| Shared graph/ranking | `src/pathway_redundancy.py` | Used independently by both routes |
| Shared IO/provenance | `src/pathway_review_io.py` | Validation and hashes |
| Dispatcher | `src/pathway_routes.py` | Calls only production A/B |

Shared QC, diagnosis, preprocessing and DA remain active. Reporting agent classes
are retained unchanged as reusable libraries; their branch-specific writer/evidence
entry points and prior artifacts are archived. No new reporting integration exists.

## D. Dependency verification

Canonical A calculates local ORA from the complete successfully tested differential
CSV plus local JSON definitions. It is **not** a copy/reader of an old scope's
results. Numerical functions have been extracted into production-only helper
modules with the same candidate scope, BH family, universe and ranking rules.

Route B's function ASTs are unchanged; only its imports point to the extracted
loader and graph helpers. Standard comparison reads only canonical A/B results,
the shared DA table and their metadata/snapshots. Its historical mode/branch
parameters and paths have been removed entirely.

An isolated installation containing only production source modules and synthetic
inputs runs A, B and comparison with neither an archive nor obsolete result folders.
Thus the archive is not a runtime dependency. Test/reference inspection of the
archive is deliberate and separate from production execution.

## E. Exhaustive search audit

`docs/CLEANUP_SEARCH_AUDIT.tsv` records every matched occurrence of the requested
terms in textual repository files and paths, with line/column location,
classification and rationale. `docs/CLEANUP_SEARCH_SUMMARY.json` records totals
and scan scope. Generated audit/report files themselves are excluded to avoid
recursive self-reporting. Regenerate with `python scripts/audit_pathway_terms.py`. Git internals, caches and binary files are excluded;
archive content and dotfiles are included.

- **HISTORICAL/ARCHIVED:** archive content/path references and immutable
  generation-time method references in saved canonical metadata. Original source
  paths referenced by those records resolve through the archive manifest.
- **TEST/REFERENCE:** tests, archive-inspection documentation and the controlled
  cleanup explanation. They expose no production choices.
- **ACTIVE:** shared DA protein-tier diagnostics (including nominal Level 4–5
  views) biological protein/pathway annotations, shared QC interpretation, and unchanged literature/dataset-strength classification. These are not retired Route A
  pathway workflows, and changing them would exceed this cleanup's scope.

Search counts: 71780 occurrences; 71554 HISTORICAL/ARCHIVED, 52 TEST/REFERENCE, 174 intentional ACTIVE occurrences.

No production A/B/comparison module contains the retired terms or imports an
archived workflow. No old branch result directories remain in active pathway
results. Saved canonical metadata is not rewritten to falsify what implementation
originally generated it; the manifest provides its source relocation context.

## F. Numerical/provenance checks

Production A, B and standard comparison result files remain byte-for-byte identical
to their before-cleanup snapshots. No production statistical command was run during
cleanup. Synthetic-input tests compare the extracted A tables, candidates, cluster
IDs, representatives, memberships, universe and levels exactly with the preserved
pre-cleanup implementation. B's complete function ASTs and frozen numerical/cluster
fixture remain identical. Comparison regenerated in temporary directories has
exactly the saved production table, and preserves inputs and user selections.

Future newly generated metadata correctly references the extracted production
modules. Existing generation-time metadata remains unchanged; original source
hashes and paths are inspectable through the archive's preserved source snapshots.

## G. Commands and locations going forward

```bash
python -m src.pathway_protein --comparison DM_vs_NDM
python -m src.pathway_routes --route pathway --comparison DM_vs_NDM \
  --gene-sets results/pathway_centric/DM_vs_NDM/gene_sets.json
python -m src.pathway_comparison --comparison DM_vs_NDM
```

All are local with the established snapshot. These are the approved commands;
they were verified through isolated tests/help, not used to overwrite existing
production statistical results during cleanup.

Protein evidence:
`results/pathway_comparison/DM_vs_NDM/pathway_protein_evidence.csv`.
Metadata: A/B `analysis_metadata.json` in their respective trees, and
`results/pathway_comparison/DM_vs_NDM/comparison_metadata.json`.
Archive: `archive/route_a_legacy/`, including `REPRODUCTION.md` and
`ARCHIVE_MANIFEST.json`.

## H. Tests

Normal collection is explicitly `tests/` via `pytest.ini`. The 26 former tests
for mixed/reference workflows are archived, not discarded. Fifteen focused
production pathway tests replace their mixed entry-point tests. Other upstream,
agent and review tests remain active; a provider-construction test deliberately
imports the archived report script as a TEST/REFERENCE check, never as a production
runtime dependency. The active suite passed: **78 tests**. Previous mixed-workflow suites are preserved in the archive. All archive hashes and production snapshot hashes
are checked by tests.
