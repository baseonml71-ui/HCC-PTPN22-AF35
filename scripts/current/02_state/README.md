# Patient-level state architecture and covariance

## Purpose
Patient-level state architecture and covariance. This is the single family review starting point.

## Primary script
None approved for public execution. The 8 source candidates below are excluded from ACTIVE_CURRENT pending complete input/source-result contracts. Their ordering is not an executable pipeline.

- [config.R](sources/HCC_ICI_project/04_scripts/config.R) — A032
- [07_functional_program_GSE206325.R](sources/HCC_ICI_project/04_scripts/phase2/07_functional_program_GSE206325.R) — A041
- [01_source_audit.py](sources/PTPN22_AF35_COUPLING_ARCHITECTURE_v1/scripts/01_source_audit.py) — A046
- [05_run_analysis.py](sources/PTPN22_AF35_COUPLING_ARCHITECTURE_v1/scripts/05_run_analysis.py) — A047
- [coupling_math.py](sources/PTPN22_AF35_COUPLING_ARCHITECTURE_v1/scripts/coupling_math.py) — A048
- [00_reuse_sources.py](sources/PTPN22_STATE_AXIS_SPECIFICITY_v1/scripts/00_reuse_sources.py) — A054
- [02_analyze.py](sources/PTPN22_STATE_AXIS_SPECIFICITY_v1/scripts/02_analyze.py) — A055
- [coupling_math.py](sources/PTPN22_STATE_AXIS_SPECIFICITY_v1/scripts/coupling_math.py) — A056

## Required inputs
See [derived input ledger](../../../release/DERIVED_INPUT_LEDGER.tsv), [dataset manifest](../../../data/SOURCE_DATA_MANIFEST.tsv) and [parameter records](../../../config/ANALYSIS_PARAMETER_MANIFEST.tsv). The dataset path example is not a completed execution configuration. Do not bypass scientific freeze checks or replace absent protocol content from memory.

## Expected outputs
See [source-result map](../../../release/PUBLIC_ANALYSIS_SOURCE_MAP.tsv). Included output tables are evidence, not proof of upstream closure.

## Environment
See [historical environment records](../../../environment/README.md). No complete locked environment has been demonstrated.

## Approximate execution order
First resolve exact acquisition and intermediate-production contracts; then verify molecular/outcome-blind freeze guards; only then perform the authorized downstream analysis. A complete family-specific runnable order cannot yet be asserted. No biological rerun is part of RC2.

## Known limitations
See [RC2 boundary](../../../docs/RC2_PACKAGING_BOUNDARY.md) and [unresolved actions](../../../release/UNRESOLVED_AUTHOR_ACTIONS.tsv). Relocation does not repair implicit working-directory or sibling-helper assumptions. Original module layout can be inspected through the code-only materializer; this does not make analyses runnable.
