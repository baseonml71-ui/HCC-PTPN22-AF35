# Clonotype context and temporal linkage

## Purpose
Clonotype context and temporal linkage. This is the single family review starting point.

## Primary script
None approved for public execution. The 10 source candidates below are excluded from ACTIVE_CURRENT pending complete input/source-result contracts. Their ordering is not an executable pipeline.

- [common.py](sources/AF35_CLONOTYPE_NEUTRALIZED_STATE_OCCUPANCY_v1/scripts/common.py) — A010
- [03_phaseA_molecular.py](sources/AF35_LONGITUDINAL_CLONAL_STATE_DYNAMICS_M2_CLEANROOM_v1/scripts/03_phaseA_molecular.py) — A016
- [07_phaseB_verify.R](sources/AF35_LONGITUDINAL_CLONAL_STATE_DYNAMICS_M2_CLEANROOM_v1/scripts/07_phaseB_verify.R) — A018
- [common.py](sources/AF35_LONGITUDINAL_CLONAL_STATE_DYNAMICS_M2_CLEANROOM_v1/scripts/common.py) — A019
- [02_phase19a_patient_analysis.py](sources/GSE235863_CLONE_STATE_DECOMPOSITION_v1/authority/04_scripts/phase19/02_phase19a_patient_analysis.py) — A026
- [01_authority.py](sources/GSE235863_CLONE_STATE_DECOMPOSITION_v1/scripts/01_authority.py) — A027
- [02_analyze.py](sources/GSE235863_CLONE_STATE_DECOMPOSITION_v1/scripts/02_analyze.py) — A028
- [01_meta_clone_compartment.R](sources/HCC_ICI_project/04_scripts/final_harmonization/01_meta_clone_compartment.R) — A033
- [03_independent_scrna_sctcr_clone_state.R](sources/HCC_ICI_project/04_scripts/phase3/03_independent_scrna_sctcr_clone_state.R) — A042
- [01_gse235863_baseline_future_clone_fate.R](sources/HCC_ICI_project/04_scripts/phase_temporal/01_gse235863_baseline_future_clone_fate.R) — A043

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
