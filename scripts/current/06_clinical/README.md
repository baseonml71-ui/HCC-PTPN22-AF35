# Pretreatment clinical association and sensitivity

## Purpose
Pretreatment clinical association and sensitivity. This is the single family review starting point.

## Primary script
None approved for public execution. The 22 source candidates below are excluded from ACTIVE_CURRENT pending complete input/source-result contracts. Their ordering is not an executable pipeline.

- [02_clinical_utf8.py](sources/AF35_AGGREGATION_COMPOSITION_DECOMPOSITION_v1/scripts/02_clinical_utf8.py) — A002
- [00_initialize.py](sources/AF35_BROAD_IMMUNE_RESIDUAL_ASSOCIATION_v1/scripts/00_initialize.py) — A004
- [03_join_and_analyze.py](sources/AF35_BROAD_IMMUNE_RESIDUAL_ASSOCIATION_v1/scripts/03_join_and_analyze.py) — A005
- [02_join.py](sources/AF35_CLINICAL_STATE_OCCUPANCY_DECOMPOSITION_v1/scripts/02_join.py) — A006
- [03_clinical.py](sources/AF35_CLINICAL_STATE_OCCUPANCY_DECOMPOSITION_v1/scripts/03_clinical.py) — A007
- [common.py](sources/AF35_CLINICAL_STATE_OCCUPANCY_DECOMPOSITION_v1/scripts/common.py) — A008
- [02_clinical.py](sources/AF35_CLONOTYPE_NEUTRALIZED_STATE_OCCUPANCY_v1/scripts/02_clinical.py) — A009
- [01_freeze.py](sources/AF35_COMPONENT_ABLATION_ASSOCIATION_ROBUSTNESS_v1/01_freeze.py) — A011
- [05_analyze.py](sources/AF35_COMPONENT_ABLATION_ASSOCIATION_ROBUSTNESS_v1/05_analyze.py) — A012
- [06_phaseB_clinical.py](sources/AF35_LONGITUDINAL_CLONAL_STATE_DYNAMICS_M2_CLEANROOM_v1/scripts/06_phaseB_clinical.py) — A017
- [02_freeze_matched_programs.py](sources/AF35_MATCHED_RANDOM_CLINICAL_NULL_v1/02_freeze_matched_programs.py) — A020
- [03_score_and_join.R](sources/AF35_MATCHED_RANDOM_CLINICAL_NULL_v1/03_score_and_join.R) — A021
- [04_summarize_and_verify.py](sources/AF35_MATCHED_RANDOM_CLINICAL_NULL_v1/04_summarize_and_verify.py) — A022
- [03_paired_and_firth.R](sources/CLINICAL_ANCHOR_READOUT_DECOMPOSITION_v1/scripts/03_paired_and_firth.R) — A025
- [01_prepare_response_blind_inputs.R](sources/HCC_ICI_project/04_scripts/phase15b_public_replication/01_prepare_response_blind_inputs.R) — A034
- [02_lock_outcome_maps.R](sources/HCC_ICI_project/04_scripts/phase15b_public_replication/02_lock_outcome_maps.R) — A035
- [03_analyze_public_replication.R](sources/HCC_ICI_project/04_scripts/phase15b_public_replication/03_analyze_public_replication.R) — A036
- [01_run_phase15c_benchmark.R](sources/HCC_ICI_project/04_scripts/phase15c_public_benchmark/01_run_phase15c_benchmark.R) — A037
- [03_phase19b_response_blind_scoring.R](sources/HCC_ICI_project/04_scripts/phase19/03_phase19b_response_blind_scoring.R) — A038
- [04_phase19b_freeze_response_blind_tables.py](sources/HCC_ICI_project/04_scripts/phase19/04_phase19b_freeze_response_blind_tables.py) — A039
- [05_phase19b_outcome_analysis.R](sources/HCC_ICI_project/04_scripts/phase19/05_phase19b_outcome_analysis.R) — A040
- [01_gse104580_tace_falsification.R](sources/HCC_ICI_project/04_scripts/pretreatment_clinical_translation/01_gse104580_tace_falsification.R) — A044

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
