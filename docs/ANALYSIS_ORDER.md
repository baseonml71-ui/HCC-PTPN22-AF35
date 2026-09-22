# Analysis order

1. Biological anchor selection using the fixed candidate framework.
2. Patient-level state architecture and PTPN22 association.
3. The fixed AF35 definition and response-blind measurement checks.
4. Exact-linked clonotype contexts, state adjustment and longitudinal/occupancy decomposition.
5. Separate GSE238264 specification-context and Wu extension analyses.
6. Cohort-specific continuous clinical associations and sensitivity analyses, after the relevant molecular locks.
7. Orthogonal perturbation and genetic-context analyses.

This conceptual sequence is not an executable pipeline. The entrypoint map is provenance/CURRENT_ANALYSIS_ENTRYPOINTS.tsv; its source references and classification corrections are recorded alongside it. Primary script directories contain ACTIVE_CURRENT sources and required helpers only. DISPLAY_ONLY and HISTORICAL entries were excluded, not presented as current analyses.

Some source scripts contain older exploratory branches, including GSE287319 in the shared clinical preparation/analysis scripts. Those branches are not current independent manuscript validation; they remain identified as mixed-source dependencies pending a bounded portability refactor. Do not run the archive indiscriminately or bypass response locks.
