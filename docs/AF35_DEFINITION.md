# AF35 definition

The authoritative repository vector is config/AF35_genes.txt: 35 unique genes in fixed order, with unsigned equal weights, excluding PTPN22. Membership was not changed during repository preparation. No clinical cutoff exists and no classifier was trained.

Platform-specific implementations are recorded in config/analysis_parameters/AF35_PLATFORM_SPECIFICATION.tsv. Single-cell scores use an equal mean of gene-wise standardized log1p(CP10K) expression. Bulk inputs use their recorded assay transformations and within-cohort standardization. Spatial scores use the fixed transformed-gene mean followed by section-specific standardization and the specified technical/immune conditioning. Scores are not numerically interchangeable across platforms.

Membership, weighting, missing-gene rules and cohort-specific reference populations must not be inferred from downstream outcomes. The broader program architecture is not a recovered candidate-gene universe. No missing historical reduction procedure is reconstructed here.
