ANALYSIS_SEED <- 20260828L
CANDIDATES <- c("MAP4K1", "PTPN22", "PIK3CG")

GENE_SETS <- list(
  T_lineage = c("CD3D", "CD3E", "TRAC"),
  NK = c("NKG7", "KLRD1", "GNLY", "FCGR3A"),
  T_naive_progenitor = c("TCF7", "CCR7", "LEF1", "IL7R", "SLAMF6"),
  T_effector = c("GZMK", "CCL5", "GZMB", "PRF1", "IFNG", "TNF", "FGFBP2"),
  T_terminal_exhaustion = c(
    "PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "ENTPD1", "CXCL13", "LAYN"
  ),
  T_proliferating = c("MKI67", "TOP2A", "STMN1"),
  Tfh_like = c("CXCL13", "IL21", "CH25H", "PDCD1", "ICOS", "TOX2"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2"),
  Myeloid = c("LST1", "TYROBP", "FCER1G", "CTSS", "CSF1R", "CD68", "APOE", "C1QC", "C1QB"),
  Myeloid_suppressive_TAM = c(
    "SPP1", "APOE", "C1QA", "C1QB", "C1QC", "TREM2", "LGALS3", "CTSD", "CTSB", "CD163", "MSR1"
  ),
  Myeloid_inflammatory = c("IL1B", "S100A8", "S100A9", "FCN1", "LILRB1", "CTSS"),
  B_cell = c("CD79A", "MS4A1", "CD37", "CD74"),
  Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "ALB", "APOA1")
)

TARGET_WEIGHTS <- c(
  tissue_state = 20,
  cross_scrna_replication = 20,
  ici_response = 25,
  sctcr_clonal_fate = 20,
  perturbation_plausibility = 10,
  druggability = 5
)

stopifnot(sum(TARGET_WEIGHTS) == 100)
