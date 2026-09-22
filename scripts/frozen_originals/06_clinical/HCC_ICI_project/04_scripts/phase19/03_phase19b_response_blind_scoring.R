options(stringsAsFactors = FALSE, width = 240)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(GSVA)
  library(jsonlite)
})

project_dir <- normalizePath("HCC_ICI_project", winslash = "/", mustWork = TRUE)
metadata_dir <- file.path(project_dir, "01_metadata", "phase19b_published_benchmark")
processed_dir <- file.path(project_dir, "03_processed_data", "phase19b_published_benchmark")
intermediate_dir <- file.path(project_dir, "10_intermediate_files", "phase19")
log_dir <- file.path(project_dir, "09_logs", "phase19b_published_benchmark")
invisible(lapply(c(metadata_dir, processed_dir, intermediate_dir, log_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

matrix_path <- file.path(project_dir, "01_metadata", "phase18_literature_landscape", "PHASE18B_DIRECT_AND_NEAR_COMPETITOR_MATRIX.tsv")
axis_path <- file.path(project_dir, "01_metadata", "clinical_translation", "PTPN22_CLINICAL_AXES_FROZEN.tsv")
erp_rds_path <- file.path(project_dir, "03_processed_data", "pretreatment_clinical_translation", "ERP117672_gene_expression_response_blind.rds")
g302_rds_path <- file.path(project_dir, "03_processed_data", "phase15b_public_replication", "GSE302495_RESPONSE_BLIND_MATRIX_AND_SCORES.rds")
erp_benchmark_path <- file.path(project_dir, "05_results", "tables", "phase15c_public_benchmark", "ERP117672_BENCHMARK_SCORES_RESPONSE_BLIND.tsv")
g302_benchmark_path <- file.path(project_dir, "05_results", "tables", "phase15c_public_benchmark", "GSE302495_BENCHMARK_SCORES_RESPONSE_BLIND.tsv")
s007_supp_path <- file.path(project_dir, "02_raw_data", "phase19_published_signature_benchmark", "PMC11590841_supp", "jitc-12-11-s001.pdf")
s006_mmc1_path <- file.path(project_dir, "02_raw_data", "phase19_published_signature_benchmark", "PMC12086051_supp", "mmc1.pdf")
s006_mmc2_path <- file.path(project_dir, "02_raw_data", "phase19_published_signature_benchmark", "PMC12086051_supp", "mmc2.xlsx")

expected_hashes <- c(
  matrix_path = "bc024912d66a9be3c103a735550b300b0c2272ed7badac01b8cb3aac88621aa5",
  axis_path = "268474ae8a9f6ca742e7b9197e4e82e71edd6297aef423619961ae723257927f",
  erp_rds_path = "ec4c1fd4c25d04e81167766c252da72047256a0352cd416920aba74847cc63f6",
  g302_rds_path = "b138295da8af3f7f8965dc68ff429aa10cf63e90f3e0ce286a1185e299322d1e",
  erp_benchmark_path = "cc87cc0bd3d6af3db1eef1e07684a65bd18812bb98640f76c88ad25adff79b91",
  g302_benchmark_path = "ab8a6acc323ceb6e3c9af69481941a2b7468490c4b84ec1d00b89d3eed037311",
  s007_supp_path = "e56ccc65198ad1cb39973054fff02e73d618d4a7aaf5285f67d65944c0475ce4",
  s006_mmc1_path = "2bca77592a2b6166890531c51753da7752edbbdf0442bd96481db3aa76d7782a",
  s006_mmc2_path = "e48e6be27015d95f1f01fb4fa2b6a8bb6ec9c2c632be25cc982b83d7b48c9ae2"
)
input_paths <- c(
  matrix_path = matrix_path, axis_path = axis_path, erp_rds_path = erp_rds_path,
  g302_rds_path = g302_rds_path, erp_benchmark_path = erp_benchmark_path,
  g302_benchmark_path = g302_benchmark_path, s007_supp_path = s007_supp_path,
  s006_mmc1_path = s006_mmc1_path, s006_mmc2_path = s006_mmc2_path
)
stopifnot(all(file.exists(input_paths)))
actual_hashes <- vapply(input_paths, digest, character(1L), algo = "sha256", file = TRUE)
stopifnot(identical(unname(actual_hashes), unname(expected_hashes)))

safe_z <- function(x) {
  x <- as.numeric(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x)) / s
}

matrix <- fread(matrix_path)
expected_ids <- c("S005", "S006", "S007", "S008", "S009", "S010", "S011", "S012", "S013", "S014", "S015", "S021", "S022", "S024", "P18R13623")
stopifnot(identical(matrix$paper_id, expected_ids))

audit <- data.table(
  paper_id = expected_ids,
  review_class = c(rep("PEER_REVIEWED", 5), "PREPRINT_CONTEXT_ONLY", rep("PEER_REVIEWED", 9)),
  hcc_ici_setting = c(
    "HCC pembrolizumab pretreatment tumor / longitudinal PBMC",
    "HCC atezolizumab plus bevacizumab pretreatment bulk with scRNA-derived states",
    "HCC sintilimab plus bevacizumab paired ascites with external bulk validation",
    "HCC neoadjuvant nivolumab longitudinal blood and tumor",
    "HCC anti-PD-1 plus lenvatinib longitudinal single-cell study",
    "HCC immunotherapy spatial preprint; context only",
    "HCC neoadjuvant nivolumab multimodal resistance study",
    "HCC atezolizumab plus bevacizumab primary refractoriness",
    "HCC combination immunotherapy longitudinal PBMC",
    "HCC neoadjuvant immunotherapy post-treatment TLS morphology",
    "HCC post-resection recurrence; not an ICI-response benchmark",
    "HCC immune-exclusion mechanism with preclinical ICB context",
    "HCC prognosis with murine anti-PD-1 context",
    "HCC immunotherapy-response radiomics/spatial TRM study",
    "HCC neoadjuvant anti-PD-1 niche and clone-state study"
  ),
  signature_name = c(
    "No frozen multigene contract", "21 scRNA-derived cell-state signatures", "S007_ISG_20",
    "No fixed transferable score", "No fixed transferable score", "SPARC model context only",
    "No transferable locked score", "Conditional myeloid refractoriness model", "No fixed transferable score",
    "TLS density/morphology", "TIMES five-marker recurrence model", "No fixed clinical score",
    "Ductular-reaction prognostic score", "Radiomics response model", "No fixed clinical score"
  ),
  biological_unit = c(
    "cytotoxic T-cell context", "21 single-cell states projected to bulk", "interferon-stimulated gene program",
    "top intratumoral clonotypes and blood-tumor sharing", "GZMK cTem/Tpex versus KIR CD8/Treg states",
    "spatial response/resistance neighborhoods", "TAM-exhausted CD8-lipid interface", "myeloid primary-refractory phenotype",
    "regimen-specific systemic immune states", "tertiary lymphoid structure morphology", "spatial recurrence risk",
    "DLL4-NOTCH3 capillary-mCAF immune exclusion", "cholic acid-NR1H4 CD8 dysfunction",
    "CD8 TRM spatial distribution", "mregDC-CXCL13 TH-progenitor exhausted CD8 niche"
  ),
  derivation_type = c(
    "mixed-or-unclear", "outcome-derived", "mixed-or-unclear", "outcome-derived", "outcome-derived",
    "outcome-derived", "mixed-or-unclear", "outcome-derived", "outcome-derived", "outcome-derived",
    "outcome-derived", "outcome-independent", "outcome-derived", "outcome-derived", "outcome-independent"
  ),
  exact_gene_list_source = c(
    rep("NOT_AVAILABLE_AS_FROZEN_TRANSFERABLE_SIGNATURE", 1),
    "Official Supplementary Table S2 in PMC12086051 mmc2.xlsx; separate up/down lists for 21 signatures",
    "Official Supplementary Table 4, jitc-12-11-s001.pdf page 44 (PDF page index 43)",
    rep("NOT_AVAILABLE_AS_FROZEN_TRANSFERABLE_SIGNATURE", 12)
  ),
  exact_weights_source = c(
    "NOT_AVAILABLE", "Separate marker directions exist, but continuous signed-combination weights are not specified",
    "Unweighted membership gene set", rep("NOT_AVAILABLE_OR_NOT_APPLICABLE", 12)
  ),
  exact_scoring_method_source = c(
    "NOT_AVAILABLE", "Supplementary Methods report ssGSEA/NTP, but the continuous combination of up/down lists is unspecified",
    "Main Methods: ssGSEA using GSVA 1.46.0 on log2(TPM+1); no outcome cutpoint used in Phase19",
    rep("NOT_AVAILABLE_AS_EXACT_BULK_GENE_EXPRESSION_RULE", 12)
  ),
  intended_platform = c(
    "bulk tumor and PBMC context", "bulk tumor RNA-seq", "bulk tumor RNA-seq",
    "TCR/transcriptomics", "single-cell RNA/TCR", "spatial multi-omics", "spatial/scRNA/bulk multi-omics",
    "pathology/IMC/RNA/scRNA conditional model", "PBMC single-cell multi-omics", "spatial histology",
    "spatial profiling / XGBoost", "spatial multi-omics / preclinical", "pathology/spatial/metabolomics",
    "radiomics/spatial profiling", "single-cell/TCR/spatial niche"
  ),
  transferable_to_bulk = c("NO", "NO", "YES", rep("NO", 12)),
  eligibility = c(
    "PUBLISHED_SIGNATURE_DEFINITION_NO_GO", "PUBLISHED_SIGNATURE_DEFINITION_NO_GO", "ELIGIBLE_TIER_A",
    rep("PUBLISHED_SIGNATURE_DEFINITION_NO_GO", 12)
  ),
  reason_if_no = c(
    "No exact frozen multigene gene list, weights and rule",
    "Exact up/down gene lists exist, but the exact continuous signed ssGSEA combination rule is not specified; proxy inference prohibited",
    "ELIGIBLE",
    "No fixed transferable bulk score", "No exact transferable bulk scoring contract", "Preprint cannot satisfy primary gate and model is not a simple bulk gene-expression signature",
    "No transferable locked score", "Conditional pathology/myeloid model is not an exact transferable bulk gene-expression signature",
    "No fixed transferable bulk score", "Post-treatment TLS morphology is not a bulk gene-expression signature",
    "Different endpoint and non-ICI; outcome-trained spatial XGBoost model is not transferable to these bulk cohorts",
    "Mechanistic/preclinical study without fixed clinical bulk score", "Prognostic/preclinical score without exact human ICI bulk rule",
    "Radiomics model is not applicable to bulk tumor expression", "No fixed transferable bulk score"
  )
)

candidate_pool <- merge(
  matrix[, .(paper_id, title, doi, pmid)], audit, by = "paper_id", sort = FALSE
)
candidate_pool <- candidate_pool[match(expected_ids, paper_id)]
candidate_pool[, phase18_frozen_pool := TRUE]

isg_genes <- c(
  "IFIT1", "IFIT2", "IFIT3", "IFIT5", "IFI35", "ISG15", "IRF1", "IRF3", "IRF7", "IRF9",
  "OASL", "ADAR", "EIF2AK2", "EPSTI1", "GBP4", "MX1", "OGFR", "SP110", "STAT1", "STAT2"
)
stopifnot(length(isg_genes) == 20L, !anyDuplicated(isg_genes))

contract <- data.table(
  paper_id = expected_ids,
  signature_name = audit$signature_name,
  exact_definition_status = audit$eligibility,
  gene_list = fifelse(expected_ids == "S007", paste(isg_genes, collapse = ";"), "NA"),
  requested_gene_n = fifelse(expected_ids == "S007", 20L, NA_integer_),
  weight_rule = fifelse(expected_ids == "S007", "UNWEIGHTED_MEMBERSHIP", "NOT_REPRODUCIBLE_OR_NOT_APPLICABLE"),
  scoring_rule = fifelse(expected_ids == "S007", "SSGSEA_GSVA_ALPHA_0.25_DEFAULT_GLOBAL_RANGE_NORMALIZATION", "NO_GO"),
  input_scale = fifelse(expected_ids == "S007", "ERP: supplied gene_expression (log2 TPM+1 convention); GSE302495: supplied log2 CPM+1 response-blind matrix", "NA"),
  missing_gene_rule = fifelse(expected_ids == "S007", ">=80% coverage if source gives no stricter rule; no outcome-informed substitution", "NA"),
  score_direction = fifelse(expected_ids == "S007", "HIGHER_EQUALS_GREATER_ISG_ENRICHMENT; NEVER_FLIPPED", "NA"),
  exact_definition_source = audit$exact_gene_list_source,
  reason_if_no_go = audit$reason_if_no,
  peer_review_gate_eligible = expected_ids != "S010",
  outcome_joined = FALSE
)

tier_assignment <- data.table(
  paper_id = expected_ids,
  signature_name = audit$signature_name,
  tier = fifelse(expected_ids == "S007", "A", fifelse(expected_ids == "S010", "CONTEXT_ONLY_PREPRINT", "NOT_ASSIGNED_DEFINITION_NO_GO")),
  tier_reason = fifelse(
    expected_ids == "S007",
    "Direct HCC ICI response-linked interferon gene-expression state, pretreatment-applicable and bulk-transferable",
    fifelse(expected_ids == "S010", "Preprint cannot satisfy primary gate", "Exact transferable signature definition not eligible")
  ),
  assigned_before_outcome = TRUE
)

axes <- fread(axis_path)
af_genes <- strsplit(axes[component == "PTPN22_activation_feedback_35", exact_features], ";", fixed = TRUE)[[1L]]
stopifnot(length(af_genes) == 35L, !("PTPN22" %chin% af_genes))

erp_object <- readRDS(erp_rds_path)
g302_object <- readRDS(g302_rds_path)
stopifnot(identical(erp_object$response_joined, FALSE), identical(g302_object$response_joined, FALSE))
erp_mat <- erp_object$gene_expression
g302_mat <- g302_object$gene_expression_log2_cpm_plus_1
stopifnot(is.matrix(erp_mat), is.matrix(g302_mat), ncol(erp_mat) == 40L, ncol(g302_mat) == 38L)

score_isg <- function(mat, dataset, patient_ids) {
  present <- intersect(isg_genes, rownames(mat))
  variable <- present[vapply(present, function(g) is.finite(sd(mat[g, ])) && sd(mat[g, ]) > 0, logical(1L))]
  stopifnot(length(present) >= 16L, length(variable) >= 16L)
  # GSVA 1.52 returns NaN when a single-set call requests internal normalization.
  # Run the identical raw ssGSEA walk and apply the documented default global
  # range normalization explicitly; this is response-blind and deterministic.
  param <- ssgseaParam(mat, list(S007_ISG_20 = variable), minSize = 1, maxSize = Inf, alpha = 0.25, normalize = FALSE)
  raw_scored <- gsva(param, verbose = FALSE)
  stopifnot(identical(colnames(raw_scored), colnames(mat)), nrow(raw_scored) == 1L, all(is.finite(raw_scored)))
  score_range <- max(raw_scored) - min(raw_scored)
  stopifnot(is.finite(score_range), score_range > 0)
  scored <- raw_scored / score_range
  table <- data.table(
    dataset = dataset, patient_id = patient_ids, S007_ISG_20 = as.numeric(scored[1L, ]),
    requested_gene_n = 20L, present_gene_n = length(present), variable_gene_n = length(variable),
    coverage_fraction = length(present) / 20, scoring_rule = "SSGSEA_GSVA_ALPHA_0.25_DEFAULT_GLOBAL_RANGE_NORMALIZATION",
    score_direction = "HIGHER_EQUALS_GREATER_ISG_ENRICHMENT", scoring_scope = "WITHIN_COHORT_RESPONSE_BLIND",
    outcome_joined = FALSE
  )
  list(table = table, present = present, variable = variable)
}

erp_scored <- score_isg(erp_mat, "ERP117672", colnames(erp_mat))
g302_scored <- score_isg(g302_mat, "GSE302495", g302_object$patient_ids)
cat(sprintf("S007 score finite counts ERP/GSE: %d/%d; %d/%d\n", sum(is.finite(erp_scored$table$S007_ISG_20)), nrow(erp_scored$table), sum(is.finite(g302_scored$table$S007_ISG_20)), nrow(g302_scored$table)))

erp_bench <- fread(erp_benchmark_path)
g302_bench <- fread(g302_benchmark_path)
stopifnot(all(erp_bench$outcome_joined == FALSE), all(g302_bench$outcome_joined == FALSE))

score_program <- function(mat, genes) {
  present <- intersect(genes, rownames(mat))
  variable <- present[vapply(present, function(g) is.finite(sd(mat[g, ])) && sd(mat[g, ]) > 0, logical(1L))]
  stopifnot(length(variable) >= ceiling(0.8 * length(genes)))
  colMeans(t(scale(t(mat[variable, , drop = FALSE]))))
}

erp_identity <- merge(
  data.table(patient_id = colnames(erp_mat), AF35 = safe_z(score_program(erp_mat, af_genes))),
  erp_bench[, .(patient_id, broad_immune, IFNg_response)], by = "patient_id", sort = FALSE
)
erp_identity <- merge(erp_identity, erp_scored$table[, .(patient_id, S007_ISG_20)], by = "patient_id", sort = FALSE)
g302_identity <- merge(
  data.table(patient_id = g302_object$patient_ids, AF35 = safe_z(g302_object$molecular_scores$AF35_score)),
  g302_bench[, .(patient_id, broad_immune, IFNg_response)], by = "patient_id", sort = FALSE
)
g302_identity <- merge(g302_identity, g302_scored$table[, .(patient_id, S007_ISG_20)], by = "patient_id", sort = FALSE)
stopifnot(nrow(erp_identity) == 40L, nrow(g302_identity) == 38L)
if (anyNA(erp_identity) || anyNA(g302_identity)) {
  erp_na <- names(erp_identity)[vapply(erp_identity, anyNA, logical(1L))]
  g302_na <- names(g302_identity)[vapply(g302_identity, anyNA, logical(1L))]
  stop(sprintf("Response-blind identity table contains NA; ERP=%s; GSE302495=%s", paste(erp_na, collapse = ","), paste(g302_na, collapse = ",")))
}

overlap <- intersect(af_genes, isg_genes)
relationships <- data.table(
  paper_id = "S007", signature_name = "S007_ISG_20", derivation_type = audit[paper_id == "S007", derivation_type],
  af35_gene_n = 35L, published_gene_n = 20L, overlap_gene_n = length(overlap), overlap_genes = paste(overlap, collapse = ";"),
  jaccard_index = length(overlap) / length(union(af_genes, isg_genes)),
  ERP117672_AF35_signature_spearman = cor(erp_identity$AF35, erp_identity$S007_ISG_20, method = "spearman"),
  GSE302495_AF35_signature_spearman = cor(g302_identity$AF35, g302_identity$S007_ISG_20, method = "spearman"),
  ERP117672_signature_broad_immune_spearman = cor(erp_identity$S007_ISG_20, erp_identity$broad_immune, method = "spearman"),
  GSE302495_signature_broad_immune_spearman = cor(g302_identity$S007_ISG_20, g302_identity$broad_immune, method = "spearman"),
  ERP117672_signature_frozen_IFNg_spearman = cor(erp_identity$S007_ISG_20, erp_identity$IFNg_response, method = "spearman"),
  GSE302495_signature_frozen_IFNg_spearman = cor(g302_identity$S007_ISG_20, g302_identity$IFNg_response, method = "spearman"),
  near_duplicate_gene_rule = length(overlap) / length(union(af_genes, isg_genes)) >= 0.8,
  near_duplicate_score_rule = abs(cor(erp_identity$AF35, erp_identity$S007_ISG_20, method = "spearman")) >= 0.95 && abs(cor(g302_identity$AF35, g302_identity$S007_ISG_20, method = "spearman")) >= 0.95,
  near_duplicate_predefined_flag = FALSE,
  outcome_joined = FALSE
)
relationships[, near_duplicate_predefined_flag := near_duplicate_gene_rule & near_duplicate_score_rule]

residualize_blind <- function(identity) {
  data.table(
    patient_id = identity$patient_id,
    AF35_residual_on_S007_ISG_20 = safe_z(residuals(lm(AF35 ~ S007_ISG_20, data = identity))),
    outcome_joined = FALSE
  )
}
erp_residual <- residualize_blind(erp_identity)
g302_residual <- residualize_blind(g302_identity)

tables <- list(
  "PHASE19B_PUBLISHED_SIGNATURE_CANDIDATE_POOL.tsv" = candidate_pool,
  "PHASE19B_EXACT_PUBLISHED_SIGNATURE_CONTRACT.tsv" = contract,
  "ERP117672_PUBLISHED_SIGNATURE_SCORES_RESPONSE_BLIND.tsv" = erp_scored$table,
  "GSE302495_PUBLISHED_SIGNATURE_SCORES_RESPONSE_BLIND.tsv" = g302_scored$table,
  "AF35_PUBLISHED_SIGNATURE_MOLECULAR_RELATIONSHIPS.tsv" = relationships,
  "PHASE19B_PUBLISHED_SIGNATURE_TIER_ASSIGNMENT.tsv" = tier_assignment
)

table_json_path <- file.path(intermediate_dir, "phase19b_response_blind_tables.json")
write_json(tables, table_json_path, dataframe = "rows", auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "string")

blind_object_path <- file.path(processed_dir, "PHASE19B_RESPONSE_BLIND_OBJECT.rds")
saveRDS(list(
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), outcome_joined = FALSE,
  isg_genes = isg_genes, af35_genes = af_genes, candidate_pool = candidate_pool, contract = contract,
  tier_assignment = tier_assignment, relationships = relationships,
  ERP117672 = list(scores = erp_scored$table, identity = erp_identity, residuals = erp_residual),
  GSE302495 = list(scores = g302_scored$table, identity = g302_identity, residuals = g302_residual)
), blind_object_path, version = 3)

vector_hash <- function(x) digest(writeBin(as.double(x), raw(), size = 8L, endian = "little"), algo = "sha256", serialize = FALSE)
freeze <- list(
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), outcome_joined = FALSE,
  candidate_pool_source_sha256 = unname(actual_hashes[["matrix_path"]]),
  candidate_n = nrow(candidate_pool), peer_reviewed_candidate_n = sum(candidate_pool$review_class == "PEER_REVIEWED"),
  eligible_signature_n = 1L, eligible_tier_a_signature_n = 1L, eligible_peer_reviewed_paper_n = 1L,
  eligible_signatures = "S007_ISG_20", primary_gate_prerequisite = "FAIL_LT_2_INDEPENDENT_PEER_REVIEWED_TIER_A",
  exact_gene_list = isg_genes, gene_list_sha256 = digest(paste(isg_genes, collapse = "\n"), algo = "sha256", serialize = FALSE),
  published_method = "ssGSEA via GSVA 1.46.0 on log2(TPM+1)",
  local_implementation = paste0("GSVA_", as.character(packageVersion("GSVA")), "_ssgseaParam_alpha_0.25_raw_then_explicit_default_global_range_normalization"),
  ERP117672_present_gene_n = length(erp_scored$present), GSE302495_present_gene_n = length(g302_scored$present),
  ERP117672_score_vector_sha256 = vector_hash(erp_scored$table$S007_ISG_20),
  GSE302495_score_vector_sha256 = vector_hash(g302_scored$table$S007_ISG_20),
  ERP117672_residual_vector_sha256 = vector_hash(erp_residual$AF35_residual_on_S007_ISG_20),
  GSE302495_residual_vector_sha256 = vector_hash(g302_residual$AF35_residual_on_S007_ISG_20),
  blind_object_sha256 = digest(blind_object_path, algo = "sha256", file = TRUE),
  response_map_read = FALSE, response_label_read = FALSE, no_cutpoint = TRUE, no_proxy = TRUE, score_direction_locked = TRUE,
  input_hashes = as.list(actual_hashes)
)
freeze_path <- file.path(metadata_dir, "PHASE19B_RESPONSE_BLIND_FREEZE.json")
write_json(freeze, freeze_path, auto_unbox = TRUE, pretty = TRUE, digits = 17)

writeLines(c(
  "phase=PHASE19B_RESPONSE_BLIND_SCORING",
  "status=RESPONSE_BLIND_FREEZE_COMPLETE",
  "outcome_joined=FALSE",
  "response_map_read=FALSE",
  "eligible_signature=S007_ISG_20",
  "eligible_peer_reviewed_tier_A_papers=1",
  "primary_gate_prerequisite=FAIL_LT_2_INDEPENDENT_PEER_REVIEWED_TIER_A",
  sprintf("blind_object_sha256=%s", freeze$blind_object_sha256),
  sprintf("response_blind_tables_json_sha256=%s", digest(table_json_path, algo = "sha256", file = TRUE)),
  "new_literature_horizon=FALSE",
  "new_signature=FALSE",
  "outcome_cutpoint=FALSE"
), file.path(log_dir, "PHASE19B_RESPONSE_BLIND_SCORING_LOG.txt"), useBytes = TRUE)

cat("Phase 19B response-blind scoring complete\n")
cat(sprintf("Eligible Tier A signatures: %d\n", freeze$eligible_tier_a_signature_n))
cat(sprintf("S007 overlap/Jaccard: %d / %.6f\n", relationships$overlap_gene_n, relationships$jaccard_index))
cat(sprintf("AF35 rho (ERP/GSE): %.6f / %.6f\n", relationships$ERP117672_AF35_signature_spearman, relationships$GSE302495_AF35_signature_spearman))
cat(sprintf("Blind object SHA256: %s\n", freeze$blind_object_sha256))
