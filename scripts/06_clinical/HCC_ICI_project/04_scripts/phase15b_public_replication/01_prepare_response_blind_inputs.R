options(stringsAsFactors = FALSE, width = 240)
invisible(try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE))

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

workspace_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
project_dir <- file.path(workspace_dir, "HCC_ICI_project")
raw_dir <- file.path(project_dir, "02_raw_data", "phase15a_horizon_access")
processed_dir <- file.path(project_dir, "03_processed_data", "phase15b_public_replication")
table_dir <- file.path(project_dir, "05_results", "tables", "phase15b_public_replication")
result_dir <- file.path(project_dir, "05_results", "phase15b_public_replication")
dir.create(processed_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

gse302_path <- file.path(raw_dir, "GSE302495_GEO_readcount_tbl.txt.gz")
gse287_path <- file.path(raw_dir, "GSE287319_Dep_n12_readcount.gct.gz")
map_path <- file.path(project_dir, "01_metadata", "pretreatment_clinical_translation", "ENSEMBL98_GENE_ID_SYMBOL_MAP.tsv")
af_path <- file.path(project_dir, "01_metadata", "spatial_phase8", "PTPN22_SPATIAL_PROGRAM_FROZEN_PROVENANCE.tsv")
comparator_path <- file.path(project_dir, "01_metadata", "PHASE3_PTPN22_SPECIFICITY_COMPARATOR_PROGRAMS.tsv")
stopifnot(all(file.exists(c(gse302_path, gse287_path, map_path, af_path, comparator_path))))

expected_hashes <- c(
  GSE302495 = "2b9a88f0f441772ffaf9561fff73f0d0ac28f71bbcbfffc6a6cad916bea8bdbd",
  GSE287319 = "db1371ce97c4023f99c256e521bf1c595427209c3fb8ad074c46566fd15ed112",
  ENSEMBL98 = "4def1e77a8c184a144d874a04bd9a8273792c7fa839a48441799b0d563b66150",
  AF35 = "e91558ff99f57a975c9114f980a265e17096482b008705b4fdb9f4f4ccec154c",
  COMPARATORS = "aa8e303efca69ea76e57e5c2c90935da04ce50d5e7555a6ca99f4e707e014434"
)
actual_hashes <- c(
  GSE302495 = digest(gse302_path, algo = "sha256", file = TRUE),
  GSE287319 = digest(gse287_path, algo = "sha256", file = TRUE),
  ENSEMBL98 = digest(map_path, algo = "sha256", file = TRUE),
  AF35 = digest(af_path, algo = "sha256", file = TRUE),
  COMPARATORS = digest(comparator_path, algo = "sha256", file = TRUE)
)
stopifnot(identical(unname(actual_hashes), unname(expected_hashes)))

af_genes <- fread(af_path)$gene
comparators <- fread(comparator_path)
broad_genes <- comparators[program == "broad_immune_abundance", gene]
stopifnot(length(af_genes) == 35L, !anyDuplicated(af_genes), length(broad_genes) == 11L, !anyDuplicated(broad_genes))

safe_z <- function(x) {
  x <- as.numeric(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x)) / s
}

score_program <- function(mat, genes, min_present, min_variable_fraction) {
  present <- intersect(genes, rownames(mat))
  variable <- present[vapply(present, function(g) is.finite(sd(mat[g, ])) && sd(mat[g, ]) > 0, logical(1L))]
  eligible <- length(present) >= min_present && length(variable) >= ceiling(min_variable_fraction * length(present))
  score <- rep(NA_real_, ncol(mat))
  if (eligible) score <- safe_z(colMeans(t(scale(t(mat[variable, , drop = FALSE])))))
  list(
    score = score, requested_n = length(genes), present_n = length(present), variable_n = length(variable),
    present = present, missing = setdiff(genes, present), eligible = eligible
  )
}

aggregate_ensembl_counts <- function(dt, id_col, sample_cols, mapping) {
  long_ids <- sub("\\..*$", "", dt[[id_col]])
  keep_map <- mapping[match(long_ids, mapping$gene_id), ]
  keep <- !is.na(keep_map$gene_symbol) & nzchar(keep_map$gene_symbol)
  count_dt <- copy(dt[keep, c(id_col, sample_cols), with = FALSE])
  count_dt[, gene_symbol := keep_map$gene_symbol[keep]]
  count_dt[, (id_col) := NULL]
  aggregate_dt <- count_dt[, lapply(.SD, sum), by = gene_symbol, .SDcols = sample_cols]
  mat <- as.matrix(aggregate_dt[, ..sample_cols])
  storage.mode(mat) <- "double"
  rownames(mat) <- aggregate_dt$gene_symbol
  mat
}

to_logcpm <- function(counts) {
  library_sizes <- colSums(counts)
  stopifnot(all(is.finite(library_sizes)), all(library_sizes > 0))
  log2(t(t(counts) / library_sizes * 1e6) + 1)
}

write_coverage <- function(dataset, af, broad, ptpn22_variable, path) {
  rows <- data.table(
    dataset = dataset,
    program = c("PTPN22", "AF35", "broad_immune_abundance"),
    requested_n = c(1L, af$requested_n, broad$requested_n),
    present_n = c(1L, af$present_n, broad$present_n),
    variable_n = c(as.integer(ptpn22_variable), af$variable_n, broad$variable_n),
    missing_genes = c("NONE", if (length(af$missing)) paste(af$missing, collapse = ";") else "NONE", if (length(broad$missing)) paste(broad$missing, collapse = ";") else "NONE"),
    eligibility = c(if (ptpn22_variable) "PASS" else "FAIL", if (af$eligible) "PASS" else "FAIL", if (broad$eligible) "PASS" else "FAIL"),
    outcome_loaded = "NO"
  )
  fwrite(rows, path, sep = "\t", quote = TRUE)
}

mapping <- fread(map_path, select = c("gene_id", "gene_symbol"))
mapping[, gene_id := sub("\\..*$", "", gene_id)]
mapping <- unique(mapping[nzchar(gene_symbol)], by = "gene_id")

# GSE302495: no outcome file is opened above or in this block.
g302 <- fread(gse302_path)
stopifnot(names(g302)[1L] == "GeneID", ncol(g302) == 71L)
g302_samples <- names(g302)[-1L]
g302_pre <- g302_samples[grepl("_PRE$", g302_samples)]
g302_post <- setdiff(g302_samples, g302_pre)
stopifnot(length(g302_pre) == 38L, length(g302_post) == 32L, !anyDuplicated(g302_pre))
g302_patient_ids <- gsub("_", "-", sub("_PRE$", "", g302_pre))
stopifnot(length(unique(g302_patient_ids)) == 38L)
g302_counts <- aggregate_ensembl_counts(g302, "GeneID", g302_pre, mapping)
g302_logcpm <- to_logcpm(g302_counts)
colnames(g302_logcpm) <- g302_patient_ids
colnames(g302_counts) <- g302_patient_ids
stopifnot("PTPN22" %in% rownames(g302_logcpm), is.finite(sd(g302_logcpm["PTPN22", ])), sd(g302_logcpm["PTPN22", ]) > 0)
g302_af <- score_program(g302_logcpm, af_genes, 28L, 0.8)
g302_broad <- score_program(g302_logcpm, broad_genes, 11L, 0.8)
stopifnot(g302_af$eligible, g302_broad$eligible)
g302_ptpn22 <- as.numeric(g302_logcpm["PTPN22", ])
g302_ptpn22_resid <- safe_z(residuals(lm(g302_ptpn22 ~ g302_broad$score)))
g302_af_resid <- safe_z(residuals(lm(g302_af$score ~ g302_broad$score)))
g302_scores <- data.table(
  patient_id = g302_patient_ids,
  PTPN22_expression = g302_ptpn22,
  AF35_score = g302_af$score,
  broad_immune_score = g302_broad$score,
  PTPN22_broad_immune_residual = g302_ptpn22_resid,
  AF35_broad_immune_residual = g302_af_resid
)
stopifnot(!anyNA(g302_scores), !anyDuplicated(g302_scores$patient_id))
g302_object <- list(
  dataset = "GSE302495", source_file = basename(gse302_path), source_sha256 = actual_hashes[["GSE302495"]],
  source_dimensions = c(genes = nrow(g302), samples = length(g302_samples)), baseline_patient_n = 38L,
  excluded_postoperation_sample_n = 32L, patient_ids = g302_patient_ids,
  gene_counts = g302_counts, gene_expression_log2_cpm_plus_1 = g302_logcpm,
  molecular_scores = g302_scores, response_joined = FALSE,
  processing = "Ensembl version stripped; frozen Ensembl98 mapping; duplicate symbol counts summed; log2(CPM+1)",
  af35_rule = ">=28/35 present and >=80% of present genes variable"
)
g302_rds <- file.path(processed_dir, "GSE302495_RESPONSE_BLIND_MATRIX_AND_SCORES.rds")
saveRDS(g302_object, g302_rds, version = 3)
g302_rds_hash <- digest(g302_rds, algo = "sha256", file = TRUE)
fwrite(g302_scores, file.path(table_dir, "GSE302495_PATIENT_LEVEL_SCORES.tsv"), sep = "\t", quote = TRUE)
write_coverage("GSE302495", g302_af, g302_broad, TRUE, file.path(processed_dir, "GSE302495_RESPONSE_BLIND_GENE_COVERAGE.tsv"))

g302_audit <- c(
  "# GSE302495 response-blind input audit", "",
  "Status: `PASS_RESPONSE_BLIND_INPUT_LOCK`  ",
  "Outcome loaded during construction: `NO`", "",
  "## Source and dimensions", "",
  sprintf("- Source: `%s`", basename(gse302_path)),
  sprintf("- Source SHA-256: `%s`", actual_hashes[["GSE302495"]]),
  sprintf("- Source dimensions: `%d` gene-ID rows × `%d` sample columns.", nrow(g302), length(g302_samples)),
  sprintf("- Pretreatment columns locked: `%d` unique patients.", length(g302_pre)),
  sprintf("- Operation/post-treatment columns excluded before outcome join: `%d`.", length(g302_post)),
  sprintf("- Ensembl98-mapped unique symbols after duplicate-symbol summation: `%d`.", nrow(g302_logcpm)), "",
  "## Frozen transformation", "",
  "Ensembl version suffixes were stripped, frozen Ensembl98 IDs were mapped to symbols, duplicate IDs mapping to one symbol were summed, and `log2(CPM+1)` was calculated within each library. No cross-cohort normalization, ComBat, cutpoint, feature selection or imputation was used.", "",
  "## Coverage lock", "",
  sprintf("- PTPN22: present and variable (`SD=%.6f`).", sd(g302_ptpn22)),
  sprintf("- AF35: `%d/35` present, `%d` variable; eligibility `%s`.", g302_af$present_n, g302_af$variable_n, ifelse(g302_af$eligible, "PASS", "FAIL")),
  sprintf("- Broad immune: `%d/11` present, `%d` variable; eligibility `%s`.", g302_broad$present_n, g302_broad$variable_n, ifelse(g302_broad$eligible, "PASS", "FAIL")),
  "- PTPN22, AF35, broad immune and both broad-immune residuals were constructed before outcome loading.", "",
  "## Immutable response-blind object", "",
  sprintf("- File: `%s`", basename(g302_rds)),
  sprintf("- SHA-256: `%s`", g302_rds_hash),
  "- Embedded flag: `response_joined=FALSE`.", "",
  "Decision: `GSE302495_RESPONSE_BLIND_INPUT_PASS`."
)
writeLines(g302_audit, file.path(result_dir, "GSE302495_RESPONSE_BLIND_INPUT_AUDIT.md"), useBytes = TRUE)

# GSE287319: no outcome file is opened above or in this block.
g287 <- fread(gse287_path, skip = 2L)
stopifnot(names(g287)[1L] == "geneName", names(g287)[2L] == "geneId", ncol(g287) == 14L)
g287_samples <- names(g287)[-(1:2)]
g287_pre <- g287_samples[grepl("pre$", g287_samples)]
g287_post <- setdiff(g287_samples, g287_pre)
stopifnot(length(g287_pre) == 9L, length(g287_post) == 3L, !anyDuplicated(g287_pre))
g287_patient_ids <- paste0("DEPARTURE_", sub("^Dep([0-9]+)pre$", "\\1", g287_pre))
stopifnot(length(unique(g287_patient_ids)) == 9L)
g287_dt <- copy(g287[, c("geneName", g287_pre), with = FALSE])
g287_dt <- g287_dt[!is.na(geneName) & nzchar(geneName)]
g287_agg <- g287_dt[, lapply(.SD, sum), by = geneName, .SDcols = g287_pre]
g287_counts <- as.matrix(g287_agg[, ..g287_pre])
storage.mode(g287_counts) <- "double"
rownames(g287_counts) <- g287_agg$geneName
colnames(g287_counts) <- g287_patient_ids
g287_logcpm <- to_logcpm(g287_counts)
stopifnot("PTPN22" %in% rownames(g287_logcpm), is.finite(sd(g287_logcpm["PTPN22", ])), sd(g287_logcpm["PTPN22", ]) > 0)
g287_af <- score_program(g287_logcpm, af_genes, 28L, 0.8)
g287_broad <- score_program(g287_logcpm, broad_genes, 11L, 0.8)
stopifnot(g287_af$eligible, g287_broad$eligible)
g287_ptpn22 <- as.numeric(g287_logcpm["PTPN22", ])
g287_scores <- data.table(
  patient_id = g287_patient_ids,
  PTPN22_expression = g287_ptpn22,
  AF35_score = g287_af$score,
  broad_immune_score = g287_broad$score,
  PTPN22_broad_immune_residual = safe_z(residuals(lm(g287_ptpn22 ~ g287_broad$score))),
  AF35_broad_immune_residual = safe_z(residuals(lm(g287_af$score ~ g287_broad$score)))
)
stopifnot(!anyNA(g287_scores), !anyDuplicated(g287_scores$patient_id))
g287_object <- list(
  dataset = "GSE287319", source_file = basename(gse287_path), source_sha256 = actual_hashes[["GSE287319"]],
  source_dimensions = c(genes = nrow(g287), samples = length(g287_samples)), baseline_patient_n = 9L,
  excluded_post_sample_n = 3L, patient_ids = g287_patient_ids,
  gene_counts = g287_counts, gene_expression_log2_cpm_plus_1 = g287_logcpm,
  molecular_scores = g287_scores, response_joined = FALSE,
  processing = "Author gene symbols; duplicate symbol counts summed; log2(CPM+1)",
  af35_rule = ">=28/35 present and >=80% of present genes variable"
)
g287_rds <- file.path(processed_dir, "GSE287319_RESPONSE_BLIND_MATRIX_AND_SCORES.rds")
saveRDS(g287_object, g287_rds, version = 3)
g287_rds_hash <- digest(g287_rds, algo = "sha256", file = TRUE)
fwrite(g287_scores, file.path(processed_dir, "GSE287319_RESPONSE_BLIND_PATIENT_SCORES.tsv"), sep = "\t", quote = TRUE)
write_coverage("GSE287319", g287_af, g287_broad, TRUE, file.path(processed_dir, "GSE287319_RESPONSE_BLIND_GENE_COVERAGE.tsv"))

g287_audit <- c(
  "# GSE287319 response-blind input audit", "",
  "Status: `PASS_RESPONSE_BLIND_INPUT_LOCK`  ",
  "Outcome loaded during construction: `NO`", "",
  sprintf("- Source SHA-256: `%s`", actual_hashes[["GSE287319"]]),
  sprintf("- Source dimensions: `%d` gene rows × `%d` sample columns.", nrow(g287), length(g287_samples)),
  sprintf("- Pretreatment profiles locked: `%d`; post-treatment profiles excluded: `%d`.", length(g287_pre), length(g287_post)),
  sprintf("- Unique symbols after duplicate-symbol summation: `%d`.", nrow(g287_logcpm)),
  sprintf("- PTPN22 present/variable; AF35 `%d/35` present and `%d` variable; broad immune `%d/11` present and `%d` variable.", g287_af$present_n, g287_af$variable_n, g287_broad$present_n, g287_broad$variable_n),
  "- All molecular variables and response-blind residuals were constructed before endpoint loading.",
  sprintf("- Immutable object SHA-256: `%s`", g287_rds_hash),
  "- Embedded flag: `response_joined=FALSE`.", "",
  "Decision: `GSE287319_RESPONSE_BLIND_INPUT_PASS`."
)
writeLines(g287_audit, file.path(result_dir, "GSE287319_RESPONSE_BLIND_INPUT_AUDIT.md"), useBytes = TRUE)

cat("GSE302495_SOURCE_ROWS=", nrow(g302), "\n", sep = "")
cat("GSE302495_BASELINE_PATIENTS=", length(g302_pre), "\n", sep = "")
cat("GSE302495_SYMBOLS=", nrow(g302_logcpm), "\n", sep = "")
cat("GSE302495_AF35=", g302_af$present_n, "/", g302_af$variable_n, "\n", sep = "")
cat("GSE302495_BROAD=", g302_broad$present_n, "/", g302_broad$variable_n, "\n", sep = "")
cat("GSE302495_RDS_SHA256=", g302_rds_hash, "\n", sep = "")
cat("GSE287319_SOURCE_ROWS=", nrow(g287), "\n", sep = "")
cat("GSE287319_BASELINE_PATIENTS=", length(g287_pre), "\n", sep = "")
cat("GSE287319_AF35=", g287_af$present_n, "/", g287_af$variable_n, "\n", sep = "")
cat("GSE287319_BROAD=", g287_broad$present_n, "/", g287_broad$variable_n, "\n", sep = "")
cat("GSE287319_RDS_SHA256=", g287_rds_hash, "\n", sep = "")
cat("OUTCOME_LOADED=NO\n")
