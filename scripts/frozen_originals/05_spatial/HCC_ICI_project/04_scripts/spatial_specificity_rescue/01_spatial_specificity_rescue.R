options(stringsAsFactors = FALSE, width = 260)
invisible(try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(FNN)
  library(rhdf5)
  library(arrow)
  library(digest)
})

set.seed(20260829L)

project_dir <- normalizePath("HCC_ICI_project", winslash = "/", mustWork = TRUE)
metadata_dir <- file.path(project_dir, "01_metadata", "spatial_specificity_rescue")
result_dir <- file.path(project_dir, "05_results", "spatial_specificity_rescue")
table_dir <- file.path(project_dir, "05_results", "tables", "spatial_specificity_rescue")
log_dir <- file.path(project_dir, "09_logs", "spatial_specificity_rescue")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, "01_spatial_specificity_rescue.log")
log_con <- file(log_file, open = "wt", encoding = "UTF-8")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("SPATIAL SPECIFICITY RESCUE\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Seed: 20260829\n")

assert_hash <- function(path, expected) {
  stopifnot(file.exists(path))
  observed <- digest(file = path, algo = "sha256")
  if (!identical(observed, expected)) stop("Hash mismatch: ", path, " observed=", observed)
  cat("HASH_PASS ", observed, " ", path, "\n", sep = "")
}

contract_file <- file.path(metadata_dir, "PTPN22_CONDITIONAL_SPATIAL_SCORE_CONTRACT.md")
input_audit_file <- file.path(result_dir, "SPATIAL_RESCUE_INPUT_AUDIT.md")
axis_file <- file.path(project_dir, "01_metadata", "clinical_translation", "PTPN22_CLINICAL_AXES_FROZEN.tsv")
comparator_file <- file.path(project_dir, "01_metadata", "PHASE3_PTPN22_SPECIFICITY_COMPARATOR_PROGRAMS.tsv")
program_provenance_file <- file.path(project_dir, "01_metadata", "spatial_phase8", "PTPN22_SPATIAL_PROGRAM_FROZEN_PROVENANCE.tsv")
wu_map_file <- file.path(project_dir, "01_metadata", "spatial_phase8", "WU_HCC_SPATIAL_PATIENT_SAMPLE_MAP.tsv")
gse_archive <- file.path(project_dir, "02_raw_data", "GSE238264", "GSE238264_RAW.tar")
gse_dir <- file.path(project_dir, "02_raw_data", "GSE238264", "extracted")
wu_h5 <- file.path(project_dir, "03_processed_data", "spatial_phase8", "Visium-ST", "visium_all.h5ad")
spot_summary_file <- file.path(project_dir, "03_processed_data", "spatial_phase8", "notebook_and_resources_v2", "Supporting_data", "spot_summary.pqt")
release_metadata_file <- file.path(project_dir, "03_processed_data", "spatial_phase8", "notebook_and_resources_v2", "metadata.csv")
codex_raw_dir <- file.path(project_dir, "03_processed_data", "spatial_phase8", "CODEX-MIF", "raw_data")
v10_patient_file <- file.path(project_dir, "05_results", "tables", "final_harmonization", "PTPN22_SPATIAL_NULL_PATIENT_METRICS.tsv")
v10_matching_file <- file.path(project_dir, "05_results", "tables", "final_harmonization", "PTPN22_SPATIAL_NULL_MATCHING_AUDIT.tsv")
v10_wu_audit_file <- file.path(project_dir, "05_results", "tables", "final_harmonization", "WU_VISIUM_V8_RECOMPUTATION_AUDIT.tsv")
v8_coreg_audit_file <- file.path(project_dir, "05_results", "tables", "spatial_phase8", "WU_SPATIAL_COREGISTRATION_REGION_AUDIT.tsv")

hash_contract <- "159309ab8f422bf8af0c59947570c553bae4b6a825670b9420b028388f50107a"
assert_hash(contract_file, hash_contract)
assert_hash(input_audit_file, "11eb790df374dba1a9a6d91814676829b1bf8ebfb6777c87132424275387728d")
assert_hash(axis_file, "268474ae8a9f6ca742e7b9197e4e82e71edd6297aef423619961ae723257927f")
assert_hash(comparator_file, "aa8e303efca69ea76e57e5c2c90935da04ce50d5e7555a6ca99f4e707e014434")
assert_hash(program_provenance_file, "e91558ff99f57a975c9114f980a265e17096482b008705b4fdb9f4f4ccec154c")
assert_hash(wu_map_file, "33a2b4cf4fefd182b274ece26c39757f834bf6ef177fd13cf8e40cae5890aa2a")
assert_hash(gse_archive, "cb64aad38f3f9d9441d93e14f98cc110ceea4b4f1fc4e8c78481d37444d2367f")
assert_hash(wu_h5, "b72efaa80712db6980fc3683244e56f6d7aa4f2f320139da5012cd7b8cf376e6")
assert_hash(spot_summary_file, "c2358e40f4f85af9370178ad840f2b1b7959a5fbbc46f5eb649f49a930dc0b78")
assert_hash(release_metadata_file, "9e87b7c374db5b06de3bc22e2a1c295232cf61e964c2b0548cf1c9fc3121b882")
assert_hash(v10_patient_file, "4b012c4425eebd347c35b7e397c027f00ded75c434e731154404c20c1f8df3eb")
assert_hash(v10_matching_file, "f0a11c9e097015a64f7b99dc271ca82b2beb6cb737535443f4554a559c9a2290")
assert_hash(v10_wu_audit_file, "d207ffef82e3958da612cc44f784162c26c030152473b2c5c361c23922141b91")
assert_hash(v8_coreg_audit_file, "acb2e9ed4e697b3555fd30a0e5fa6453f1d7b3c3f1c09d0f8f808df0482bb3df")

axes <- fread(axis_file)
signature <- strsplit(axes[component == "PTPN22_activation_feedback_35", exact_features], ";", fixed = TRUE)[[1L]]
comparators <- fread(comparator_file)
broad_immune <- comparators[program == "broad_immune_abundance" & !startsWith(gene, "DERIVED:"), gene]
stopifnot(length(signature) == 35L, !"PTPN22" %chin% signature, length(broad_immune) == 11L)

n_random <- 500L
n_nearest_candidates <- 200L
n_boot <- 20000L

z_safe <- function(x) {
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  as.numeric((x - mean(x)) / s)
}

conditional_residuals <- function(score_matrix, broad, log_total, log_detected) {
  y <- apply(as.matrix(score_matrix), 2L, z_safe)
  if (is.null(dim(y))) y <- matrix(y, ncol = 1L)
  covariates <- cbind(
    broad_immune = z_safe(broad),
    log_total_counts = z_safe(log_total),
    log_detected_genes = z_safe(log_detected)
  )
  keep <- apply(covariates, 2L, sd) > 0
  x <- cbind(intercept = 1, covariates[, keep, drop = FALSE])
  residual <- qr.resid(qr(x), y)
  residual <- apply(residual, 2L, z_safe)
  if (is.null(dim(residual))) residual <- matrix(residual, ncol = 1L)
  list(residual = residual, covariates = covariates, kept_covariates = colnames(covariates)[keep])
}

moran_matrix <- function(scores, nn_index) {
  scores <- as.matrix(scores)
  centered <- sweep(scores, 2L, colMeans(scores), "-")
  denom <- colSums(centered^2)
  neighbor_sum <- matrix(0, nrow = nrow(centered), ncol = ncol(centered))
  for (j in seq_len(ncol(nn_index))) neighbor_sum <- neighbor_sum + centered[nn_index[, j], , drop = FALSE]
  numerator <- colSums(centered * (neighbor_sum / ncol(nn_index)))
  out <- numerator / denom
  out[!is.finite(out)] <- NA_real_
  out
}

csr_library_sizes <- function(values, row_index, n_rows) {
  out <- numeric(n_rows)
  if (length(values)) {
    grouped <- rowsum(matrix(values, ncol = 1L), group = row_index, reorder = FALSE)
    out[as.integer(rownames(grouped))] <- grouped[, 1L]
    stopifnot(all(out[row_index] > 0))
  }
  out
}

make_random_sets <- function(gene_stats, target_genes, excluded_genes, cohort_label) {
  stopifnot(all(target_genes %chin% gene_stats$gene))
  gene_stats <- copy(gene_stats)
  gene_stats[, log_mean := log10(mean_log1p_cp10k + 1e-8)]
  clipped_detection <- pmin(pmax(gene_stats$detection_fraction, 1e-6), 1 - 1e-6)
  gene_stats[, logit_detection := qlogis(clipped_detection)]
  gene_stats[, mean_z := as.numeric(scale(log_mean))]
  gene_stats[, detection_z := as.numeric(scale(logit_detection))]
  candidate <- gene_stats[
    is.finite(mean_z) & is.finite(detection_z) & detection_fraction > 0 &
      !gene %chin% excluded_genes & !grepl("^(MT-|RPL|RPS)", gene)
  ]
  targets <- gene_stats[match(target_genes, gene)]
  nearest <- lapply(seq_len(nrow(targets)), function(i) {
    distance <- (candidate$mean_z - targets$mean_z[i])^2 +
      (candidate$detection_z - targets$detection_z[i])^2
    candidate$gene[head(order(distance), n_nearest_candidates)]
  })
  sets <- matrix(NA_character_, nrow = length(target_genes), ncol = n_random)
  for (b in seq_len(n_random)) {
    used <- character()
    for (i in sample(seq_along(target_genes))) {
      available <- setdiff(nearest[[i]], used)
      if (!length(available)) stop("No unused matched candidate for ", cohort_label, "/", target_genes[i])
      chosen <- sample(available, 1L)
      sets[i, b] <- chosen
      used <- c(used, chosen)
    }
  }
  stopifnot(all(apply(sets, 2L, function(z) length(unique(z)) == length(z))))
  matched_rows <- gene_stats[match(as.vector(sets), gene)]
  matched_array_mean <- matrix(matched_rows$mean_log1p_cp10k, nrow = length(target_genes), ncol = n_random)
  matched_array_detection <- matrix(matched_rows$detection_fraction, nrow = length(target_genes), ncol = n_random)
  audit <- data.table(
    cohort = cohort_label,
    target_gene = target_genes,
    target_mean_log1p_cp10k = targets$mean_log1p_cp10k,
    matched_mean_log1p_cp10k_mean = rowMeans(matched_array_mean),
    target_detection_fraction = targets$detection_fraction,
    matched_detection_fraction_mean = rowMeans(matched_array_detection),
    n_random_sets = n_random,
    nearest_candidate_pool_per_target = n_nearest_candidates,
    matching_is_response_blinded = TRUE
  )
  list(sets = sets, audit = audit)
}

make_membership <- function(genes, random_sets, broad_genes) {
  set_names <- c("TARGET_35", "BROAD_IMMUNE", sprintf("RANDOM_%03d", seq_len(ncol(random_sets))))
  set_genes <- c(list(signature, broad_genes), lapply(seq_len(ncol(random_sets)), function(j) random_sets[, j]))
  i <- unlist(lapply(set_genes, function(gs) match(gs, genes)), use.names = FALSE)
  if (anyNA(i)) stop("Membership gene absent")
  j <- rep(seq_along(set_genes), lengths(set_genes))
  weights <- unlist(lapply(set_genes, function(gs) rep(1 / length(gs), length(gs))), use.names = FALSE)
  sparseMatrix(i = i, j = j, x = weights, dims = c(length(genes), length(set_genes)), dimnames = list(genes, set_names))
}

empirical_two_sided <- function(observed, null) {
  centered_observed <- observed - mean(null)
  (1 + sum(abs(null - mean(null)) >= abs(centered_observed))) / (length(null) + 1)
}

group_effect <- function(values, response) {
  keep <- is.finite(values) & response %chin% c("Responder", "Nonresponder")
  values <- values[keep]
  response <- response[keep]
  r <- values[response == "Responder"]
  nr <- values[response == "Nonresponder"]
  stopifnot(length(r) > 0L, length(nr) > 0L)
  observed <- mean(r) - mean(nr)
  boot <- replicate(n_boot, mean(sample(r, length(r), replace = TRUE)) - mean(sample(nr, length(nr), replace = TRUE)))
  combined <- c(r, nr)
  allocations <- combn(seq_along(combined), length(r))
  perm <- apply(allocations, 2L, function(ind) mean(combined[ind]) - mean(combined[-ind]))
  loo <- vapply(seq_along(values), function(j) {
    v <- values[-j]
    g <- response[-j]
    if (!any(g == "Responder") || !any(g == "Nonresponder")) return(NA_real_)
    mean(v[g == "Responder"]) - mean(v[g == "Nonresponder"])
  }, numeric(1L))
  data.table(
    n_responder = length(r), n_nonresponder = length(nr), responder_mean = mean(r), nonresponder_mean = mean(nr),
    effect = observed, bootstrap_ci_low = unname(quantile(boot, 0.025)), bootstrap_ci_high = unname(quantile(boot, 0.975)),
    exact_permutation_p_two_sided = mean(abs(perm) >= abs(observed) - 1e-12),
    loo_min = min(loo, na.rm = TRUE), loo_max = max(loo, na.rm = TRUE),
    loo_same_direction_fraction = mean(sign(loo) == sign(observed), na.rm = TRUE)
  )
}

one_sample_summary <- function(values) {
  values <- values[is.finite(values)]
  stopifnot(length(values) >= 2L)
  observed <- mean(values)
  boot <- replicate(n_boot, mean(sample(values, length(values), replace = TRUE)))
  loo <- vapply(seq_along(values), function(j) mean(values[-j]), numeric(1L))
  data.table(
    n_patients = length(values), mean = observed, median = median(values),
    bootstrap_ci_low = unname(quantile(boot, 0.025)), bootstrap_ci_high = unname(quantile(boot, 0.975)),
    loo_min = min(loo), loo_max = max(loo), loo_same_direction_fraction = mean(sign(loo) == sign(observed))
  )
}

match_regions <- function(conditional_score, raw_program, broad, log_total, log_detected, nn_index) {
  q <- quantile(conditional_score, c(0.25, 0.75), names = FALSE, type = 7)
  low <- which(conditional_score <= q[1L])
  high <- which(conditional_score >= q[2L])
  covariates <- cbind(z_safe(broad), z_safe(log_total), z_safe(log_detected))
  d_broad <- abs(outer(covariates[high, 1L], covariates[low, 1L], "-"))
  d_sq <- matrix(0, nrow = length(high), ncol = length(low))
  for (j in seq_len(ncol(covariates))) d_sq <- d_sq + outer(covariates[high, j], covariates[low, j], "-")^2
  eligible <- which(d_broad <= 0.50 & sqrt(d_sq) <= 1.00, arr.ind = TRUE)
  if (!nrow(eligible)) return(data.table(status = "NOT_ESTIMABLE", n_pairs = 0L))
  candidates <- data.table(
    high = high[eligible[, 1L]], low = low[eligible[, 2L]],
    distance = sqrt(d_sq[eligible]), broad_distance = d_broad[eligible]
  )
  setorder(candidates, distance, high, low)
  used_high <- logical(length(conditional_score))
  used_low <- logical(length(conditional_score))
  chosen <- integer()
  for (i in seq_len(nrow(candidates))) {
    h <- candidates$high[i]
    l <- candidates$low[i]
    if (!used_high[h] && !used_low[l]) {
      chosen <- c(chosen, i)
      used_high[h] <- TRUE
      used_low[l] <- TRUE
    }
  }
  pairs <- candidates[chosen]
  if (nrow(pairs) < 20L) return(data.table(status = "NOT_ESTIMABLE", n_pairs = nrow(pairs)))
  local_score <- rowMeans(matrix(conditional_score[nn_index], nrow = nrow(nn_index), ncol = ncol(nn_index)))
  data.table(
    status = "ESTIMATED", n_pairs = nrow(pairs),
    mean_pair_distance = mean(pairs$distance), max_pair_distance = max(pairs$distance),
    conditional_score_high_minus_low = mean(conditional_score[pairs$high] - conditional_score[pairs$low]),
    local_conditional_neighbor_high_minus_low = mean(local_score[pairs$high] - local_score[pairs$low]),
    raw_program_high_minus_low = mean(raw_program[pairs$high] - raw_program[pairs$low]),
    broad_immune_high_minus_low = mean(broad[pairs$high] - broad[pairs$low]),
    log_total_high_minus_low = mean(log_total[pairs$high] - log_total[pairs$low]),
    log_detected_high_minus_low = mean(log_detected[pairs$high] - log_detected[pairs$low])
  )
}

# -----------------------------------------------------------------------------
# GSE238264: reproduce v10 inputs, then calculate conditional organization.
# -----------------------------------------------------------------------------

gse_map <- data.table(
  sample_id = paste0("HCC", 1:7, c(rep("R", 4), rep("NR", 3))),
  patient_id = paste0("patient", 1:7),
  response = c(rep("Responder", 4), rep("Nonresponder", 3))
)

read_gse_counts <- function(sample_id) {
  matrix_dir <- file.path(gse_dir, sample_id, "filtered_feature_bc_matrix")
  features <- fread(file.path(matrix_dir, "features.tsv.gz"), header = FALSE)
  barcodes <- readLines(gzfile(file.path(matrix_dir, "barcodes.tsv.gz")))
  counts <- readMM(gzfile(file.path(matrix_dir, "matrix.mtx.gz")))
  genes <- make.unique(as.character(features[[2L]]))
  rownames(counts) <- genes
  colnames(counts) <- barcodes
  keep <- Matrix::colSums(counts) > 0
  counts[, keep, drop = FALSE]
}

cat("GSE238264_STATS_START\n")
gse_stats <- NULL
gse_total_spots <- 0L
for (sid in gse_map$sample_id) {
  counts <- read_gse_counts(sid)
  lib <- Matrix::colSums(counts)
  lognorm <- t(t(counts) * (10000 / lib))
  lognorm@x <- log1p(lognorm@x)
  if (is.null(gse_stats)) {
    gse_stats <- data.table(gene = rownames(lognorm), detected = 0, expression_sum = 0)
  } else stopifnot(identical(gse_stats$gene, rownames(lognorm)))
  gse_stats[, detected := detected + as.numeric(Matrix::rowSums(counts > 0))]
  gse_stats[, expression_sum := expression_sum + as.numeric(Matrix::rowSums(lognorm))]
  gse_total_spots <- gse_total_spots + ncol(lognorm)
  cat("GSE_STATS ", sid, " spots=", ncol(lognorm), "\n", sep = "")
}
gse_stats[, `:=`(
  detection_fraction = detected / gse_total_spots,
  mean_log1p_cp10k = expression_sum / gse_total_spots
)]
stopifnot(all(signature %chin% gse_stats$gene))
broad_gse <- intersect(broad_immune, gse_stats$gene)
stopifnot(identical(broad_gse, broad_immune))
gse_random <- make_random_sets(gse_stats, signature, unique(c(signature, broad_immune, "PTPN22")), "GSE238264")
gse_membership <- make_membership(gse_stats$gene, gse_random$sets, broad_gse)

gse_patient_rows <- vector("list", nrow(gse_map))
gse_regional_rows <- vector("list", nrow(gse_map))
gse_random_conditional_morans <- matrix(NA_real_, nrow = nrow(gse_map), ncol = n_random)

for (i in seq_len(nrow(gse_map))) {
  sid <- gse_map$sample_id[i]
  counts <- read_gse_counts(sid)
  total_counts <- as.numeric(Matrix::colSums(counts))
  detected_genes <- as.numeric(Matrix::colSums(counts > 0))
  lognorm <- t(t(counts) * (10000 / total_counts))
  lognorm@x <- log1p(lognorm@x)
  scores <- as.matrix(Matrix::t(lognorm) %*% gse_membership)
  positions <- fread(file.path(gse_dir, sid, "spatial", "tissue_positions_list.csv"), header = FALSE)
  setnames(positions, c("barcode", "in_tissue", "array_row", "array_col", "pixel_row", "pixel_col"))
  positions <- positions[match(colnames(lognorm), barcode)]
  stopifnot(!anyNA(positions$barcode))
  nn <- get.knn(as.matrix(positions[, .(array_row, array_col)]), k = 6L)$nn.index
  raw_morans <- moran_matrix(scores[, 1:2, drop = FALSE], nn)
  conditional <- conditional_residuals(
    scores[, c(1L, 3:ncol(scores)), drop = FALSE], scores[, 2L],
    log1p(total_counts), log1p(detected_genes)
  )
  conditional_morans <- moran_matrix(conditional$residual, nn)
  gse_random_conditional_morans[i, ] <- conditional_morans[-1L]
  gse_patient_rows[[i]] <- data.table(
    cohort = "GSE238264", record_type = "PATIENT", sample_id = sid,
    patient_id = gse_map$patient_id[i], response = gse_map$response[i], timing = "post_treatment",
    n_spots = nrow(scores), conditional_moran = conditional_morans[1L],
    raw_target_moran = raw_morans[1L], broad_immune_moran = raw_morans[2L],
    raw_target_minus_broad_moran = raw_morans[1L] - raw_morans[2L],
    conditional_random_moran_mean = mean(conditional_morans[-1L]),
    conditional_random_moran_sd = sd(conditional_morans[-1L]),
    conditional_random_null_z = (conditional_morans[1L] - mean(conditional_morans[-1L])) / sd(conditional_morans[-1L]),
    geometry = "within-slide directed k=6 nearest-neighbor graph", inference_unit = "patient"
  )
  regional <- match_regions(
    conditional$residual[, 1L], scores[, 1L], scores[, 2L],
    log1p(total_counts), log1p(detected_genes), nn
  )
  regional[, `:=`(
    cohort = "GSE238264", record_type = "PATIENT", sample_id = sid,
    patient_id = gse_map$patient_id[i], response = gse_map$response[i], timing = "post_treatment"
  )]
  gse_regional_rows[[i]] <- regional
  cat("GSE_SCORED ", sid, " raw=", format(raw_morans[1L], digits = 7),
      " broad=", format(raw_morans[2L], digits = 7),
      " conditional=", format(conditional_morans[1L], digits = 7),
      " matched_pairs=", regional$n_pairs, "\n", sep = "")
}
gse_patient <- rbindlist(gse_patient_rows, fill = TRUE)
gse_regional <- rbindlist(gse_regional_rows, fill = TRUE)

# -----------------------------------------------------------------------------
# Wu: reconstruct v10-corrected rows exactly and retain all samples for CODEX.
# -----------------------------------------------------------------------------

decode_categorical <- function(file, field) {
  categories <- h5read(file, paste0("/obs/", field, "/categories"))
  codes <- h5read(file, paste0("/obs/", field, "/codes"))
  categories[codes + 1L]
}

wu_obs <- data.table(
  spot_id = as.character(h5read(wu_h5, "/obs/index")),
  sample_id = as.character(decode_categorical(wu_h5, "sample_id")),
  state = as.character(decode_categorical(wu_h5, "diagnosis")),
  response_raw = as.character(decode_categorical(wu_h5, "Response"))
)
wu_obs[, patient_id := sub("^cytassist_([0-9]+)_(pre|post)$", "\\1", sample_id)]
wu_obs[, response := fifelse(response_raw == "Responder", "Responder", "Nonresponder")]
spatial <- h5read(wu_h5, "/obsm/spatial")
wu_obs[, `:=`(x = as.numeric(spatial[1L, ]), y = as.numeric(spatial[2L, ]))]
wu_genes <- as.character(h5read(wu_h5, "/var/_index"))
wu_indptr <- h5read(wu_h5, "/X/indptr")
stopifnot(length(wu_indptr) == nrow(wu_obs) + 1L, !anyDuplicated(wu_genes), nrow(wu_obs) == 104418L)
stopifnot(sum(diff(wu_indptr) == 0L) == 401L)

wu_map <- fread(wu_map_file)
stopifnot(uniqueN(wu_obs$sample_id) == 22L, sum(wu_obs$state == "Pre") > 0L)

cache_file <- file.path(project_dir, "03_processed_data", "final_harmonization", "WU_SPATIAL_GENE_MATCHING_STATS.rds")
stopifnot(file.exists(cache_file))
wu_stats <- readRDS(cache_file)
stopifnot(identical(wu_stats$gene, wu_genes), nrow(wu_stats) == length(wu_genes), all(signature %chin% wu_genes))
broad_wu <- intersect(broad_immune, wu_genes)
stopifnot(length(broad_wu) == 10L, identical(setdiff(broad_immune, broad_wu), "HLA-DRA"))
wu_random <- make_random_sets(wu_stats, signature, unique(c(signature, broad_immune, "PTPN22")), "Wu_Zenodo_19123188")
wu_membership <- make_membership(wu_genes, wu_random$sets, broad_wu)

matching_now <- rbindlist(list(gse_random$audit, wu_random$audit))
matching_v10 <- fread(v10_matching_file)
setorder(matching_now, cohort, target_gene)
setorder(matching_v10, cohort, target_gene)
matching_columns <- c("target_mean_log1p_cp10k", "matched_mean_log1p_cp10k_mean", "target_detection_fraction", "matched_detection_fraction_mean")
stopifnot(
  identical(matching_now[, .(cohort, target_gene)], matching_v10[, .(cohort, target_gene)]),
  max(abs(as.matrix(matching_now[, ..matching_columns]) - as.matrix(matching_v10[, ..matching_columns]))) < 1e-12
)
cat("V10_RANDOM_MATCHING_REPRODUCED\n")

sample_runs <- rle(wu_obs$sample_id)
stopifnot(length(sample_runs$values) == uniqueN(wu_obs$sample_id))
run_ends <- cumsum(sample_runs$lengths)
run_starts <- c(1L, head(run_ends, -1L) + 1L)
pre_run_indices <- which(wu_obs$state[run_starts] == "Pre")

wu_patient_rows <- list()
wu_regional_rows <- list()
wu_random_conditional_morans <- matrix(NA_real_, nrow = length(pre_run_indices), ncol = n_random)
wu_spot_rows <- vector("list", length(pre_run_indices))
pre_counter <- 0L

for (jj in seq_along(pre_run_indices)) {
  ii <- pre_run_indices[jj]
  first_row <- run_starts[ii]
  last_row <- run_ends[ii]
  ptr <- wu_indptr[first_row:(last_row + 1L)]
  first_entry <- ptr[1L] + 1L
  last_entry <- ptr[length(ptr)]
  if (last_entry >= first_entry) {
    idx <- as.integer(h5read(wu_h5, "/X/indices", index = list(first_entry:last_entry)))
    val <- as.numeric(h5read(wu_h5, "/X/data", index = list(first_entry:last_entry)))
    row_local <- rep.int(seq_len(last_row - first_row + 1L), diff(ptr))
  } else {
    idx <- integer()
    val <- numeric()
    row_local <- integer()
  }
  n_rows <- last_row - first_row + 1L
  total_counts <- csr_library_sizes(val, row_local, n_rows)
  detected_genes <- tabulate(row_local, nbins = n_rows)
  logval <- if (length(val)) log1p(val * 10000 / total_counts[row_local]) else numeric()
  selected <- wu_membership@i + 1L
  keep <- (idx + 1L) %in% selected
  block <- sparseMatrix(
    i = row_local[keep], j = idx[keep] + 1L, x = logval[keep],
    dims = c(n_rows, length(wu_genes))
  )
  scores <- as.matrix(block %*% wu_membership)
  conditional <- conditional_residuals(
    scores[, c(1L, 3:ncol(scores)), drop = FALSE], scores[, 2L],
    log1p(total_counts), log1p(detected_genes)
  )
  obs_block <- wu_obs[first_row:last_row]
  wu_spot_rows[[jj]] <- data.table(
    spot_id = obs_block$spot_id, sample_id = obs_block$sample_id,
    patient_id = obs_block$patient_id, state = obs_block$state, response = obs_block$response,
    conditional_score = conditional$residual[, 1L], raw_target_score = scores[, 1L],
    broad_immune_score = scores[, 2L], log_total_counts = log1p(total_counts),
    log_detected_genes = log1p(detected_genes), zero_count_spot = total_counts == 0
  )
  if (obs_block$state[1L] == "Pre") {
    pre_counter <- pre_counter + 1L
    nn <- get.knn(as.matrix(obs_block[, .(x, y)]), k = 6L)$nn.index
    raw_morans <- moran_matrix(scores[, 1:2, drop = FALSE], nn)
    conditional_morans <- moran_matrix(conditional$residual, nn)
    wu_random_conditional_morans[pre_counter, ] <- conditional_morans[-1L]
    wu_patient_rows[[pre_counter]] <- data.table(
      cohort = "Wu_Zenodo_19123188", record_type = "PATIENT", sample_id = obs_block$sample_id[1L],
      patient_id = obs_block$patient_id[1L], response = obs_block$response[1L], timing = "pretreatment",
      n_spots = n_rows, conditional_moran = conditional_morans[1L],
      raw_target_moran = raw_morans[1L], broad_immune_moran = raw_morans[2L],
      raw_target_minus_broad_moran = raw_morans[1L] - raw_morans[2L],
      conditional_random_moran_mean = mean(conditional_morans[-1L]),
      conditional_random_moran_sd = sd(conditional_morans[-1L]),
      conditional_random_null_z = (conditional_morans[1L] - mean(conditional_morans[-1L])) / sd(conditional_morans[-1L]),
      geometry = "within-sample directed k=6 nearest-neighbor graph", inference_unit = "patient"
    )
    regional <- match_regions(
      conditional$residual[, 1L], scores[, 1L], scores[, 2L],
      log1p(total_counts), log1p(detected_genes), nn
    )
    regional[, `:=`(
      cohort = "Wu_Zenodo_19123188", record_type = "PATIENT", sample_id = obs_block$sample_id[1L],
      patient_id = obs_block$patient_id[1L], response = obs_block$response[1L], timing = "pretreatment"
    )]
    wu_regional_rows[[pre_counter]] <- regional
    cat("WU_PRE_SCORED ", obs_block$sample_id[1L], " raw=", format(raw_morans[1L], digits = 7),
        " broad=", format(raw_morans[2L], digits = 7),
        " conditional=", format(conditional_morans[1L], digits = 7),
        " matched_pairs=", regional$n_pairs, "\n", sep = "")
  }
}

wu_patient <- rbindlist(wu_patient_rows, fill = TRUE)
wu_regional <- rbindlist(wu_regional_rows, fill = TRUE)
wu_spots <- rbindlist(wu_spot_rows, fill = TRUE)
stopifnot(nrow(wu_patient) == 11L, sum(wu_spots$zero_count_spot) == 397L)

# Runtime reproduction check: raw target/broad Moran values must equal frozen v10.
v10_patient <- fread(v10_patient_file)
current_raw <- rbindlist(list(gse_patient, wu_patient))
raw_check <- merge(
  current_raw[, .(cohort, sample_id, raw_target_moran, broad_immune_moran)],
  v10_patient[, .(cohort, sample_id, target_moran, broad_moran)],
  by = c("cohort", "sample_id"), all = TRUE
)
stopifnot(
  nrow(raw_check) == 18L, !anyNA(raw_check),
  max(abs(raw_check$raw_target_moran - raw_check$target_moran)) < 1e-12,
  max(abs(raw_check$broad_immune_moran - raw_check$broad_moran)) < 1e-12
)
cat("V10_RAW_SPATIAL_REPRODUCTION_PASS rows=18\n")

# -----------------------------------------------------------------------------
# Patient-level spatial summaries and immune-conditioned regional summaries.
# -----------------------------------------------------------------------------

summarize_spatial <- function(patient, random_morans) {
  conditional_effect <- group_effect(patient$conditional_moran, patient$response)
  raw_effect <- mean(patient$response == "Responder")
  is_r <- patient$response == "Responder"
  is_nr <- patient$response == "Nonresponder"
  random_response <- colMeans(random_morans[is_r, , drop = FALSE]) - colMeans(random_morans[is_nr, , drop = FALSE])
  random_cohort_mean <- colMeans(random_morans)
  response_null_z <- (conditional_effect$effect - mean(random_response)) / sd(random_response)
  cohort_mean <- mean(patient$conditional_moran)
  cohort_mean_null_z <- (cohort_mean - mean(random_cohort_mean)) / sd(random_cohort_mean)
  raw_target_effect <- mean(patient$raw_target_moran[is_r]) - mean(patient$raw_target_moran[is_nr])
  broad_effect <- mean(patient$broad_immune_moran[is_r]) - mean(patient$broad_immune_moran[is_nr])
  support <- conditional_effect$effect > 0 && response_null_z > 0 && conditional_effect$loo_same_direction_fraction == 1
  data.table(
    cohort = patient$cohort[1L], record_type = "COHORT_SUMMARY", sample_id = NA_character_, patient_id = NA_character_,
    response = NA_character_, timing = patient$timing[1L], n_spots = sum(patient$n_spots),
    n_patients = nrow(patient), n_responder = conditional_effect$n_responder, n_nonresponder = conditional_effect$n_nonresponder,
    conditional_moran = cohort_mean, raw_target_moran = mean(patient$raw_target_moran),
    broad_immune_moran = mean(patient$broad_immune_moran),
    raw_target_minus_broad_moran = mean(patient$raw_target_minus_broad_moran),
    conditional_R_minus_NR = conditional_effect$effect,
    conditional_R_mean = conditional_effect$responder_mean,
    conditional_NR_mean = conditional_effect$nonresponder_mean,
    bootstrap_ci_low = conditional_effect$bootstrap_ci_low,
    bootstrap_ci_high = conditional_effect$bootstrap_ci_high,
    exact_permutation_p_two_sided = conditional_effect$exact_permutation_p_two_sided,
    loo_min = conditional_effect$loo_min, loo_max = conditional_effect$loo_max,
    loo_same_direction_fraction = conditional_effect$loo_same_direction_fraction,
    raw_target_R_minus_NR = raw_target_effect, broad_immune_R_minus_NR = broad_effect,
    matched_random_response_null_mean = mean(random_response),
    matched_random_response_null_sd = sd(random_response),
    matched_random_response_null_z = response_null_z,
    matched_random_response_empirical_p = empirical_two_sided(conditional_effect$effect, random_response),
    conditional_cohort_mean = cohort_mean,
    matched_random_cohort_mean_null_mean = mean(random_cohort_mean),
    matched_random_cohort_mean_null_sd = sd(random_cohort_mean),
    matched_random_cohort_mean_null_z = cohort_mean_null_z,
    matched_random_cohort_mean_empirical_p = empirical_two_sided(cohort_mean, random_cohort_mean),
    conditional_support = support,
    geometry = patient$geometry[1L], inference_unit = "patient"
  )
}

gse_summary <- summarize_spatial(gse_patient, gse_random_conditional_morans)
wu_summary <- summarize_spatial(wu_patient, wu_random_conditional_morans)
spatial_effects <- rbindlist(list(gse_patient, wu_patient, gse_summary, wu_summary), fill = TRUE, use.names = TRUE)

summarize_regional <- function(regional) {
  estimated <- regional[status == "ESTIMATED"]
  if (!nrow(estimated)) {
    return(data.table(
      cohort = regional$cohort[1L], record_type = "COHORT_SUMMARY", status = "NOT_ESTIMABLE",
      n_patients = nrow(regional), n_patients_estimated = 0L
    ))
  }
  effect <- group_effect(estimated$local_conditional_neighbor_high_minus_low, estimated$response)
  data.table(
    cohort = regional$cohort[1L], record_type = "COHORT_SUMMARY", status = "ESTIMATED",
    n_patients = nrow(regional), n_patients_estimated = nrow(estimated), n_pairs = sum(estimated$n_pairs),
    mean_pair_distance = mean(estimated$mean_pair_distance), max_pair_distance = max(estimated$max_pair_distance),
    conditional_score_high_minus_low = mean(estimated$conditional_score_high_minus_low),
    local_conditional_neighbor_high_minus_low = mean(estimated$local_conditional_neighbor_high_minus_low),
    raw_program_high_minus_low = mean(estimated$raw_program_high_minus_low),
    broad_immune_high_minus_low = mean(estimated$broad_immune_high_minus_low),
    log_total_high_minus_low = mean(estimated$log_total_high_minus_low),
    log_detected_high_minus_low = mean(estimated$log_detected_high_minus_low),
    regional_local_effect_R_minus_NR = effect$effect,
    bootstrap_ci_low = effect$bootstrap_ci_low, bootstrap_ci_high = effect$bootstrap_ci_high,
    exact_permutation_p_two_sided = effect$exact_permutation_p_two_sided,
    loo_min = effect$loo_min, loo_max = effect$loo_max,
    loo_same_direction_fraction = effect$loo_same_direction_fraction,
    inference_unit = "patient"
  )
}

regional_effects <- rbindlist(list(
  gse_regional, wu_regional,
  summarize_regional(gse_regional), summarize_regional(wu_regional)
), fill = TRUE, use.names = TRUE)

fwrite(spatial_effects, file.path(table_dir, "PTPN22_CONDITIONAL_SPATIAL_EFFECTS.tsv"), sep = "\t", na = "NA")
fwrite(regional_effects, file.path(table_dir, "PTPN22_IMMUNE_MATCHED_REGIONAL_EFFECTS.tsv"), sep = "\t", na = "NA")

cat("\nCONDITIONAL SPATIAL SUMMARIES\n")
print(spatial_effects[record_type == "COHORT_SUMMARY"])
cat("\nIMMUNE-CONDITIONED REGIONAL SUMMARIES\n")
print(regional_effects[record_type == "COHORT_SUMMARY"])

# -----------------------------------------------------------------------------
# Pretreatment Wu ST-CODEX conditional neighborhood analysis.
# -----------------------------------------------------------------------------

spot_summary <- as.data.table(read_parquet(spot_summary_file))
spot_summary[, patient_id := as.character(patient_id)]
spot_summary <- spot_summary[state == "Pre"]
release_metadata <- fread(release_metadata_file)
release_metadata[, patient_id := as.character(patient_id)]
region_meta <- unique(release_metadata[, .(region_id, region_label, patient_id, state, response, treatment)])
coreg_regions <- unique(spot_summary$region_id)
stopifnot(all(coreg_regions %chin% region_meta$region_id))

cell_type_levels <- c(
  "B cells", "CD4 T cells", "CD8 T cells", "Dendritic cells", "Endothelial cells", "Fibroblasts", "INFg+",
  "Macrophages", "Macrophages M2-like", "Neutrophils", "Stroma Uncharacterized",
  "Epithelium (INOS-)", "Epithelium (INOS+)", "Unknown"
)
t_cell_types <- c("CD4 T cells", "CD8 T cells")
myeloid_types <- c("Dendritic cells", "Macrophages", "Macrophages M2-like", "Neutrophils")
file_for <- function(region, suffix) file.path(codex_raw_dir, paste0(region, suffix))

codex_spot_rows <- vector("list", length(coreg_regions))
for (ri in seq_along(coreg_regions)) {
  region <- coreg_regions[ri]
  cell_data <- fread(file_for(region, ".cell_data.csv"))
  cell_types <- fread(file_for(region, ".cell_types.csv"))
  expression <- fread(file_for(region, ".expression.csv"))
  stopifnot(
    identical(cell_data$CELL_ID, cell_types$CELL_ID),
    identical(cell_data$CELL_ID, expression$CELL_ID),
    all(cell_types$CELL_TYPE %chin% cell_type_levels),
    all(c("LAG3", "TOX") %chin% names(expression))
  )
  sub <- spot_summary[region_id == region]
  score_match <- match(sub$spot_id, wu_spots$spot_id)
  stopifnot(!anyNA(score_match))
  mapped_ids <- lapply(sub$cell_ids_expanded, function(x) as.integer(unlist(x)))
  n_mapped <- lengths(mapped_ids)
  all_mapped <- unlist(mapped_ids, use.names = FALSE)
  cell_row <- match(all_mapped, cell_data$CELL_ID)
  stopifnot(!anyNA(cell_row))
  spot_rep <- rep.int(seq_len(nrow(sub)), n_mapped)
  incidence <- sparseMatrix(i = spot_rep, j = cell_row, x = 1, dims = c(nrow(sub), nrow(cell_data)))
  t_indicator <- as.numeric(cell_types$CELL_TYPE %chin% t_cell_types)
  myeloid_indicator <- as.numeric(cell_types$CELL_TYPE %chin% myeloid_types)
  t_fraction <- as.numeric(incidence %*% t_indicator) / pmax(n_mapped, 1L)
  myeloid_fraction <- as.numeric(incidence %*% myeloid_indicator) / pmax(n_mapped, 1L)
  marker_means <- as.matrix(incidence %*% as.matrix(expression[, .(LAG3, TOX)])) / pmax(n_mapped, 1L)
  codex_spot_rows[[ri]] <- data.table(
    region_id = region, spot_id = sub$spot_id,
    patient_id = wu_spots$patient_id[score_match], state = wu_spots$state[score_match],
    response = wu_spots$response[score_match], sample_id = wu_spots$sample_id[score_match],
    conditional_score = wu_spots$conditional_score[score_match],
    broad_immune_score = wu_spots$broad_immune_score[score_match],
    log_total_counts = wu_spots$log_total_counts[score_match],
    log_detected_genes = wu_spots$log_detected_genes[score_match],
    n_mapped_cells = n_mapped, eligible_primary = n_mapped >= 5L,
    T_cells = t_fraction, myeloid = myeloid_fraction,
    LAG3 = marker_means[, 1L], TOX = marker_means[, 2L]
  )
  meta <- region_meta[region_id == region]
  expected_response <- fifelse(meta$response[1L] == "Responder", "Responder", "Nonresponder")
  stopifnot(
    uniqueN(codex_spot_rows[[ri]]$patient_id) == 1L,
    codex_spot_rows[[ri]]$patient_id[1L] == as.character(meta$patient_id[1L]),
    uniqueN(codex_spot_rows[[ri]]$response) == 1L,
    codex_spot_rows[[ri]]$response[1L] == expected_response
  )
  cat("CODEX_REGION ", region, " spots=", nrow(sub), " eligible=", sum(n_mapped >= 5L), "\n", sep = "")
}
codex_spots <- rbindlist(codex_spot_rows, fill = TRUE)
codex_pre <- codex_spots[state == "Pre" & eligible_primary == TRUE]
stopifnot(uniqueN(codex_pre$patient_id) == 10L, all(c("Responder", "Nonresponder") %chin% codex_pre$response))

fit_patient_feature <- function(d, feature) {
  required <- c("conditional_score", "broad_immune_score", "n_mapped_cells", feature)
  if (feature != "T_cells") required <- c(required, "T_cells")
  keep <- complete.cases(d[, ..required])
  d <- d[keep]
  if (nrow(d) < 20L || sd(d[[feature]]) == 0 || sd(d$conditional_score) == 0) {
    return(data.table(status = "NOT_ESTIMABLE", beta_conditional = NA_real_, n_spots = nrow(d), n_regions = uniqueN(d$region_id)))
  }
  model_data <- data.frame(
    y = z_safe(d[[feature]]), cond = z_safe(d$conditional_score), broad = z_safe(d$broad_immune_score),
    log_mapped = z_safe(log1p(d$n_mapped_cells)), region = factor(d$region_id)
  )
  has_multiple_regions <- nlevels(model_data$region) > 1L
  if (feature == "T_cells" && has_multiple_regions) {
    fit <- lm(y ~ cond + broad + log_mapped + region, data = model_data)
  } else if (feature == "T_cells") {
    fit <- lm(y ~ cond + broad + log_mapped, data = model_data)
  } else {
    model_data$tcell <- z_safe(d$T_cells)
    if (has_multiple_regions) {
      fit <- lm(y ~ cond + broad + tcell + log_mapped + region, data = model_data)
    } else {
      fit <- lm(y ~ cond + broad + tcell + log_mapped, data = model_data)
    }
  }
  beta <- unname(coef(fit)["cond"])
  data.table(
    status = ifelse(is.finite(beta), "ESTIMATED", "NOT_ESTIMABLE"), beta_conditional = beta,
    n_spots = nrow(d), n_regions = uniqueN(d$region_id), model_rank = fit$rank, residual_df = df.residual(fit)
  )
}

features <- c("LAG3", "TOX", "myeloid", "T_cells")
codex_patient_rows <- list()
row_counter <- 0L
for (patient in sort(unique(codex_pre$patient_id))) {
  d <- codex_pre[patient_id == patient]
  for (feature in features) {
    row_counter <- row_counter + 1L
    fit_row <- fit_patient_feature(d, feature)
    fit_row[, `:=`(
      record_type = "PATIENT", patient_id = patient, response = d$response[1L],
      feature = feature,
      feature_type = fifelse(feature %chin% c("LAG3", "TOX"), "released_marker_mean", "cell_type_fraction"),
      inference_unit = "patient; region fixed effects within patient"
    )]
    codex_patient_rows[[row_counter]] <- fit_row
  }
}
codex_patient <- rbindlist(codex_patient_rows, fill = TRUE)

codex_summary_rows <- list()
summary_counter <- 0L
for (feature_name in features) {
  d <- codex_patient[feature == feature_name & status == "ESTIMATED"]
  all_summary <- one_sample_summary(d$beta_conditional)
  summary_counter <- summary_counter + 1L
  codex_summary_rows[[summary_counter]] <- data.table(
    record_type = "ALL_PATIENT_SUMMARY", feature = feature_name,
    feature_type = d$feature_type[1L], n_patients = all_summary$n_patients,
    mean_beta_conditional = all_summary$mean, median_beta_conditional = all_summary$median,
    bootstrap_ci_low = all_summary$bootstrap_ci_low, bootstrap_ci_high = all_summary$bootstrap_ci_high,
    loo_min = all_summary$loo_min, loo_max = all_summary$loo_max,
    loo_same_direction_fraction = all_summary$loo_same_direction_fraction,
    inference_unit = "patient"
  )
  response_summary <- group_effect(d$beta_conditional, d$response)
  summary_counter <- summary_counter + 1L
  codex_summary_rows[[summary_counter]] <- data.table(
    record_type = "RESPONSE_SUMMARY", feature = feature_name,
    feature_type = d$feature_type[1L], n_patients = nrow(d),
    n_responder = response_summary$n_responder, n_nonresponder = response_summary$n_nonresponder,
    responder_mean_beta = response_summary$responder_mean,
    nonresponder_mean_beta = response_summary$nonresponder_mean,
    responder_minus_nonresponder_beta = response_summary$effect,
    bootstrap_ci_low = response_summary$bootstrap_ci_low, bootstrap_ci_high = response_summary$bootstrap_ci_high,
    exact_permutation_p_two_sided = response_summary$exact_permutation_p_two_sided,
    loo_min = response_summary$loo_min, loo_max = response_summary$loo_max,
    loo_same_direction_fraction = response_summary$loo_same_direction_fraction,
    inference_unit = "patient"
  )
}
codex_effects <- rbindlist(c(list(codex_patient), codex_summary_rows), fill = TRUE, use.names = TRUE)
fwrite(codex_effects, file.path(table_dir, "PTPN22_CODEX_CONDITIONAL_EFFECTS.tsv"), sep = "\t", na = "NA")

safe_spearman <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 4L || sd(x[keep]) == 0 || sd(y[keep]) == 0) return(NA_real_)
  suppressWarnings(cor(x[keep], y[keep], method = "spearman"))
}

crossmodal_rows <- list()
for (i in seq_along(features)) {
  feature_name <- features[i]
  d <- merge(
    wu_patient[, .(patient_id, response, conditional_moran)],
    codex_patient[feature == feature_name & status == "ESTIMATED", .(patient_id, beta_conditional)],
    by = "patient_id"
  )
  rho <- safe_spearman(d$conditional_moran, d$beta_conditional)
  boot <- replicate(n_boot, {
    ind <- sample(seq_len(nrow(d)), nrow(d), replace = TRUE)
    safe_spearman(d$conditional_moran[ind], d$beta_conditional[ind])
  })
  boot <- boot[is.finite(boot)]
  loo <- vapply(seq_len(nrow(d)), function(j) safe_spearman(d$conditional_moran[-j], d$beta_conditional[-j]), numeric(1L))
  crossmodal_rows[[i]] <- data.table(
    feature = feature_name, n_patients = nrow(d), n_responder = sum(d$response == "Responder"),
    n_nonresponder = sum(d$response == "Nonresponder"), spearman_rho = rho,
    bootstrap_ci_low = unname(quantile(boot, 0.025)), bootstrap_ci_high = unname(quantile(boot, 0.975)),
    loo_min = min(loo, na.rm = TRUE), loo_max = max(loo, na.rm = TRUE),
    loo_same_direction_fraction = mean(sign(loo) == sign(rho), na.rm = TRUE),
    inference_unit = "patient", caveat = "small-n descriptive cross-modal concordance; adjacent sections"
  )
}
crossmodal <- rbindlist(crossmodal_rows)
fwrite(crossmodal, file.path(table_dir, "PTPN22_ST_CODEX_CROSSMODAL_CONCORDANCE.tsv"), sep = "\t", na = "NA")

cat("\nCONDITIONAL CODEX ALL-PATIENT SUMMARIES\n")
print(codex_effects[record_type == "ALL_PATIENT_SUMMARY"])
cat("\nCONDITIONAL CODEX RESPONSE SUMMARIES\n")
print(codex_effects[record_type == "RESPONSE_SUMMARY"])
cat("\nST-CODEX CROSSMODAL CONCORDANCE\n")
print(crossmodal)

# -----------------------------------------------------------------------------
# Prespecified gate inputs and falsification ledger.
# -----------------------------------------------------------------------------

spatial_gate_rows <- spatial_effects[record_type == "COHORT_SUMMARY", .(
  evidence_type = "SPATIAL_COHORT", evidence_id = cohort,
  observed_effect = conditional_R_minus_NR,
  matched_random_null_z = matched_random_response_null_z,
  loo_min, loo_max, loo_same_direction_fraction,
  positive_loo_stable = conditional_R_minus_NR > 0 & loo_same_direction_fraction == 1,
  counts_for_upgrade = conditional_support,
  comparator = paste0("raw_target_R_minus_NR=", format(raw_target_R_minus_NR, digits = 8),
    ";broad_immune_R_minus_NR=", format(broad_immune_R_minus_NR, digits = 8))
)]

codex_primary <- codex_effects[
  record_type == "ALL_PATIENT_SUMMARY" & feature %chin% c("LAG3", "TOX", "myeloid"),
  .(
    evidence_type = "CONDITIONAL_CODEX", evidence_id = feature,
    observed_effect = mean_beta_conditional,
    matched_random_null_z = NA_real_, loo_min, loo_max, loo_same_direction_fraction,
    positive_loo_stable = mean_beta_conditional > 0 & loo_min > 0,
    counts_for_upgrade = mean_beta_conditional > 0 & loo_min > 0,
    comparator = "adjusted for broad immune, total T-cell fraction, technical covariates and region"
  )
]

tcell_falsification <- codex_effects[
  record_type == "ALL_PATIENT_SUMMARY" & feature == "T_cells",
  .(
    evidence_type = "T_CELL_FALSIFICATION", evidence_id = feature,
    observed_effect = mean_beta_conditional,
    matched_random_null_z = NA_real_, loo_min, loo_max, loo_same_direction_fraction,
    positive_loo_stable = mean_beta_conditional > 0 & loo_min > 0,
    counts_for_upgrade = FALSE,
    comparator = "total T-cell fraction cannot count toward rescue support"
  )
]

gate_inputs <- rbindlist(list(spatial_gate_rows, codex_primary, tcell_falsification), fill = TRUE)
n_spatial_support <- sum(spatial_gate_rows$counts_for_upgrade)
n_codex_support <- sum(codex_primary$counts_for_upgrade)
if (n_spatial_support == 2L || (n_spatial_support == 1L && n_codex_support >= 2L)) {
  rescue_gate <- "SPATIAL_SPECIFICITY_UPGRADED"
} else if (n_spatial_support == 0L && n_codex_support == 0L) {
  rescue_gate <- "SPATIAL_SPECIFICITY_NOT_ADDITIONAL"
} else {
  rescue_gate <- "SPATIAL_SPECIFICITY_PARTIAL_CONFIRMED"
}
gate_inputs[, `:=`(
  n_spatial_cohorts_supporting = n_spatial_support,
  n_primary_codex_features_supporting = n_codex_support,
  rescue_gate = rescue_gate,
  criteria_frozen_before_results = TRUE
)]
fwrite(gate_inputs, file.path(table_dir, "PTPN22_SPATIAL_RESCUE_GATE_INPUTS.tsv"), sep = "\t", na = "NA")

falsification <- rbindlist(list(
  spatial_effects[record_type == "COHORT_SUMMARY", .(
    control = "broad_immune_and_matched_random_programs", cohort,
    target_conditional_R_minus_NR = conditional_R_minus_NR,
    raw_target_R_minus_NR, broad_immune_R_minus_NR,
    matched_random_response_null_z, matched_random_response_empirical_p,
    matched_random_cohort_mean_null_z, matched_random_cohort_mean_empirical_p,
    patient_loo_same_direction_fraction = loo_same_direction_fraction,
    response_blind_score_construction = TRUE, patient_level_inference = TRUE
  )],
  codex_effects[record_type == "ALL_PATIENT_SUMMARY", .(
    control = fifelse(feature == "T_cells", "total_T_cell_fraction_comparator", "conditional_CODEX_primary_feature"),
    cohort = "Wu_Zenodo_19123188", feature,
    conditional_mean_beta = mean_beta_conditional,
    bootstrap_ci_low, bootstrap_ci_high,
    patient_loo_same_direction_fraction = loo_same_direction_fraction,
    response_blind_score_construction = TRUE, patient_level_inference = TRUE
  )]
), fill = TRUE, use.names = TRUE)
fwrite(falsification, file.path(table_dir, "PTPN22_SPATIAL_RESCUE_FALSIFICATION_INPUTS.tsv"), sep = "\t", na = "NA")

cat("\nPRESPECIFIED GATE INPUTS\n")
print(gate_inputs)
cat("\nRESCUE_GATE: ", rescue_gate, "\n", sep = "")

output_files <- c(
  file.path(table_dir, "PTPN22_CONDITIONAL_SPATIAL_EFFECTS.tsv"),
  file.path(table_dir, "PTPN22_IMMUNE_MATCHED_REGIONAL_EFFECTS.tsv"),
  file.path(table_dir, "PTPN22_CODEX_CONDITIONAL_EFFECTS.tsv"),
  file.path(table_dir, "PTPN22_ST_CODEX_CROSSMODAL_CONCORDANCE.tsv"),
  file.path(table_dir, "PTPN22_SPATIAL_RESCUE_GATE_INPUTS.tsv"),
  file.path(table_dir, "PTPN22_SPATIAL_RESCUE_FALSIFICATION_INPUTS.tsv")
)
cat("\nOUTPUT HASHES\n")
for (path in output_files) cat(digest(file = path, algo = "sha256"), " ", path, "\n", sep = "")
cat("\nSESSION INFO\n")
print(sessionInfo())
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
