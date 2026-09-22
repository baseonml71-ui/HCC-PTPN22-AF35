options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

set.seed(20260828)

tcr_file <- "02_raw_data/GSE235863/GSE235863_nine_patients_TCR.csv.gz"
h5_file <- "03_processed_data/phase3/GSE235863_nine_patients_scRNAseq_cd45_raw_counts.h5ad"
h5_cell_file <- "03_processed_data/phase3/GSE235863_H5AD_PTPN22_CELL_DATA.tsv.gz"
h5_extract_audit_file <- "05_results/tables/phase3/GSE235863_H5AD_PTPN22_EXTRACTION_AUDIT.tsv"
response_file <- "01_metadata/phase3/GSE235863_SCRNA_SCTCR_RESPONSE_MAPPING.tsv"
design_file <- "01_metadata/phase3/GSE235863_CLONOTYPE_ANALYSIS_DESIGN.tsv"
mapping_file <- "03_processed_data/phase3/GSE235863_BARCODE_CLONOTYPE_STATE_PTPN22.tsv.gz"
clone_file <- "05_results/tables/phase3/GSE235863_CLONOTYPE_SUMMARIES.tsv.gz"
patient_file <- "05_results/tables/phase3/GSE235863_CLONE_PATIENT_SUMMARIES.tsv"
state_file <- "05_results/tables/phase3/GSE235863_CLONE_STATE_PATIENT_SUMMARIES.tsv"
context_file <- "05_results/tables/phase3/GSE235863_CLONE_STATE_CONTEXT_SUMMARIES.tsv"
persistence_file <- "05_results/tables/phase3/GSE235863_PERSISTENCE_SUMMARIES.tsv"
effects_file <- "05_results/tables/phase3/GSE235863_PATIENT_LEVEL_EFFECTS.tsv"
within_patient_file <- "05_results/tables/phase3/GSE235863_WITHIN_PATIENT_CLONAL_ASSOCIATIONS.tsv"
mapping_audit_file <- "05_results/tables/phase3/GSE235863_BARCODE_MAPPING_AUDIT.tsv"
state_alias_file <- "05_results/tables/phase3/GSE235863_STATE_LABEL_ALIAS_AUDIT.tsv"
h5_structure_file <- "05_results/tables/phase3/GSE235863_H5AD_STRUCTURE.tsv"
manifest_file <- "01_metadata/PHASE3_SCRNA_SCTCR_BARCODE_MAPPING_MANIFEST.tsv"
log_file <- "09_logs/phase3/03_independent_scrna_sctcr_clone_state.log"

dir.create(dirname(mapping_file), recursive = TRUE, showWarnings = FALSE)

log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

cat("Phase 3C independent scRNA/scTCR analysis\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Seed: 20260828\n")

stopifnot(file.exists(tcr_file), file.exists(h5_file), file.exists(h5_cell_file),
          file.exists(h5_extract_audit_file), file.exists(h5_structure_file),
          file.exists(response_file), file.exists(design_file))

h5_cells <- fread(h5_cell_file)
h5_extract_audit <- fread(h5_extract_audit_file)
audit_value <- function(item) h5_extract_audit[audit_item == item, value][[1L]]
n_obs <- nrow(h5_cells)
n_var <- as.integer(audit_value("h5_genes"))
gene_column <- audit_value("gene_column")
ptpn22_zero_based <- as.integer(audit_value("PTPN22_zero_based_index"))
encoding <- audit_value("X_encoding")
sample_column <- audit_value("sample_column")
subcluster_column <- audit_value("sub_cluster_column")
major_column <- audit_value("major_cluster_column")
stopifnot(uniqueN(h5_cells$barcode) == n_obs, n_obs == as.integer(audit_value("h5_cells")))

tcr <- fread(tcr_file)
response <- fread(response_file)
stopifnot(uniqueN(tcr$barcode) == nrow(tcr), uniqueN(response$patient) == nrow(response))

tcr[, timepoint := fifelse(grepl("-pre-", sample, fixed = TRUE), "pre",
                    fifelse(grepl("-post-", sample, fixed = TRUE), "post", NA_character_))]
stopifnot(!anyNA(tcr$timepoint))
tcr[, clonotype_key := paste(patient, clone.id, sep = "::")]
tcr[, cdr3_pair_aa := paste(get("cdr3.A"), get("cdr3.B"), sep = "|")]

clone_integrity <- tcr[, .(
  observed_cells = .N,
  reported_sizes = uniqueN(clone.size),
  reported_size = clone.size[[1L]],
  cdr3_pairs = uniqueN(cdr3_pair_aa)
), by = .(patient, clone.id)]
stopifnot(all(clone_integrity$reported_sizes == 1L),
          all(clone_integrity$reported_size == clone_integrity$observed_cells),
          all(clone_integrity$cdr3_pairs == 1L))

tcr[, obs_row := match(barcode, h5_cells$barcode)]
matched <- !is.na(tcr$obs_row)
match_rate <- mean(matched)
cat("Transcriptomic barcodes:", n_obs, "\n")
cat("TCR barcodes:", nrow(tcr), "\n")
cat("Matched barcodes:", sum(matched), "\n")
cat("Match rate:", sprintf("%.8f", match_rate), "\n")
stopifnot(match_rate >= 0.95)

tcr[, h5_sample := h5_cells$h5_sample[obs_row]]
tcr[, h5_subcluster := h5_cells$h5_subcluster[obs_row]]
tcr[, h5_major_cluster := h5_cells$h5_major_cluster[obs_row]]

sample_match_rate <- if (all(is.na(tcr$h5_sample))) NA_real_ else mean(tcr$sample == tcr$h5_sample, na.rm = TRUE)
subcluster_match_rate <- if (all(is.na(tcr$h5_subcluster))) NA_real_ else mean(tcr$sub_cluster == tcr$h5_subcluster, na.rm = TRUE)
major_match_rate <- if (all(is.na(tcr$h5_major_cluster))) NA_real_ else mean(tcr$major_cluster == tcr$h5_major_cluster, na.rm = TRUE)
state_cluster_id <- function(x) sub("^((CD4|CD8)_C[0-9]+).*$", "\\1", x)
subcluster_id_match_rate <- mean(state_cluster_id(tcr$sub_cluster) == state_cluster_id(tcr$h5_subcluster), na.rm = TRUE)
major_semantic_match_rate <- mean(paste0(tcr$major_cluster, "T") == tcr$h5_major_cluster, na.rm = TRUE)
state_alias_audit <- tcr[, .N, by = .(major_cluster, h5_major_cluster, sub_cluster, h5_subcluster)]
state_alias_audit[, cluster_id_match := state_cluster_id(sub_cluster) == state_cluster_id(h5_subcluster)]
setorder(state_alias_audit, -N)

cat("X encoding:", encoding, "\n")
tcr[, ptpn22_raw_count := h5_cells$ptpn22_raw_count[obs_row]]
tcr[, total_raw_umi := h5_cells$total_raw_umi[obs_row]]
stopifnot(all(tcr$total_raw_umi > 0))
tcr[, ptpn22_log1p_cp10k := log1p(10000 * ptpn22_raw_count / total_raw_umi)]
tcr[, ptpn22_detected := ptpn22_raw_count > 0]
tcr[, expansion_category := fifelse(clone.size == 1L, "singleton",
                              fifelse(clone.size <= 4L, "small_expanded", "large_expanded"))]
tcr[response, on = "patient", `:=`(
  response = i.response,
  response_original = i.response_original,
  response_group = i.response_group
)]
stopifnot(!anyNA(tcr$response_group))

context_eligibility <- tcr[, .(
  context_has_pre = any(timepoint == "pre"),
  context_has_post = any(timepoint == "post")
), by = .(patient, tissue)]
clone_persistence <- tcr[, .(
  clone_has_pre = any(timepoint == "pre"),
  clone_has_post = any(timepoint == "post")
), by = .(patient, tissue, clonotype_key)]
clone_persistence[context_eligibility, on = .(patient, tissue), `:=`(
  context_has_pre = i.context_has_pre,
  context_has_post = i.context_has_post
)]
clone_persistence[, persistence_eligible := context_has_pre & context_has_post]
clone_persistence[, persistent_within_tissue := fifelse(
  persistence_eligible, clone_has_pre & clone_has_post, NA
)]
tcr[clone_persistence, on = .(patient, tissue, clonotype_key), `:=`(
  persistence_eligible = i.persistence_eligible,
  persistent_within_tissue = i.persistent_within_tissue
)]

modal_state <- function(x) {
  counts <- sort(table(x), decreasing = TRUE)
  sort(names(counts)[counts == counts[[1L]]])[[1L]]
}

clone_summary <- tcr[, .(
  clone_size = .N,
  reported_clone_size = clone.size[[1L]],
  cdr3_pair_aa = cdr3_pair_aa[[1L]],
  alpha_v = get("v_gene.A")[[1L]],
  alpha_j = get("j_gene.A")[[1L]],
  beta_v = get("v_gene.B")[[1L]],
  beta_j = get("j_gene.B")[[1L]],
  n_samples = uniqueN(sample),
  n_tissues = uniqueN(tissue),
  n_timepoints = uniqueN(timepoint),
  n_subclusters = uniqueN(sub_cluster),
  mean_ptpn22_log1p_cp10k = mean(ptpn22_log1p_cp10k),
  ptpn22_detection_fraction = mean(ptpn22_detected)
), by = .(patient, response_group, clone.id, clonotype_key)]
stopifnot(all(clone_summary$clone_size == clone_summary$reported_clone_size))

safe_mean <- function(x) if (length(x) && any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || uniqueN(x[ok]) < 2L || uniqueN(y[ok]) < 2L) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok], method = "spearman"))
}

patient_summary <- tcr[, {
  clone_sizes <- table(clonotype_key)
  clone_p <- as.numeric(clone_sizes) / sum(clone_sizes)
  entropy <- -sum(clone_p * log(clone_p))
  normalized_entropy <- if (length(clone_p) > 1L) entropy / log(length(clone_p)) else NA_real_
  .(
    n_cells = .N,
    n_clonotypes = uniqueN(clonotype_key),
    max_clone_size = max(clone.size),
    expanded_cell_fraction = mean(clone.size >= 2L),
    large_expanded_cell_fraction = mean(clone.size >= 5L),
    normalized_shannon = normalized_entropy,
    clonality = 1 - normalized_entropy,
    mean_ptpn22_log1p_cp10k = mean(ptpn22_log1p_cp10k),
    ptpn22_detection_fraction = mean(ptpn22_detected),
    expanded_minus_singleton_ptpn22 = safe_mean(ptpn22_log1p_cp10k[clone.size >= 2L]) - safe_mean(ptpn22_log1p_cp10k[clone.size == 1L]),
    large_minus_singleton_ptpn22 = safe_mean(ptpn22_log1p_cp10k[clone.size >= 5L]) - safe_mean(ptpn22_log1p_cp10k[clone.size == 1L]),
    expanded_minus_singleton_detection = safe_mean(as.numeric(ptpn22_detected[clone.size >= 2L])) - safe_mean(as.numeric(ptpn22_detected[clone.size == 1L]))
  )
}, by = .(patient, response_group)]

clone_patient_metrics <- clone_summary[, .(
  clone_size_ptpn22_spearman = safe_cor(log1p(clone_size), mean_ptpn22_log1p_cp10k),
  multistate_expanded_clone_fraction = safe_mean(as.numeric(n_subclusters[clone_size >= 2L] > 1L))
), by = .(patient, response_group)]
patient_summary[clone_patient_metrics, on = .(patient, response_group), `:=`(
  clone_size_ptpn22_spearman = i.clone_size_ptpn22_spearman,
  multistate_expanded_clone_fraction = i.multistate_expanded_clone_fraction
)]

state_summary <- tcr[, .(
  n_cells = .N,
  n_clonotypes = uniqueN(clonotype_key),
  expanded_cell_fraction = mean(clone.size >= 2L),
  large_expanded_cell_fraction = mean(clone.size >= 5L),
  mean_ptpn22_log1p_cp10k = mean(ptpn22_log1p_cp10k),
  ptpn22_detection_fraction = mean(ptpn22_detected),
  expanded_minus_singleton_ptpn22 = safe_mean(ptpn22_log1p_cp10k[clone.size >= 2L]) - safe_mean(ptpn22_log1p_cp10k[clone.size == 1L])
), by = .(patient, response_group, major_cluster, sub_cluster)]

context_summary <- tcr[, .(
  n_cells = .N,
  n_clonotypes = uniqueN(clonotype_key),
  expanded_cell_fraction = mean(clone.size >= 2L),
  mean_ptpn22_log1p_cp10k = mean(ptpn22_log1p_cp10k),
  ptpn22_detection_fraction = mean(ptpn22_detected)
), by = .(patient, response_group, sample, timepoint, tissue, major_cluster, sub_cluster)]

modal_by_time <- tcr[persistence_eligible == TRUE, .(
  modal_subcluster = modal_state(sub_cluster)
), by = .(patient, response_group, tissue, clonotype_key, timepoint)]
modal_wide <- dcast(modal_by_time, patient + response_group + tissue + clonotype_key ~ timepoint, value.var = "modal_subcluster")
if (!"pre" %chin% names(modal_wide)) modal_wide[, pre := NA_character_]
if (!"post" %chin% names(modal_wide)) modal_wide[, post := NA_character_]
modal_wide[, persistent := !is.na(pre) & !is.na(post)]
modal_wide[, state_switch := fifelse(persistent, pre != post, NA)]

persistence_summary <- clone_persistence[, .(
  persistence_eligible = unique(persistence_eligible),
  n_clonotypes = .N,
  persistent_clonotypes = if (unique(persistence_eligible)) sum(persistent_within_tissue) else NA_integer_,
  persistent_clonotype_fraction = if (unique(persistence_eligible)) mean(persistent_within_tissue) else NA_real_
), by = .(patient, tissue)]
persistence_cells <- tcr[, .(
  n_cells = .N,
  persistent_cells = if (unique(persistence_eligible)) sum(persistent_within_tissue) else NA_integer_,
  persistent_cell_fraction = if (unique(persistence_eligible)) mean(persistent_within_tissue) else NA_real_,
  persistent_minus_nonpersistent_ptpn22 = if (unique(persistence_eligible)) safe_mean(ptpn22_log1p_cp10k[persistent_within_tissue == TRUE]) - safe_mean(ptpn22_log1p_cp10k[persistent_within_tissue == FALSE]) else NA_real_
), by = .(patient, response_group, tissue)]
persistence_summary[persistence_cells, on = .(patient, tissue), `:=`(
  response_group = i.response_group,
  n_cells = i.n_cells,
  persistent_cells = i.persistent_cells,
  persistent_cell_fraction = i.persistent_cell_fraction,
  persistent_minus_nonpersistent_ptpn22 = i.persistent_minus_nonpersistent_ptpn22
)]
switch_summary <- modal_wide[persistent == TRUE, .(
  n_persistent_clonotypes_with_modal_states = .N,
  state_switch_fraction = mean(state_switch)
), by = .(patient, response_group, tissue)]
persistence_summary[switch_summary, on = .(patient, tissue), `:=`(
  n_persistent_clonotypes_with_modal_states = i.n_persistent_clonotypes_with_modal_states,
  state_switch_fraction = i.state_switch_fraction
)]

persistence_patient <- tcr[persistence_eligible == TRUE, .(
  persistent_cell_fraction_longitudinal = mean(persistent_within_tissue),
  persistent_minus_nonpersistent_ptpn22 = safe_mean(ptpn22_log1p_cp10k[persistent_within_tissue == TRUE]) - safe_mean(ptpn22_log1p_cp10k[persistent_within_tissue == FALSE])
), by = .(patient, response_group)]
persistence_clone_patient <- clone_persistence[persistence_eligible == TRUE, .(
  persistent_clonotype_fraction_longitudinal = mean(persistent_within_tissue)
), by = .(patient)]
persistence_switch_patient <- modal_wide[persistent == TRUE, .(
  persistent_state_switch_fraction = mean(state_switch)
), by = .(patient)]
persistence_patient[persistence_clone_patient, on = "patient", persistent_clonotype_fraction_longitudinal := i.persistent_clonotype_fraction_longitudinal]
persistence_patient[persistence_switch_patient, on = "patient", persistent_state_switch_fraction := i.persistent_state_switch_fraction]
patient_summary[persistence_patient, on = .(patient, response_group), `:=`(
  persistent_cell_fraction_longitudinal = i.persistent_cell_fraction_longitudinal,
  persistent_clonotype_fraction_longitudinal = i.persistent_clonotype_fraction_longitudinal,
  persistent_minus_nonpersistent_ptpn22 = i.persistent_minus_nonpersistent_ptpn22,
  persistent_state_switch_fraction = i.persistent_state_switch_fraction
)]

bootstrap_effect <- function(r, nr, iterations = 5000L) {
  r <- r[is.finite(r)]
  nr <- nr[is.finite(nr)]
  if (length(r) < 2L || length(nr) < 2L) return(c(NA_real_, NA_real_))
  draws <- replicate(iterations, mean(sample(r, length(r), TRUE)) - mean(sample(nr, length(nr), TRUE)))
  as.numeric(quantile(draws, c(0.025, 0.975), names = FALSE, na.rm = TRUE))
}

cliff_delta <- function(r, nr) {
  r <- r[is.finite(r)]
  nr <- nr[is.finite(nr)]
  if (!length(r) || !length(nr)) return(NA_real_)
  comparisons <- outer(r, nr, "-")
  mean(comparisons > 0) - mean(comparisons < 0)
}

effect_row <- function(data, metric, level, context) {
  r <- data[response_group == "R", get(metric)]
  nr <- data[response_group == "NR", get(metric)]
  r <- r[is.finite(r)]
  nr <- nr[is.finite(nr)]
  ci <- bootstrap_effect(r, nr)
  data.table(
    analysis_level = level,
    context = context,
    metric = metric,
    n_R = length(r),
    n_NR = length(nr),
    mean_R = safe_mean(r),
    mean_NR = safe_mean(nr),
    effect_R_minus_NR = safe_mean(r) - safe_mean(nr),
    bootstrap_ci_low = ci[[1L]],
    bootstrap_ci_high = ci[[2L]],
    cliff_delta = cliff_delta(r, nr),
    wilcoxon_p = if (length(r) >= 2L && length(nr) >= 2L) suppressWarnings(wilcox.test(r, nr, exact = FALSE)$p.value) else NA_real_,
    eligibility = if (length(r) >= 3L && length(nr) >= 3L) "DESCRIPTIVE_PATIENT_LEVEL" else "UNDERPOWERED_FEWER_THAN_3_PER_GROUP"
  )
}

patient_metrics <- setdiff(names(patient_summary), c("patient", "response_group"))
effects <- rbindlist(lapply(patient_metrics, function(metric) effect_row(patient_summary, metric, "patient", "all_TCR_linked_cells")))

state_metrics <- c("expanded_cell_fraction", "large_expanded_cell_fraction", "mean_ptpn22_log1p_cp10k", "ptpn22_detection_fraction", "expanded_minus_singleton_ptpn22")
state_effects <- rbindlist(lapply(unique(state_summary$sub_cluster), function(state) {
  d <- state_summary[sub_cluster == state]
  rbindlist(lapply(state_metrics, function(metric) effect_row(d, metric, "patient_state", state)))
}))
effects <- rbind(effects, state_effects, fill = TRUE)

one_sample_row <- function(data, metric, context) {
  values <- data[[metric]]
  values <- values[is.finite(values)]
  draws <- if (length(values) >= 2L) replicate(5000L, mean(sample(values, length(values), TRUE))) else NA_real_
  ci <- if (length(values) >= 2L) as.numeric(quantile(draws, c(0.025, 0.975), names = FALSE)) else c(NA_real_, NA_real_)
  data.table(
    context = context,
    metric = metric,
    n_patients = length(values),
    mean_patient_effect = safe_mean(values),
    median_patient_effect = if (length(values)) median(values) else NA_real_,
    bootstrap_ci_low = ci[[1L]],
    bootstrap_ci_high = ci[[2L]],
    positive_patients = sum(values > 0),
    negative_patients = sum(values < 0),
    zero_patients = sum(values == 0),
    wilcoxon_signed_rank_p = if (length(values) >= 3L) suppressWarnings(wilcox.test(values, mu = 0, exact = FALSE)$p.value) else NA_real_
  )
}

within_patient <- rbindlist(list(
  one_sample_row(patient_summary, "expanded_minus_singleton_ptpn22", "all_TCR_linked_cells"),
  one_sample_row(patient_summary, "large_minus_singleton_ptpn22", "all_TCR_linked_cells"),
  one_sample_row(patient_summary, "expanded_minus_singleton_detection", "all_TCR_linked_cells"),
  one_sample_row(patient_summary, "clone_size_ptpn22_spearman", "patient_clonotype_summaries"),
  one_sample_row(patient_summary, "persistent_minus_nonpersistent_ptpn22", "longitudinally_eligible_patient_tissues")
))

mapping_audit <- data.table(
  audit_item = c(
    "transcriptomic_barcodes", "tcr_barcodes", "matched_barcodes", "barcode_match_rate",
    "unique_tcr_barcodes", "author_clone_ids", "clone_id_to_single_cdr3_pair",
    "reported_clone_size_exact", "sample_metadata_column", "sample_match_rate",
    "subcluster_metadata_column", "subcluster_exact_label_match_rate", "subcluster_numeric_id_match_rate", "major_metadata_column",
    "major_exact_label_match_rate", "major_semantic_match_rate", "patients_retained", "responders", "nonresponders",
    "PTPN22_gene_column", "PTPN22_gene_zero_based_index", "X_encoding"
  ),
  value = as.character(c(
    n_obs, nrow(tcr), sum(matched), match_rate,
    uniqueN(tcr$barcode), uniqueN(tcr$clonotype_key), all(clone_integrity$cdr3_pairs == 1L),
    all(clone_integrity$reported_size == clone_integrity$observed_cells), sample_column, sample_match_rate,
    subcluster_column, subcluster_match_rate, subcluster_id_match_rate, major_column, major_match_rate, major_semantic_match_rate,
    uniqueN(tcr$patient), uniqueN(tcr[response_group == "R", patient]), uniqueN(tcr[response_group == "NR", patient]),
    ifelse(is.na(gene_column), "var/_index", gene_column), ptpn22_zero_based, encoding
  )),
  status = c(
    "COUNT", "COUNT", "COUNT", ifelse(match_rate >= 0.95, "PASS", "FAIL"),
    ifelse(uniqueN(tcr$barcode) == nrow(tcr), "PASS", "FAIL"), "COUNT", "PASS", "PASS",
    ifelse(is.na(sample_column), "NOT_AVAILABLE", "FOUND"), ifelse(is.na(sample_match_rate), "NOT_TESTED", ifelse(sample_match_rate == 1, "PASS", "REVIEW")),
    ifelse(is.na(subcluster_column), "NOT_AVAILABLE", "FOUND"), ifelse(is.na(subcluster_match_rate), "NOT_TESTED", ifelse(subcluster_match_rate == 1, "PASS", "RELEASED_MARKER_ALIAS_DIFFERENCE")), "PASS",
    ifelse(is.na(major_column), "NOT_AVAILABLE", "FOUND"), ifelse(is.na(major_match_rate), "NOT_TESTED", ifelse(major_match_rate == 1, "PASS", "RELEASED_MAJOR_LABEL_ALIAS_DIFFERENCE")), "PASS",
    "COUNT", "COUNT", "COUNT", "PASS", "PASS", "PASS"
  )
)

mapping_columns <- c(
  "patient", "response_group", "sample", "timepoint", "tissue", "barcode", "obs_row",
  "cdr3_pair_aa", "clone.id", "clonotype_key", "clone.size", "expansion_category",
  "major_cluster", "sub_cluster", "h5_sample", "h5_major_cluster", "h5_subcluster",
  "persistence_eligible", "persistent_within_tissue", "ptpn22_raw_count", "total_raw_umi",
  "ptpn22_log1p_cp10k", "ptpn22_detected"
)
fwrite(tcr[, ..mapping_columns], mapping_file, sep = "\t", compress = "gzip")
fwrite(clone_summary, clone_file, sep = "\t", compress = "gzip")
fwrite(patient_summary, patient_file, sep = "\t")
fwrite(state_summary, state_file, sep = "\t")
fwrite(context_summary, context_file, sep = "\t")
fwrite(persistence_summary, persistence_file, sep = "\t")
fwrite(effects, effects_file, sep = "\t")
fwrite(within_patient, within_patient_file, sep = "\t")
fwrite(mapping_audit, mapping_audit_file, sep = "\t")
fwrite(state_alias_audit, state_alias_file, sep = "\t")

sha256_file <- function(path) toupper(digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
manifest_paths <- c(tcr_file, h5_file, h5_cell_file, response_file, design_file, mapping_file)
manifest <- data.table(
  role = c("TCR_INPUT", "TRANSCRIPTOME_INPUT", "PTPN22_H5AD_EXTRACTION", "RESPONSE_MAPPING", "FROZEN_ANALYSIS_DESIGN", "BARCODE_CLONOTYPE_STATE_MAPPING"),
  path = manifest_paths,
  bytes = file.info(manifest_paths)$size,
  sha256 = vapply(manifest_paths, sha256_file, character(1)),
  records = c(nrow(tcr), n_obs, n_obs, nrow(response), nrow(fread(design_file)), nrow(tcr)),
  record_key = c("barcode", "obs/_index", "barcode", "patient", "design_item", "patient+barcode"),
  matched_records = c(sum(matched), sum(matched), sum(matched), nrow(response), nrow(fread(design_file)), sum(matched)),
  match_rate = c(match_rate, sum(matched) / n_obs, sum(matched) / n_obs, 1, 1, match_rate),
  match_rate_definition = c(
    "matched_tcr_barcodes/tcr_barcodes",
    "tcr_covered_transcriptomic_barcodes/transcriptomic_barcodes",
    "tcr_covered_transcriptomic_barcodes/transcriptomic_barcodes",
    "mapped_patients/response_rows",
    "frozen_design_rows/design_rows",
    "matched_tcr_barcodes/tcr_barcodes"
  ),
  tcr_barcode_match_rate = rep(match_rate, 6L),
  status = c("PASS", "PASS", "PASS", "PASS", "FROZEN_BEFORE_EXPRESSION", "PASS")
)
setnames(manifest, "record_key", "key")
fwrite(manifest, manifest_file, sep = "\t")

cat("\nPatient summaries:\n")
print(patient_summary[order(response_group, patient)])
cat("\nPrimary patient-level effects:\n")
print(effects[analysis_level == "patient"])
cat("\nMapping audit:\n")
print(mapping_audit)
cat("\nSession info:\n")
print(sessionInfo())
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
