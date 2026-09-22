options(stringsAsFactors = FALSE, scipen = 999)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

set.seed(20260828)

mapping_file <- "03_processed_data/phase3/GSE235863_BARCODE_CLONOTYPE_STATE_PTPN22.tsv.gz"
soft_file <- "02_raw_data/GSE235863/GSE235863_family.soft.gz"
response_file <- "01_metadata/phase3/GSE235863_SCRNA_SCTCR_RESPONSE_MAPPING.tsv"
design_file <- "01_metadata/phase3/GSE235863_CLONOTYPE_ANALYSIS_DESIGN.tsv"
manifest_file <- "01_metadata/phase_temporal/GSE235863_LONGITUDINAL_CLONE_FATE_MANIFEST.tsv"
clone_output_file <- "05_results/tables/temporal/GSE235863_BASELINE_CLONE_FATE_TABLE.tsv"
patient_output_file <- "05_results/tables/temporal/PTPN22_PATIENT_LEVEL_TEMPORAL_EFFECTS.tsv"
expansion_output_file <- "05_results/tables/temporal/PTPN22_FUTURE_EXPANSION_ASSOCIATIONS.tsv"
sensitivity_output_file <- "05_results/tables/temporal/PTPN22_TEMPORAL_SENSITIVITY.tsv"
log_file <- "09_logs/temporal/01_gse235863_baseline_future_clone_fate.log"

for (path in c(manifest_file, clone_output_file, patient_output_file,
               expansion_output_file, sensitivity_output_file, log_file)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}

log_con <- file(log_file, open = "wt", encoding = "UTF-8")
sink(log_con, split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

cat("GSE235863 baseline-to-future clone-fate analysis\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Seed: 20260828\n")
cat("Primary tissue: blood (released tissue code P)\n")
cat("Secondary tissue sensitivity: liver tumor (released tissue code T)\n")
cat("Primary follow-up rule: the single explicitly titled Post-treatment sample within the same patient and tissue\n")
cat("Common follow-up downsampling depth: minimum primary-blood follow-up TCR-linked depth\n")

stopifnot(file.exists(mapping_file), file.exists(soft_file),
          file.exists(response_file), file.exists(design_file))

x <- fread(mapping_file)
response <- fread(response_file)
design <- fread(design_file)

required_columns <- c(
  "patient", "response_group", "sample", "timepoint", "tissue", "barcode",
  "cdr3_pair_aa", "clone.id", "clonotype_key", "major_cluster", "sub_cluster",
  "ptpn22_raw_count", "total_raw_umi", "ptpn22_log1p_cp10k", "ptpn22_detected"
)
stopifnot(all(required_columns %chin% names(x)))
stopifnot(nrow(x) == 58872L, uniqueN(x$barcode) == nrow(x))
stopifnot(all(x$timepoint %chin% c("pre", "post")))
recomputed_ptpn22 <- log1p(10000 * x$ptpn22_raw_count / x$total_raw_umi)
stopifnot(max(abs(x$ptpn22_log1p_cp10k - recomputed_ptpn22)) < 1e-12)

identity_rule <- design[design_item == "clonotype_identity", prespecified_rule]
expression_rule <- design[design_item == "PTPN22_expression", prespecified_rule]
stopifnot(length(identity_rule) == 1L, length(expression_rule) == 1L)
stopifnot(grepl("patient plus released author clone.id", identity_rule, fixed = TRUE))
stopifnot(grepl("log1p(10000 times raw PTPN22 UMI", expression_rule, fixed = TRUE))

# Exact author states are preserved. These label families are prespecified sensitivity
# summaries only and do not replace or relabel the author sub_cluster field.
proliferating_states <- c("CD4_C09_MKI67", "CD8_C11_MKI67")
effector_like_states <- c(
  "CD4_C04_GZMK", "CD4_C05_GNLY", "CD8_C03_CX3CR1",
  "CD8_C07_KIR2DL4", "CD8_C08_GZMK"
)
exhausted_context_states <- c("CD4_C06_CXCL13", "CD4_C07_CTLA4", "CD8_C09_LAYN")

modal_lexical <- function(values) {
  counts <- table(values)
  sort(names(counts)[counts == max(counts)])[[1L]]
}

safe_mean <- function(values) {
  values <- values[is.finite(values)]
  if (length(values)) mean(values) else NA_real_
}

safe_median <- function(values) {
  values <- values[is.finite(values)]
  if (length(values)) median(values) else NA_real_
}

safe_spearman <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < 3L || uniqueN(a[ok]) < 2L || uniqueN(b[ok]) < 2L) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok], method = "spearman"))
}

bootstrap_stat_ci <- function(values, statistic = median, iterations = 10000L) {
  values <- values[is.finite(values)]
  if (length(values) < 2L) return(c(NA_real_, NA_real_))
  draws <- replicate(iterations, statistic(sample(values, length(values), replace = TRUE)))
  as.numeric(quantile(draws, c(0.025, 0.975), names = FALSE, na.rm = TRUE))
}

parse_soft_samples <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE, encoding = "UTF-8")
  starts <- grep("^\\^SAMPLE = ", lines)
  ends <- c(starts[-1L] - 1L, length(lines))
  rows <- lapply(seq_along(starts), function(index) {
    block <- lines[starts[[index]]:ends[[index]]]
    first_value <- function(prefix) {
      hit <- block[startsWith(block, prefix)]
      if (!length(hit)) return(NA_character_)
      sub(prefix, "", hit[[1L]], fixed = TRUE)
    }
    characteristics <- block[startsWith(block, "!Sample_characteristics_ch1 = ")]
    characteristic <- function(key) {
      prefix <- paste0("!Sample_characteristics_ch1 = ", key, ": ")
      hit <- characteristics[startsWith(characteristics, prefix)]
      if (!length(hit)) return(NA_character_)
      sub(prefix, "", hit[[1L]], fixed = TRUE)
    }
    data.table(
      gsm = first_value("^SAMPLE = "),
      title = first_value("!Sample_title = "),
      source_name = first_value("!Sample_source_name_ch1 = "),
      characteristic_patient = characteristic("patient"),
      characteristic_tissue = characteristic("tissue"),
      characteristic_response = characteristic("response"),
      characteristic_library_type = characteristic("library type"),
      characteristic_treatment = characteristic("treatment")
    )
  })
  rbindlist(rows, fill = TRUE)
}

soft <- parse_soft_samples(soft_file)

context <- x[, .(
  has_baseline = any(timepoint == "pre"),
  has_followup = any(timepoint == "post"),
  baseline_sample_count = uniqueN(sample[timepoint == "pre"]),
  followup_sample_count = uniqueN(sample[timepoint == "post"]),
  baseline_sample_id = paste(sort(unique(sample[timepoint == "pre"])), collapse = "|"),
  followup_sample_id = paste(sort(unique(sample[timepoint == "post"])), collapse = "|"),
  baseline_tcr_linked_cells = sum(timepoint == "pre"),
  followup_tcr_linked_cells = sum(timepoint == "post"),
  baseline_clonotypes = uniqueN(clonotype_key[timepoint == "pre"]),
  followup_clonotypes = uniqueN(clonotype_key[timepoint == "post"])
), by = .(patient, response_group, tissue)]

eligible_context <- context[
  has_baseline & has_followup & baseline_sample_count == 1L & followup_sample_count == 1L
]
stopifnot(nrow(eligible_context[tissue == "P"]) == 7L)
stopifnot(nrow(eligible_context[tissue == "T"]) == 5L)
stopifnot(all(eligible_context$patient %chin% response[longitudinal_eligibility != "post-treatment only; excluded from persistence analyses", patient]))

baseline_cells <- x[timepoint == "pre"][
  eligible_context[, .(patient, tissue)], on = .(patient, tissue), nomatch = 0L
]
followup_cells <- x[timepoint == "post"][
  eligible_context[, .(patient, tissue)], on = .(patient, tissue), nomatch = 0L
]

baseline_context_totals <- baseline_cells[, .(
  baseline_tcr_linked_cells = .N
), by = .(patient, response_group, tissue, baseline_sample_id = sample)]

baseline_clones <- baseline_cells[, .(
  clone_id = unique(get("clone.id"))[[1L]],
  cdr3_pair_aa = unique(cdr3_pair_aa)[[1L]],
  baseline_mean_ptpn22_log1p_cp10k = mean(ptpn22_log1p_cp10k),
  baseline_median_ptpn22_log1p_cp10k = median(ptpn22_log1p_cp10k),
  baseline_ptpn22_detection_fraction = mean(ptpn22_detected),
  baseline_ptpn22_raw_umi = sum(ptpn22_raw_count),
  baseline_clone_cells = .N,
  baseline_singleton_flag = .N == 1L,
  baseline_dominant_major_cluster = modal_lexical(major_cluster),
  baseline_dominant_sub_cluster = modal_lexical(sub_cluster),
  baseline_author_state_count = uniqueN(sub_cluster),
  baseline_proliferating_state_fraction = mean(sub_cluster %chin% proliferating_states),
  baseline_effector_like_state_fraction = mean(sub_cluster %chin% effector_like_states),
  baseline_exhausted_context_state_fraction = mean(sub_cluster %chin% exhausted_context_states)
), by = .(patient, response_group, tissue, baseline_sample_id = sample, clonotype_key)]

baseline_clones[baseline_context_totals, on = .(patient, response_group, tissue, baseline_sample_id),
                baseline_tcr_linked_cells := i.baseline_tcr_linked_cells]
baseline_clones[, baseline_clone_fraction := baseline_clone_cells / baseline_tcr_linked_cells]
baseline_clones[, baseline_clone_size_category := fifelse(
  baseline_clone_cells == 1L, "singleton_1",
  fifelse(baseline_clone_cells <= 4L, "small_2_4", "large_ge5")
)]
baseline_clones[, baseline_clone_size_category := factor(
  baseline_clone_size_category,
  levels = c("singleton_1", "small_2_4", "large_ge5")
)]

followup_context_totals <- followup_cells[, .(
  followup_tcr_linked_cells = .N
), by = .(patient, response_group, tissue, followup_sample_id = sample)]
followup_clones <- followup_cells[, .(
  followup_clone_cells = .N
), by = .(patient, response_group, tissue, followup_sample_id = sample, clonotype_key)]
followup_clones[followup_context_totals, on = .(patient, response_group, tissue, followup_sample_id),
                followup_tcr_linked_cells := i.followup_tcr_linked_cells]

clone_fate <- copy(baseline_clones)
clone_fate[followup_clones, on = .(patient, response_group, tissue, clonotype_key), `:=`(
  followup_sample_id = i.followup_sample_id,
  followup_clone_cells = i.followup_clone_cells,
  followup_tcr_linked_cells = i.followup_tcr_linked_cells
)]
clone_fate[followup_context_totals, on = .(patient, response_group, tissue), `:=`(
  followup_sample_id = fifelse(is.na(followup_sample_id), i.followup_sample_id, followup_sample_id),
  followup_tcr_linked_cells = fifelse(is.na(followup_tcr_linked_cells), i.followup_tcr_linked_cells, followup_tcr_linked_cells)
)]
clone_fate[is.na(followup_clone_cells), followup_clone_cells := 0L]
clone_fate[, future_persistent_observed := followup_clone_cells > 0L]
clone_fate[, followup_clone_fraction_observed := followup_clone_cells / followup_tcr_linked_cells]
clone_fate[, abundance_change_all_observed := followup_clone_fraction_observed - baseline_clone_fraction]
clone_fate[, future_abundance_change_persistent := fifelse(
  future_persistent_observed, abundance_change_all_observed, NA_real_
)]
clone_fate[, future_log2_fraction_ratio_persistent := fifelse(
  future_persistent_observed,
  log2(followup_clone_fraction_observed / baseline_clone_fraction),
  NA_real_
)]
clone_fate[, future_fate := fifelse(
  !future_persistent_observed, "observed_disappeared",
  fifelse(future_abundance_change_persistent > 0, "persistent_increased_fraction",
          fifelse(future_abundance_change_persistent < 0, "persistent_decreased_fraction",
                  "persistent_unchanged_fraction"))
)]
clone_fate[, analysis_role := fifelse(tissue == "P", "PRIMARY_BLOOD", "SECONDARY_TUMOR_SENSITIVITY")]
clone_fate[, tissue_source := fifelse(tissue == "P", "blood", "liver tumor")]

stopifnot(all(clone_fate$baseline_clone_cells > 0L))
stopifnot(all(clone_fate$followup_tcr_linked_cells > 0L))
stopifnot(!anyNA(clone_fate$future_persistent_observed))
stopifnot(all(clone_fate[future_persistent_observed == TRUE, followup_clone_cells] > 0L))
stopifnot(all(clone_fate[future_persistent_observed == FALSE, followup_clone_cells] == 0L))

# Parse exact public wording and retain source inconsistencies as explicit ambiguity.
find_soft_record <- function(patient_id, tissue_code, timepoint_value) {
  time_word <- if (timepoint_value == "pre") "Pre-treatment" else "Post-treatment"
  tissue_word <- if (tissue_code == "P") "blood" else "liver tumor"
  prefix <- paste0(time_word, " ", tissue_word, " of patient ", patient_id, " (")
  hits <- soft[startsWith(title, prefix)]
  if (nrow(hits) != 1L) {
    stop(sprintf("Expected one SOFT title for %s/%s/%s; found %d", patient_id, tissue_code, timepoint_value, nrow(hits)))
  }
  hits
}

manifest <- rbindlist(lapply(seq_len(nrow(eligible_context)), function(index) {
  row <- eligible_context[index]
  pre_soft <- find_soft_record(row$patient, row$tissue, "pre")
  post_soft <- find_soft_record(row$patient, row$tissue, "post")
  fate_context <- clone_fate[patient == row$patient & tissue == row$tissue]
  ambiguity_items <- character()
  if (!identical(pre_soft$characteristic_patient, row$patient)) {
    ambiguity_items <- c(ambiguity_items, sprintf(
      "baseline SOFT characteristics patient=%s conflicts with title/released sample patient=%s",
      pre_soft$characteristic_patient, row$patient
    ))
  }
  if (!identical(post_soft$characteristic_patient, row$patient)) {
    ambiguity_items <- c(ambiguity_items, sprintf(
      "follow-up SOFT characteristics patient=%s conflicts with title/released sample patient=%s",
      post_soft$characteristic_patient, row$patient
    ))
  }
  data.table(
    patient_id = row$patient,
    response_group = row$response_group,
    response = response[patient == row$patient, response][[1L]],
    baseline_sample_id = row$baseline_sample_id,
    followup_sample_id = row$followup_sample_id,
    baseline_gsm = pre_soft$gsm,
    followup_gsm = post_soft$gsm,
    exact_baseline_timepoint_wording = pre_soft$title,
    exact_followup_timepoint_wording = post_soft$title,
    tissue_code = row$tissue,
    tissue_source = if (row$tissue == "P") "blood" else "liver tumor",
    analysis_role = if (row$tissue == "P") "PRIMARY_BLOOD" else "SECONDARY_TUMOR_SENSITIVITY",
    baseline_tcr_linked_t_cell_count = row$baseline_tcr_linked_cells,
    followup_tcr_linked_t_cell_count = row$followup_tcr_linked_cells,
    baseline_clonotypes = row$baseline_clonotypes,
    persistent_at_followup = sum(fate_context$future_persistent_observed),
    absent_at_followup = sum(!fate_context$future_persistent_observed),
    primary_followup_rule = "single explicitly titled Post-treatment sample within the same patient and tissue; no sample-number inference",
    mapping_provenance = "frozen 58,872/58,872 exact barcode mapping plus released sample token plus NCBI GEO family SOFT exact title",
    ambiguity_flag = length(ambiguity_items) > 0L,
    ambiguity_detail = if (length(ambiguity_items)) paste(ambiguity_items, collapse = "; ") else "NONE"
  )
}))
setorder(manifest, analysis_role, patient_id)

bootstrap_clone_effect <- function(data, iterations = 5000L) {
  persistent <- data[future_persistent_observed == TRUE, baseline_mean_ptpn22_log1p_cp10k]
  disappeared <- data[future_persistent_observed == FALSE, baseline_mean_ptpn22_log1p_cp10k]
  if (!length(persistent) || !length(disappeared)) {
    return(c(low = NA_real_, high = NA_real_, positive_fraction = NA_real_))
  }
  draws <- replicate(iterations,
    mean(sample(persistent, length(persistent), replace = TRUE)) -
      mean(sample(disappeared, length(disappeared), replace = TRUE))
  )
  interval <- quantile(draws, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
  c(low = interval[[1L]], high = interval[[2L]], positive_fraction = mean(draws > 0))
}

patient_effects <- rbindlist(lapply(split(clone_fate, by = c("patient", "tissue"), keep.by = TRUE), function(data) {
  boot <- bootstrap_clone_effect(data)
  persistent <- data[future_persistent_observed == TRUE]
  disappeared <- data[future_persistent_observed == FALSE]
  data.table(
    patient = data$patient[[1L]],
    response_group = data$response_group[[1L]],
    tissue = data$tissue[[1L]],
    tissue_source = data$tissue_source[[1L]],
    analysis_role = data$analysis_role[[1L]],
    baseline_sample_id = data$baseline_sample_id[[1L]],
    followup_sample_id = data$followup_sample_id[[1L]],
    baseline_tcr_linked_cells = data$baseline_tcr_linked_cells[[1L]],
    followup_tcr_linked_cells = data$followup_tcr_linked_cells[[1L]],
    baseline_clonotypes = nrow(data),
    future_persistent_clonotypes = nrow(persistent),
    future_disappearing_clonotypes = nrow(disappeared),
    persistent_fraction_observed = mean(data$future_persistent_observed),
    mean_baseline_ptpn22_future_persistent = safe_mean(persistent$baseline_mean_ptpn22_log1p_cp10k),
    mean_baseline_ptpn22_future_disappearing = safe_mean(disappeared$baseline_mean_ptpn22_log1p_cp10k),
    mean_difference_persistent_minus_disappearing =
      safe_mean(persistent$baseline_mean_ptpn22_log1p_cp10k) - safe_mean(disappeared$baseline_mean_ptpn22_log1p_cp10k),
    median_baseline_ptpn22_future_persistent = safe_median(persistent$baseline_mean_ptpn22_log1p_cp10k),
    median_baseline_ptpn22_future_disappearing = safe_median(disappeared$baseline_mean_ptpn22_log1p_cp10k),
    median_difference_persistent_minus_disappearing =
      safe_median(persistent$baseline_mean_ptpn22_log1p_cp10k) - safe_median(disappeared$baseline_mean_ptpn22_log1p_cp10k),
    detection_fraction_difference_persistent_minus_disappearing =
      safe_mean(persistent$baseline_ptpn22_detection_fraction) - safe_mean(disappeared$baseline_ptpn22_detection_fraction),
    clone_bootstrap_ci_low = boot[["low"]],
    clone_bootstrap_ci_high = boot[["high"]],
    clone_bootstrap_positive_fraction = boot[["positive_fraction"]],
    effect_direction = fifelse(
      is.na(safe_mean(persistent$baseline_mean_ptpn22_log1p_cp10k)), "NOT_ESTIMABLE",
      fifelse(safe_mean(persistent$baseline_mean_ptpn22_log1p_cp10k) - safe_mean(disappeared$baseline_mean_ptpn22_log1p_cp10k) > 0,
              "POSITIVE", fifelse(safe_mean(persistent$baseline_mean_ptpn22_log1p_cp10k) - safe_mean(disappeared$baseline_mean_ptpn22_log1p_cp10k) < 0,
                                  "NEGATIVE", "ZERO"))
    ),
    inference_unit = "PATIENT; clone bootstrap is within-patient distributional robustness only"
  )
}), fill = TRUE)
setorder(patient_effects, analysis_role, patient)

expansion_associations <- clone_fate[future_persistent_observed == TRUE, .(
  persistent_clonotypes = .N,
  spearman_baseline_ptpn22_vs_fraction_change = safe_spearman(
    baseline_mean_ptpn22_log1p_cp10k, future_abundance_change_persistent
  ),
  spearman_baseline_ptpn22_vs_log2_fraction_ratio = safe_spearman(
    baseline_mean_ptpn22_log1p_cp10k, future_log2_fraction_ratio_persistent
  ),
  non_singleton_persistent_clonotypes = sum(baseline_clone_cells >= 2L),
  spearman_excluding_baseline_singletons = safe_spearman(
    baseline_mean_ptpn22_log1p_cp10k[baseline_clone_cells >= 2L],
    future_abundance_change_persistent[baseline_clone_cells >= 2L]
  )
), by = .(patient, response_group, tissue, tissue_source, analysis_role)]
expansion_associations[, direction := fifelse(
  spearman_baseline_ptpn22_vs_fraction_change > 0, "POSITIVE",
  fifelse(spearman_baseline_ptpn22_vs_fraction_change < 0, "NEGATIVE", "ZERO_OR_NOT_ESTIMABLE")
)]
setorder(expansion_associations, analysis_role, patient)

new_sensitivity_row <- function(
    analysis_family, analysis_id, scope, patient = NA_character_,
    response_group = NA_character_, tissue = NA_character_, n_patients = NA_integer_,
    n_clones = NA_integer_, estimate_name, estimate, ci_low = NA_real_, ci_high = NA_real_,
    secondary_metric = NA_character_, secondary_value = NA_real_, positive_units = NA_integer_,
    negative_units = NA_integer_, p_value_descriptive = NA_real_, interpretation) {
  data.table(
    analysis_family, analysis_id, scope, patient, response_group, tissue,
    n_patients, n_clones, estimate_name, estimate, ci_low, ci_high,
    secondary_metric, secondary_value, positive_units, negative_units,
    p_value_descriptive, interpretation
  )
}

sensitivity_rows <- list()
add_sensitivity <- function(row) {
  sensitivity_rows[[length(sensitivity_rows) + 1L]] <<- row
}

primary_patient <- patient_effects[analysis_role == "PRIMARY_BLOOD" & is.finite(mean_difference_persistent_minus_disappearing)]
primary_effect_values <- primary_patient$mean_difference_persistent_minus_disappearing
primary_ci <- bootstrap_stat_ci(primary_effect_values, median)
primary_positive <- sum(primary_effect_values > 0)
primary_negative <- sum(primary_effect_values < 0)
primary_nonzero <- primary_positive + primary_negative
primary_sign_p <- if (primary_nonzero) binom.test(primary_positive, primary_nonzero, 0.5, alternative = "greater")$p.value else NA_real_
add_sensitivity(new_sensitivity_row(
  "PRIMARY_PATIENT_SUMMARY", "blood_unadjusted", "7 longitudinal blood pairs",
  n_patients = nrow(primary_patient), n_clones = sum(primary_patient$baseline_clonotypes),
  estimate_name = "median_patient_mean_difference", estimate = median(primary_effect_values),
  ci_low = primary_ci[[1L]], ci_high = primary_ci[[2L]],
  secondary_metric = "mean_patient_mean_difference", secondary_value = mean(primary_effect_values),
  positive_units = primary_positive, negative_units = primary_negative,
  p_value_descriptive = primary_sign_p,
  interpretation = "PATIENT_IS_BIOLOGICAL_REPLICATE; across-patient bootstrap and exact sign test are descriptive with n=7"
))

# Frozen baseline-size-category stratification.
size_strata <- clone_fate[, {
  persistent <- baseline_mean_ptpn22_log1p_cp10k[future_persistent_observed == TRUE]
  disappeared <- baseline_mean_ptpn22_log1p_cp10k[future_persistent_observed == FALSE]
  .(
    persistent_clones = length(persistent),
    disappearing_clones = length(disappeared),
    stratum_effect = safe_mean(persistent) - safe_mean(disappeared),
    estimable = length(persistent) > 0L & length(disappeared) > 0L
  )
}, by = .(patient, response_group, tissue, analysis_role, baseline_clone_size_category)]

for (index in seq_len(nrow(size_strata))) {
  row <- size_strata[index]
  add_sensitivity(new_sensitivity_row(
    "BASELINE_CLONE_SIZE_STRATIFICATION", paste0(row$patient, "_", row$tissue, "_", row$baseline_clone_size_category),
    as.character(row$baseline_clone_size_category), patient = row$patient,
    response_group = row$response_group, tissue = row$tissue,
    n_clones = row$persistent_clones + row$disappearing_clones,
    estimate_name = "persistent_minus_disappearing_within_size_stratum",
    estimate = row$stratum_effect,
    secondary_metric = "persistent_clones", secondary_value = row$persistent_clones,
    interpretation = if (row$estimable) "FROZEN_1_2TO4_GE5_SIZE_STRATUM" else "NOT_ESTIMABLE_ONE_FATE_GROUP_ABSENT"
  ))
}

size_adjusted_patient <- size_strata[estimable == TRUE, .(
  size_stratified_weighted_effect = weighted.mean(
    stratum_effect, persistent_clones + disappearing_clones
  ),
  estimable_size_strata = .N,
  represented_clones = sum(persistent_clones + disappearing_clones)
), by = .(patient, response_group, tissue, analysis_role)]

for (index in seq_len(nrow(size_adjusted_patient))) {
  row <- size_adjusted_patient[index]
  add_sensitivity(new_sensitivity_row(
    "BASELINE_CLONE_SIZE_STRATIFICATION", paste0(row$patient, "_", row$tissue, "_weighted"),
    "weighted across estimable frozen size strata", patient = row$patient,
    response_group = row$response_group, tissue = row$tissue,
    n_clones = row$represented_clones,
    estimate_name = "size_stratified_weighted_patient_effect", estimate = row$size_stratified_weighted_effect,
    secondary_metric = "estimable_size_strata", secondary_value = row$estimable_size_strata,
    interpretation = "PATIENT_LEVEL_SIZE_STRATIFIED_SENSITIVITY"
  ))
}

primary_clones <- clone_fate[analysis_role == "PRIMARY_BLOOD"]
primary_clones[, patient_factor := factor(patient)]
primary_clones[, baseline_clone_log_fraction := log1p(10000 * baseline_clone_fraction)]
primary_clones[, persistent_numeric := as.integer(future_persistent_observed)]

size_only_model <- glm(
  persistent_numeric ~ baseline_clone_log_fraction + patient_factor,
  data = primary_clones, family = binomial()
)
size_ptpn22_model <- glm(
  persistent_numeric ~ baseline_mean_ptpn22_log1p_cp10k + baseline_clone_log_fraction + patient_factor,
  data = primary_clones, family = binomial()
)
state_adjusted_model <- glm(
  persistent_numeric ~ baseline_mean_ptpn22_log1p_cp10k + baseline_clone_log_fraction +
    patient_factor + baseline_dominant_major_cluster + baseline_proliferating_state_fraction +
    baseline_effector_like_state_fraction + baseline_exhausted_context_state_fraction,
  data = primary_clones, family = binomial()
)

model_row <- function(model, id, interpretation) {
  coefficients <- summary(model)$coefficients
  coefficient <- coefficients["baseline_mean_ptpn22_log1p_cp10k", ]
  new_sensitivity_row(
    "FIXED_PATIENT_LOGISTIC_SENSITIVITY", id, "primary blood baseline clonotypes",
    n_patients = uniqueN(primary_clones$patient), n_clones = nrow(primary_clones),
    estimate_name = "PTPN22_log_odds_coefficient", estimate = coefficient[["Estimate"]],
    ci_low = coefficient[["Estimate"]] - 1.96 * coefficient[["Std. Error"]],
    ci_high = coefficient[["Estimate"]] + 1.96 * coefficient[["Std. Error"]],
    secondary_metric = "AIC", secondary_value = AIC(model),
    p_value_descriptive = coefficient[["Pr(>|z|)"]], interpretation = interpretation
  )
}
add_sensitivity(model_row(
  size_ptpn22_model, "patient_fixed_effect_plus_clone_fraction_plus_PTPN22",
  "CLONE_LEVEL_MODEL_SENSITIVITY_NOT_INDEPENDENT_PATIENT_REPLICATION"
))
add_sensitivity(model_row(
  state_adjusted_model, "patient_fixed_effect_plus_clone_fraction_plus_state_composition_plus_PTPN22",
  "STATE_FAMILIES_ARE_SENSITIVITY_COVARIATES; EXACT_AUTHOR_STATES_RETAINED"
))

model_comparison <- anova(size_only_model, size_ptpn22_model, test = "Chisq")
add_sensitivity(new_sensitivity_row(
  "FALSIFICATION_MODEL_COMPARISON", "clone_size_only_vs_clone_size_plus_PTPN22",
  "primary blood baseline clonotypes", n_patients = uniqueN(primary_clones$patient),
  n_clones = nrow(primary_clones), estimate_name = "AIC_change_size_plus_PTPN22_minus_size_only",
  estimate = AIC(size_ptpn22_model) - AIC(size_only_model),
  secondary_metric = "deviance_improvement", secondary_value = deviance(size_only_model) - deviance(size_ptpn22_model),
  p_value_descriptive = model_comparison$`Pr(>Chi)`[[2L]],
  interpretation = "NEGATIVE_CONTROL_COMPARISON; clone-level likelihood is descriptive and does not change biological N"
))

# Patient-specific size-adjusted coefficients preserve the patient as the summarization unit.
for (patient_id in sort(unique(primary_clones$patient))) {
  data <- primary_clones[patient == patient_id]
  fit <- try(glm(
    persistent_numeric ~ baseline_mean_ptpn22_log1p_cp10k + baseline_clone_log_fraction,
    data = data, family = binomial()
  ), silent = TRUE)
  if (inherits(fit, "try-error") || !"baseline_mean_ptpn22_log1p_cp10k" %chin% rownames(summary(fit)$coefficients)) {
    estimate <- low <- high <- p_value <- NA_real_
    status <- "MODEL_NOT_ESTIMABLE"
  } else {
    coefficient <- summary(fit)$coefficients["baseline_mean_ptpn22_log1p_cp10k", ]
    estimate <- coefficient[["Estimate"]]
    low <- estimate - 1.96 * coefficient[["Std. Error"]]
    high <- estimate + 1.96 * coefficient[["Std. Error"]]
    p_value <- coefficient[["Pr(>|z|)"]]
    status <- "PATIENT_SPECIFIC_SIZE_ADJUSTED_DIRECTION_SENSITIVITY"
  }
  add_sensitivity(new_sensitivity_row(
    "PATIENT_SPECIFIC_SIZE_MODEL", paste0(patient_id, "_size_adjusted"),
    "primary blood", patient = patient_id,
    response_group = data$response_group[[1L]], tissue = "P", n_clones = nrow(data),
    estimate_name = "PTPN22_log_odds_coefficient", estimate = estimate,
    ci_low = low, ci_high = high, p_value_descriptive = p_value,
    interpretation = status
  ))
}

# Exact author-state and broad-compartment stratification.
state_strata <- primary_clones[, {
  persistent <- baseline_mean_ptpn22_log1p_cp10k[future_persistent_observed == TRUE]
  disappeared <- baseline_mean_ptpn22_log1p_cp10k[future_persistent_observed == FALSE]
  .(
    persistent_clones = length(persistent),
    disappearing_clones = length(disappeared),
    effect = safe_mean(persistent) - safe_mean(disappeared),
    estimable = length(persistent) >= 2L & length(disappeared) >= 2L
  )
}, by = .(patient, response_group, baseline_dominant_sub_cluster)]

state_adjusted_patient <- state_strata[estimable == TRUE, .(
  exact_state_stratified_weighted_effect = weighted.mean(
    effect, persistent_clones + disappearing_clones
  ),
  estimable_exact_states = .N,
  represented_clones = sum(persistent_clones + disappearing_clones)
), by = .(patient, response_group)]

priority_state_family <- data.table(
  baseline_dominant_sub_cluster = c(proliferating_states, effector_like_states, exhausted_context_states),
  priority_family = c(
    rep("proliferating", length(proliferating_states)),
    rep("effector_like", length(effector_like_states)),
    rep("exhausted_context", length(exhausted_context_states))
  )
)
priority_state_rows <- state_strata[
  priority_state_family, on = "baseline_dominant_sub_cluster", nomatch = 0L
]
for (index in seq_len(nrow(priority_state_rows))) {
  row <- priority_state_rows[index]
  add_sensitivity(new_sensitivity_row(
    "PRIORITY_EXACT_AUTHOR_STATE", paste0(row$patient, "_", row$baseline_dominant_sub_cluster),
    paste0(row$priority_family, ":", row$baseline_dominant_sub_cluster),
    patient = row$patient, response_group = row$response_group, tissue = "P",
    n_clones = row$persistent_clones + row$disappearing_clones,
    estimate_name = "persistent_minus_disappearing_within_exact_author_state",
    estimate = row$effect,
    secondary_metric = "persistent_clones", secondary_value = row$persistent_clones,
    interpretation = if (row$estimable) {
      "EXACT_AUTHOR_STATE_ESTIMABLE_AT_LEAST_2_CLONES_PER_FATE_GROUP"
    } else {
      "INSUFFICIENT_WITHIN_STATE_FATE_COUNTS; retained without promotion"
    }
  ))
}

priority_state_summary <- priority_state_rows[estimable == TRUE & is.finite(effect), .(
  n_patients = .N,
  median_patient_effect = median(effect),
  positive_patients = sum(effect > 0),
  negative_patients = sum(effect < 0),
  represented_clones = sum(persistent_clones + disappearing_clones)
), by = .(priority_family, baseline_dominant_sub_cluster)]
for (index in seq_len(nrow(priority_state_summary))) {
  row <- priority_state_summary[index]
  add_sensitivity(new_sensitivity_row(
    "PRIORITY_EXACT_AUTHOR_STATE_SUMMARY", row$baseline_dominant_sub_cluster,
    paste0(row$priority_family, ":", row$baseline_dominant_sub_cluster),
    tissue = "P", n_patients = row$n_patients, n_clones = row$represented_clones,
    estimate_name = "median_estimable_patient_effect", estimate = row$median_patient_effect,
    positive_units = row$positive_patients, negative_units = row$negative_patients,
    interpretation = "DESCRIPTIVE_ONLY; exact author state; sparse states are not promoted"
  ))
}

for (index in seq_len(nrow(state_adjusted_patient))) {
  row <- state_adjusted_patient[index]
  add_sensitivity(new_sensitivity_row(
    "EXACT_AUTHOR_STATE_STRATIFICATION", paste0(row$patient, "_exact_state_weighted"),
    "primary blood exact dominant author sub_cluster", patient = row$patient,
    response_group = row$response_group, tissue = "P", n_clones = row$represented_clones,
    estimate_name = "exact_state_stratified_weighted_patient_effect",
    estimate = row$exact_state_stratified_weighted_effect,
    secondary_metric = "estimable_exact_states", secondary_value = row$estimable_exact_states,
    interpretation = "NO_STATE_RELABELING; states with fewer than 2 clones per future-fate group excluded"
  ))
}

major_strata <- primary_clones[, {
  persistent <- baseline_mean_ptpn22_log1p_cp10k[future_persistent_observed == TRUE]
  disappeared <- baseline_mean_ptpn22_log1p_cp10k[future_persistent_observed == FALSE]
  .(
    persistent_clones = length(persistent), disappearing_clones = length(disappeared),
    effect = safe_mean(persistent) - safe_mean(disappeared),
    estimable = length(persistent) >= 2L & length(disappeared) >= 2L
  )
}, by = .(patient, response_group, baseline_dominant_major_cluster)]
for (index in seq_len(nrow(major_strata))) {
  row <- major_strata[index]
  add_sensitivity(new_sensitivity_row(
    "BROAD_CD4_CD8_STRATIFICATION", paste0(row$patient, "_", row$baseline_dominant_major_cluster),
    row$baseline_dominant_major_cluster, patient = row$patient,
    response_group = row$response_group, tissue = "P",
    n_clones = row$persistent_clones + row$disappearing_clones,
    estimate_name = "persistent_minus_disappearing_within_broad_compartment", estimate = row$effect,
    secondary_metric = "persistent_clones", secondary_value = row$persistent_clones,
    interpretation = if (row$estimable) "AUTHOR_MAJOR_CLUSTER_STRATIFICATION" else "INSUFFICIENT_WITHIN_COMPARTMENT_FATE_COUNTS"
  ))
}

# Response is secondary context only.
response_summary <- primary_patient[, .(
  n_patients = .N,
  median_effect = median(mean_difference_persistent_minus_disappearing),
  mean_effect = mean(mean_difference_persistent_minus_disappearing),
  positive_patients = sum(mean_difference_persistent_minus_disappearing > 0),
  negative_patients = sum(mean_difference_persistent_minus_disappearing < 0)
), by = response_group]
for (index in seq_len(nrow(response_summary))) {
  row <- response_summary[index]
  add_sensitivity(new_sensitivity_row(
    "RESPONSE_SECONDARY_CONTEXT", paste0("response_", row$response_group),
    "primary blood patient effects", response_group = row$response_group,
    n_patients = row$n_patients, estimate_name = "median_patient_effect", estimate = row$median_effect,
    secondary_metric = "mean_patient_effect", secondary_value = row$mean_effect,
    positive_units = row$positive_patients, negative_units = row$negative_patients,
    interpretation = "SECONDARY_DIRECTION_ONLY; no credible differential-prediction claim with 4 R and 3 NR"
  ))
}

# Within-patient permutation, constrained within frozen baseline-size categories.
permutation_iterations <- 5000L
primary_by_patient <- split(primary_clones, by = "patient", keep.by = TRUE)
observed_permutation_statistic <- median(primary_effect_values)
permuted_statistics <- numeric(permutation_iterations)
permuted_positive_patients <- integer(permutation_iterations)
for (iteration in seq_len(permutation_iterations)) {
  patient_draws <- vapply(primary_by_patient, function(data) {
    permuted <- data$baseline_mean_ptpn22_log1p_cp10k
    for (category in levels(data$baseline_clone_size_category)) {
      indices <- which(data$baseline_clone_size_category == category)
      if (length(indices) > 1L) permuted[indices] <- sample(permuted[indices], length(indices), replace = FALSE)
    }
    safe_mean(permuted[data$future_persistent_observed]) - safe_mean(permuted[!data$future_persistent_observed])
  }, numeric(1))
  permuted_statistics[[iteration]] <- median(patient_draws)
  permuted_positive_patients[[iteration]] <- sum(patient_draws > 0)
}
permutation_p <- (1 + sum(permuted_statistics >= observed_permutation_statistic)) / (permutation_iterations + 1)
add_sensitivity(new_sensitivity_row(
  "WITHIN_PATIENT_PERMUTATION", "PTPN22_permuted_within_patient_and_size_category",
  "primary blood; median patient effect statistic", n_patients = nrow(primary_patient),
  n_clones = nrow(primary_clones), estimate_name = "observed_median_patient_effect",
  estimate = observed_permutation_statistic,
  ci_low = quantile(permuted_statistics, 0.025), ci_high = quantile(permuted_statistics, 0.975),
  secondary_metric = "null_mean_median_effect", secondary_value = mean(permuted_statistics),
  positive_units = primary_positive, p_value_descriptive = permutation_p,
  interpretation = "ONE_SIDED_EMPIRICAL_P; null interval shown; permutation preserves patient and frozen clone-size category"
))

# Prespecified depth and baseline-clone-size sensitivities.
adequate_patients <- manifest[
  analysis_role == "PRIMARY_BLOOD" &
    baseline_tcr_linked_t_cell_count >= 1000L & followup_tcr_linked_t_cell_count >= 1000L,
  patient_id
]
adequate_effects <- primary_patient[patient %chin% adequate_patients, mean_difference_persistent_minus_disappearing]
adequate_ci <- bootstrap_stat_ci(adequate_effects, median)
add_sensitivity(new_sensitivity_row(
  "SAMPLING_DEPTH_RESTRICTION", "both_timepoints_ge1000_TCR_linked_cells",
  "primary blood", n_patients = length(adequate_effects),
  n_clones = sum(primary_patient[patient %chin% adequate_patients, baseline_clonotypes]),
  estimate_name = "median_patient_effect", estimate = median(adequate_effects),
  ci_low = adequate_ci[[1L]], ci_high = adequate_ci[[2L]],
  positive_units = sum(adequate_effects > 0), negative_units = sum(adequate_effects < 0),
  interpretation = "PRESPECIFIED_ADEQUATE_DEPTH_THRESHOLD; observational persistence remains detection-limited"
))

non_singleton_effects <- primary_clones[baseline_clone_cells >= 2L, .(
  n_clones = .N,
  persistent_clones = sum(future_persistent_observed),
  effect = safe_mean(baseline_mean_ptpn22_log1p_cp10k[future_persistent_observed]) -
    safe_mean(baseline_mean_ptpn22_log1p_cp10k[!future_persistent_observed])
), by = .(patient, response_group)]
non_singleton_values <- non_singleton_effects[is.finite(effect), effect]
non_singleton_ci <- bootstrap_stat_ci(non_singleton_values, median)
add_sensitivity(new_sensitivity_row(
  "BASELINE_MINIMAL_ABUNDANCE", "baseline_clone_cells_ge2",
  "primary blood", n_patients = length(non_singleton_values),
  n_clones = sum(non_singleton_effects$n_clones), estimate_name = "median_patient_effect",
  estimate = median(non_singleton_values), ci_low = non_singleton_ci[[1L]], ci_high = non_singleton_ci[[2L]],
  positive_units = sum(non_singleton_values > 0), negative_units = sum(non_singleton_values < 0),
  interpretation = "SINGLETONS_EXCLUDED_ONLY_AS_SENSITIVITY; primary analysis retains them"
))

tumor_effects <- patient_effects[
  analysis_role == "SECONDARY_TUMOR_SENSITIVITY" & is.finite(mean_difference_persistent_minus_disappearing),
  mean_difference_persistent_minus_disappearing
]
tumor_ci <- bootstrap_stat_ci(tumor_effects, median)
add_sensitivity(new_sensitivity_row(
  "TISSUE_SENSITIVITY", "complete_pre_post_liver_tumor_pairs",
  "secondary tumor contexts; same patients are not added as independent replication",
  n_patients = length(tumor_effects),
  n_clones = sum(patient_effects[analysis_role == "SECONDARY_TUMOR_SENSITIVITY", baseline_clonotypes]),
  estimate_name = "median_patient_effect", estimate = median(tumor_effects),
  ci_low = tumor_ci[[1L]], ci_high = tumor_ci[[2L]],
  positive_units = sum(tumor_effects > 0), negative_units = sum(tumor_effects < 0),
  interpretation = "TUMOR_IS_SECONDARY_TISSUE_SENSITIVITY_NOT_EXTRA_BIOLOGICAL_N"
))

# Equalize all primary follow-up blood samples to the minimum observed follow-up depth.
common_followup_depth <- min(manifest[analysis_role == "PRIMARY_BLOOD", followup_tcr_linked_t_cell_count])
downsample_iterations <- 1000L
primary_post_by_patient <- split(followup_cells[tissue == "P"], by = "patient", keep.by = TRUE)
primary_baseline_by_patient <- split(primary_clones, by = "patient", keep.by = TRUE)
downsample_matrix <- matrix(
  NA_real_, nrow = downsample_iterations, ncol = length(primary_baseline_by_patient),
  dimnames = list(NULL, names(primary_baseline_by_patient))
)
for (iteration in seq_len(downsample_iterations)) {
  for (patient_id in names(primary_baseline_by_patient)) {
    post_data <- primary_post_by_patient[[patient_id]]
    baseline_data <- primary_baseline_by_patient[[patient_id]]
    sampled_indices <- if (nrow(post_data) > common_followup_depth) {
      sample.int(nrow(post_data), common_followup_depth, replace = FALSE)
    } else {
      seq_len(nrow(post_data))
    }
    detected <- unique(post_data$clonotype_key[sampled_indices])
    persistent_sampled <- baseline_data$clonotype_key %chin% detected
    downsample_matrix[iteration, patient_id] <-
      safe_mean(baseline_data$baseline_mean_ptpn22_log1p_cp10k[persistent_sampled]) -
      safe_mean(baseline_data$baseline_mean_ptpn22_log1p_cp10k[!persistent_sampled])
  }
}
for (patient_id in colnames(downsample_matrix)) {
  draws <- downsample_matrix[, patient_id]
  data <- primary_clones[patient == patient_id]
  add_sensitivity(new_sensitivity_row(
    "FOLLOWUP_TECHNICAL_DOWNSAMPLING", paste0(patient_id, "_followup_depth_", common_followup_depth),
    "primary blood", patient = patient_id, response_group = data$response_group[[1L]], tissue = "P",
    n_clones = nrow(data), estimate_name = "median_downsampled_patient_effect", estimate = median(draws, na.rm = TRUE),
    ci_low = quantile(draws, 0.025, na.rm = TRUE), ci_high = quantile(draws, 0.975, na.rm = TRUE),
    secondary_metric = "positive_downsample_fraction", secondary_value = mean(draws > 0, na.rm = TRUE),
    interpretation = "1000_WITHIN_FOLLOWUP_TECHNICAL_DOWNSAMPLES; does not create biological replication"
  ))
}
downsample_overall <- apply(downsample_matrix, 1L, median, na.rm = TRUE)
add_sensitivity(new_sensitivity_row(
  "FOLLOWUP_TECHNICAL_DOWNSAMPLING", paste0("overall_followup_depth_", common_followup_depth),
  "median across 7 primary-blood patient effects per technical draw",
  n_patients = ncol(downsample_matrix), n_clones = nrow(primary_clones),
  estimate_name = "median_of_downsampled_median_patient_effects", estimate = median(downsample_overall),
  ci_low = quantile(downsample_overall, 0.025), ci_high = quantile(downsample_overall, 0.975),
  secondary_metric = "positive_overall_draw_fraction", secondary_value = mean(downsample_overall > 0),
  interpretation = "COMMON_FOLLOWUP_DEPTH_580; technical sampling sensitivity only"
))

# Across-patient expansion summaries.
primary_expansion <- expansion_associations[analysis_role == "PRIMARY_BLOOD"]
for (metric in c("spearman_baseline_ptpn22_vs_fraction_change", "spearman_excluding_baseline_singletons")) {
  values <- primary_expansion[[metric]]
  values <- values[is.finite(values)]
  interval <- bootstrap_stat_ci(values, median)
  add_sensitivity(new_sensitivity_row(
    "FUTURE_ABUNDANCE_CHANGE_SUMMARY", metric,
    "persistent primary-blood clonotypes; patient-level Spearman rho",
    n_patients = length(values), n_clones = sum(primary_expansion$persistent_clonotypes),
    estimate_name = "median_patient_rho", estimate = median(values),
    ci_low = interval[[1L]], ci_high = interval[[2L]],
    positive_units = sum(values > 0), negative_units = sum(values < 0),
    interpretation = "PATIENT_LEVEL_CORRELATION_SUMMARY; no pooled unstratified clone correlation"
  ))
}

sensitivity <- rbindlist(sensitivity_rows, fill = TRUE)

clone_output_columns <- c(
  "patient", "response_group", "tissue", "tissue_source", "analysis_role",
  "baseline_sample_id", "followup_sample_id", "clone_id", "clonotype_key", "cdr3_pair_aa",
  "baseline_mean_ptpn22_log1p_cp10k", "baseline_median_ptpn22_log1p_cp10k",
  "baseline_ptpn22_detection_fraction", "baseline_ptpn22_raw_umi",
  "baseline_clone_cells", "baseline_tcr_linked_cells", "baseline_clone_fraction",
  "baseline_clone_size_category", "baseline_singleton_flag",
  "baseline_dominant_major_cluster", "baseline_dominant_sub_cluster", "baseline_author_state_count",
  "baseline_proliferating_state_fraction", "baseline_effector_like_state_fraction",
  "baseline_exhausted_context_state_fraction", "followup_clone_cells",
  "followup_tcr_linked_cells", "followup_clone_fraction_observed", "future_persistent_observed",
  "abundance_change_all_observed", "future_abundance_change_persistent",
  "future_log2_fraction_ratio_persistent", "future_fate"
)

fwrite(manifest, manifest_file, sep = "\t")
fwrite(clone_fate[, ..clone_output_columns], clone_output_file, sep = "\t")
fwrite(patient_effects, patient_output_file, sep = "\t")
fwrite(expansion_associations, expansion_output_file, sep = "\t")
fwrite(sensitivity, sensitivity_output_file, sep = "\t")

sha256_file <- function(path) {
  digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

cat("\nInput integrity\n")
cat("Mapping rows:", nrow(x), "\n")
cat("Mapping SHA256:", sha256_file(mapping_file), "\n")
cat("SOFT SHA256:", sha256_file(soft_file), "\n")
cat("Frozen identity rule:", identity_rule, "\n")
cat("Frozen expression rule:", expression_rule, "\n")
cat("\nEligible contexts\n")
print(manifest)
cat("\nPrimary patient effects\n")
print(patient_effects[analysis_role == "PRIMARY_BLOOD"])
cat("\nPrimary expansion associations\n")
print(primary_expansion)
cat("\nKey sensitivity rows\n")
print(sensitivity[analysis_family %chin% c(
  "PRIMARY_PATIENT_SUMMARY", "FIXED_PATIENT_LOGISTIC_SENSITIVITY",
  "FALSIFICATION_MODEL_COMPARISON", "WITHIN_PATIENT_PERMUTATION",
  "SAMPLING_DEPTH_RESTRICTION", "BASELINE_MINIMAL_ABUNDANCE",
  "TISSUE_SENSITIVITY", "FOLLOWUP_TECHNICAL_DOWNSAMPLING",
  "FUTURE_ABUNDANCE_CHANGE_SUMMARY"
)])
cat("\nOutput hashes\n")
for (path in c(manifest_file, clone_output_file, patient_output_file,
               expansion_output_file, sensitivity_output_file)) {
  cat(path, sha256_file(path), "\n")
}
cat("\nSession info\n")
print(sessionInfo())
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
