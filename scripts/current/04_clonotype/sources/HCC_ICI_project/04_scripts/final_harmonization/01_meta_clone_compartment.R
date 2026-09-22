options(stringsAsFactors = FALSE, width = 260)
invisible(try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE))

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

set.seed(20260829L)

project_dir <- normalizePath("HCC_ICI_project", winslash = "/", mustWork = TRUE)
metadata_dir <- file.path(project_dir, "01_metadata", "final_harmonization")
table_dir <- file.path(project_dir, "05_results", "tables", "final_harmonization")
log_dir <- file.path(project_dir, "09_logs", "final_harmonization")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

role_file <- file.path(metadata_dir, "HCC_FINAL_COHORT_ROLE_TABLE.tsv")
estimand_file <- file.path(metadata_dir, "HCC_PRIMARY_ICI_ESTIMAND_CONTRACT.md")
effect_file <- file.path(project_dir, "05_results", "tables", "clinical_translation", "PTPN22_CONTINUOUS_AXIS_EFFECTS.tsv")
mapping_file <- file.path(project_dir, "03_processed_data", "phase3", "GSE235863_BARCODE_CLONOTYPE_STATE_PTPN22.tsv.gz")
v6_effect_file <- file.path(project_dir, "05_results", "tables", "temporal", "PTPN22_PATIENT_LEVEL_TEMPORAL_EFFECTS.tsv")
soft_file <- file.path(project_dir, "02_raw_data", "GSE235863", "GSE235863_family.soft.gz")

required_files <- c(role_file, estimand_file, effect_file, mapping_file, v6_effect_file, soft_file)
stopifnot(all(file.exists(required_files)))
stopifnot(
  digest(file = role_file, algo = "sha256") == "91721aa3aaef8cb6b25ef2f57e43f211ed8071d506f4069743f5984d4dfc4778",
  digest(file = estimand_file, algo = "sha256") == "754384aa3d7eaa0dd5c5b95fbd49b299833b6c5940303d7777bb0e3835d5b995"
)

log_file <- file.path(log_dir, "01_meta_clone_compartment.log")
log_con <- file(log_file, open = "wt", encoding = "UTF-8")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("FINAL STATISTICAL HARMONIZATION: META + CLONE ID + COMPARTMENT\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("Seed: 20260829\n")
cat("Role table SHA256:", digest(file = role_file, algo = "sha256"), "\n")
cat("Estimand contract SHA256:", digest(file = estimand_file, algo = "sha256"), "\n")

# -----------------------------------------------------------------------------
# Final ICI meta-analysis. All yi/vi values come from the frozen Phase 9 effects.
# -----------------------------------------------------------------------------

roles <- fread(role_file)
effects <- fread(effect_file)
effects <- effects[axis_or_metric == "PTPN22_anchor" & is.finite(hedges_g)]
effects[, vi := ((hedges_g_ci_high - hedges_g_ci_low) / (2 * 1.96))^2]
stopifnot(nrow(effects) == 6L, !anyDuplicated(effects$dataset), all(effects$vi > 0))

design <- data.table(
  dataset = c("GSE206325", "GSE202069", "GSE238264", "GSE215011", "GSE235863", "GSE279750"),
  timing = c("post_treatment", "unknown", "post_treatment", "unknown", "post_treatment", "unknown"),
  modality = c("state_specific_scRNA_pseudobulk", "bulk_whole_tumor", "spatial_whole_slide_pseudobulk", "bulk_whole_tumor", "bulk_whole_tumor", "bulk_whole_tumor"),
  endpoint_type = c("pathologic_or_clinical", "RECIST", "pathologic", "GEO_clinical", "RECIST_CRPR_vs_SDPD", "GEO_clinical"),
  therapy_class = c("monotherapy", "monotherapy", "combination", "monotherapy", "combination", "unknown"),
  target_selection_involved = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  primary_compatible = c(FALSE, TRUE, FALSE, TRUE, TRUE, TRUE)
)
effects <- merge(effects, design, by = "dataset", all.x = TRUE, sort = FALSE)
stopifnot(!anyNA(effects$timing))

tau2_reml <- function(yi, vi) {
  if (length(yi) < 2L) return(NA_real_)
  nll <- function(tau2) {
    w <- 1 / (vi + tau2)
    mu <- sum(w * yi) / sum(w)
    0.5 * (sum(log(vi + tau2)) + log(sum(w)) + sum(w * (yi - mu)^2))
  }
  upper <- max(10, stats::var(yi) * 100, max(vi) * 100)
  fit <- optimize(nll, interval = c(0, upper), tol = 1e-12)
  if (nll(0) <= fit$objective + 1e-10) 0 else fit$minimum
}

tau2_dl <- function(yi, vi) {
  w <- 1 / vi
  mu <- sum(w * yi) / sum(w)
  q <- sum(w * (yi - mu)^2)
  c_value <- sum(w) - sum(w^2) / sum(w)
  max(0, (q - (length(yi) - 1L)) / c_value)
}

not_estimable_row <- function(model_id, family, method, cohort_note, reason) {
  data.table(
    model_id = model_id, analysis_family = family, method = method,
    status = "NOT_ESTIMABLE", k = 0L, included_cohorts = cohort_note,
    omitted_cohort = "NONE", pooled_hedges_g = NA_real_, standard_error = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_, tau2 = NA_real_,
    Q = NA_real_, Q_df = NA_integer_, Q_p = NA_real_, I2_percent = NA_real_,
    hksj_scale = NA_real_, prediction_interval_low = NA_real_, prediction_interval_high = NA_real_,
    positive_cohorts = NA_integer_, negative_cohorts = NA_integer_, reason = reason
  )
}

fit_meta <- function(d, model_id, family, method = c("REML_HKSJ", "DL_NORMAL"), omitted = "NONE", reason = "") {
  method <- match.arg(method)
  d <- d[is.finite(hedges_g) & is.finite(vi)]
  if (nrow(d) < 2L) {
    return(not_estimable_row(model_id, family, method, paste(d$dataset, collapse = ";"),
      if (nzchar(reason)) reason else "Fewer than two eligible cohorts"))
  }
  yi <- d$hedges_g
  vi <- d$vi
  k <- length(yi)
  fixed_w <- 1 / vi
  fixed_mu <- sum(fixed_w * yi) / sum(fixed_w)
  Q <- sum(fixed_w * (yi - fixed_mu)^2)
  Q_df <- k - 1L
  Q_p <- pchisq(Q, df = Q_df, lower.tail = FALSE)
  I2 <- if (Q > 0) max(0, (Q - Q_df) / Q) * 100 else 0
  tau2 <- if (method == "REML_HKSJ") tau2_reml(yi, vi) else tau2_dl(yi, vi)
  w <- 1 / (vi + tau2)
  pooled <- sum(w * yi) / sum(w)
  if (method == "REML_HKSJ") {
    hk_scale <- sum(w * (yi - pooled)^2) / (k - 1L)
    se <- sqrt(hk_scale / sum(w))
    critical <- qt(0.975, df = k - 1L)
    p_value <- if (se > 0) 2 * pt(-abs(pooled / se), df = k - 1L) else 0
  } else {
    hk_scale <- NA_real_
    se <- sqrt(1 / sum(w))
    critical <- qnorm(0.975)
    p_value <- 2 * pnorm(-abs(pooled / se))
  }
  pred_critical <- if (k >= 3L) qt(0.975, df = k - 2L) else NA_real_
  pred_half <- if (is.finite(pred_critical)) pred_critical * sqrt(tau2 + se^2) else NA_real_
  data.table(
    model_id = model_id, analysis_family = family, method = method,
    status = "ESTIMATED", k = k, included_cohorts = paste(d$dataset, collapse = ";"),
    omitted_cohort = omitted, pooled_hedges_g = pooled, standard_error = se,
    ci_low = pooled - critical * se, ci_high = pooled + critical * se, p_value = p_value,
    tau2 = tau2, Q = Q, Q_df = Q_df, Q_p = Q_p, I2_percent = I2,
    hksj_scale = hk_scale,
    prediction_interval_low = if (is.finite(pred_half)) pooled - pred_half else NA_real_,
    prediction_interval_high = if (is.finite(pred_half)) pooled + pred_half else NA_real_,
    positive_cohorts = sum(yi > 0), negative_cohorts = sum(yi < 0), reason = reason
  )
}

primary <- effects[primary_compatible == TRUE]
external_all_modalities <- effects[target_selection_involved == FALSE]
meta_rows <- list(
  fit_meta(primary, "PRIMARY_BULK_TUMOR_REML_HKSJ", "PRIMARY", "REML_HKSJ",
    reason = "Frozen cross-sectional patient-level whole-tumor bulk estimand"),
  fit_meta(primary, "PRIMARY_BULK_TUMOR_DL", "METHOD_SENSITIVITY", "DL_NORMAL",
    reason = "DerSimonian-Laird sensitivity on the identical primary cohort set"),
  fit_meta(effects, "LEGACY_ALL_SIX_CONTEXTS_REML_HKSJ", "MODALITY_SENSITIVITY", "REML_HKSJ",
    reason = "Legacy context-inclusive sensitivity; not the manuscript primary estimand"),
  fit_meta(effects, "LEGACY_ALL_SIX_CONTEXTS_DL", "METHOD_SENSITIVITY", "DL_NORMAL",
    reason = "DL comparison for the legacy mixed-context set"),
  fit_meta(external_all_modalities, "EXTERNAL_ONLY_ALL_MODALITIES", "SELECTION_SENSITIVITY", "REML_HKSJ",
    reason = "Excludes GSE206325, which materially participated in target selection"),
  fit_meta(primary, "EXTERNAL_ONLY_PRIMARY_BULK", "SELECTION_SENSITIVITY", "REML_HKSJ",
    reason = "All primary-compatible cohorts are external to target selection; composition equals primary"),
  not_estimable_row("PRETREATMENT_ONLY", "TIMING_SENSITIVITY", "REML_HKSJ", "NONE",
    "No confirmed pretreatment whole-tumor PTPN22 response cohort is available; unknown timing is not relabeled as baseline"),
  fit_meta(primary, "TUMOR_WHOLE_TUMOR_COMPATIBLE_ONLY", "MODALITY_SENSITIVITY", "REML_HKSJ",
    reason = "Same composition as the frozen primary estimand"),
  fit_meta(effects[dataset != "GSE238264"], "EXCLUDE_SPATIAL_PSEUDOBULK", "MODALITY_SENSITIVITY", "REML_HKSJ",
    reason = "Retains state-specific GSE206325 and four bulk cohorts"),
  fit_meta(effects[dataset != "GSE206325"], "EXCLUDE_STATE_SPECIFIC_SINGLE_CELL", "MODALITY_SENSITIVITY", "REML_HKSJ",
    reason = "Retains spatial pseudobulk and four bulk cohorts"),
  fit_meta(effects[timing == "post_treatment"], "KNOWN_POST_TREATMENT_ONLY", "TIMING_SENSITIVITY", "REML_HKSJ",
    reason = "Known post-treatment effects across state-specific, spatial and bulk modalities"),
  fit_meta(effects[timing == "unknown"], "TIMING_UNKNOWN_BULK_ONLY", "TIMING_SENSITIVITY", "REML_HKSJ",
    reason = "Unknown-timepoint bulk cohorts; not interpreted as pretreatment"),
  fit_meta(effects[endpoint_type %chin% c("RECIST", "GEO_clinical", "RECIST_CRPR_vs_SDPD")],
    "RADIOGRAPHIC_OR_GEO_CLINICAL_ENDPOINT", "ENDPOINT_SENSITIVITY", "REML_HKSJ",
    reason = "Clinical/radiographic response labels; endpoint definitions remain heterogeneous"),
  fit_meta(effects[endpoint_type %chin% c("pathologic", "pathologic_or_clinical")],
    "PATHOLOGIC_ENDPOINT_CONTEXT", "ENDPOINT_SENSITIVITY", "REML_HKSJ",
    reason = "Only two cohorts and different modalities; descriptive small-k sensitivity"),
  fit_meta(effects[therapy_class == "monotherapy"], "ANTI_PD1_MONOTHERAPY_CONTEXT", "REGIMEN_SENSITIVITY", "REML_HKSJ",
    reason = "Three monotherapy cohorts with heterogeneous timing/modalities"),
  fit_meta(effects[therapy_class == "combination"], "COMBINATION_THERAPY_CONTEXT", "REGIMEN_SENSITIVITY", "REML_HKSJ",
    reason = "Only two combination cohorts; prediction interval not estimable")
)

for (omitted in primary$dataset) {
  meta_rows[[length(meta_rows) + 1L]] <- fit_meta(
    primary[dataset != omitted], paste0("PRIMARY_LOO_OMIT_", omitted), "LEAVE_ONE_STUDY_OUT",
    "REML_HKSJ", omitted = omitted,
    reason = "Leave-one-study-out refit of the frozen primary bulk estimand"
  )
}

meta_models <- rbindlist(meta_rows, fill = TRUE)
fwrite(meta_models, file.path(table_dir, "PTPN22_FINAL_META_MODELS.tsv"), sep = "\t", na = "NA")

# -----------------------------------------------------------------------------
# GSE235863 clonotype identity sensitivity.
# -----------------------------------------------------------------------------

x <- fread(mapping_file)
required_mapping_columns <- c(
  "patient", "response_group", "sample", "timepoint", "tissue", "barcode",
  "cdr3_pair_aa", "clone.id", "ptpn22_log1p_cp10k"
)
stopifnot(all(required_mapping_columns %chin% names(x)))
stopifnot(nrow(x) == 58872L, uniqueN(x$barcode) == nrow(x))
stopifnot(!anyNA(x$cdr3_pair_aa), all(nzchar(x$cdr3_pair_aa)))

primary_patients <- c("P1", "P5", "P11", "P15", "P18", "P26", "P27")
primary_cells <- x[patient %chin% primary_patients & tissue == "P" & timepoint %chin% c("pre", "post")]
primary_cells[, author_key := paste(patient, clone.id, sep = "::")]
primary_cells[, sequence_key := paste(patient, cdr3_pair_aa, sep = "::")]

identity_pairs <- unique(primary_cells[, .(patient, author_key, sequence_key)])
identity_pairs[, n_sequence_per_author := uniqueN(sequence_key), by = author_key]
identity_pairs[, n_author_per_sequence := uniqueN(author_key), by = sequence_key]
identity_pairs[, one_to_one := n_sequence_per_author == 1L & n_author_per_sequence == 1L]
primary_cells[identity_pairs, on = .(author_key, sequence_key), one_to_one := i.one_to_one]
one_to_one_cell_fraction <- mean(primary_cells$one_to_one)
one_to_one_pair_fraction <- mean(identity_pairs$one_to_one)

compute_identity_effects <- function(key_column, identity_label) {
  baseline <- primary_cells[timepoint == "pre", .(
    baseline_mean_ptpn22 = mean(ptpn22_log1p_cp10k),
    baseline_clone_cells = .N
  ), by = c("patient", "response_group", key_column)]
  setnames(baseline, key_column, "identity_key")
  future <- unique(primary_cells[timepoint == "post", c("patient", key_column), with = FALSE])
  setnames(future, key_column, "identity_key")
  future[, future_persistent_observed := TRUE]
  baseline[future, on = .(patient, identity_key), future_persistent_observed := i.future_persistent_observed]
  baseline[is.na(future_persistent_observed), future_persistent_observed := FALSE]
  out <- baseline[, .(
    n_baseline_clonotypes = .N,
    n_future_persistent = sum(future_persistent_observed),
    n_future_disappearing = sum(!future_persistent_observed),
    mean_baseline_ptpn22_persistent = mean(baseline_mean_ptpn22[future_persistent_observed]),
    mean_baseline_ptpn22_disappearing = mean(baseline_mean_ptpn22[!future_persistent_observed]),
    effect_persistent_minus_disappearing =
      mean(baseline_mean_ptpn22[future_persistent_observed]) - mean(baseline_mean_ptpn22[!future_persistent_observed])
  ), by = .(patient, response_group)]
  out[, identity_definition := identity_label]
  out[, effect_direction := fifelse(effect_persistent_minus_disappearing > 0, "POSITIVE",
    fifelse(effect_persistent_minus_disappearing < 0, "NEGATIVE", "ZERO"))]
  out[]
}

author_effects <- compute_identity_effects("author_key", "patient_plus_author_clone.id")
sequence_effects <- compute_identity_effects("sequence_key", "patient_plus_paired_CDR3alpha_CDR3beta")
v6 <- fread(v6_effect_file)[analysis_role == "PRIMARY_BLOOD"]
author_effects[v6, on = "patient", frozen_v6_effect := i.mean_difference_persistent_minus_disappearing]
author_effects[, recomputed_minus_frozen := effect_persistent_minus_disappearing - frozen_v6_effect]
stopifnot(max(abs(author_effects$recomputed_minus_frozen)) < 1e-12)

comparison <- merge(
  author_effects[, .(patient, author_effect = effect_persistent_minus_disappearing)],
  sequence_effects[, .(patient, sequence_effect = effect_persistent_minus_disappearing)],
  by = "patient"
)
comparison[, sequence_minus_author := sequence_effect - author_effect]
comparison[, same_direction := sign(sequence_effect) == sign(author_effect)]

same_direction_n <- sum(comparison$same_direction)
sequence_positive_n <- sum(comparison$sequence_effect > 0)
sequence_median <- median(comparison$sequence_effect)
median_abs_shift <- median(abs(comparison$sequence_minus_author))
clone_gate <- if (
  same_direction_n == 7L && sequence_positive_n == 7L && sequence_median > 0 && one_to_one_cell_fraction >= 0.95
) {
  "CLONE_ID_ROBUST"
} else if (sequence_positive_n >= 6L && sequence_median > 0) {
  "CLONE_ID_PARTIAL"
} else {
  "CLONE_ID_SENSITIVE"
}

clone_table <- rbindlist(list(author_effects, sequence_effects), fill = TRUE)
clone_table[comparison, on = "patient", `:=`(
  author_effect = i.author_effect,
  sequence_effect = i.sequence_effect,
  sequence_minus_author = i.sequence_minus_author,
  same_direction_between_definitions = i.same_direction
)]
clone_table[, one_to_one_cell_fraction := one_to_one_cell_fraction]
clone_table[, one_to_one_pair_fraction := one_to_one_pair_fraction]
clone_table[, gate := clone_gate]
fwrite(clone_table, file.path(table_dir, "PTPN22_CLONE_IDENTITY_SENSITIVITY.tsv"), sep = "\t", na = "NA")

clone_summary <- data.table(
  gate = clone_gate,
  n_primary_patients = 7L,
  author_positive_patients = sum(author_effects$effect_persistent_minus_disappearing > 0),
  sequence_positive_patients = sequence_positive_n,
  same_direction_patients = same_direction_n,
  author_median_effect = median(author_effects$effect_persistent_minus_disappearing),
  sequence_median_effect = sequence_median,
  median_absolute_patient_effect_shift = median_abs_shift,
  one_to_one_cell_fraction = one_to_one_cell_fraction,
  one_to_one_pair_fraction = one_to_one_pair_fraction,
  n_unique_author_identities = uniqueN(primary_cells$author_key),
  n_unique_sequence_identities = uniqueN(primary_cells$sequence_key)
)
fwrite(clone_summary, file.path(table_dir, "PTPN22_CLONE_IDENTITY_SENSITIVITY_SUMMARY.tsv"), sep = "\t", na = "NA")

# -----------------------------------------------------------------------------
# Blood-versus-tumor label audit for every released longitudinal TCR sample.
# -----------------------------------------------------------------------------

parse_soft_samples <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE, encoding = "UTF-8")
  starts <- grep("^\\^SAMPLE = ", lines)
  ends <- c(starts[-1L] - 1L, length(lines))
  rbindlist(lapply(seq_along(starts), function(index) {
    block <- lines[starts[[index]]:ends[[index]]]
    first_value <- function(prefix) {
      hit <- block[startsWith(block, prefix)]
      if (!length(hit)) return(NA_character_)
      sub(prefix, "", hit[[1L]], fixed = TRUE)
    }
    data.table(gsm = first_value("^SAMPLE = "), title = first_value("!Sample_title = "))
  }))
}

soft <- parse_soft_samples(soft_file)
sample_audit <- unique(x[, .(sample, patient, timepoint, tissue, response_group)])
sample_audit[, sample_phase := fifelse(timepoint == "pre", "baseline", "follow-up")]
sample_audit[, tissue_source := fifelse(tissue == "P", "peripheral blood", "liver tumor")]
sample_audit[, expected_title_prefix := paste0(
  fifelse(timepoint == "pre", "Pre-treatment", "Post-treatment"), " ",
  fifelse(tissue == "P", "blood", "liver tumor"), " of patient ", patient, " ("
)]
sample_audit[, gsm := vapply(expected_title_prefix, function(prefix) {
  hit <- soft[startsWith(title, prefix), gsm]
  if (length(hit) == 1L) hit else NA_character_
}, character(1L))]
sample_audit[, exact_geo_title := vapply(expected_title_prefix, function(prefix) {
  hit <- soft[startsWith(title, prefix), title]
  if (length(hit) == 1L) hit else NA_character_
}, character(1L))]

complete_pairs <- sample_audit[, .(
  has_baseline = any(timepoint == "pre"),
  has_followup = any(timepoint == "post")
), by = .(patient, tissue)]
sample_audit[complete_pairs, on = .(patient, tissue), `:=`(
  has_complete_pair = i.has_baseline & i.has_followup
)]
sample_audit[, inclusion := fcase(
  patient %chin% primary_patients & tissue == "P" & has_complete_pair, "PRIMARY_BLOOD_INCLUDED",
  patient %chin% primary_patients & tissue == "T" & has_complete_pair, "SECONDARY_TUMOR_SENSITIVITY_INCLUDED",
  tissue == "T" & !has_complete_pair, "EXCLUDED_INCOMPLETE_TUMOR_PAIR",
  !has_complete_pair, "EXCLUDED_NO_BASELINE_FOLLOWUP_PAIR",
  default = "EXCLUDED_NOT_IN_FROZEN_ELIGIBLE_PATIENT_SET"
)]
sample_audit[, provenance := "frozen 58,872/58,872 exact barcode-linked table + released sample token + NCBI GEO family SOFT exact title"]
sample_audit[, ambiguity_flag := patient == "P27" & sample == "P27-pre-P"]
sample_audit[, ambiguity_detail := fifelse(ambiguity_flag,
  "Released sample token and GEO title identify P27; one GEO characteristics field says P18; title/token retained as authority",
  "NONE")]
setorder(sample_audit, patient, tissue, timepoint)
fwrite(sample_audit[, .(
  patient, sample_phase, sample_id = sample, gsm, exact_geo_title, tissue_code = tissue,
  tissue_source, response_group, inclusion, provenance, ambiguity_flag, ambiguity_detail
)], file.path(metadata_dir, "GSE235863_TEMPORAL_COMPARTMENT_AUDIT.tsv"), sep = "\t", na = "NA")

# -----------------------------------------------------------------------------
# Manuscript-ready cohort characteristics table (Supplementary Table 1 input).
# -----------------------------------------------------------------------------

manuscript_table <- roles[, .(
  cohort = cohort_id,
  accession,
  N,
  `R/NR` = R_NR,
  treatment,
  setting = treatment_setting_line_if_known,
  `line of therapy if known` = treatment_setting_line_if_known,
  `sample timing` = sample_timing,
  tissue,
  `assay/platform` = modality,
  endpoint = response_definition,
  `response definition` = response_definition,
  role,
  `primary meta yes/no` = enters_primary_meta,
  `sensitivity meta yes/no` = enters_sensitivity_meta,
  `major limitation` = major_limitation
)]
fwrite(manuscript_table, file.path(table_dir, "HCC_MANUSCRIPT_COHORT_CHARACTERISTICS.tsv"), sep = "\t", na = "NA")

cat("\nMETA MODELS\n")
print(meta_models)
cat("\nCLONE SUMMARY\n")
print(clone_summary)
cat("\nCOMPARTMENT AUDIT COUNTS\n")
print(sample_audit[, .N, by = .(tissue_source, inclusion)])
cat("\nOUTPUT HASHES\n")
for (path in c(
  file.path(table_dir, "PTPN22_FINAL_META_MODELS.tsv"),
  file.path(table_dir, "PTPN22_CLONE_IDENTITY_SENSITIVITY.tsv"),
  file.path(table_dir, "PTPN22_CLONE_IDENTITY_SENSITIVITY_SUMMARY.tsv"),
  file.path(metadata_dir, "GSE235863_TEMPORAL_COMPARTMENT_AUDIT.tsv"),
  file.path(table_dir, "HCC_MANUSCRIPT_COHORT_CHARACTERISTICS.tsv")
)) cat(path, digest(file = path, algo = "sha256"), "\n")
cat("\nSESSION INFO\n")
print(sessionInfo())
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
