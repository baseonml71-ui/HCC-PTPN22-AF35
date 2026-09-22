options(stringsAsFactors = FALSE, width = 240)

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(jsonlite)
  library(logistf)
})

project_dir <- normalizePath("HCC_ICI_project", winslash = "/", mustWork = TRUE)
metadata_dir <- file.path(project_dir, "01_metadata", "phase19b_published_benchmark")
processed_dir <- file.path(project_dir, "03_processed_data", "phase19b_published_benchmark")
intermediate_dir <- file.path(project_dir, "10_intermediate_files", "phase19")
log_dir <- file.path(project_dir, "09_logs", "phase19b_published_benchmark")
invisible(lapply(c(processed_dir, intermediate_dir, log_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

lock_path <- file.path(metadata_dir, "PHASE19B_RESPONSE_BLIND_TABLE_LOCK.json")
blind_object_path <- file.path(processed_dir, "PHASE19B_RESPONSE_BLIND_OBJECT.rds")
erp_map_path <- file.path(project_dir, "01_metadata", "pretreatment_clinical_translation", "ERP117672_PATIENT_RESPONSE_MAP.tsv")
g302_map_path <- file.path(project_dir, "01_metadata", "phase15b_public_replication", "GSE302495_PATIENT_RESPONSE_MAP.tsv")

stopifnot(
  digest(lock_path, algo = "sha256", file = TRUE) == "12abef9bff69b79b860943682adc7dfbe08807993ea85b254b13043114746f20",
  digest(erp_map_path, algo = "sha256", file = TRUE) == "4f28ae48fabf509f1248ac808d2c07b1a6beaf37b4f8f14fbeb999431f7082b2",
  digest(g302_map_path, algo = "sha256", file = TRUE) == "fd2d8101719bd178e0b28680645e86319a30e485f210dd8f27fae04dc1c83b8d"
)

lock <- fromJSON(lock_path, simplifyVector = FALSE)
stopifnot(
  identical(lock$status, "PHASE19B_RESPONSE_BLIND_TABLES_LOCKED_BEFORE_OUTCOME_JOIN"),
  identical(lock$outcome_joined, FALSE),
  identical(lock$eligible_signature_n, 1L),
  identical(lock$eligible_signature, "S007_ISG_20"),
  identical(lock$eligible_peer_reviewed_tier_a_paper_n, 1L),
  identical(lock$gate_prerequisite, "FAIL_LT_2_INDEPENDENT_PEER_REVIEWED_TIER_A")
)
for (item in names(lock$file_sha256)) {
  candidates <- list(
    candidate_pool = file.path(metadata_dir, "PHASE19B_PUBLISHED_SIGNATURE_CANDIDATE_POOL.tsv"),
    contract = file.path(metadata_dir, "PHASE19B_EXACT_PUBLISHED_SIGNATURE_CONTRACT.tsv"),
    tier = file.path(metadata_dir, "PHASE19B_PUBLISHED_SIGNATURE_TIER_ASSIGNMENT.tsv"),
    erp_scores = file.path(project_dir, "05_results", "tables", "phase19b_published_benchmark", "ERP117672_PUBLISHED_SIGNATURE_SCORES_RESPONSE_BLIND.tsv"),
    g302_scores = file.path(project_dir, "05_results", "tables", "phase19b_published_benchmark", "GSE302495_PUBLISHED_SIGNATURE_SCORES_RESPONSE_BLIND.tsv"),
    relationships = file.path(project_dir, "05_results", "tables", "phase19b_published_benchmark", "AF35_PUBLISHED_SIGNATURE_MOLECULAR_RELATIONSHIPS.tsv"),
    blind_object = blind_object_path,
    blind_freeze = file.path(metadata_dir, "PHASE19B_RESPONSE_BLIND_FREEZE.json")
  )
  stopifnot(digest(candidates[[item]], algo = "sha256", file = TRUE) == lock$file_sha256[[item]])
}

blind <- readRDS(blind_object_path)
stopifnot(identical(blind$outcome_joined, FALSE))

safe_z <- function(x) {
  x <- as.numeric(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x)) / s
}

hedges_g <- function(values, groups) {
  r <- values[groups == "R"]
  nr <- values[groups == "NR"]
  pooled_sd <- sqrt(((length(r) - 1L) * var(r) + (length(nr) - 1L) * var(nr)) / (length(values) - 2L))
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)
  correction <- 1 - 3 / (4 * (length(values) - 2L) - 1)
  correction * (mean(r) - mean(nr)) / pooled_sd
}

combo_cache <- new.env(parent = emptyenv())
permutation_p <- function(values, groups, draws, seed, exact_limit = 2000000L) {
  n <- length(values)
  n_r <- sum(groups == "R")
  observed <- mean(values[groups == "R"]) - mean(values[groups == "NR"])
  total <- choose(n, n_r)
  if (is.finite(total) && total <= exact_limit) {
    key <- paste(n, n_r, sep = "_")
    if (!exists(key, envir = combo_cache, inherits = FALSE)) assign(key, combn(n, n_r), envir = combo_cache)
    combos <- get(key, envir = combo_cache, inherits = FALSE)
    responder_sums <- colSums(matrix(values[as.vector(combos)], nrow = n_r))
    diffs <- responder_sums / n_r - (sum(values) - responder_sums) / (n - n_r)
    return(list(p = mean(abs(diffs) >= abs(observed) - 1e-12), method = "EXACT_ALL_LABEL_ALLOCATIONS", allocations = as.integer(total)))
  }
  set.seed(seed)
  permuted <- replicate(draws, {
    idx <- sample.int(n, n_r, replace = FALSE)
    mean(values[idx]) - mean(values[-idx])
  })
  list(p = (1 + sum(abs(permuted) >= abs(observed) - 1e-12)) / (draws + 1), method = "FIXED_SEED_MONTE_CARLO_LABEL_PERMUTATION", allocations = draws)
}

fit_univariate_firth <- function(values, groups) {
  predictor <- safe_z(values)
  fit <- logistf(I(groups == "R") ~ predictor, firth = TRUE, pl = TRUE)
  list(
    or = exp(unname(coef(fit)["predictor"])), low = exp(unname(fit$ci.lower["predictor"])),
    high = exp(unname(fit$ci.upper["predictor"])), p = unname(fit$prob["predictor"]),
    convergence_max = max(abs(fit$conv)), converged = max(abs(fit$conv)) < 1e-3
  )
}

effect_summary <- function(values, groups, ids, dataset, variable, role, seed, draws = 20000L) {
  keep <- is.finite(values) & groups %chin% c("R", "NR")
  values <- as.numeric(values[keep])
  groups <- groups[keep]
  ids <- ids[keep]
  r <- values[groups == "R"]
  nr <- values[groups == "NR"]
  observed_g <- hedges_g(values, groups)
  set.seed(seed)
  bootstrap_g <- replicate(draws, {
    rb <- sample(r, replace = TRUE)
    nrb <- sample(nr, replace = TRUE)
    hedges_g(c(rb, nrb), c(rep("R", length(rb)), rep("NR", length(nrb))))
  })
  bootstrap_g <- bootstrap_g[is.finite(bootstrap_g)]
  stopifnot(length(bootstrap_g) > 1000L)
  firth <- fit_univariate_firth(values, groups)
  permutation <- permutation_p(values, groups, draws, seed + 1L)
  loo <- vapply(seq_along(values), function(i) hedges_g(values[-i], groups[-i]), numeric(1L))
  data.table(
    dataset = dataset, variable = variable, variable_role = role,
    n = length(values), n_R = length(r), n_NR = length(nr), R_mean = mean(r), NR_mean = mean(nr),
    R_minus_NR = mean(r) - mean(nr), hedges_g = observed_g,
    hedges_g_ci_low = unname(quantile(bootstrap_g, 0.025)), hedges_g_ci_high = unname(quantile(bootstrap_g, 0.975)),
    bootstrap_draws = draws, firth_or_per_sd = firth$or, firth_ci_low = firth$low, firth_ci_high = firth$high,
    firth_p = firth$p, firth_converged = firth$converged, firth_convergence_max = firth$convergence_max,
    permutation_p_two_sided = permutation$p, permutation_method = permutation$method,
    permutation_allocations_or_draws = permutation$allocations,
    loo_hedges_g_min = min(loo, na.rm = TRUE), loo_hedges_g_max = max(loo, na.rm = TRUE),
    loo_positive_fraction = mean(loo > 0, na.rm = TRUE),
    loo_same_direction_fraction = mean(sign(loo) == sign(observed_g), na.rm = TRUE),
    biological_replicate = "patient", seed = seed
  )
}

assemble <- function(dataset, identity, residuals, map) {
  scores <- merge(
    identity[, .(patient_id, AF35, S007_ISG_20)],
    residuals[, .(patient_id, AF35_residual_on_S007_ISG_20)], by = "patient_id", sort = FALSE
  )
  d <- merge(map, scores, by = "patient_id", all.x = TRUE, sort = FALSE)
  d[, dataset := dataset]
  setcolorder(d, c("dataset", "patient_id", "group", "AF35", "S007_ISG_20", "AF35_residual_on_S007_ISG_20"))
  stopifnot(!anyNA(d), !anyDuplicated(d$patient_id))
  d
}

# Response labels are first loaded here, after the response-blind table lock was verified.
erp_map <- fread(erp_map_path)[primary_included == TRUE, .(patient_id, group = fifelse(primary_group == "Responder", "R", "NR"))]
g302_map <- fread(g302_map_path)[, .(patient_id, group = response_group)]
stopifnot(nrow(erp_map) == 35L, sum(erp_map$group == "R") == 6L, sum(erp_map$group == "NR") == 29L)
stopifnot(nrow(g302_map) == 38L, sum(g302_map$group == "R") == 12L, sum(g302_map$group == "NR") == 26L)

erp_data <- assemble("ERP117672", blind$ERP117672$identity, blind$ERP117672$residuals, erp_map)
g302_data <- assemble("GSE302495", blind$GSE302495$identity, blind$GSE302495$residuals, g302_map)

run_univariate <- function(data, seed_base) {
  result <- rbindlist(list(
    effect_summary(data$AF35, data$group, data$patient_id, unique(data$dataset), "AF35", "FROZEN_CLINICAL_STATE_UNIT", seed_base + 1L),
    effect_summary(data$S007_ISG_20, data$group, data$patient_id, unique(data$dataset), "S007_ISG_20", "ELIGIBLE_PUBLISHED_TIER_A", seed_base + 2L)
  ))
  result[, `:=`(
    permutation_p_BH_published_family = fifelse(variable == "S007_ISG_20", p.adjust(permutation_p_two_sided[variable == "S007_ISG_20"], method = "BH"), NA_real_),
    BH_family = "ELIGIBLE_PUBLISHED_SIGNATURES_ONLY_N1",
    BH_included = variable == "S007_ISG_20"
  )]
  result
}

conditional_row <- function(data) {
  d <- data.table(y = as.integer(data$group == "R"), b = safe_z(data$S007_ISG_20), af35 = safe_z(data$AF35))
  base <- logistf(y ~ b, data = d, firth = TRUE, pl = FALSE)
  full <- logistf(y ~ b + af35, data = d, firth = TRUE, pl = TRUE)
  beta <- unname(coef(full)["af35"])
  loo_beta <- vapply(seq_len(nrow(d)), function(i) {
    fit <- tryCatch(logistf(y ~ b + af35, data = d[-i], firth = TRUE, pl = FALSE), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    unname(coef(fit)["af35"])
  }, numeric(1L))
  events <- sum(d$y == 1L)
  epp <- events / 2
  warning <- character()
  if (epp < 10) warning <- c(warning, sprintf("LOW_EVENT_PER_PARAMETER:%.1f", epp))
  if (sum(is.finite(loo_beta)) < nrow(d)) warning <- c(warning, sprintf("LOO_FIT_FAILURE:%d/%d", sum(!is.finite(loo_beta)), nrow(d)))
  if (max(abs(full$conv)) >= 1e-3) warning <- c(warning, "FULL_MODEL_CONVERGENCE_WARNING")
  data.table(
    dataset = unique(data$dataset), comparator_paper_id = "S007", comparator = "S007_ISG_20", candidate = "AF35",
    base_formula = "response ~ S007_ISG_20", candidate_formula = "response ~ S007_ISG_20 + AF35",
    n = nrow(d), n_R = events, n_NR = nrow(d) - events, parameters_excluding_intercept = 2L,
    events_per_parameter = epp, af35_conditional_beta = beta, af35_conditional_or_per_sd = exp(beta),
    conditional_ci_low = exp(unname(full$ci.lower["af35"])), conditional_ci_high = exp(unname(full$ci.upper["af35"])),
    conditional_profile_p = unname(full$prob["af35"]), coefficient_direction = fifelse(beta > 0, "POSITIVE", fifelse(beta < 0, "NEGATIVE", "ZERO")),
    convergence = fifelse(max(abs(full$conv)) < 1e-3, "CONVERGED", "CONVERGENCE_WARNING"),
    convergence_max = max(abs(full$conv)), loo_af35_coefficient_positive_fraction = mean(loo_beta > 0, na.rm = TRUE),
    loo_af35_coefficient_same_direction_fraction = mean(sign(loo_beta) == sign(beta), na.rm = TRUE),
    loo_fit_success_n = sum(is.finite(loo_beta)), loo_total_n = length(loo_beta),
    model_warning = fifelse(length(warning) > 0, paste(warning, collapse = ";"), "NONE"),
    model_scope = "ONE_PUBLISHED_COMPARATOR_PLUS_AF35", biological_replicate = "patient"
  )
}

residual_row <- function(data, seed) {
  row <- effect_summary(
    data$AF35_residual_on_S007_ISG_20, data$group, data$patient_id, unique(data$dataset),
    "AF35_residual_on_S007_ISG_20", "RESPONSE_BLIND_RESIDUAL_SPECIFICITY", seed
  )
  row[, `:=`(
    comparator_paper_id = "S007", comparator = "S007_ISG_20", candidate = "AF35",
    residualization_timing = "BEFORE_OUTCOME_JOIN_ALL_RESPONSE_BLIND_SAMPLES",
    residualization_role = "SPECIFICITY_FALSIFICATION_ONLY"
  )]
  setcolorder(row, c("dataset", "comparator_paper_id", "comparator", "candidate", setdiff(names(row), c("dataset", "comparator_paper_id", "comparator", "candidate"))))
  row
}

erp_effects <- run_univariate(erp_data, 1900000000L)
g302_effects <- run_univariate(g302_data, 1900000100L)
conditional <- rbindlist(list(conditional_row(erp_data), conditional_row(g302_data)))
residual <- rbindlist(list(residual_row(erp_data, 1900000200L), residual_row(g302_data, 1900000300L)))

gate <- "AF35_PUBLISHED_SIGNATURE_BENCHMARK_NOT_EVALUABLE"
summary <- list(
  created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), response_blind_lock_sha256 = digest(lock_path, algo = "sha256", file = TRUE),
  eligible_peer_reviewed_tier_a_signature_n = 1L, eligible_peer_reviewed_paper_n = 1L,
  eligible_signature = "S007_ISG_20", gate = gate,
  gate_reason = "Fewer than 2 independent peer-reviewed Tier A signatures have exact reconstructable definitions",
  ERP117672 = list(
    n = nrow(erp_data), n_R = sum(erp_data$group == "R"), n_NR = sum(erp_data$group == "NR"),
    af35_univariate_hedges_g = erp_effects[variable == "AF35", hedges_g],
    s007_univariate_hedges_g = erp_effects[variable == "S007_ISG_20", hedges_g],
    af35_conditional_beta = conditional[dataset == "ERP117672", af35_conditional_beta],
    af35_conditional_loo_positive_fraction = conditional[dataset == "ERP117672", loo_af35_coefficient_positive_fraction],
    af35_residual_hedges_g = residual[dataset == "ERP117672", hedges_g]
  ),
  GSE302495 = list(
    n = nrow(g302_data), n_R = sum(g302_data$group == "R"), n_NR = sum(g302_data$group == "NR"),
    af35_univariate_hedges_g = g302_effects[variable == "AF35", hedges_g],
    s007_univariate_hedges_g = g302_effects[variable == "S007_ISG_20", hedges_g],
    af35_conditional_beta = conditional[dataset == "GSE302495", af35_conditional_beta],
    af35_conditional_loo_positive_fraction = conditional[dataset == "GSE302495", loo_af35_coefficient_positive_fraction],
    af35_residual_hedges_g = residual[dataset == "GSE302495", hedges_g]
  ),
  no_cutoff = TRUE, no_roc = TRUE, no_ml = TRUE, no_meta_analysis = TRUE, cohorts_not_pooled = TRUE
)

tables <- list(
  "ERP117672_PUBLISHED_SIGNATURE_EFFECTS.tsv" = erp_effects,
  "GSE302495_PUBLISHED_SIGNATURE_EFFECTS.tsv" = g302_effects,
  "PHASE19B_AF35_CONDITIONAL_PUBLISHED_SIGNATURE_MODELS.tsv" = conditional,
  "PHASE19B_AF35_PUBLISHED_SIGNATURE_RESIDUAL_SENSITIVITY.tsv" = residual
)
write_json(tables, file.path(intermediate_dir, "phase19b_outcome_tables.json"), dataframe = "rows", auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "string")
write_json(summary, file.path(intermediate_dir, "phase19b_outcome_summary.json"), auto_unbox = TRUE, pretty = TRUE, digits = 17, na = "string")
saveRDS(list(summary = summary, ERP117672_effects = erp_effects, GSE302495_effects = g302_effects, conditional = conditional, residual = residual), file.path(processed_dir, "PHASE19B_OUTCOME_ANALYSIS_RESULTS.rds"), version = 3)

writeLines(c(
  "phase=PHASE19B_OUTCOME_ANALYSIS",
  "status=COMPLETE",
  sprintf("response_blind_lock_sha256=%s", summary$response_blind_lock_sha256),
  "eligible_signature=S007_ISG_20",
  "eligible_peer_reviewed_tier_A_papers=1",
  sprintf("gate=%s", gate),
  "cohorts_pooled=FALSE",
  "cutoff_optimization=FALSE",
  "roc=FALSE",
  "machine_learning=FALSE",
  "survival_analysis=FALSE"
), file.path(log_dir, "PHASE19B_OUTCOME_ANALYSIS_LOG.txt"), useBytes = TRUE)

cat("Phase 19B outcome analysis complete\n")
cat(sprintf("Gate: %s\n", gate))
cat(sprintf("ERP AF35 conditional beta/LOO+: %.6f / %.3f\n", summary$ERP117672$af35_conditional_beta, summary$ERP117672$af35_conditional_loo_positive_fraction))
cat(sprintf("GSE AF35 conditional beta/LOO+: %.6f / %.3f\n", summary$GSE302495$af35_conditional_beta, summary$GSE302495$af35_conditional_loo_positive_fraction))
