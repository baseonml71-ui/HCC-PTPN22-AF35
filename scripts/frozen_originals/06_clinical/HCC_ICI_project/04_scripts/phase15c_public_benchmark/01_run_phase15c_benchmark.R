options(stringsAsFactors = FALSE, width = 240)
invisible(try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE))

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(logistf)
})

project_dir <- normalizePath("HCC_ICI_project", winslash = "/", mustWork = TRUE)
metadata_dir <- file.path(project_dir, "01_metadata", "phase15c_public_benchmark")
processed_dir <- file.path(project_dir, "03_processed_data", "phase15c_public_benchmark")
table_dir <- file.path(project_dir, "05_results", "tables", "phase15c_public_benchmark")
result_dir <- file.path(project_dir, "05_results", "phase15c_public_benchmark")
log_dir <- file.path(project_dir, "09_logs", "phase15c_public_benchmark")
invisible(lapply(c(metadata_dir, processed_dir, table_dir, result_dir, log_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

contract_path <- file.path(metadata_dir, "PHASE15C_BENCHMARK_CONTRACT.tsv")
contract_sha256 <- "bf643cc7959c09da6c1765ce5cb6005c33566b02d4e2e8ce6d8e44949f6c844f"
erp_rds_path <- file.path(project_dir, "03_processed_data", "pretreatment_clinical_translation", "ERP117672_gene_expression_response_blind.rds")
g302_rds_path <- file.path(project_dir, "03_processed_data", "phase15b_public_replication", "GSE302495_RESPONSE_BLIND_MATRIX_AND_SCORES.rds")
erp_map_path <- file.path(project_dir, "01_metadata", "pretreatment_clinical_translation", "ERP117672_PATIENT_RESPONSE_MAP.tsv")
g302_map_path <- file.path(project_dir, "01_metadata", "phase15b_public_replication", "GSE302495_PATIENT_RESPONSE_MAP.tsv")
axis_path <- file.path(project_dir, "01_metadata", "clinical_translation", "PTPN22_CLINICAL_AXES_FROZEN.tsv")
comparator_path <- file.path(project_dir, "01_metadata", "PHASE3_PTPN22_SPECIFICITY_COMPARATOR_PROGRAMS.tsv")

stopifnot(all(file.exists(c(contract_path, erp_rds_path, g302_rds_path, erp_map_path, g302_map_path, axis_path, comparator_path))))
stopifnot(
  digest(contract_path, algo = "sha256", file = TRUE) == contract_sha256,
  digest(erp_rds_path, algo = "sha256", file = TRUE) == "ec4c1fd4c25d04e81167766c252da72047256a0352cd416920aba74847cc63f6",
  digest(g302_rds_path, algo = "sha256", file = TRUE) == "b138295da8af3f7f8965dc68ff429aa10cf63e90f3e0ce286a1185e299322d1e",
  digest(erp_map_path, algo = "sha256", file = TRUE) == "4f28ae48fabf509f1248ac808d2c07b1a6beaf37b4f8f14fbeb999431f7082b2",
  digest(g302_map_path, algo = "sha256", file = TRUE) == "fd2d8101719bd178e0b28680645e86319a30e485f210dd8f27fae04dc1c83b8d"
)

safe_z <- function(x) {
  x <- as.numeric(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x)) / s
}

score_program <- function(mat, genes, minimum_fraction = 0.8) {
  present <- intersect(genes, rownames(mat))
  variable <- present[vapply(present, function(g) is.finite(sd(mat[g, ])) && sd(mat[g, ]) > 0, logical(1L))]
  required <- ceiling(minimum_fraction * length(genes))
  eligible <- length(present) >= required && length(variable) >= required
  score <- rep(NA_real_, ncol(mat))
  if (eligible) score <- colMeans(t(scale(t(mat[variable, , drop = FALSE]))))
  list(
    score = score, eligible = eligible, requested_n = length(genes), present_n = length(present),
    variable_n = length(variable), missing = setdiff(genes, present)
  )
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
  list(
    p = (1 + sum(abs(permuted) >= abs(observed) - 1e-12)) / (draws + 1),
    method = "FIXED_SEED_MONTE_CARLO_LABEL_PERMUTATION", allocations = draws
  )
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

effect_summary <- function(values, groups, ids, dataset, variable, variable_role, seed, draws = 20000L) {
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
    dataset = dataset, variable = variable, variable_role = variable_role,
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

make_blind_residual_object <- function(dataset, candidates, benchmarks, eligible_ids, tier_a) {
  blind <- merge(candidates, benchmarks, by = "patient_id", all = FALSE, sort = FALSE)
  blind <- blind[match(eligible_ids, patient_id)]
  stopifnot(nrow(blind) == length(eligible_ids), !anyNA(blind))
  full <- data.table(patient_id = blind$patient_id)
  loo <- list()
  for (candidate in c("PTPN22", "AF35")) {
    for (comparator in tier_a) {
      key <- paste(candidate, comparator, sep = "__residual_on__")
      full[[key]] <- safe_z(residuals(lm(blind[[candidate]] ~ blind[[comparator]])))
      loo[[key]] <- setNames(lapply(seq_len(nrow(blind)), function(i) {
        residual_i <- safe_z(residuals(lm(blind[[candidate]][-i] ~ blind[[comparator]][-i])))
        setNames(residual_i, blind$patient_id[-i])
      }), blind$patient_id)
    }
  }
  list(
    dataset = dataset, patient_ids = blind$patient_id, full_residuals = full, loo_residuals = loo,
    tier_a = tier_a, residualization_scope = "ENDPOINT_ELIGIBLE_IDS_WITHOUT_OUTCOME_LABELS",
    outcome_joined = FALSE
  )
}

incremental_row <- function(data, comparator, candidate) {
  d <- data.table(y = as.integer(data$group == "R"), b = safe_z(data[[comparator]]), candidate_z = safe_z(data[[candidate]]))
  base <- tryCatch(logistf(y ~ b, data = d, firth = TRUE, pl = FALSE), error = function(e) NULL)
  full <- tryCatch(logistf(y ~ b + candidate_z, data = d, firth = TRUE, pl = TRUE), error = function(e) NULL)
  if (is.null(base) || is.null(full)) stop(sprintf("Firth model failed for %s / %s / %s", unique(data$dataset), comparator, candidate))
  beta <- unname(coef(full)["candidate_z"])
  nested_stat <- 2 * (unname(full$loglik["full"]) - unname(base$loglik["full"]))
  nested_p <- pchisq(nested_stat, df = 1, lower.tail = FALSE)
  loo_beta <- vapply(seq_len(nrow(d)), function(i) {
    fit <- tryCatch(logistf(y ~ b + candidate_z, data = d[-i], firth = TRUE, pl = FALSE), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    unname(coef(fit)["candidate_z"])
  }, numeric(1L))
  loo_success <- sum(is.finite(loo_beta))
  events <- sum(d$y == 1L)
  epp <- events / 2
  warnings <- character()
  if (epp < 10) warnings <- c(warnings, sprintf("LOW_EVENT_PER_PARAMETER:%.1f", epp))
  if (loo_success < nrow(d)) warnings <- c(warnings, sprintf("LOO_FIT_FAILURE:%d/%d", nrow(d) - loo_success, nrow(d)))
  if (max(abs(full$conv)) >= 1e-3) warnings <- c(warnings, "FULL_MODEL_CONVERGENCE_WARNING")
  data.table(
    dataset = unique(data$dataset), comparator = comparator, candidate = candidate,
    n = nrow(d), n_R = events, n_NR = nrow(d) - events, parameters_excluding_intercept = 2L,
    events_per_parameter = epp, candidate_conditional_coefficient = beta,
    candidate_conditional_or_per_sd = exp(beta), conditional_ci_low = exp(unname(full$ci.lower["candidate_z"])),
    conditional_ci_high = exp(unname(full$ci.upper["candidate_z"])), conditional_profile_p = unname(full$prob["candidate_z"]),
    coefficient_direction = fifelse(beta > 0, "POSITIVE", fifelse(beta < 0, "NEGATIVE", "ZERO")),
    convergence = ifelse(max(abs(full$conv)) < 1e-3, "CONVERGED", "CONVERGENCE_WARNING"),
    convergence_max = max(abs(full$conv)), nested_penalized_lr_statistic = nested_stat,
    nested_penalized_lr_df = 1L, nested_penalized_lr_p = nested_p,
    loo_candidate_coefficient_positive_fraction = sum(loo_beta > 0, na.rm = TRUE) / length(loo_beta),
    loo_candidate_coefficient_same_direction_fraction = sum(sign(loo_beta) == sign(beta), na.rm = TRUE) / length(loo_beta),
    loo_fit_success_n = loo_success, loo_total_n = length(loo_beta),
    conditional_uncertainty_nonnull = (exp(unname(full$ci.lower["candidate_z"])) > 1 || exp(unname(full$ci.upper["candidate_z"])) < 1 || nested_p < 0.05),
    model_warning = ifelse(length(warnings), paste(warnings, collapse = ";"), "NONE"),
    model_scope = "ONE_COMPARATOR_PLUS_ONE_CANDIDATE", biological_replicate = "patient"
  )
}

auroc_rank <- function(y, probability) {
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  (sum(rank(probability, ties.method = "average")[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

average_precision <- function(y, probability) {
  order_index <- order(probability, decreasing = TRUE)
  y_ordered <- y[order_index]
  precision <- cumsum(y_ordered) / seq_along(y_ordered)
  mean(precision[y_ordered == 1L])
}

descriptive_metric_row <- function(data, comparator, candidate = NA_character_) {
  d <- data.table(y = as.integer(data$group == "R"), b = safe_z(data[[comparator]]))
  if (!is.na(candidate)) d[, candidate_z := safe_z(data[[candidate]])]
  formula <- if (is.na(candidate)) y ~ b else y ~ b + candidate_z
  fit <- logistf(formula, data = d, firth = TRUE, pl = FALSE)
  apparent <- as.numeric(predict(fit, newdata = d, type = "response"))
  loo_probability <- vapply(seq_len(nrow(d)), function(i) {
    loo_fit <- tryCatch(logistf(formula, data = d[-i], firth = TRUE, pl = FALSE), error = function(e) NULL)
    if (is.null(loo_fit)) return(NA_real_)
    tryCatch(as.numeric(predict(loo_fit, newdata = d[i], type = "response")), error = function(e) NA_real_)
  }, numeric(1L))
  clipped <- pmin(pmax(loo_probability, 1e-6), 1 - 1e-6)
  available <- sum(is.finite(loo_probability))
  data.table(
    dataset = unique(data$dataset), comparator = comparator,
    model = ifelse(is.na(candidate), "COMPARATOR_ALONE", paste0("COMPARATOR_PLUS_", candidate)),
    candidate = ifelse(is.na(candidate), "NONE", candidate), n = nrow(d), n_R = sum(d$y), n_NR = sum(d$y == 0L),
    apparent_AUROC = auroc_rank(d$y, apparent), apparent_PR_AUC = average_precision(d$y, apparent),
    PR_AUC_method = "AVERAGE_PRECISION", LOO_Brier = mean((loo_probability - d$y)^2, na.rm = TRUE),
    LOO_log_loss = -mean(d$y * log(clipped) + (1 - d$y) * log(1 - clipped), na.rm = TRUE),
    LOO_predictions_available = available, LOO_predictions_expected = nrow(d),
    interpretation_status = "NOT_STABLE_FOR_INTERPRETATION",
    stability_reason = sprintf("SMALL_EVENT_COUNT_AND_LOW_EVENTS_PER_PARAMETER:n_R=%d", sum(d$y)),
    threshold_used = FALSE, clinical_prediction_claim = FALSE, biological_replicate = "patient"
  )
}

# Phase 15C contract has already been frozen and hashed. Everything above this point is definition-only.
contract <- fread(contract_path)
stopifnot(nrow(contract) == 11L, all(contract[eligible_yes_no == "YES", benchmark_name] == c("broad_immune", "CD274", "IFNg_response", "cytotoxicity", "TCR_activation", "checkpoint_exhaustion")))
tier_a <- contract[tier_A_or_B == "A" & eligible_yes_no == "YES", benchmark_name]
eligible_benchmarks <- contract[eligible_yes_no == "YES", benchmark_name]
stopifnot(identical(tier_a, c("broad_immune", "CD274", "IFNg_response", "cytotoxicity")))

axes <- fread(axis_path)
af_genes <- strsplit(axes[component == "PTPN22_activation_feedback_35", exact_features], ";", fixed = TRUE)[[1L]]
stopifnot(length(af_genes) == 35L)

erp_object <- readRDS(erp_rds_path)
g302_object <- readRDS(g302_rds_path)
stopifnot(identical(erp_object$response_joined, FALSE), identical(g302_object$response_joined, FALSE))
erp_mat <- erp_object$gene_expression
g302_mat <- g302_object$gene_expression_log2_cpm_plus_1
stopifnot(is.matrix(erp_mat), is.matrix(g302_mat), ncol(erp_mat) == 40L, ncol(g302_mat) == 38L)

program_ids <- c(
  broad_immune = "broad_immune_abundance", IFNg_response = "IFN_response", cytotoxicity = "cytotoxicity",
  TCR_activation = "TCR_activation", checkpoint_exhaustion = "checkpoint_exhaustion"
)
program_definitions <- fread(comparator_path)[program %chin% unname(program_ids) & !startsWith(gene, "DERIVED:")]

build_benchmark_scores <- function(dataset, mat, patient_ids) {
  scores <- data.table(dataset = dataset, patient_id = patient_ids)
  coverage <- list()
  for (benchmark in names(program_ids)) {
    genes <- program_definitions[program == program_ids[[benchmark]], gene]
    scored <- score_program(mat, genes, 0.8)
    stopifnot(scored$eligible)
    scores[[benchmark]] <- scored$score
    coverage[[benchmark]] <- scored
  }
  stopifnot("CD274" %chin% rownames(mat), is.finite(sd(mat["CD274", ])), sd(mat["CD274", ]) > 0)
  scores[, CD274 := as.numeric(mat["CD274", ])]
  setcolorder(scores, c("dataset", "patient_id", eligible_benchmarks))
  scores[, `:=`(scoring_scope = "WITHIN_COHORT_RESPONSE_BLIND", outcome_joined = FALSE)]
  list(scores = scores, coverage = coverage)
}

erp_benchmark <- build_benchmark_scores("ERP117672", erp_mat, colnames(erp_mat))
g302_benchmark <- build_benchmark_scores("GSE302495", g302_mat, g302_object$patient_ids)

erp_af <- score_program(erp_mat, af_genes, 0.8)
g302_af <- score_program(g302_mat, af_genes, 0.8)
stopifnot(erp_af$eligible, g302_af$eligible)
erp_candidates <- data.table(patient_id = colnames(erp_mat), PTPN22 = safe_z(erp_mat["PTPN22", ]), AF35 = safe_z(erp_af$score))
g302_candidates <- data.table(patient_id = g302_object$patient_ids, PTPN22 = g302_object$molecular_scores$PTPN22_expression, AF35 = g302_object$molecular_scores$AF35_score)

# Existing frozen scores are checked without outcomes; only positive affine rescaling is permitted by the frozen rule.
stopifnot(
  max(abs(safe_z(erp_benchmark$scores$broad_immune) - safe_z(fread(file.path(project_dir, "05_results", "tables", "pretreatment_clinical_translation", "ERP117672_PRETREATMENT_PATIENT_SCORES.tsv"))$broad_immune_z))) < 1e-12,
  max(abs(safe_z(g302_benchmark$scores$broad_immune) - safe_z(g302_object$molecular_scores$broad_immune_score))) < 1e-12,
  max(abs(safe_z(g302_af$score) - safe_z(g302_object$molecular_scores$AF35_score))) < 1e-12
)

erp_score_path <- file.path(table_dir, "ERP117672_BENCHMARK_SCORES_RESPONSE_BLIND.tsv")
g302_score_path <- file.path(table_dir, "GSE302495_BENCHMARK_SCORES_RESPONSE_BLIND.tsv")
fwrite(erp_benchmark$scores, erp_score_path, sep = "\t", quote = FALSE, na = "NA")
fwrite(g302_benchmark$scores, g302_score_path, sep = "\t", quote = FALSE, na = "NA")

# Endpoint eligibility is loaded without response labels so residuals and every LOO residual fold remain response-blind.
erp_blind_eligibility <- fread(erp_map_path, select = c("patient_id", "primary_included"))[primary_included == TRUE, patient_id]
g302_blind_eligibility <- g302_object$patient_ids
erp_residual_object <- make_blind_residual_object("ERP117672", erp_candidates, erp_benchmark$scores[, c("patient_id", eligible_benchmarks), with = FALSE], erp_blind_eligibility, tier_a)
g302_residual_object <- make_blind_residual_object("GSE302495", g302_candidates, g302_benchmark$scores[, c("patient_id", eligible_benchmarks), with = FALSE], g302_blind_eligibility, tier_a)
erp_residual_rds <- file.path(processed_dir, "ERP117672_TIERA_RESIDUALS_RESPONSE_BLIND.rds")
g302_residual_rds <- file.path(processed_dir, "GSE302495_TIERA_RESIDUALS_RESPONSE_BLIND.rds")
saveRDS(erp_residual_object, erp_residual_rds, version = 3)
saveRDS(g302_residual_object, g302_residual_rds, version = 3)

score_lock_lines <- c(
  "# Phase 15C response-blind score and residual lock", "",
  sprintf("Contract SHA256: `%s`", contract_sha256),
  sprintf("ERP117672 score SHA256: `%s`", digest(erp_score_path, algo = "sha256", file = TRUE)),
  sprintf("GSE302495 score SHA256: `%s`", digest(g302_score_path, algo = "sha256", file = TRUE)),
  sprintf("ERP117672 residual-object SHA256: `%s`", digest(erp_residual_rds, algo = "sha256", file = TRUE)),
  sprintf("GSE302495 residual-object SHA256: `%s`", digest(g302_residual_rds, algo = "sha256", file = TRUE)), "",
  "Both score tables and both residual objects have `outcome_joined=FALSE`.",
  "ERP117672 endpoint eligibility was loaded without response labels; GSE302495 uses all 38 frozen baseline patients.",
  "This lock was written before response labels were loaded into the analysis."
)
writeLines(score_lock_lines, file.path(result_dir, "PHASE15C_RESPONSE_BLIND_SCORE_LOCK.md"), useBytes = TRUE)

# Response labels are loaded only after the score/residual lock above.
erp_map <- fread(erp_map_path)[primary_included == TRUE, .(patient_id, group = fifelse(primary_group == "Responder", "R", "NR"))]
g302_map <- fread(g302_map_path)[, .(patient_id, group = response_group)]
stopifnot(nrow(erp_map) == 35L, sum(erp_map$group == "R") == 6L, sum(erp_map$group == "NR") == 29L)
stopifnot(nrow(g302_map) == 38L, sum(g302_map$group == "R") == 12L, sum(g302_map$group == "NR") == 26L)

assemble_analysis_data <- function(dataset, map, candidates, benchmarks) {
  d <- merge(map, candidates, by = "patient_id", all.x = TRUE, sort = FALSE)
  d <- merge(d, benchmarks[, c("patient_id", eligible_benchmarks), with = FALSE], by = "patient_id", all.x = TRUE, sort = FALSE)
  d[, dataset := dataset]
  setcolorder(d, c("dataset", "patient_id", "group", "PTPN22", "AF35", eligible_benchmarks))
  stopifnot(!anyNA(d), !anyDuplicated(d$patient_id))
  d
}

erp_data <- assemble_analysis_data("ERP117672", erp_map, erp_candidates, erp_benchmark$scores)
g302_data <- assemble_analysis_data("GSE302495", g302_map, g302_candidates, g302_benchmark$scores)

run_univariate <- function(data, seed_base) {
  variables <- c("PTPN22", "AF35", eligible_benchmarks)
  roles <- c(PTPN22 = "CANDIDATE", AF35 = "CANDIDATE", setNames(rep("BENCHMARK", length(eligible_benchmarks)), eligible_benchmarks))
  rows <- lapply(seq_along(variables), function(i) {
    variable <- variables[[i]]
    effect_summary(data[[variable]], data$group, data$patient_id, unique(data$dataset), variable, roles[[variable]], seed_base + i)
  })
  result <- rbindlist(rows)
  result[, permutation_p_BH := p.adjust(permutation_p_two_sided, method = "BH")]
  result[, BH_family := "PTPN22_AF35_AND_ALL_ELIGIBLE_BENCHMARKS"]
  result
}

erp_univariate <- run_univariate(erp_data, 1500000000L)
g302_univariate <- run_univariate(g302_data, 1500000100L)
fwrite(erp_univariate, file.path(table_dir, "ERP117672_UNIVARIATE_BENCHMARK_EFFECTS.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(g302_univariate, file.path(table_dir, "GSE302495_UNIVARIATE_BENCHMARK_EFFECTS.tsv"), sep = "\t", quote = FALSE, na = "NA")

run_incremental <- function(data) {
  rbindlist(lapply(tier_a, function(comparator) {
    rbindlist(lapply(c("PTPN22", "AF35"), function(candidate) incremental_row(data, comparator, candidate)))
  }))
}

erp_incremental <- run_incremental(erp_data)
g302_incremental <- run_incremental(g302_data)
fwrite(erp_incremental, file.path(table_dir, "ERP117672_INCREMENTAL_FIRTH_MODELS.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(g302_incremental, file.path(table_dir, "GSE302495_INCREMENTAL_FIRTH_MODELS.tsv"), sep = "\t", quote = FALSE, na = "NA")

run_residual_sensitivity <- function(data, residual_object, seed_base) {
  rows <- list()
  index <- 0L
  for (candidate in c("PTPN22", "AF35")) {
    for (comparator in tier_a) {
      index <- index + 1L
      key <- paste(candidate, comparator, sep = "__residual_on__")
      full <- residual_object$full_residuals
      values <- full[[key]][match(data$patient_id, full$patient_id)]
      row <- effect_summary(values, data$group, data$patient_id, unique(data$dataset), key, "RESPONSE_BLIND_RESIDUAL", seed_base + index)
      loo_g <- vapply(data$patient_id, function(omitted_id) {
        residual_i <- residual_object$loo_residuals[[key]][[omitted_id]]
        group_i <- data$group[match(names(residual_i), data$patient_id)]
        hedges_g(as.numeric(residual_i), group_i)
      }, numeric(1L))
      row[, `:=`(
        candidate = candidate, comparator = comparator,
        loo_hedges_g_min = min(loo_g, na.rm = TRUE), loo_hedges_g_max = max(loo_g, na.rm = TRUE),
        loo_positive_fraction = mean(loo_g > 0, na.rm = TRUE),
        loo_same_direction_fraction = mean(sign(loo_g) == sign(hedges_g), na.rm = TRUE),
        residualization_timing = "BEFORE_OUTCOME_JOIN", residualization_role = "SPECIFICITY_FALSIFICATION_ONLY"
      )]
      rows[[key]] <- row
    }
  }
  result <- rbindlist(rows)
  setcolorder(result, c("dataset", "comparator", "candidate", setdiff(names(result), c("dataset", "comparator", "candidate"))))
  result
}

erp_residual <- run_residual_sensitivity(erp_data, erp_residual_object, 1500000200L)
g302_residual <- run_residual_sensitivity(g302_data, g302_residual_object, 1500000300L)

# Deterministic broad-immune effects must reproduce the frozen v14/v15B analyses.
erp_old <- fread(file.path(project_dir, "05_results", "tables", "pretreatment_clinical_translation", "ERP117672_BROAD_IMMUNE_RESPONSE_EFFECTS.tsv"))
g302_old <- fread(file.path(project_dir, "05_results", "tables", "phase15b_public_replication", "GSE302495_ALL_RECIST_CONTINUOUS_EFFECTS.tsv"))
erp_expected <- rbindlist(list(
  erp_old[population == "PRIMARY_RESPONSE_EVALUABLE" & metric == "PTPN22_residual_after_broad_immune",
    .(candidate = "PTPN22", old_hedges_g = hedges_g, old_firth_or = firth_or_per_sd, old_loo_positive = loo_positive_fraction)],
  erp_old[population == "PRIMARY_RESPONSE_EVALUABLE" & metric == "AF_program_residual_after_broad_immune",
    .(candidate = "AF35", old_hedges_g = hedges_g, old_firth_or = firth_or_per_sd, old_loo_positive = loo_positive_fraction)]
))
g302_expected <- rbindlist(list(
  g302_old[endpoint == "RECIST" & metric == "PTPN22_broad_immune_residual",
    .(candidate = "PTPN22", old_hedges_g = hedges_g, old_firth_or = firth_or_per_sd, old_loo_positive = loo_positive_fraction)],
  g302_old[endpoint == "RECIST" & metric == "AF35_broad_immune_residual",
    .(candidate = "AF35", old_hedges_g = hedges_g, old_firth_or = firth_or_per_sd, old_loo_positive = loo_positive_fraction)]
))

check_reproduction <- function(current, expected, dataset) {
  checks <- lapply(c("PTPN22", "AF35"), function(candidate_name) {
    new <- current[comparator == "broad_immune" & candidate == candidate_name]
    old <- expected[candidate == candidate_name]
    stopifnot(nrow(new) == 1L, nrow(old) == 1L)
    data.table(
      dataset = dataset, candidate = candidate_name,
      hedges_g_abs_diff = abs(new$hedges_g - old$old_hedges_g),
      firth_or_abs_diff = abs(new$firth_or_per_sd - old$old_firth_or),
      loo_positive_fraction_abs_diff = abs(new$loo_positive_fraction - old$old_loo_positive)
    )
  })
  result <- rbindlist(checks)
  result[, reproduction_pass := hedges_g_abs_diff < 1e-10 & firth_or_abs_diff < 1e-8 & loo_positive_fraction_abs_diff < 1e-12]
  result
}

reproduction <- rbindlist(list(
  check_reproduction(erp_residual, erp_expected, "ERP117672"),
  check_reproduction(g302_residual, g302_expected, "GSE302495")
))
if (!all(reproduction$reproduction_pass)) stop("BROAD_IMMUNE_RESIDUAL_REPRODUCTION_FAILED")
erp_residual <- merge(erp_residual, reproduction[dataset == "ERP117672", .(candidate, broad_residual_reproduction_pass = reproduction_pass)], by = "candidate", all.x = TRUE, sort = FALSE)
g302_residual <- merge(g302_residual, reproduction[dataset == "GSE302495", .(candidate, broad_residual_reproduction_pass = reproduction_pass)], by = "candidate", all.x = TRUE, sort = FALSE)
erp_residual[comparator != "broad_immune", broad_residual_reproduction_pass := NA]
g302_residual[comparator != "broad_immune", broad_residual_reproduction_pass := NA]
fwrite(erp_residual, file.path(table_dir, "ERP117672_TIERA_RESIDUAL_SENSITIVITY.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(g302_residual, file.path(table_dir, "GSE302495_TIERA_RESIDUAL_SENSITIVITY.tsv"), sep = "\t", quote = FALSE, na = "NA")

run_descriptive_metrics <- function(data) {
  rbindlist(lapply(tier_a, function(comparator) {
    rbindlist(list(
      descriptive_metric_row(data, comparator),
      descriptive_metric_row(data, comparator, "PTPN22"),
      descriptive_metric_row(data, comparator, "AF35")
    ))
  }))
}

descriptive_metrics <- rbindlist(list(run_descriptive_metrics(erp_data), run_descriptive_metrics(g302_data)))
fwrite(descriptive_metrics, file.path(table_dir, "PHASE15C_DESCRIPTIVE_MODEL_METRICS.tsv"), sep = "\t", quote = FALSE, na = "NA")

saveRDS(list(
  contract_sha256 = contract_sha256, eligible_benchmarks = eligible_benchmarks, tier_a = tier_a,
  ERP117672 = list(univariate = erp_univariate, incremental = erp_incremental, residual = erp_residual),
  GSE302495 = list(univariate = g302_univariate, incremental = g302_incremental, residual = g302_residual),
  reproduction = reproduction, descriptive_metrics = descriptive_metrics
), file.path(processed_dir, "PHASE15C_ANALYSIS_RESULTS.rds"), version = 3)

analysis_log <- c(
  "phase=PHASE15C_PUBLIC_CLINICAL_BENCHMARK_INCREMENTAL",
  "status=ANALYSIS_TABLES_COMPLETE",
  sprintf("contract_sha256=%s", contract_sha256),
  sprintf("eligible_tier_A=%s", paste(tier_a, collapse = ";")),
  sprintf("eligible_tier_B=%s", paste(setdiff(eligible_benchmarks, tier_a), collapse = ";")),
  "response_blind_scores_locked_before_outcome=TRUE",
  "response_blind_residuals_locked_before_outcome=TRUE",
  "broad_residual_reproduction=PASS",
  "large_model_fit=FALSE",
  "raw_fastq_used=FALSE",
  "survival_analysis=FALSE",
  "meta_analysis=FALSE",
  "new_cohort=FALSE"
)
writeLines(analysis_log, file.path(log_dir, "PHASE15C_ANALYSIS_RUN_LOG.txt"), useBytes = TRUE)
