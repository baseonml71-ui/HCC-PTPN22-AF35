options(stringsAsFactors = FALSE, width = 240)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(logistf)
})

workspace_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
project_dir <- file.path(workspace_dir, "HCC_ICI_project")
processed_dir <- file.path(project_dir, "03_processed_data", "phase15b_public_replication")
metadata_dir <- file.path(project_dir, "01_metadata", "phase15b_public_replication")
table_dir <- file.path(project_dir, "05_results", "tables", "phase15b_public_replication")
result_dir <- file.path(project_dir, "05_results", "phase15b_public_replication")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

g302_rds <- file.path(processed_dir, "GSE302495_RESPONSE_BLIND_MATRIX_AND_SCORES.rds")
g287_rds <- file.path(processed_dir, "GSE287319_RESPONSE_BLIND_MATRIX_AND_SCORES.rds")
g302_map_path <- file.path(metadata_dir, "GSE302495_PATIENT_RESPONSE_MAP.tsv")
g287_map_path <- file.path(metadata_dir, "GSE287319_PATIENT_ENDPOINT_MAP.tsv")
contract_path <- file.path(metadata_dir, "PHASE15B_ANALYSIS_CONTRACT.md")
stopifnot(all(file.exists(c(g302_rds, g287_rds, g302_map_path, g287_map_path, contract_path))))
stopifnot(
  digest(g302_rds, algo = "sha256", file = TRUE) == "b138295da8af3f7f8965dc68ff429aa10cf63e90f3e0ce286a1185e299322d1e",
  digest(g287_rds, algo = "sha256", file = TRUE) == "07994a1ada18f4aeb8e499299a2d088f1ab9a23ab448f67e7e70998723e51034",
  digest(g302_map_path, algo = "sha256", file = TRUE) == "fd2d8101719bd178e0b28680645e86319a30e485f210dd8f27fae04dc1c83b8d",
  digest(g287_map_path, algo = "sha256", file = TRUE) == "6e27b0493342772830fd853235ac5d9d1f7864c47d95eccc8e58788d8f35093f"
)

g302_object <- readRDS(g302_rds)
g287_object <- readRDS(g287_rds)
stopifnot(identical(g302_object$response_joined, FALSE), identical(g287_object$response_joined, FALSE))

safe_z <- function(x) {
  x <- as.numeric(x)
  s <- sd(x)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x)) / s
}

hedges_g <- function(values, groups, group1, group2) {
  x1 <- values[groups == group1]
  x2 <- values[groups == group2]
  n <- length(x1) + length(x2)
  pooled_sd <- sqrt(((length(x1) - 1L) * var(x1) + (length(x2) - 1L) * var(x2)) / (n - 2L))
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)
  correction <- 1 - 3 / (4 * (n - 2L) - 1)
  correction * (mean(x1) - mean(x2)) / pooled_sd
}

combo_cache <- new.env(parent = emptyenv())
permutation_p <- function(values, groups, group1, group2, draws, seed, exact_limit = 2000000L) {
  n <- length(values)
  n1 <- sum(groups == group1)
  observed <- mean(values[groups == group1]) - mean(values[groups == group2])
  total <- choose(n, n1)
  if (is.finite(total) && total <= exact_limit) {
    key <- paste(n, n1, sep = "_")
    if (!exists(key, envir = combo_cache, inherits = FALSE)) assign(key, combn(n, n1), envir = combo_cache)
    combos <- get(key, envir = combo_cache, inherits = FALSE)
    sums <- colSums(matrix(values[as.vector(combos)], nrow = n1))
    diffs <- sums / n1 - (sum(values) - sums) / (n - n1)
    return(list(p = mean(abs(diffs) >= abs(observed) - 1e-12), method = "EXACT_ALL_LABEL_ALLOCATIONS", allocations = as.integer(total)))
  }
  set.seed(seed)
  perm <- replicate(draws, {
    idx <- sample.int(n, n1, replace = FALSE)
    mean(values[idx]) - mean(values[-idx])
  })
  list(p = (1 + sum(abs(perm) >= abs(observed) - 1e-12)) / (draws + 1), method = "FIXED_SEED_MONTE_CARLO_LABEL_PERMUTATION", allocations = draws)
}

effect_summary <- function(values, groups, ids, dataset, endpoint, metric, group1, group2, seed, firth = TRUE, draws = 20000L) {
  keep <- is.finite(values) & groups %chin% c(group1, group2)
  values <- values[keep]
  groups <- groups[keep]
  ids <- ids[keep]
  x1 <- values[groups == group1]
  x2 <- values[groups == group2]
  observed_g <- hedges_g(values, groups, group1, group2)
  set.seed(seed)
  boot_g <- replicate(draws, {
    b1 <- sample(x1, replace = TRUE)
    b2 <- sample(x2, replace = TRUE)
    hedges_g(c(b1, b2), c(rep(group1, length(b1)), rep(group2, length(b2))), group1, group2)
  })
  finite_boot <- boot_g[is.finite(boot_g)]
  stopifnot(length(finite_boot) > 1000L)
  perm <- permutation_p(values, groups, group1, group2, draws, seed + 1L)
  loo <- vapply(seq_along(values), function(i) hedges_g(values[-i], groups[-i], group1, group2), numeric(1L))
  fit_values <- list(or = NA_real_, low = NA_real_, high = NA_real_, p = NA_real_)
  if (firth) {
    predictor <- safe_z(values)
    fit <- logistf(I(groups == group1) ~ predictor, firth = TRUE, pl = TRUE)
    fit_values <- list(
      or = exp(unname(coef(fit)["predictor"])),
      low = exp(unname(fit$ci.lower["predictor"])),
      high = exp(unname(fit$ci.upper["predictor"])),
      p = unname(fit$prob["predictor"])
    )
  }
  summary <- data.table(
    dataset = dataset, endpoint = endpoint, analysis_label = ifelse(endpoint == "PATHOLOGICAL_RESPONSE", "SECONDARY_EXPLORATORY_PATHOLOGICAL_RESPONSE", "PRESPECIFIED"),
    metric = metric, group1 = group1, group2 = group2, n = length(values), n_group1 = length(x1), n_group2 = length(x2),
    group1_mean = mean(x1), group2_mean = mean(x2), group1_median = median(x1), group2_median = median(x2),
    group1_minus_group2 = mean(x1) - mean(x2), hedges_g = observed_g,
    hedges_g_ci_low = unname(quantile(finite_boot, 0.025)), hedges_g_ci_high = unname(quantile(finite_boot, 0.975)),
    bootstrap_draws_requested = draws, bootstrap_draws_finite = length(finite_boot),
    firth_or_per_sd = fit_values$or, firth_ci_low = fit_values$low, firth_ci_high = fit_values$high, firth_p = fit_values$p,
    permutation_p_two_sided = perm$p, permutation_method = perm$method, permutation_allocations_or_draws = perm$allocations,
    loo_hedges_g_min = min(loo, na.rm = TRUE), loo_hedges_g_max = max(loo, na.rm = TRUE),
    loo_positive_fraction = mean(loo > 0, na.rm = TRUE), loo_same_direction_fraction = mean(sign(loo) == sign(observed_g), na.rm = TRUE),
    biological_replicate = "patient", seed = seed
  )
  loo_table <- data.table(
    dataset = dataset, endpoint = endpoint, metric = metric, omitted_patient_id = ids,
    loo_hedges_g = loo, loo_positive = loo > 0, same_direction_as_full = sign(loo) == sign(observed_g)
  )
  list(summary = summary, loo = loo_table)
}

residual_effect_summary <- function(raw_values, broad, groups, ids, dataset, endpoint, metric, group1, group2, seed, firth = TRUE, draws = 20000L) {
  residual <- safe_z(residuals(lm(raw_values ~ broad)))
  result <- effect_summary(residual, groups, ids, dataset, endpoint, metric, group1, group2, seed, firth, draws)
  loo <- vapply(seq_along(raw_values), function(i) {
    residual_i <- safe_z(residuals(lm(raw_values[-i] ~ broad[-i])))
    hedges_g(residual_i, groups[-i], group1, group2)
  }, numeric(1L))
  result$summary[, `:=`(
    loo_hedges_g_min = min(loo, na.rm = TRUE), loo_hedges_g_max = max(loo, na.rm = TRUE),
    loo_positive_fraction = mean(loo > 0, na.rm = TRUE),
    loo_same_direction_fraction = mean(sign(loo) == sign(hedges_g), na.rm = TRUE)
  )]
  result$loo[, `:=`(
    loo_hedges_g = loo, loo_positive = loo > 0,
    same_direction_as_full = sign(loo) == sign(result$summary$hedges_g)
  )]
  result
}

fmt <- function(x, digits = 3L) ifelse(is.finite(x), formatC(x, format = "f", digits = digits), "NA")

# GSE302495 join and primary/secondary analyses.
g302_scores <- copy(g302_object$molecular_scores)
g302_map <- fread(g302_map_path)
setkey(g302_scores, patient_id)
setkey(g302_map, patient_id)
g302 <- g302_map[g302_scores]
stopifnot(nrow(g302) == 38L, !anyNA(g302$response_group), sum(g302$response_group == "R") == 12L, sum(g302$response_group == "NR") == 26L)
g302[, recist_group := fifelse(response_group == "R", "Responder", "Nonresponder")]

recist_specs <- list(
  list(metric = "PTPN22_expression", seed = 202608301L, residual = FALSE, raw = "PTPN22_expression"),
  list(metric = "AF35_score", seed = 202608302L, residual = FALSE, raw = "AF35_score"),
  list(metric = "broad_immune_score", seed = 202608303L, residual = FALSE, raw = "broad_immune_score"),
  list(metric = "PTPN22_broad_immune_residual", seed = 202608304L, residual = TRUE, raw = "PTPN22_expression"),
  list(metric = "AF35_broad_immune_residual", seed = 202608305L, residual = TRUE, raw = "AF35_score")
)
recist_results <- vector("list", length(recist_specs))
for (i in seq_along(recist_specs)) {
  spec <- recist_specs[[i]]
  if (spec$residual) {
    recist_results[[i]] <- residual_effect_summary(g302[[spec$raw]], g302$broad_immune_score, g302$recist_group, g302$patient_id, "GSE302495", "RECIST", spec$metric, "Responder", "Nonresponder", spec$seed)
  } else {
    recist_results[[i]] <- effect_summary(g302[[spec$metric]], g302$recist_group, g302$patient_id, "GSE302495", "RECIST", spec$metric, "Responder", "Nonresponder", spec$seed)
  }
}
recist_effects <- rbindlist(lapply(recist_results, `[[`, "summary"))
recist_loo <- rbindlist(lapply(recist_results, `[[`, "loo"))
fwrite(recist_effects, file.path(table_dir, "GSE302495_ALL_RECIST_CONTINUOUS_EFFECTS.tsv"), sep = "\t", quote = TRUE, na = "NA")
fwrite(recist_effects[metric == "PTPN22_expression"], file.path(table_dir, "GSE302495_PTPN22_RECIST_ASSOCIATION.tsv"), sep = "\t", quote = TRUE, na = "NA")
fwrite(recist_effects[metric == "AF35_score"], file.path(table_dir, "GSE302495_AF35_RECIST_ASSOCIATION.tsv"), sep = "\t", quote = TRUE, na = "NA")
fwrite(recist_loo, file.path(table_dir, "GSE302495_RECIST_LOO_AUDIT.tsv"), sep = "\t", quote = TRUE, na = "NA")

path_data <- g302[pathological_response_group %chin% c("MPR", "NON_MPR")]
stopifnot(nrow(path_data) == 22L, sum(path_data$pathological_response_group == "MPR") == 7L, sum(path_data$pathological_response_group == "NON_MPR") == 15L)
path_specs <- list(
  list(metric = "PTPN22_expression", seed = 202608311L),
  list(metric = "AF35_score", seed = 202608312L),
  list(metric = "broad_immune_score", seed = 202608313L),
  list(metric = "PTPN22_broad_immune_residual", seed = 202608314L),
  list(metric = "AF35_broad_immune_residual", seed = 202608315L)
)
path_results <- lapply(path_specs, function(spec) {
  effect_summary(path_data[[spec$metric]], path_data$pathological_response_group, path_data$patient_id, "GSE302495", "PATHOLOGICAL_RESPONSE", spec$metric, "MPR", "NON_MPR", spec$seed)
})
path_effects <- rbindlist(lapply(path_results, `[[`, "summary"))
path_loo <- rbindlist(lapply(path_results, `[[`, "loo"))
fwrite(path_effects, file.path(table_dir, "GSE302495_PATHOLOGICAL_RESPONSE_ASSOCIATION.tsv"), sep = "\t", quote = TRUE, na = "NA")
fwrite(path_loo, file.path(table_dir, "GSE302495_PATHOLOGICAL_RESPONSE_LOO_AUDIT.tsv"), sep = "\t", quote = TRUE, na = "NA")

rho_ptpn22 <- cor(g302$PTPN22_expression, g302$broad_immune_score, method = "spearman")
rho_af35 <- cor(g302$AF35_score, g302$broad_immune_score, method = "spearman")
p <- recist_effects[metric == "PTPN22_expression"]
a <- recist_effects[metric == "AF35_score"]
b <- recist_effects[metric == "broad_immune_score"]
pr <- recist_effects[metric == "PTPN22_broad_immune_residual"]
ar <- recist_effects[metric == "AF35_broad_immune_residual"]

broad_report <- c(
  "# GSE302495 broad-immune falsification/context report", "",
  "Population: 38 pretreatment tumor patients (12 RECIST responders, 26 nonresponders)  ",
  "Biological replicate: patient  ",
  "Interpretation: context/falsification only; residualization does not prove mechanistic independence.", "",
  sprintf("- Broad immune response effect: Hedges' g `%s` (bootstrap 95%% CI `%s` to `%s`), permutation P=`%s`, LOO-positive fraction `%s`.", fmt(b$hedges_g), fmt(b$hedges_g_ci_low), fmt(b$hedges_g_ci_high), fmt(b$permutation_p_two_sided, 4L), fmt(b$loo_positive_fraction, 2L)),
  sprintf("- PTPN22 versus broad immune Spearman rho: `%s`.", fmt(rho_ptpn22)),
  sprintf("- AF35 versus broad immune Spearman rho: `%s`.", fmt(rho_af35)),
  sprintf("- PTPN22 residual effect: Hedges' g `%s` (`%s` to `%s`), Firth OR/SD `%s` (`%s` to `%s`), permutation P=`%s`, LOO-positive fraction `%s`.", fmt(pr$hedges_g), fmt(pr$hedges_g_ci_low), fmt(pr$hedges_g_ci_high), fmt(pr$firth_or_per_sd), fmt(pr$firth_ci_low), fmt(pr$firth_ci_high), fmt(pr$permutation_p_two_sided, 4L), fmt(pr$loo_positive_fraction, 2L)),
  sprintf("- AF35 residual effect: Hedges' g `%s` (`%s` to `%s`), Firth OR/SD `%s` (`%s` to `%s`), permutation P=`%s`, LOO-positive fraction `%s`.", fmt(ar$hedges_g), fmt(ar$hedges_g_ci_low), fmt(ar$hedges_g_ci_high), fmt(ar$firth_or_per_sd), fmt(ar$firth_ci_low), fmt(ar$firth_ci_high), fmt(ar$permutation_p_two_sided, 4L), fmt(ar$loo_positive_fraction, 2L)), "",
  sprintf("Raw-to-residual PTPN22 g ratio: `%s`.", fmt(pr$hedges_g / p$hedges_g)),
  "No clinical covariate, cutoff, feature selection, cross-cohort normalization or rescue model was added."
)
writeLines(broad_report, file.path(result_dir, "GSE302495_BROAD_IMMUNE_FALSIFICATION.md"), useBytes = TRUE)

strong_attenuation <- p$hedges_g > 0 && (pr$hedges_g <= 0 || pr$hedges_g < 0.5 * p$hedges_g)
g302_gate <- if (p$hedges_g < 0) {
  "PUBLIC_PRETREATMENT_REPLICATION_CONFLICTING_OR_NULL"
} else if (p$hedges_g > 0 && p$loo_positive_fraction < 0.9) {
  "PUBLIC_PRETREATMENT_REPLICATION_PARTIAL"
} else if (p$hedges_g > 0 && p$loo_positive_fraction >= 0.9 && a$hedges_g <= 0) {
  "PUBLIC_PRETREATMENT_REPLICATION_PARTIAL"
} else if (strong_attenuation) {
  "PUBLIC_PRETREATMENT_REPLICATION_PARTIAL"
} else if (p$hedges_g > 0 && p$loo_positive_fraction >= 0.9 && a$hedges_g > 0) {
  "PUBLIC_PRETREATMENT_REPLICATION_SUPPORTED"
} else {
  "PUBLIC_PRETREATMENT_REPLICATION_CONFLICTING_OR_NULL"
}

gate_report <- c(
  "# GSE302495 public pretreatment replication gate", "",
  sprintf("Decision: `%s`", g302_gate), "",
  "## Frozen gate inputs", "",
  "- Exact map: PASS; 38 baseline patients, 12 R and 26 NR.",
  sprintf("- PTPN22 direction: g `%s`; positive LOO fraction `%s`.", fmt(p$hedges_g), fmt(p$loo_positive_fraction, 2L)),
  sprintf("- AF35 direction: g `%s`; positive LOO fraction `%s`.", fmt(a$hedges_g), fmt(a$loo_positive_fraction, 2L)),
  sprintf("- PTPN22 residual g `%s`; raw-to-residual ratio `%s`; frozen strong-attenuation flag `%s`.", fmt(pr$hedges_g), fmt(pr$hedges_g / p$hedges_g), ifelse(strong_attenuation, "YES", "NO")),
  "- Processing and gate rules were not changed after outcome loading.", "",
  "## Primary estimates", "",
  sprintf("PTPN22 R–NR difference `%s`; Hedges' g `%s` (95%% CI `%s` to `%s`); Firth OR/SD `%s` (`%s` to `%s`); permutation P=`%s`; LOO g `%s` to `%s`.", fmt(p$group1_minus_group2), fmt(p$hedges_g), fmt(p$hedges_g_ci_low), fmt(p$hedges_g_ci_high), fmt(p$firth_or_per_sd), fmt(p$firth_ci_low), fmt(p$firth_ci_high), fmt(p$permutation_p_two_sided, 4L), fmt(p$loo_hedges_g_min), fmt(p$loo_hedges_g_max)),
  sprintf("AF35 R–NR difference `%s`; Hedges' g `%s` (95%% CI `%s` to `%s`); Firth OR/SD `%s` (`%s` to `%s`); permutation P=`%s`; LOO g `%s` to `%s`.", fmt(a$group1_minus_group2), fmt(a$hedges_g), fmt(a$hedges_g_ci_low), fmt(a$hedges_g_ci_high), fmt(a$firth_or_per_sd), fmt(a$firth_ci_low), fmt(a$firth_ci_high), fmt(a$permutation_p_two_sided, 4L), fmt(a$loo_hedges_g_min), fmt(a$loo_hedges_g_max)), "",
  "No post-gate rescue analysis was performed."
)
writeLines(gate_report, file.path(result_dir, "GSE302495_PUBLIC_REPLICATION_GATE.md"), useBytes = TRUE)

# GSE287319 exploratory compatibility.
g287_scores <- copy(g287_object$molecular_scores)
g287_map <- fread(g287_map_path)
setkey(g287_scores, patient_id)
setkey(g287_map, patient_id)
g287_all <- g287_map[g287_scores]
stopifnot(nrow(g287_all) == 9L, sum(g287_all$endpoint_eligible) == 8L)
g287 <- g287_all[endpoint_eligible == TRUE]
stopifnot(nrow(g287) == 8L, sum(g287$endpoint_group == "REDUCTION") == 4L, sum(g287$endpoint_group == "PROGRESSION") == 4L)
# Refit residuals within the prespecified eight eligible patients without using group labels in the model.
g287[, PTPN22_broad_immune_residual := safe_z(residuals(lm(PTPN22_expression ~ broad_immune_score)))]
g287[, AF35_broad_immune_residual := safe_z(residuals(lm(AF35_score ~ broad_immune_score)))]
g287_specs <- list(
  list(metric = "PTPN22_expression", seed = 202608321L),
  list(metric = "AF35_score", seed = 202608322L),
  list(metric = "broad_immune_score", seed = 202608323L),
  list(metric = "PTPN22_broad_immune_residual", seed = 202608324L),
  list(metric = "AF35_broad_immune_residual", seed = 202608325L)
)
g287_results <- lapply(g287_specs, function(spec) {
  if (spec$metric == "PTPN22_broad_immune_residual") {
    residual_effect_summary(g287$PTPN22_expression, g287$broad_immune_score, g287$endpoint_group, g287$patient_id, "GSE287319", "NONIRRADIATED_LESION_DYNAMICS", spec$metric, "REDUCTION", "PROGRESSION", spec$seed, firth = FALSE)
  } else if (spec$metric == "AF35_broad_immune_residual") {
    residual_effect_summary(g287$AF35_score, g287$broad_immune_score, g287$endpoint_group, g287$patient_id, "GSE287319", "NONIRRADIATED_LESION_DYNAMICS", spec$metric, "REDUCTION", "PROGRESSION", spec$seed, firth = FALSE)
  } else {
    effect_summary(g287[[spec$metric]], g287$endpoint_group, g287$patient_id, "GSE287319", "NONIRRADIATED_LESION_DYNAMICS", spec$metric, "REDUCTION", "PROGRESSION", spec$seed, firth = FALSE)
  }
})
g287_effects <- rbindlist(lapply(g287_results, `[[`, "summary"))
g287_loo <- rbindlist(lapply(g287_results, `[[`, "loo"))
fwrite(g287_effects, file.path(table_dir, "GSE287319_EXPLORATORY_CONTINUOUS_EFFECTS.tsv"), sep = "\t", quote = TRUE, na = "NA")
fwrite(g287_loo, file.path(table_dir, "GSE287319_EXPLORATORY_LOO_AUDIT.tsv"), sep = "\t", quote = TRUE, na = "NA")
g287_p <- g287_effects[metric == "PTPN22_expression"]
g287_a <- g287_effects[metric == "AF35_score"]
g287_b <- g287_effects[metric == "broad_immune_score"]
g287_decision <- if (g287_p$hedges_g > 0 && g287_a$hedges_g > 0) {
  "GSE287319_EXPLORATORY_CONCORDANT"
} else if (g287_p$hedges_g < 0 && g287_a$hedges_g < 0) {
  "GSE287319_EXPLORATORY_DISCORDANT"
} else {
  "GSE287319_EXPLORATORY_INCONCLUSIVE"
}
g287_report <- c(
  "# GSE287319 exploratory compatibility report", "",
  sprintf("Decision: `%s`", g287_decision), "",
  "Role: `EXPLORATORY_ADVANCED_HCC_COMPATIBILITY`  ",
  "Endpoint: nonirradiated-lesion reduction versus progression; not RECIST.  ",
  "Population: eight classified pretreatment patients (4/4); one unclassifiable baseline patient and three post-treatment profiles excluded.  ",
  "Biological replicate: patient.", "",
  sprintf("- PTPN22: reduction mean `%s`, progression mean `%s`, difference `%s`, Hedges' g `%s` with unstable small-n bootstrap 95%% CI `%s` to `%s`, exact 70-allocation P=`%s`, LOO g `%s` to `%s`, positive fraction `%s`.", fmt(g287_p$group1_mean), fmt(g287_p$group2_mean), fmt(g287_p$group1_minus_group2), fmt(g287_p$hedges_g), fmt(g287_p$hedges_g_ci_low), fmt(g287_p$hedges_g_ci_high), fmt(g287_p$permutation_p_two_sided, 4L), fmt(g287_p$loo_hedges_g_min), fmt(g287_p$loo_hedges_g_max), fmt(g287_p$loo_positive_fraction, 2L)),
  sprintf("- AF35: reduction mean `%s`, progression mean `%s`, difference `%s`, Hedges' g `%s` with unstable small-n bootstrap 95%% CI `%s` to `%s`, exact P=`%s`, LOO positive fraction `%s`.", fmt(g287_a$group1_mean), fmt(g287_a$group2_mean), fmt(g287_a$group1_minus_group2), fmt(g287_a$hedges_g), fmt(g287_a$hedges_g_ci_low), fmt(g287_a$hedges_g_ci_high), fmt(g287_a$permutation_p_two_sided, 4L), fmt(g287_a$loo_positive_fraction, 2L)),
  sprintf("- Broad immune descriptive context: Hedges' g `%s` (`%s` to `%s`), exact P=`%s`.", fmt(g287_b$hedges_g), fmt(g287_b$hedges_g_ci_low), fmt(g287_b$hedges_g_ci_high), fmt(g287_b$permutation_p_two_sided, 4L)), "",
  "Residual effects are retained in the continuous-effects table as response-blind context. No Firth model, multivariable model, survival analysis, subgrouping, pooling or gate rescue was performed. The n=8 interval is intrinsically unstable and cannot upgrade the GSE302495 gate."
)
writeLines(g287_report, file.path(result_dir, "GSE287319_EXPLORATORY_COMPATIBILITY_REPORT.md"), useBytes = TRUE)

cat("GSE302495_GATE=", g302_gate, "\n", sep = "")
cat("GSE302495_PTPN22_G=", p$hedges_g, " P=", p$permutation_p_two_sided, " LOO_POS=", p$loo_positive_fraction, "\n", sep = "")
cat("GSE302495_AF35_G=", a$hedges_g, " P=", a$permutation_p_two_sided, " LOO_POS=", a$loo_positive_fraction, "\n", sep = "")
cat("GSE302495_BROAD_G=", b$hedges_g, " P=", b$permutation_p_two_sided, "\n", sep = "")
cat("GSE302495_PTPN22_RESID_G=", pr$hedges_g, " P=", pr$permutation_p_two_sided, "\n", sep = "")
cat("GSE302495_AF35_RESID_G=", ar$hedges_g, " P=", ar$permutation_p_two_sided, "\n", sep = "")
cat("GSE287319_DECISION=", g287_decision, "\n", sep = "")
cat("GSE287319_PTPN22_G=", g287_p$hedges_g, " P=", g287_p$permutation_p_two_sided, " LOO_POS=", g287_p$loo_positive_fraction, "\n", sep = "")
cat("GSE287319_AF35_G=", g287_a$hedges_g, " P=", g287_a$permutation_p_two_sided, " LOO_POS=", g287_a$loo_positive_fraction, "\n", sep = "")
