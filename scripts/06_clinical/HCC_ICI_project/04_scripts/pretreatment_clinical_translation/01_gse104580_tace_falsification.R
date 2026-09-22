options(stringsAsFactors = FALSE, width = 240)
invisible(try(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"), silent = TRUE))

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
  library(logistf)
})

set.seed(20260829L)
project_dir <- normalizePath("HCC_ICI_project", winslash = "/", mustWork = TRUE)
raw_dir <- file.path(project_dir, "02_raw_data", "pretreatment_clinical_translation", "GSE104580")
metadata_dir <- file.path(project_dir, "01_metadata", "pretreatment_clinical_translation")
table_dir <- file.path(project_dir, "05_results", "tables", "pretreatment_clinical_translation")
result_dir <- file.path(project_dir, "05_results", "pretreatment_clinical_translation")
log_dir <- file.path(project_dir, "09_logs", "pretreatment_clinical_translation")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

matrix_path <- file.path(raw_dir, "GSE104580_series_matrix.txt.gz")
annotation_path <- file.path(raw_dir, "GPL570.annot.gz")
map_path <- file.path(metadata_dir, "GSE104580_PATIENT_RESPONSE_MAP.tsv")
axis_path <- file.path(project_dir, "01_metadata", "clinical_translation", "PTPN22_CLINICAL_AXES_FROZEN.tsv")
comparator_path <- file.path(project_dir, "01_metadata", "PHASE3_PTPN22_SPECIFICITY_COMPARATOR_PROGRAMS.tsv")
contract_path <- file.path(metadata_dir, "PRETREATMENT_CLINICAL_ANALYSIS_CONTRACT.md")
stopifnot(all(file.exists(c(matrix_path, annotation_path, map_path, axis_path, comparator_path, contract_path))))

unquote <- function(x) sub('^"(.*)"$', '\\1', x)
read_geo_matrix <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE, encoding = "UTF-8")
  meta <- function(prefix) {
    hit <- lines[startsWith(lines, prefix)]
    stopifnot(length(hit) == 1L)
    unquote(strsplit(hit, "\t", fixed = TRUE)[[1L]][-1L])
  }
  begin <- match("!series_matrix_table_begin", lines)
  end <- match("!series_matrix_table_end", lines)
  stopifnot(is.finite(begin), is.finite(end), end > begin + 1L)
  list(
    titles = meta("!Sample_title"),
    gsm = meta("!Sample_geo_accession"),
    matrix = fread(text = paste(lines[(begin + 1L):(end - 1L)], collapse = "\n"), sep = "\t", check.names = FALSE)
  )
}
read_geo_platform <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE, encoding = "UTF-8")
  begin <- match("!platform_table_begin", lines)
  end <- match("!platform_table_end", lines)
  stopifnot(is.finite(begin), is.finite(end), end > begin + 1L)
  fread(text = paste(lines[(begin + 1L):(end - 1L)], collapse = "\n"), sep = "\t", check.names = FALSE, fill = TRUE)
}
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
hedges_g <- function(values, groups) {
  r <- values[groups == "Responder"]
  nr <- values[groups == "Nonresponder"]
  pooled_sd <- sqrt(((length(r) - 1L) * var(r) + (length(nr) - 1L) * var(nr)) / (length(values) - 2L))
  if (!is.finite(pooled_sd) || pooled_sd == 0) return(NA_real_)
  correction <- 1 - 3 / (4 * (length(values) - 2L) - 1)
  correction * (mean(r) - mean(nr)) / pooled_sd
}
label_permutation_p <- function(values, groups, draws = 20000L) {
  observed <- mean(values[groups == "Responder"]) - mean(values[groups == "Nonresponder"])
  n_r <- sum(groups == "Responder")
  perm <- replicate(draws, {
    idx <- sample.int(length(values), n_r, replace = FALSE)
    mean(values[idx]) - mean(values[-idx])
  })
  (1 + sum(abs(perm) >= abs(observed) - 1e-12)) / (draws + 1)
}
effect_summary <- function(values, groups, ids, draws = 20000L) {
  keep <- is.finite(values) & groups %chin% c("Responder", "Nonresponder")
  values <- values[keep]
  groups <- groups[keep]
  ids <- ids[keep]
  r <- values[groups == "Responder"]
  nr <- values[groups == "Nonresponder"]
  observed_g <- hedges_g(values, groups)
  boot_g <- replicate(draws, {
    rb <- sample(r, replace = TRUE)
    nrb <- sample(nr, replace = TRUE)
    hedges_g(c(rb, nrb), c(rep("Responder", length(rb)), rep("Nonresponder", length(nrb))))
  })
  boot_g <- boot_g[is.finite(boot_g)]
  predictor <- safe_z(values)
  fit <- logistf(I(groups == "Responder") ~ predictor, firth = TRUE, pl = TRUE)
  loo <- vapply(seq_along(values), function(i) hedges_g(values[-i], groups[-i]), numeric(1L))
  full_sign <- sign(observed_g)
  data.table(
    n = length(values), n_responder = length(r), n_nonresponder = length(nr),
    responder_mean = mean(r), nonresponder_mean = mean(nr), responder_minus_nonresponder = mean(r) - mean(nr),
    hedges_g = observed_g, hedges_g_ci_low = unname(quantile(boot_g, 0.025)),
    hedges_g_ci_high = unname(quantile(boot_g, 0.975)), bootstrap_draws = draws,
    firth_or_per_sd = exp(unname(coef(fit)["predictor"])),
    firth_ci_low = exp(unname(fit$ci.lower["predictor"])),
    firth_ci_high = exp(unname(fit$ci.upper["predictor"])),
    firth_p = unname(fit$prob["predictor"]),
    permutation_p_two_sided = label_permutation_p(values, groups, draws), permutation_draws = draws,
    loo_hedges_g_min = min(loo, na.rm = TRUE), loo_hedges_g_max = max(loo, na.rm = TRUE),
    loo_same_direction_fraction = mean(sign(loo) == full_sign, na.rm = TRUE), biological_replicate = "patient", seed = 20260829L
  )
}

# Response-blind platform mapping and probe collapse occur before the response map is read.
geo <- read_geo_matrix(matrix_path)
platform <- read_geo_platform(annotation_path)
stopifnot(length(geo$titles) == 147L, ncol(geo$matrix) == 148L)
expr <- as.matrix(geo$matrix[, -1L, with = FALSE])
storage.mode(expr) <- "double"
rownames(expr) <- as.character(geo$matrix[[1L]])
colnames(expr) <- geo$titles
stopifnot(all(is.finite(expr)), max(expr) < 25, quantile(expr, 0.99) < 15)

symbols <- trimws(as.character(platform[["Gene symbol"]]))
unambiguous <- nzchar(symbols) & symbols != "---" & !grepl("///", symbols, fixed = TRUE)
probe_symbol <- data.table(probe_id = as.character(platform$ID[unambiguous]), gene_symbol = symbols[unambiguous])
probe_symbol <- probe_symbol[probe_id %chin% rownames(expr)]
stopifnot(!anyDuplicated(probe_symbol$probe_id))

gene_expr <- rowsum(expr[probe_symbol$probe_id, , drop = FALSE], group = probe_symbol$gene_symbol, reorder = FALSE)
probe_n <- table(factor(probe_symbol$gene_symbol, levels = rownames(gene_expr)))
gene_expr <- gene_expr / as.numeric(probe_n)
stopifnot("PTPN22" %chin% rownames(gene_expr), length(probe_symbol[gene_symbol == "PTPN22", probe_id]) == 4L)

axes <- fread(axis_path)
af_genes <- strsplit(axes[component == "PTPN22_activation_feedback_35", exact_features], ";", fixed = TRUE)[[1L]]
comparators <- fread(comparator_path)
broad_genes <- comparators[program == "broad_immune_abundance" & !startsWith(gene, "DERIVED:"), gene]
stopifnot(length(af_genes) == 35L, length(broad_genes) == 11L)

af <- score_program(gene_expr, af_genes, 28L, 0.8)
broad <- score_program(gene_expr, broad_genes, 11L, 0.8)
ptpn22 <- safe_z(gene_expr["PTPN22", ])
stopifnot(af$eligible, broad$eligible, all(is.finite(ptpn22)), all(is.finite(af$score)), all(is.finite(broad$score)))

scores <- data.table(
  sample_title = colnames(gene_expr), PTPN22_z = ptpn22,
  AF_program_z = af$score, broad_immune_z = broad$score
)
scores[, PTPN22_broad_residual_z := safe_z(residuals(lm(PTPN22_z ~ broad_immune_z, data = scores)))]
scores[, AF_broad_residual_z := safe_z(residuals(lm(AF_program_z ~ broad_immune_z, data = scores)))]

sample_map <- fread(map_path)
scores <- merge(scores, sample_map, by = "sample_title", all.x = TRUE, sort = FALSE)
stopifnot(nrow(scores) == 147L, !anyNA(scores$response), !anyDuplicated(scores$patient_id))

metric_map <- c(
  PTPN22_z = "PTPN22",
  AF_program_z = "PTPN22_activation_feedback_35",
  broad_immune_z = "broad_immune_abundance",
  PTPN22_broad_residual_z = "PTPN22_residual_after_broad_immune",
  AF_broad_residual_z = "AF_program_residual_after_broad_immune"
)
effect_rows <- lapply(names(metric_map), function(metric) {
  out <- effect_summary(scores[[metric]], scores$response, scores$patient_id)
  out[, `:=`(dataset = "GSE104580", metric = unname(metric_map[[metric]]), metric_column = metric)]
  out
})
effects <- rbindlist(effect_rows, fill = TRUE)
setcolorder(effects, c("dataset", "metric", "metric_column", setdiff(names(effects), c("dataset", "metric", "metric_column"))))

correlations <- data.table(
  comparison = c("PTPN22_vs_broad_immune", "AF_program_vs_broad_immune", "PTPN22_vs_AF_program"),
  spearman_rho = c(
    cor(scores$PTPN22_z, scores$broad_immune_z, method = "spearman"),
    cor(scores$AF_program_z, scores$broad_immune_z, method = "spearman"),
    cor(scores$PTPN22_z, scores$AF_program_z, method = "spearman")
  ),
  biological_replicate = "patient"
)

coverage <- data.table(
  component = c("PTPN22", "PTPN22_activation_feedback_35", "broad_immune_abundance"),
  requested_n = c(1L, af$requested_n, broad$requested_n),
  present_n = c(1L, af$present_n, broad$present_n),
  variable_n = c(1L, af$variable_n, broad$variable_n),
  missing_genes = c("", paste(af$missing, collapse = ";"), paste(broad$missing, collapse = ";")),
  status = c("AVAILABLE_4_PROBE_RESPONSE_BLIND_MEAN", if (af$eligible) "SCORE_AVAILABLE" else "SCORE_UNAVAILABLE", if (broad$eligible) "SCORE_AVAILABLE" else "SCORE_UNAVAILABLE")
)

fwrite(scores, file.path(table_dir, "GSE104580_TACE_PATIENT_SCORES.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(effects, file.path(table_dir, "GSE104580_TACE_CONTINUOUS_EFFECTS.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(correlations, file.path(table_dir, "GSE104580_TACE_BROAD_IMMUNE_CORRELATIONS.tsv"), sep = "\t", quote = FALSE, na = "NA")
fwrite(coverage, file.path(table_dir, "GSE104580_TACE_PROGRAM_COVERAGE.tsv"), sep = "\t", quote = FALSE, na = "NA")

fmt <- function(x, d = 3L) ifelse(is.finite(x), formatC(x, format = "f", digits = d), "NA")
row_for <- function(metric_name) effects[metric == metric_name][1L]
p <- row_for("PTPN22")
a <- row_for("PTPN22_activation_feedback_35")
b <- row_for("broad_immune_abundance")
pr <- row_for("PTPN22_residual_after_broad_immune")
ar <- row_for("AF_program_residual_after_broad_immune")
report <- c(
  "# GSE104580 TACE context falsification",
  "",
  "Date: 2026-08-29  ",
  "Role: non-ICI treatment-response negative control  ",
  "Biological replicate: patient",
  "",
  "## Cohort and response-blind preprocessing",
  "",
  "The official GEO series matrix contains 147 unique pretreatment HCC tumor biopsies: 81 TACE responders and 66 nonresponders. All sample titles, GSMs and title-derived patient identifiers are unique. No sample or probe was counted as a patient.",
  "",
  sprintf("GPL570 maps PTPN22 unambiguously to four probes (`%s`). Only probes with one unambiguous gene symbol were retained, and multiple probes were collapsed by their all-patient arithmetic mean before response labels were joined.", paste(probe_symbol[gene_symbol == "PTPN22", probe_id], collapse = ";")),
  sprintf("The frozen 35-gene program passed coverage (%d/35 present; %d variable). The exact broad-immune comparator passed coverage (%d/11 present; %d variable).", af$present_n, af$variable_n, broad$present_n, broad$variable_n),
  "",
  "The released values were already on a bounded log-like microarray scale (overall range -2.403 to 18.639; 99th percentile 13.191); no outcome-informed transformation or renormalization was applied. Scores use the frozen gene-wise-z, mean, re-z rule.",
  "",
  "## Patient-level continuous effects",
  "",
  sprintf("- PTPN22: responder-minus-nonresponder %s SD units; Hedges' g %s (patient bootstrap 95%% CI %s to %s); Firth OR per SD %s (%s to %s); permutation P=%s; LOO g range %s to %s, same-direction fraction %s.", fmt(p$responder_minus_nonresponder), fmt(p$hedges_g), fmt(p$hedges_g_ci_low), fmt(p$hedges_g_ci_high), fmt(p$firth_or_per_sd), fmt(p$firth_ci_low), fmt(p$firth_ci_high), fmt(p$permutation_p_two_sided, 4L), fmt(p$loo_hedges_g_min), fmt(p$loo_hedges_g_max), fmt(p$loo_same_direction_fraction, 2L)),
  sprintf("- Frozen activation-feedback program: difference %s; Hedges' g %s (%s to %s); Firth OR %s (%s to %s); permutation P=%s; LOO same-direction fraction %s.", fmt(a$responder_minus_nonresponder), fmt(a$hedges_g), fmt(a$hedges_g_ci_low), fmt(a$hedges_g_ci_high), fmt(a$firth_or_per_sd), fmt(a$firth_ci_low), fmt(a$firth_ci_high), fmt(a$permutation_p_two_sided, 4L), fmt(a$loo_same_direction_fraction, 2L)),
  sprintf("- Broad immune comparator: difference %s; Hedges' g %s (%s to %s); permutation P=%s.", fmt(b$responder_minus_nonresponder), fmt(b$hedges_g), fmt(b$hedges_g_ci_low), fmt(b$hedges_g_ci_high), fmt(b$permutation_p_two_sided, 4L)),
  "",
  "## Broad-immune conditional context",
  "",
  sprintf("PTPN22 correlated with broad immune at Spearman rho %s; the activation-feedback program correlated with broad immune at rho %s.", fmt(correlations[comparison == "PTPN22_vs_broad_immune", spearman_rho]), fmt(correlations[comparison == "AF_program_vs_broad_immune", spearman_rho])),
  sprintf("After response-blind broad-immune residualization, PTPN22 Hedges' g was %s (%s to %s; permutation P=%s), and the program residual Hedges' g was %s (%s to %s; P=%s).", fmt(pr$hedges_g), fmt(pr$hedges_g_ci_low), fmt(pr$hedges_g_ci_high), fmt(pr$permutation_p_two_sided, 4L), fmt(ar$hedges_g), fmt(ar$hedges_g_ci_low), fmt(ar$hedges_g_ci_high), fmt(ar$permutation_p_two_sided, 4L)),
  "",
  "## Decision",
  "",
  "The final categorical TACE decision is assigned only after the independent ERP117672 pretreatment ICI direction is known. GSE104580 is not pooled with ICI data and is not used to build a TACE classifier.",
  "",
  "No cutpoint, feature selection, covariate search, model tuning or treatment-interaction claim was used."
)
writeLines(report, file.path(result_dir, "GSE104580_TACE_CONTEXT_FALSIFICATION.md"), useBytes = TRUE)

log_path <- file.path(log_dir, "02_gse104580_tace_falsification.log")
log_con <- file(log_path, open = "wt", encoding = "UTF-8")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({ sink(type = "message"); sink(type = "output"); close(log_con) }, add = TRUE)
cat("GSE104580 TACE FALSIFICATION\n")
cat("Seed: 20260829; bootstrap/permutation draws: 20000\n")
cat("Contract SHA256:", digest(file = contract_path, algo = "sha256", serialize = FALSE), "\n")
cat("Matrix SHA256:", digest(file = matrix_path, algo = "sha256", serialize = FALSE), "\n")
cat("GPL570 annotation SHA256:", digest(file = annotation_path, algo = "sha256", serialize = FALSE), "\n")
cat("Matrix dimensions:", nrow(expr), "probes x", ncol(expr), "patients\n")
print(coverage)
print(effects)
print(correlations)
cat("SESSION_INFO\n")
print(sessionInfo())
