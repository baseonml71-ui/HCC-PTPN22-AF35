options(stringsAsFactors = FALSE, width = 240)

source("HCC_ICI_project/04_scripts/config.R")
set.seed(ANALYSIS_SEED)

project_dir <- normalizePath("HCC_ICI_project", winslash = "/", mustWork = TRUE)
table_dir <- file.path(project_dir, "05_results", "tables")
stats_dir <- file.path(project_dir, "05_results", "statistics")
log_file <- file.path(project_dir, "09_logs", "target_ranking.log")

patient_expression <- read.delim(file.path(table_dir, "GSE206325_candidate_expression_by_patient.tsv"), check.names = FALSE)
tcr_metrics <- read.delim(file.path(table_dir, "GSE206325_TCR_repertoire_metrics_by_sample.tsv"), check.names = FALSE)
tcr_persistence <- read.delim(file.path(table_dir, "GSE206325_TCR_persistence_by_patient.tsv"), check.names = FALSE)

cliffs_delta <- function(r, nr) mean(sign(outer(r, nr, "-")))

response_loo_rows <- list()
i <- 0L
for (candidate in CANDIDATES) {
  x <- patient_expression[patient_expression$candidate == candidate, ]
  for (metric in c("mean_normalized", "detected_fraction")) {
    for (excluded_patient in x$patient_ID) {
      y <- x[x$patient_ID != excluded_patient, ]
      r <- y[[metric]][y$response == "anti-PD1_R"]
      nr <- y[[metric]][y$response == "anti-PD1_NR"]
      i <- i + 1L
      response_loo_rows[[i]] <- data.frame(
        candidate = candidate,
        metric = metric,
        excluded_patient = excluded_patient,
        median_difference_r_minus_nr = median(r) - median(nr),
        cliffs_delta = cliffs_delta(r, nr),
        stringsAsFactors = FALSE
      )
    }
  }
}
response_loo <- do.call(rbind, response_loo_rows)

response_loo_summary <- do.call(rbind, lapply(split(response_loo, interaction(response_loo$candidate, response_loo$metric)), function(x) {
  data.frame(
    candidate = x$candidate[1L],
    metric = x$metric[1L],
    n_loo_runs = nrow(x),
    minimum_median_difference = min(x$median_difference_r_minus_nr),
    maximum_median_difference = max(x$median_difference_r_minus_nr),
    fraction_same_sign_as_full = mean(sign(x$median_difference_r_minus_nr) == sign(median(x$median_difference_r_minus_nr))),
    minimum_cliffs_delta = min(x$cliffs_delta),
    maximum_cliffs_delta = max(x$cliffs_delta),
    stringsAsFactors = FALSE
  )
}))

bootstrap_spearman <- function(x, y, n_boot = 5000L) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  n <- length(x)
  if (n < 4L) return(c(n = n, rho = NA, lower = NA, upper = NA, p = NA))
  rho <- suppressWarnings(stats::cor(x, y, method = "spearman"))
  draws <- replicate(n_boot, {
    idx <- sample.int(n, n, replace = TRUE)
    suppressWarnings(stats::cor(x[idx], y[idx], method = "spearman"))
  })
  draws <- draws[is.finite(draws)]
  p <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE)$p.value)
  c(n = n, rho = rho, lower = stats::quantile(draws, 0.025, names = FALSE), upper = stats::quantile(draws, 0.975, names = FALSE), p = p)
}

patient_tcr <- aggregate(
  cbind(normalized_shannon_clonality, simpson_concentration, top_clone_proportion) ~ patient_ID,
  data = tcr_metrics[!is.na(tcr_metrics$response) & nzchar(tcr_metrics$response), ],
  FUN = mean
)
patient_persistence <- aggregate(
  cbind(jaccard_overlap, baseline_mass_persistent, post_mass_from_persistent) ~ patient_ID,
  data = tcr_persistence,
  FUN = mean
)

tcr_association_rows <- list()
i <- 0L
for (candidate in CANDIDATES) {
  candidate_data <- patient_expression[patient_expression$candidate == candidate, ]
  for (endpoint in c("sample_repertoire", "paired_persistence")) {
    y <- if (endpoint == "sample_repertoire") patient_tcr else patient_persistence
    merged <- merge(candidate_data, y, by = "patient_ID")
    metrics <- setdiff(names(y), "patient_ID")
    for (metric in metrics) {
      estimate <- bootstrap_spearman(merged$mean_normalized, merged[[metric]])
      i <- i + 1L
      tcr_association_rows[[i]] <- data.frame(
        dataset = "GSE206325",
        candidate = candidate,
        endpoint = endpoint,
        metric = metric,
        n_patients = estimate[["n"]],
        spearman_rho = estimate[["rho"]],
        bootstrap_ci_lower = estimate[["lower"]],
        bootstrap_ci_upper = estimate[["upper"]],
        spearman_p = estimate[["p"]],
        stringsAsFactors = FALSE
      )
    }
  }
}
tcr_associations <- do.call(rbind, tcr_association_rows)

grades <- data.frame(
  candidate = rep(CANDIDATES, each = 6L),
  dimension = rep(names(TARGET_WEIGHTS), times = length(CANDIDATES)),
  grade = c(
    4, 4, 1, 1, 4, 4,
    4, 4, 2, 2, 3, 2,
    1, 2, 3, 1, 4, 4
  ),
  evidence = c(
    "Discovery tumor: positive terminal-exhaustion and Tfh-like associations with patient-level bootstrap CIs above zero.",
    "Independent GSE326201 reproduces positive terminal-exhaustion, Tfh-like, Treg and proliferating associations.",
    "Responder mean difference is positive but uncertain; detection-rate effect is null and relevant author states are not response-discriminatory.",
    "Public SCE lacks cell-level clonotype; sample-level MAP4K1 associations with clonality are weak and directionally mixed.",
    "Published genetic, inhibitor and degrader perturbations support enhanced T-cell function after HPK1 loss/inhibition (PMID 32860752).",
    "HPK1 has genetic, inhibitor and degrader evidence and active clinical development; no docking was run.",
    "Discovery tumor shows positive effector, terminal-exhaustion, Tfh-like and Treg associations, several with CIs above zero.",
    "GSE326201 reproduces effector, terminal-exhaustion and Tfh-like directions with CIs above zero.",
    "Responder mean expression is higher with moderate Cliff effect, but bootstrap median-difference CI crosses zero and detection-rate evidence is weak.",
    "No cell-level fate link; however, patient-level PTPN22 expression is coherently associated with multiple repertoire clonality metrics.",
    "Knockout and small-molecule inhibition support antitumor activity, with chronic-stimulation exhaustion as a caveat (PMIDs 34283806, 38056892).",
    "A reported small-molecule inhibitor exists, but phosphatase selectivity/translation is less mature than kinase targets; no docking was run.",
    "Discovery inflammatory-myeloid association is weak and suppressive-TAM direction is inconsistent.",
    "Inflammatory association strengthens in GSE326201, while suppressive-TAM direction reverses relative to discovery.",
    "PIK3CG detection is lower in responders with a bootstrap CI below zero, consistent with a resistance-associated myeloid target.",
    "TCR evidence is indirect for a myeloid target and candidate/repertoire correlations are weak or mixed.",
    "PI3K-gamma inhibition has preclinical macrophage activation and checkpoint-combination support (PMID 33785652).",
    "Selective PI3K-gamma inhibitors provide strong tractability evidence; no docking was run."
  ),
  stringsAsFactors = FALSE
)
grades$weight <- TARGET_WEIGHTS[grades$dimension]
grades$weighted_points <- grades$grade / 4 * grades$weight

ranking <- aggregate(weighted_points ~ candidate, data = grades, sum)
ranking$biology_points_excluding_druggability <- vapply(ranking$candidate, function(candidate) {
  sum(grades$weighted_points[grades$candidate == candidate & grades$dimension != "druggability"])
}, numeric(1L))
ranking$gate3_ici_grade <- grades$grade[match(paste(ranking$candidate, "ici_response"), paste(grades$candidate, grades$dimension))]
ranking$gate3_sctcr_grade <- grades$grade[match(paste(ranking$candidate, "sctcr_clonal_fate"), paste(grades$candidate, grades$dimension))]
ranking$decision <- c(MAP4K1 = "NO-GO as lead; backup", PTPN22 = "CONDITIONAL GO; selected lead", PIK3CG = "NO-GO for current paper")[ranking$candidate]
ranking <- ranking[order(-ranking$weighted_points, -ranking$biology_points_excluding_druggability, ranking$candidate), ]
ranking$rank <- seq_len(nrow(ranking))
ranking <- ranking[c("rank", "candidate", "weighted_points", "biology_points_excluding_druggability", "gate3_ici_grade", "gate3_sctcr_grade", "decision")]

write.table(response_loo, file.path(stats_dir, "GSE206325_candidate_response_leave_one_patient_out.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(response_loo_summary, file.path(stats_dir, "GSE206325_candidate_response_LOO_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tcr_associations, file.path(stats_dir, "GSE206325_candidate_TCR_associations_with_CI.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(grades, file.path(table_dir, "TARGET_RANKING_EVIDENCE_GRADES_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ranking, file.path(table_dir, "TARGET_RANKING_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)
cat("RESPONSE_LOO_SUMMARY\n")
print(response_loo_summary)
cat("\nTCR_ASSOCIATIONS_WITH_CI\n")
print(tcr_associations)
cat("\nEVIDENCE_GRADES\n")
print(grades)
cat("\nTARGET_RANKING\n")
print(ranking)
cat("\nSESSION_INFO\n")
print(sessionInfo())
