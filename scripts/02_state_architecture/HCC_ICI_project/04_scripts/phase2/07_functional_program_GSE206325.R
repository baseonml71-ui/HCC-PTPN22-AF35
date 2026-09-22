options(stringsAsFactors = FALSE, width = 240)
invisible(Sys.setlocale("LC_ALL", "Chinese (Simplified)_China.utf8"))

suppressPackageStartupMessages({
  library(Matrix)
  library(edgeR)
  library(limma)
  library(fgsea)
  library(msigdbr)
})

source("HCC_ICI_project/04_scripts/config.R")
set.seed(ANALYSIS_SEED)

project_dir <- normalizePath("HCC_ICI_project", winslash = "/", mustWork = TRUE)
input_file <- file.path(project_dir, "03_processed_data", "phase2", "GSE206325_phase2_patient_state_pseudobulk.rds")
table_dir <- file.path(project_dir, "05_results", "tables", "phase2")
log_dir <- file.path(project_dir, "09_logs", "phase2")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

pb <- readRDS(input_file)
counts <- pb$counts
meta <- pb$metadata
rownames(meta) <- meta$pseudobulk_ID
if (!identical(colnames(counts), rownames(meta))) stop("Pseudobulk metadata order mismatch")

program_sets <- list(
  TCR_activation = c("CD3D", "CD3E", "CD3G", "CD247", "LCK", "FYN", "ZAP70", "LAT", "LCP2", "TRAT1"),
  NFAT = c("NFATC1", "NFATC2", "NFATC3", "RCAN1", "EGR1", "EGR2", "FOS", "JUN", "IL2"),
  AP1 = c("JUN", "JUNB", "JUND", "FOS", "FOSB", "FOSL1", "FOSL2", "ATF3", "BATF"),
  NFkB = c("NFKB1", "NFKB2", "REL", "RELA", "RELB", "NFKBIA", "TNFAIP3", "BCL3"),
  IFN_response = c("IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "STAT1", "IRF1", "IFI6", "IFI44L"),
  IFNG_TNF = c("IFNG", "TNF", "CXCL9", "CXCL10", "CXCL11", "CCL3", "CCL4", "CCL5", "TNFRSF9"),
  cytotoxicity = c("NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "CTSW", "CCL5", "KLRD1"),
  exhaustion = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "ENTPD1", "CXCL13", "LAYN", "BATF"),
  checkpoint = c("PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "CD274", "PDCD1LG2", "TNFRSF9"),
  proliferation = c("MKI67", "TOP2A", "STMN1", "TUBA1B", "TYMS", "PCNA", "UBE2C", "CENPF"),
  effector_differentiation = c("GZMK", "CCL5", "IFNG", "TNF", "PRF1", "GZMB", "FGFBP2", "TBX21", "EOMES"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "TNFRSF18", "IKZF2", "TIGIT", "BATF")
)

hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
msig <- rbind(hallmark[, c("gs_name", "gene_symbol")], reactome[, c("gs_name", "gene_symbol")])
msig_pathways <- split(msig$gene_symbol, msig$gs_name)

gene_results <- list()
program_results <- list()
score_rows <- list()
fgsea_results <- list()
gi <- pi <- si <- fi <- 0L

for (state_name in unique(meta$state)) {
  use <- which(meta$state == state_name & meta$n_cells >= 20L)
  state_meta <- meta[use, , drop = FALSE]
  state_counts <- counts[, use, drop = FALSE]
  if (nrow(state_meta) < 10L || length(unique(state_meta$response)) < 2L) next

  dge_all <- DGEList(state_counts)
  dge_all <- calcNormFactors(dge_all)
  logcpm_all <- cpm(dge_all, log = TRUE, prior.count = 1)
  if (!"PTPN22" %in% rownames(logcpm_all)) stop("PTPN22 missing")
  ptpn22 <- as.numeric(logcpm_all["PTPN22", ])
  ptpn22_z <- as.numeric(scale(ptpn22))
  response_R <- as.integer(state_meta$response == "anti-PD1_R")
  design <- model.matrix(~ response_R + ptpn22_z)

  keep_genes <- filterByExpr(dge_all, group = state_meta$response, min.count = 5)
  dge <- dge_all[keep_genes, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  voomed <- voom(dge, design, plot = FALSE)
  fit <- eBayes(lmFit(voomed, design), robust = TRUE)
  tt <- topTable(fit, coef = "ptpn22_z", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$state <- state_name
  tt$n_patient_states <- nrow(state_meta)
  gi <- gi + 1L
  gene_results[[gi]] <- tt[, c("state", "n_patient_states", "gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

  for (program_name in names(program_sets)) {
    genes <- intersect(program_sets[[program_name]], rownames(logcpm_all))
    genes <- genes[apply(logcpm_all[genes, , drop = FALSE], 1L, sd) > 0]
    if (length(genes) < 3L) next
    z <- t(scale(t(logcpm_all[genes, , drop = FALSE])))
    score <- colMeans(z, na.rm = TRUE)
    score_fit <- lm(score ~ response_R + ptpn22_z)
    sm <- summary(score_fit)$coefficients
    ci <- confint(score_fit, "ptpn22_z")
    pi <- pi + 1L
    program_results[[pi]] <- data.frame(
      state = state_name,
      program = program_name,
      n_patient_states = length(score),
      n_genes = length(genes),
      genes_used = paste(genes, collapse = ";"),
      beta_per_1SD_PTPN22 = sm["ptpn22_z", "Estimate"],
      ci_lower = ci[1L],
      ci_upper = ci[2L],
      p_value = sm["ptpn22_z", "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    si <- si + 1L
    score_rows[[si]] <- data.frame(
      pseudobulk_ID = state_meta$pseudobulk_ID,
      patient_ID = state_meta$patient_ID,
      response = state_meta$response,
      state = state_name,
      program = program_name,
      PTPN22_logCPM = ptpn22,
      score = score,
      stringsAsFactors = FALSE
    )
  }

  ranks <- tt$t
  names(ranks) <- tt$gene
  ranks <- sort(ranks[is.finite(ranks) & names(ranks) != "PTPN22"], decreasing = TRUE)
  fg <- suppressWarnings(fgseaMultilevel(
    pathways = msig_pathways,
    stats = ranks,
    minSize = 10,
    maxSize = 500,
    eps = 1e-10,
    nproc = 1
  ))
  fg <- as.data.frame(fg)
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1L))
  fg$state <- state_name
  if (!"nMoreExtreme" %in% names(fg)) fg$nMoreExtreme <- NA_integer_
  fi <- fi + 1L
  fgsea_results[[fi]] <- fg[, c("state", "pathway", "pval", "padj", "ES", "NES", "nMoreExtreme", "size", "leadingEdge")]
}

gene_results <- do.call(rbind, gene_results)
program_results <- do.call(rbind, program_results)
score_rows <- do.call(rbind, score_rows)
fgsea_results <- do.call(rbind, fgsea_results)

write.table(gene_results, file.path(table_dir, "PTPN22_FUNCTIONAL_GENE_ASSOCIATIONS.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(program_results, file.path(table_dir, "PTPN22_PATIENT_AWARE_PROGRAM_ASSOCIATIONS.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(score_rows, file.path(table_dir, "PTPN22_PATIENT_PROGRAM_SCORES.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fgsea_results, file.path(table_dir, "PTPN22_FUNCTIONAL_PATHWAY_ENRICHMENT.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE)

log_file <- file.path(log_dir, "07_functional_program_GSE206325.log")
sink(log_file, split = TRUE)
cat("PROGRAM_ASSOCIATIONS\n")
print(program_results[order(program_results$state, program_results$p_value), ])
cat("\nTOP_FGSEA_PER_STATE\n")
top_fg <- do.call(rbind, lapply(split(fgsea_results, fgsea_results$state), function(x) {
  head(x[order(x$padj, -abs(x$NES)), ], 15L)
}))
print(top_fg[, c("state", "pathway", "NES", "padj", "size")])
cat("\nSESSION_INFO\n")
print(sessionInfo())
sink()
