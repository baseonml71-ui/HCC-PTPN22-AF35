#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(rhdf5)
  library(png)
})

project <- "./HCC_ICI_project"
stage2b <- "./HCC_ICI_MAIN_FIGURE_FIGUREABILITY_STAGE2B_v1"
out <- "./HCC_ICI_BIOLOGICAL_RESOLUTION_BR1_v1"
data_dir <- file.path(out, "data")
proof_dir <- file.path(out, "proofs", "figure5")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proof_dir, recursive = TRUE, showWarnings = FALSE)

stop_exactly_one <- function(x, label) {
  if (length(x) != 1L) stop(sprintf("Expected exactly one %s, found %d", label, length(x)))
  x
}

read_stage2b <- function(path, expected_n) {
  x <- as.data.table(readRDS(path))
  if (nrow(x) != expected_n) stop(sprintf("Unexpected rows in %s: %d", basename(path), nrow(x)))
  x
}

wu_rds <- file.path(stage2b, "Figure5_STAGE2B_WU_CONDITIONAL_RESIDUAL_SPOTS.rds")
gse_rds <- file.path(stage2b, "Figure5_STAGE2B_GSE238264_CONDITIONAL_RESIDUAL_SPOTS.rds")
wu <- read_stage2b(wu_rds, 11637L)
gse <- read_stage2b(gse_rds, 17292L)

# -----------------------------------------------------------------------------
# GSE238264: deterministic join of source numeric clusters to frozen Stage 2B
# spot values. No biological names are introduced.
# -----------------------------------------------------------------------------
gse_cluster_parts <- vector("list", length(unique(gse$sample_id)))
for (i in seq_along(sort(unique(gse$sample_id)))) {
  sample_id <- sort(unique(gse$sample_id))[i]
  sample_dir <- file.path(project, "02_raw_data", "GSE238264", "extracted", sample_id)
  cluster_path <- stop_exactly_one(
    list.files(sample_dir, pattern = "^data_SME_.*_identity\\.csv$", full.names = TRUE),
    sprintf("source numeric cluster file for %s", sample_id)
  )
  part <- fread(cluster_path)
  setnames(part, names(part)[1], "barcode")
  if (!"X_pca_kmeans" %in% names(part)) stop(sprintf("Missing X_pca_kmeans in %s", cluster_path))
  part <- part[, .(
    sample_id = sample_id,
    barcode = as.character(barcode),
    source_numeric_cluster = as.integer(X_pca_kmeans),
    source_numeric_cluster_file = normalizePath(cluster_path, winslash = "/", mustWork = TRUE)
  )]
  if (anyDuplicated(part$barcode)) stop(sprintf("Duplicated barcodes in %s", cluster_path))
  gse_cluster_parts[[i]] <- part
}
gse_clusters <- rbindlist(gse_cluster_parts)
setkey(gse_clusters, sample_id, barcode)
setkey(gse, sample_id, barcode)
gse_plot <- gse_clusters[gse]
if (anyNA(gse_plot$source_numeric_cluster)) stop("Not all frozen GSE238264 spots mapped to source numeric clusters")
if (nrow(gse_plot) != 17292L) stop("GSE238264 join changed row count")

gse_plot[, `:=`(
  source_numeric_cluster_label = paste0("source cluster ", source_numeric_cluster),
  coordinate_status = "released tissue_positions_list.csv pixel coordinates; no re-embedding",
  biological_unit = "patient; spots nested",
  semantic_boundary = "source numeric cluster only; no ecology/tumor/stroma/cell-type interpretation"
)]

gse_plot[, `:=`(
  histology_path = NA_character_,
  tissue_hires_scalef = NA_real_,
  x_hires = NA_real_,
  y_hires = NA_real_
)]
for (sid in unique(gse_plot$sample_id)) {
  sample_dir <- file.path(project, "02_raw_data", "GSE238264", "extracted", sid)
  scale_path <- file.path(sample_dir, "spatial", "scalefactors_json.json")
  hist_path <- file.path(sample_dir, "spatial", "tissue_hires_image.png")
  scale <- as.numeric(fromJSON(scale_path)$tissue_hires_scalef)
  gse_plot[sample_id == sid, `:=`(
    histology_path = normalizePath(hist_path, winslash = "/", mustWork = TRUE),
    tissue_hires_scalef = scale,
    x_hires = x * scale,
    y_hires = y * scale
  )]
}
setorder(gse_plot, sample_id, barcode)
fwrite(gse_plot, file.path(data_dir, "GSE238264_REPLICATION_MAP_PLOT_READY.tsv"), sep = "\t", na = "NA")

# -----------------------------------------------------------------------------
# Wu: exact Stage 2B residuals joined to the existing adjacent-section ST-CODEX
# mapping. Region-spot rows are retained and are not treated as replicates.
# -----------------------------------------------------------------------------
spot_map_path <- file.path(project, "05_results", "tables", "spatial_phase8", "PTPN22_CODEX_SPOT_NEIGHBORHOOD_FEATURES.tsv")
region_audit_path <- file.path(project, "05_results", "tables", "spatial_phase8", "WU_SPATIAL_COREGISTRATION_REGION_AUDIT.tsv")
spot_map <- fread(spot_map_path)
spot_map[, eligible_primary := as.logical(eligible_primary)]
spot_map_pre <- spot_map[state == "Pre" & eligible_primary == TRUE]
if (nrow(spot_map_pre) != 5205L) stop(sprintf("Expected 5205 eligible pretreatment region-spot rows, got %d", nrow(spot_map_pre)))
setkey(wu, spot_id)
setkey(spot_map_pre, spot_id)
wu_context <- wu[spot_map_pre, allow.cartesian = TRUE]
if (nrow(wu_context) != 5205L) stop(sprintf("Wu context join changed row count: %d", nrow(wu_context)))
if (anyNA(wu_context$conditional_score)) stop("Wu mapping rows missing frozen Stage 2B residuals")
if (any(wu_context$patient_id != wu_context$i.patient_id)) stop("Patient mismatch across frozen ST and ST-CODEX sources")

wu_context[, `:=`(
  mapped_patient_id = i.patient_id,
  mapped_response = response,
  mapped_state = state,
  adjacent_section_status = "released adjacent-section ST-CODEX mapping; not identical-cell colocalization",
  biological_unit = "patient; region-spot rows nested",
  inference_status = "descriptive only; no response test; no new neighborhood inference"
)]
wu_context[, "i.patient_id" := NULL]
setorder(wu_context, patient_id, region_id, spot_id)
fwrite(wu_context, file.path(data_dir, "WU_COREGISTERED_TISSUE_CONTEXT_PLOT_READY.tsv"), sep = "\t", na = "NA")

# Freeze and record the response-blind technical representative-region rule.
region_audit <- fread(region_audit_path)
selection_pool <- region_audit[state == "Pre"]
setorder(selection_pool, -n_spots_primary_eligible_ge5, region_id)
selected <- selection_pool[1]
selected_region <- selected$region_id
selected_context <- wu_context[region_id == selected_region]
selected_sample <- stop_exactly_one(unique(selected_context$sample_id), "Visium sample for selected Wu region")
selection_out <- selected[, .(
  selection_rule = "pretreatment regions; decreasing n_spots_primary_eligible_ge5; lexical region_id tie-break",
  selection_variables = "state + mapping coverage + region_id only; response and AF35 not used",
  region_id,
  region_label,
  patient_id,
  state,
  response_annotation_only = response,
  n_cells,
  n_spots_primary_eligible_ge5,
  visium_sample_id = selected_sample,
  adjacent_section_status = "released adjacent-section ST-CODEX mapping; not identical-cell colocalization"
)]
fwrite(selection_out, file.path(data_dir, "WU_MODULE_D_PROOF_REGION_SELECTION.tsv"), sep = "\t", na = "NA")

# Resolve and export the exact source H&E image and its scale factor.
h5_path <- file.path(project, "03_processed_data", "spatial_phase8", "Visium-ST", "visium_all.h5ad")
h5_listing <- h5ls(h5_path, recursive = TRUE)
sample_root <- paste0("/uns/spatial/", selected_sample, "/")
hires_rows <- h5_listing[startsWith(h5_listing$group, sample_root) & h5_listing$name == "hires", ]
if (nrow(hires_rows) != 1L) stop(sprintf("Expected one H&E hires dataset for %s, got %d", selected_sample, nrow(hires_rows)))
hires_h5_path <- paste0(hires_rows$group, "/", hires_rows$name)
image_array <- h5read(h5_path, hires_h5_path)
if (length(dim(image_array)) != 3L) stop("Unexpected H&E array dimensionality")
if (dim(image_array)[1] %in% c(3L, 4L)) {
  image_array <- aperm(image_array, c(3, 2, 1))
} else if (!dim(image_array)[3] %in% c(3L, 4L)) {
  stop(sprintf("Unexpected H&E channel dimensions: %s", paste(dim(image_array), collapse = "x")))
}
he_source_path <- file.path(proof_dir, "WU_SELECTED_REGION_SOURCE_HE.png")
writePNG(image_array, he_source_path)

library_rows <- unique(sub("/images$", "", hires_rows$group))
if (length(library_rows) != 1L) stop("Unable to resolve unique Wu spatial library path")
scale_h5_path <- paste0(library_rows, "/scalefactors/tissue_hires_scalef")
hires_scale <- as.numeric(h5read(h5_path, scale_h5_path))
if (length(hires_scale) != 1L || !is.finite(hires_scale)) stop("Invalid Wu H&E hires scale factor")

selected_st <- copy(wu[sample_id == selected_sample])
selected_st[, `:=`(
  x_hires = x * hires_scale,
  y_hires = y * hires_scale,
  tissue_hires_scalef = hires_scale,
  histology_path = normalizePath(he_source_path, winslash = "/", mustWork = TRUE),
  selected_region_id = selected_region,
  coordinate_status = "Visium full-resolution coordinates scaled to released hires H&E",
  adjacent_section_status = "CODEX context is from a mapped adjacent section; no same-cell claim"
)]
selected_ann <- unique(selected_context[, .(spot_id, HE_Annotation, n_mapped_cells)])
if (anyDuplicated(selected_ann$spot_id)) stop("Selected region has duplicate annotations for a spot")
setkey(selected_st, spot_id)
setkey(selected_ann, spot_id)
selected_st <- selected_ann[selected_st]
setorder(selected_st, spot_id)
fwrite(selected_st, file.path(data_dir, "WU_SELECTED_SAMPLE_ST_LAYERS_PLOT_READY.tsv"), sep = "\t", na = "NA")

# Existing bounded ST x CODEX matrices are copied into an explicit plot-ready
# semantic wrapper. No coefficients, intervals, or tests are recomputed.
patient_coeff_path <- file.path(project, "05_results", "tables", "codex_activation_feedback_panel", "CODEX_AF_PANEL_PATIENT_COEFFICIENTS.tsv")
coeff <- fread(patient_coeff_path)
coeff[, `:=`(
  record_status = "existing frozen patient coefficient",
  biological_unit = "patient",
  adjacent_section_status = "released adjacent-section ST-CODEX tissue-context association",
  semantic_boundary = "contextual coefficient; not identical-cell colocalization; no BR1 response test"
)]
fwrite(coeff, file.path(data_dir, "WU_ST_CODEX_CONTEXT_MATRIX_PLOT_READY.tsv"), sep = "\t", na = "NA")

crossmodal_path <- file.path(project, "05_results", "tables", "spatial_specificity_rescue", "PTPN22_ST_CODEX_CROSSMODAL_CONCORDANCE.tsv")
crossmodal <- fread(crossmodal_path)
crossmodal[, `:=`(
  record_status = "existing frozen cross-modal summary",
  adjacent_section_status = "released adjacent-section ST-CODEX tissue-context association",
  semantic_boundary = "bounded four-feature context matrix; no BR1 recomputation"
)]
fwrite(crossmodal, file.path(data_dir, "WU_ST_CODEX_EXISTING_CONCORDANCE_PLOT_READY.tsv"), sep = "\t", na = "NA")

cat(sprintf("MODULE_D_SPATIAL_TABLES_PASS gse_rows=%d wu_context_rows=%d selected_region=%s selected_sample=%s\n",
            nrow(gse_plot), nrow(wu_context), selected_region, selected_sample))
