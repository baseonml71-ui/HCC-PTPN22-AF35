options(stringsAsFactors = FALSE, width = 240)
suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

workspace_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
project_dir <- file.path(workspace_dir, "HCC_ICI_project")
phase15a_dir <- file.path(project_dir, "03_processed_data", "phase15a_horizon_access")
phase15b_processed <- file.path(project_dir, "03_processed_data", "phase15b_public_replication")
metadata_dir <- file.path(project_dir, "01_metadata", "phase15b_public_replication")
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)

g302_rds <- file.path(phase15b_processed, "GSE302495_RESPONSE_BLIND_MATRIX_AND_SCORES.rds")
g287_rds <- file.path(phase15b_processed, "GSE287319_RESPONSE_BLIND_MATRIX_AND_SCORES.rds")
g302_source <- file.path(phase15a_dir, "HCC_PHASE15A_GSE302495_PUBLIC_ENDPOINT_CROSSWALK.tsv")
g287_source <- file.path(phase15a_dir, "HCC_PHASE15A_GSE287319_PUBLIC_ENDPOINT_CROSSWALK.tsv")
stopifnot(all(file.exists(c(g302_rds, g287_rds, g302_source, g287_source))))

stopifnot(
  digest(g302_rds, algo = "sha256", file = TRUE) == "b138295da8af3f7f8965dc68ff429aa10cf63e90f3e0ce286a1185e299322d1e",
  digest(g287_rds, algo = "sha256", file = TRUE) == "07994a1ada18f4aeb8e499299a2d088f1ab9a23ab448f67e7e70998723e51034",
  digest(g302_source, algo = "sha256", file = TRUE) == "ef4f21ffaaeba58e9f56e5d0c1611c7e7bb35ba92be31d1ca0d87b334ccbd2f8",
  digest(g287_source, algo = "sha256", file = TRUE) == "d1c185371068c28eb7842825506a65aa380614ee6038e21729884f5eb8a5d767"
)

g302_object <- readRDS(g302_rds)
g287_object <- readRDS(g287_rds)
stopifnot(identical(g302_object$response_joined, FALSE), identical(g287_object$response_joined, FALSE))

g302 <- fread(g302_source)[timepoint == "PRETREATMENT"]
stopifnot(nrow(g302) == 38L, !anyDuplicated(g302$patient_id))
g302[, response_group := fifelse(frozen_response_group == "Responder", "R", fifelse(frozen_response_group == "Non-responder", "NR", NA_character_))]
g302[, pathological_response_group := fifelse(major_pathological_response == "Yes", "MPR", fifelse(major_pathological_response == "No", "NON_MPR", "NOT_EVALUABLE"))]
stopifnot(sum(g302$response_group == "R") == 12L, sum(g302$response_group == "NR") == 26L)
stopifnot(sum(g302$pathological_response_group == "MPR") == 7L, sum(g302$pathological_response_group == "NON_MPR") == 15L, sum(g302$pathological_response_group == "NOT_EVALUABLE") == 16L)
stopifnot(setequal(g302$patient_id, g302_object$patient_ids))
g302_map <- g302[, .(
  patient_id, sample_accession, sample_title, timepoint,
  recist_category, tumor_size_change_percent, response_group,
  surgery, major_pathological_response, pathological_response_group,
  endpoint_source, mapping_status
)]
setorder(g302_map, patient_id)
g302_out <- file.path(metadata_dir, "GSE302495_PATIENT_RESPONSE_MAP.tsv")
fwrite(g302_map, g302_out, sep = "\t", quote = TRUE)

g287 <- fread(g287_source)[timepoint == "PRETREATMENT"]
stopifnot(nrow(g287) == 9L, !anyDuplicated(g287$public_patient_id))
g287[, endpoint_eligible := public_lesion_dynamics_group %chin% c("REDUCTION", "PROGRESSION")]
stopifnot(sum(g287$public_lesion_dynamics_group == "REDUCTION") == 4L, sum(g287$public_lesion_dynamics_group == "PROGRESSION") == 4L, sum(!g287$endpoint_eligible) == 1L)
stopifnot(setequal(g287$public_patient_id, g287_object$patient_ids))
g287_map <- g287[, .(
  patient_id = public_patient_id, geo_sample_accession, geo_sample_title, timepoint,
  endpoint_group = public_lesion_dynamics_group, endpoint_eligible,
  endpoint_definition, mapping_source, mapping_status,
  recist_response_available, individual_pfs_time_available, individual_os_time_available
)]
setorder(g287_map, patient_id)
g287_out <- file.path(metadata_dir, "GSE287319_PATIENT_ENDPOINT_MAP.tsv")
fwrite(g287_map, g287_out, sep = "\t", quote = TRUE)

cat("GSE302495_MAP_ROWS=", nrow(g302_map), "\n", sep = "")
cat("GSE302495_RECIST_R_NR=", sum(g302_map$response_group == "R"), "/", sum(g302_map$response_group == "NR"), "\n", sep = "")
cat("GSE302495_PATH_MPR_NONMPR_NE=", sum(g302_map$pathological_response_group == "MPR"), "/", sum(g302_map$pathological_response_group == "NON_MPR"), "/", sum(g302_map$pathological_response_group == "NOT_EVALUABLE"), "\n", sep = "")
cat("GSE302495_MAP_SHA256=", digest(g302_out, algo = "sha256", file = TRUE), "\n", sep = "")
cat("GSE287319_MAP_ROWS=", nrow(g287_map), "\n", sep = "")
cat("GSE287319_ELIGIBLE_REDUCTION_PROGRESSION=", sum(g287_map$endpoint_eligible), "/", sum(g287_map$endpoint_group == "REDUCTION"), "/", sum(g287_map$endpoint_group == "PROGRESSION"), "\n", sep = "")
cat("GSE287319_MAP_SHA256=", digest(g287_out, algo = "sha256", file = TRUE), "\n", sep = "")
