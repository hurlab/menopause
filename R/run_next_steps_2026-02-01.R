################################################################################
## Driver: Next Steps Plan Implementation (2026-02-01)
##
## Runs the “next steps” analyses in a fixed order and writes logs/results into
## dedicated dated subfolders for easier sharing with collaborators.
################################################################################

scripts <- c(
  "R/next_steps_2026-02-01_covariate_sensitivity.R",
  "R/next_steps_2026-02-01_composition_markers.R",
  "R/next_steps_2026-02-01_adrenergic_modules.R",
  "R/next_steps_2026-02-01_signature_refinement.R",
  "R/next_steps_2026-02-01_external_dataset_search_and_triage.R",
  "R/external_validation_GSE44000_postmenopausal_adipocytes_obese_vs_lean.R"
)

for (s in scripts) {
  cat("\n============================================================\n")
  cat("Running:", s, "\n")
  cat("============================================================\n\n")
  if (!file.exists(s)) {
    stop("Missing script: ", s)
  }
  status <- system2("Rscript", c("--vanilla", s))
  if (!identical(status, 0L)) {
    stop("Script failed (non-zero exit): ", s)
  }
}

cat("\nAll next-steps scripts completed.\n")
