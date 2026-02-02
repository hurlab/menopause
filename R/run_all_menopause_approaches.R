################################################################################
## Menopause Signature Discovery - Orchestration Script
##
## Runs all four approaches sequentially to discover menopause gene signatures
##
## Approaches:
## 1. Improved age proxy (30-45 vs 55-70, excluding transition zone)
## 2. Estrogen-responsive gene filter
## 3. GSE86244 dataset analysis (with actual menopause status)
## 4. Literature review and signature comparison
################################################################################
# Ensure a default CRAN mirror is set for non-interactive runs
.repos <- getOption("repos")
.cran <- unname(.repos["CRAN"])
if (is.null(.cran) || is.na(.cran) || .cran == "@CRAN@" || .cran == "") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
rm(.repos, .cran)


cat("\n===========================================\n")
cat("MENOPAUSE SIGNATURE DISCOVERY\n")
cat("Running All Four Approaches\n")
cat("===========================================\n\n")

## List of R scripts to run
scripts <- c(
  "R/menopause_signature_v1_improved_age_proxy.R",
  "R/menopause_signature_v2_estrogen_filter.R",
  "R/menopause_signature_v3_gse86244.R",
  "R/menopause_signature_v4_literature.R"
)

approach_names <- c(
  "Approach 1: Improved Age Proxy",
  "Approach 2: Estrogen-Responsive Filter",
  "Approach 3: GSE86244 Dataset",
  "Approach 4: Literature Review"
)

## Create output directories
out_dir <- "GTEx_v10_AT_analysis_out"
ensure_dirs <- function(paths) {
  for (pth in paths) {
    dir.create(pth, showWarnings = FALSE, recursive = TRUE)
  }
}

## Function to run a single approach
run_approach <- function(script_path, approach_name) {
  cat(sprintf("\n===========================================\n"))
  cat(sprintf("RUNNING: %s\n", approach_name))
  cat(sprintf("Script: %s\n", script_path))
  cat("===========================================\n\n")

  result <- tryCatch({
    output <- system2("Rscript", args = c(script_path), stdout = TRUE, stderr = TRUE)
    list(status = "SUCCESS", output = output)
  }, error = function(e) {
    list(status = "ERROR", error = e$message, output = "")
  })

  return(list(
    approach = approach_name,
    script = script_path,
    result = result
  ))
}

## Run all approaches sequentially
cat("Launching all approaches...\n\n")

results <- list()
for (i in seq_along(scripts)) {
  res <- run_approach(scripts[i], approach_names[i])
  results[[i]] <- res
}

## -----------------------------------------------------------------------------
## Collect and summarize results
## -----------------------------------------------------------------------------
cat("\n\n===========================================\n")
cat("SUMMARY OF ALL APPROACHES\n")
cat("===========================================\n\n")

for (i in seq_along(results)) {
  res <- results[[i]]
  cat(sprintf("\n%s\n", res$approach))
  cat(sprintf("Status: %s\n", res$result$status))
  if (res$result$status == "ERROR") {
    cat(sprintf("Error: %s\n", res$result$error))
  }
}

## -----------------------------------------------------------------------------
## Check for generated signature files
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("GENERATED SIGNATURE FILES\n")
cat("===========================================\n\n")

# Find all signature gene lists
sig_files <- list.files(
  file.path(out_dir, "tables"),
  # list.files() takes a regex pattern, not a glob.
  pattern = "signature.*\\.txt$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(sig_files) > 0) {
  for (f in sig_files) {
    n_lines <- length(readLines(f))
    cat(sprintf("%s (%d genes)\n", basename(f), n_lines))
  }
} else {
  cat("No signature files generated yet.\n")
}

## -----------------------------------------------------------------------------
## Create consolidated summary
## -----------------------------------------------------------------------------
summary_file <- file.path(out_dir, "tables", "menopause_discovery_summary.md")

summary_content <- paste0(
  "# Menopause Signature Discovery - Summary\n\n",
  "## Approaches Completed\n\n",
  paste0("- ", approach_names, "\n"),
  "\n## Output Locations\n\n",
  "```\n",
  out_dir, "/\n",
  "├── figs/\n",
  "│   ├── approach1_improved_age_proxy/\n",
  "│   ├── approach2_estrogen_filter/\n",
  "│   ├── approach3_gse86244/\n",
  "│   └── approach4_literature/\n",
  "└── tables/\n",
  "    ├── approach1_improved_age_proxy/\n",
  "    ├── approach2_estrogen_filter/\n",
  "    ├── approach3_gse86244/\n",
  "    └── approach4_literature/\n",
  "```\n\n",
  "## Next Steps\n\n",
  "1. Review signature gene lists from each approach\n",
  "2. Identify overlapping genes across approaches\n",
  "3. Create consensus menopause signature\n",
  "4. Validate against external datasets\n\n",
  "---\n",
  "Generated: ", Sys.time(), "\n"
)

writeLines(summary_content, summary_file)

cat(sprintf("\nSaved summary to: %s\n", summary_file))

cat("\n===========================================\n")
cat("ORCHESTRATION COMPLETE\n")
cat("===========================================\n\n")
