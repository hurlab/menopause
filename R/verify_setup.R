################################################################################
## R Setup Verification Script
## Checks if all required R packages are installed and accessible
################################################################################

cat("\n===========================================\n")
cat("R SETUP VERIFICATION\n")
cat("===========================================\n\n")

## Check if pacman is installed, if not install it
# Ensure a default CRAN mirror is set for non-interactive runs
.repos <- getOption("repos")
.cran <- unname(.repos["CRAN"])
if (is.null(.cran) || is.na(.cran) || .cran == "@CRAN@" || .cran == "") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
rm(.repos, .cran)

if (!requireNamespace("pacman", quietly = TRUE)) {
  cat("Installing pacman package...\n")
  install.packages("pacman")
}

## Load required packages
cat("Checking required packages...\n")

# For R 4.5, use Bioconductor 3.22
if (requireNamespace("BiocManager", quietly = TRUE)) {
  BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)
}

required_packages <- c(
  "data.table", "readr", "dplyr", "tidyr", "stringr",
  "tibble", "purrr", "glue", "magrittr",
  "ggplot2", "patchwork", "vsn",
  "Matrix", "matrixStats"
)

## Load packages using pacman
pacman::p_load(char = required_packages, character.only = TRUE)

## Load Bioconductor packages
cat("Loading Bioconductor packages...\n")
bioc_packages <- c("DESeq2", "AnnotationDbi", "org.Hs.eg.db")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s from Bioconductor...\n", pkg))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  library(pkg, character.only = TRUE)
  cat(sprintf("  [OK] %s\n", pkg))
}

cat("\nAll required packages loaded successfully!\n\n")

## Check data files
cat("Checking data files...\n")
data_dir <- "GTExDatav10"
required_files <- c(
  "gene_reads_v10_adipose_subcutaneous.gct.gz",
  "gene_reads_v10_adipose_visceral_omentum.gct.gz",
  "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt",
  "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt"
)

for (f in required_files) {
  path <- file.path(data_dir, f)
  if (file.exists(path)) {
    size <- file.size(path) / 1024 / 1024  # MB
    cat(sprintf("  [OK] %s (%.1f MB)\n", f, size))
  } else {
    cat(sprintf("  [MISSING] %s\n", f))
  }
}

## Check output directories
cat("\nChecking output directories...\n")
out_dir <- "GTEx_v10_AT_analysis_out"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "figs"), showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(out_dir, "rds"), showWarnings = FALSE)

cat(sprintf("  [OK] Output directory: %s\n", out_dir))

cat("\n===========================================\n")
cat("SETUP VERIFICATION COMPLETE\n")
cat("===========================================\n\n")

cat("You can now run the menopause signature discovery:\n")
cat("  Rscript R/run_all_menopause_approaches.R\n\n")

cat("Optional: run the original exploratory GTEx PCA/DE script (archived):\n")
cat("  Rscript R/archive/01_GTEx_analysis_v0.2.R\n\n")
cat("Or check for batch confounding:\n")
cat("  Rscript R/check_batch_confound.R\n\n")
