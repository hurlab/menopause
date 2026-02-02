################################################################################
## APPROACH 3: GSE86244 Dataset Analysis
##
## Strategy: Analyze GSE86244 dataset with ACTUAL menopause status
## Reference: Shan et al. - Adipose-derived stem cells from pre/post-menopausal women
##
## Note: This script provides a template for analyzing GSE86244
## The actual data needs to be downloaded from GEO
################################################################################
# Ensure a default CRAN mirror is set for non-interactive runs
.repos <- getOption("repos")
.cran <- unname(.repos["CRAN"])
if (is.null(.cran) || is.na(.cran) || .cran == "@CRAN@" || .cran == "") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
rm(.repos, .cran)


if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  data.table, readr, dplyr, tidyr, stringr, tibble, purrr, glue, magrittr,
  ggplot2, patchwork, DESeq2, vsn, Matrix, matrixStats, AnnotationDbi, org.Hs.eg.db,
  GEOquery, limma, hgu133plus.db
)

## Utility functions (ensure_dirs, etc.)
source("R/utils.R")

## -----------------------------------------------------------------------------
## Setup
## -----------------------------------------------------------------------------
out_dir  <- "GTEx_v10_AT_analysis_out"
fig_dir  <- file.path(out_dir, "figs", "approach3_gse86244")
tab_dir  <- file.path(out_dir, "tables", "approach3_gse86244")
ensure_dirs(c(out_dir, fig_dir, tab_dir))

cat("\n===========================================\n")
cat("APPROACH 3: GSE86244 ANALYSIS\n")
cat("===========================================\n\n")

cat("Dataset: GSE86244 (Shan et al.)\n")
cat("Tissue: Adipose-derived stem cells\n")
cat("Groups: Premenopausal (<45, n=12) vs Postmenopausal (>55, n=3)\n\n")

## -----------------------------------------------------------------------------
## Download GSE86244 data from GEO
## -----------------------------------------------------------------------------
cat("Checking for local GSE86244 data...\n")

local_data_dir <- "GSE86244_data"
if (dir.exists(local_data_dir)) {
  cat(sprintf("Found local data directory: %s\n", local_data_dir))
} else {
  cat("Local data not found. Would you like to download from GEO?\n")
  cat("To download manually:\n")
  cat("  1. Visit https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86244\n")
  cat("  2. Download series matrix and processed data\n")
  cat("  3. Save to: ", local_data_dir, "\n\n")

  ensure_dirs(local_data_dir)

  ## Attempt automated download
  cat("Attempting automated download via GEOquery...\n")
  tryCatch({
    gse <- getGEO("GSE86244", GSEMatrix = FALSE)
    cat("Successfully downloaded GSE86244 metadata\n")

    ## Save metadata
    saveRDS(gse, file.path(local_data_dir, "GSE86244_metadata.rds"))
  }, error = function(e) {
    cat("Automated download failed:\n")
    cat(sprintf("  Error: %s\n", e$message))
    cat("\nPlease download data manually from GEO.\n")

    ## Create placeholder for template code
    cat("\nCreating template analysis script for when data is available...\n")
  })
}

## -----------------------------------------------------------------------------
## Template analysis (when data is available)
## -----------------------------------------------------------------------------
template_script <- '
# Template for analyzing GSE86244 when data is downloaded
# Replace "GSE86244_data/" with actual data path

library(limma)
library(GEOquery)

# Load expression data
# Options:
# 1. Load processed data from GEO
gse <- getGEO("GSE86244", GSEMatrix = TRUE)
expr_data <- exprs(gse[[1]])
pheno_data <- pData(gse[[1]])

# 2. Load local data if downloaded
# expr_data <- read.table("GSE86244_data/expression_matrix.txt", header=TRUE, row.names=1)
# pheno_data <- read.table("GSE86244_data/metadata.txt", header=TRUE, sep="\\t")

# Filter for female samples only (all should be female in this dataset)
# This dataset only has females: 12 premenopausal (<45), 3 postmenopausal (>55)

# Create design matrix
group <- factor(c(rep("Pre", 12), rep("Post", 3)))
design <- model.matrix(~ group)

# Fit linear model
fit <- lmFit(expr_data, design)
fit <- eBayes(fit)

# Get results
results <- topTable(fit, coef=2, number=Inf, adjust="fdr")

# Filter significant genes
sig_genes <- results[results$adj.P.Val < 0.05, ]

# Save results
write.table(results, "GSE86244_DE_results.tsv", sep="\\t", quote=FALSE)
write.table(sig_genes, "GSE86244_significant_genes.tsv", sep="\\t", quote=FALSE)

# Create gene list for signature
sig_genes <- rownames(sig_genes)
writeLines(sig_genes, "menopause_signature_gse86244.txt")

cat(sprintf("Found %d significant genes\\n", nrow(sig_genes)))
'

## Save template script
writeLines(template_script, file.path(tab_dir, "GSE86244_analysis_template.R"))

cat("\n===========================================\n")
cat("APPROACH 3: SUMMARY\n")
cat("===========================================\n\n")

cat("Status: Template created\n")
cat("\nNext steps:\n")
cat("  1. Download GSE86244 data from GEO\n")
cat("  2. Place in: ", local_data_dir, "\n")
cat("  3. Run: Rscript ", file.path(tab_dir, "GSE86244_analysis_template.R"), "\n\n")

cat("Key advantages of GSE86244:\n")
cat("  - Actual menopause status (not age proxy)\n")
cat("  - Same tissue type (adipose-derived stem cells)\n")
cat("  - Published menopause signature candidates\n")
cat("  - Small but well-defined groups\n\n")

cat("Limitations:\n")
cat("  - Small sample size (n=15 total)\n")
cat("  - Different cell type (ASC vs whole adipose)\n")
cat("  - Requires external download\n")

cat("\nOutput files:\n")
cat(sprintf("  Template script: %s\n", file.path(tab_dir, "GSE86244_analysis_template.R")))
cat("===========================================\n\n")
