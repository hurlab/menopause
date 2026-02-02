################################################################################
## Batch Correction Strategy Comparison for GTEx Analysis
##
## Two Approaches:
## 1. Include batch as covariate in DESeq2 design (RECOMMENDED for DE)
## 2. Apply ComBat correction before analysis (USEFUL for visualization)
################################################################################
# Ensure a default CRAN mirror is set for non-interactive runs
.repos <- getOption("repos")
.cran <- unname(.repos["CRAN"])
if (is.null(.cran) || is.na(.cran) || .cran == "@CRAN@" || .cran == "") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
rm(.repos, .cran)


if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(data.table, dplyr, tidyr, ggplot2, readr, DESeq2, vsn, sva)

cat("\n===========================================\n")
cat("BATCH CORRECTION STRATEGY GUIDE\n")
cat("===========================================\n\n")

cat("APPROACH 1: Include SMCENTER in DESeq2 Design\n")
cat("------------------------------------------------\n")
cat("Pros:\n")
cat("  + Statistically rigorous - DESeq2 handles batch in GLM\n")
cat("  + Preserves count distribution for testing\n")
cat("  + No risk of over-correcting biological signal\n")
cat("  + Can extract batch effect size from results\n")
cat("\nCons:\n")
cat("  - PCA/clustering still shows batch effects\n")
cat("  - More complex design formula\n")
cat("\nDesign formula:\n")
cat("  design = ~ SMCENTER + sex + age_mid + DEPOT + age_group\n\n")

cat("\nAPPROACH 2: ComBat Correction on VST Data\n")
cat("------------------------------------------\n")
cat("Pros:\n")
cat("  + Removes batch effects for visualization\n")
cat("  + Cleaner PCA plots after correction\n")
cat("  + Can use corrected data for clustering\n")
cat("\nCons:\n")
cat("  - Risk of removing biological signal if batch is confounded\n")
cat("  - Cannot use ComBat-corrected data for DESeq2 (breaks count model)\n")
cat("  - Must re-run DESeq2 on original counts\n")
cat("\nRecommended workflow:\n")
cat("  1. VST transform counts\n")
cat("  2. Apply ComBat to VST data (batch = SMCENTER)\n")
cat("  3. Use ComBat-corrected VST for PCA/visualization\n")
cat("  4. Use ORIGINAL counts for DESeq2 with batch in design\n\n")

cat("\nRECOMMENDATION:\n")
cat("---------------\n")
cat("For your GTEx analysis with age/sex/batch confounding:\n\n")

cat("PRIMARY APPROACH (DE analysis):\n")
cat("  Use DESeq2 with SMCENTER as covariate:\n")
cat("  dds <- DESeqDataSetFromMatrix(countData, colData, \n")
cat("                                design = ~ SMCENTER + sex + age_mid + DEPOT)\n")
cat("  dds <- DESeq(dds)\n")
cat("  results(dds, contrast=c('age_group', '50-59', '40-49'))\n\n")

cat("SECONDARY APPROACH (Visualization):\n")
cat("  Apply ComBat to VST data for PCA plots:\n")
cat("  vsd <- vst(dds, blind=FALSE)\n")
cat("  assay(vsd) <- ComBat_seq(counts(dds), batch=colData$SMCENTER)\n\n")

cat("BEST PRACTICE:\n")
cat("  1. Run BOTH approaches\n")
cat("  2. Compare DE results with and without batch correction\n")
cat("  3. If results differ substantially, batch correction is critical\n")
cat("  4. Report both corrected and uncorrected results\n\n")

cat("\n===========================================\n\n")
