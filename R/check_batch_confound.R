################################################################################
## Batch Confounding Diagnostic Check
## Checks if age groups (40-49 vs 50-59) are confounded with batch variables
################################################################################
# Ensure a default CRAN mirror is set for non-interactive runs
.repos <- getOption("repos")
.cran <- unname(.repos["CRAN"])
if (is.null(.cran) || is.na(.cran) || .cran == "@CRAN@" || .cran == "") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
rm(.repos, .cran)


if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(data.table, dplyr, tidyr, ggplot2, readr)

## Load metadata
meta <- readr::read_tsv("GTEx_v10_AT_analysis_out/tables/sample_metadata_full.tsv",
                       show_col_types = FALSE)

## Clean age bin labels
meta$age_bin_label <- gsub("\\s*–\\s*|\\s*-\\s*", "-", meta$AGE)

## Filter to relevant samples (40-49 and 50-59)
meta_filt <- meta %>%
  filter(age_bin_label %in% c("40-49", "50-59"),
         DEPOT %in% c("Subcutaneous", "Visceral"))

cat("\n===========================================\n")
cat("BATCH CONFOUNDING DIAGNOSTIC CHECK\n")
cat("===========================================\n\n")

cat("Total samples in analysis:", nrow(meta_filt), "\n\n")

## Check 1: SMCENTER vs age_bin_label
cat("1. SMCENTER vs Age Bin Cross-Tabulation:\n")
cat("   --------------------------------------\n")
tbl_center <- table(meta_filt$age_bin_label, meta_filt$SMCENTER)
print(tbl_center)

# Chi-square test
if (nrow(tbl_center) > 1 && ncol(tbl_center) > 1) {
  chi_center <- chisq.test(tbl_center)
  cat("\n   Chi-square p-value:", format.pval(chi_center$p.value), "\n")
  if (chi_center$p.value < 0.05) {
    cat("   *** WARNING: Significant association detected! ***\n")
  }
}
cat("\n")

## Check 2: SMNABTCH vs age_bin_label
cat("2. SMNABTCH (NA Batch) vs Age Bin Cross-Tabulation:\n")
cat("   -----------------------------------------------\n")
tbl_na <- table(meta_filt$age_bin_label, meta_filt$SMNABTCH)
print(tbl_na)

if (nrow(tbl_na) > 1 && ncol(tbl_na) > 1) {
  chi_na <- chisq.test(tbl_na)
  cat("\n   Chi-square p-value:", format.pval(chi_na$p.value), "\n")
  if (chi_na$p.value < 0.05) {
    cat("   *** WARNING: Significant association detected! ***\n")
  }
}
cat("\n")

## Check 3: SMGEBTCH vs age_bin_label
cat("3. SMGEBTCH (GE Batch) vs Age Bin Cross-Tabulation:\n")
cat("   -----------------------------------------------\n")
tbl_ge <- table(meta_filt$age_bin_label, meta_filt$SMGEBTCH)
print(tbl_ge)

if (nrow(tbl_ge) > 1 && ncol(tbl_ge) > 1) {
  chi_ge <- chisq.test(tbl_ge)
  cat("\n   Chi-square p-value:", format.pval(chi_ge$p.value), "\n")
  if (chi_ge$p.value < 0.05) {
    cat("   *** WARNING: Significant association detected! ***\n")
  }
}
cat("\n")

## Check 4: Separate by tissue type
cat("4. Analysis by Tissue (DEPOT):\n")
cat("   ===========================\n\n")

for (depot in c("Subcutaneous", "Visceral")) {
  sub <- meta_filt %>% filter(DEPOT == depot)

  cat(sprintf("   %s (n=%d):\n", depot, nrow(sub)))

  # SMCENTER
  tbl <- table(sub$age_bin_label, sub$SMCENTER)
  if (nrow(tbl) > 0 && ncol(tbl) > 0) {
    cat(sprintf("     SMCENTER distribution:\n"))
    print(tbl)
    if (nrow(tbl) > 1 && ncol(tbl) > 1) {
      chi <- chisq.test(tbl)
      cat(sprintf("     Chi-square p-value: %s\n", format.pval(chi$p.value)))
    }
  }

  # SMNABTCH (show top batches only if many)
  tbl_na <- table(sub$age_bin_label, sub$SMNABTCH)
  if (ncol(tbl_na) > 10) {
    # Show summary instead
    cat(sprintf("\n     SMNABTCH: %d unique batches (showing counts by age)\n",
                ncol(tbl_na)))
    cat("       40-49 samples:", sum(sub$age_bin_label == "40-49"), "\n")
    cat("       50-59 samples:", sum(sub$age_bin_label == "50-59"), "\n")
  } else {
    cat(sprintf("\n     SMNABTCH distribution:\n"))
    print(tbl_na)
    if (nrow(tbl_na) > 1 && ncol(tbl_na) > 1) {
      chi <- chisq.test(tbl_na)
      cat(sprintf("     Chi-square p-value: %s\n", format.pval(chi$p.value)))
    }
  }
  cat("\n")
}

## Check 5: Visualization data
cat("5. Creating visualization...\n")

# Create a summary for plotting
summary_data <- meta_filt %>%
  group_by(age_bin_label, SMCENTER) %>%
  summarise(n = n(), .groups = "drop")

p1 <- ggplot(summary_data, aes(x = SMCENTER, y = n, fill = age_bin_label)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Sample Count by SMCENTER and Age Bin",
       x = "SMCENTER", y = "Sample Count", fill = "Age Bin") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("GTEx_v10_AT_analysis_out/figs/batch_confound_SMCENTER.png",
       p1, width = 10, height = 6, dpi = 300)

# By depot
for (depot in c("Subcutaneous", "Visceral")) {
  sub <- meta_filt %>% filter(DEPOT == depot)
  summary_data <- sub %>%
    group_by(age_bin_label, SMCENTER) %>%
    summarise(n = n(), .groups = "drop")

  p <- ggplot(summary_data, aes(x = SMCENTER, y = n, fill = age_bin_label)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste("Sample Count by SMCENTER and Age Bin -", depot),
         x = "SMCENTER", y = "Sample Count", fill = "Age Bin") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(sprintf("GTEx_v10_AT_analysis_out/figs/batch_confound_SMCENTER_%s.png",
                gsub(" ", "_", depot)),
         p, width = 10, height = 6, dpi = 300)
}

cat("   Saved visualizations to GTEx_v10_AT_analysis_out/figs/\n\n")

## Summary
cat("===========================================\n")
cat("SUMMARY:\n")
cat("===========================================\n")
cat("If any chi-square p-value < 0.05, this indicates BATCH CONFOUNDING.\n")
cat("This means age groups are unevenly distributed across batches,\n")
cat("and batch correction is CRITICAL for valid DE analysis.\n")
cat("===========================================\n\n")
