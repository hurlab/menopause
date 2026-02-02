################################################################################
## Batch Correction Impact Analysis
##
## Runs DE analysis TWO ways and compares results:
## 1. WITHOUT batch correction (naive)
## 2. WITH SMCENTER as covariate (corrected)
##
## Quantifies the impact of batch correction on DE results
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
  ggplot2, patchwork, DESeq2, vsn, Matrix, matrixStats, AnnotationDbi, org.Hs.eg.db
)

## Load utility functions
source("R/utils.R")

## -----------------------------------------------------------------------------
## Paths
## -----------------------------------------------------------------------------
in_dir   <- "GTExDatav10"
gct_sc   <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
sattr    <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph   <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

out_dir  <- "GTEx_v10_AT_analysis_out"
fig_dir  <- file.path(out_dir, "figs")
tab_dir  <- file.path(out_dir, "tables")
rds_dir  <- file.path(out_dir, "rds")
ensure_dirs(c(out_dir, fig_dir, tab_dir, rds_dir))

cat("\n===========================================\n")
cat("BATCH CORRECTION IMPACT ANALYSIS\n")
cat("===========================================\n\n")
cat("This script compares DE results:\n")
cat("  1. WITHOUT batch correction (naive model)\n")
cat("  2. WITH SMCENTER as covariate (corrected model)\n\n")

## -----------------------------------------------------------------------------
## Load and prepare data
## -----------------------------------------------------------------------------
cat("Loading data...\n")

## Load sample attributes (tab-separated)
sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  # Clean up column names
  rename_with(~ gsub("^X", "", .), everything())

## Load subject phenotypes (tab-separated, contains AGE and SEX)
subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  rename_with(~ gsub("^X", "", .), everything())

## Extract subject ID from SAMPID (format: GTEX-XXXXX-XXXX-SM-XXXXX)
## The subject ID in GTEx is just GTEX-XXXXX (first two parts)
sample_attr <- sample_attr %>%
  mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])

## Merge sample attributes with subject phenotypes
meta_full <- sample_attr %>%
  left_join(subject_pheno, by = "SUBJID")

## Load gene expression data
gct_list_sc <- read_gct_v12(gct_sc)

## Subset to common samples
common_samples <- intersect(colnames(gct_list_sc$counts), meta_full$SAMPID)

## Filter to adipose tissue samples
meta <- meta_full %>%
  filter(SAMPID %in% common_samples, SMTS == "Adipose Tissue") %>%
  mutate(
    DEPOT = ifelse(SMTSD == "Adipose - Subcutaneous", "Subcutaneous",
                   ifelse(SMTSD == "Adipose - Visceral (Omentum)", "Visceral", "Other")),
    sex = ifelse(SEX == 1, "Male", ifelse(SEX == 2, "Female", NA))
  ) %>%
  filter(DEPOT == "Subcutaneous") %>%
  mutate(age_bin_label = AGE) %>%
  bind_cols(parse_age_bin(.$AGE)) %>%
  derive_age_numeric()

cat(sprintf("Loaded %d subcutaneous samples\n", nrow(meta)))

## -----------------------------------------------------------------------------
## Filter for DE analysis: Females, age 40-49 vs 50-59
## -----------------------------------------------------------------------------
meta_de <- meta %>%
  filter(
    sex == "Female",
    age_bin_label %in% c("40-49", "50-59")
  ) %>%
  mutate(
    age_group = factor(age_bin_label, levels = c("40-49", "50-59")),
    SMCENTER = factor(SMCENTER)
  )

cat(sprintf("\nDE analysis cohort:\n"))
cat(sprintf("  Females 40-49: %d samples\n", sum(meta_de$age_group == "40-49")))
cat(sprintf("  Females 50-59: %d samples\n", sum(meta_de$age_group == "50-59")))
cat(sprintf("  B1 center: %d samples\n", sum(meta_de$SMCENTER == "B1")))
cat(sprintf("  C1 center: %d samples\n", sum(meta_de$SMCENTER == "C1")))

## Prepare count matrix
counts_de <- gct_list_sc$counts[, meta_de$SAMPID, drop = FALSE]
counts_de <- collapse_dupe_rowsum(counts_de)

cat(sprintf("\nCount matrix: %d genes x %d samples\n", nrow(counts_de), ncol(counts_de)))

## -----------------------------------------------------------------------------
## ANALYSIS 1: WITHOUT batch correction (Naive)
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("ANALYSIS 1: WITHOUT Batch Correction (Naive)\n")
cat("===========================================\n\n")

cat("Design formula: design = ~ age_group\n\n")

## Create DESeq2 dataset WITHOUT batch
dds_naive <- DESeqDataSetFromMatrix(
  countData = counts_de,
  colData = meta_de,
  design = ~ age_group
)

## Pre-filtering
keep <- rowSums(counts(dds_naive)) >= 10
dds_naive <- dds_naive[keep, ]
cat(sprintf("After filtering: %d genes\n", nrow(dds_naive)))

## Run DESeq
cat("Running DESeq...\n")
dds_naive <- DESeq(dds_naive)

## Get results
res_naive <- results(dds_naive, contrast = c("age_group", "50-59", "40-49"))
# lfcShrink may fail with some coefficient name formats, use results directly if needed
tryCatch({
  res_naive <- lfcShrink(dds_naive, coef = "age_group_50.59_vs_40.49", type = "apeglm")
}, error = function(e) {
  cat("Note: lfcShrink failed, using unshrunk results\n")
})

## Summary
cat("\nNaive Model Results Summary:\n")
summary(res_naive)

## Convert to data frame
df_naive <- as.data.frame(res_naive) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    gene_symbol = gct_list_sc$gene_symbol[match(gene_id, gct_list_sc$gene_id)],
    model = "naive",
    significant = padj < 0.05 & !is.na(padj)
  )

## Count significant genes
n_sig_naive <- sum(df_naive$significant, na.rm = TRUE)
n_up_naive <- sum(df_naive$significant & df_naive$log2FoldChange > 0, na.rm = TRUE)
n_down_naive <- sum(df_naive$significant & df_naive$log2FoldChange < 0, na.rm = TRUE)

cat(sprintf("\nSignificant genes (padj < 0.05): %d\n", n_sig_naive))
cat(sprintf("  Up in 50-59: %d\n", n_up_naive))
cat(sprintf("  Down in 50-59: %d\n", n_down_naive))

## -----------------------------------------------------------------------------
## ANALYSIS 2: WITH batch correction (SMCENTER covariate)
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("ANALYSIS 2: WITH Batch Correction (SMCENTER Covariate)\n")
cat("===========================================\n\n")

cat("Design formula: design = ~ SMCENTER + age_group\n\n")

## Create DESeq2 dataset WITH batch
dds_corrected <- DESeqDataSetFromMatrix(
  countData = counts_de,
  colData = meta_de,
  design = ~ SMCENTER + age_group
)

## Use same filtering
dds_corrected <- dds_corrected[keep, ]
cat(sprintf("After filtering: %d genes\n", nrow(dds_corrected)))

## Run DESeq
cat("Running DESeq...\n")
dds_corrected <- DESeq(dds_corrected)

## Get results
res_corrected <- results(dds_corrected, contrast = c("age_group", "50-59", "40-49"))
# lfcShrink may fail with some coefficient name formats, use results directly if needed
tryCatch({
  res_corrected <- lfcShrink(dds_corrected, coef = "age_group_50.59_vs_40.49", type = "apeglm")
}, error = function(e) {
  cat("Note: lfcShrink failed, using unshrunk results\n")
})

## Summary
cat("\nCorrected Model Results Summary:\n")
summary(res_corrected)

## Convert to data frame
df_corrected <- as.data.frame(res_corrected) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    gene_symbol = gct_list_sc$gene_symbol[match(gene_id, gct_list_sc$gene_id)],
    model = "corrected",
    significant = padj < 0.05 & !is.na(padj)
  )

## Count significant genes
n_sig_corrected <- sum(df_corrected$significant, na.rm = TRUE)
n_up_corrected <- sum(df_corrected$significant & df_corrected$log2FoldChange > 0, na.rm = TRUE)
n_down_corrected <- sum(df_corrected$significant & df_corrected$log2FoldChange < 0, na.rm = TRUE)

cat(sprintf("\nSignificant genes (padj < 0.05): %d\n", n_sig_corrected))
cat(sprintf("  Up in 50-59: %d\n", n_up_corrected))
cat(sprintf("  Down in 50-59: %d\n", n_down_corrected))

## -----------------------------------------------------------------------------
## COMPARISON ANALYSIS
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("COMPARISON: Batch Correction Impact\n")
cat("===========================================\n\n")

## Merge results
df_combined <- bind_rows(df_naive, df_corrected) %>%
  dplyr::select(gene_id, gene_symbol, model, log2FoldChange, lfcSE, pvalue, padj, significant) %>%
  pivot_wider(
    names_from = model,
    values_from = c(log2FoldChange, lfcSE, pvalue, padj, significant),
    names_glue = "{model}_{.value}"
  ) %>%
  as_tibble()

## Classify genes by their behavior
df_comparison <- df_combined %>%
  mutate(
    sig_naive = coalesce(naive_significant, FALSE),
    sig_corrected = coalesce(corrected_significant, FALSE),
    batch_category = case_when(
      sig_naive & !sig_corrected ~ "Lost (false positive)",
      !sig_naive & sig_corrected ~ "Gained (masked by batch)",
      sig_naive & sig_corrected ~ "Stable (true signal)",
      !sig_naive & !sig_corrected ~ "Not significant",
      TRUE ~ "Other"
    ),
    lfc_change = corrected_log2FoldChange - naive_log2FoldChange,
    pct_lfc_change = (lfc_change / abs(naive_log2FoldChange)) * 100
  )

## Summary statistics
cat("Significant Gene Comparison:\n")
cat(sprintf("  Naive model: %d genes\n", n_sig_naive))
cat(sprintf("  Corrected model: %d genes\n", n_sig_corrected))

# Overlap
both_sig <- sum(df_comparison$sig_naive & df_comparison$sig_corrected)
only_naive <- sum(df_comparison$sig_naive & !df_comparison$sig_corrected)
only_corrected <- sum(!df_comparison$sig_naive & df_comparison$sig_corrected)

cat(sprintf("\nOverlap analysis:\n"))
cat(sprintf("  Significant in BOTH models: %d genes\n", both_sig))
cat(sprintf("  Only in NAIVE (likely false positives): %d genes\n", only_naive))
cat(sprintf("  Only in CORRECTED (masked by batch): %d genes\n", only_corrected))

## Venn diagram-style summary
venn_summary <- df_comparison %>%
  group_by(batch_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

cat("\nGene category breakdown:\n")
print(venn_summary)

## -----------------------------------------------------------------------------
## VISUALIZATIONS
## -----------------------------------------------------------------------------
cat("\nCreating visualizations...\n")

## Figure 1: Side-by-side volcano plots
df_naive_plot <- df_naive %>%
  mutate(
    neglog10p = -log10(pvalue),
    label = ifelse(significant & abs(log2FoldChange) > 1, gene_symbol, "")
  )

df_corrected_plot <- df_corrected %>%
  mutate(
    neglog10p = -log10(pvalue),
    label = ifelse(significant & abs(log2FoldChange) > 1, gene_symbol, "")
  )

p_volc_naive <- ggplot(df_naive_plot, aes(x = log2FoldChange, y = neglog10p)) +
  geom_point(aes(color = significant), alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70")) +
  labs(
    title = "WITHOUT Batch Correction",
    x = "log2 Fold Change (50-59 vs 40-49)",
    y = "-log10(p-value)",
    color = "Significant"
  ) +
  theme_bw()

p_volc_corrected <- ggplot(df_corrected_plot, aes(x = log2FoldChange, y = neglog10p)) +
  geom_point(aes(color = significant), alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70")) +
  labs(
    title = "WITH Batch Correction",
    x = "log2 Fold Change (50-59 vs 40-49)",
    y = "-log10(p-value)",
    color = "Significant"
  ) +
  theme_bw()

p_volc_combined <- p_volc_naive | p_volc_corrected
ggsave(file.path(fig_dir, "volcano_comparison_batch_correction.png"), p_volc_combined, width = 12, height = 5, dpi = 300)

## Figure 2: MA plots side by side
## plotMA() produces base graphics, so save via a graphics device instead of ggsave().
grDevices::png(file.path(fig_dir, "MA_plot_naive.png"), width = 8, height = 6, units = "in", res = 150)
plotMA(res_naive, main = "MA Plot: WITHOUT Batch Correction", ylim = c(-5, 5))
grDevices::dev.off()

grDevices::png(file.path(fig_dir, "MA_plot_corrected.png"), width = 8, height = 6, units = "in", res = 150)
plotMA(res_corrected, main = "MA Plot: WITH Batch Correction", ylim = c(-5, 5))
grDevices::dev.off()

## Figure 3: Scatter plot comparing logFC between models
p_scatter <- ggplot(df_comparison, aes(x = naive_log2FoldChange, y = corrected_log2FoldChange)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Log2 Fold Change Comparison",
    subtitle = "Naive vs Batch-Corrected Model",
    x = "Log2 FC (Naive Model)",
    y = "Log2 FC (Corrected Model)"
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "logFC_scatter_comparison.png"), p_scatter, width = 7, height = 6, dpi = 300)

## Figure 4: Bar plot of gene categories
p_category <- ggplot(venn_summary, aes(x = reorder(batch_category, n), y = n, fill = batch_category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "Impact of Batch Correction on DE Results",
    x = "Category",
    y = "Number of Genes",
    fill = "Category"
  ) +
  theme_bw() +
  scale_fill_brewer(type = "qual", palette = "Set2")

ggsave(file.path(fig_dir, "batch_correction_category_barplot.png"), p_category, width = 8, height = 5, dpi = 300)

## Figure 5: Summary panel
summary_panel <- (p_volc_combined / p_scatter) | p_category
ggsave(file.path(fig_dir, "batch_correction_summary_panel.png"), summary_panel, width = 14, height = 12, dpi = 300)

## -----------------------------------------------------------------------------
## SAVE RESULTS
## -----------------------------------------------------------------------------
cat("\nSaving results...\n")

## Full comparison table
readr::write_tsv(
  df_comparison %>%
    arrange(naive_padj) %>%
    dplyr::select(gene_id, gene_symbol, batch_category,
           naive_log2FoldChange, naive_padj, naive_significant,
           corrected_log2FoldChange, corrected_padj, corrected_significant,
           lfc_change, pct_lfc_change),
  file.path(tab_dir, "batch_correction_comparison_full.tsv")
)

## Top genes lost (false positives in naive)
lost_genes <- df_comparison %>%
  filter(batch_category == "Lost (false positive)") %>%
  arrange(naive_padj) %>%
  dplyr::select(gene_symbol, naive_log2FoldChange, naive_padj, corrected_padj, lfc_change)

readr::write_tsv(lost_genes, file.path(tab_dir, "genes_lost_after_correction.tsv"))

## Top genes gained (masked by batch)
gained_genes <- df_comparison %>%
  filter(batch_category == "Gained (masked by batch)") %>%
  arrange(corrected_padj) %>%
  dplyr::select(gene_symbol, corrected_log2FoldChange, naive_padj, corrected_padj, lfc_change)

readr::write_tsv(gained_genes, file.path(tab_dir, "genes_gained_after_correction.tsv"))

## Stable significant genes
stable_genes <- df_comparison %>%
  filter(batch_category == "Stable (true signal)") %>%
  arrange(corrected_padj) %>%
  dplyr::select(gene_symbol, naive_log2FoldChange, corrected_log2FoldChange, naive_padj, corrected_padj)

readr::write_tsv(stable_genes, file.path(tab_dir, "genes_stable_both_models.tsv"))

## -----------------------------------------------------------------------------
## FINAL SUMMARY
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("FINAL SUMMARY\n")
cat("===========================================\n\n")

cat("Results saved to:\n")
cat(sprintf("  Figures: %s\n", fig_dir))
cat(sprintf("  Tables: %s\n\n", tab_dir))

cat("Key Outputs:\n")
cat("  1. batch_correction_comparison_full.tsv - All genes compared\n")
cat(sprintf("  2. genes_lost_after_correction.tsv - Likely false positives (%d genes)\n", only_naive))
cat(sprintf("  3. genes_gained_after_correction.tsv - Batch-masked signals (%d genes)\n", only_corrected))
cat(sprintf("  4. genes_stable_both_models.tsv - Robust DE genes (%d genes)\n", both_sig))

cat("\nGenerated Figures:\n")
cat("  1. volcano_comparison_batch_correction.png - Side-by-side volcano plots\n")
cat("  2. logFC_scatter_comparison.png - Correlation between models\n")
cat("  3. batch_correction_category_barplot.png - Gene category breakdown\n")
cat("  4. batch_correction_summary_panel.png - Combined summary\n\n")

cat("INTERPRETATION GUIDE:\n")
cat(sprintf("  - %d genes were significant WITHOUT batch correction\n", n_sig_naive))
cat(sprintf("  - %d genes were significant WITH batch correction\n", n_sig_corrected))
cat(sprintf("  - %d genes (%.1f%%) were likely false positives\n",
            only_naive, 100*only_naive/max(1, n_sig_naive)))
cat(sprintf("  - %d genes were masked by batch effects\n", only_corrected))
cat(sprintf("  - %d genes are robust (significant in both)\n\n", both_sig))

cat("RECOMMENDATION:\n")
if (only_naive > n_sig_naive * 0.2) {
  cat("  *** BATCH CORRECTION IS CRITICAL ***\n")
  cat("  More than 20%% of 'significant' genes are false positives!\n")
  cat("  Use the corrected model for all biological interpretations.\n\n")
} else if (only_naive > 50) {
  cat("  ** Batch correction has MODERATE impact **\n")
  cat("  Use corrected model but verify results manually.\n\n")
} else {
  cat("  * Batch correction has MINOR impact *\n")
  cat("  Results are relatively robust to batch effects.\n\n")
}

cat("===========================================\n\n")
