################################################################################
## Wider Age Group Comparison
##
## Tests DE between wider age bins to maximize contrast:
## - 30-39 vs 60-69 (younger vs older)
## - 20-29 vs 70-79 (youngest vs oldest, if sample size permits)
##
## Includes SMCENTER batch correction in all models
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
cat("WIDER AGE GROUP COMPARISON\n")
cat("===========================================\n\n")

## -----------------------------------------------------------------------------
## Load and prepare data
## -----------------------------------------------------------------------------
cat("Loading data...\n")

sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  rename_with(~ gsub("^X", "", .), everything())

subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  rename_with(~ gsub("^X", "", .), everything())

sample_attr <- sample_attr %>%
  mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])

meta_full <- sample_attr %>%
  left_join(subject_pheno, by = "SUBJID")

gct_list_sc <- read_gct_v12(gct_sc)

common_samples <- intersect(colnames(gct_list_sc$counts), meta_full$SAMPID)

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
## Define age group comparisons to test
## -----------------------------------------------------------------------------
comparisons <- list(
  list(name = "30-39 vs 60-69", young = "30-39", old = "60-69"),
  list(name = "20-29 vs 70-79", young = "20-29", old = "70-79"),
  list(name = "30-39 vs 50-59", young = "30-39", old = "50-59"),
  list(name = "40-49 vs 60-69", young = "40-49", old = "60-69")
)

## -----------------------------------------------------------------------------
## Function to run DE analysis for a given comparison
## -----------------------------------------------------------------------------
run_de_analysis <- function(meta, counts, young_age, old_age, comparison_name) {
  cat("\n===========================================\n")
  cat(sprintf("COMPARISON: %s\n", comparison_name))
  cat("===========================================\n\n")

  ## Filter samples
  meta_de <- meta %>%
    filter(
      sex == "Female",
      age_bin_label %in% c(young_age, old_age)
    ) %>%
    mutate(
      age_group = factor(age_bin_label, levels = c(young_age, old_age)),
      SMCENTER = factor(SMCENTER)
    )

  ## Check sample sizes
  n_young <- sum(meta_de$age_group == young_age)
  n_old <- sum(meta_de$age_group == old_age)

  cat(sprintf("Sample sizes:\n"))
  cat(sprintf("  %s: %d samples\n", young_age, n_young))
  cat(sprintf("  %s: %d samples\n", old_age, n_old))
  cat(sprintf("  Total: %d samples\n", nrow(meta_de)))

  ## Check batch distribution
  cat("\nBatch distribution:\n")
  batch_table <- table(meta_de$age_group, meta_de$SMCENTER)
  print(batch_table)

  chi_test <- tryCatch(chisq.test(batch_table), error = function(e) list(p.value = NA))
  if (!is.na(chi_test$p.value)) {
    cat(sprintf("\nChi-square p-value: %g\n", chi_test$p.value))
  }

  ## Check minimum sample size
  if (n_young < 10 || n_old < 10) {
    cat("\n*** WARNING: Sample size too small (< 10 per group) ***\n")
    cat("Skipping this comparison.\n\n")
    return(NULL)
  }

  ## Prepare count matrix
  samples_de <- meta_de$SAMPID
  counts_de <- counts[, samples_de, drop = FALSE]

  cat(sprintf("\nCount matrix: %d genes x %d samples\n", nrow(counts_de), ncol(counts_de)))

  ## Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = counts_de,
    colData = meta_de,
    design = ~ SMCENTER + age_group
  )

  ## Pre-filtering
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  cat(sprintf("After filtering: %d genes\n", nrow(dds)))

  ## Run DESeq
  cat("\nRunning DESeq...\n")
  dds <- DESeq(dds)

  ## Get results
  res <- results(dds, contrast = c("age_group", old_age, young_age))

  ## Summary
  cat("\nResults Summary:\n")
  summary(res)

  ## Convert to data frame
  df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    mutate(
      gene_symbol = gct_list_sc$gene_symbol[match(gene_id, gct_list_sc$gene_id)],
      comparison = comparison_name,
      young_age = young_age,
      old_age = old_age,
      significant = padj < 0.05 & !is.na(padj)
    )

  ## Count significant genes
  n_sig <- sum(df$significant, na.rm = TRUE)
  n_up <- sum(df$significant & df$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(df$significant & df$log2FoldChange < 0, na.rm = TRUE)

  cat(sprintf("\nSignificant genes (padj < 0.05): %d\n", n_sig))
  cat(sprintf("  Up in %s: %d\n", old_age, n_up))
  cat(sprintf("  Down in %s: %d\n", old_age, n_down))

  return(list(
    df = df,
    n_sig = n_sig,
    n_up = n_up,
    n_down = n_down,
    n_young = n_young,
    n_old = n_old,
    chi_p = chi_test$p.value
  ))
}

## -----------------------------------------------------------------------------
## Run all comparisons
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("RUNNING ALL AGE GROUP COMPARISONS\n")
cat("===========================================\n\n")

## Get count matrix (collapse duplicates once)
counts_sc <- gct_list_sc$counts[, meta$SAMPID, drop = FALSE]
counts_sc <- collapse_dupe_rowsum(counts_sc)

## Run each comparison
results_list <- list()

for (comp in comparisons) {
  result <- run_de_analysis(
    meta = meta,
    counts = counts_sc,
    young_age = comp$young,
    old_age = comp$old,
    comparison_name = comp$name
  )

  if (!is.null(result)) {
    results_list[[comp$name]] <- result
  }
}

## -----------------------------------------------------------------------------
## Combine and summarize results
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("SUMMARY OF ALL COMPARISONS\n")
cat("===========================================\n\n")

summary_table <- bind_rows(lapply(names(results_list), function(comp_name) {
  res <- results_list[[comp_name]]
  tibble::tibble(
    comparison = comp_name,
    n_young = res$n_young,
    n_old = res$n_old,
    n_total = res$n_young + res$n_old,
    n_sig = res$n_sig,
    n_up = res$n_up,
    n_down = res$n_down,
    batch_p = res$chi_p
  )
}))

cat("Comparison Summary:\n")
print(summary_table)

## -----------------------------------------------------------------------------
## Save detailed results
## -----------------------------------------------------------------------------
cat("\nSaving results...\n")

## Combine all results
all_results <- bind_rows(lapply(names(results_list), function(comp_name) {
  results_list[[comp_name]]$df
}))

readr::write_tsv(
  all_results %>%
    arrange(padj) %>%
    dplyr::select(gene_id, gene_symbol, comparison, young_age, old_age,
           baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, significant),
  file.path(tab_dir, "DE_wider_age_groups_with_batch_correction.tsv")
)

cat(sprintf("Saved to: %s\n", file.path(tab_dir, "DE_wider_age_groups_with_batch_correction.tsv")))

## Save summary table
readr::write_tsv(
  summary_table,
  file.path(tab_dir, "DE_comparison_summary_wider_ages.tsv")
)

cat(sprintf("Saved summary to: %s\n", file.path(tab_dir, "DE_comparison_summary_wider_ages.tsv")))

## -----------------------------------------------------------------------------
## Get top genes across all comparisons
## -----------------------------------------------------------------------------
if (sum(summary_table$n_sig) > 0) {
  cat("\n===========================================\n")
  cat("TOP SIGNIFICANT GENES (across all comparisons)\n")
  cat("===========================================\n\n")

  top_genes <- all_results %>%
    filter(significant) %>%
    arrange(padj) %>%
    head(20) %>%
    dplyr::select(gene_symbol, comparison, log2FoldChange, padj)

  print(top_genes)

  ## Save top genes
  readr::write_tsv(
    top_genes,
    file.path(tab_dir, "top_DE_genes_wider_ages.tsv")
  )

  cat(sprintf("\nSaved top genes to: %s\n", file.path(tab_dir, "top_DE_genes_wider_ages.tsv")))
} else {
  cat("\n===========================================\n")
  cat("NO SIGNIFICANT GENES FOUND IN ANY COMPARISON\n")
  cat("===========================================\n\n")
}

## -----------------------------------------------------------------------------
## Create visualization
## -----------------------------------------------------------------------------
cat("\nCreating visualizations...\n")

## Bar plot of significant genes by comparison
p_bar <- ggplot(summary_table, aes(x = reorder(comparison, n_sig), y = n_sig, fill = comparison)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n_sig), vjust = -0.5) +
  coord_flip() +
  labs(
    title = "Differential Expression: Wider Age Group Comparisons",
    subtitle = "Females, Subcutaneous Adipose, with SMCENTER batch correction",
    x = "Age Group Comparison",
    y = "Number of Significant Genes (padj < 0.05)"
  ) +
  theme_bw() +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "DE_wider_age_groups_summary.png"),
  p_bar,
  width = 10,
  height = 6,
  dpi = 300
)

## Sample size comparison
p_samples <- ggplot(summary_table, aes(x = comparison, y = n_total, fill = comparison)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0("n=", n_total)), vjust = -0.5) +
  labs(
    title = "Sample Sizes for Age Group Comparisons",
    x = "Age Group Comparison",
    y = "Total Sample Size"
  ) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(fig_dir, "DE_wider_age_groups_sample_sizes.png"),
  p_samples,
  width = 10,
  height = 6,
  dpi = 300
)

## Combined panel
p_combined <- p_bar / p_samples
ggsave(
  file.path(fig_dir, "DE_wider_age_groups_panel.png"),
  p_combined,
  width = 10,
  height = 12,
  dpi = 300
)

cat(sprintf("Saved figures to: %s\n", fig_dir))

## -----------------------------------------------------------------------------
## Final summary
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("FINAL SUMMARY\n")
cat("===========================================\n\n")

cat(sprintf("Total comparisons run: %d\n", length(results_list)))
cat(sprintf("Comparisons with significant genes: %d\n", sum(summary_table$n_sig > 0)))

if (any(summary_table$n_sig > 0)) {
  best_comp <- summary_table %>%
    filter(n_sig > 0) %>%
    arrange(desc(n_sig)) %>%
    head(1)

  cat(sprintf("\nMost fruitful comparison: %s\n", best_comp$comparison))
  cat(sprintf("  Significant genes: %d\n", best_comp$n_sig))
  cat(sprintf("  Sample sizes: %d vs %d\n", best_comp$n_young, best_comp$n_old))
} else {
  cat("\nNo significant genes found in any comparison.\n")
  cat("This suggests either:\n")
  cat("  1. No strong age-related expression changes in subcutaneous adipose\n")
  cat("  2. Sample sizes are still insufficient\n")
  cat("  3. Need to consider other approaches (continuous age, different depots)\n")
}

cat("\nFiles saved:\n")
cat(sprintf("  Tables: %s\n", tab_dir))
cat(sprintf("  Figures: %s\n", fig_dir))
cat("===========================================\n\n")
