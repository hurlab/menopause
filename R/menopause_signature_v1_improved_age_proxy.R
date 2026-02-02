################################################################################
## APPROACH 1: Improved Age Proxy for Menopause Signature
##
## Strategy: Use age 30-45 vs 55-70, excluding the 45-55 transition zone
## to improve classification accuracy for pre/post-menopausal status
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
## Setup
## -----------------------------------------------------------------------------
in_dir   <- "GTExDatav10"
gct_sc   <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
sattr    <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph   <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

out_dir  <- "GTEx_v10_AT_analysis_out"
fig_dir  <- file.path(out_dir, "figs", "approach1_improved_age_proxy")
tab_dir  <- file.path(out_dir, "tables", "approach1_improved_age_proxy")
ensure_dirs(c(out_dir, fig_dir, tab_dir))

cat("\n===========================================\n")
cat("APPROACH 1: IMPROVED AGE PROXY FOR MENOPAUSE\n")
cat("===========================================\n\n")
cat("Strategy:\n")
cat("  - Pre-menopausal proxy: age 30-45\n")
cat("  - Post-menopausal proxy: age 55-70\n")
cat("  - EXCLUDE transition zone (45-55)\n\n")

## -----------------------------------------------------------------------------
## Load data
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
## Define menopause proxy categories
## -----------------------------------------------------------------------------
meta <- meta %>%
  mutate(
    menopause_proxy = case_when(
      sex != "Female" ~ NA_character_,
      age_mid < 45 ~ "Pre",
      age_mid >= 45 & age_mid <= 55 ~ "Transition",
      age_mid > 55 ~ "Post",
      TRUE ~ NA_character_
    ),
    menopause_proxy = factor(menopause_proxy, levels = c("Pre", "Post"))
  )

## Show sample distribution
cat("\nMenopause Proxy Distribution:\n")
proxy_table <- table(meta$menopause_proxy, useNA = "ifany")
print(proxy_table)

## Filter to Pre and Post only (exclude transition zone)
meta_de <- meta %>%
  filter(!is.na(menopause_proxy), menopause_proxy %in% c("Pre", "Post"))

cat(sprintf("\nAnalysis cohort (excluding transition zone):\n"))
cat(sprintf("  Pre-menopausal (30-45): %d samples\n", sum(meta_de$menopause_proxy == "Pre")))
cat(sprintf("  Post-menopausal (55-70): %d samples\n", sum(meta_de$menopause_proxy == "Post")))

## Check batch distribution
cat("\nBatch distribution:\n")
batch_table <- table(meta_de$menopause_proxy, meta_de$SMCENTER)
print(batch_table)

chi_test <- chisq.test(batch_table)
cat(sprintf("\nChi-square p-value: %g\n", chi_test$p.value))

## -----------------------------------------------------------------------------
## Prepare count matrix
## -----------------------------------------------------------------------------
counts_sc <- gct_list_sc$counts[, meta_de$SAMPID, drop = FALSE]
counts_sc <- collapse_dupe_rowsum(counts_sc)

cat(sprintf("\nCount matrix: %d genes x %d samples\n", nrow(counts_sc), ncol(counts_sc)))

## -----------------------------------------------------------------------------
## DESeq2 analysis with batch correction
## -----------------------------------------------------------------------------
cat("\nRunning DESeq2 analysis...\n")

meta_de$SMCENTER <- factor(meta_de$SMCENTER)

dds <- DESeqDataSetFromMatrix(
  countData = counts_sc,
  colData = meta_de,
  design = ~ SMCENTER + menopause_proxy
)

## Pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(sprintf("After filtering: %d genes\n", nrow(dds)))

## Run DESeq
dds <- DESeq(dds)

## Get results (Post vs Pre)
res <- results(dds, contrast = c("menopause_proxy", "Post", "Pre"))

## Summary
cat("\nResults Summary:\n")
summary(res)

## Convert to data frame
df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    gene_symbol = gct_list_sc$gene_symbol[match(gene_id, gct_list_sc$gene_id)],
    significant = padj < 0.05 & !is.na(padj),
    significant_relaxed = padj < 0.1 & !is.na(padj)
  )

## Count significant genes
n_sig <- sum(df$significant, na.rm = TRUE)
n_sig_relax <- sum(df$significant_relaxed, na.rm = TRUE)
n_up <- sum(df$significant & df$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(df$significant & df$log2FoldChange < 0, na.rm = TRUE)

cat(sprintf("\nSignificant genes (padj < 0.05): %d\n", n_sig))
cat(sprintf("  Up in Post: %d\n", n_up))
cat(sprintf("  Down in Post: %d\n", n_down))
cat(sprintf("Significant genes (padj < 0.1): %d\n", n_sig_relax))

## -----------------------------------------------------------------------------
## Save results
## -----------------------------------------------------------------------------
cat("\nSaving results...\n")

readr::write_tsv(
  df %>%
    arrange(padj) %>%
    dplyr::select(gene_id, gene_symbol, baseMean, log2FoldChange, lfcSE,
           stat, pvalue, padj, significant, significant_relaxed),
  file.path(tab_dir, "DE_menopause_proxy_improved.tsv")
)

## Top significant genes
top_genes <- df %>%
  filter(significant) %>%
  arrange(padj) %>%
  dplyr::select(gene_id, gene_symbol, log2FoldChange, padj) %>%
  head(100)

readr::write_tsv(
  top_genes,
  file.path(tab_dir, "menopause_signature_candidate_genes.tsv")
)

cat(sprintf("\nSaved results to: %s\n", tab_dir))

## -----------------------------------------------------------------------------
## Create menopause signature gene set
## -----------------------------------------------------------------------------
cat("\nCreating menopause signature gene set...\n")

# Define signature genes (strict criteria)
signature_genes <- df %>%
  filter(
    significant,
    abs(log2FoldChange) > 1  # Minimum 2-fold change
  ) %>%
  arrange(padj) %>%
  dplyr::select(gene_symbol, log2FoldChange, padj)

# Save as gene list
writeLines(
  signature_genes$gene_symbol,
  file.path(tab_dir, "menopause_signature_gene_list.txt")
)

cat(sprintf("Signature size: %d genes\n", nrow(signature_genes)))
cat(sprintf("Saved gene list to: %s\n", file.path(tab_dir, "menopause_signature_gene_list.txt")))

## -----------------------------------------------------------------------------
## Visualizations
## -----------------------------------------------------------------------------
cat("\nCreating visualizations...\n")

## Volcano plot
df$neglog10p <- -log10(df$pvalue)
df$label <- ifelse(df$significant & abs(df$log2FoldChange) > 1,
                  df$gene_symbol, "")

p_volc <- ggplot(df, aes(x = log2FoldChange, y = neglog10p)) +
  geom_point(aes(color = significant), alpha = 0.4) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_text(aes(label = label), size = 3, vjust = -0.5, hjust = 0.5, check_overlap = TRUE) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70")) +
  labs(
    title = "Menopause Proxy DE Analysis: Improved Age Groups",
    subtitle = "Post-menopausal (55-70) vs Pre-menopausal (30-45), excluding 45-55",
    x = "log2 Fold Change (Post vs Pre)",
    y = "-log10(p-value)",
    color = "Significant"
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "volcano_menopause_proxy_improved.png"), p_volc, width = 8, height = 6, dpi = 300)

## Sample distribution plot
sample_dist <- meta_de %>%
  group_by(menopause_proxy, age_bin_label) %>%
  summarise(n = n(), .groups = "drop")

p_samples <- ggplot(sample_dist, aes(x = age_bin_label, y = n, fill = menopause_proxy)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Sample Distribution by Age and Menopause Proxy",
    x = "Age Group",
    y = "Sample Count",
    fill = "Menopause Proxy"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "sample_distribution.png"), p_samples, width = 8, height = 5, dpi = 300)

cat(sprintf("Saved figures to: %s\n", fig_dir))

## -----------------------------------------------------------------------------
## Summary
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("APPROACH 1: SUMMARY\n")
cat("===========================================\n\n")

cat(sprintf("Menopause signature candidates: %d genes (padj < 0.05)\n", n_sig))
cat(sprintf("  - Relaxed criteria (padj < 0.1): %d genes\n", n_sig_relax))
cat(sprintf("  - Signature genes (|FC| > 2): %d genes\n", nrow(signature_genes)))

cat("\nKey characteristics:\n")
cat("  - Uses improved age proxy (excludes 45-55 transition)\n")
cat("  - Includes SMCENTER batch correction\n")
cat("  - Focus on female subcutaneous adipose\n")

cat("\nOutput files:\n")
cat(sprintf("  Results: %s\n", file.path(tab_dir, "DE_menopause_proxy_improved.tsv")))
cat(sprintf("  Signature genes: %s\n", file.path(tab_dir, "menopause_signature_candidate_genes.tsv")))
cat(sprintf("  Gene list: %s\n", file.path(tab_dir, "menopause_signature_gene_list.txt")))
cat("===========================================\n\n")
