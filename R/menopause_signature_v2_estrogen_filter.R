################################################################################
## APPROACH 2: Estrogen-Responsive Gene Filter
##
## Strategy: Filter DE results for estrogen-responsive genes first
## to focus on biologically relevant menopause pathways
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
fig_dir  <- file.path(out_dir, "figs", "approach2_estrogen_filter")
tab_dir  <- file.path(out_dir, "tables", "approach2_estrogen_filter")
ensure_dirs(c(out_dir, fig_dir, tab_dir))

cat("\n===========================================\n")
cat("APPROACH 2: ESTROGEN-RESPONSIVE GENE FILTER\n")
cat("===========================================\n\n")
cat("Strategy:\n")
cat("  1. Load known estrogen-responsive gene sets\n")
cat("  2. Run DE on wider age groups (30-39 vs 60-69)\n")
cat("  3. Filter results for estrogen-responsive genes\n\n")

## -----------------------------------------------------------------------------
## Define estrogen-responsive gene sets
## -----------------------------------------------------------------------------
cat("Loading estrogen-responsive gene sets...\n")

## Known estrogen-responsive genes from literature
estrogen_genes_known <- c(
  "ESR1", "ESR2", "PGR", "GREB1", "TFF1", "TFF3", "CTSD", "FOXA1",
  "GATA3", "MYC", "CCND1", "c-MYC", "BCL2", "BAX", "FOS", "JUN",
  "EGFR", "ERBB2", "IGF1R", "VEGFA", "MMP1", "MMP2", "MMP9",
  "CDH1", "CDH2", "VIM", "SNAI1", "SNAI2", "TWIST1", "ZEB1",
  "KI67", "MKI67", "PCNA", "TOP2A", "CENPF", "AURKA", "PLK1"
)

## Get ESR1/ESR2 target genes from org.Hs.eg.db
cat("Querying ESR target genes from database...\n")
tryCatch({
  esr_targets <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = c("ESR1", "ESR2"),
    columns = c("SYMBOL", "ENTREZID", "ENSEMBL"),
    keytype = "SYMBOL"
  )

  # Get genes with estrogen response elements (EREs)
  # This is a simplified approach - ideally would use curated sets
  ere_genes <- c(
    "TFF1", "TFF3", "GREB1", "PGR", "CTSD", "CXCL12", "CXCL14",
    "IGFBP4", "IGFBP5", "CCND1", "c-MYC", "FOS", "JUN", "EGR1",
    "EGR3", "C3", "CDH1", "KRT19", "KRT8", "GATA3", "FOXA1"
  )

  cat(sprintf("Found %d known estrogen-responsive genes\n", length(estrogen_genes_known)))
}, error = function(e) {
  cat("Note: Could not query ESR targets, using known genes only\n")
  ere_genes <- character(0)
})

## Combine all estrogen-responsive genes
all_estrogen_genes <- unique(c(estrogen_genes_known, ere_genes))

## -----------------------------------------------------------------------------
## Load GTEx data
## -----------------------------------------------------------------------------
cat("\nLoading GTEx data...\n")

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
## Filter for young vs old (wider age groups)
## -----------------------------------------------------------------------------
meta_de <- meta %>%
  filter(
    sex == "Female",
    age_bin_label %in% c("30-39", "60-69")
  ) %>%
  mutate(
    age_group = factor(age_bin_label, levels = c("30-39", "60-69")),
    SMCENTER = factor(SMCENTER)
  )

cat(sprintf("\nAnalysis cohort:\n"))
cat(sprintf("  30-39: %d samples\n", sum(meta_de$age_group == "30-39")))
cat(sprintf("  60-69: %d samples\n", sum(meta_de$age_group == "60-69")))

## -----------------------------------------------------------------------------
## DESeq2 analysis
## -----------------------------------------------------------------------------
cat("\nRunning DESeq2...\n")

counts_sc <- gct_list_sc$counts[, meta_de$SAMPID, drop = FALSE]
counts_sc <- collapse_dupe_rowsum(counts_sc)

dds <- DESeqDataSetFromMatrix(
  countData = counts_sc,
  colData = meta_de,
  design = ~ SMCENTER + age_group
)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)
res <- results(dds, contrast = c("age_group", "60-69", "30-39"))

## -----------------------------------------------------------------------------
## Filter for estrogen-responsive genes
## -----------------------------------------------------------------------------
cat("\nFiltering for estrogen-responsive genes...\n")

df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    gene_symbol = gct_list_sc$gene_symbol[match(gene_id, gct_list_sc$gene_id)],
    estrogen_responsive = gene_symbol %in% all_estrogen_genes
  )

## Summary statistics
cat("\nAll genes:\n")
cat(sprintf("  Significant (padj < 0.05): %d\n", sum(df$padj < 0.05 & !is.na(df$padj))))

cat("\nEstrogen-responsive genes only:\n")
estrogen_df <- df %>%
  filter(estrogen_responsive)

cat(sprintf("  Total estrogen-responsive genes: %d\n", nrow(estrogen_df)))
cat(sprintf("  Significant (padj < 0.05): %d\n", sum(estrogen_df$padj < 0.05 & !is.na(estrogen_df$padj))))

## Save all results
readr::write_tsv(
  df %>%
    arrange(padj) %>%
    dplyr::select(gene_id, gene_symbol, estrogen_responsive, baseMean,
           log2FoldChange, lfcSE, stat, pvalue, padj),
  file.path(tab_dir, "DE_all_genes_with_estrogen_flag.tsv")
)

## Save estrogen-responsive genes
estrogen_sig <- estrogen_df %>%
  arrange(padj) %>%
  mutate(
    significant = padj < 0.05,
    significant_relaxed = padj < 0.1
  )

readr::write_tsv(
  estrogen_sig,
  file.path(tab_dir, "DE_estrogen_responsive_genes.tsv")
)

## Top estrogen-responsive genes
top_estrogen <- estrogen_df %>%
  arrange(padj) %>%
  head(50) %>%
  dplyr::select(gene_id, gene_symbol, log2FoldChange, padj)

readr::write_tsv(
  top_estrogen,
  file.path(tab_dir, "top_estrogen_responsive_genes.tsv")
)

## -----------------------------------------------------------------------------
## Create menopause signature from estrogen-responsive genes
## -----------------------------------------------------------------------------
cat("\nCreating estrogen-focused menopause signature...\n")

estrogen_signature <- estrogen_df %>%
  filter(
    padj < 0.05,
    abs(log2FoldChange) > 0.5  # Minimum 1.4-fold change
  ) %>%
  arrange(padj) %>%
  dplyr::select(gene_symbol, log2FoldChange, padj)

if (nrow(estrogen_signature) > 0) {
  writeLines(
    estrogen_signature$gene_symbol,
    file.path(tab_dir, "estrogen_menopause_signature.txt")
  )

  cat(sprintf("Estrogen-focused signature: %d genes\n", nrow(estrogen_signature)))
  cat(sprintf("Saved to: %s\n", file.path(tab_dir, "estrogen_menopause_signature.txt")))
} else {
  cat("No significant estrogen-responsive genes found.\n")
}

## -----------------------------------------------------------------------------
## Visualizations
## -----------------------------------------------------------------------------
cat("\nCreating visualizations...\n")

## Volcano plot colored by estrogen responsiveness
df$neglog10p <- -log10(df$pvalue)
df$is_sig <- df$padj < 0.05 & !is.na(df$padj)

p_volc <- ggplot(df, aes(x = log2FoldChange, y = neglog10p)) +
  geom_point(aes(color = estrogen_responsive, shape = is_sig), alpha = 0.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey70"), name = "Estrogen") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), name = "Significant") +
  labs(
    title = "DE Analysis with Estrogen-Responsive Genes Highlighted",
    subtitle = "Females 30-39 vs 60-69, Subcutaneous Adipose",
    x = "log2 Fold Change",
    y = "-log10(p-value)"
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "volcano_estrogen_highlight.png"), p_volc, width = 9, height = 6, dpi = 300)

## Bar plot of top estrogen genes
if (nrow(estrogen_sig) > 0) {
  top_n <- min(20, nrow(estrogen_sig))

  top_bar <- estrogen_sig %>%
    arrange(padj) %>%
    head(top_n) %>%
    mutate(gene_symbol = reorder(gene_symbol, log2FoldChange))

  p_bar <- ggplot(top_bar, aes(x = gene_symbol, y = log2FoldChange, fill = log2FoldChange > 0)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue"), name = "Up in 60-69") +
    labs(
      title = paste0("Top ", top_n, " Estrogen-Responsive DE Genes"),
      x = "Gene",
      y = "log2 Fold Change"
    ) +
    theme_bw()

  ggsave(file.path(fig_dir, "barplot_estrogen_genes.png"), p_bar, width = 8, height = 10, dpi = 300)
}

cat(sprintf("Saved figures to: %s\n", fig_dir))

## -----------------------------------------------------------------------------
## Summary
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("APPROACH 2: SUMMARY\n")
cat("===========================================\n\n")

cat(sprintf("Estrogen-responsive genes tested: %d\n", length(all_estrogen_genes)))
cat(sprintf("  Significant DE estrogen genes: %d\n", sum(estrogen_df$padj < 0.05 & !is.na(estrogen_df$padj))))
if (nrow(estrogen_signature) > 0) {
  cat(sprintf("  Estrogen-focused signature: %d genes\n", nrow(estrogen_signature)))
}

cat("\nKey characteristics:\n")
cat("  - Focuses on biologically relevant pathways\n")
cat("  - Uses estrogen-responsive gene filter\n")
cat("  - May reduce false positives from age-related noise\n")

cat("\nOutput files:\n")
cat(sprintf("  All results: %s\n", file.path(tab_dir, "DE_all_genes_with_estrogen_flag.tsv")))
cat(sprintf("  Estrogen genes: %s\n", file.path(tab_dir, "DE_estrogen_responsive_genes.tsv")))
cat("===========================================\n\n")
