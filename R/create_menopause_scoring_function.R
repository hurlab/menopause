################################################################################
## Menopause Signature Scoring Function (GSVA-like)
##
## Creates a scoring function to calculate menopause signature scores
## Can be applied to any GTEx sample or external dataset
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
  ggplot2, patchwork, GSVA, limma, genefilter
)

## -----------------------------------------------------------------------------
## Source utility functions
## -----------------------------------------------------------------------------
source("R/utils.R")

## -----------------------------------------------------------------------------
## Setup
## -----------------------------------------------------------------------------
in_dir   <- "GTExDatav10"
gct_sc   <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
sattr    <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph   <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

out_dir  <- "GTEx_v10_AT_analysis_out"
fig_dir  <- file.path(out_dir, "figs", "scoring_function")
tab_dir  <- file.path(out_dir, "tables", "scoring_function")
ensure_dirs(c(out_dir, fig_dir, tab_dir))

cat("\n===========================================\n")
cat("MENOPAUSE SIGNATURE SCORING FUNCTION\n")
cat("===========================================\n\n")

## -----------------------------------------------------------------------------
## Load menopause signature
## -----------------------------------------------------------------------------
cat("Loading menopause signature genes...\n")

sig_file <- file.path("GTEx_v10_AT_analysis_out/tables/approach1_improved_age_proxy/menopause_signature_gene_list.txt")
if (!file.exists(sig_file)) {
  stop("Signature file not found. Run Approach 1 first.")
}

signature_genes <- read_lines(sig_file)
cat(sprintf("Loaded %d signature genes\n", length(signature_genes)))

## -----------------------------------------------------------------------------
## Function to calculate menopause signature score
## -----------------------------------------------------------------------------
calculate_menopause_score <- function(counts, signature_genes, method = "zscore") {
  #' Calculate menopause signature score for each sample
  #'
  #' @param counts Gene expression matrix (genes x samples)
  #' @param signature_genes Vector of gene symbols to use for signature
  #' @param method Scoring method: "zscore", "gsva", "mean", "ssgsea"
  #' @return Named vector of scores (one per sample)

  ## Ensure gene symbols match
  rownames(counts) <- toupper(rownames(counts))

  ## For GSVA, we need Ensembll IDs - map gene symbols to Ensembl
  ## For simplicity, we'll use intersection
  common_genes <- intersect(toupper(signature_genes), rownames(counts))

  cat(sprintf("Found %d signature genes in data\n", length(common_genes)))

  if (length(common_genes) < 10) {
    stop("Too few signature genes found in data")
  }

  signature_expr <- counts[common_genes, , drop = FALSE]

  if (method == "zscore") {
    ## Z-score method (like SenMayo)
    gene_means <- rowMeans(signature_expr, na.rm = TRUE)
    gene_sds <- apply(signature_expr, 1, function(x) sd(x, na.rm = TRUE))

    ## Handle zero SD
    gene_sds[gene_sds == 0] <- 1

    z_scores <- (signature_expr - gene_means) / gene_sds
    sample_scores <- colMeans(z_scores, na.rm = TRUE)

  } else if (method == "mean") {
    ## Simple mean expression
    sample_scores <- colMeans(signature_expr, na.rm = TRUE)

  } else if (method == "gsva") {
    ## GSVA-like method
    gsva_mat <- gsva(
      as.matrix(signature_expr),
      list(signature = common_genes),
      method = "gsva",
      kcdf = "Gaussian",
      abs.ranking = FALSE
    )
    # GSVA returns a geneset x sample matrix.
    sample_scores <- as.numeric(gsva_mat["signature", ])
    names(sample_scores) <- colnames(gsva_mat)

  } else if (method == "ssgsea") {
    ## ssGSEA-like enrichment score
    sample_scores <- apply(signature_expr, 2, function(sample_expr) {
      ## Calculate rank-based enrichment
      ranks <- rank(sample_expr)
      max_rank <- max(ranks)

      ## Weighted average (higher rank = higher contribution)
      weighted_score <- sum(ranks * (max_rank - ranks + 1)) / sum(max_rank - ranks + 1)
      weighted_score / length(common_genes)
    })
  }

  return(sample_scores)
}

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

## Prepare counts
counts_sc <- gct_list_sc$counts[, meta$SAMPID, drop = FALSE]
counts_sc <- collapse_dupe_rowsum(counts_sc)

## Convert to gene symbols (and collapse duplicates by mean to avoid inflating signal)
sym <- toupper(gct_list_sc$gene_symbol[match(rownames(counts_sc), gct_list_sc$gene_id)])
keep_sym <- !is.na(sym) & nzchar(sym)
counts_sc <- counts_sc[keep_sym, , drop = FALSE]
sym <- sym[keep_sym]
if (any(duplicated(sym))) {
  summed <- rowsum(counts_sc, group = sym, reorder = FALSE)
  denom <- as.numeric(table(sym)[rownames(summed)])
  counts_sc <- summed / denom
  rownames(counts_sc) <- rownames(summed)
} else {
  rownames(counts_sc) <- sym
}

## Subset to female samples only
meta_f <- meta %>% filter(sex == "Female")
counts_f <- counts_sc[, meta_f$SAMPID, drop = FALSE]

cat(sprintf("Using %d female samples\n", ncol(counts_f)))

## -----------------------------------------------------------------------------
## Calculate menopause scores for all samples
## -----------------------------------------------------------------------------
cat("\nCalculating menopause signature scores...\n")

# Z-score method (primary - like SenMayo)
scores_zscore <- calculate_menopause_score(counts_f, signature_genes, method = "zscore")

# Add sample names
names(scores_zscore) <- meta_f$SAMPID

# Order by score
scores_ordered <- sort(scores_zscore, decreasing = TRUE)

## Create results data frame
score_df <- tibble::tibble(
  SAMPID = names(scores_zscore),
  menopause_score = scores_zscore,
  age_bin_label = meta_f$age_bin_label[match(names(scores_zscore), meta_f$SAMPID)],
  age_mid = meta_f$age_mid[match(names(scores_zscore), meta_f$SAMPID)]
)

## -----------------------------------------------------------------------------
## Analyze score vs age
## -----------------------------------------------------------------------------
cat("\nAnalyzing score-age correlation...\n")

age_score_cor <- cor.test(score_df$menopause_score, score_df$age_mid, use = "complete.obs")

cat(sprintf("Correlation with age: r = %.3f, p = %g\n", age_score_cor$estimate, age_score_cor$p.value))

## Score by age group
age_group_scores <- score_df %>%
  group_by(age_bin_label) %>%
  summarise(
    n = n(),
  mean_score = mean(menopause_score),
  sd_score = sd(menopause_score),
    se_score = sd(menopause_score) / sqrt(n()),
    .groups = "drop"
  )

cat("\nMean scores by age group:\n")
print(age_group_scores)

## -----------------------------------------------------------------------------
## Save results
## -----------------------------------------------------------------------------
cat("\nSaving results...\n")

## Save full scores
readr::write_tsv(
  score_df %>%
    arrange(desc(menopause_score)),
  file.path(tab_dir, "menopause_signature_scores_all_samples.tsv")
)

## Save summary by age
readr::write_tsv(
  age_group_scores,
  file.path(tab_dir, "menopause_score_by_age_group.tsv")
)

## Save signature genes for reference
sig_info <- tibble::tibble(
  gene_symbol = signature_genes,
  in_data = toupper(signature_genes) %in% rownames(counts_sc),
  note = "Menopause signature from Approach 1 (Improved Age Proxy)"
)
readr::write_tsv(sig_info, file.path(tab_dir, "menopause_signature_gene_info.tsv"))

## -----------------------------------------------------------------------------
## Visualizations
## -----------------------------------------------------------------------------
cat("\nCreating visualizations...\n")

## Score vs age scatter
p_scatter <- ggplot(score_df, aes(x = age_mid, y = menopause_score)) +
  geom_point(aes(color = age_bin_label), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(
    title = "Menopause Signature Score vs Age",
    subtitle = "Females, Subcutaneous Adipose",
    x = "Age (midpoint)",
    y = "Menopause Signature Score (Z-score)",
    color = "Age Group"
  ) +
  theme_bw() +
  labs(color = "Age Group")

ggsave(file.path(fig_dir, "menopause_score_vs_age.png"), p_scatter, width = 8, height = 6, dpi = 300)

## Score by age group boxplot
p_boxplot <- ggplot(score_df, aes(x = age_bin_label, y = menopause_score, fill = age_bin_label)) +
  geom_boxplot(outlier.alpha = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(
    title = "Menopause Signature Score by Age Group",
    x = "Age Group",
    y = "Menopause Signature Score (Z-score)",
    fill = "Age Group"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "menopause_score_by_age_group.png"), p_boxplot, width = 8, height = 6, dpi = 300)

## Distribution plot
p_dist <- ggplot(score_df, aes(x = menopause_score)) +
  geom_histogram(aes(fill = ..density..), bins = 30, alpha = 0.7) +
  geom_density(aes(color = "blue"), size = 1) +
  labs(
    title = "Distribution of Menopause Signature Scores",
    x = "Menopause Signature Score (Z-score)",
    y = "Density"
  ) +
  theme_bw()

ggsave(file.path(fig_dir, "menopause_score_distribution.png"), p_dist, width = 8, height = 5, dpi = 300)

cat(sprintf("Saved figures to: %s\n", fig_dir))

## -----------------------------------------------------------------------------
## Summary
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("SCORING FUNCTION SUMMARY\n")
cat("===========================================\n\n")

cat(sprintf("Signature size: %d genes\n", length(signature_genes)))
cat(sprintf("Samples scored: %d female samples\n", length(scores_zscore)))
cat(sprintf("Score-age correlation: r = %.3f (p = %g)\n", age_score_cor$estimate, age_score_cor$p.value))

cat("\nScoring function created:\n")
cat("  - Input: Gene expression matrix (genes x samples)\n")
cat("  - Signature: 81 genes from Approach 1\n")
cat("  - Method: Z-score (like SenMayo)\n")
cat("  - Output: Menopause score per sample\n")

cat("\nOutput files:\n")
cat(sprintf("  All scores: %s\n", file.path(tab_dir, "menopause_signature_scores_all_samples.tsv")))
cat(sprintf("  By age group: %s\n", file.path(tab_dir, "menopause_score_by_age_group.tsv")))
cat(sprintf("  Gene info: %s\n", file.path(tab_dir, "menopause_signature_gene_info.tsv")))

cat("\nUsage example:\n")
cat("  # Calculate scores for new data:\n")
cat("  scores <- calculate_menopause_score(counts_matrix, signature_genes)\n")
cat("  # Apply to GTEx other tissues:\n")
cat("  # scores_blood <- calculate_menopause_score(counts_blood, signature_genes)\n\n")

cat("===========================================\n\n")
