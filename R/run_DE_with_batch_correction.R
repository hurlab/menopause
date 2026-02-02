################################################################################
## Differential Expression with Proper Batch Correction
##
## Strategy:
## 1. DESeq2 with SMCENTER as covariate (for DE results)
## 2. ComBat on VST data (for visualization only)
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
  ggplot2, patchwork, DESeq2, apeglm, vsn, sva, Matrix, matrixStats, AnnotationDbi, org.Hs.eg.db
)

## Load utility functions
source("R/utils.R")

## -----------------------------------------------------------------------------
## Paths
## -----------------------------------------------------------------------------
in_dir   <- "GTExDatav10"
gct_sc   <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
gct_vis  <- file.path(in_dir, "gene_reads_v10_adipose_visceral_omentum.gct.gz")
sattr    <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph   <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

out_dir  <- "GTEx_v10_AT_analysis_out"
fig_dir  <- file.path(out_dir, "figs")
tab_dir  <- file.path(out_dir, "tables")
rds_dir  <- file.path(out_dir, "rds")
ensure_dirs(c(out_dir, fig_dir, tab_dir, rds_dir))

cat("\n===========================================\n")
cat("DIFFERENTIAL EXPRESSION WITH BATCH CORRECTION\n")
cat("===========================================\n\n")

## -----------------------------------------------------------------------------
## Load and prepare data
## -----------------------------------------------------------------------------
cat("Loading data...\n")

## Load sample attributes
sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything())

## Load subject phenotypes (contains SEX and AGE bins)
subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything())

## Extract subject ID from SAMPID (format: GTEX-XXXXX-XXXX-SM-XXXXX)
sample_attr <- sample_attr %>%
  dplyr::mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])

## Merge sample attributes with subject phenotypes
meta_full <- sample_attr %>%
  dplyr::left_join(subject_pheno, by = "SUBJID")

## Load gene expression data
gct_list_sc <- read_gct_v12(gct_sc)
gct_list_vis <- read_gct_v12(gct_vis)

## Samples for each depot matrix (these are typically disjoint between depots)
sc_samples_all <- intersect(colnames(gct_list_sc$counts), sample_attr$SAMPID)
vis_samples_all <- intersect(colnames(gct_list_vis$counts), sample_attr$SAMPID)
adipose_samples <- union(sc_samples_all, vis_samples_all)

cat(sprintf("Found %d subcutaneous samples and %d visceral samples\n",
            length(sc_samples_all), length(vis_samples_all)))

## Filter to adipose tissue samples
meta <- meta_full %>%
  filter(SAMPID %in% adipose_samples, SMTS == "Adipose Tissue") %>%
  mutate(
    DEPOT = ifelse(SMTSD == "Adipose - Subcutaneous", "Subcutaneous",
                   ifelse(SMTSD == "Adipose - Visceral (Omentum)", "Visceral", "Other")),
    sex = ifelse(SEX == 1, "Male", ifelse(SEX == 2, "Female", NA))
  ) %>%
  filter(DEPOT %in% c("Subcutaneous", "Visceral"))

## Add age variables
meta <- meta %>%
  mutate(age_bin_label = AGE) %>%
  bind_cols(parse_age_bin(.$AGE)) %>%
  derive_age_numeric()

cat(sprintf("Filtered to %d adipose samples\n", nrow(meta)))

cat("\nSex distribution (adipose samples):\n")
print(table(meta$sex, useNA = "ifany"))
cat("\nAge-bin distribution (adipose samples):\n")
print(table(meta$AGE, useNA = "ifany"))

## -----------------------------------------------------------------------------
## Filter samples: Females, age 40-49 vs 50-59, Subcutaneous only
## -----------------------------------------------------------------------------
cat("\nFiltering for DE analysis...\n")

meta_de <- meta %>%
  filter(
    sex == "Female",
    DEPOT == "Subcutaneous",
    age_bin_label %in% c("40-49", "50-59")
  )

cat(sprintf("  Females 40-49: %d samples\n", sum(meta_de$age_bin_label == "40-49")))
cat(sprintf("  Females 50-59: %d samples\n", sum(meta_de$age_bin_label == "50-59")))

## Check batch distribution
cat("\nBatch distribution in DE cohort:\n")
batch_table <- table(meta_de$age_bin_label, meta_de$SMCENTER)
print(batch_table)

if (sum(batch_table) > 0 && nrow(batch_table) > 1 && ncol(batch_table) > 1) {
  chi_test <- chisq.test(batch_table)
  cat(sprintf("\nChi-square p-value: %g\n", chi_test$p.value))
} else {
  cat("\nChi-square test skipped (insufficient samples in DE cohort).\n")
}

## -----------------------------------------------------------------------------
## Prepare count matrices
## -----------------------------------------------------------------------------
cat("\nPreparing count matrices...\n")

## Get samples for each depot
sc_samples <- meta %>% filter(DEPOT == "Subcutaneous") %>% pull(SAMPID)
vis_samples <- meta %>% filter(DEPOT == "Visceral") %>% pull(SAMPID)

## Subset count matrices
counts_sc <- gct_list_sc$counts[, intersect(sc_samples, colnames(gct_list_sc$counts)), drop = FALSE]
counts_vis <- gct_list_vis$counts[, intersect(vis_samples, colnames(gct_list_vis$counts)), drop = FALSE]

## Collapse duplicate genes
counts_sc <- collapse_dupe_rowsum(counts_sc)
counts_vis <- collapse_dupe_rowsum(counts_vis)

cat(sprintf("  Subcutaneous: %d genes, %d samples\n", nrow(counts_sc), ncol(counts_sc)))
cat(sprintf("  Visceral: %d genes, %d samples\n", nrow(counts_vis), ncol(counts_vis)))

## -----------------------------------------------------------------------------
## ANALYSIS 1: DESeq2 with SMCENTER as covariate (PRIMARY)
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("ANALYSIS 1: DESeq2 with Batch Covariate\n")
cat("===========================================\n\n")

## Prepare DESeq2 dataset for subcutaneous depot
meta_sc <- meta %>% filter(DEPOT == "Subcutaneous")

## Create age group factor
meta_sc$age_group <- factor(meta_sc$age_bin_label, levels = c("40-49", "50-59"))
meta_sc$SMCENTER <- factor(meta_sc$SMCENTER)

## Subset to females with age 40-59
keep_rows <- meta_sc$SAMPID %in% meta_de$SAMPID
meta_sc_de <- meta_sc[keep_rows, ]
counts_sc_de <- counts_sc[, meta_sc_de$SAMPID, drop = FALSE]

cat(sprintf("DESeq2 analysis: %d genes, %d samples\n", nrow(counts_sc_de), ncol(counts_sc_de)))

## Show batch distribution
cat("\nDesign formula with batch covariate:\n")
cat("  design = ~ SMCENTER + age_group\n\n")

## Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts_sc_de,
  colData = meta_sc_de,
  design = ~ SMCENTER + age_group
)

## Pre-filtering: keep genes with >= 10 counts total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat(sprintf("After filtering: %d genes\n", nrow(dds)))

## Run DESeq
cat("\nRunning DESeq (this may take a few minutes)...\n")
dds <- DESeq(dds)

## Get results
res <- results(dds, contrast = c("age_group", "50-59", "40-49"))
coef_name <- grep("^age_group", resultsNames(dds), value = TRUE)
res <- tryCatch(
  {
    if (length(coef_name) < 1) {
      stop("No age_group coefficient found in resultsNames(dds)")
    }
    lfcShrink(dds, coef = coef_name[[1]], type = "apeglm")
  },
  error = function(e) {
    message("Note: lfcShrink failed, using unshrunk results. Reason: ", e$message)
    res
  }
)

## Summary
cat("\nDESeq2 Results Summary:\n")
summary(res)

## Save results
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    gene_symbol = gct_list_sc$gene_symbol[match(gene_id, gct_list_sc$gene_id)],
    significant = padj < 0.05 & !is.na(padj)
  )

## Write results
readr::write_tsv(
  res_df %>%
    arrange(padj) %>%
    dplyr::select(dplyr::any_of(c(
      "gene_id", "gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "significant"
    ))),
  file.path(tab_dir, "DE_results_with_batch_covariate.tsv")
)

cat(sprintf("\nSaved results to: %s\n", file.path(tab_dir, "DE_results_with_batch_covariate.tsv")))

## Count significant genes
n_sig <- sum(res_df$significant, na.rm = TRUE)
n_up <- sum(res_df$significant & res_df$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(res_df$significant & res_df$log2FoldChange < 0, na.rm = TRUE)

cat(sprintf("\nSignificant genes (padj < 0.05): %d\n", n_sig))
cat(sprintf("  Up in 50-59: %d\n", n_up))
cat(sprintf("  Down in 50-59: %d\n", n_down))

## -----------------------------------------------------------------------------
## ANALYSIS 2: ComBat on VST data (for visualization only)
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("ANALYSIS 2: ComBat on VST Data (Visualization)\n")
cat("===========================================\n\n")

## Create VST for ALL subcutaneous samples (not just DE cohort)
meta_sc_all <- meta %>% filter(DEPOT == "Subcutaneous")
meta_sc_all$SMCENTER <- factor(meta_sc_all$SMCENTER)

dds_vst <- DESeqDataSetFromMatrix(
  countData = counts_sc,
  colData = meta_sc_all,
  design = ~ SMCENTER  # Minimal design for VST
)

## Filter
keep <- rowSums(counts(dds_vst)) >= 10
dds_vst <- dds_vst[keep, ]

## VST transform
cat("Running VST transformation...\n")
vsd <- vst(dds_vst, blind = FALSE)

## Extract VST matrix
vst_mat <- assay(vsd)

## Apply ComBat
cat("Applying ComBat correction...\n")
vst_combat <- ComBat(
  dat = vst_mat,
  batch = meta_sc_all$SMCENTER,
  mod = NULL,
  par.prior = TRUE,
  prior.plots = FALSE
)

## Create comparison PCA plots
cat("\nCreating PCA comparison plots...\n")

## PCA on original VST (before ComBat)
pca_orig <- prcomp(t(vst_mat), center = TRUE, scale. = FALSE)
df_orig <- pca_scores_df(pca_orig, meta_sc_all, c("DEPOT", "sex", "age_mid", "SMCENTER"))

varpct <- (pca_orig$sdev^2 / sum(pca_orig$sdev^2)) * 100

p_orig <- ggplot(df_orig, aes(PC1, PC2, color = SMCENTER, shape = DEPOT)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "PCA BEFORE ComBat Correction",
    x = glue::glue("PC1 ({round(varpct[1], 1)}%)"),
    y = glue::glue("PC2 ({round(varpct[2], 1)}%)"),
    color = "Sequencing Center"
  ) +
  theme_bw()

## PCA on ComBat-corrected VST
pca_combat <- prcomp(t(vst_combat), center = TRUE, scale. = FALSE)
df_combat <- pca_scores_df(pca_combat, meta_sc_all, c("DEPOT", "sex", "age_mid", "SMCENTER"))

varpct_combat <- (pca_combat$sdev^2 / sum(pca_combat$sdev^2)) * 100

p_combat <- ggplot(df_combat, aes(PC1, PC2, color = age_mid, shape = DEPOT)) +
  geom_point(size = 2, alpha = 0.7) +
  labs(
    title = "PCA AFTER ComBat Correction (colored by age)",
    x = glue::glue("PC1 ({round(varpct_combat[1], 1)}%)"),
    y = glue::glue("PC2 ({round(varpct_combat[2], 1)}%)"),
    color = "Age (midpoint)"
  ) +
  theme_bw()

## Combined plot
p_combined <- p_orig / p_combat

ggsave(
  file.path(fig_dir, "PCA_combat_correction_comparison.png"),
  p_combined,
  width = 10,
  height = 12,
  dpi = 300
)

cat(sprintf("Saved PCA comparison to: %s\n", file.path(fig_dir, "PCA_combat_correction_comparison.png")))

## -----------------------------------------------------------------------------
## SUMMARY
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("SUMMARY\n")
cat("===========================================\n\n")

cat("Analysis complete!\n\n")

cat("Results saved:\n")
cat(sprintf("  1. DE results (with batch covariate): %s\n", file.path(tab_dir, "DE_results_with_batch_covariate.tsv")))
cat(sprintf("  2. PCA comparison: %s\n", file.path(fig_dir, "PCA_combat_correction_comparison.png")))

cat("\nKey findings:\n")
cat(sprintf("  - Significant DE genes: %d\n", n_sig))
cat(sprintf("  - Up in 50-59: %d\n", n_up))
cat(sprintf("  - Down in 50-59: %d\n", n_down))

cat("\nInterpretation:\n")
cat("  - DESeq2 results account for batch effects (SMCENTER) in the model\n")
cat("  - ComBat-corrected PCA shows data after batch removal\n")
cat("  - Compare the two PCA plots to visualize batch effect impact\n\n")

cat("RECOMMENDED NEXT STEPS:\n")
cat("  1. Examine top DE genes for biological relevance\n")
cat("  2. Compare with results WITHOUT batch correction\n")
cat("  3. If results differ substantially, batch correction is critical\n")
cat("===========================================\n\n")
