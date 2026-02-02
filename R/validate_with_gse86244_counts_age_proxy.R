################################################################################
## Validate GTEx Menopause Signature with GSE86244 (Counts + Age Proxy)
##
## Context:
## - GSE86244 is an RNA-seq series with processed count matrices provided as GEO
##   supplementary files (not a simple microarray ExpressionSet).
## - Menopause status is not explicitly encoded; we use the age proxy described
##   in literature: pre (<45) vs post (>55), excluding 45–55.
##
## Output:
## - GTEx_v10_AT_analysis_out/tables/gse86244_validation_counts_age_proxy/
## - GTEx_v10_AT_analysis_out/figs/gse86244_validation_counts_age_proxy/
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
  readr, dplyr, tidyr, stringr, tibble, purrr,
  ggplot2,
  DESeq2, apeglm,
  GEOquery
)

source("R/utils.R")

## -----------------------------------------------------------------------------
## Setup
## -----------------------------------------------------------------------------
out_dir <- "GTEx_v10_AT_analysis_out"
fig_dir <- file.path(out_dir, "figs", "gse86244_validation_counts_age_proxy")
tab_dir <- file.path(out_dir, "tables", "gse86244_validation_counts_age_proxy")
ensure_dirs(c(fig_dir, tab_dir))

cat("\n===========================================\n")
cat("GSE86244 VALIDATION (COUNTS + AGE PROXY)\n")
cat("===========================================\n\n")

## Load GTEx signature (Approach 1)
sig_file <- file.path(out_dir, "tables", "approach1_improved_age_proxy", "menopause_signature_gene_list.txt")
if (!file.exists(sig_file)) {
  stop("Signature file not found. Run Approach 1 first: ", sig_file)
}
signature_genes <- readr::read_lines(sig_file) %>% stringr::str_trim() %>% toupper()
signature_genes <- signature_genes[nzchar(signature_genes)]
cat(sprintf("Loaded GTEx signature: %d genes\n", length(signature_genes)))

## -----------------------------------------------------------------------------
## Download + load GSE86244 supplementary count matrix
## -----------------------------------------------------------------------------
gse_dir <- file.path(out_dir, "GSE86244_data")
ensure_dirs(gse_dir)

counts_file <- file.path(
  gse_dir,
  "GSE86244_FINAL_master_list_of_genes_counts_MIN.RNA_Pooling_Experiment_29_samples.highExp.txt.gz"
)

if (!file.exists(counts_file)) {
  cat("Downloading GSE86244 supplementary count matrices...\n")
  GEOquery::getGEOSuppFiles("GSE86244", baseDir = gse_dir, makeDirectory = FALSE)
}
if (!file.exists(counts_file)) {
  stop("Missing expected supplementary file: ", counts_file)
}

cat("Reading: ", counts_file, "\n", sep = "")
df <- readr::read_tsv(counts_file, show_col_types = FALSE)

if (!all(c("id", "geneSymbol") %in% names(df))) {
  stop("Unexpected columns in counts file; expected at least: id, geneSymbol")
}

sample_cols <- names(df)[stringr::str_detect(names(df), "^[^-]+-[0-9]+-")]
if (!length(sample_cols)) {
  stop("Could not identify sample columns from header names.")
}

ages <- stringr::str_match(sample_cols, "^[^-]+-([0-9]+)-")[, 2]
ages <- as.numeric(ages)
if (anyNA(ages)) {
  stop("Failed to parse ages from sample IDs in counts file.")
}

group <- dplyr::case_when(
  ages < 45 ~ "pre",
  ages > 55 ~ "post",
  TRUE ~ NA_character_
)

meta <- tibble::tibble(
  sample_id = sample_cols,
  age = ages,
  menopause_proxy = group
)

meta_use <- meta %>% dplyr::filter(!is.na(.data$menopause_proxy)) %>%
  dplyr::mutate(menopause_proxy = factor(.data$menopause_proxy, levels = c("pre", "post")))

cat(sprintf("GSE86244 samples in count file: %d\n", nrow(meta)))
cat("Menopause proxy counts (age):\n")
print(table(meta_use$menopause_proxy, useNA = "ifany"))

counts <- as.matrix(df[, meta_use$sample_id, drop = FALSE])
storage.mode(counts) <- "integer"

gene_id <- sub("^gene:", "", as.character(df$id))
rownames(counts) <- gene_id
counts <- collapse_dupe_rowsum(counts)

gene_symbol_lookup <- toupper(as.character(df$geneSymbol))
names(gene_symbol_lookup) <- gene_id

## -----------------------------------------------------------------------------
## DESeq2 (post vs pre)
## -----------------------------------------------------------------------------
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts,
  colData = data.frame(menopause_proxy = meta_use$menopause_proxy, row.names = meta_use$sample_id),
  design = ~ menopause_proxy
)

keep <- rowSums(DESeq2::counts(dds)) >= 10
dds <- dds[keep, ]

dds <- DESeq2::DESeq(dds)
res <- DESeq2::results(dds, contrast = c("menopause_proxy", "post", "pre"))
res_shr <- tryCatch(
  DESeq2::lfcShrink(dds, coef = "menopause_proxy_post_vs_pre", res = res, type = "apeglm"),
  error = function(e) res
)

de <- as.data.frame(res_shr) %>%
  tibble::rownames_to_column("gene_id") %>%
  dplyr::mutate(
    gene_symbol = gene_symbol_lookup[sub("\\..*$", "", .data$gene_id)],
    gene_symbol = ifelse(is.na(.data$gene_symbol) | !nzchar(.data$gene_symbol), NA_character_, .data$gene_symbol)
  ) %>%
  dplyr::arrange(.data$padj)

readr::write_tsv(de, file.path(tab_dir, "GSE86244_DESeq2_post_vs_pre_all.tsv"))

de_sig <- de %>%
  dplyr::filter(!is.na(.data$padj), .data$padj < 0.05) %>%
  dplyr::mutate(gene_symbol = toupper(.data$gene_symbol))

readr::write_tsv(de_sig, file.path(tab_dir, "GSE86244_DESeq2_post_vs_pre_significant.tsv"))

de_up <- de_sig %>% dplyr::filter(.data$log2FoldChange > 0, !is.na(.data$gene_symbol))
de_down <- de_sig %>% dplyr::filter(.data$log2FoldChange < 0, !is.na(.data$gene_symbol))

writeLines(unique(de_up$gene_symbol), file.path(tab_dir, "GSE86244_post_up_gene_symbols.txt"))
writeLines(unique(de_down$gene_symbol), file.path(tab_dir, "GSE86244_post_down_gene_symbols.txt"))

cat(sprintf("DESig (padj<0.05): %d genes\n", nrow(de_sig)))
cat(sprintf("  Up in post: %d\n", nrow(de_up)))
cat(sprintf("  Down in post: %d\n", nrow(de_down)))

## -----------------------------------------------------------------------------
## Overlap vs GTEx signature
## -----------------------------------------------------------------------------
gse_sig_syms <- unique(na.omit(de_sig$gene_symbol))
overlap <- intersect(signature_genes, gse_sig_syms)

summary_tbl <- tibble::tibble(
  gtex_signature_n = length(signature_genes),
  gse_sig_n = length(gse_sig_syms),
  overlap_n = length(overlap),
  overlap_frac_of_gtex = length(overlap) / length(signature_genes)
)
readr::write_tsv(summary_tbl, file.path(tab_dir, "GTEx_vs_GSE86244_overlap_summary.tsv"))

overlap_tbl <- tibble::tibble(gene_symbol = overlap) %>%
  dplyr::left_join(de_sig %>% dplyr::select(gene_symbol, log2FoldChange, padj), by = "gene_symbol") %>%
  dplyr::arrange(.data$padj)
readr::write_tsv(overlap_tbl, file.path(tab_dir, "GTEx_vs_GSE86244_overlap_genes.tsv"))

## Fisher enrichment (background = tested genes with non-missing symbol)
bg_syms <- unique(na.omit(toupper(de$gene_symbol)))
a <- length(overlap)
b <- length(signature_genes) - a
c <- length(gse_sig_syms) - a
d <- length(setdiff(bg_syms, union(signature_genes, gse_sig_syms)))
cont <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
rownames(cont) <- c("InGTExSig", "NotInGTExSig")
colnames(cont) <- c("InGSE86244Sig", "NotInGSE86244Sig")
fisher <- fisher.test(cont)

fisher_tbl <- tibble::tibble(
  odds_ratio = unname(fisher$estimate),
  p_value = fisher$p.value
)
readr::write_tsv(fisher_tbl, file.path(tab_dir, "GTEx_vs_GSE86244_fisher.tsv"))

## -----------------------------------------------------------------------------
## Simple visualization
## -----------------------------------------------------------------------------
plot_df <- tibble::tibble(
  set = c("GTEx signature", "GSE86244 sig (age proxy)", "Overlap"),
  n = c(length(signature_genes), length(gse_sig_syms), length(overlap))
)

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = set, y = n)) +
  ggplot2::geom_col() +
  ggplot2::theme_bw() +
  ggplot2::labs(
    title = "GTEx menopause signature overlap with GSE86244 (age proxy)",
    subtitle = sprintf("Fisher OR=%.2f, p=%.3g", fisher_tbl$odds_ratio, fisher_tbl$p_value),
    x = NULL,
    y = "Gene count"
  )
ggplot2::ggsave(file.path(fig_dir, "overlap_counts.png"), p, width = 7, height = 4.5, dpi = 200)

cat("\nWrote tables under: ", tab_dir, "\n", sep = "")
cat("Wrote figures under: ", fig_dir, "\n", sep = "")

