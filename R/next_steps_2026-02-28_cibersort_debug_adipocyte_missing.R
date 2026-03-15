################################################################################
## Debug: Why is adipocyte fraction ~0 in GTEx adipose CIBERSORT results?
##
## This script:
## - Summarizes adipocyte fraction distribution.
## - Checks whether key adipocyte markers are present in the mixture/signature.
## - Tests a minimal intervention: remove MALAT1 (extreme, ubiquitous) and re-run
##   CIBERSORT to see if adipocyte fractions recover.
################################################################################

options(stringsAsFactors = FALSE)
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readr, dplyr, tidyr, stringr, tibble, glue, ggplot2, IOBR)

sig_path <- "external_data/GSE176171/derived/signature_topN50_cpm_noversion.tsv.gz"
feat_path <- "external_data/GSE176171/GSE176171_Hs10X.counts.features.tsv.gz"

# Input from the previous CIBERSORT runs (already on the same gene ID space)
run_dir <- "GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions"
sub_sig_used <- file.path(run_dir, "cibersort_subcutaneous_perm0_topN50_minCells200__signature_used.tsv.gz")
sub_mix_used <- file.path(run_dir, "cibersort_subcutaneous_perm0_topN50_minCells200__mixture_used.tsv.gz")
sub_frac_meta <- file.path(run_dir, "cibersort_subcutaneous_perm0_topN50_minCells200__fractions_with_meta.tsv.gz")

vis_sig_used <- file.path(run_dir, "cibersort_visceral_perm0_topN50_minCells200__signature_used.tsv.gz")
vis_mix_used <- file.path(run_dir, "cibersort_visceral_perm0_topN50_minCells200__mixture_used.tsv.gz")
vis_frac_meta <- file.path(run_dir, "cibersort_visceral_perm0_topN50_minCells200__fractions_with_meta.tsv.gz")

out_dir <- "GTEx_v10_AT_analysis_out"
tab_out <- file.path(out_dir, "tables", "next_steps_2026-02-28_cibersort_debug")
fig_out <- file.path(out_dir, "figs", "next_steps_2026-02-28_cibersort_debug")
dir.create(tab_out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(p) if (!file.exists(p)) stop("Missing: ", p)
lapply(c(sig_path, feat_path, sub_sig_used, sub_mix_used, sub_frac_meta, vis_sig_used, vis_mix_used, vis_frac_meta), stop_if_missing)

features <- readr::read_tsv(feat_path, col_names = c("feature_id", "symbol", "feature_type"), show_col_types = FALSE) %>%
  dplyr::filter(feature_type == "Gene Expression") %>%
  dplyr::mutate(id_nv = sub("\\..*$", "", feature_id))

id_to_symbol <- setNames(features$symbol, features$id_nv)

summarize_fractions <- function(df, label) {
  ct_cols <- intersect(colnames(df), c(
    "adipocyte","aspc","endothelial","lec","pericyte","smc","mesothelium",
    "macrophage","monocyte","dendritic_cell","mast_cell","t_cell","nk_cell","b_cell"
  ))
  s <- df %>%
    dplyr::summarise(
      depot = label,
      n = dplyr::n(),
      adipocyte_mean = mean(adipocyte, na.rm = TRUE),
      adipocyte_median = median(adipocyte, na.rm = TRUE),
      adipocyte_max = max(adipocyte, na.rm = TRUE),
      adipocyte_pct_zero = mean(adipocyte == 0, na.rm = TRUE),
      frac_sum_mean = mean(rowSums(dplyr::across(dplyr::all_of(ct_cols))), na.rm = TRUE)
    )
  s
}

sub_frac <- readr::read_tsv(sub_frac_meta, show_col_types = FALSE)
vis_frac <- readr::read_tsv(vis_frac_meta, show_col_types = FALSE)

sum_tbl <- dplyr::bind_rows(
  summarize_fractions(sub_frac, "subcutaneous"),
  summarize_fractions(vis_frac, "visceral")
)
readr::write_tsv(sum_tbl, file.path(tab_out, "fraction_summary_adipocyte.tsv"))

# Check whether canonical adipocyte genes are present among the genes used for the run
canonical_symbols <- c("ADIPOQ","PLIN1","LPL","FABP4","PPARG","LEP")
canon_ids <- names(id_to_symbol)[id_to_symbol %in% canonical_symbols]
canon_tbl <- tibble::tibble(ensg = canon_ids, symbol = id_to_symbol[canon_ids])
readr::write_tsv(canon_tbl, file.path(tab_out, "canonical_adipocyte_genes_ensg.tsv"))

sub_mix <- readr::read_tsv(sub_mix_used, show_col_types = FALSE)
sub_sig <- readr::read_tsv(sub_sig_used, show_col_types = FALSE)

present_in_mix <- canon_ids[canon_ids %in% sub_mix$Gene]
present_in_sig <- canon_ids[canon_ids %in% sub_sig$Gene]

present_tbl <- tibble::tibble(
  symbol = canonical_symbols,
  ensg = canon_ids[match(canonical_symbols, id_to_symbol[canon_ids])],
  in_sub_mix_used = canonical_symbols %in% id_to_symbol[present_in_mix],
  in_sub_sig_used = canonical_symbols %in% id_to_symbol[present_in_sig]
)
readr::write_tsv(present_tbl, file.path(tab_out, "canonical_presence_in_used_matrices.tsv"))

# Minimal intervention: remove MALAT1 (ENSG00000251562) from signature+mixture and re-run
rerun_no_malat1 <- function(sig_used, mix_used, label) {
  sig <- readr::read_tsv(sig_used, show_col_types = FALSE)
  mix <- readr::read_tsv(mix_used, show_col_types = FALSE)

  drop_id <- names(id_to_symbol)[id_to_symbol == "MALAT1"]
  if (length(drop_id) == 0) drop_id <- "ENSG00000251562"

  sig2 <- sig %>% dplyr::filter(Gene != drop_id)
  mix2 <- mix %>% dplyr::filter(Gene != drop_id)

  # Convert to matrices as expected by IOBR::CIBERSORT
  sig_mat <- as.data.frame(sig2[, -1, drop = FALSE])
  rownames(sig_mat) <- sig2$Gene
  mix_mat <- as.data.frame(mix2[, -1, drop = FALSE])
  rownames(mix_mat) <- mix2$Gene

  est <- IOBR::CIBERSORT(sig_matrix = sig_mat, mixture_file = mix_mat, perm = 0, QN = FALSE, absolute = FALSE)
  est_df <- as.data.frame(est)
  est_df$Mixture <- rownames(est_df)

  # Summarize adipocyte
  out_sum <- est_df %>%
    dplyr::summarise(
      depot = label,
      n = dplyr::n(),
      adipocyte_mean = mean(adipocyte, na.rm = TRUE),
      adipocyte_median = median(adipocyte, na.rm = TRUE),
      adipocyte_max = max(adipocyte, na.rm = TRUE),
      adipocyte_pct_zero = mean(adipocyte == 0, na.rm = TRUE)
    )

  readr::write_tsv(out_sum, file.path(tab_out, glue::glue("rerun_no_malat1__summary_{label}.tsv")))
  readr::write_tsv(est_df, file.path(tab_out, glue::glue("rerun_no_malat1__fractions_{label}.tsv.gz")))

  # Quick histogram
  p <- ggplot2::ggplot(est_df, ggplot2::aes(x = adipocyte)) +
    ggplot2::geom_histogram(bins = 60, fill = "#4C78A8") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = glue::glue("Adipocyte fraction after dropping MALAT1 ({label})"), x = "adipocyte fraction", y = "# samples")
  ggplot2::ggsave(file.path(fig_out, glue::glue("rerun_no_malat1__adipocyte_hist_{label}.png")), p, width = 7, height = 4, dpi = 250)

  out_sum
}

sub_r <- rerun_no_malat1(sub_sig_used, sub_mix_used, "subcutaneous")
vis_r <- rerun_no_malat1(vis_sig_used, vis_mix_used, "visceral")

cat("Original fractions summary written to: ", file.path(tab_out, "fraction_summary_adipocyte.tsv"), "\n", sep="")
cat("Re-run (no MALAT1) summaries:\n")
print(dplyr::bind_rows(sub_r, vis_r))
