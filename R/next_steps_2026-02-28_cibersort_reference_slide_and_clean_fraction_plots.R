################################################################################
## Next steps (2026-02-28): CIBERSORT reference summary + clean fraction plots
##
## Goal:
## - Generate a reference-data summary plot for GSE176171 (cell counts per type).
## - Regenerate fraction plots excluding non-cell-type columns like P-value.
##
## Inputs:
## - external_data/GSE176171/derived/reference_cell_type_counts.tsv
## - external_data/GSE176171/derived/signature_topN50_cpm_noversion.tsv.gz
## - GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/*__fractions_with_meta.tsv.gz
##
## Outputs:
## - GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_cibersort_fractions/
################################################################################

.repos <- getOption("repos")
.cran <- unname(.repos["CRAN"])
if (is.null(.cran) || is.na(.cran) || .cran == "@CRAN@" || .cran == "") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
rm(.repos, .cran)

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readr, dplyr, tidyr, stringr, tibble, glue, ggplot2)

ref_counts <- "external_data/GSE176171/derived/reference_cell_type_counts.tsv"
sig_path <- "external_data/GSE176171/derived/signature_topN50_cpm_noversion.tsv.gz"

out_dir <- "GTEx_v10_AT_analysis_out"
fig_dir <- file.path(out_dir, "figs", "next_steps_2026-02-28_cibersort_fractions")
ensure_dir <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
ensure_dir(fig_dir)

# 1) Reference summary plot
if (file.exists(ref_counts)) {
  dfc <- readr::read_tsv(ref_counts, show_col_types = FALSE) %>%
    dplyr::arrange(dplyr::desc(n_cells))

  p_ref <- ggplot2::ggplot(dfc, ggplot2::aes(x = reorder(cell_type, n_cells), y = n_cells)) +
    ggplot2::geom_col(fill = "#4C78A8") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(
      title = "GSE176171 reference (human): cell counts by broad cell type",
      subtitle = "Used to build a compact signature matrix for CIBERSORT (first pass)",
      x = NULL,
      y = "# cells"
    )

  ggplot2::ggsave(file.path(fig_dir, "gse176171_reference_celltype_counts.png"), p_ref, width = 8.5, height = 5.0, dpi = 250)
}

# Cell-type columns should match signature matrix columns (preferred) OR fraction-like columns.
cell_types <- NULL
if (file.exists(sig_path)) {
  sig <- readr::read_tsv(sig_path, show_col_types = FALSE)
  cell_types <- colnames(sig)[-1]
}

# 2) Clean fraction plots for both depots
frac_paths <- c(
  subcutaneous = "GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/cibersort_subcutaneous_perm0_topN50_minCells200__fractions_with_meta.tsv.gz",
  visceral = "GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/cibersort_visceral_perm0_topN50_minCells200__fractions_with_meta.tsv.gz"
)

for (nm in names(frac_paths)) {
  path <- frac_paths[[nm]]
  if (!file.exists(path)) {
    warning("Missing: ", path)
    next
  }

  est2 <- readr::read_tsv(path, show_col_types = FALSE)

  if (!is.null(cell_types)) {
    frac_cols <- intersect(colnames(est2), cell_types)
  } else {
    # Fallback: drop known non-fraction columns.
    frac_cols <- setdiff(colnames(est2), c(
      "SAMPID", "Mixture", "P.value", "P-value", "P_value", "Correlation", "RMSE",
      "depot", "sex", "age_bin_label", "SMCENTER"
    ))
  }

  if (!length(frac_cols)) stop("No fraction columns found after filtering for: ", nm)

  sum_tbl <- est2 %>%
    dplyr::group_by(depot, sex, age_bin_label) %>%
    dplyr::summarise(
      n = dplyr::n(),
      dplyr::across(dplyr::all_of(frac_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )

  sum_long <- sum_tbl %>%
    tidyr::pivot_longer(cols = dplyr::all_of(frac_cols), names_to = "cell_type", values_to = "mean_fraction")

  p <- ggplot2::ggplot(sum_long, ggplot2::aes(x = age_bin_label, y = mean_fraction, fill = cell_type)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~ sex, nrow = 1) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = glue::glue("CIBERSORT fractions (mean) by sex x age bin | {nm}"),
      subtitle = "Cleaned: only cell-type fraction columns (no P-value/Correlation/RMSE)",
      x = "GTEx age bin",
      y = "Mean estimated fraction",
      fill = "Cell type"
    )

  ggplot2::ggsave(
    file.path(fig_dir, glue::glue("cibersort_{nm}_perm0_topN50_minCells200__fractions_mean_by_sex_agebin_clean.png")),
    p,
    width = 12.5,
    height = 4.2,
    dpi = 250
  )
}

cat("Wrote figures under: ", fig_dir, "\n", sep = "")
