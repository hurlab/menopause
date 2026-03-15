################################################################################
## Plot helper: build mean-by-sex-agebin fraction stacked bars from an existing
## __fractions_with_meta.tsv.gz file.
################################################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readr, dplyr, tidyr, stringr, tibble, glue, ggplot2)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript --vanilla R/next_steps_2026-02-28_cibersort_plot_from_fractions_with_meta.R <fractions_with_meta.tsv.gz> <out_png>")
}

in_path <- args[[1]]
out_png <- args[[2]]

df <- readr::read_tsv(in_path, show_col_types = FALSE)

meta_cols <- c("SAMPID","depot","sex","age_bin_label","SMCENTER")
stat_cols <- c("P.value","P-value","Correlation","RMSE")
frac_cols <- setdiff(colnames(df), c(meta_cols, stat_cols))
# Keep only numeric columns
frac_cols <- frac_cols[sapply(df[frac_cols], is.numeric)]

sum_tbl <- df %>%
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
    title = glue::glue("CIBERSORT fractions (mean) by sex x age bin | {unique(df$depot)}"),
    subtitle = glue::glue("Input: {basename(in_path)}"),
    x = "GTEx age bin",
    y = "Mean estimated fraction",
    fill = "Cell type"
  )

ggplot2::ggsave(out_png, p, width = 12.5, height = 4.2, dpi = 250)
