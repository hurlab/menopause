################################################################################
## Compare reference cell-count proportions vs CIBERSORT mean fractions
## (using corrected noMALAT1 runs)
################################################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readr, dplyr, tidyr, stringr, tibble, glue, ggplot2)

ref_counts <- "external_data/GSE176171/derived/reference_cell_type_counts.tsv"
sub_frac <- "GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/cibersort_subcutaneous_perm0_topN50_minCells200_noMALAT1__fractions_with_meta.tsv.gz"
vis_frac <- "GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/cibersort_visceral_perm0_topN50_minCells200_noMALAT1__fractions_with_meta.tsv.gz"

out_dir <- "GTEx_v10_AT_analysis_out"
fig_dir <- file.path(out_dir, "figs", "next_steps_2026-02-28_cibersort_debug")
tab_dir <- file.path(out_dir, "tables", "next_steps_2026-02-28_cibersort_debug")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

ref <- readr::read_tsv(ref_counts, show_col_types = FALSE) %>%
  dplyr::mutate(ref_pct = n_cells / sum(n_cells)) %>%
  dplyr::select(cell_type, ref_pct)

load_mean <- function(path, depot_label) {
  df <- readr::read_tsv(path, show_col_types = FALSE)
  meta_cols <- c("SAMPID","depot","sex","age_bin_label","SMCENTER")
  stat_cols <- c("P.value","P-value","Correlation","RMSE")
  frac_cols <- setdiff(colnames(df), c(meta_cols, stat_cols))
  frac_cols <- frac_cols[sapply(df[frac_cols], is.numeric)]

  means <- df %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(frac_cols), ~ mean(.x, na.rm = TRUE))) %>%
    tidyr::pivot_longer(cols = dplyr::everything(), names_to = "cell_type", values_to = "cib_mean") %>%
    dplyr::mutate(depot = depot_label)
  means
}

sub_m <- load_mean(sub_frac, "Subcutaneous")
vis_m <- load_mean(vis_frac, "Visceral")

cmp <- dplyr::bind_rows(sub_m, vis_m) %>%
  dplyr::left_join(ref, by = "cell_type") %>%
  dplyr::mutate(ref_pct = ifelse(is.na(ref_pct), 0, ref_pct))

# Correlation per depot
cors <- cmp %>%
  dplyr::group_by(depot) %>%
  dplyr::summarise(r = suppressWarnings(stats::cor(ref_pct, cib_mean)), .groups = "drop")

readr::write_tsv(cmp, file.path(tab_dir, "ref_pct_vs_cibersort_mean_noMALAT1.tsv"))
readr::write_tsv(cors, file.path(tab_dir, "ref_pct_vs_cibersort_mean_noMALAT1_cor.tsv"))

# Scatter plot with labels
p <- ggplot2::ggplot(cmp, ggplot2::aes(x = ref_pct, y = cib_mean, label = cell_type)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, color = "#999999", linewidth = 0.4) +
  ggplot2::geom_point(size = 2.0, color = "#4C78A8", alpha = 0.9) +
  ggplot2::geom_text(size = 3.0, vjust = -0.6, check_overlap = TRUE) +
  ggplot2::facet_wrap(~ depot, nrow = 1) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
  ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  ggplot2::labs(
    title = "Reference cell-count % vs CIBERSORT mean fraction (noMALAT1)",
    subtitle = paste0(
      "GSE176171 reference (human cells kept) vs GTEx adipose bulk deconvolution. ",
      paste(cors$depot, sprintf("r=%.2f", cors$r), collapse = "; ")
    ),
    x = "Reference proportion (cell counts)",
    y = "CIBERSORT mean fraction"
  )

ggplot2::ggsave(file.path(fig_dir, "ref_pct_vs_cibersort_mean_noMALAT1.png"), p, width = 12.0, height = 4.5, dpi = 250)
