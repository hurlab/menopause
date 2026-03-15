################################################################################
## Next steps (2026-02-28): Kenichi panel with F/M interleaved groups
##
## Goal:
## - Use already-generated per-sample score tables and re-plot so that
##   x-axis order is: F-pre, M-pre, F-peri, M-peri, F-post, M-post.
## - Produce one figure per depot (subq + visceral), matching Slide 7–8 style.
##
## Inputs:
## - GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_kenichi_panel_by_sex/
##
## Outputs:
## - GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_kenichi_panel_by_sex/
################################################################################

.repos <- getOption("repos")
.cran <- unname(.repos["CRAN"])
if (is.null(.cran) || is.na(.cran) || .cran == "@CRAN@" || .cran == "") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
rm(.repos, .cran)

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readr, dplyr, tidyr, stringr, tibble, glue, ggplot2, patchwork)

out_dir <- "GTEx_v10_AT_analysis_out"
tab_dir <- file.path(out_dir, "tables", "next_steps_2026-02-28_kenichi_panel_by_sex")
fig_dir <- file.path(out_dir, "figs", "next_steps_2026-02-28_kenichi_panel_by_sex")

dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

plot_one <- function(path_tsv, suffix, depot_label) {
  df0 <- readr::read_tsv(path_tsv, show_col_types = FALSE) %>%
    dplyr::filter(!is.na(sex), !is.na(menopause_group_3)) %>%
    dplyr::mutate(
      sex_short = dplyr::case_when(sex == "Female" ~ "F", sex == "Male" ~ "M", TRUE ~ as.character(sex)),
      group_fm = paste0(sex_short, "-", as.character(menopause_group_3)),
      group_fm = factor(group_fm, levels = c("F-pre", "M-pre", "F-peri", "M-peri", "F-post", "M-post"))
    )

  # Keep only score columns we know exist
  score_cols <- intersect(c("score__lipolysis_core", "score__thermogenesis_program"), colnames(df0))
  if (!length(score_cols)) stop("No expected score columns found in: ", path_tsv)

  make_panel <- function(col, title) {
    ggplot2::ggplot(df0, ggplot2::aes(x = group_fm, y = .data[[col]], color = sex_short)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.25) +
      ggplot2::geom_jitter(width = 0.18, alpha = 0.35, size = 0.9) +
      ggplot2::stat_summary(
        fun = median,
        geom = "text",
        ggplot2::aes(label = sprintf("%.2f", after_stat(y))),
        color = "black",
        size = 3.0,
        vjust = -0.7
      ) +
      ggplot2::theme_bw() +
      ggplot2::scale_color_manual(values = c(F = "#e41a1c", M = "#377eb8")) +
      ggplot2::labs(
        title = title,
        subtitle = depot_label,
        x = NULL,
        y = "Mean z-score across genes"
      ) +
      ggplot2::theme(
        legend.position = "none",
        axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5)
      )
  }

  p1 <- make_panel("score__lipolysis_core", "lipolysis_core (VST z-score activity)")
  p2 <- make_panel("score__thermogenesis_program", "thermogenesis_program (VST z-score activity)")
  fig <- p1 / p2

  out_png <- file.path(fig_dir, paste0("figure_candidate__lipolysis_thermogenesis__F_M_interleaved_with_medians", suffix, ".png"))
  ggplot2::ggsave(out_png, fig, width = 11.0, height = 8.2, dpi = 250)
  out_png
}

p_sub <- file.path(tab_dir, "kenichi_panel_scores_per_sample_by_sex.tsv")
p_vis <- file.path(tab_dir, "kenichi_panel_scores_per_sample_by_sex_visceral.tsv")

cat("Writing interleaved-sex figures...\n")
cat("  - ", plot_one(p_sub, "", "Adipose - Subcutaneous"), "\n", sep = "")
cat("  - ", plot_one(p_vis, "_visceral", "Adipose - Visceral (Omentum)"), "\n", sep = "")
