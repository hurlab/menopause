################################################################################
## Make simple table graphics for the enhanced PPTX (2026-02-01)
##
## Outputs:
## - GTEx_v10_AT_analysis_out/figs/ppt_tables_2026-02-01/
################################################################################

# Ensure a default CRAN mirror is set for non-interactive runs
.repos <- getOption("repos")
.cran <- unname(.repos["CRAN"])
if (is.null(.cran) || is.na(.cran) || .cran == "@CRAN@" || .cran == "") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
rm(.repos, .cran)

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readr, dplyr, tidyr, stringr, tibble, ggplot2, glue)

source("R/utils.R")

out_dir <- "GTEx_v10_AT_analysis_out"
fig_dir <- file.path(out_dir, "figs", "ppt_tables_2026-02-01")
ensure_dirs(fig_dir)

table_plot <- function(df, title, out_path, font_size = 4.2) {
  stopifnot(all(c("row", "col", "value") %in% names(df)))
  df <- df %>%
    dplyr::mutate(
      row = factor(.data$row, levels = rev(unique(.data$row))),
      col = factor(.data$col, levels = unique(.data$col))
    )
  p <- ggplot2::ggplot(df, ggplot2::aes(x = col, y = row, label = value)) +
    ggplot2::geom_tile(fill = "white", color = "grey80") +
    ggplot2::geom_text(size = font_size, family = "sans") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0, size = 14),
      axis.text.x = ggplot2::element_text(size = 10, face = "bold", margin = ggplot2::margin(t = 4)),
      axis.text.y = ggplot2::element_text(size = 10, face = "bold", margin = ggplot2::margin(r = 4)),
      axis.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    ) +
    ggplot2::labs(title = title)
  ggplot2::ggsave(out_path, p, width = 10.5, height = 5.5, dpi = 200)
}

## -----------------------------------------------------------------------------
## Table 1: subcutaneous Kenichi set effects (post vs pre) with composition adj
## -----------------------------------------------------------------------------
ken_path <- file.path(
  out_dir,
  "tables",
  "next_steps_2026-02-01_composition",
  "kenichi_gene_set_effects_with_composition_adjustment.tsv"
)
ken <- readr::read_tsv(ken_path, show_col_types = FALSE) %>%
  dplyr::filter(contrast == "post_vs_pre") %>%
  dplyr::filter(model %in% c("base", "adj_comp", "adj_comp_tech")) %>%
  dplyr::mutate(
    est = sprintf("%.3f", estimate),
    ptxt = ifelse(p < 1e-3, format(p, scientific = TRUE, digits = 2), sprintf("%.3f", p)),
    val = paste0(est, "\n(p=", ptxt, ")")
  ) %>%
  dplyr::select(set_name, model, val) %>%
  tidyr::pivot_wider(names_from = model, values_from = val) %>%
  dplyr::arrange(set_name)

ken_df <- ken %>%
  tidyr::pivot_longer(cols = c("base", "adj_comp", "adj_comp_tech"), names_to = "col", values_to = "value") %>%
  dplyr::mutate(
    row = set_name,
    col = dplyr::recode(col, base = "Base\n(SMCENTER)", adj_comp = "+Composition", adj_comp_tech = "+Comp+Tech")
  ) %>%
  dplyr::select(row, col, value)

table_plot(
  ken_df,
  title = "Subcutaneous: Kenichi gene-set effect (post vs pre)\nBase vs composition-adjusted models",
  out_path = file.path(fig_dir, "table_subq_kenichi_sets_composition_adjustment.png"),
  font_size = 4.0
)

## -----------------------------------------------------------------------------
## Table 2: senMayoLike v2 gene list
## -----------------------------------------------------------------------------
sig2 <- readr::read_tsv("references/menopause_signature_senMayoLike_v2.tsv", show_col_types = FALSE) %>%
  dplyr::mutate(
    log2FC = sprintf("%.2f", log2FoldChange),
    padj_txt = ifelse(padj < 1e-3, format(padj, scientific = TRUE, digits = 2), sprintf("%.3f", padj))
  ) %>%
  dplyr::select(gene_symbol, direction_in_post, log2FC, padj_txt) %>%
  dplyr::arrange(direction_in_post, gene_symbol)

sig2_grid <- sig2 %>%
  dplyr::mutate(row = gene_symbol) %>%
  tidyr::pivot_longer(cols = c("direction_in_post", "log2FC", "padj_txt"), names_to = "col", values_to = "value") %>%
  dplyr::mutate(
    col = dplyr::recode(
      col,
      direction_in_post = "Direction\n(in post)",
      log2FC = "Approach1\nlog2FC",
      padj_txt = "Approach1\npadj"
    )
  ) %>%
  dplyr::select(row, col, value)

table_plot(
  sig2_grid,
  title = "senMayoLike v2 menopause signature (curated)\nGenes and direction from Approach 1",
  out_path = file.path(fig_dir, "table_senMayoLike_v2_gene_list.png"),
  font_size = 3.8
)

cat("Wrote PPT table images to: ", fig_dir, "\n", sep = "")

