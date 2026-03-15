################################################################################
## Next steps (2026-02-28): Adrenergic receptor gene expression heatmaps
##
## Goal:
## - Visualize expression levels of key adrenergic receptor genes as heatmaps
##   with depot / sex / age-bin context.
##
## Output:
## - GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_receptor_heatmaps/
## - GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_receptor_heatmaps/
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
  readr, dplyr, tidyr, stringr, tibble, purrr, glue,
  ggplot2,
  DESeq2
)

source("R/utils.R")

in_dir <- "GTExDatav10"
gct_sc <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
gct_vis <- file.path(in_dir, "gene_reads_v10_adipose_visceral_omentum.gct.gz")
sattr <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

out_dir <- "GTEx_v10_AT_analysis_out"
tab_dir <- file.path(out_dir, "tables", "next_steps_2026-02-28_receptor_heatmaps")
fig_dir <- file.path(out_dir, "figs", "next_steps_2026-02-28_receptor_heatmaps")
log_dir <- file.path(out_dir, "logs")
ensure_dirs(c(tab_dir, fig_dir, log_dir))

log_path <- file.path(log_dir, "next_steps_2026-02-28_receptor_heatmaps.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

depots <- list(
  list(key = "subcutaneous", label = "Adipose - Subcutaneous", gct_path = gct_sc, suffix = ""),
  list(key = "visceral", label = "Adipose - Visceral (Omentum)", gct_path = gct_vis, suffix = "_visceral")
)

genes_beta <- c("ADRB1", "ADRB2", "ADRB3")
genes_alpha2 <- c("ADRA2A", "ADRA2B", "ADRA2C")

age_levels <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")

cat("Loading metadata...\n")
sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything()) %>%
  dplyr::mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])
subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything())
meta_full <- dplyr::left_join(sample_attr, subject_pheno, by = "SUBJID")

all_long <- list()

for (dpt in depots) {
  depot_label <- dpt$label
  depot_key <- dpt$key

  cat("\n-------------------------------------------\n")
  cat("Depot: ", depot_label, " (", depot_key, ")\n", sep = "")
  cat("-------------------------------------------\n\n")

  gct <- read_gct_v12(dpt$gct_path)
  common_samples <- intersect(colnames(gct$counts), meta_full$SAMPID)

  meta <- meta_full %>%
    dplyr::filter(SAMPID %in% common_samples, SMTS == "Adipose Tissue", SMTSD == depot_label) %>%
    dplyr::mutate(
      sex = dplyr::case_when(SEX == 1 ~ "Male", SEX == 2 ~ "Female", TRUE ~ NA_character_),
      age_bin_label = AGE
    ) %>%
    dplyr::filter(!is.na(sex), !is.na(age_bin_label), !is.na(SMCENTER)) %>%
    dplyr::bind_cols(parse_age_bin(.$AGE)) %>%
    derive_age_numeric() %>%
    dplyr::mutate(
      sex = factor(sex, levels = c("Male", "Female")),
      age_bin_label = factor(as.character(age_bin_label), levels = age_levels),
      depot = depot_key
    ) %>%
    dplyr::filter(!is.na(age_bin_label))

  counts <- gct$counts[, meta$SAMPID, drop = FALSE]
  counts <- collapse_dupe_rowsum(counts)
  counts <- round(counts)
  storage.mode(counts) <- "integer"

  cat("Computing VST (this depot, both sexes)...\n")
  dds0 <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = data.frame(row.names = colnames(counts)),
    design = ~ 1
  )
  vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds0, blind = TRUE))
  expr_sym <- collapse_to_symbol_mean(vst_mat, gct$gene_id, gct$gene_symbol, uppercase = TRUE)

  genes <- unique(c(genes_beta, genes_alpha2))
  genes_present <- intersect(genes, rownames(expr_sym))
  if (length(genes_present) == 0) stop("None of the requested genes were found in expr_sym for depot=", depot_key)

  long <- as.data.frame(t(expr_sym[genes_present, meta$SAMPID, drop = FALSE])) %>%
    tibble::rownames_to_column("SAMPID") %>%
    tidyr::pivot_longer(cols = -SAMPID, names_to = "gene", values_to = "vst") %>%
    dplyr::left_join(meta %>% dplyr::select(SAMPID, depot, sex, age_bin_label, age_mid, age_years, SMCENTER), by = "SAMPID") %>%
    dplyr::mutate(
      gene = factor(gene, levels = genes),
      family = dplyr::case_when(
        gene %in% genes_beta ~ "beta",
        gene %in% genes_alpha2 ~ "alpha2",
        TRUE ~ "other"
      ),
      family = factor(family, levels = c("beta", "alpha2"))
    )

  all_long[[length(all_long) + 1]] <- long
}

expr_long <- dplyr::bind_rows(all_long) %>%
  dplyr::filter(!is.na(vst))

summary_tbl <- expr_long %>%
  dplyr::group_by(depot, sex, age_bin_label, family, gene) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean_vst = mean(vst, na.rm = TRUE),
    sd_vst = sd(vst, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(summary_tbl, file.path(tab_dir, "adrenergic_receptor_expression_group_means_vst.tsv"))

plot_heatmap <- function(df, family_label, out_png) {
  df <- df %>% dplyr::filter(family == family_label)
  if (!nrow(df)) stop("No rows for family=", family_label)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = age_bin_label, y = gene, fill = mean_vst)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::facet_grid(sex ~ depot, scales = "free_x") +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = median(df$mean_vst, na.rm = TRUE),
      name = "Mean VST"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = glue::glue("Adrenergic receptor expression (VST; group means)"),
      subtitle = glue::glue("Family: {family_label}; facets: sex x depot; x-axis: GTEx AGE bins"),
      x = "GTEx age bin",
      y = "Gene"
    )
  ggplot2::ggsave(out_png, p, width = 11.0, height = 4.8, dpi = 200)
  invisible(p)
}

plot_heatmap_with_values <- function(df, family_label, out_png) {
  df <- df %>% dplyr::filter(family == family_label)
  if (!nrow(df)) stop("No rows for family=", family_label)
  df <- df %>% dplyr::mutate(label = formatC(mean_vst, format = "f", digits = 2))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = age_bin_label, y = gene, fill = mean_vst)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 2.6) +
    ggplot2::facet_grid(sex ~ depot, scales = "free_x") +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = median(df$mean_vst, na.rm = TRUE),
      name = "Mean VST"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = glue::glue("Adrenergic receptor expression (VST; group means; values shown)"),
      subtitle = glue::glue("Family: {family_label}; facets: sex x depot; x-axis: GTEx AGE bins; cells: mean VST (2 decimals)"),
      x = "GTEx age bin",
      y = "Gene"
    )
  ggplot2::ggsave(out_png, p, width = 11.0, height = 4.8, dpi = 250)
  invisible(p)
}

plot_per_sample_heatmap <- function(df_long, family_label, out_png) {
  df <- df_long %>%
    dplyr::filter(family == family_label) %>%
    dplyr::filter(!is.na(vst), !is.na(age_bin_label), !is.na(sex), !is.na(depot))
  if (!nrow(df)) stop("No rows for family=", family_label)

  df <- df %>%
    dplyr::group_by(depot, sex) %>%
    dplyr::arrange(age_mid, age_years, age_bin_label, SAMPID, .by_group = TRUE) %>%
    dplyr::mutate(sample_order = dplyr::row_number()) %>%
    dplyr::ungroup()

  # Use a numeric y scale so we can place age-bin labels above the heatmap.
  gene_levels <- levels(df$gene)
  df <- df %>%
    dplyr::mutate(
      gene_i = length(gene_levels) - as.numeric(factor(gene, levels = gene_levels)) + 1
    )

  # Bin boundary vertical lines (within each facet).
  bounds <- df %>%
    dplyr::group_by(depot, sex, age_bin_label) %>%
    dplyr::summarise(max_order = max(sample_order), .groups = "drop") %>%
    dplyr::group_by(depot, sex) %>%
    dplyr::arrange(max_order, .by_group = TRUE) %>%
    dplyr::mutate(xint = max_order + 0.5) %>%
    dplyr::mutate(is_last = dplyr::row_number() == dplyr::n()) %>%
    dplyr::filter(!is_last) %>%  # don't draw the last line
    dplyr::select(-is_last) %>%
    dplyr::ungroup()

  bins <- df %>%
    dplyr::group_by(depot, sex, age_bin_label) %>%
    dplyr::summarise(
      xmin = min(sample_order),
      xmax = max(sample_order),
      xmid = (xmin + xmax) / 2,
      .groups = "drop"
    ) %>%
    dplyr::mutate(y = max(df$gene_i) + 0.7)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = sample_order, y = gene_i, fill = vst)) +
    ggplot2::geom_tile() +
    ggplot2::geom_vline(data = bounds, ggplot2::aes(xintercept = xint), inherit.aes = FALSE,
                        color = "black", alpha = 0.15, linewidth = 0.2) +
    ggplot2::geom_text(
      data = bins,
      ggplot2::aes(x = xmid, y = y, label = as.character(age_bin_label)),
      inherit.aes = FALSE,
      size = 2.3
    ) +
    ggplot2::facet_grid(sex ~ depot, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_gradient2(
      low = "#2166ac",
      mid = "white",
      high = "#b2182b",
      midpoint = median(df$vst, na.rm = TRUE),
      name = "VST"
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_len(length(gene_levels)),
      labels = rev(gene_levels),
      expand = ggplot2::expansion(mult = c(0.02, 0.20))
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = glue::glue("Adrenergic receptor expression (VST; per-sample heatmap)"),
      subtitle = glue::glue("Family: {family_label}; columns are individual samples ordered by age bin (labels shown above heatmap) within each (sex x depot) facet"),
      x = "Samples (ordered by age bin; labels omitted)",
      y = "Gene"
    )
  ggplot2::ggsave(out_png, p, width = 12.0, height = 5.2, dpi = 250)
  invisible(p)
}

plot_heatmap(summary_tbl, "beta", file.path(fig_dir, "heatmap_ADRB1_ADRB2_ADRB3_meanVST_by_depot_sex_agebin.png"))
plot_heatmap(summary_tbl, "alpha2", file.path(fig_dir, "heatmap_ADRA2A_ADRA2B_ADRA2C_meanVST_by_depot_sex_agebin.png"))

plot_heatmap_with_values(summary_tbl, "beta", file.path(fig_dir, "heatmap_ADRB1_ADRB2_ADRB3_meanVST_by_depot_sex_agebin_with_values.png"))
plot_heatmap_with_values(summary_tbl, "alpha2", file.path(fig_dir, "heatmap_ADRA2A_ADRA2B_ADRA2C_meanVST_by_depot_sex_agebin_with_values.png"))

plot_per_sample_heatmap(expr_long, "beta", file.path(fig_dir, "heatmap_ADRB1_ADRB2_ADRB3_VST_per_sample_ordered_by_agebin.png"))
plot_per_sample_heatmap(expr_long, "alpha2", file.path(fig_dir, "heatmap_ADRA2A_ADRA2B_ADRA2C_VST_per_sample_ordered_by_agebin.png"))

cat("\nWrote:\n")
cat("  - ", file.path(tab_dir, "adrenergic_receptor_expression_group_means_vst.tsv"), "\n", sep = "")
cat("  - ", file.path(fig_dir, "heatmap_ADRB1_ADRB2_ADRB3_meanVST_by_depot_sex_agebin.png"), "\n", sep = "")
cat("  - ", file.path(fig_dir, "heatmap_ADRA2A_ADRA2B_ADRA2C_meanVST_by_depot_sex_agebin.png"), "\n", sep = "")
cat("  - ", file.path(fig_dir, "heatmap_ADRB1_ADRB2_ADRB3_meanVST_by_depot_sex_agebin_with_values.png"), "\n", sep = "")
cat("  - ", file.path(fig_dir, "heatmap_ADRA2A_ADRA2B_ADRA2C_meanVST_by_depot_sex_agebin_with_values.png"), "\n", sep = "")
cat("  - ", file.path(fig_dir, "heatmap_ADRB1_ADRB2_ADRB3_VST_per_sample_ordered_by_agebin.png"), "\n", sep = "")
cat("  - ", file.path(fig_dir, "heatmap_ADRA2A_ADRA2B_ADRA2C_VST_per_sample_ordered_by_agebin.png"), "\n", sep = "")
cat("\nLog: ", log_path, "\n", sep = "")
