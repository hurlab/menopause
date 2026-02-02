################################################################################
## Next steps (2026-02-02): Sex + age-bin trend plots for key adipose scores
##
## Motivation (Kenichi email):
## - Male vs female comparisons can provide context but are imperfect controls.
## - Create simple, interpretable plots of score trajectories across GTEx age bins
##   by sex, separately for subcutaneous and visceral adipose.
##
## Outputs:
## - GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-02_sex_interaction_scores/
## - GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-02_sex_interaction_scores/
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
tab_dir <- file.path(out_dir, "tables", "next_steps_2026-02-02_sex_interaction_scores")
fig_dir <- file.path(out_dir, "figs", "next_steps_2026-02-02_sex_interaction_scores")
log_dir <- file.path(out_dir, "logs")
ensure_dirs(c(tab_dir, fig_dir, log_dir))

log_path <- file.path(log_dir, "next_steps_2026-02-02_sex_age_trends_plots.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

depots <- list(
  list(key = "subcutaneous", label = "Adipose - Subcutaneous", gct_path = gct_sc, suffix = ""),
  list(key = "visceral", label = "Adipose - Visceral (Omentum)", gct_path = gct_vis, suffix = "_visceral")
)

kenichi_sets <- list(
  acute_beta_adrenergic = c("NR4A1", "NR4A2", "NR4A3", "FOS", "JUN", "EGR1"),
  lipolysis_core = c("PNPLA2", "LIPE", "MGLL", "ABHD5", "PLIN1", "G0S2", "LPL", "FABP4"),
  thermogenesis_program = c("UCP1", "PPARGC1A", "PPARGC1B", "CIDEA", "DIO2", "VEGFA")
)
kenichi_sets <- purrr::map(kenichi_sets, ~ unique(toupper(.x)))

composition_sets <- list(
  adipocyte = c("ADIPOQ", "PLIN1", "FABP4"),
  macrophage_mono = c("LST1", "C1QC", "TYROBP")
)
composition_sets <- purrr::map(composition_sets, ~ unique(toupper(.x)))

adrenergic_modules <- list(
  camp_core = c("ADCY3", "ADCY6", "GNAS", "PRKACA", "PRKACB", "PRKAR1A", "PRKAR2B")
)
adrenergic_modules <- purrr::map(adrenergic_modules, ~ unique(toupper(.x)))

sig_v2 <- readr::read_tsv("references/menopause_signature_senMayoLike_v2.tsv", show_col_types = FALSE) %>%
  dplyr::mutate(gene_symbol = toupper(as.character(gene_symbol)))
sig2_up <- sig_v2 %>% dplyr::filter(direction_in_post == "up") %>% dplyr::pull(gene_symbol)
sig2_down <- sig_v2 %>% dplyr::filter(direction_in_post == "down") %>% dplyr::pull(gene_symbol)

score_set_z <- function(expr_sym, genes) {
  genes <- intersect(unique(toupper(genes)), rownames(expr_sym))
  if (length(genes) < 2) {
    return(rep(NA_real_, ncol(expr_sym)))
  }
  m <- expr_sym[genes, , drop = FALSE]
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  z <- sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
  colMeans(z, na.rm = TRUE)
}

add_scores <- function(tbl, expr_sym, sets) {
  for (nm in names(sets)) {
    tbl[[nm]] <- score_set_z(expr_sym, sets[[nm]])
  }
  tbl
}

cat("Loading metadata...\n")
sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything()) %>%
  dplyr::mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])
subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything())
meta_full <- dplyr::left_join(sample_attr, subject_pheno, by = "SUBJID")

all_score_rows <- list()

for (dpt in depots) {
  depot_label <- dpt$label
  depot_key <- dpt$key
  suffix <- dpt$suffix

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
    dplyr::mutate(
      sex = factor(sex, levels = c("Male", "Female")),
      SMCENTER = factor(SMCENTER),
      age_bin_label = factor(age_bin_label, levels = sort(unique(as.character(age_bin_label))))
    )

  counts <- gct$counts[, meta$SAMPID, drop = FALSE]
  counts <- collapse_dupe_rowsum(counts)
  counts <- round(counts)
  storage.mode(counts) <- "integer"

  dds0 <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = data.frame(row.names = colnames(counts)),
    design = ~ 1
  )
  vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds0, blind = TRUE))
  expr_sym <- collapse_to_symbol_mean(vst_mat, gct$gene_id, gct$gene_symbol, uppercase = TRUE)

  score_tbl <- tibble::tibble(SAMPID = colnames(expr_sym)) %>%
    dplyr::mutate(
      senMayoLike_v2_z_signed = score_set_z(expr_sym, sig2_up) - score_set_z(expr_sym, sig2_down)
    ) %>%
    add_scores(expr_sym, kenichi_sets) %>%
    add_scores(expr_sym, composition_sets) %>%
    add_scores(expr_sym, adrenergic_modules) %>%
    dplyr::mutate(depot = depot_key)

  meta2 <- meta %>%
    dplyr::left_join(score_tbl, by = "SAMPID") %>%
    dplyr::mutate(depot = depot_key)

  # Long-form
  score_cols <- c(
    "senMayoLike_v2_z_signed",
    names(kenichi_sets),
    names(composition_sets),
    names(adrenergic_modules)
  )

  long <- meta2 %>%
    dplyr::select(SAMPID, depot, sex, age_bin_label, dplyr::all_of(score_cols)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(score_cols), names_to = "score_name", values_to = "score") %>%
    dplyr::filter(!is.na(score))

  all_score_rows[[length(all_score_rows) + 1]] <- long

  # Summary by sex x age bin for key scores (for quick tables)
  sum_tbl <- long %>%
    dplyr::group_by(depot, score_name, sex, age_bin_label) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(score, na.rm = TRUE),
      sd = sd(score, na.rm = TRUE),
      se = sd / sqrt(n),
      .groups = "drop"
    )

  readr::write_tsv(sum_tbl, file.path(tab_dir, glue::glue("scores_by_sex_agebin_summary{suffix}.tsv")))
}

scores_long <- dplyr::bind_rows(all_score_rows)
readr::write_tsv(scores_long, file.path(tab_dir, "scores_by_sex_agebin_long.tsv"))

## -----------------------------------------------------------------------------
## Plot: trajectories across age bins by sex (facet by depot and score)
## -----------------------------------------------------------------------------
plot_scores <- c(
  "lipolysis_core", "thermogenesis_program", "acute_beta_adrenergic",
  "adipocyte", "macrophage_mono",
  "camp_core",
  "senMayoLike_v2_z_signed"
)

plot_df <- scores_long %>%
  dplyr::filter(score_name %in% plot_scores) %>%
  dplyr::group_by(depot, score_name, sex, age_bin_label) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(score, na.rm = TRUE),
    se = sd(score, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  ) %>%
  dplyr::filter(n >= 5)

# order age bins in a biologically sensible way (20-29, 30-39, ..., 80+)
age_levels <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
plot_df <- plot_df %>%
  dplyr::mutate(age_bin_label = factor(as.character(age_bin_label), levels = age_levels)) %>%
  dplyr::filter(!is.na(age_bin_label))

pretty_names <- c(
  lipolysis_core = "Lipolysis core",
  thermogenesis_program = "Thermogenesis program",
  acute_beta_adrenergic = "Acute beta-adrenergic (IEG)",
  adipocyte = "Adipocyte proxy",
  macrophage_mono = "Macrophage/mono proxy",
  camp_core = "cAMP core module",
  senMayoLike_v2_z_signed = "senMayoLike v2 (signed z)"
)

plot_df <- plot_df %>%
  dplyr::mutate(score_label = dplyr::recode(score_name, !!!pretty_names))

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = age_bin_label, y = mean, color = sex, group = sex)) +
  ggplot2::geom_line(linewidth = 0.6, alpha = 0.9) +
  ggplot2::geom_point(size = 1.7) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), width = 0.15, linewidth = 0.4) +
  ggplot2::facet_grid(score_label ~ depot, scales = "free_y") +
  ggplot2::theme_bw(base_size = 9) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  ggplot2::labs(
    title = "Adipose score trajectories across GTEx age bins (Male vs Female)",
    subtitle = "Scores are VST-based gene-set mean z-scores within each depot; error bars show 95% CI of the mean",
    x = "GTEx age bin",
    y = "Score (relative; z-based)"
  )

ggplot2::ggsave(file.path(fig_dir, "sex_agebin_score_trajectories_key_scores.png"), p, width = 10.5, height = 8.5, dpi = 150)

cat("\nDone.\n")
cat("Wrote per-sample long scores: ", file.path(tab_dir, "scores_by_sex_agebin_long.tsv"), "\n", sep = "")
cat("Wrote key plot: ", file.path(fig_dir, "sex_agebin_score_trajectories_key_scores.png"), "\n", sep = "")
cat("Log: ", log_path, "\n", sep = "")

