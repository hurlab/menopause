################################################################################
## Next steps (2026-02-01): Expanded Adrenergic / cAMP / Desensitization Modules
##
## Motivation:
## - Kenichi's caveat: baseline GTEx adipose may not show acute SNS markers.
## - Chronic SNS stimulation can produce compensatory feedback/desensitization
##   footprints (e.g., PDEs, DUSPs, RGS genes) and catecholamine resistance.
##
## Here we score small, interpretable gene modules intended to reflect:
## - cAMP/PKA/CREB and MAPK feedback regulators
## - receptor desensitization machinery
## - anti-lipolytic α2 axis / PDE3B axis
##
## Outputs (dated subfolders):
## - GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_adrenergic_modules/
## - GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_adrenergic_modules/
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
  ggplot2, patchwork,
  DESeq2
)

source("R/utils.R")

## -----------------------------------------------------------------------------
## Paths + output folders
## -----------------------------------------------------------------------------
in_dir <- "GTExDatav10"
gct_sc <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
gct_vis <- file.path(in_dir, "gene_reads_v10_adipose_visceral_omentum.gct.gz")
sattr <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

out_dir <- "GTEx_v10_AT_analysis_out"
tab_dir <- file.path(out_dir, "tables", "next_steps_2026-02-01_adrenergic_modules")
fig_dir <- file.path(out_dir, "figs", "next_steps_2026-02-01_adrenergic_modules")
log_dir <- file.path(out_dir, "logs")
ensure_dirs(c(tab_dir, fig_dir, log_dir))

log_path <- file.path(log_dir, "next_steps_2026-02-01_adrenergic_modules.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

depots <- list(
  list(
    key = "subcutaneous",
    label = "Adipose - Subcutaneous",
    gct_path = gct_sc,
    suffix = ""
  ),
  list(
    key = "visceral",
    label = "Adipose - Visceral (Omentum)",
    gct_path = gct_vis,
    suffix = "_visceral"
  )
)

## -----------------------------------------------------------------------------
## Gene modules (compact, interpretable; can be revised as needed)
## -----------------------------------------------------------------------------
## NOTE: these are intended as *proxies* in baseline tissue and are not definitive
## markers of SNS tone.
gene_modules <- list(
  camp_pka_feedback = c("PDE4D", "PDE4B", "PDE3B", "RGS2", "DUSP1", "DUSP2", "ATF3"),
  receptor_desensitization = c("ADRBK1", "ADRBK2", "ARRB1", "ARRB2", "RGS2"),
  alpha2_antilipolytic_axis = c("ADRA2A", "ADRA2B", "ADRA2C", "PDE3B"),
  camp_core = c("ADCY3", "ADCY6", "GNAS", "PRKACA", "PRKACB", "PRKAR1A", "PRKAR2B"),
  immediate_early = c("NR4A1", "NR4A2", "NR4A3", "FOS", "JUN", "EGR1", "DUSP1", "ATF3")
)
gene_modules <- purrr::map(gene_modules, ~ unique(toupper(.x)))

## Composition proxies to include in an adjusted model
composition_markers <- list(
  adipocyte = c("ADIPOQ", "PLIN1", "FABP4"),
  macrophage_mono = c("LST1", "C1QC", "TYROBP"),
  fibroblast_ecm = c("COL1A1", "COL3A1", "DCN")
)
composition_markers <- purrr::map(composition_markers, ~ unique(toupper(.x)))

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------
score_set_z <- function(expr_sym, genes) {
  genes <- intersect(genes, rownames(expr_sym))
  if (length(genes) < 3) {
    return(rep(NA_real_, ncol(expr_sym)))
  }
  m <- expr_sym[genes, , drop = FALSE]
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  z <- sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
  colMeans(z, na.rm = TRUE)
}

lm_contrast <- function(fit, term_a, term_b) {
  cf <- stats::coef(fit)
  vc <- stats::vcov(fit)
  cn <- names(cf)
  v <- rep(0, length(cf))
  names(v) <- cn
  if (!is.null(term_a) && nzchar(term_a) && term_a %in% cn) v[term_a] <- v[term_a] + 1
  if (!is.null(term_b) && nzchar(term_b) && term_b %in% cn) v[term_b] <- v[term_b] - 1
  est <- sum(v * cf)
  se <- sqrt(as.numeric(t(v) %*% vc %*% v))
  tval <- est / se
  df <- fit$df.residual
  p <- 2 * stats::pt(abs(tval), df = df, lower.tail = FALSE)
  tibble::tibble(estimate = est, se = se, t = tval, p = p)
}

## -----------------------------------------------------------------------------
## Load metadata
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("ADRIVERGIC MODULE SCORING (GTEx v10)\n")
cat("===========================================\n\n")

sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything()) %>%
  dplyr::mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])
subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything())
meta_full <- dplyr::left_join(sample_attr, subject_pheno, by = "SUBJID")

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
    dplyr::filter(
      SAMPID %in% common_samples,
      SMTS == "Adipose Tissue",
      SMTSD == depot_label
    ) %>%
    dplyr::mutate(
      sex = dplyr::case_when(SEX == 1 ~ "Male", SEX == 2 ~ "Female", TRUE ~ NA_character_),
      age_bin_label = AGE
    ) %>%
    dplyr::filter(!is.na(sex), !is.na(age_bin_label), !is.na(SMCENTER)) %>%
    dplyr::bind_cols(parse_age_bin(.$AGE)) %>%
    derive_age_numeric() %>%
    dplyr::mutate(SMCENTER = factor(SMCENTER))

  meta_f <- meta %>%
    dplyr::filter(sex == "Female") %>%
    dplyr::mutate(
      menopause_group_3 = dplyr::case_when(
        is.na(age_mid) ~ NA_character_,
        age_mid < 45 ~ "pre",
        age_mid <= 55 ~ "peri",
        age_mid > 55 ~ "post",
        TRUE ~ NA_character_
      ),
      menopause_group_3 = factor(menopause_group_3, levels = c("pre", "peri", "post"))
    ) %>%
    dplyr::filter(!is.na(menopause_group_3))

  counts <- gct$counts[, meta_f$SAMPID, drop = FALSE]
  counts <- collapse_dupe_rowsum(counts)
  counts <- round(counts)
  storage.mode(counts) <- "integer"

  cat("Computing VST...\n")
  dds0 <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = data.frame(row.names = colnames(counts)),
    design = ~ 1
  )
  vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds0, blind = TRUE))
  expr_sym <- collapse_to_symbol_mean(vst_mat, gct$gene_id, gct$gene_symbol, uppercase = TRUE)

  ## Score modules + composition proxies
  module_scores <- purrr::imap_dfc(gene_modules, ~ tibble::tibble(!!.y := score_set_z(expr_sym, .x)))
  rownames(module_scores) <- colnames(expr_sym)

  comp_scores <- purrr::imap_dfc(composition_markers, ~ tibble::tibble(!!.y := score_set_z(expr_sym, .x)))
  rownames(comp_scores) <- colnames(expr_sym)

  score_df <- meta_f %>%
    dplyr::select(SAMPID, SMCENTER, SMRIN, SMTSISCH, DTHHRDY, menopause_group_3) %>%
    dplyr::mutate(SAMPID = as.character(SAMPID)) %>%
    dplyr::bind_cols(module_scores[.$SAMPID, , drop = FALSE]) %>%
    dplyr::bind_cols(comp_scores[.$SAMPID, , drop = FALSE]) %>%
    dplyr::rename(group = menopause_group_3)

  readr::write_tsv(score_df, file.path(tab_dir, glue::glue("adrenergic_module_scores_per_sample{suffix}.tsv")))

  ## Stats: base vs composition-adjusted vs composition+tech
  stat_rows <- list()
  for (mod in names(gene_modules)) {
    df <- score_df %>% dplyr::mutate(score = .data[[mod]])
    fits <- list(
      base = stats::lm(score ~ SMCENTER + group, data = df),
      adj_comp = stats::lm(score ~ SMCENTER + group + adipocyte + macrophage_mono + fibroblast_ecm, data = df),
      adj_comp_tech = stats::lm(score ~ SMCENTER + SMRIN + SMTSISCH + group + adipocyte + macrophage_mono + fibroblast_ecm, data = df)
    )
    for (nm in names(fits)) {
      fit <- fits[[nm]]
      peri <- lm_contrast(fit, "groupperi", "") %>% dplyr::mutate(contrast = "peri_vs_pre")
      post <- lm_contrast(fit, "grouppost", "") %>% dplyr::mutate(contrast = "post_vs_pre")
      stat_rows[[length(stat_rows) + 1]] <- dplyr::bind_rows(peri, post) %>%
        dplyr::mutate(depot = depot_key, module = mod, model = nm, n = nrow(df))
    }
  }
  stats <- dplyr::bind_rows(stat_rows) %>%
    dplyr::mutate(conf_low = estimate - 1.96 * se, conf_high = estimate + 1.96 * se) %>%
    dplyr::relocate(depot, module, model, contrast, n)
  readr::write_tsv(stats, file.path(tab_dir, glue::glue("adrenergic_module_effect_stats{suffix}.tsv")))

  ## Plot: module scores by group
  plot_df <- score_df %>%
    tidyr::pivot_longer(cols = names(gene_modules), names_to = "module", values_to = "score")
  p1 <- ggplot2::ggplot(plot_df, ggplot2::aes(x = group, y = score)) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::facet_wrap(~ module, scales = "free_y", ncol = 3) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      title = paste0("Adrenergic module scores (", depot_label, ", female; baseline tissue)"),
      x = "Menopause proxy group",
      y = "Module score (mean z)"
    )

  ## Plot: effect estimates across models (post vs pre)
  p2 <- stats %>%
    dplyr::filter(contrast == "post_vs_pre") %>%
    dplyr::mutate(model = factor(model, levels = c("base", "adj_comp", "adj_comp_tech"))) %>%
    ggplot2::ggplot(ggplot2::aes(x = model, y = estimate, ymin = conf_low, ymax = conf_high)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey70") +
    ggplot2::geom_pointrange(size = 0.25) +
    ggplot2::facet_wrap(~ module, scales = "free_y", ncol = 3) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = paste0("Post vs pre effect sensitivity (", depot_label, ", female)"),
      x = "Model",
      y = "Effect estimate"
    )

  ggplot2::ggsave(
    filename = file.path(fig_dir, glue::glue("adrenergic_module_scores_by_group{suffix}.png")),
    plot = p1,
    width = 12,
    height = 6,
    dpi = 150
  )
  ggplot2::ggsave(
    filename = file.path(fig_dir, glue::glue("adrenergic_module_effect_sensitivity_post_vs_pre{suffix}.png")),
    plot = p2,
    width = 12,
    height = 5,
    dpi = 150
  )

  md <- c(
    glue::glue("# Adrenergic module scoring summary ({depot_label}; female-only)"),
    "",
    "This folder scores small adrenergic/cAMP/desensitization modules intended to better match baseline tissue and chronic SNS adaptation caveats.",
    "",
    "## Modules (v1)",
    paste0("- ", names(gene_modules), ": ", vapply(gene_modules, function(g) paste(g, collapse = ", "), character(1))),
    "",
    "## Outputs",
    glue::glue("- Per-sample module scores: `adrenergic_module_scores_per_sample{suffix}.tsv`"),
    glue::glue("- Module effect stats (base vs adjusted models): `adrenergic_module_effect_stats{suffix}.tsv`"),
    glue::glue("- Module boxplots: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_adrenergic_modules/adrenergic_module_scores_by_group{suffix}.png`"),
    glue::glue("- Effect sensitivity plot: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_adrenergic_modules/adrenergic_module_effect_sensitivity_post_vs_pre{suffix}.png`"),
    "",
    "Generated:",
    as.character(Sys.time())
  )
  writeLines(md, file.path(tab_dir, glue::glue("README_adrenergic_modules{suffix}.md")))
}

cat("\nDone. Log written to: ", log_path, "\n", sep = "")

