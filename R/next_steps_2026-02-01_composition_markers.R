################################################################################
## Next steps (2026-02-01): Cell-type Composition Proxy Analysis (Marker Scores)
##
## Motivation:
## - Bulk adipose transcriptional changes across age/menopause proxies can be
##   driven by shifts in adipocyte fraction, immune infiltration, and fibrosis.
## - We quantify simple marker-based scores and test whether adjusting for these
##   scores changes the key Kenichi gene-set results.
##
## Outputs (dated subfolders):
## - GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_composition/
## - GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/
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

## -----------------------------------------------------------------------------
## Paths + output folders
## -----------------------------------------------------------------------------
in_dir <- "GTExDatav10"
gct_sc <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
gct_vis <- file.path(in_dir, "gene_reads_v10_adipose_visceral_omentum.gct.gz")
sattr <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

out_dir <- "GTEx_v10_AT_analysis_out"
tab_dir <- file.path(out_dir, "tables", "next_steps_2026-02-01_composition")
fig_dir <- file.path(out_dir, "figs", "next_steps_2026-02-01_composition")
log_dir <- file.path(out_dir, "logs")
ensure_dirs(c(tab_dir, fig_dir, log_dir))

log_path <- file.path(log_dir, "next_steps_2026-02-01_composition_markers.log")
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

## Kenichi gene sets (reuse; we will test adjustment sensitivity)
kenichi_gene_sets <- list(
  acute_beta_adrenergic = c("NR4A1", "NR4A2", "NR4A3", "FOS", "JUN", "EGR1"),
  adrenergic_receptors = c("ADRB1", "ADRB2", "ADRB3", "ADRA2A", "ADRA2B", "ADRA2C"),
  lipolysis_core = c("PNPLA2", "LIPE", "MGLL", "ABHD5", "PLIN1", "G0S2", "LPL", "FABP4"),
  thermogenesis_program = c("UCP1", "PPARGC1A", "PPARGC1B", "CIDEA", "DIO2", "VEGFA")
)
kenichi_gene_sets <- purrr::map(kenichi_gene_sets, ~ unique(toupper(.x)))

## Marker panels (intentionally small, robust, interpretable)
marker_sets <- list(
  adipocyte = c("ADIPOQ", "PLIN1", "FABP4"),
  macrophage_mono = c("LST1", "C1QC", "TYROBP"),
  endothelial = c("PECAM1", "VWF", "KDR"),
  fibroblast_ecm = c("COL1A1", "COL3A1", "DCN"),
  t_cell = c("TRAC", "CD3D", "CD3E")
)
marker_sets <- purrr::map(marker_sets, ~ unique(toupper(.x)))

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------
score_set_z <- function(expr_sym, genes) {
  genes <- intersect(genes, rownames(expr_sym))
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
cat("COMPOSITION MARKER SCORES (GTEx v10)\n")
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
      menopause_group_3 = factor(menopause_group_3, levels = c("pre", "peri", "post")),
      menopause_group_strict2 = dplyr::case_when(
        is.na(age_mid) ~ NA_character_,
        age_mid < 45 ~ "pre",
        age_mid > 55 ~ "post",
        TRUE ~ NA_character_
      ),
      menopause_group_strict2 = factor(menopause_group_strict2, levels = c("pre", "post"))
    )

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

  ## Marker scores
  marker_scores <- purrr::imap_dfc(marker_sets, ~ tibble::tibble(!!.y := score_set_z(expr_sym, .x)))
  rownames(marker_scores) <- colnames(expr_sym)

  meta_s <- meta_f %>%
    dplyr::select(SAMPID, SMCENTER, SMRIN, SMTSISCH, DTHHRDY, menopause_group_3, menopause_group_strict2) %>%
    dplyr::mutate(SAMPID = as.character(SAMPID))
  ms <- dplyr::bind_cols(meta_s, marker_scores[meta_s$SAMPID, , drop = FALSE])
  readr::write_tsv(ms, file.path(tab_dir, glue::glue("marker_scores_per_sample{suffix}.tsv")))

  ## Marker score stats (3-group; SMCENTER-adjusted only + sensitivity with RIN/ischemic)
  stat_rows <- list()
  for (mname in names(marker_sets)) {
    df <- ms %>%
      dplyr::filter(!is.na(menopause_group_3)) %>%
      dplyr::mutate(group = menopause_group_3, score = .data[[mname]])
    fits <- list(
      base = stats::lm(score ~ SMCENTER + group, data = df),
      rin_isch = stats::lm(score ~ SMCENTER + SMRIN + SMTSISCH + group, data = df)
    )
    for (nm in names(fits)) {
      fit <- fits[[nm]]
      peri <- lm_contrast(fit, "groupperi", "")
      post <- lm_contrast(fit, "grouppost", "")
      stat_rows[[length(stat_rows) + 1]] <- dplyr::bind_rows(
        peri %>% dplyr::mutate(contrast = "peri_vs_pre"),
        post %>% dplyr::mutate(contrast = "post_vs_pre")
      ) %>%
        dplyr::mutate(depot = depot_key, marker = mname, model = nm, n = nrow(df))
    }
  }
  marker_stats <- dplyr::bind_rows(stat_rows) %>%
    dplyr::mutate(conf_low = estimate - 1.96 * se, conf_high = estimate + 1.96 * se) %>%
    dplyr::relocate(depot, marker, model, contrast, n)
  readr::write_tsv(marker_stats, file.path(tab_dir, glue::glue("marker_score_stats{suffix}.tsv")))

  ## Plot: marker scores across menopause proxy groups
  plot_df <- ms %>%
    dplyr::filter(!is.na(menopause_group_3)) %>%
    tidyr::pivot_longer(cols = names(marker_sets), names_to = "marker", values_to = "score")
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = menopause_group_3, y = score)) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::facet_wrap(~ marker, scales = "free_y") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      title = paste0("Composition proxy marker scores (", depot_label, ", female)"),
      x = "Menopause proxy group",
      y = "Marker score (mean z)"
    )
  ggplot2::ggsave(
    filename = file.path(fig_dir, glue::glue("marker_scores_by_group{suffix}.png")),
    plot = p,
    width = 10,
    height = 5,
    dpi = 150
  )

  ## ---------------------------------------------------------------------------
  ## Does adjusting for composition proxies change the Kenichi gene-set effects?
  ## ---------------------------------------------------------------------------
  kenichi_scores <- purrr::imap_dfc(kenichi_gene_sets, ~ tibble::tibble(!!.y := score_set_z(expr_sym, .x)))
  rownames(kenichi_scores) <- colnames(expr_sym)
  df0 <- dplyr::bind_cols(meta_s, kenichi_scores[meta_s$SAMPID, , drop = FALSE], marker_scores[meta_s$SAMPID, , drop = FALSE]) %>%
    dplyr::filter(!is.na(menopause_group_3)) %>%
    dplyr::mutate(group = menopause_group_3)

  adjust_rows <- list()
  for (set_name in names(kenichi_gene_sets)) {
    df <- df0 %>% dplyr::mutate(score = .data[[set_name]])

    fit_base <- stats::lm(score ~ SMCENTER + group, data = df)
    fit_adj <- stats::lm(score ~ SMCENTER + group + adipocyte + macrophage_mono + fibroblast_ecm, data = df)
    # optional: add RIN/ischemic too
    fit_adj_tech <- stats::lm(score ~ SMCENTER + SMRIN + SMTSISCH + group + adipocyte + macrophage_mono + fibroblast_ecm, data = df)

    for (nm in c("base", "adj_comp", "adj_comp_tech")) {
      fit <- switch(nm, base = fit_base, adj_comp = fit_adj, adj_comp_tech = fit_adj_tech)
      peri <- lm_contrast(fit, "groupperi", "") %>% dplyr::mutate(contrast = "peri_vs_pre")
      post <- lm_contrast(fit, "grouppost", "") %>% dplyr::mutate(contrast = "post_vs_pre")
      adjust_rows[[length(adjust_rows) + 1]] <- dplyr::bind_rows(peri, post) %>%
        dplyr::mutate(depot = depot_key, set_name = set_name, model = nm, n = nrow(df))
    }
  }
  adj_stats <- dplyr::bind_rows(adjust_rows) %>%
    dplyr::mutate(conf_low = estimate - 1.96 * se, conf_high = estimate + 1.96 * se) %>%
    dplyr::relocate(depot, set_name, model, contrast, n)
  readr::write_tsv(adj_stats, file.path(tab_dir, glue::glue("kenichi_gene_set_effects_with_composition_adjustment{suffix}.tsv")))

  p2 <- adj_stats %>%
    dplyr::filter(contrast %in% c("peri_vs_pre", "post_vs_pre")) %>%
    dplyr::mutate(model = factor(model, levels = c("base", "adj_comp", "adj_comp_tech"))) %>%
    ggplot2::ggplot(ggplot2::aes(x = model, y = estimate, ymin = conf_low, ymax = conf_high)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey70") +
    ggplot2::geom_pointrange(size = 0.25) +
    ggplot2::facet_grid(contrast ~ set_name, scales = "free_y") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = paste0("Kenichi gene-set effects with composition adjustment (", depot_label, ", female)"),
      x = "Model",
      y = "Effect estimate (vs pre)"
    )
  ggplot2::ggsave(
    filename = file.path(fig_dir, glue::glue("kenichi_gene_set_effects_with_composition_adjustment{suffix}.png")),
    plot = p2,
    width = 12,
    height = 5,
    dpi = 150
  )

  ## README
  md <- c(
    glue::glue("# Composition proxy summary ({depot_label}; female-only)"),
    "",
    "This folder quantifies simple marker-based composition proxy scores and tests whether adjusting for these proxies changes the key Kenichi gene-set effects.",
    "",
    "## Marker panels",
    "- adipocyte: ADIPOQ, PLIN1, FABP4",
    "- macrophage_mono: LST1, C1QC, TYROBP",
    "- endothelial: PECAM1, VWF, KDR",
    "- fibroblast_ecm: COL1A1, COL3A1, DCN",
    "- t_cell: TRAC, CD3D, CD3E",
    "",
    "## Outputs",
    glue::glue("- Per-sample marker scores: `marker_scores_per_sample{suffix}.tsv`"),
    glue::glue("- Marker score group stats: `marker_score_stats{suffix}.tsv`"),
    glue::glue("- Kenichi gene-set effects (base vs composition-adjusted): `kenichi_gene_set_effects_with_composition_adjustment{suffix}.tsv`"),
    glue::glue("- Marker score plot: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/marker_scores_by_group{suffix}.png`"),
    glue::glue("- Kenichi effect sensitivity plot: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/kenichi_gene_set_effects_with_composition_adjustment{suffix}.png`"),
    "",
    "Generated:",
    as.character(Sys.time())
  )
  writeLines(md, file.path(tab_dir, glue::glue("README_composition_markers{suffix}.md")))
}

cat("\nDone. Log written to: ", log_path, "\n", sep = "")

