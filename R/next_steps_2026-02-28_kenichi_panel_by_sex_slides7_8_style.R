################################################################################
## Next steps (2026-02-28): Slides 7-8 style Kenichi panel, but with Male+Female
##
## Slide context:
## - Slides 7-8 in docs/menopause_sns_gtex_enhanced_2026-02-02_v4.1.2.pdf show
##   female-only "Kenichi panel" baseline gene-set activity (lipolysis + thermo)
##   within each depot, grouped into pre/peri/post (age-bin proxy).
##
## Here:
## - We reproduce the same plot style, but facet by sex (Male vs Female).
## - Scores are computed as in the deck: VST, per-gene z across samples within
##   the plotted cohort (depot+sex), then mean(z) across genes in the set.
## - Stats: score ~ SMCENTER + group, with contrasts peri-pre, post-pre, post-peri.
##
## Outputs:
## - GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_kenichi_panel_by_sex/
## - GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_kenichi_panel_by_sex/
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

in_dir <- "GTExDatav10"
gct_sc <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
gct_vis <- file.path(in_dir, "gene_reads_v10_adipose_visceral_omentum.gct.gz")
sattr <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

out_dir <- "GTEx_v10_AT_analysis_out"
tab_dir <- file.path(out_dir, "tables", "next_steps_2026-02-28_kenichi_panel_by_sex")
fig_dir <- file.path(out_dir, "figs", "next_steps_2026-02-28_kenichi_panel_by_sex")
log_dir <- file.path(out_dir, "logs")
ensure_dirs(c(tab_dir, fig_dir, log_dir))

log_path <- file.path(log_dir, "next_steps_2026-02-28_kenichi_panel_by_sex.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

depots <- list(
  list(key = "subcutaneous", label = "Adipose - Subcutaneous", gct_path = gct_sc, suffix = ""),
  list(key = "visceral", label = "Adipose - Visceral (Omentum)", gct_path = gct_vis, suffix = "_visceral")
)

gene_sets <- list(
  lipolysis_core = c("PNPLA2", "LIPE", "MGLL", "ABHD5", "PLIN1", "G0S2", "LPL", "FABP4"),
  thermogenesis_program = c("UCP1", "PPARGC1A", "PPARGC1B", "CIDEA", "DIO2", "VEGFA")
)
gene_sets <- purrr::map(gene_sets, ~ unique(toupper(.x)))

score_set_z <- function(expr_sym, genes) {
  genes <- intersect(genes, rownames(expr_sym))
  if (length(genes) < 3) return(rep(NA_real_, ncol(expr_sym)))
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

cat("Loading metadata...\n")
sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything()) %>%
  dplyr::mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])
subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything())
meta_full <- dplyr::left_join(sample_attr, subject_pheno, by = "SUBJID")

stats_rows <- list()

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
    dplyr::bind_cols(parse_age_bin(.$AGE)) %>%
    derive_age_numeric() %>%
    dplyr::mutate(
      sex = factor(sex, levels = c("Male", "Female")),
      SMCENTER = factor(SMCENTER),
      menopause_group_3 = dplyr::case_when(
        is.na(age_mid) ~ NA_character_,
        age_mid < 45 ~ "pre",         # maps to 20-49
        age_mid <= 55 ~ "peri",       # maps to 50-59
        age_mid > 55 ~ "post",        # maps to 60+
        TRUE ~ NA_character_
      ),
      menopause_group_3 = factor(menopause_group_3, levels = c("pre", "peri", "post"))
    ) %>%
    dplyr::filter(!is.na(menopause_group_3))

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
  expr_sym_all <- collapse_to_symbol_mean(vst_mat, gct$gene_id, gct$gene_symbol, uppercase = TRUE)

  score_rows <- list()
  for (sx in levels(meta$sex)) {
    meta_sx <- meta %>% dplyr::filter(sex == sx)
    expr_sym <- expr_sym_all[, meta_sx$SAMPID, drop = FALSE]
    score_tbl <- tibble::tibble(SAMPID = colnames(expr_sym))
    for (nm in names(gene_sets)) {
      score_tbl[[paste0("score__", nm)]] <- score_set_z(expr_sym, gene_sets[[nm]])
    }
    score_tbl$sex <- sx
    score_rows[[length(score_rows) + 1]] <- score_tbl
  }
  score_tbl_all <- dplyr::bind_rows(score_rows)

  df0 <- meta %>%
    dplyr::left_join(score_tbl_all, by = c("SAMPID", "sex")) %>%
    dplyr::mutate(sex = factor(as.character(sex), levels = c("Male", "Female")))

  # Save per-sample table
  out_tsv <- file.path(tab_dir, paste0("kenichi_panel_scores_per_sample_by_sex", suffix, ".tsv"))
  readr::write_tsv(df0 %>%
                     dplyr::select(SAMPID, depot = SMTSD, sex, AGE, age_mid, menopause_group_3, SMCENTER,
                                   dplyr::starts_with("score__")),
                   out_tsv)

  # Fit models + plots
  plot_panels <- list()
  for (nm in names(gene_sets)) {
    col <- paste0("score__", nm)

    stat_by_sex <- list()
    for (sx in levels(df0$sex)) {
      df <- df0 %>%
        dplyr::filter(sex == sx, !is.na(.data[[col]]), !is.na(menopause_group_3))
      fit <- stats::lm(stats::as.formula(paste0(col, " ~ SMCENTER + menopause_group_3")), data = df)
      term_peri <- "menopause_group_3peri"
      term_post <- "menopause_group_3post"
      peri_vs_pre <- lm_contrast(fit, term_peri, "")
      post_vs_pre <- lm_contrast(fit, term_post, "")
      post_vs_peri <- lm_contrast(fit, term_post, term_peri)
      stat_by_sex[[sx]] <- tibble::tibble(
        depot = depot_key,
        sex = sx,
        gene_set = nm,
        n = nrow(df),
        peri_vs_pre_p = peri_vs_pre$p,
        post_vs_pre_p = post_vs_pre$p,
        post_vs_peri_p = post_vs_peri$p
      )
    }
    stats_rows[[length(stats_rows) + 1]] <- dplyr::bind_rows(stat_by_sex)

    df_plot <- df0 %>%
      dplyr::filter(!is.na(.data[[col]]), !is.na(menopause_group_3)) %>%
      dplyr::mutate(sex = factor(as.character(sex), levels = c("Male", "Female")))

    # Add p-values per facet from stats computed above
    stat_df <- dplyr::bind_rows(stat_by_sex) %>%
      dplyr::mutate(
        subtitle = paste0(
          "LM adj (SMCENTER): peri vs pre p=", signif(peri_vs_pre_p, 3),
          "; post vs pre p=", signif(post_vs_pre_p, 3),
          "; post vs peri p=", signif(post_vs_peri_p, 3)
        )
      )

    # Build panel plot and facet by sex
    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = menopause_group_3, y = .data[[col]], color = menopause_group_3)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.25) +
      ggplot2::geom_jitter(width = 0.18, alpha = 0.35, size = 0.9) +
      ggplot2::facet_wrap(~ sex, nrow = 1) +
      ggplot2::theme_bw() +
      ggplot2::scale_color_manual(values = c(pre = "#1b9e77", peri = "#d95f02", post = "#7570b3")) +
      ggplot2::labs(
        title = glue::glue("{nm} (VST z-score activity)"),
        x = NULL,
        y = "Mean z-score across genes"
      ) +
      ggplot2::theme(legend.position = "none")

    # Add per-sex subtitle text into each facet
    p <- p + ggplot2::geom_text(
      data = stat_df,
      ggplot2::aes(x = "peri", y = Inf, label = subtitle),
      inherit.aes = FALSE,
      vjust = 1.2,
      size = 2.7
    )

    plot_panels[[nm]] <- p
    ggplot2::ggsave(file.path(fig_dir, glue::glue("by_sex__gene_set_activity__{nm}{suffix}.png")), p,
                    width = 11.0, height = 4.2, dpi = 250)
  }

  if (all(c("lipolysis_core", "thermogenesis_program") %in% names(plot_panels))) {
    fig_main <- plot_panels[["lipolysis_core"]] / plot_panels[["thermogenesis_program"]]
    ggplot2::ggsave(
      filename = file.path(fig_dir, paste0("figure_candidate__lipolysis_thermogenesis__by_sex", suffix, ".png")),
      plot = fig_main,
      width = 11.0,
      height = 8.2,
      dpi = 250
    )
  }
}

stats_tbl <- dplyr::bind_rows(stats_rows) %>%
  dplyr::arrange(depot, gene_set, sex)
readr::write_tsv(stats_tbl, file.path(tab_dir, "kenichi_panel_lm_pvalues_by_sex.tsv"))

cat("\nWrote:\n")
cat("  - ", file.path(tab_dir, "kenichi_panel_lm_pvalues_by_sex.tsv"), "\n", sep = "")
cat("  - Figures under: ", fig_dir, "\n", sep = "")
cat("\nLog: ", log_path, "\n", sep = "")
