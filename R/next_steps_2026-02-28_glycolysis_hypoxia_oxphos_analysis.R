################################################################################
## Next steps (2026-02-28): Glycolysis-focused hypotheses support (GTEx v10)
##
## Motivation:
## - Revisit "glycolysis" interpretation in baseline GTEx adipose in the same
##   spirit as the lipolysis/thermogenesis analyses:
##   - depot-specific
##   - sex-aware
##   - composition-aware (using first-pass CIBERSORT fractions when available)
##
## What this script produces:
## - Per-sample gene-set scores (VST-based mean gene-wise z) for MSigDB Hallmarks:
##   - HALLMARK_GLYCOLYSIS
##   - HALLMARK_OXIDATIVE_PHOSPHORYLATION
##   - HALLMARK_HYPOXIA
##   - HALLMARK_TNFA_SIGNALING_VIA_NFKB
##   - HALLMARK_INFLAMMATORY_RESPONSE
##   - HALLMARK_ADIPOGENESIS (context)
## - Female-only "menopause proxy" comparisons (pre/peri/post) with models:
##   - base: score ~ SMCENTER + group
##   - +cibersort: base + selected CIBERSORT fractions (if present for the depot)
## - Sex-stratified age-bin trajectories for key scores.
##
## Outputs:
## - GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_glycolysis/
## - GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_glycolysis/
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
  DESeq2,
  msigdbr
)

source("R/utils.R")

## -----------------------------------------------------------------------------
## Config
## -----------------------------------------------------------------------------
in_dir <- "GTExDatav10"
gct_sc <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
gct_vis <- file.path(in_dir, "gene_reads_v10_adipose_visceral_omentum.gct.gz")
sattr <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

out_dir <- "GTEx_v10_AT_analysis_out"
tab_dir <- file.path(out_dir, "tables", "next_steps_2026-02-28_glycolysis")
fig_dir <- file.path(out_dir, "figs", "next_steps_2026-02-28_glycolysis")
log_dir <- file.path(out_dir, "logs")
ensure_dirs(c(tab_dir, fig_dir, log_dir))

log_path <- file.path(log_dir, "next_steps_2026-02-28_glycolysis.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

age_levels <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")

depots <- list(
  list(key = "subcutaneous", label = "Adipose - Subcutaneous", gct_path = gct_sc, suffix = ""),
  list(key = "visceral", label = "Adipose - Visceral (Omentum)", gct_path = gct_vis, suffix = "_visceral")
)

hallmarks <- c(
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_ADIPOGENESIS"
)

## Prefer the most recent first-pass CIBERSORT outputs from this repo (TopN50, perm0).
cibersort_dir <- file.path(out_dir, "tables", "next_steps_2026-02-28_cibersort_fractions")
cibersort_paths <- list(
  subcutaneous = file.path(
    cibersort_dir,
    "cibersort_subcutaneous_perm0_topN50_minCells200__fractions_with_meta.tsv.gz"
  ),
  visceral = file.path(
    cibersort_dir,
    "cibersort_visceral_perm0_topN50_minCells200__fractions_with_meta.tsv.gz"
  )
)

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------
score_set_z <- function(expr_sym, genes) {
  genes <- intersect(unique(toupper(genes)), rownames(expr_sym))
  if (length(genes) < 5) return(rep(NA_real_, ncol(expr_sym)))
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

fmt_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 1e-3) return(format(x, scientific = TRUE, digits = 2))
  format(round(x, 4), nsmall = 4, trim = TRUE)
}

## -----------------------------------------------------------------------------
## MSigDB sets
## -----------------------------------------------------------------------------
cat("Loading MSigDB hallmarks via msigdbr...\n")
msig <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::filter(.data$gs_name %in% hallmarks) %>%
  dplyr::mutate(gene_symbol = toupper(.data$gene_symbol))

missing_sets <- setdiff(hallmarks, unique(msig$gs_name))
if (length(missing_sets)) warning("Missing hallmarks in msigdbr: ", paste(missing_sets, collapse = ", "))

gene_sets <- msig %>%
  dplyr::group_by(.data$gs_name) %>%
  dplyr::summarise(genes = list(unique(.data$gene_symbol)), .groups = "drop")
gene_sets <- stats::setNames(gene_sets$genes, gene_sets$gs_name)

## -----------------------------------------------------------------------------
## GTEx metadata
## -----------------------------------------------------------------------------
cat("Loading GTEx metadata...\n")
sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything()) %>%
  dplyr::mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])
subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything())
meta_full <- dplyr::left_join(sample_attr, subject_pheno, by = "SUBJID")

all_score_rows <- list()
all_stats_rows <- list()

for (dpt in depots) {
  depot_key <- dpt$key
  depot_label <- dpt$label
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
      depot = depot_key,
      sex = factor(sex, levels = c("Male", "Female")),
      age_bin_label = factor(as.character(age_bin_label), levels = age_levels),
      SMCENTER = factor(SMCENTER)
    ) %>%
    dplyr::filter(!is.na(age_bin_label))

  # VST for scoring (all samples in this depot)
  counts <- gct$counts[, meta$SAMPID, drop = FALSE]
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
  expr_sym_all <- collapse_to_symbol_mean(vst_mat, gct$gene_id, gct$gene_symbol, uppercase = TRUE)

  # Optional: CIBERSORT fractions for this depot (first pass)
  frac_tbl <- NULL
  frac_path <- cibersort_paths[[depot_key]]
  if (file.exists(frac_path)) {
    # The fractions-with-meta file includes both fractions and duplicated GTEx
    # metadata. Select only fraction columns to avoid join column collisions.
    frac_cols_possible <- c(
      "adipocyte", "aspc", "macrophage", "monocyte", "endothelial", "lec",
      "t_cell", "b_cell", "nk_cell", "dendritic_cell", "mast_cell",
      "pericyte", "smc", "mesothelium", "neutrophil", "endometrium"
    )
    frac_tbl <- readr::read_tsv(frac_path, show_col_types = FALSE) %>%
      dplyr::select(SAMPID, dplyr::any_of(frac_cols_possible))
    cat("Loaded CIBERSORT fractions: ", frac_path, "\n", sep = "")
  } else {
    cat("Note: CIBERSORT fractions not found for depot=", depot_key, " at: ", frac_path, "\n", sep = "")
  }

  # Compute scores within sex (z across samples in depot+sex)
  score_rows <- list()
  for (sx in levels(meta$sex)) {
    meta_sx <- meta %>% dplyr::filter(sex == sx)
    expr_sym <- expr_sym_all[, meta_sx$SAMPID, drop = FALSE]

    score_tbl <- tibble::tibble(SAMPID = colnames(expr_sym))
    for (gs in names(gene_sets)) {
      score_tbl[[gs]] <- score_set_z(expr_sym, gene_sets[[gs]])
    }
    score_tbl$sex <- sx
    score_rows[[length(score_rows) + 1]] <- score_tbl
  }
  score_tbl <- dplyr::bind_rows(score_rows) %>%
    dplyr::mutate(sex = factor(as.character(sex), levels = c("Male", "Female")))

  df0 <- meta %>%
    dplyr::left_join(score_tbl, by = c("SAMPID", "sex")) %>%
    dplyr::mutate(
      menopause_group_3 = dplyr::case_when(
        is.na(age_mid) ~ NA_character_,
        age_mid < 45 ~ "pre",
        age_mid <= 55 ~ "peri",
        age_mid > 55 ~ "post",
        TRUE ~ NA_character_
      ),
      menopause_group_3 = factor(menopause_group_3, levels = c("pre", "peri", "post"))
    )

  # Persist per-sample scores
  out_scores <- file.path(tab_dir, glue::glue("hallmark_scores_vst_z_by_sex_per_sample{suffix}.tsv.gz"))
  readr::write_tsv(df0 %>%
                     dplyr::select(SAMPID, depot, sex, age_bin_label, age_mid, menopause_group_3, SMCENTER, dplyr::all_of(names(gene_sets))),
                   out_scores)

  # Female-only: base vs +CIBERSORT models (pre/peri/post)
  df_f <- df0 %>% dplyr::filter(sex == "Female", !is.na(menopause_group_3))
  df_f2 <- df_f
  frac_cols_use <- character(0)
  if (!is.null(frac_tbl)) {
    df_f2 <- df_f2 %>% dplyr::left_join(frac_tbl, by = "SAMPID")
    frac_cols <- setdiff(colnames(frac_tbl), "SAMPID")
    # Keep a compact set of broad fraction covariates if present.
    # Use a small, interpretable subset to reduce multicollinearity (CIBERSORT
    # fractions sum to ~1 if you include too many).
    frac_cols_use <- intersect(frac_cols, c("adipocyte", "aspc", "macrophage", "endothelial", "t_cell"))
  }

  model_defs <- list(
    base = "score ~ SMCENTER + menopause_group_3"
  )
  if (length(frac_cols_use)) {
    model_defs[["adj_cibersort"]] <- paste0("score ~ SMCENTER + menopause_group_3 + ", paste(frac_cols_use, collapse = " + "))
  }

  stats_rows <- list()
  for (gs in names(gene_sets)) {
    for (model_nm in names(model_defs)) {
      df <- df_f2 %>% dplyr::mutate(score = .data[[gs]]) %>% dplyr::filter(!is.na(score))
      fit <- tryCatch(stats::lm(stats::as.formula(model_defs[[model_nm]]), data = df), error = function(e) NULL)
      if (is.null(fit)) next
      peri <- lm_contrast(fit, "menopause_group_3peri", "")
      post <- lm_contrast(fit, "menopause_group_3post", "")
      post_peri <- lm_contrast(fit, "menopause_group_3post", "menopause_group_3peri")
      stats_rows[[length(stats_rows) + 1]] <- tibble::tibble(
        depot = depot_key,
        gene_set = gs,
        model = model_nm,
        n = nrow(df),
        peri_vs_pre_est = peri$estimate, peri_vs_pre_p = peri$p,
        post_vs_pre_est = post$estimate, post_vs_pre_p = post$p,
        post_vs_peri_est = post_peri$estimate, post_vs_peri_p = post_peri$p,
        cibersort_covariates = if (model_nm == "adj_cibersort") paste(frac_cols_use, collapse = ",") else ""
      )
    }
  }
  stats_tbl <- dplyr::bind_rows(stats_rows) %>%
    dplyr::mutate(
      peri_vs_pre_p_fmt = vapply(.data$peri_vs_pre_p, fmt_p, character(1)),
      post_vs_pre_p_fmt = vapply(.data$post_vs_pre_p, fmt_p, character(1)),
      post_vs_peri_p_fmt = vapply(.data$post_vs_peri_p, fmt_p, character(1))
    )
  all_stats_rows[[length(all_stats_rows) + 1]] <- stats_tbl
  readr::write_tsv(stats_tbl, file.path(tab_dir, glue::glue("female_pre_peri_post_lm_stats_hallmarks{suffix}.tsv")))

  # Plot (female-only, by group) for glycolysis and a couple context sets
  to_plot <- intersect(
    c("HALLMARK_GLYCOLYSIS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_HYPOXIA", "HALLMARK_ADIPOGENESIS"),
    names(gene_sets)
  )
  df_long <- df_f %>%
    dplyr::select(SAMPID, menopause_group_3, dplyr::all_of(to_plot)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(to_plot), names_to = "gene_set", values_to = "score") %>%
    dplyr::filter(!is.na(score))

  p_box <- ggplot2::ggplot(df_long, ggplot2::aes(x = menopause_group_3, y = score, color = menopause_group_3)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.25) +
    ggplot2::geom_jitter(width = 0.18, alpha = 0.25, size = 0.8) +
    ggplot2::facet_wrap(~ gene_set, scales = "free_y", ncol = 2) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::scale_color_manual(values = c(pre = "#1b9e77", peri = "#d95f02", post = "#7570b3")) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(
      title = glue::glue("Hallmark score shifts across menopause proxy groups (Female; {depot_label})"),
      subtitle = "Scores are mean gene-wise z across set genes; z computed within depot+sex on VST expression",
      x = NULL,
      y = "Mean z-score across genes"
    )
  ggplot2::ggsave(file.path(fig_dir, glue::glue("female_menopause_groups_hallmarks_boxplots{suffix}.png")),
                  p_box, width = 10.5, height = 6.2, dpi = 200)

  # Sex trajectories across age bins (means with 95% CI)
  traj_sets <- intersect(
    c("HALLMARK_GLYCOLYSIS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_HYPOXIA"),
    names(gene_sets)
  )
  traj_long <- df0 %>%
    dplyr::select(SAMPID, depot, sex, age_bin_label, dplyr::all_of(traj_sets)) %>%
    tidyr::pivot_longer(cols = dplyr::all_of(traj_sets), names_to = "gene_set", values_to = "score") %>%
    dplyr::filter(!is.na(score))

  traj_sum <- traj_long %>%
    dplyr::group_by(depot, gene_set, sex, age_bin_label) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(score, na.rm = TRUE),
      se = sd(score, na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    ) %>%
    dplyr::filter(n >= 8)

  readr::write_tsv(traj_sum, file.path(tab_dir, glue::glue("scores_by_sex_agebin_summary_hallmarks{suffix}.tsv")))

  p_traj <- ggplot2::ggplot(traj_sum, ggplot2::aes(x = age_bin_label, y = mean, color = sex, group = sex)) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::geom_point(size = 1.6) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), width = 0.15, linewidth = 0.4) +
    ggplot2::facet_wrap(~ gene_set, scales = "free_y", ncol = 1) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "top") +
    ggplot2::labs(
      title = glue::glue("Hallmark trajectories across GTEx age bins (Male vs Female; {depot_label})"),
      subtitle = "Scores computed within depot+sex (mean gene-wise z on VST); error bars are 95% CI of mean",
      x = "GTEx age bin",
      y = "Mean z-score"
    )
  ggplot2::ggsave(file.path(fig_dir, glue::glue("sex_agebin_trajectories_hallmarks{suffix}.png")),
                  p_traj, width = 8.8, height = 7.6, dpi = 200)

  all_score_rows[[length(all_score_rows) + 1]] <- df0 %>%
    dplyr::select(SAMPID, depot, sex, age_bin_label, age_mid, menopause_group_3, SMCENTER, dplyr::all_of(names(gene_sets)))
}

scores_all <- dplyr::bind_rows(all_score_rows)
readr::write_tsv(scores_all, file.path(tab_dir, "hallmark_scores_vst_z_by_sex_all_depots.tsv.gz"))

stats_all <- dplyr::bind_rows(all_stats_rows)
readr::write_tsv(stats_all, file.path(tab_dir, "female_pre_peri_post_lm_stats_hallmarks_ALL.tsv"))

cat("\nWrote outputs:\n")
cat("  - ", tab_dir, "\n", sep = "")
cat("  - ", fig_dir, "\n", sep = "")
cat("Log: ", log_path, "\n", sep = "")
