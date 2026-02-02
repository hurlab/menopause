################################################################################
## Next steps (2026-02-02): Sex interaction analysis on gene-set/module scores
##
## Motivation (Kenichi email):
## - Male vs female comparisons can provide context but are imperfect controls.
## - We quantify sex-by-age interactions for *scores* (Kenichi sets, composition
##   proxies, adrenergic modules, senMayoLike signatures) in both depots.
##
## Model family:
## - score ~ SMCENTER + sex + age_group + sex:age_group
## - Optional sensitivity: add SMRIN + SMTSISCH + composition proxies where relevant
##
## Age-group comparisons (by GTEx age bins):
## - 40-49 vs 60-69 (strongest “menopause-adjacent” contrast)
## - 30-39 vs 60-69 (wider contrast)
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

log_path <- file.path(log_dir, "next_steps_2026-02-02_sex_interaction_scores.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

depots <- list(
  list(key = "subcutaneous", label = "Adipose - Subcutaneous", gct_path = gct_sc, suffix = ""),
  list(key = "visceral", label = "Adipose - Visceral (Omentum)", gct_path = gct_vis, suffix = "_visceral")
)

comparisons <- list(
  list(name = "40-49_vs_60-69", young = "40-49", old = "60-69"),
  list(name = "30-39_vs_60-69", young = "30-39", old = "60-69")
)

min_n_per_sex_group <- 8

## -----------------------------------------------------------------------------
## Define score sets (reusing prior work)
## -----------------------------------------------------------------------------
kenichi_sets <- list(
  acute_beta_adrenergic = c("NR4A1", "NR4A2", "NR4A3", "FOS", "JUN", "EGR1"),
  adrenergic_receptors = c("ADRB1", "ADRB2", "ADRB3", "ADRA2A", "ADRA2B", "ADRA2C"),
  lipolysis_core = c("PNPLA2", "LIPE", "MGLL", "ABHD5", "PLIN1", "G0S2", "LPL", "FABP4"),
  thermogenesis_program = c("UCP1", "PPARGC1A", "PPARGC1B", "CIDEA", "DIO2", "VEGFA")
)
kenichi_sets <- purrr::map(kenichi_sets, ~ unique(toupper(.x)))

composition_sets <- list(
  adipocyte = c("ADIPOQ", "PLIN1", "FABP4"),
  macrophage_mono = c("LST1", "C1QC", "TYROBP"),
  fibroblast_ecm = c("COL1A1", "COL3A1", "DCN")
)
composition_sets <- purrr::map(composition_sets, ~ unique(toupper(.x)))

adrenergic_modules <- list(
  camp_pka_feedback = c("PDE4D", "PDE4B", "PDE3B", "RGS2", "DUSP1", "DUSP2", "ATF3"),
  receptor_desensitization = c("ADRBK1", "ADRBK2", "ARRB1", "ARRB2", "RGS2"),
  alpha2_antilipolytic_axis = c("ADRA2A", "ADRA2B", "ADRA2C", "PDE3B"),
  camp_core = c("ADCY3", "ADCY6", "GNAS", "PRKACA", "PRKACB", "PRKAR1A", "PRKAR2B"),
  immediate_early = c("NR4A1", "NR4A2", "NR4A3", "FOS", "JUN", "EGR1", "DUSP1", "ATF3")
)
adrenergic_modules <- purrr::map(adrenergic_modules, ~ unique(toupper(.x)))

sig_v2 <- readr::read_tsv("references/menopause_signature_senMayoLike_v2.tsv", show_col_types = FALSE) %>%
  dplyr::mutate(gene_symbol = toupper(as.character(gene_symbol)))
sig2_up <- sig_v2 %>% dplyr::filter(direction_in_post == "up") %>% dplyr::pull(gene_symbol)
sig2_down <- sig_v2 %>% dplyr::filter(direction_in_post == "down") %>% dplyr::pull(gene_symbol)

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

fit_and_extract <- function(df, formula_txt, effect_label, comparison_id, depot_key, score_name) {
  fit <- tryCatch(stats::lm(stats::as.formula(formula_txt), data = df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  cf <- summary(fit)$coefficients

  # Coeff naming: baseline sex=Male, age_group=young
  # - age_groupold: male aging effect
  # - sexFemale: female-male difference at young
  # - sexFemale:age_groupold: interaction (female aging difference vs male)
  term_age <- "age_groupold"
  term_int <- "sexFemale:age_groupold"

  get_row <- function(term) {
    if (!(term %in% rownames(cf))) return(NULL)
    tibble::tibble(
      estimate = unname(cf[term, "Estimate"]),
      se = unname(cf[term, "Std. Error"]),
      t = unname(cf[term, "t value"]),
      p = unname(cf[term, "Pr(>|t|)"])
    )
  }

  male <- get_row(term_age)
  inter <- get_row(term_int)
  if (is.null(male) || is.null(inter)) return(NULL)

  # Female aging effect = male + interaction
  est_f <- male$estimate + inter$estimate
  # Approximate SE for sum using vcov
  vc <- stats::vcov(fit)
  if (!all(c(term_age, term_int) %in% rownames(vc))) return(NULL)
  se_f <- sqrt(vc[term_age, term_age] + vc[term_int, term_int] + 2 * vc[term_age, term_int])
  t_f <- est_f / se_f
  p_f <- 2 * stats::pt(abs(t_f), df = fit$df.residual, lower.tail = FALSE)

  dplyr::bind_rows(
    male %>% dplyr::mutate(effect = "male_old_vs_young"),
    tibble::tibble(estimate = est_f, se = se_f, t = t_f, p = p_f) %>% dplyr::mutate(effect = "female_old_vs_young"),
    inter %>% dplyr::mutate(effect = "sex_interaction_female_minus_male")
  ) %>%
    dplyr::mutate(
      depot = depot_key,
      comparison = comparison_id,
      score_name = score_name,
      model = effect_label,
      n = nrow(df)
    ) %>%
    dplyr::relocate(depot, comparison, score_name, model, effect, n)
}

## -----------------------------------------------------------------------------
## Metadata
## -----------------------------------------------------------------------------
sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything()) %>%
  dplyr::mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])
subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything())
meta_full <- dplyr::left_join(sample_attr, subject_pheno, by = "SUBJID")

all_rows <- list()

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
      DTHHRDY = factor(DTHHRDY)
    )

  counts <- gct$counts[, meta$SAMPID, drop = FALSE]
  counts <- collapse_dupe_rowsum(counts)
  counts <- round(counts)
  storage.mode(counts) <- "integer"

  # VST for scoring
  dds0 <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = data.frame(row.names = colnames(counts)),
    design = ~ 1
  )
  vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds0, blind = TRUE))
  expr_sym <- collapse_to_symbol_mean(vst_mat, gct$gene_id, gct$gene_symbol, uppercase = TRUE)

  # Compute all scores for all samples in this depot
  score_tbl <- tibble::tibble(SAMPID = colnames(expr_sym)) %>%
    dplyr::mutate(
      senMayoLike_v2_z_signed = score_set_z(expr_sym, sig2_up) - score_set_z(expr_sym, sig2_down)
    )

  add_scores <- function(tbl, sets) {
    for (nm in names(sets)) {
      tbl[[nm]] <- score_set_z(expr_sym, sets[[nm]])
    }
    tbl
  }
  score_tbl <- score_tbl %>%
    add_scores(kenichi_sets) %>%
    add_scores(composition_sets) %>%
    add_scores(adrenergic_modules)

  meta2 <- meta %>%
    dplyr::left_join(score_tbl, by = c("SAMPID" = "SAMPID"))

  for (cmp in comparisons) {
    comparison_id <- cmp$name
    young_age <- cmp$young
    old_age <- cmp$old

    df <- meta2 %>%
      dplyr::filter(age_bin_label %in% c(young_age, old_age)) %>%
      dplyr::mutate(
        age_group = factor(ifelse(age_bin_label == young_age, "young", "old"), levels = c("young", "old"))
      )

    n_tbl <- df %>%
      dplyr::count(sex, age_group, name = "n") %>%
      tidyr::complete(sex, age_group, fill = list(n = 0L))

    cat(glue::glue("Comparison: {comparison_id} ({young_age} vs {old_age})\n"))
    print(n_tbl)

    if (any(n_tbl$n < min_n_per_sex_group)) {
      cat("Skipping (insufficient n per sex x age_group)\n\n")
      next
    }

    # model variants
    model_defs <- list(
      base = "score ~ SMCENTER + sex + age_group + sex:age_group",
      tech = "score ~ SMCENTER + SMRIN + SMTSISCH + sex + age_group + sex:age_group",
      tech_comp = "score ~ SMCENTER + SMRIN + SMTSISCH + adipocyte + macrophage_mono + fibroblast_ecm + sex + age_group + sex:age_group"
    )

    score_names <- c(
      names(kenichi_sets),
      names(composition_sets),
      names(adrenergic_modules),
      "senMayoLike_v2_z_signed"
    )

    for (sc in score_names) {
      df_sc <- df %>% dplyr::mutate(score = .data[[sc]])
      if (all(is.na(df_sc$score))) next

      for (mdl in names(model_defs)) {
        # require non-missing covariates for the model
        req_cols <- c("score", "SMCENTER", "sex", "age_group")
        if (mdl %in% c("tech", "tech_comp")) req_cols <- c(req_cols, "SMRIN", "SMTSISCH")
        if (mdl == "tech_comp") req_cols <- c(req_cols, "adipocyte", "macrophage_mono", "fibroblast_ecm")
        df_fit <- df_sc %>% dplyr::filter(if_all(dplyr::all_of(req_cols), ~ !is.na(.x)))
        if (nrow(df_fit) < 2 * min_n_per_sex_group * 2) next

        res <- fit_and_extract(df_fit, model_defs[[mdl]], mdl, comparison_id, depot_key, sc)
        if (!is.null(res)) all_rows[[length(all_rows) + 1]] <- res
      }
    }
  }
}

res_all <- dplyr::bind_rows(all_rows) %>%
  dplyr::mutate(padj_within = ave(p, depot, comparison, model, effect, FUN = p.adjust, method = "BH"))

readr::write_tsv(res_all, file.path(tab_dir, "sex_interaction_score_effects.tsv"))

## -----------------------------------------------------------------------------
## Simple “heatmap-like” plot: interaction p-values for key scores (subq)
## -----------------------------------------------------------------------------
plot_scores <- c(
  "lipolysis_core", "thermogenesis_program",
  "adipocyte", "macrophage_mono", "fibroblast_ecm",
  "camp_core", "receptor_desensitization",
  "senMayoLike_v2_z_signed"
)

plot_df <- res_all %>%
  dplyr::filter(score_name %in% plot_scores) %>%
  dplyr::filter(effect == "sex_interaction_female_minus_male") %>%
  dplyr::filter(model == "base") %>%
  dplyr::mutate(neglog10p = -log10(pmax(p, 1e-300)))

p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = comparison, y = score_name, fill = neglog10p)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::facet_wrap(~ depot, ncol = 1) +
  ggplot2::scale_fill_viridis_c(option = "C", direction = 1) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    title = "Sex × age-bin interaction strength (score-level; base model)",
    subtitle = "Fill = -log10(p) for interaction term (female age-effect minus male age-effect)",
    x = "Age-bin comparison (young vs old)",
    y = "Score"
  )

ggplot2::ggsave(file.path(fig_dir, "sex_interaction_scores_heatmap_base.png"), p, width = 10, height = 6, dpi = 150)

cat("Done. Wrote: ", file.path(tab_dir, "sex_interaction_score_effects.tsv"), "\n", sep = "")
cat("Log: ", log_path, "\n", sep = "")

