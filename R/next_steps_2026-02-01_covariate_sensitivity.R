################################################################################
## Next steps (2026-02-01): Covariate Sensitivity Ladder
##
## Motivation (see A quick summary.docx):
## - PC-level variation in GTEx adipose is strongly associated with ischemic time
##   (SMTSISCH) and RNA integrity (SMRIN).
## - Existing Kenichi panel analyses adjust for SMCENTER; here we stress-test the
##   menopause-proxy associations against additional technical covariates.
##
## Outputs (dated subfolders):
## - GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_covariate_sensitivity/
## - GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_covariate_sensitivity/
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

options(error = function() {
  traceback(2)
  quit(status = 1)
})

## -----------------------------------------------------------------------------
## Paths + output folders
## -----------------------------------------------------------------------------
in_dir <- "GTExDatav10"
gct_sc <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
gct_vis <- file.path(in_dir, "gene_reads_v10_adipose_visceral_omentum.gct.gz")
sattr <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

out_dir <- "GTEx_v10_AT_analysis_out"
tab_dir <- file.path(out_dir, "tables", "next_steps_2026-02-01_covariate_sensitivity")
fig_dir <- file.path(out_dir, "figs", "next_steps_2026-02-01_covariate_sensitivity")
log_dir <- file.path(out_dir, "logs")
ensure_dirs(c(tab_dir, fig_dir, log_dir))

log_path <- file.path(log_dir, "next_steps_2026-02-01_covariate_sensitivity.log")
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

## Reuse the Kenichi panel gene sets as the main stress-test targets
gene_sets <- list(
  acute_beta_adrenergic = c("NR4A1", "NR4A2", "NR4A3", "FOS", "JUN", "EGR1"),
  adrenergic_receptors = c("ADRB1", "ADRB2", "ADRB3", "ADRA2A", "ADRA2B", "ADRA2C"),
  lipolysis_core = c("PNPLA2", "LIPE", "MGLL", "ABHD5", "PLIN1", "G0S2", "LPL", "FABP4"),
  thermogenesis_program = c("UCP1", "PPARGC1A", "PPARGC1B", "CIDEA", "DIO2", "VEGFA")
)
gene_sets <- purrr::map(gene_sets, ~ unique(toupper(.x)))
panel_symbols <- unique(unlist(gene_sets, use.names = FALSE))

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

fit_lm_safe <- function(formula, data) {
  tryCatch(stats::lm(formula, data = data), error = function(e) NULL)
}

tidy_group_contrasts <- function(fit, group_var) {
  if (is.null(fit)) return(tibble::tibble())
  cn <- names(stats::coef(fit))
  # treatment coding: baseline is first factor level
  # terms like group_varperi, group_varpost
  peri_term <- paste0(group_var, "peri")
  post_term <- paste0(group_var, "post")
  has_peri <- peri_term %in% cn
  has_post <- post_term %in% cn
  out <- list()
  if (has_peri) {
    out$peri_vs_pre <- lm_contrast(fit, peri_term, "")
  }
  if (has_post) {
    out$post_vs_pre <- lm_contrast(fit, post_term, "")
  }
  if (has_peri && has_post) {
    out$post_vs_peri <- lm_contrast(fit, post_term, peri_term)
  }
  if (!length(out)) return(tibble::tibble())
  dplyr::bind_rows(out, .id = "contrast")
}

## -----------------------------------------------------------------------------
## Load metadata
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("COVARIATE SENSITIVITY LADDER (GTEx v10)\n")
cat("===========================================\n\n")

cat("Loading metadata...\n")
sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything()) %>%
  dplyr::mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])
subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything())
meta_full <- dplyr::left_join(sample_attr, subject_pheno, by = "SUBJID")

## -----------------------------------------------------------------------------
## Model ladder definitions
## -----------------------------------------------------------------------------
## Note: BMI is not present in the v10 SubjectPhenotypesDS in this repo.
## We stress-test technical covariates that exist here: SMRIN, SMTSISCH, DTHHRDY, SMGEBTCH.
model_ladder <- list(
  base = "score ~ SMCENTER + group",
  rin = "score ~ SMCENTER + SMRIN + group",
  rin_isch = "score ~ SMCENTER + SMRIN + SMTSISCH + group",
  rin_isch_hardy = "score ~ SMCENTER + SMRIN + SMTSISCH + DTHHRDY + group",
  rin_isch_hardy_ge = "score ~ SMCENTER + SMRIN + SMTSISCH + DTHHRDY + SMGEBTCH + group"
)

## -----------------------------------------------------------------------------
## Main loop (both depots)
## -----------------------------------------------------------------------------
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
    dplyr::mutate(
      SMCENTER = factor(SMCENTER),
      DTHHRDY = factor(DTHHRDY),
      SMGEBTCH = factor(SMGEBTCH)
    )

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

  cat(sprintf("Female samples (all age bins): %d\n", nrow(meta_f)))
  cat("Menopause proxy (3-group) counts:\n")
  print(table(meta_f$menopause_group_3, useNA = "ifany"))
  cat("Menopause proxy (strict2) counts:\n")
  print(table(meta_f$menopause_group_strict2, useNA = "ifany"))

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

  ## Per-sample gene-set scores
  scores <- purrr::imap_dfc(gene_sets, ~ tibble::tibble(!!.y := score_set_z(expr_sym, .x)))
  rownames(scores) <- colnames(expr_sym)

  ## Join scores to metadata
  meta_s <- meta_f %>%
    dplyr::select(SAMPID, SMCENTER, SMRIN, SMTSISCH, DTHHRDY, SMGEBTCH, menopause_group_3, menopause_group_strict2) %>%
    dplyr::mutate(SAMPID = as.character(SAMPID))
  score_df <- dplyr::bind_cols(meta_s, scores[meta_s$SAMPID, , drop = FALSE])

  ## -----------------------
  ## Gene-set sensitivity (3-group)
  ## -----------------------
  gene_set_rows <- list()
  for (set_name in names(gene_sets)) {
    df <- score_df %>%
      dplyr::filter(!is.na(menopause_group_3)) %>%
      dplyr::mutate(group = menopause_group_3, score = .data[[set_name]])

    for (mdl in names(model_ladder)) {
      ftxt <- model_ladder[[mdl]]
      fml <- stats::as.formula(ftxt)
      fit <- fit_lm_safe(fml, df)
      ct <- tidy_group_contrasts(fit, group_var = "group")
      if (!nrow(ct)) next
      gene_set_rows[[length(gene_set_rows) + 1]] <- ct %>%
        dplyr::mutate(
          depot = depot_key,
          menopause_def = "3group",
          set_name = set_name,
          model = mdl,
          n = nrow(df)
        )
    }
  }
  gene_set_stats_3 <- dplyr::bind_rows(gene_set_rows) %>%
    dplyr::relocate(depot, menopause_def, set_name, model, contrast, n)

  ## -----------------------
  ## Gene-set sensitivity (strict2)
  ## -----------------------
  gene_set_rows <- list()
  model_ladder_strict <- list(
    base = "score ~ SMCENTER + group",
    rin_isch = "score ~ SMCENTER + SMRIN + SMTSISCH + group",
    rin_isch_hardy = "score ~ SMCENTER + SMRIN + SMTSISCH + DTHHRDY + group",
    rin_isch_hardy_ge = "score ~ SMCENTER + SMRIN + SMTSISCH + DTHHRDY + SMGEBTCH + group"
  )
  for (set_name in names(gene_sets)) {
    df <- score_df %>%
      dplyr::filter(!is.na(menopause_group_strict2)) %>%
      dplyr::mutate(group = menopause_group_strict2, score = .data[[set_name]])

    for (mdl in names(model_ladder_strict)) {
      ftxt <- model_ladder_strict[[mdl]]
      fml <- stats::as.formula(ftxt)
      fit <- fit_lm_safe(fml, df)
      if (is.null(fit)) next
      # in strict2, term is grouppost
      ct <- lm_contrast(fit, "grouppost", "")
      gene_set_rows[[length(gene_set_rows) + 1]] <- ct %>%
        dplyr::mutate(
          depot = depot_key,
          menopause_def = "strict2",
          set_name = set_name,
          model = mdl,
          contrast = "post_vs_pre",
          n = nrow(df)
        )
    }
  }
  gene_set_stats_strict <- dplyr::bind_rows(gene_set_rows) %>%
    dplyr::relocate(depot, menopause_def, set_name, model, contrast, n)

  gene_set_stats_all <- dplyr::bind_rows(gene_set_stats_3, gene_set_stats_strict) %>%
    dplyr::mutate(
      conf_low = estimate - 1.96 * se,
      conf_high = estimate + 1.96 * se
    )

  out_tsv <- file.path(tab_dir, glue::glue("gene_set_covariate_sensitivity{suffix}.tsv"))
  readr::write_tsv(gene_set_stats_all, out_tsv)

  ## Plot: effect estimate across model ladder for key contrasts
  plot_df <- gene_set_stats_all %>%
    dplyr::filter(contrast %in% c("peri_vs_pre", "post_vs_pre")) %>%
    dplyr::mutate(model = factor(model, levels = unique(model)))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = model, y = estimate, ymin = conf_low, ymax = conf_high)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey70") +
    ggplot2::geom_pointrange(size = 0.25) +
    ggplot2::facet_grid(menopause_def ~ set_name, scales = "free_y") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = paste0("Covariate sensitivity: menopause proxy effects (", depot_label, ", female)"),
      x = "Model ladder",
      y = "Effect estimate (vs baseline 'pre')"
    )
  ggplot2::ggsave(
    filename = file.path(fig_dir, glue::glue("gene_set_covariate_sensitivity{suffix}.png")),
    plot = p,
    width = 12,
    height = 6,
    dpi = 150
  )

  ## -----------------------
  ## Gene-level sensitivity (panel genes; VST linear model)
  ## -----------------------
  cat("Gene-level sensitivity (panel genes; VST lm)...\n")
  panel_present <- intersect(panel_symbols, rownames(expr_sym))
  expr_panel <- expr_sym[panel_present, , drop = FALSE]

  strict_df <- meta_f %>%
    dplyr::filter(!is.na(menopause_group_strict2)) %>%
    dplyr::select(SAMPID, SMCENTER, SMRIN, SMTSISCH, DTHHRDY, SMGEBTCH, menopause_group_strict2) %>%
    dplyr::mutate(group = menopause_group_strict2)
  rownames(strict_df) <- strict_df$SAMPID

  gene_rows <- list()
  for (g in rownames(expr_panel)) {
    df <- strict_df
    df$expr <- as.numeric(expr_panel[g, df$SAMPID])
    # base vs full covariate model (VST-scale)
    fits <- list(
      base = fit_lm_safe(expr ~ SMCENTER + group, df),
      rin_isch_hardy_ge = fit_lm_safe(expr ~ SMCENTER + SMRIN + SMTSISCH + DTHHRDY + SMGEBTCH + group, df)
    )
    for (nm in names(fits)) {
      fit <- fits[[nm]]
      if (is.null(fit)) next
      ct <- lm_contrast(fit, "grouppost", "")
      gene_rows[[length(gene_rows) + 1]] <- ct %>%
        dplyr::mutate(
          depot = depot_key,
          gene_symbol = g,
          model = nm,
          contrast = "post_vs_pre",
          n = nrow(df)
        )
    }
  }
  gene_stats <- dplyr::bind_rows(gene_rows) %>%
    dplyr::mutate(conf_low = estimate - 1.96 * se, conf_high = estimate + 1.96 * se) %>%
    dplyr::arrange(model, p) %>%
    dplyr::relocate(depot, gene_symbol, model, contrast, n)
  readr::write_tsv(gene_stats, file.path(tab_dir, glue::glue("panel_gene_vst_lm_sensitivity{suffix}.tsv")))

  ## -----------------------
  ## DESeq2 sensitivity (strict2; base vs expanded covariates)
  ## -----------------------
  cat("DESeq2 sensitivity (strict2; base vs expanded)...\n")
  meta_strict <- meta_f %>%
    dplyr::filter(!is.na(menopause_group_strict2)) %>%
    dplyr::mutate(group = menopause_group_strict2) %>%
    dplyr::filter(!is.na(SMRIN), !is.na(SMTSISCH), !is.na(DTHHRDY), !is.na(SMGEBTCH)) %>%
    dplyr::mutate(SAMPID = as.character(SAMPID)) %>%
    dplyr::filter(SAMPID %in% colnames(counts))

  counts_strict <- counts[, as.character(meta_strict$SAMPID), drop = FALSE]
  stopifnot(all(colnames(counts_strict) == meta_strict$SAMPID))
  coldata <- meta_strict %>%
    dplyr::select(SAMPID, SMCENTER, SMRIN, SMTSISCH, DTHHRDY, SMGEBTCH, group) %>%
    dplyr::mutate(
      SMCENTER = factor(SMCENTER),
      DTHHRDY = factor(DTHHRDY),
      SMGEBTCH = factor(SMGEBTCH),
      group = factor(group, levels = c("pre", "post")),
      SMRIN_z = as.numeric(scale(SMRIN)),
      SMTSISCH_z = as.numeric(scale(SMTSISCH))
    ) %>%
    as.data.frame()
  rownames(coldata) <- coldata$SAMPID

  ## Reduce the gene universe for faster sensitivity runs:
  ## - retain top-N expressed genes (for stable dispersion estimation)
  ## - always include panel genes (mapped by symbol)
  N_TOP <- 5000
  gene_id_by_symbol <- gct$gene_id[match(toupper(gct$gene_symbol), toupper(gct$gene_symbol))]
  panel_gene_ids <- gct$gene_id[toupper(gct$gene_symbol) %in% panel_symbols]
  means <- rowMeans(counts_strict, na.rm = TRUE)
  top_ids <- names(sort(means, decreasing = TRUE))[seq_len(min(N_TOP, length(means)))]
  keep_ids <- unique(c(top_ids, panel_gene_ids))
  keep_ids <- keep_ids[keep_ids %in% rownames(counts_strict)]
  counts_strict_red <- counts_strict[keep_ids, , drop = FALSE]

  run_deseq <- function(design_formula) {
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = counts_strict_red,
      colData = coldata,
      design = design_formula
    )
    dds <- DESeq2::DESeq(dds, quiet = TRUE, minReplicatesForReplace = Inf)
    res <- DESeq2::results(dds, contrast = c("group", "post", "pre"))
    as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      dplyr::mutate(
        gene_symbol = toupper(gct$gene_symbol[match(.data$gene_id, gct$gene_id)]),
        design = deparse(design_formula)
      ) %>%
      dplyr::relocate(design, gene_id, gene_symbol)
  }

  safe_run_deseq <- function(formulas) {
    for (f in formulas) {
      out <- tryCatch(run_deseq(f), error = function(e) NULL)
      if (!is.null(out)) return(out)
    }
    stop("All DESeq2 designs failed (rank-deficient or other errors).")
  }

  res_base <- run_deseq(~ SMCENTER + group)
  # Full model can be rank-deficient in some depots after filtering; fall back gracefully.
  res_full <- safe_run_deseq(list(
    ~ SMCENTER + SMRIN_z + SMTSISCH_z + DTHHRDY + SMGEBTCH + group,
    ~ SMCENTER + SMRIN_z + SMTSISCH_z + DTHHRDY + group,
    ~ SMCENTER + SMRIN_z + SMTSISCH_z + group
  ))

  out_all <- dplyr::bind_rows(res_base, res_full)
  readr::write_tsv(out_all, file.path(tab_dir, glue::glue("DESeq2_strict2_base_vs_full_all_genes{suffix}.tsv")))

  out_panel <- out_all %>%
    dplyr::filter(!is.na(gene_symbol), gene_symbol %in% panel_symbols) %>%
    dplyr::arrange(design, padj)
  readr::write_tsv(out_panel, file.path(tab_dir, glue::glue("DESeq2_strict2_base_vs_full_PANEL{suffix}.tsv")))

  ## Quick markdown summary (per depot)
  md <- c(
    glue::glue("# Covariate Sensitivity Summary ({depot_label}; female-only)"),
    "",
    "This folder stress-tests menopause-proxy effects against key GTEx technical covariates available in this repo.",
    "",
    "## Inputs",
    glue::glue("- Counts: `GTExDatav10/{basename(dpt$gct_path)}`"),
    "- Sample attributes: `GTExDatav10/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt`",
    "- Subject phenotypes: `GTExDatav10/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt`",
    "",
    "## Covariates tested (availability-limited)",
    "- SMCENTER (factor)",
    "- SMRIN (numeric)",
    "- SMTSISCH (numeric; ischemic time)",
    "- DTHHRDY (factor; Hardy scale)",
    "- SMGEBTCH (factor; gene expression batch)",
    "",
    "## Outputs",
    glue::glue("- Gene-set sensitivity ladder: `gene_set_covariate_sensitivity{suffix}.tsv`"),
    glue::glue("- Gene-set effect plot: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_covariate_sensitivity/gene_set_covariate_sensitivity{suffix}.png`"),
    glue::glue("- Panel-gene VST lm sensitivity (base vs full): `panel_gene_vst_lm_sensitivity{suffix}.tsv`"),
    glue::glue("- DESeq2 strict2 base vs full (all genes + panel subset): `DESeq2_strict2_base_vs_full_*{suffix}.tsv`"),
    "",
    "Generated: ",
    as.character(Sys.time())
  )
  writeLines(md, file.path(tab_dir, glue::glue("README_covariate_sensitivity{suffix}.md")))
}

cat("\nDone. Log written to: ", log_path, "\n", sep = "")
