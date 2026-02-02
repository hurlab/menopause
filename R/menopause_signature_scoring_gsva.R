################################################################################
## Menopause Signature Scoring (GSVA-like; SenMayo-style)
##
## Goal:
## - Compute per-sample menopause signature activity scores in GTEx v10
##   subcutaneous adipose (female-only), using VST expression and multiple
##   GSVA package methods.
## - Use a curated/edited signature list in references/ so gene refinement is
##   explicit and versioned.
##
## Inputs:
## - GTExDatav10/ gene_reads_v10_adipose_subcutaneous.gct.gz
## - GTExDatav10/ GTEx_Analysis_v10_Annotations_*DS.txt
## - references/menopause_signature_curated.tsv
## - (optional) Approach 1 DE table for direction:
##   GTEx_v10_AT_analysis_out/tables/approach1_improved_age_proxy/DE_menopause_proxy_improved.tsv
##
## Outputs:
## - GTEx_v10_AT_analysis_out/tables/menopause_signature_scoring_gsva/
## - GTEx_v10_AT_analysis_out/figs/menopause_signature_scoring_gsva/
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
  readr, dplyr, tidyr, stringr, tibble, purrr,
  ggplot2,
  DESeq2,
  GSVA
)

source("R/utils.R")

## -----------------------------------------------------------------------------
## Paths
## -----------------------------------------------------------------------------
in_dir <- "GTExDatav10"
gct_sc <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
gct_vis <- file.path(in_dir, "gene_reads_v10_adipose_visceral_omentum.gct.gz")
sattr <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

out_dir <- "GTEx_v10_AT_analysis_out"
fig_dir <- file.path(out_dir, "figs", "menopause_signature_scoring_gsva")
tab_dir <- file.path(out_dir, "tables", "menopause_signature_scoring_gsva")
ensure_dirs(c(fig_dir, tab_dir))

cur_sig_path <- file.path("references", "menopause_signature_curated.tsv")
approach1_de_path <- file.path(out_dir, "tables", "approach1_improved_age_proxy", "DE_menopause_proxy_improved.tsv")

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
## Signature loading + direction (up/down)
## -----------------------------------------------------------------------------
if (!file.exists(cur_sig_path)) {
  stop("Missing curated signature file: ", cur_sig_path)
}

sig_df <- readr::read_tsv(cur_sig_path, show_col_types = FALSE)
if (!all(c("gene_symbol", "keep") %in% names(sig_df))) {
  stop("Curated signature TSV must contain columns: gene_symbol, keep")
}

keep_vec <- sig_df$keep
keep_mask <- if (is.logical(keep_vec)) {
  !is.na(keep_vec) & keep_vec
} else {
  toupper(as.character(keep_vec)) %in% c("TRUE", "T", "1", "YES", "Y")
}

sig_keep <- sig_df %>%
  dplyr::mutate(gene_symbol = toupper(as.character(.data$gene_symbol))) %>%
  dplyr::filter(keep_mask) %>%
  dplyr::pull(.data$gene_symbol) %>%
  unique()

if (length(sig_keep) < 10) {
  stop("Too few signature genes after curation (keep==TRUE): ", length(sig_keep))
}

sig_up <- character(0)
sig_down <- character(0)
if (file.exists(approach1_de_path)) {
  de <- readr::read_tsv(approach1_de_path, show_col_types = FALSE) %>%
    dplyr::mutate(gene_symbol = toupper(as.character(.data$gene_symbol))) %>%
    dplyr::filter(.data$gene_symbol %in% sig_keep) %>%
    dplyr::filter(!is.na(.data$log2FoldChange))

  sig_up <- de %>% dplyr::filter(.data$log2FoldChange > 0) %>% dplyr::pull(.data$gene_symbol) %>% unique()
  sig_down <- de %>% dplyr::filter(.data$log2FoldChange < 0) %>% dplyr::pull(.data$gene_symbol) %>% unique()
}

gene_sets <- list(
  menopause_all = sig_keep,
  menopause_up = sig_up,
  menopause_down = sig_down
)
gene_sets <- gene_sets[lengths(gene_sets) > 0]

## -----------------------------------------------------------------------------
## Load GTEx female subcutaneous data + VST
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("MENOPAUSE SIGNATURE SCORING (GSVA; Female-only)\n")
cat("===========================================\n\n")

sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything()) %>%
  dplyr::mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])

subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything())

meta_full <- dplyr::left_join(sample_attr, subject_pheno, by = "SUBJID")

## -----------------------------------------------------------------------------
## Scoring helpers
## -----------------------------------------------------------------------------
score_set_z <- function(expr, genes) {
  genes <- intersect(unique(toupper(genes)), rownames(expr))
  if (length(genes) < 3) {
    return(rep(NA_real_, ncol(expr)))
  }
  m <- expr[genes, , drop = FALSE]
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

run_gsva <- function(expr, sets, method) {
  expr <- as.matrix(expr)
  sets <- lapply(sets, function(g) intersect(unique(toupper(g)), rownames(expr)))
  sets <- sets[lengths(sets) >= 5]
  if (!length(sets)) stop("No gene sets remain after filtering (>=5 genes).")

  res <- tryCatch(
    GSVA::gsva(
      expr,
      sets,
      method = method,
      kcdf = "Gaussian",
      abs.ranking = FALSE,
      verbose = FALSE
    ),
    error = function(e) NULL
  )
  if (!is.null(res)) return(res)

  if (exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    param_fun <- switch(
      method,
      gsva = GSVA::gsvaParam,
      ssgsea = GSVA::ssgseaParam,
      plage = GSVA::plageParam,
      zscore = GSVA::zscoreParam,
      stop("Unsupported GSVA method: ", method)
    )
    param <- switch(
      method,
      gsva = param_fun(expr, sets, kcdf = "Gaussian", absRanking = FALSE),
      ssgsea = param_fun(expr, sets),
      plage = param_fun(expr, sets),
      zscore = param_fun(expr, sets)
    )
    return(GSVA::gsva(param, verbose = FALSE))
  }

  stop("GSVA call failed; could not find a compatible API for method: ", method)
}

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
    dplyr::filter(sex == "Female")

  counts <- gct$counts[, meta_f$SAMPID, drop = FALSE]
  counts <- collapse_dupe_rowsum(counts)
  counts <- round(counts)
  storage.mode(counts) <- "integer"

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = data.frame(row.names = colnames(counts)),
    design = ~ 1
  )
  vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = TRUE))
  expr_sym <- collapse_to_symbol_mean(vst_mat, gct$gene_id, gct$gene_symbol, uppercase = TRUE)

  ## Define menopause groups
  meta_f <- meta_f %>%
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

  ## -----------------------------------------------------------------------------
  ## Compute scores
  ## -----------------------------------------------------------------------------
  z_scores <- tibble::tibble(
    SAMPID = colnames(expr_sym),
    z_menopause_all = score_set_z(expr_sym, gene_sets$menopause_all)
  )

  if (!is.null(gene_sets$menopause_up) && !is.null(gene_sets$menopause_down)) {
    z_scores <- z_scores %>%
      dplyr::mutate(
        z_menopause_up = score_set_z(expr_sym, gene_sets$menopause_up),
        z_menopause_down = score_set_z(expr_sym, gene_sets$menopause_down),
        z_menopause_signed = .data$z_menopause_up - .data$z_menopause_down
      )
  }

  methods <- c("gsva", "ssgsea", "plage", "zscore")
  gsva_scores <- purrr::map(methods, function(m) {
    mat <- run_gsva(expr_sym, gene_sets, method = m)
    df <- as.data.frame(t(mat)) %>%
      tibble::rownames_to_column("SAMPID") %>%
      dplyr::rename_with(~ paste0(m, "__", .x), -SAMPID)
    df
  })
  gsva_scores <- Reduce(function(x, y) dplyr::left_join(x, y, by = "SAMPID"), gsva_scores)

  score_df <- meta_f %>%
    dplyr::select(SAMPID, SMCENTER, age_mid, age_bin_label, menopause_group_3) %>%
    dplyr::left_join(z_scores, by = "SAMPID") %>%
    dplyr::left_join(gsva_scores, by = "SAMPID")

  readr::write_tsv(score_df, file.path(tab_dir, paste0("menopause_signature_scores_vst", suffix, ".tsv")))

  ## -----------------------------------------------------------------------------
  ## Group stats (SMCENTER-adjusted)
  ## -----------------------------------------------------------------------------
  score_cols <- names(score_df)
  score_cols <- score_cols[grepl("^(z_|gsva__|ssgsea__|plage__|zscore__)", score_cols)]

  stats_tbl <- purrr::map_dfr(score_cols, function(col_nm) {
    df <- score_df %>% dplyr::filter(!is.na(.data$menopause_group_3), !is.na(.data[[col_nm]]))
    if (nrow(df) < 30) return(NULL)
    fit <- stats::lm(stats::as.formula(paste(col_nm, "~ SMCENTER + menopause_group_3")), data = df)
    a <- lm_contrast(fit, "menopause_group_3peri", NULL)
    b <- lm_contrast(fit, "menopause_group_3post", NULL)
    c <- lm_contrast(fit, "menopause_group_3post", "menopause_group_3peri")
    tibble::tibble(
      score = col_nm,
      n = nrow(df),
      peri_vs_pre_est = a$estimate, peri_vs_pre_p = a$p,
      post_vs_pre_est = b$estimate, post_vs_pre_p = b$p,
      post_vs_peri_est = c$estimate, post_vs_peri_p = c$p,
      n_pre = sum(df$menopause_group_3 == "pre"),
      n_peri = sum(df$menopause_group_3 == "peri"),
      n_post = sum(df$menopause_group_3 == "post")
    )
  })

  readr::write_tsv(stats_tbl, file.path(tab_dir, paste0("menopause_signature_score_stats", suffix, ".tsv")))

  ## -----------------------------------------------------------------------------
  ## One summary plot (signed z-score if available, else z_all)
  ## -----------------------------------------------------------------------------
  plot_score <- if ("z_menopause_signed" %in% names(score_df)) "z_menopause_signed" else "z_menopause_all"
  score_df %>%
    dplyr::filter(!is.na(.data$menopause_group_3), !is.na(.data[[plot_score]])) %>%
    ggplot2::ggplot(ggplot2::aes(x = menopause_group_3, y = .data[[plot_score]])) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.35, size = 1) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste0("Menopause signature score (", depot_key, "; VST)"),
      subtitle = paste0("Score: ", plot_score, " | groups by age_mid (<45, 45-55, >55)"),
      x = "Menopause proxy group",
      y = "Score"
    )
  ggplot2::ggsave(
    filename = file.path(fig_dir, paste0("menopause_signature_score__", plot_score, suffix, ".png")),
    width = 8,
    height = 5,
    dpi = 200
  )

  cat("\nWrote (", depot_key, "):\n", sep = "")
  cat("  - ", file.path(tab_dir, paste0("menopause_signature_scores_vst", suffix, ".tsv")), "\n", sep = "")
  cat("  - ", file.path(tab_dir, paste0("menopause_signature_score_stats", suffix, ".tsv")), "\n", sep = "")
  cat("  - ", file.path(fig_dir, paste0("menopause_signature_score__", plot_score, suffix, ".png")), "\n", sep = "")
}
