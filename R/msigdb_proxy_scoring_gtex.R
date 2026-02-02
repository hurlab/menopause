################################################################################
## Score MSigDB estrogen/menopause-proxy gene sets in GTEx adipose (female)
##
## Purpose:
## - If MSigDB lacks a direct "menopause" signature, estrogen-response gene sets
##   can be used as a biologically motivated proxy axis to evaluate alongside
##   data-driven menopause signatures.
##
## Inputs:
## - references/msigdb_menopause_proxy_set_names.txt
## - GTExDatav10/ (subcutaneous adipose + annotations)
##
## Outputs:
## - GTEx_v10_AT_analysis_out/tables/msigdb_proxy_scoring/
## - GTEx_v10_AT_analysis_out/figs/msigdb_proxy_scoring/
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
  GSVA,
  msigdbr
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
fig_dir <- file.path(out_dir, "figs", "msigdb_proxy_scoring")
tab_dir <- file.path(out_dir, "tables", "msigdb_proxy_scoring")
ensure_dirs(c(fig_dir, tab_dir))

set_names_path <- file.path("references", "msigdb_menopause_proxy_set_names.txt")
if (!file.exists(set_names_path)) stop("Missing: ", set_names_path)
set_names <- readr::read_lines(set_names_path) %>%
  stringr::str_trim() %>%
  discard(~ !nzchar(.x)) %>%
  unique()

## -----------------------------------------------------------------------------
## MSigDB gene sets
## -----------------------------------------------------------------------------
msig <- msigdbr::msigdbr(species = "Homo sapiens") %>%
  dplyr::filter(.data$gs_name %in% set_names)

missing_sets <- setdiff(set_names, unique(msig$gs_name))
if (length(missing_sets)) {
  warning("Some requested set names not found in msigdbr: ", paste(missing_sets, collapse = ", "))
}

gene_sets <- msig %>%
  dplyr::group_by(.data$gs_name) %>%
  dplyr::summarise(genes = list(unique(toupper(.data$gene_symbol))), .groups = "drop")

gene_sets <- stats::setNames(gene_sets$genes, gene_sets$gs_name)

if (!length(gene_sets)) stop("No MSigDB gene sets loaded; check set names.")

## -----------------------------------------------------------------------------
## GTEx female adipose VST (subcutaneous + visceral)
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("MSIGDB PROXY SCORING (GTEx v10, Female-only)\n")
cat("===========================================\n\n")

sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything()) %>%
  dplyr::mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])

subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything())

meta_full <- dplyr::left_join(sample_attr, subject_pheno, by = "SUBJID")

depots <- list(
  list(key = "subcutaneous", label = "Adipose - Subcutaneous", gct_path = gct_sc, suffix = ""),
  list(key = "visceral", label = "Adipose - Visceral (Omentum)", gct_path = gct_vis, suffix = "_visceral")
)

## -----------------------------------------------------------------------------
## GSVA scoring + stats
## -----------------------------------------------------------------------------
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
    GSVA::gsva(expr, sets, method = method, kcdf = "Gaussian", abs.ranking = FALSE, verbose = FALSE),
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
  stop("GSVA call failed for method: ", method)
}

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

  meta_f <- meta %>% dplyr::filter(sex == "Female")

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

  ## Menopause proxy groups (same as Kenichi panel script)
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

  gsva_mat <- run_gsva(expr_sym, gene_sets, method = "ssgsea")
  score_df <- as.data.frame(t(gsva_mat)) %>%
    tibble::rownames_to_column("SAMPID") %>%
    dplyr::left_join(meta_f %>% dplyr::select(SAMPID, SMCENTER, age_mid, age_bin_label, menopause_group_3), by = "SAMPID")

  readr::write_tsv(score_df, file.path(tab_dir, paste0("msigdb_proxy_scores_ssgsea_vst", suffix, ".tsv")))

  stats_tbl <- purrr::map_dfr(names(gene_sets), function(gs) {
    df <- score_df %>% dplyr::filter(!is.na(.data$menopause_group_3), !is.na(.data[[gs]]))
    fit <- stats::lm(stats::as.formula(paste0("`", gs, "` ~ SMCENTER + menopause_group_3")), data = df)
    a <- lm_contrast(fit, "menopause_group_3peri", NULL)
    b <- lm_contrast(fit, "menopause_group_3post", NULL)
    c <- lm_contrast(fit, "menopause_group_3post", "menopause_group_3peri")
    tibble::tibble(
      gene_set = gs,
      n = nrow(df),
      peri_vs_pre_est = a$estimate, peri_vs_pre_p = a$p,
      post_vs_pre_est = b$estimate, post_vs_pre_p = b$p,
      post_vs_peri_est = c$estimate, post_vs_peri_p = c$p
    )
  })

  readr::write_tsv(stats_tbl, file.path(tab_dir, paste0("msigdb_proxy_score_stats_ssgsea_vst", suffix, ".tsv")))

  ## Simple plot for hallmark sets if present
  to_plot <- intersect(
    c("HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_ESTROGEN_RESPONSE_LATE"),
    names(gene_sets)
  )
  if (length(to_plot)) {
    df_long <- score_df %>%
      dplyr::select(SAMPID, menopause_group_3, dplyr::all_of(to_plot)) %>%
      tidyr::pivot_longer(cols = dplyr::all_of(to_plot), names_to = "gene_set", values_to = "score") %>%
      dplyr::filter(!is.na(.data$menopause_group_3))

    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = menopause_group_3, y = score)) +
      ggplot2::geom_boxplot(outlier.shape = NA) +
      ggplot2::geom_jitter(width = 0.15, alpha = 0.25, size = 0.9) +
      ggplot2::facet_wrap(~ gene_set, scales = "free_y") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste0("MSigDB estrogen-response proxy scores (", depot_key, "; ssGSEA; VST)"),
        x = "Menopause proxy group",
        y = "Score"
      )
    ggplot2::ggsave(file.path(fig_dir, paste0("msigdb_estrogen_hallmarks_ssgsea", suffix, ".png")), p, width = 10, height = 4.8, dpi = 200)
  }

  cat("Wrote (", depot_key, "):\n", sep = "")
  cat("  - ", file.path(tab_dir, paste0("msigdb_proxy_scores_ssgsea_vst", suffix, ".tsv")), "\n", sep = "")
  cat("  - ", file.path(tab_dir, paste0("msigdb_proxy_score_stats_ssgsea_vst", suffix, ".tsv")), "\n", sep = "")
  cat("  - ", fig_dir, "\n", sep = "")
}
