################################################################################
## Kenichi Manuscript Support: SNS / Adrenergic / Lipolysis Signals in GTEx
##
## Subcutaneous adipose only.
## Primary focus: women-only contrasts (with batch/center covariate).
##
## We run two complementary analyses:
## 1) 3-group model (includes transition age range):
##    pre  : age_mid < 45
##    peri : 45 <= age_mid <= 55
##    post : age_mid > 55
##
## 2) "Strict" 2-group model (excludes transition range, like Approach 1):
##    pre  : age_mid < 45
##    post : age_mid > 55
##
## Outputs:
## - DESeq2 results (all genes + panel subset) for key contrasts
## - Gene-set z-score activity plots + covariate-adjusted stats
## - Markdown report suitable for collaborator review
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
  DESeq2, apeglm
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
fig_dir <- file.path(out_dir, "figs", "kenichi_sns_panel")
tab_dir <- file.path(out_dir, "tables", "kenichi_sns_panel")
ensure_dirs(c(fig_dir, tab_dir))

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
## Curated panels (from Kenichi email + common adipose readouts)
## -----------------------------------------------------------------------------
gene_sets <- list(
  acute_beta_adrenergic = c("NR4A1", "NR4A2", "NR4A3", "FOS", "JUN", "EGR1"),
  adrenergic_receptors = c("ADRB1", "ADRB2", "ADRB3", "ADRA2A", "ADRA2B", "ADRA2C"),
  lipolysis_core = c("PNPLA2", "LIPE", "MGLL", "ABHD5", "PLIN1", "G0S2", "LPL", "FABP4"),
  thermogenesis_program = c("UCP1", "PPARGC1A", "PPARGC1B", "CIDEA", "DIO2", "VEGFA"),
  myers2009_context = c(
    "ACVR1", "AKT2", "CBL", "CEBPA", "FASN", "GDF10", "GYS1",
    "HDAC5", "HDAC7", "IFNG", "IL15", "IRS1", "LEP", "PPARGC1B",
    "PYGM", "RAB18", "SLC2A2", "SLC2A4", "SLC2A8", "SMAD2", "STAT5B", "UCP3"
  )
)
gene_sets <- purrr::map(gene_sets, ~ unique(toupper(.x)))

## Optional: merge in extra gene sets (e.g., Christoph/Kenichi curated lists)
## Format: TSV with columns: set_name<TAB>gene_symbol (header required).
extra_gene_sets_path <- file.path("references", "kenichi_extra_gene_sets.tsv")
if (file.exists(extra_gene_sets_path)) {
  extra_df <- readr::read_tsv(extra_gene_sets_path, show_col_types = FALSE, comment = "#") %>%
    dplyr::filter(!is.na(.data$set_name), !is.na(.data$gene_symbol)) %>%
    dplyr::mutate(
      set_name = as.character(.data$set_name),
      gene_symbol = toupper(as.character(.data$gene_symbol))
    )
  if (nrow(extra_df)) {
    extra_sets <- extra_df %>%
      dplyr::group_by(.data$set_name) %>%
      dplyr::summarise(genes = list(unique(.data$gene_symbol)), .groups = "drop")
    for (i in seq_len(nrow(extra_sets))) {
      nm <- extra_sets$set_name[[i]]
      genes <- extra_sets$genes[[i]]
      if (nm %in% names(gene_sets)) {
        gene_sets[[nm]] <- unique(c(gene_sets[[nm]], genes))
      } else {
        gene_sets[[nm]] <- unique(genes)
      }
    }
    message("Loaded extra gene sets from: ", extra_gene_sets_path)
  }
}

panel_symbols <- unique(unlist(gene_sets, use.names = FALSE))

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------
lfc_shrink_safe <- function(dds, coef_name, res_obj) {
  tryCatch(
    DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm"),
    error = function(e) {
      message("Note: lfcShrink failed for ", coef_name, ", using unshrunk results. Reason: ", e$message)
      res_obj
    }
  )
}

res_to_df <- function(res_obj, gene_symbol_lookup, contrast_label, model_label) {
  as.data.frame(res_obj) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::mutate(
      gene_symbol = gene_symbol_lookup[match(.data$gene_id, names(gene_symbol_lookup))],
      model = model_label,
      contrast = contrast_label
    ) %>%
    dplyr::relocate(model, contrast, gene_id, gene_symbol)
}

subset_panel <- function(df, panel_syms) {
  df %>%
    dplyr::mutate(gene_symbol_upper = toupper(.data$gene_symbol)) %>%
    dplyr::filter(!is.na(.data$gene_symbol_upper), .data$gene_symbol_upper %in% panel_syms) %>%
    dplyr::select(-gene_symbol_upper) %>%
    dplyr::arrange(.data$padj)
}

make_symbol_expr_vst <- function(vst_mat, gene_id, gene_symbol) {
  sym <- toupper(gene_symbol[match(rownames(vst_mat), gene_id)])
  keep <- !is.na(sym) & nzchar(sym)
  vst_mat <- vst_mat[keep, , drop = FALSE]
  sym <- sym[keep]
  if (any(duplicated(sym))) {
    summed <- rowsum(vst_mat, group = sym, reorder = FALSE)
    denom <- as.numeric(table(sym)[rownames(summed)])
    vst_mat <- summed / denom
    rownames(vst_mat) <- rownames(summed)
  } else {
    rownames(vst_mat) <- sym
  }
  vst_mat
}

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
  ## Compare two factor levels a vs b in a model with treatment contrasts.
  ## Returns (estimate, se, t, p).
  cf <- stats::coef(fit)
  vc <- stats::vcov(fit)
  cn <- names(cf)
  # baseline is the first level; treatment terms are e.g. menopause_groupperi, menopause_grouppost
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
## Load GTEx female subcutaneous data
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("KENICHI SNS PANEL (GTEx v10, Female-only)\n")
cat("===========================================\n\n")

cat("Loading metadata...\n")
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

  cat("Loading counts: ", dpt$gct_path, "\n", sep = "")
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
    derive_age_numeric()

  meta_f <- meta %>%
    dplyr::filter(sex == "Female") %>%
    dplyr::mutate(SMCENTER = factor(SMCENTER))

  counts <- gct$counts[, meta_f$SAMPID, drop = FALSE]
  counts <- collapse_dupe_rowsum(counts)

  gene_symbol_lookup <- gct$gene_symbol
  names(gene_symbol_lookup) <- gct$gene_id

  cat(sprintf("Female %s samples: %d\n", depot_key, nrow(meta_f)))
  cat("Age-bin counts (Female):\n")
  print(table(meta_f$age_bin_label, useNA = "ifany"))

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

  meta_strict <- meta_f %>%
    dplyr::filter(!is.na(menopause_group_3), menopause_group_3 %in% c("pre", "post")) %>%
    dplyr::mutate(menopause_group_2 = factor(as.character(menopause_group_3), levels = c("pre", "post")))

  cat("\nMenopause groups (3-level, Female):\n")
  print(table(meta_f$menopause_group_3, useNA = "ifany"))
  cat("\nMenopause groups (strict 2-level, Female):\n")
  print(table(meta_strict$menopause_group_2, useNA = "ifany"))

## -----------------------------------------------------------------------------
## DESeq2: 3-group model (includes peri)
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("DESeq2: 3-GROUP MODEL (pre/peri/post)\n")
cat("design = ~ SMCENTER + menopause_group_3\n")
cat("===========================================\n\n")

meta_3 <- meta_f %>%
  dplyr::filter(!is.na(menopause_group_3)) %>%
  dplyr::select(SAMPID, SMCENTER, menopause_group_3, age_mid, age_bin_label)

counts_3 <- counts[, meta_3$SAMPID, drop = FALSE]

dds3 <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_3,
  colData = meta_3,
  design = ~ SMCENTER + menopause_group_3
)
dds3 <- dds3[rowSums(DESeq2::counts(dds3)) >= 10, ]
dds3 <- DESeq2::DESeq(dds3)

rn3 <- DESeq2::resultsNames(dds3)
cat("resultsNames(dds3):\n")
print(rn3)

coef_peri <- grep("^menopause_group_3_peri_vs_pre$", rn3, value = TRUE)
coef_post <- grep("^menopause_group_3_post_vs_pre$", rn3, value = TRUE)
if (length(coef_peri) != 1 || length(coef_post) != 1) {
  stop("Unexpected menopause_group_3 coefficient names: ", paste(rn3, collapse = ", "))
}

res3_post_vs_pre <- DESeq2::results(dds3, name = coef_post[[1]])
res3_post_vs_pre <- lfc_shrink_safe(dds3, coef_post[[1]], res3_post_vs_pre)

res3_peri_vs_pre <- DESeq2::results(dds3, name = coef_peri[[1]])
res3_peri_vs_pre <- lfc_shrink_safe(dds3, coef_peri[[1]], res3_peri_vs_pre)

## post vs peri = (post vs pre) - (peri vs pre)
res3_post_vs_peri <- DESeq2::results(dds3, contrast = c("menopause_group_3", "post", "peri"))
## Note: apeglm shrinkage in DESeq2 is coef-based; for this non-baseline contrast we
## keep the unshrunk estimate and use shrunk estimates for the baseline contrasts.

df3_post_vs_pre <- res_to_df(res3_post_vs_pre, gene_symbol_lookup, "post_vs_pre", "3group") %>% dplyr::arrange(.data$padj)
df3_peri_vs_pre <- res_to_df(res3_peri_vs_pre, gene_symbol_lookup, "peri_vs_pre", "3group") %>% dplyr::arrange(.data$padj)
df3_post_vs_peri <- res_to_df(res3_post_vs_peri, gene_symbol_lookup, "post_vs_peri", "3group") %>% dplyr::arrange(.data$padj)

readr::write_tsv(df3_post_vs_pre, file.path(tab_dir, paste0("DE_3group_post_vs_pre_all_genes", suffix, ".tsv")))
readr::write_tsv(df3_peri_vs_pre, file.path(tab_dir, paste0("DE_3group_peri_vs_pre_all_genes", suffix, ".tsv")))
readr::write_tsv(df3_post_vs_peri, file.path(tab_dir, paste0("DE_3group_post_vs_peri_all_genes", suffix, ".tsv")))

readr::write_tsv(subset_panel(df3_post_vs_pre, panel_symbols), file.path(tab_dir, paste0("DE_3group_post_vs_pre_PANEL", suffix, ".tsv")))
readr::write_tsv(subset_panel(df3_peri_vs_pre, panel_symbols), file.path(tab_dir, paste0("DE_3group_peri_vs_pre_PANEL", suffix, ".tsv")))
readr::write_tsv(subset_panel(df3_post_vs_peri, panel_symbols), file.path(tab_dir, paste0("DE_3group_post_vs_peri_PANEL", suffix, ".tsv")))

## -----------------------------------------------------------------------------
## DESeq2: strict 2-group model (excludes peri)
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("DESeq2: STRICT 2-GROUP MODEL (pre vs post)\n")
cat("design = ~ SMCENTER + menopause_group_2\n")
cat("===========================================\n\n")

meta_2 <- meta_strict %>%
  dplyr::select(SAMPID, SMCENTER, menopause_group_2, age_mid, age_bin_label)

counts_2 <- counts[, meta_2$SAMPID, drop = FALSE]

dds2 <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_2,
  colData = meta_2,
  design = ~ SMCENTER + menopause_group_2
)
dds2 <- dds2[rowSums(DESeq2::counts(dds2)) >= 10, ]
dds2 <- DESeq2::DESeq(dds2)

res2_post_vs_pre <- DESeq2::results(dds2, contrast = c("menopause_group_2", "post", "pre"))
rn2 <- DESeq2::resultsNames(dds2)
coef2 <- grep("^menopause_group_2_post_vs_pre$", rn2, value = TRUE)
if (length(coef2) != 1) {
  message("Note: unexpected coef naming for menopause_group_2, skipping shrink. resultsNames(dds2) = ", paste(rn2, collapse = ", "))
} else {
  res2_post_vs_pre <- lfc_shrink_safe(dds2, coef2[[1]], res2_post_vs_pre)
}

df2_post_vs_pre <- res_to_df(res2_post_vs_pre, gene_symbol_lookup, "post_vs_pre", "strict2") %>%
  dplyr::arrange(.data$padj)

readr::write_tsv(df2_post_vs_pre, file.path(tab_dir, paste0("DE_strict2_post_vs_pre_all_genes", suffix, ".tsv")))
readr::write_tsv(subset_panel(df2_post_vs_pre, panel_symbols), file.path(tab_dir, paste0("DE_strict2_post_vs_pre_PANEL", suffix, ".tsv")))

## -----------------------------------------------------------------------------
## DESeq2: explicit menopause-adjacent bins (40-49 vs 50-59) in women
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("DESeq2: AGE BIN 40-49 vs 50-59 (Female only)\n")
cat("design = ~ SMCENTER + age_group\n")
cat("===========================================\n\n")

meta_bin <- meta_f %>%
  dplyr::filter(age_bin_label %in% c("40-49", "50-59")) %>%
  dplyr::mutate(
    age_group = factor(age_bin_label, levels = c("40-49", "50-59")),
    SMCENTER = factor(SMCENTER)
  ) %>%
  dplyr::select(SAMPID, SMCENTER, age_group, age_mid, age_bin_label)

cat(sprintf("Samples: 40-49=%d, 50-59=%d\n",
            sum(meta_bin$age_group == "40-49"), sum(meta_bin$age_group == "50-59")))

counts_bin <- counts[, meta_bin$SAMPID, drop = FALSE]

dds_bin <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts_bin,
  colData = meta_bin,
  design = ~ SMCENTER + age_group
)
dds_bin <- dds_bin[rowSums(DESeq2::counts(dds_bin)) >= 10, ]
dds_bin <- DESeq2::DESeq(dds_bin)

res_bin <- DESeq2::results(dds_bin, contrast = c("age_group", "50-59", "40-49"))
rn_bin <- DESeq2::resultsNames(dds_bin)
coef_bin <- grep("^age_group", rn_bin, value = TRUE)
if (length(coef_bin) == 1) {
  res_bin <- lfc_shrink_safe(dds_bin, coef_bin[[1]], res_bin)
} else {
  message("Note: unexpected age_group coefficient naming, skipping shrink. resultsNames(dds_bin) = ", paste(rn_bin, collapse = ", "))
}

df_bin <- res_to_df(res_bin, gene_symbol_lookup, "50-59_vs_40-49", "bins_40-49_50-59") %>%
  dplyr::arrange(.data$padj)

readr::write_tsv(df_bin, file.path(tab_dir, paste0("DE_bins_40-49_vs_50-59_all_genes", suffix, ".tsv")))
readr::write_tsv(subset_panel(df_bin, panel_symbols), file.path(tab_dir, paste0("DE_bins_40-49_vs_50-59_PANEL", suffix, ".tsv")))

## -----------------------------------------------------------------------------
## Gene-set activity (z-score across genes; VST expression)
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("Gene-set activity (VST z-score)\n")
cat("===========================================\n\n")

vst3 <- DESeq2::vst(dds3, blind = FALSE)
expr_vst <- SummarizedExperiment::assay(vst3)
expr_sym <- make_symbol_expr_vst(expr_vst, gct$gene_id, gct$gene_symbol)

score_tbl <- meta_3 %>%
  dplyr::mutate(SAMPID = as.character(SAMPID)) %>%
  dplyr::select(SAMPID, SMCENTER, menopause_group_3, age_mid, age_bin_label)

for (nm in names(gene_sets)) {
  sc <- score_set_z(expr_sym, gene_sets[[nm]])
  score_tbl[[paste0("score__", nm)]] <- sc[match(score_tbl$SAMPID, colnames(expr_sym))]
}

readr::write_tsv(score_tbl, file.path(tab_dir, paste0("gene_set_scores_vst_z", suffix, ".tsv")))

## Stats + plots
stats_rows <- list()
plot_list <- list()

for (nm in names(gene_sets)) {
  col <- paste0("score__", nm)
  df <- score_tbl %>%
    dplyr::filter(!is.na(.data[[col]]), !is.na(.data$menopause_group_3)) %>%
    dplyr::mutate(menopause_group_3 = droplevels(menopause_group_3))

  fit <- stats::lm(stats::as.formula(paste0(col, " ~ SMCENTER + menopause_group_3")), data = df)
  cf <- summary(fit)$coefficients

  ## Baseline = pre; terms are menopause_group_3peri and menopause_group_3post
  term_peri <- "menopause_group_3peri"
  term_post <- "menopause_group_3post"
  peri_vs_pre <- lm_contrast(fit, term_peri, "")
  post_vs_pre <- lm_contrast(fit, term_post, "")
  post_vs_peri <- lm_contrast(fit, term_post, term_peri)

  means <- df %>%
    dplyr::group_by(menopause_group_3) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(.data[[col]], na.rm = TRUE),
      sd = sd(.data[[col]], na.rm = TRUE),
      .groups = "drop"
    )

  stats_rows[[nm]] <- tibble::tibble(
    gene_set = nm,
    n = nrow(df),
    peri_vs_pre_est = peri_vs_pre$estimate, peri_vs_pre_p = peri_vs_pre$p,
    post_vs_pre_est = post_vs_pre$estimate, post_vs_pre_p = post_vs_pre$p,
    post_vs_peri_est = post_vs_peri$estimate, post_vs_peri_p = post_vs_peri$p
  )
  means_wide <- means %>% tidyr::pivot_wider(names_from = menopause_group_3, values_from = c(n, mean, sd))
  stats_rows[[nm]] <- dplyr::bind_cols(stats_rows[[nm]], means_wide)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = menopause_group_3, y = .data[[col]], color = menopause_group_3)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.25) +
    ggplot2::geom_jitter(width = 0.18, alpha = 0.35, size = 0.9) +
    ggplot2::theme_bw() +
    ggplot2::scale_color_manual(values = c(pre = "#1b9e77", peri = "#d95f02", post = "#7570b3")) +
    ggplot2::labs(
      title = glue::glue("{nm} (VST z-score activity)"),
      subtitle = glue::glue("LM adj: peri vs pre p={signif(peri_vs_pre$p, 3)}; post vs pre p={signif(post_vs_pre$p, 3)}; post vs peri p={signif(post_vs_peri$p, 3)}"),
      x = NULL,
      y = "Mean z-score across genes"
    ) +
    ggplot2::theme(legend.position = "none")

  plot_list[[nm]] <- p
  ggplot2::ggsave(file.path(fig_dir, glue::glue("gene_set_activity__{nm}{suffix}.png")), p, width = 7.2, height = 4.2, dpi = 300)
}

stats_tbl <- dplyr::bind_rows(stats_rows) %>% dplyr::arrange(.data$gene_set)
readr::write_tsv(stats_tbl, file.path(tab_dir, paste0("gene_set_activity_stats", suffix, ".tsv")))

## Additional: gene-set activity for menopause-adjacent bins only (40-49 vs 50-59)
bin_stats <- list()
for (nm in names(gene_sets)) {
  col <- paste0("score__", nm)
  df <- score_tbl %>%
    dplyr::filter(age_bin_label %in% c("40-49", "50-59")) %>%
    dplyr::mutate(age_group = factor(age_bin_label, levels = c("40-49", "50-59"))) %>%
    dplyr::filter(!is.na(.data[[col]]), !is.na(.data$age_group))

  fit <- stats::lm(stats::as.formula(paste0(col, " ~ SMCENTER + age_group")), data = df)
  # treatment contrast: age_group50-59 is 50-59 vs 40-49
  term <- "age_group50-59"
  ctr <- lm_contrast(fit, term, "")
  bin_stats[[nm]] <- tibble::tibble(
    gene_set = nm,
    n = nrow(df),
    est_50_59_vs_40_49 = ctr$estimate,
    p_50_59_vs_40_49 = ctr$p,
    mean_40_49 = mean(df[[col]][df$age_group == "40-49"], na.rm = TRUE),
    mean_50_59 = mean(df[[col]][df$age_group == "50-59"], na.rm = TRUE)
  )
}
bin_stats_tbl <- dplyr::bind_rows(bin_stats) %>% dplyr::arrange(.data$gene_set)
readr::write_tsv(bin_stats_tbl, file.path(tab_dir, paste0("gene_set_activity_stats_40-49_vs_50-59", suffix, ".tsv")))

## Candidate "manuscript-support" figure (stacked panels)
if (all(c("lipolysis_core", "thermogenesis_program") %in% names(plot_list))) {
  fig_main <- plot_list[["lipolysis_core"]] / plot_list[["thermogenesis_program"]]
  ggplot2::ggsave(
    filename = file.path(fig_dir, paste0("figure_candidate__lipolysis_thermogenesis", suffix, ".png")),
    plot = fig_main,
    width = 7.2,
    height = 8.2,
    dpi = 300
  )
}
if ("acute_beta_adrenergic" %in% names(plot_list)) {
  ggplot2::ggsave(
    filename = file.path(fig_dir, paste0("figure_candidate__acute_beta_adrenergic", suffix, ".png")),
    plot = plot_list[["acute_beta_adrenergic"]],
    width = 7.2,
    height = 4.2,
    dpi = 300
  )
}

## -----------------------------------------------------------------------------
## Write a short markdown report
## -----------------------------------------------------------------------------
report_path <- file.path(tab_dir, paste0("kenichi_sns_panel_report", suffix, ".md"))

fmt_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 1e-3) return(format(x, scientific = TRUE, digits = 2))
  format(round(x, 4), nsmall = 4, trim = TRUE)
}

highlight_sets <- c("acute_beta_adrenergic", "lipolysis_core", "thermogenesis_program", "adrenergic_receptors")
stats_small <- stats_tbl %>%
  dplyr::filter(.data$gene_set %in% highlight_sets) %>%
  dplyr::mutate(
    peri_vs_pre_p_fmt = vapply(.data$peri_vs_pre_p, fmt_p, character(1)),
    post_vs_pre_p_fmt = vapply(.data$post_vs_pre_p, fmt_p, character(1)),
    post_vs_peri_p_fmt = vapply(.data$post_vs_peri_p, fmt_p, character(1))
  )

top_strict_panel <- subset_panel(df2_post_vs_pre, panel_symbols) %>%
  dplyr::filter(!is.na(.data$padj)) %>%
  dplyr::arrange(.data$padj) %>%
  dplyr::slice_head(n = 8) %>%
  dplyr::mutate(
    padj_fmt = vapply(.data$padj, fmt_p, character(1)),
    line = glue::glue("{gene_symbol}: log2FC={round(log2FoldChange, 3)}, padj={padj_fmt}")
  ) %>%
  dplyr::pull(.data$line)

readme <- c(
  glue::glue("# Kenichi SNS/Lipolysis Panel (GTEx v10, {depot_label}, Women-only)"),
  "",
  "This report provides targeted, hypothesis-driven readouts requested for the estrogen deficiency project.",
  "",
  "## Key caveats (from collaborator emails)",
  "- GTEx adipose is baseline/unstimulated; acute SNS-responsive transcripts may be weak or transient.",
  "- Chronic SNS stimulation can reduce ADRB3 (catecholamine resistance); directionality is not a simple proxy for SNS tone.",
  "- \"Men as a negative control\" is imperfect (men have hormonal aging and psychosocial confounds); male comparisons are supportive only.",
  "",
  "## Definitions",
  glue::glue("- Tissue: {depot_label} only."),
  "- Cohort: Female samples only.",
  "- 3-group model (includes transition range): pre (<45), peri (45-55), post (>55) using age-bin midpoints.",
  "- Strict 2-group model: pre (<45) vs post (>55), excluding peri.",
  "",
  "## Snapshot (quick read)",
  "### Gene-set activity (adjusted for SMCENTER; VST z-score)",
  paste0(
    "- ", stats_small$gene_set,
    " | peri vs pre p=", stats_small$peri_vs_pre_p_fmt,
    " | post vs pre p=", stats_small$post_vs_pre_p_fmt,
    " | post vs peri p=", stats_small$post_vs_peri_p_fmt
  ),
  "",
  "### Strict pre vs post (panel genes; top hits by padj)",
  if (length(top_strict_panel)) paste0("- ", top_strict_panel) else "- (no panel genes available)",
  "",
  "## Outputs",
  "- DESeq2 (3-group):",
  glue::glue("  - All genes: `DE_3group_*_all_genes{suffix}.tsv`"),
  glue::glue("  - Panel only: `DE_3group_*_PANEL{suffix}.tsv`"),
  "- DESeq2 (strict 2-group):",
  glue::glue("  - All genes: `DE_strict2_post_vs_pre_all_genes{suffix}.tsv`"),
  glue::glue("  - Panel only: `DE_strict2_post_vs_pre_PANEL{suffix}.tsv`"),
  "- DESeq2 (menopause-adjacent bins):",
  glue::glue("  - All genes: `DE_bins_40-49_vs_50-59_all_genes{suffix}.tsv`"),
  glue::glue("  - Panel only: `DE_bins_40-49_vs_50-59_PANEL{suffix}.tsv`"),
  "- Gene-set activity (VST z-score):",
  glue::glue("  - Per-sample scores: `gene_set_scores_vst_z{suffix}.tsv`"),
  glue::glue("  - Adjusted stats: `gene_set_activity_stats{suffix}.tsv`"),
  glue::glue("  - 40-49 vs 50-59 stats: `gene_set_activity_stats_40-49_vs_50-59{suffix}.tsv`"),
  "  - Plots: see figs in `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/`",
  "",
  "## Gene sets",
  paste0("- ", names(gene_sets), collapse = "\n"),
  "",
  "Generated: ",
  as.character(Sys.time())
)

writeLines(readme, report_path)
cat("Wrote report: ", report_path, "\n", sep = "")

## -----------------------------------------------------------------------------
## Draft message to Kenichi (editable)
## -----------------------------------------------------------------------------
draft_path <- file.path(tab_dir, paste0("kenichi_draft_message", suffix, ".md"))

draft <- c(
  "# Draft: Message to Kenichi (GTEx human adipose support)",
  "",
  "Tissue / cohort / model:",
  glue::glue("- GTEx v10 {depot_label}, Female only."),
  "- Adjusted for collection center (SMCENTER) in all models.",
  "- Menopause proxy groups by age-bin midpoint:",
  "  - pre: <45",
  "  - peri (transition): 45-55",
  "  - post: >55",
  "",
  "## Main points (what we can and cannot claim)",
  "- We do **not** see robust upregulation of acute beta-adrenergic immediate-early markers at baseline (NR4A1/2/3, FOS/JUN/EGR1) when comparing post vs pre; any shift is small and not statistically strong at the gene-set level.",
  "- Program-level gene-set shifts (SMCENTER-adjusted; VST z-score) are summarized below (estimates and p-values are from linear models):",
  glue::glue("  - lipolysis_core: peri vs pre est={round(stats_tbl$peri_vs_pre_est[stats_tbl$gene_set==\"lipolysis_core\"],3)} p={fmt_p(stats_tbl$peri_vs_pre_p[stats_tbl$gene_set==\"lipolysis_core\"])}; post vs pre est={round(stats_tbl$post_vs_pre_est[stats_tbl$gene_set==\"lipolysis_core\"],3)} p={fmt_p(stats_tbl$post_vs_pre_p[stats_tbl$gene_set==\"lipolysis_core\"])}."),
  glue::glue("  - thermogenesis_program: peri vs pre est={round(stats_tbl$peri_vs_pre_est[stats_tbl$gene_set==\"thermogenesis_program\"],3)} p={fmt_p(stats_tbl$peri_vs_pre_p[stats_tbl$gene_set==\"thermogenesis_program\"])}; post vs pre est={round(stats_tbl$post_vs_pre_est[stats_tbl$gene_set==\"thermogenesis_program\"],3)} p={fmt_p(stats_tbl$post_vs_pre_p[stats_tbl$gene_set==\"thermogenesis_program\"])}."),
  glue::glue("  - acute_beta_adrenergic: peri vs pre est={round(stats_tbl$peri_vs_pre_est[stats_tbl$gene_set==\"acute_beta_adrenergic\"],3)} p={fmt_p(stats_tbl$peri_vs_pre_p[stats_tbl$gene_set==\"acute_beta_adrenergic\"])}; post vs pre est={round(stats_tbl$post_vs_pre_est[stats_tbl$gene_set==\"acute_beta_adrenergic\"],3)} p={fmt_p(stats_tbl$post_vs_pre_p[stats_tbl$gene_set==\"acute_beta_adrenergic\"])}."),
  "",
  "Interpretation suggestion (aligned with your caveats):",
  "- In baseline human adipose, the data are more consistent with **program-level metabolic remodeling** across menopause proxy groups, rather than clear evidence for **acute SNS activation**.",
  "- This does not rule out increased sympathetic outflow; chronic stimulation can produce catecholamine resistance and does not necessarily yield higher expression of adrenergic receptors or acute-response transcripts in unstimulated tissue.",
  "",
  "## Files to review",
  glue::glue("- Gene-set activity stats (pre/peri/post): `gene_set_activity_stats{suffix}.tsv`"),
  glue::glue("- Gene-set activity stats (40-49 vs 50-59): `gene_set_activity_stats_40-49_vs_50-59{suffix}.tsv`"),
  glue::glue("- Panel DE (strict pre vs post): `DE_strict2_post_vs_pre_PANEL{suffix}.tsv`"),
  "- Candidate figure panels:",
  glue::glue("  - `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__lipolysis_thermogenesis{suffix}.png`"),
  glue::glue("  - `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__acute_beta_adrenergic{suffix}.png`"),
  "",
  "## Next refinement",
  "- If you share Christoph's exact SNS/adrenergic/lipolysis gene lists, we can add them via `references/kenichi_extra_gene_sets.tsv` and regenerate the same tables/figures.",
  "",
  paste0("Generated: ", Sys.time())
)

writeLines(draft, draft_path)
cat("Wrote draft message: ", draft_path, "\n", sep = "")

} # end depot loop
