################################################################################
## Next steps (2026-03-04): Kenichi follow-up gene sets (inflammation/fibrosis)
## and literature-linked catecholamine resistance mechanisms.
##
## Goal:
## - Score Kenichi-provided gene sets (inflammation, fibrosis, senescence, DNL).
## - Test literature-linked catecholamine resistance mechanism axes (TGFb/ALK7,
##   TNF->NFkB->IKKε/TBK1->PDE3B).
## - Provide plots stratified by depot / sex / age bin.
## - Provide correlations to ADR receptors and CIBERSORT (noMALAT1) fractions.
##
## Inputs:
## - GTEx v10 adipose GCTs
## - references/kenichi_gene_sets_immune_fibrosis_2026-03-04.tsv
## - CIBERSORT fractions tables (noMALAT1)
##
## Output:
## - GTEx_v10_AT_analysis_out/{tables,figs}/next_steps_2026-03-04_kenichi_followup/
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
base_key <- "next_steps_2026-03-04_kenichi_followup"
tab_dir <- file.path(out_dir, "tables", base_key)
fig_dir <- file.path(out_dir, "figs", base_key)
log_dir <- file.path(out_dir, "logs")
ensure_dirs(c(tab_dir, fig_dir, log_dir))

log_path <- file.path(log_dir, paste0(base_key, ".log"))
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

age_levels <- c("20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")

depots <- list(
  list(key = "subcutaneous", label = "Adipose - Subcutaneous", gct_path = gct_sc, suffix = ""),
  list(key = "visceral", label = "Adipose - Visceral (Omentum)", gct_path = gct_vis, suffix = "_visceral")
)

# Kenichi gene sets + review-derived axes
sets_df <- readr::read_tsv(
  "references/kenichi_gene_sets_immune_fibrosis_2026-03-04.tsv",
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    set_name = as.character(set_name),
    gene_symbol = toupper(as.character(gene_symbol))
  )

gene_sets <- sets_df %>%
  dplyr::group_by(set_name) %>%
  dplyr::summarise(genes = list(unique(gene_symbol)), n_genes = dplyr::n_distinct(gene_symbol), .groups = "drop")

readr::write_tsv(gene_sets %>% dplyr::mutate(genes = purrr::map_chr(genes, ~ paste(.x, collapse = ","))),
                 file.path(tab_dir, "kenichi_and_review_gene_sets.tsv"))

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

add_scores <- function(tbl, expr_sym, sets_tbl) {
  for (i in seq_len(nrow(sets_tbl))) {
    nm <- sets_tbl$set_name[i]
    tbl[[nm]] <- score_set_z(expr_sym, sets_tbl$genes[[i]])
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

read_cibersort <- function(depot_key) {
  base <- file.path(
    out_dir,
    "tables",
    "next_steps_2026-02-28_cibersort_fractions",
    paste0("cibersort_", depot_key, "_perm0_topN50_minCells200_noMALAT1__fractions_with_meta.tsv.gz")
  )
  if (!file.exists(base)) {
    warning("Missing CIBERSORT fractions file: ", base)
    return(NULL)
  }
  df <- readr::read_tsv(base, show_col_types = FALSE)
  df
}

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
    dplyr::bind_cols(parse_age_bin(.$AGE)) %>%
    derive_age_numeric() %>%
    dplyr::mutate(
      sex = factor(sex, levels = c("Male", "Female")),
      SMCENTER = factor(SMCENTER),
      age_bin_label = factor(as.character(age_bin_label), levels = age_levels),
      depot = depot_key
    ) %>%
    dplyr::filter(!is.na(age_bin_label))

  readr::write_tsv(
    meta %>% dplyr::count(sex, age_bin_label) %>% dplyr::arrange(sex, age_bin_label),
    file.path(tab_dir, paste0("sample_counts_by_sex_agebin__", depot_key, ".tsv"))
  )

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
  expr_sym <- collapse_to_symbol_mean(vst_mat, gct$gene_id, gct$gene_symbol, uppercase = TRUE)

  # Build score table
  score_tbl <- tibble::tibble(SAMPID = colnames(expr_sym))
  score_tbl <- add_scores(score_tbl, expr_sym, gene_sets)

  # Add selected mechanism genes as raw VST (for directionality intuition)
  key_genes <- c("ADRB1", "ADRB2", "ADRB3", "ADRA2A", "ADRA2B", "ADRA2C",
                 "TGFB1", "ACVR1C", "TNF", "NFKB1", "RELA", "IKBKE", "TBK1", "PDE3B",
                 "IL6", "IL1B", "CCL2", "NOS2")
  key_genes <- intersect(key_genes, rownames(expr_sym))
  if (length(key_genes)) {
    key_tbl <- as.data.frame(t(expr_sym[key_genes, , drop = FALSE])) %>%
      tibble::rownames_to_column("SAMPID")
    score_tbl <- dplyr::left_join(score_tbl, key_tbl, by = "SAMPID")
  }

  cib <- read_cibersort(depot_key)
  if (!is.null(cib)) {
    # Keep only fractions, not CIBERSORT stats columns.
    stats_cols <- c("P-value", "Correlation", "RMSE")
    frac_cols <- setdiff(names(cib), c("SAMPID", stats_cols, "depot", "sex", "age_bin_label", "SMCENTER"))
    cib_small <- cib %>% dplyr::select(SAMPID, dplyr::all_of(frac_cols))
    score_tbl <- dplyr::left_join(score_tbl, cib_small, by = "SAMPID")
  }

  df <- meta %>%
    dplyr::select(SAMPID, depot, sex, age_bin_label, age_mid, age_years, SMCENTER) %>%
    dplyr::left_join(score_tbl, by = "SAMPID")

  readr::write_tsv(df, file.path(tab_dir, paste0("kenichi_followup_scores_and_key_genes_vst__", depot_key, ".tsv")))

  all_rows[[length(all_rows) + 1]] <- df
}

all_df <- dplyr::bind_rows(all_rows)
readr::write_tsv(all_df, file.path(tab_dir, "kenichi_followup_scores_and_key_genes_vst__both_depots.tsv"))

# --- Plot: gene set trends by age bin & sex ---
plot_sets <- c(
  "kenichi__inflammation_cytokines",
  "kenichi__macrophage_markers",
  "kenichi__fibrosis",
  "kenichi__senescence",
  "review2017__tgfb_alk7_axis",
  "review2017__tnf_nfkB_ikkepsilon_tbk1_pde3b_axis"
)
plot_sets <- intersect(plot_sets, names(all_df))

trend_long <- all_df %>%
  dplyr::select(SAMPID, depot, sex, age_bin_label, dplyr::any_of(plot_sets)) %>%
  tidyr::pivot_longer(cols = dplyr::any_of(plot_sets), names_to = "set_name", values_to = "score") %>%
  dplyr::filter(!is.na(score)) %>%
  dplyr::mutate(
    set_label = dplyr::recode(
      set_name,
      kenichi__inflammation_cytokines = "Inflammation (cytokines)",
      kenichi__macrophage_markers = "Macrophage markers",
      kenichi__fibrosis = "Fibrosis/ECM",
      kenichi__senescence = "Senescence",
      review2017__tgfb_alk7_axis = "TGFb/ALK7 axis",
      review2017__tnf_nfkB_ikkepsilon_tbk1_pde3b_axis = "TNF/NFkB/IKKε/TBK1/PDE3B",
      .default = set_name
    ),
    set_label = factor(set_label, levels = unique(set_label))
  )

trend_sum <- trend_long %>%
  dplyr::group_by(depot, sex, age_bin_label, set_label) %>%
  dplyr::summarise(
    n = dplyr::n(),
    mean = mean(score, na.rm = TRUE),
    se = sd(score, na.rm = TRUE) / sqrt(n),
    .groups = "drop"
  )

p_trend <- ggplot2::ggplot(trend_sum, ggplot2::aes(x = age_bin_label, y = mean, group = sex, color = sex)) +
  ggplot2::geom_line(linewidth = 0.6) +
  ggplot2::geom_point(size = 1.3) +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - 1.96 * se, ymax = mean + 1.96 * se), width = 0.2, linewidth = 0.3, alpha = 0.7) +
  ggplot2::facet_grid(set_label ~ depot, scales = "free_y") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    panel.grid.minor = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    title = "Kenichi gene sets + catecholamine resistance mechanism axes (z-scores)",
    subtitle = "Per depot: VST -> per-gene z within depot -> mean z per set; points=mean, bars=~95% CI",
    x = "GTEx age bin",
    y = "Set score (mean z)"
  )

ggplot2::ggsave(file.path(fig_dir, "kenichi_followup_gene_set_trends_by_age_sex__both_depots.png"), p_trend,
                width = 11.5, height = 8.0, dpi = 220)

# --- Plot: correlations with receptors and CIBERSORT fractions ---
make_cor_plot <- function(df, depot_key) {
  sub <- df %>% dplyr::filter(depot == depot_key)

  vars <- c(
    plot_sets,
    c("ADRB2", "ADRB3", "ACVR1C", "PDE3B", "TBK1", "IKBKE", "TGFB1", "TNF"),
    c("adipocyte", "aspc", "endothelial", "macrophage", "t_cell", "nk_cell", "mast_cell")
  )
  vars <- unique(vars)
  vars <- vars[vars %in% names(sub)]
  mat <- sub %>% dplyr::select(dplyr::all_of(vars))

  # Correlation (pairwise complete).
  cmat <- stats::cor(as.matrix(mat), use = "pairwise.complete.obs", method = "pearson")

  keep_order <- c(
    plot_sets,
    c("ADRB2", "ADRB3"),
    c("TGFB1", "ACVR1C", "TNF", "NFKB1", "RELA", "IKBKE", "TBK1", "PDE3B"),
    c("adipocyte", "aspc", "endothelial", "macrophage", "t_cell", "nk_cell", "mast_cell")
  )
  keep_order <- keep_order[keep_order %in% colnames(cmat)]
  cmat <- cmat[keep_order, keep_order, drop = FALSE]

  cor_long <- as.data.frame(cmat) %>%
    tibble::rownames_to_column("var1") %>%
    tidyr::pivot_longer(cols = -var1, names_to = "var2", values_to = "r")

  pretty <- function(x) {
    x <- gsub("^kenichi__", "", x)
    x <- gsub("^review2017__", "review:", x)
    x
  }

  cor_long <- cor_long %>%
    dplyr::mutate(
      var1 = factor(pretty(var1), levels = unique(pretty(keep_order))),
      var2 = factor(pretty(var2), levels = unique(pretty(keep_order)))
    )

  p <- ggplot2::ggplot(cor_long, ggplot2::aes(x = var2, y = var1, fill = r)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.15) +
    ggplot2::scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-1, 1)) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = glue::glue("Correlations (Pearson r) in {depot_key} depot"),
      subtitle = "Variables: gene set scores, key mechanism genes (VST), and selected CIBERSORT fractions (noMALAT1)",
      x = NULL,
      y = NULL,
      fill = "r"
    )

  out_png <- file.path(fig_dir, paste0("kenichi_followup_correlations__", depot_key, ".png"))
  ggplot2::ggsave(out_png, p, width = 12.5, height = 9.2, dpi = 220)

  # Save matrix too.
  readr::write_tsv(
    cor_long %>% dplyr::mutate(var1 = as.character(var1), var2 = as.character(var2)),
    file.path(tab_dir, paste0("kenichi_followup_correlations__", depot_key, ".tsv"))
  )

  invisible(out_png)
}

for (dpt in c("subcutaneous", "visceral")) {
  make_cor_plot(all_df, dpt)
}

cat("\nOutputs written to:\n")
cat("  - ", tab_dir, "\n", sep = "")
cat("  - ", fig_dir, "\n", sep = "")
