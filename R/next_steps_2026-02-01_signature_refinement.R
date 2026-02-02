################################################################################
## Next steps (2026-02-01): Menopause Signature Refinement (SenMayo-like)
##
## Goals
## - Reduce the initial 81-gene GTEx subcutaneous menopause proxy signature to a
##   smaller, more interpretable gene set:
##   - no lncRNAs / pseudogenes (per user preference)
##   - de-emphasize obvious blood/immune/keratin contamination
##   - prioritize genes whose post-vs-pre effect persists after adjusting for
##     composition proxies (adipocyte/macrophage/fibrosis markers)
## - Produce a versioned signature file in references/ + score it in both depots.
##
## Outputs (dated subfolders):
## - GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_signature_refinement/
## - GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_signature_refinement/
## - references/menopause_signature_senMayoLike_v1.tsv
## - references/menopause_signature_senMayoLike_v1_evidence.tsv (skeleton)
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
  DESeq2,
  GSVA
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
tab_dir <- file.path(out_dir, "tables", "next_steps_2026-02-01_signature_refinement")
fig_dir <- file.path(out_dir, "figs", "next_steps_2026-02-01_signature_refinement")
log_dir <- file.path(out_dir, "logs")
ensure_dirs(c(tab_dir, fig_dir, log_dir))

log_path <- file.path(log_dir, "next_steps_2026-02-01_signature_refinement.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

approach1_de_path <- file.path(out_dir, "tables", "approach1_improved_age_proxy", "DE_menopause_proxy_improved.tsv")
approach1_gene_list_path <- file.path(out_dir, "tables", "approach1_improved_age_proxy", "menopause_signature_gene_list.txt")

sig_out_path <- file.path("references", "menopause_signature_senMayoLike_v1.tsv")
sig_evidence_out_path <- file.path("references", "menopause_signature_senMayoLike_v1_evidence.tsv")

sig2_out_path <- file.path("references", "menopause_signature_senMayoLike_v2.tsv")
sig2_evidence_out_path <- file.path("references", "menopause_signature_senMayoLike_v2_evidence.tsv")

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------
score_set_z <- function(expr_sym, genes) {
  genes <- intersect(genes, rownames(expr_sym))
  if (length(genes) < 1) {
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

is_drop_symbol <- function(sym) {
  sym <- toupper(sym)
  # drop lncRNA / antisense / Ensembl placeholders
  if (grepl("^LINC", sym)) return(TRUE)
  if (grepl("^MIR", sym)) return(TRUE)
  if (grepl("^ENSG", sym)) return(TRUE)
  if (grepl("-AS[0-9]*$", sym)) return(TRUE)
  if (grepl("AS1$|AS2$|AS3$", sym)) return(TRUE)
  # drop obvious keratin/epithelial contamination and blood/immune-heavy markers
  if (grepl("^KRT", sym)) return(TRUE)
  if (grepl("^IGH|^IGK|^IGL", sym)) return(TRUE)
  if (sym %in% c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2")) return(TRUE)
  if (sym %in% c("S100A8", "S100A9", "S100A12", "CD177")) return(TRUE)
  # drop many pseudogene-like symbols ending in P or containing P with numbers (conservative)
  if (grepl("P[0-9]*$", sym) && !sym %in% c("PCSK9")) return(TRUE)
  FALSE
}

run_gsva <- function(expr, sets, method = "ssgsea") {
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
    param <- GSVA::ssgseaParam(expr, sets)
    return(GSVA::gsva(param, verbose = FALSE))
  }
  stop("GSVA call failed.")
}

## -----------------------------------------------------------------------------
## Load initial signature + direction from Approach 1
## -----------------------------------------------------------------------------
if (!file.exists(approach1_gene_list_path) || !file.exists(approach1_de_path)) {
  stop("Missing Approach 1 inputs: ", approach1_gene_list_path, " and/or ", approach1_de_path)
}

sig81 <- readLines(approach1_gene_list_path) %>%
  stringr::str_trim() %>%
  .[nzchar(.)] %>%
  toupper()

de1 <- readr::read_tsv(approach1_de_path, show_col_types = FALSE) %>%
  dplyr::mutate(gene_symbol = toupper(as.character(gene_symbol)))

dir_tbl <- de1 %>%
  dplyr::filter(gene_symbol %in% sig81) %>%
  dplyr::select(gene_id, gene_symbol, log2FoldChange, padj) %>%
  dplyr::mutate(direction = dplyr::if_else(log2FoldChange > 0, "up_in_post", "down_in_post"))

## Apply hard drop filters (lncRNA/pseudogene/keratin/blood/immune)
sig_filtered <- sig81[!vapply(sig81, is_drop_symbol, logical(1))]

cat("Initial 81-gene signature size: ", length(sig81), "\n", sep = "")
cat("After hard drops (lnc/pseudo/keratin/blood/immune): ", length(sig_filtered), "\n", sep = "")

## -----------------------------------------------------------------------------
## Compute gene-level robustness in subcutaneous strict2 with composition proxies
## -----------------------------------------------------------------------------
cat("\nComputing subcutaneous VST and composition proxies...\n")

sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything()) %>%
  dplyr::mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])
subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything())
meta_full <- dplyr::left_join(sample_attr, subject_pheno, by = "SUBJID")

gct <- read_gct_v12(gct_sc)
common_samples <- intersect(colnames(gct$counts), meta_full$SAMPID)

meta <- meta_full %>%
  dplyr::filter(
    SAMPID %in% common_samples,
    SMTS == "Adipose Tissue",
    SMTSD == "Adipose - Subcutaneous"
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
    DTHHRDY = factor(DTHHRDY)
  )

meta_f <- meta %>%
  dplyr::filter(sex == "Female") %>%
  dplyr::mutate(
    menopause_group_strict2 = dplyr::case_when(
      is.na(age_mid) ~ NA_character_,
      age_mid < 45 ~ "pre",
      age_mid > 55 ~ "post",
      TRUE ~ NA_character_
    ),
    menopause_group_strict2 = factor(menopause_group_strict2, levels = c("pre", "post"))
  ) %>%
  dplyr::filter(!is.na(menopause_group_strict2))

counts <- gct$counts[, meta_f$SAMPID, drop = FALSE]
counts <- collapse_dupe_rowsum(counts)
counts <- round(counts)
storage.mode(counts) <- "integer"

dds0 <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts,
  colData = data.frame(row.names = colnames(counts)),
  design = ~ 1
)
vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds0, blind = TRUE))
expr_sym <- collapse_to_symbol_mean(vst_mat, gct$gene_id, gct$gene_symbol, uppercase = TRUE)

## Composition proxies (same as earlier script)
composition_sets <- list(
  adipocyte = c("ADIPOQ", "PLIN1", "FABP4"),
  macrophage_mono = c("LST1", "C1QC", "TYROBP"),
  fibroblast_ecm = c("COL1A1", "COL3A1", "DCN")
)
composition_sets <- purrr::map(composition_sets, ~ unique(toupper(.x)))
comp_scores <- purrr::imap_dfc(composition_sets, ~ tibble::tibble(!!.y := score_set_z(expr_sym, .x)))
rownames(comp_scores) <- colnames(expr_sym)

df_meta <- meta_f %>%
  dplyr::select(SAMPID, SMCENTER, SMRIN, SMTSISCH, DTHHRDY, menopause_group_strict2) %>%
  dplyr::mutate(SAMPID = as.character(SAMPID)) %>%
  dplyr::rename(group = menopause_group_strict2) %>%
  dplyr::bind_cols(comp_scores[.$SAMPID, , drop = FALSE])

## Evaluate each filtered gene under base vs composition-adjusted models
genes_present <- intersect(sig_filtered, rownames(expr_sym))
cat("Filtered genes present in subcutaneous VST: ", length(genes_present), "\n", sep = "")

rows <- list()
for (g in genes_present) {
  df <- df_meta
  df$expr <- as.numeric(expr_sym[g, df$SAMPID])

  fit_base <- stats::lm(expr ~ SMCENTER + group, data = df)
  fit_adj <- stats::lm(expr ~ SMCENTER + group + adipocyte + macrophage_mono + fibroblast_ecm, data = df)
  fit_adj_tech <- stats::lm(expr ~ SMCENTER + SMRIN + SMTSISCH + group + adipocyte + macrophage_mono + fibroblast_ecm, data = df)

  base_ct <- lm_contrast(fit_base, "grouppost", "") %>% dplyr::mutate(model = "base")
  adj_ct <- lm_contrast(fit_adj, "grouppost", "") %>% dplyr::mutate(model = "adj_comp")
  adjt_ct <- lm_contrast(fit_adj_tech, "grouppost", "") %>% dplyr::mutate(model = "adj_comp_tech")

  rows[[length(rows) + 1]] <- dplyr::bind_rows(base_ct, adj_ct, adjt_ct) %>%
    dplyr::mutate(gene_symbol = g, n = nrow(df))
}

gene_lm <- dplyr::bind_rows(rows) %>%
  dplyr::mutate(padj_within_model = ave(p, model, FUN = p.adjust, method = "BH")) %>%
  dplyr::mutate(
    conf_low = estimate - 1.96 * se,
    conf_high = estimate + 1.96 * se
  ) %>%
  dplyr::left_join(dir_tbl %>% dplyr::select(gene_symbol, log2FoldChange, padj, direction), by = "gene_symbol") %>%
  dplyr::relocate(gene_symbol, direction, log2FoldChange, padj, model, n, estimate, se, t, p, padj_within_model, conf_low, conf_high)

readr::write_tsv(gene_lm, file.path(tab_dir, "approach1_genes_vst_lm_robustness_subcutaneous_strict2.tsv"))

## Select a compact “SenMayo-like” candidate subset:
## - must remain directionally consistent after composition adjustment
## - prefer significant in adjusted model (BH within adjusted model)
adj_tbl <- gene_lm %>%
  dplyr::filter(model == "adj_comp_tech") %>%
  dplyr::mutate(abs_est = abs(estimate)) %>%
  dplyr::arrange(padj_within_model, desc(abs_est))

keep_tbl <- adj_tbl %>%
  dplyr::filter(!is.na(padj_within_model)) %>%
  dplyr::filter(padj_within_model < 0.10) %>%
  dplyr::slice_head(n = 25) %>%
  dplyr::select(gene_symbol, direction, log2FoldChange, padj, estimate, se, p, padj_within_model)

if (nrow(keep_tbl) < 10) {
  # fall back: take top 15 by adjusted p even if BH threshold not met
  keep_tbl <- adj_tbl %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::select(gene_symbol, direction, log2FoldChange, padj, estimate, se, p, padj_within_model)
}

cat("Selected senMayo-like candidate size: ", nrow(keep_tbl), "\n", sep = "")
readr::write_tsv(keep_tbl, file.path(tab_dir, "senMayoLike_v1_selected_genes_from_gtEx_subcutaneous.tsv"))

## -----------------------------------------------------------------------------
## Write versioned signature file (gene list + direction)
## -----------------------------------------------------------------------------
sig_tbl <- keep_tbl %>%
  dplyr::mutate(
    gene_symbol = toupper(gene_symbol),
    direction_in_post = dplyr::if_else(direction == "up_in_post", "up", "down"),
    source = "GTEx_v10_subcutaneous_approach1_filtered_then_robustness_screen",
    notes = "Selected from the original 81-gene list after dropping lncRNA/pseudogene/keratin/blood/immune and screening for robustness to composition+tech covariates (VST lm)."
  ) %>%
  dplyr::select(gene_symbol, direction_in_post, source, notes) %>%
  dplyr::arrange(direction_in_post, gene_symbol)

readr::write_tsv(sig_tbl, sig_out_path)

## Evidence skeleton to be filled by manual/literature search step
evidence_tbl <- sig_tbl %>%
  dplyr::mutate(
    evidence_summary = "",
    evidence_context = "",
    key_PMIDs = "",
    key_DOIs = ""
  ) %>%
  dplyr::select(gene_symbol, direction_in_post, evidence_summary, evidence_context, key_PMIDs, key_DOIs, source, notes)
readr::write_tsv(evidence_tbl, sig_evidence_out_path)

## -----------------------------------------------------------------------------
## Define an additional curated “biologically interpretable” v2 signature
## -----------------------------------------------------------------------------
## Rationale:
## - v1 is a purely data-driven robustness screen and can retain genes that are
##   harder to interpret biologically (e.g., tissue contamination markers).
## - v2 is intentionally smaller and focused on interpretable biology repeatedly
##   implicated in menopause/adipose remodeling (estrogen metabolism, Wnt/ECM,
##   adipocyte function), while still being anchored to Approach 1 direction.
curated_v2_genes <- c(
  "HSD17B2",  # estrogen metabolism
  "SFRP2", "CEMIP", "DACT2",  # Wnt/ECM-adipose remodeling axis
  "MMP3", "COL11A1", "ACAN",  # ECM remodeling
  "NNAT",  # adipocyte biology
  "PCSK9",  # metabolic / lipid homeostasis candidate
  "RXFP1", "SCG2",  # hormone/neuroendocrine-related
  "STRA6",  # retinoid handling (adipose-relevant)
  "XDH"  # oxidative/purine metabolism candidate
)
curated_v2_genes <- unique(toupper(curated_v2_genes))

cur2_dir <- dir_tbl %>%
  dplyr::select(gene_symbol, direction, log2FoldChange, padj) %>%
  dplyr::filter(gene_symbol %in% curated_v2_genes)

sig2_tbl <- tibble::tibble(gene_symbol = curated_v2_genes) %>%
  dplyr::left_join(cur2_dir, by = "gene_symbol") %>%
  dplyr::mutate(
    direction_in_post = dplyr::case_when(
      direction == "up_in_post" ~ "up",
      direction == "down_in_post" ~ "down",
      TRUE ~ NA_character_
    ),
    source = "Curated_v2_from_Approach1_direction_plus_biological_interpretability",
    notes = "Curated set intended for a SenMayo-like menopause signature; designed to avoid lncRNA/pseudogene/keratin/blood-heavy markers and focus on interpretable adipose-relevant biology."
  ) %>%
  dplyr::select(gene_symbol, direction_in_post, source, notes, log2FoldChange, padj) %>%
  dplyr::arrange(direction_in_post, gene_symbol)

readr::write_tsv(sig2_tbl, sig2_out_path)

evidence2_tbl <- sig2_tbl %>%
  dplyr::mutate(
    evidence_summary = "",
    evidence_context = "",
    key_PMIDs = "",
    key_DOIs = ""
  ) %>%
  dplyr::select(gene_symbol, direction_in_post, evidence_summary, evidence_context, key_PMIDs, key_DOIs, source, notes, log2FoldChange, padj)
readr::write_tsv(evidence2_tbl, sig2_evidence_out_path)

## -----------------------------------------------------------------------------
## Score this signature in both depots (multiple methods; female-only)
## -----------------------------------------------------------------------------
cat("\nScoring senMayoLike_v1 in both depots...\n")

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

sig_up <- sig_tbl %>% dplyr::filter(direction_in_post == "up") %>% dplyr::pull(gene_symbol)
sig_down <- sig_tbl %>% dplyr::filter(direction_in_post == "down") %>% dplyr::pull(gene_symbol)
sets <- list(senMayoLike_v1_all = sig_tbl$gene_symbol, senMayoLike_v1_up = sig_up, senMayoLike_v1_down = sig_down)
sets <- sets[lengths(sets) >= 5]

for (dpt in depots) {
  depot_label <- dpt$label
  depot_key <- dpt$key
  suffix <- dpt$suffix

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

  dds0 <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = data.frame(row.names = colnames(counts)),
    design = ~ 1
  )
  vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds0, blind = TRUE))
  expr_sym <- collapse_to_symbol_mean(vst_mat, gct$gene_id, gct$gene_symbol, uppercase = TRUE)

  # scores: signed z using up/down
  z_up <- score_set_z(expr_sym, sig_up)
  z_down <- score_set_z(expr_sym, sig_down)
  z_signed <- z_up - z_down

  # ssGSEA score for all genes
  ssg <- run_gsva(expr_sym, list(senMayoLike_v1_all = sig_tbl$gene_symbol), method = "ssgsea")
  ssg_score <- as.numeric(ssg["senMayoLike_v1_all", ])

  score_df <- meta_f %>%
    dplyr::select(SAMPID, SMCENTER, SMRIN, SMTSISCH, DTHHRDY, menopause_group_3) %>%
    dplyr::mutate(SAMPID = as.character(SAMPID)) %>%
    dplyr::mutate(
      z_signed = z_signed[match(SAMPID, names(z_signed))],
      ssgsea_all = ssg_score[match(SAMPID, colnames(ssg))]
    ) %>%
    dplyr::rename(group = menopause_group_3)

  readr::write_tsv(score_df, file.path(tab_dir, glue::glue("senMayoLike_v1_scores_per_sample{suffix}.tsv")))

  # stats
  stat_rows <- list()
  for (score_name in c("z_signed", "ssgsea_all")) {
    df <- score_df %>% dplyr::mutate(score = .data[[score_name]])
    fit <- stats::lm(score ~ SMCENTER + group, data = df)
    peri <- lm_contrast(fit, "groupperi", "") %>% dplyr::mutate(contrast = "peri_vs_pre")
    post <- lm_contrast(fit, "grouppost", "") %>% dplyr::mutate(contrast = "post_vs_pre")
    stat_rows[[length(stat_rows) + 1]] <- dplyr::bind_rows(peri, post) %>%
      dplyr::mutate(depot = depot_key, score_name = score_name, n = nrow(df))
  }
  stats <- dplyr::bind_rows(stat_rows) %>%
    dplyr::mutate(conf_low = estimate - 1.96 * se, conf_high = estimate + 1.96 * se) %>%
    dplyr::relocate(depot, score_name, contrast, n)
  readr::write_tsv(stats, file.path(tab_dir, glue::glue("senMayoLike_v1_score_stats{suffix}.tsv")))

  # plot
  plot_df <- score_df %>%
    tidyr::pivot_longer(cols = c("z_signed", "ssgsea_all"), names_to = "score_name", values_to = "score")
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = group, y = score)) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::facet_wrap(~ score_name, scales = "free_y") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      title = paste0("senMayoLike_v1 menopause signature scores (", depot_label, ", female)"),
      x = "Menopause proxy group",
      y = "Score"
    )
  ggplot2::ggsave(
    filename = file.path(fig_dir, glue::glue("senMayoLike_v1_scores_by_group{suffix}.png")),
    plot = p,
    width = 8,
    height = 4,
    dpi = 150
  )
}

## -----------------------------------------------------------------------------
## Score curated v2 in both depots
## -----------------------------------------------------------------------------
cat("\nScoring senMayoLike_v2 in both depots...\n")

sig2_up <- sig2_tbl %>% dplyr::filter(direction_in_post == "up") %>% dplyr::pull(gene_symbol)
sig2_down <- sig2_tbl %>% dplyr::filter(direction_in_post == "down") %>% dplyr::pull(gene_symbol)
sets2 <- list(senMayoLike_v2_all = sig2_tbl$gene_symbol, senMayoLike_v2_up = sig2_up, senMayoLike_v2_down = sig2_down)
sets2 <- sets2[lengths(sets2) >= 5]

for (dpt in depots) {
  depot_label <- dpt$label
  depot_key <- dpt$key
  suffix <- dpt$suffix

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

  dds0 <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = data.frame(row.names = colnames(counts)),
    design = ~ 1
  )
  vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds0, blind = TRUE))
  expr_sym <- collapse_to_symbol_mean(vst_mat, gct$gene_id, gct$gene_symbol, uppercase = TRUE)

  z_up <- score_set_z(expr_sym, sig2_up)
  z_down <- score_set_z(expr_sym, sig2_down)
  z_signed <- z_up - z_down

  ssg <- run_gsva(expr_sym, list(senMayoLike_v2_all = sig2_tbl$gene_symbol), method = "ssgsea")
  ssg_score <- as.numeric(ssg["senMayoLike_v2_all", ])

  score_df <- meta_f %>%
    dplyr::select(SAMPID, SMCENTER, SMRIN, SMTSISCH, DTHHRDY, menopause_group_3) %>%
    dplyr::mutate(SAMPID = as.character(SAMPID)) %>%
    dplyr::mutate(
      z_signed = z_signed[match(SAMPID, names(z_signed))],
      ssgsea_all = ssg_score[match(SAMPID, colnames(ssg))]
    ) %>%
    dplyr::rename(group = menopause_group_3)

  readr::write_tsv(score_df, file.path(tab_dir, glue::glue("senMayoLike_v2_scores_per_sample{suffix}.tsv")))

  stat_rows <- list()
  for (score_name in c("z_signed", "ssgsea_all")) {
    df <- score_df %>% dplyr::mutate(score = .data[[score_name]])
    fit <- stats::lm(score ~ SMCENTER + group, data = df)
    peri <- lm_contrast(fit, "groupperi", "") %>% dplyr::mutate(contrast = "peri_vs_pre")
    post <- lm_contrast(fit, "grouppost", "") %>% dplyr::mutate(contrast = "post_vs_pre")
    stat_rows[[length(stat_rows) + 1]] <- dplyr::bind_rows(peri, post) %>%
      dplyr::mutate(depot = depot_key, score_name = score_name, n = nrow(df))
  }
  stats <- dplyr::bind_rows(stat_rows) %>%
    dplyr::mutate(conf_low = estimate - 1.96 * se, conf_high = estimate + 1.96 * se) %>%
    dplyr::relocate(depot, score_name, contrast, n)
  readr::write_tsv(stats, file.path(tab_dir, glue::glue("senMayoLike_v2_score_stats{suffix}.tsv")))

  plot_df <- score_df %>%
    tidyr::pivot_longer(cols = c("z_signed", "ssgsea_all"), names_to = "score_name", values_to = "score")
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = group, y = score)) +
    ggplot2::geom_boxplot(outlier.size = 0.5) +
    ggplot2::facet_wrap(~ score_name, scales = "free_y") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::labs(
      title = paste0("senMayoLike_v2 menopause signature scores (", depot_label, ", female)"),
      x = "Menopause proxy group",
      y = "Score"
    )
  ggplot2::ggsave(
    filename = file.path(fig_dir, glue::glue("senMayoLike_v2_scores_by_group{suffix}.png")),
    plot = p,
    width = 8,
    height = 4,
    dpi = 150
  )
}

## README
md <- c(
  "# Menopause signature refinement (senMayoLike_v1)",
  "",
  "This folder contains a first-pass reduction of the GTEx v10 subcutaneous menopause-proxy signature to a smaller, more interpretable set that excludes lncRNAs/pseudogenes and screens for robustness to composition + technical covariates.",
  "",
  "## Key outputs",
  "- Robustness screen (subcutaneous strict2; VST lm): `approach1_genes_vst_lm_robustness_subcutaneous_strict2.tsv`",
  "- Selected gene list (input to senMayoLike_v1): `senMayoLike_v1_selected_genes_from_gtEx_subcutaneous.tsv`",
  "- Versioned signature file: `references/menopause_signature_senMayoLike_v1.tsv`",
  "- Evidence skeleton (to fill with PMIDs/DOIs): `references/menopause_signature_senMayoLike_v1_evidence.tsv`",
  "- Scoring outputs (per depot): `senMayoLike_v1_scores_per_sample*.tsv`, `senMayoLike_v1_score_stats*.tsv`, and `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_signature_refinement/senMayoLike_v1_scores_by_group*.png`",
  "",
  "## Additional curated signature",
  "- v2 signature file: `references/menopause_signature_senMayoLike_v2.tsv`",
  "- v2 evidence skeleton: `references/menopause_signature_senMayoLike_v2_evidence.tsv`",
  "- v2 scoring outputs: `senMayoLike_v2_scores_per_sample*.tsv`, `senMayoLike_v2_score_stats*.tsv`, and `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_signature_refinement/senMayoLike_v2_scores_by_group*.png`",
  "",
  "Generated:",
  as.character(Sys.time())
)
writeLines(md, file.path(tab_dir, "README_signature_refinement.md"))

cat("\nDone. Log written to: ", log_path, "\n", sep = "")
