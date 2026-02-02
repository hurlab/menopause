################################################################################
## External dataset analysis: GSE44000
##
## Study: "Transcriptome profile of subcutaneous adipocytes isolated from obese
## vs. lean postmenopausal women" (Agilent microarray).
##
## Purpose in this project:
## - Not a menopause-status validation (all subjects are postmenopausal).
## - Provides an adipocyte-enriched context to interpret adipocyte programs
##   (lipolysis/thermogenesis/adrenergic modules) without bulk composition.
##
## Inputs (downloaded if missing):
## - Series matrix: GSE44000_series_matrix.txt.gz
## - Platform annotation: GPL6480.annot.gz
##
## Outputs:
## - GTEx_v10_AT_analysis_out/tables/external_validation_GSE44000/
## - GTEx_v10_AT_analysis_out/figs/external_validation_GSE44000/
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
  limma
)

source("R/utils.R")

out_dir <- "GTEx_v10_AT_analysis_out"
tab_dir <- file.path(out_dir, "tables", "external_validation_GSE44000")
fig_dir <- file.path(out_dir, "figs", "external_validation_GSE44000")
log_dir <- file.path(out_dir, "logs")
ensure_dirs(c(tab_dir, fig_dir, log_dir))

log_path <- file.path(log_dir, "external_validation_GSE44000.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

data_dir <- file.path("external_data", "GSE44000")
ensure_dirs(data_dir)

series_gz <- file.path(data_dir, "GSE44000_series_matrix.txt.gz")
annot_gz <- file.path(data_dir, "GPL6480.annot.gz")

if (!file.exists(series_gz)) {
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44000/matrix/GSE44000_series_matrix.txt.gz",
    destfile = series_gz,
    mode = "wb",
    quiet = TRUE
  )
}
if (!file.exists(annot_gz)) {
  download.file(
    "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6480/annot/GPL6480.annot.gz",
    destfile = annot_gz,
    mode = "wb",
    quiet = TRUE
  )
}

## -----------------------------------------------------------------------------
## Read series matrix header (metadata) and expression table
## -----------------------------------------------------------------------------
cat("Reading series matrix...\n")
lines <- readr::read_lines(series_gz)
begin_idx <- which(lines == "!series_matrix_table_begin")
end_idx <- which(lines == "!series_matrix_table_end")
stopifnot(length(begin_idx) == 1, length(end_idx) == 1, begin_idx < end_idx)

header_lines <- lines[seq_len(begin_idx - 1)]
table_lines <- lines[(begin_idx + 1):(end_idx - 1)]

tmp_path <- file.path(tempdir(), "GSE44000_matrix_table.tsv")
writeLines(table_lines, tmp_path)

expr_df <- readr::read_tsv(tmp_path, show_col_types = FALSE)
stopifnot("ID_REF" %in% names(expr_df))

## Sample titles to derive group
sample_titles_line <- header_lines[stringr::str_starts(header_lines, "!Sample_title")]
stopifnot(length(sample_titles_line) == 1)
sample_titles <- readr::read_tsv(
  I(gsub("^!Sample_title\\t", "", sample_titles_line)),
  col_names = FALSE,
  show_col_types = FALSE,
  quote = "\""
)[1, ] %>%
  unlist(use.names = FALSE) %>%
  as.character()

sample_ids_line <- header_lines[stringr::str_starts(header_lines, "!Sample_geo_accession")]
stopifnot(length(sample_ids_line) == 1)
sample_ids <- readr::read_tsv(
  I(gsub("^!Sample_geo_accession\\t", "", sample_ids_line)),
  col_names = FALSE,
  show_col_types = FALSE,
  quote = "\""
)[1, ] %>%
  unlist(use.names = FALSE) %>%
  as.character()

stopifnot(length(sample_ids) == (ncol(expr_df) - 1))
colnames(expr_df) <- c("ID_REF", sample_ids)

meta <- tibble::tibble(
  sample_id = sample_ids,
  sample_title = sample_titles
) %>%
  dplyr::mutate(
    group = dplyr::case_when(
      stringr::str_detect(sample_title, "_Obese_") ~ "Obese",
      stringr::str_detect(sample_title, "_Lean_") ~ "Lean",
      TRUE ~ NA_character_
    ),
    group = factor(group, levels = c("Lean", "Obese"))
  )
stopifnot(!any(is.na(meta$group)))

## -----------------------------------------------------------------------------
## Read GPL annotation and map probe -> gene symbol
## -----------------------------------------------------------------------------
cat("Reading platform annotation...\n")
annot_lines <- readr::read_lines(annot_gz, n_max = 2000)
begin_idx <- which(annot_lines == "!platform_table_begin")
stopifnot(length(begin_idx) == 1)

annot_df <- readr::read_tsv(
  annot_gz,
  skip = begin_idx,
  show_col_types = FALSE
)
stopifnot(all(c("ID", "Gene symbol") %in% names(annot_df)))

probe_to_symbol <- annot_df %>%
  dplyr::select(probe_id = ID, gene_symbol = `Gene symbol`) %>%
  dplyr::mutate(
    probe_id = as.character(probe_id),
    gene_symbol = toupper(as.character(gene_symbol))
  ) %>%
  dplyr::filter(!is.na(probe_id), nzchar(probe_id), !is.na(gene_symbol), nzchar(gene_symbol))

## -----------------------------------------------------------------------------
## Build expression matrix collapsed to gene symbols (mean across probes)
## -----------------------------------------------------------------------------
expr_mat <- as.matrix(expr_df[, -1, drop = FALSE])
storage.mode(expr_mat) <- "double"
rownames(expr_mat) <- expr_df$ID_REF
expr_mat <- expr_mat[rownames(expr_mat) %in% probe_to_symbol$probe_id, , drop = FALSE]

sym <- probe_to_symbol$gene_symbol[match(rownames(expr_mat), probe_to_symbol$probe_id)]
keep <- !is.na(sym) & nzchar(sym)
expr_mat <- expr_mat[keep, , drop = FALSE]
sym <- sym[keep]

collapsed <- rowsum(expr_mat, group = sym, reorder = FALSE)
denom <- as.numeric(table(sym)[rownames(collapsed)])
expr_sym <- collapsed / denom
rownames(expr_sym) <- rownames(collapsed)

## -----------------------------------------------------------------------------
## Define gene sets to score
## -----------------------------------------------------------------------------
kenichi_sets <- list(
  acute_beta_adrenergic = c("NR4A1", "NR4A2", "NR4A3", "FOS", "JUN", "EGR1"),
  adrenergic_receptors = c("ADRB1", "ADRB2", "ADRB3", "ADRA2A", "ADRA2B", "ADRA2C"),
  lipolysis_core = c("PNPLA2", "LIPE", "MGLL", "ABHD5", "PLIN1", "G0S2", "LPL", "FABP4"),
  thermogenesis_program = c("UCP1", "PPARGC1A", "PPARGC1B", "CIDEA", "DIO2", "VEGFA")
)
kenichi_sets <- purrr::map(kenichi_sets, ~ unique(toupper(.x)))

sig_v2 <- readr::read_tsv("references/menopause_signature_senMayoLike_v2.tsv", show_col_types = FALSE) %>%
  dplyr::mutate(gene_symbol = toupper(as.character(gene_symbol)))
sig2_up <- sig_v2 %>% dplyr::filter(direction_in_post == "up") %>% dplyr::pull(gene_symbol)
sig2_down <- sig_v2 %>% dplyr::filter(direction_in_post == "down") %>% dplyr::pull(gene_symbol)

score_set_z <- function(expr, genes) {
  genes <- intersect(genes, rownames(expr))
  if (!length(genes)) return(rep(NA_real_, ncol(expr)))
  m <- expr[genes, , drop = FALSE]
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  z <- sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
  colMeans(z, na.rm = TRUE)
}

scores <- purrr::imap_dfc(kenichi_sets, ~ tibble::tibble(!!.y := score_set_z(expr_sym, .x)))
scores$senMayoLike_v2_z_signed <- score_set_z(expr_sym, sig2_up) - score_set_z(expr_sym, sig2_down)
rownames(scores) <- colnames(expr_sym)

score_df <- meta %>%
  dplyr::bind_cols(scores[meta$sample_id, , drop = FALSE])

readr::write_tsv(score_df, file.path(tab_dir, "scores_per_sample.tsv"))

## -----------------------------------------------------------------------------
## Differential (Obese vs Lean) for scores (limma; tiny n)
## -----------------------------------------------------------------------------
fit_rows <- list()
for (nm in setdiff(names(score_df), c("sample_id", "sample_title", "group"))) {
  y <- score_df[[nm]]
  fit <- stats::lm(y ~ group, data = score_df)
  sm <- summary(fit)$coefficients
  if (!("groupObese" %in% rownames(sm))) next
  fit_rows[[length(fit_rows) + 1]] <- tibble::tibble(
    feature = nm,
    estimate = unname(sm["groupObese", "Estimate"]),
    t = unname(sm["groupObese", "t value"]),
    p = unname(sm["groupObese", "Pr(>|t|)"])
  )
}
stats <- dplyr::bind_rows(fit_rows) %>%
  dplyr::mutate(padj = p.adjust(p, method = "BH")) %>%
  dplyr::arrange(p)
readr::write_tsv(stats, file.path(tab_dir, "score_differential_obese_vs_lean.tsv"))

## Plot: key scores by group
plot_df <- score_df %>%
  tidyr::pivot_longer(
    cols = setdiff(names(score_df), c("sample_id", "sample_title", "group")),
    names_to = "score_name",
    values_to = "score"
  )
p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = group, y = score)) +
  ggplot2::geom_boxplot(outlier.size = 0.7) +
  ggplot2::facet_wrap(~ score_name, scales = "free_y", ncol = 3) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::labs(
    title = "GSE44000 (postmenopausal adipocytes): obese vs lean score comparison",
    x = "Group",
    y = "Score (mean z)"
  )
ggplot2::ggsave(file.path(fig_dir, "scores_by_group.png"), p, width = 12, height = 7, dpi = 150)

## README
md <- c(
  "# External dataset: GSE44000 (postmenopausal adipocytes; obese vs lean)",
  "",
  "This analysis uses the GEO series matrix (processed expression) and GPL annotation to score adipocyte gene programs.",
  "",
  "Important: this dataset does not provide a pre- vs post-menopause comparison (all are postmenopausal). It is used as an adipocyte-enriched context to interpret adipocyte programs without bulk composition confounding.",
  "",
  "## Outputs",
  "- Per-sample scores: `scores_per_sample.tsv`",
  "- Obese vs lean score stats: `score_differential_obese_vs_lean.tsv`",
  "- Figure: `GTEx_v10_AT_analysis_out/figs/external_validation_GSE44000/scores_by_group.png`",
  "",
  "Generated:",
  as.character(Sys.time())
)
writeLines(md, file.path(tab_dir, "README_GSE44000.md"))

cat("Done. Log written to: ", log_path, "\n", sep = "")
