################################################################################
## Next steps (2026-02-28): Cell-fraction deconvolution (CIBERSORT; first pass)
##
## Reference:
## - GSE176171 (human+mouse WAT sc/snRNA; provides cell-level metadata with
##   broad cell types and a combined Hs10X count matrix in MatrixMarket format).
##
## Plan:
## 1) Build a compact signature matrix from the sc reference:
##    - aggregate raw counts to pseudo-bulk per cell_type__custom
##    - convert to CPM per cell type
##    - select Top-N marker genes per cell type by "margin" (top - second)
## 2) Compute GTEx bulk mixture CPM (per depot) and run IOBR::CIBERSORT.
##
## Outputs:
## - external_data/GSE176171/derived/signature_topN{N}_cpm.tsv.gz
## - GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/
## - GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_cibersort_fractions/
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
  Matrix, matrixStats,
  DESeq2,
  IOBR
)

source("R/utils.R")

## -----------------------------------------------------------------------------
## Config
## -----------------------------------------------------------------------------
ref_dir <- file.path("external_data", "GSE176171")
ref_meta_gz <- file.path(ref_dir, "GSE176171_cell_metadata.tsv.gz")
ref_mtx_gz <- file.path(ref_dir, "GSE176171_Hs10X.counts.mtx.gz")
ref_feat_gz <- file.path(ref_dir, "GSE176171_Hs10X.counts.features.tsv.gz")
ref_bar_gz <- file.path(ref_dir, "GSE176171_Hs10X.counts.barcodes.tsv.gz")

top_n <- suppressWarnings(as.integer(Sys.getenv("CIBERSORT_SIGNATURE_TOP_N", "50")))
if (is.na(top_n) || top_n <= 0) top_n <- 50L
min_cells_per_type <- suppressWarnings(as.integer(Sys.getenv("CIBERSORT_MIN_CELLS_PER_TYPE", "200")))
if (is.na(min_cells_per_type) || min_cells_per_type <= 0) min_cells_per_type <- 200L

exclude_cell_types <- c("neutrophil", "endometrium")  # keep first-pass reference compact

# Genes to exclude from both signature+mixture. Default drops MALAT1, which can
# dominate scaling in IOBR::CIBERSORT and collapse adipocyte fractions to ~0.
exclude_ensg <- Sys.getenv("CIBERSORT_EXCLUDE_ENSG", "ENSG00000251562")
exclude_ensg <- unlist(strsplit(exclude_ensg, ",", fixed = TRUE))
exclude_ensg <- trimws(exclude_ensg)
exclude_ensg <- exclude_ensg[nzchar(exclude_ensg)]

perm_n <- suppressWarnings(as.integer(Sys.getenv("CIBERSORT_PERM", "0")))
if (is.na(perm_n) || perm_n < 0) perm_n <- 0L

qn_flag <- tolower(Sys.getenv("CIBERSORT_QN", "false")) %in% c("1", "true", "t", "yes", "y")

# IOBR::CIBERSORT uses parallel::mclapply() internally in some paths; keep it
# deterministic and less fragile by default.
mc_cores <- suppressWarnings(as.integer(Sys.getenv("CIBERSORT_CORES", "1")))
if (is.na(mc_cores) || mc_cores <= 0) mc_cores <- 1L
options(mc.cores = mc_cores)

out_dir <- "GTEx_v10_AT_analysis_out"
tab_dir <- file.path(out_dir, "tables", "next_steps_2026-02-28_cibersort_fractions")
fig_dir <- file.path(out_dir, "figs", "next_steps_2026-02-28_cibersort_fractions")
log_dir <- file.path(out_dir, "logs")
ensure_dirs(c(tab_dir, fig_dir, log_dir, file.path(ref_dir, "derived")))

log_path <- file.path(log_dir, "next_steps_2026-02-28_cibersort_fractions.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------
stop_if_missing <- function(path, label) {
  if (!file.exists(path)) stop(label, " not found: ", path)
}

read_mm_gz <- function(path_gz) {
  # Matrix::readMM can read from a gzfile() connection.
  Matrix::readMM(gzfile(path_gz))
}

strip_ver <- function(x) sub("\\..*$", "", as.character(x))

compute_cpm <- function(count_mat) {
  lib <- colSums(count_mat)
  lib[lib == 0] <- 1
  sweep(count_mat, 2, lib / 1e6, "/")
}

prune_signature_topN_margin <- function(signature_mat, top_n) {
  mat <- as.matrix(signature_mat)
  if (ncol(mat) < 2) return(list(sig = signature_mat, selected = rownames(signature_mat)))
  top_idx <- max.col(mat, ties.method = "first")
  top1 <- matrixStats::rowMaxs(mat, na.rm = TRUE)
  mat2 <- mat
  mat2[cbind(seq_len(nrow(mat2)), top_idx)] <- -Inf
  top2 <- matrixStats::rowMaxs(mat2, na.rm = TRUE)
  margin <- top1 - top2
  df <- tibble::tibble(
    gene = rownames(mat),
    top_group = colnames(mat)[top_idx],
    margin = margin
  ) %>%
    dplyr::filter(is.finite(.data$margin)) %>%
    dplyr::group_by(.data$top_group) %>%
    dplyr::arrange(dplyr::desc(.data$margin)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
  selected <- unique(df$gene)
  list(sig = signature_mat[selected, , drop = FALSE], selected = selected, marker_tbl = df)
}

run_cibersort_one <- function(signature_mat, mixture_mat, out_prefix) {
  common <- intersect(rownames(signature_mat), rownames(mixture_mat))
  if (length(common) < 200) stop("Too few overlapping genes between signature and mixture (n=", length(common), ").")
  sig <- signature_mat[common, , drop = FALSE]
  mix <- mixture_mat[common, , drop = FALSE]

  cat("Running CIBERSORT: genes=", length(common), " types=", ncol(sig), " samples=", ncol(mix),
      " perm=", perm_n, " QN=", qn_flag, "\n", sep = "")

  est <- IOBR::CIBERSORT(
    sig_matrix = as.data.frame(sig),
    mixture_file = as.data.frame(mix),
    perm = perm_n,
    QN = qn_flag,
    absolute = FALSE
  )

  # Save everything needed to reproduce the run.
  readr::write_tsv(
    tibble::tibble(Gene = rownames(sig)) %>% dplyr::bind_cols(as.data.frame(sig)),
    paste0(out_prefix, "__signature_used.tsv.gz")
  )
  readr::write_tsv(
    tibble::tibble(Gene = rownames(mix)) %>% dplyr::bind_cols(as.data.frame(mix)),
    paste0(out_prefix, "__mixture_used.tsv.gz")
  )
  est_df <- as.data.frame(est)
  est_df <- tibble::rownames_to_column(est_df, "Mixture")
  readr::write_csv(est_df, paste0(out_prefix, "__cibersort_fractions.csv"))
  est_df
}

## -----------------------------------------------------------------------------
## 1) Build (or load) reference signature
## -----------------------------------------------------------------------------
stop_if_missing(ref_meta_gz, "Reference metadata")
stop_if_missing(ref_mtx_gz, "Reference matrix (Hs10X mtx.gz)")
stop_if_missing(ref_feat_gz, "Reference features")
stop_if_missing(ref_bar_gz, "Reference barcodes")

sig_out_base <- file.path(ref_dir, "derived", glue::glue("signature_topN{top_n}_cpm_noversion"))
sig_tsv_gz <- paste0(sig_out_base, ".tsv.gz")
marker_tsv_gz <- paste0(sig_out_base, "__marker_table.tsv.gz")

if (!file.exists(sig_tsv_gz)) {
  # Prefer using an existing cached full reference CPM (fast), if available.
  cached_full_rds <- file.path(ref_dir, "derived", "signature_topN200_cpm_noversion.rds")
  if (file.exists(cached_full_rds)) {
    cat("\n[Ref] Building signature from cached ref CPM: ", cached_full_rds, "\n", sep = "")
    cached <- readRDS(cached_full_rds)
    ref_cpm <- cached$ref_cpm_full
    pruned <- prune_signature_topN_margin(ref_cpm, top_n = top_n)
    sig <- pruned$sig
    markers <- pruned$marker_tbl

    cat("[Ref] Writing signature + marker table...\n")
    readr::write_tsv(
      tibble::tibble(Gene = rownames(sig)) %>% dplyr::bind_cols(as.data.frame(sig)),
      sig_tsv_gz
    )
    readr::write_tsv(markers, marker_tsv_gz)

    saveRDS(
      list(
        signature_cpm = sig,
        ref_cpm_full = ref_cpm,
        selected_genes = pruned$selected,
        cell_type_counts = cached$cell_type_counts,
        kept_cell_types = colnames(ref_cpm),
        exclude_cell_types = cached$exclude_cell_types,
        top_n = top_n,
        min_cells_per_type = cached$min_cells_per_type,
        ref_paths = cached$ref_paths
      ),
      file.path(ref_dir, "derived", glue::glue("signature_topN{top_n}_cpm_noversion.rds"))
    )
  } else {
  cat("\n[Ref] Loading cell metadata...\n")
  meta <- readr::read_tsv(ref_meta_gz, show_col_types = FALSE) %>%
    dplyr::filter(species__ontology_label == "Homo sapiens") %>%
    dplyr::mutate(
      cell_id = as.character(cell_id),
      cell_type = tolower(as.character(cell_type__custom))
    ) %>%
    dplyr::filter(!is.na(cell_type), cell_type != "", !(cell_type %in% exclude_cell_types))

  cat("[Ref] Loading features + barcodes...\n")
  feats <- readr::read_tsv(ref_feat_gz, col_names = c("feature_id", "feature_name", "feature_type"), show_col_types = FALSE)
  bars <- readr::read_tsv(ref_bar_gz, col_names = c("cell_id"), show_col_types = FALSE)

  # Keep only gene expression rows
  keep_feat <- feats$feature_type == "Gene Expression"
  feats <- feats[keep_feat, , drop = FALSE]

  cat("[Ref] Reading sparse matrix (this can take a bit)...\n")
  m <- read_mm_gz(ref_mtx_gz)
  # Matrix now deprecates dgT -> dgC direct coercion; use CsparseMatrix.
  m <- methods::as(m, "CsparseMatrix")

  n_feat_all <- length(keep_feat)
  if (nrow(m) != n_feat_all) {
    stop("Unexpected ref matrix row count. m=", nrow(m), " feats(all)=", n_feat_all)
  }
  # Keep only "Gene Expression" features (no-op when keep_feat is all TRUE).
  m <- m[keep_feat, , drop = FALSE]
  if (ncol(m) != nrow(bars)) stop("Unexpected ref matrix col count. m=", ncol(m), " barcodes=", nrow(bars))

  rownames(m) <- feats$feature_id
  colnames(m) <- bars$cell_id

  # Align to metadata cells (intersection); keep order to match colnames(m)
  meta <- meta %>% dplyr::filter(cell_id %in% colnames(m))
  meta <- meta %>% dplyr::distinct(cell_id, .keep_all = TRUE)

  # Drop cells missing from metadata
  keep_cells <- intersect(colnames(m), meta$cell_id)
  m <- m[, keep_cells, drop = FALSE]
  meta <- meta %>% dplyr::filter(cell_id %in% keep_cells)
  meta <- meta[match(colnames(m), meta$cell_id), , drop = FALSE]
  stopifnot(all(meta$cell_id == colnames(m)))

  # Filter cell types with enough cells
  ct_counts <- meta %>%
    dplyr::count(cell_type, name = "n_cells") %>%
    dplyr::arrange(dplyr::desc(n_cells))
  readr::write_tsv(ct_counts, file.path(ref_dir, "derived", "reference_cell_type_counts.tsv"))

  keep_types <- ct_counts %>%
    dplyr::filter(n_cells >= min_cells_per_type) %>%
    dplyr::pull(cell_type)
  if (length(keep_types) < 5) stop("Too few cell types passed min_cells_per_type=", min_cells_per_type)

  meta <- meta %>% dplyr::filter(cell_type %in% keep_types)
  m <- m[, meta$cell_id, drop = FALSE]

  cat("[Ref] Aggregating to pseudo-bulk counts per cell_type (types=", length(keep_types), ")...\n", sep = "")
  types <- sort(unique(meta$cell_type))
  idx <- split(seq_len(nrow(meta)), meta$cell_type)
  counts_by_type <- matrix(0, nrow = nrow(m), ncol = length(types))
  rownames(counts_by_type) <- rownames(m)
  colnames(counts_by_type) <- types
  for (k in types) {
    counts_by_type[, k] <- Matrix::rowSums(m[, idx[[k]], drop = FALSE])
  }

  # Strip Ensembl version suffixes to increase overlap with GTEx, then collapse duplicates.
  rownames(counts_by_type) <- strip_ver(rownames(counts_by_type))
  if (any(duplicated(rownames(counts_by_type)))) {
    counts_by_type <- rowsum(counts_by_type, group = rownames(counts_by_type), reorder = FALSE)
  }

  ref_cpm <- compute_cpm(counts_by_type)
  pruned <- prune_signature_topN_margin(ref_cpm, top_n = top_n)
  sig <- pruned$sig
  markers <- pruned$marker_tbl

  cat("[Ref] Writing signature + marker table...\n")
  readr::write_tsv(
    tibble::tibble(Gene = rownames(sig)) %>% dplyr::bind_cols(as.data.frame(sig)),
    sig_tsv_gz
  )
  readr::write_tsv(markers, marker_tsv_gz)

  saveRDS(
    list(
      signature_cpm = sig,
      ref_cpm_full = ref_cpm,
      selected_genes = pruned$selected,
      cell_type_counts = ct_counts,
      kept_cell_types = types,
      exclude_cell_types = exclude_cell_types,
      top_n = top_n,
      min_cells_per_type = min_cells_per_type,
      ref_paths = list(mtx = ref_mtx_gz, features = ref_feat_gz, barcodes = ref_bar_gz, meta = ref_meta_gz)
    ),
    file.path(ref_dir, "derived", glue::glue("signature_topN{top_n}_cpm_noversion.rds"))
  )
  }
} else {
  cat("[Ref] Using existing signature: ", sig_tsv_gz, "\n", sep = "")
}

sig_df <- readr::read_tsv(sig_tsv_gz, show_col_types = FALSE)
sig_mat <- as.matrix(sig_df[, -1, drop = FALSE])
rownames(sig_mat) <- sig_df$Gene
if (length(exclude_ensg)) {
  sig_mat <- sig_mat[setdiff(rownames(sig_mat), exclude_ensg), , drop = FALSE]
}

## -----------------------------------------------------------------------------
## 2) GTEx mixture CPM + run CIBERSORT (subcutaneous + visceral)
## -----------------------------------------------------------------------------
in_dir <- "GTExDatav10"
gct_sc <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
gct_vis <- file.path(in_dir, "gene_reads_v10_adipose_visceral_omentum.gct.gz")
sattr <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

depots <- list(
  list(key = "subcutaneous", label = "Adipose - Subcutaneous", gct_path = gct_sc, suffix = ""),
  list(key = "visceral", label = "Adipose - Visceral (Omentum)", gct_path = gct_vis, suffix = "_visceral")
)

cat("\n[GTEx] Loading metadata...\n")
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
  cat("[GTEx] Depot: ", depot_label, " (", depot_key, ")\n", sep = "")
  cat("-------------------------------------------\n")

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
      depot = depot_key
    )

  counts <- gct$counts[, meta$SAMPID, drop = FALSE]
  counts <- collapse_dupe_rowsum(counts)
  rownames(counts) <- strip_ver(rownames(counts))
  if (any(duplicated(rownames(counts)))) {
    counts <- rowsum(counts, group = rownames(counts), reorder = FALSE)
  }

  # Restrict to genes in the signature
  common_genes <- intersect(rownames(counts), rownames(sig_mat))
  if (length(exclude_ensg)) common_genes <- setdiff(common_genes, exclude_ensg)
  counts2 <- counts[common_genes, , drop = FALSE]
  cat("[GTEx] Using signature-overlap genes: ", length(common_genes), "\n", sep = "")

  mixture_cpm <- compute_cpm(counts2)

  # Run CIBERSORT
  excl_tag <- ""
  if (length(exclude_ensg)) {
    if (identical(sort(exclude_ensg), "ENSG00000251562")) {
      excl_tag <- "_noMALAT1"
    } else {
      excl_tag <- "_exclCustom"
    }
  }
  run_base <- glue::glue("cibersort_{depot_key}_perm{perm_n}_topN{top_n}_minCells{min_cells_per_type}{excl_tag}")
  run_prefix_tab <- file.path(tab_dir, run_base)
  run_prefix_fig <- file.path(fig_dir, run_base)
  est <- run_cibersort_one(sig_mat, mixture_cpm, run_prefix_tab)

  # Join back sample metadata for downstream summaries
  est2 <- est %>%
    dplyr::rename(SAMPID = Mixture) %>%
    dplyr::left_join(meta %>% dplyr::select(SAMPID, depot, sex, age_bin_label, SMCENTER), by = "SAMPID")

  readr::write_tsv(est2, paste0(run_prefix_tab, "__fractions_with_meta.tsv.gz"))

  # Summarize mean fractions by sex x age bin (quick, interpretable)
  frac_cols <- setdiff(colnames(est), c("Mixture", "P.value", "P-value", "Correlation", "RMSE"))
  sum_tbl <- est2 %>%
    dplyr::group_by(depot, sex, age_bin_label) %>%
    dplyr::summarise(
      n = dplyr::n(),
      dplyr::across(dplyr::all_of(frac_cols), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    )
  readr::write_tsv(sum_tbl, paste0(run_prefix_tab, "__fractions_mean_by_sex_agebin.tsv"))

  sum_long <- sum_tbl %>%
    tidyr::pivot_longer(cols = dplyr::all_of(frac_cols), names_to = "cell_type", values_to = "mean_fraction")

  p <- ggplot2::ggplot(sum_long, ggplot2::aes(x = age_bin_label, y = mean_fraction, fill = cell_type)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~ sex, nrow = 1) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = glue::glue("CIBERSORT fractions (mean) by sex x age bin | {depot_label}"),
      subtitle = glue::glue("Reference: GSE176171; signature TopN={top_n} (min cells/type={min_cells_per_type}); perm={perm_n}; QN={qn_flag}"),
      x = "GTEx age bin",
      y = "Mean estimated fraction",
      fill = "Cell type"
    )

  ggplot2::ggsave(paste0(run_prefix_fig, "__fractions_mean_by_sex_agebin.png"), p, width = 12.5, height = 4.2, dpi = 200)
}

cat("\nDone.\n")
cat("Signature: ", sig_tsv_gz, "\n", sep = "")
cat("Outputs: ", tab_dir, "\n", sep = "")
cat("Figures: ", fig_dir, "\n", sep = "")
cat("Log: ", log_path, "\n", sep = "")
