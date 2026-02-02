#' Ensure that a set of directories exist.
#'
#' This is a thin wrapper around [dir.create()] that quietly creates each
#' directory if missing and returns the normalized paths invisibly.
#'
#' @param paths Character vector of directory paths.
#' @return Invisibly returns the normalized paths.
#' @examples
#' ensure_dirs(tempdir())
# Ensure a default CRAN mirror is set for non-interactive runs
.repos <- getOption("repos")
.cran <- unname(.repos["CRAN"])
if (is.null(.cran) || is.na(.cran) || .cran == "@CRAN@" || .cran == "") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
rm(.repos, .cran)

ensure_dirs <- function(paths) {
  stopifnot(is.character(paths))
  paths <- paths[nzchar(paths)]
  for (pth in paths) {
    dir.create(pth, showWarnings = FALSE, recursive = TRUE)
  }
  invisible(normalizePath(paths, winslash = "/", mustWork = FALSE))
}


#' Read a GTEx v12 style GCT file into a matrix.
#'
#' The file is expected to have the first two rows as header metadata that
#' specify the number of rows and columns. The function reads the remaining
#' rows using `readr::read_tsv()` and returns both the expression matrix and
#' annotation vectors to keep gene identifiers aligned.
#'
#' @param path_gct Path to a GCT file, optionally gzipped.
#' @return A list with elements `counts`, `gene_id`, `gene_symbol`, and
#'   `sample_ids`.
#' @examples
#' \dontrun{
#' gct <- read_gct_v12("GTExDatav10/example.gct.gz")
#' str(gct$counts)
#' }
read_gct_v12 <- function(path_gct) {
  stopifnot(file.exists(path_gct))
  header_lines <- readr::read_lines(path_gct, n_max = 2)
  if (length(header_lines) < 2) {
    stop("GCT header missing or unreadable: ", path_gct)
  }
  df <- readr::read_tsv(
    path_gct,
    skip = 2,
    col_types = readr::cols(
      .default = readr::col_double(),
      Name = readr::col_character(),
      Description = readr::col_character()
    )
  )
  stopifnot(all(c("Name", "Description") %in% names(df)))
  gene_id <- df$Name
  gene_symbol <- df$Description
  counts <- as.matrix(df[, -(1:2)])
  storage.mode(counts) <- "double"
  rownames(counts) <- gene_id
  list(
    counts = counts,
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    sample_ids = colnames(counts)
  )
}


#' Collapse duplicated matrix rows using row-wise sums.
#'
#' Many GTEx inputs contain duplicated Ensembl IDs due to alternative mappings.
#' This helper collapses them by summing rows with identical row names.
#'
#' @param mat Numeric matrix with row names.
#' @return Matrix with duplicated row names collapsed and ordered.
#' @examples
#' mat <- matrix(1:6, nrow = 3, dimnames = list(c("A", "A", "B"), NULL))
#' collapse_dupe_rowsum(mat)
collapse_dupe_rowsum <- function(mat) {
  stopifnot(!is.null(rownames(mat)))
  if (!any(duplicated(rownames(mat)))) {
    return(mat)
  }
  aggregated <- rowsum(mat, group = rownames(mat))
  aggregated[order(rownames(aggregated)), , drop = FALSE]
}

## -----------------------------------------------------------------------------
## Expression helpers
## -----------------------------------------------------------------------------

#' Map Ensembl IDs to gene symbols and collapse duplicates by mean.
#'
#' GTEx GCT inputs provide Ensembl IDs (Name) and a "Description" column that is
#' typically the gene symbol. This helper maps matrix rownames (gene IDs) to
#' symbols using the provided vectors and collapses duplicated symbols by mean
#' (sum divided by the number of rows per symbol).
#'
#' @param mat Numeric matrix (genes x samples) with rownames as gene IDs.
#' @param gene_id Character vector of gene IDs aligned to `gene_symbol`.
#' @param gene_symbol Character vector of gene symbols aligned to `gene_id`.
#' @param uppercase Logical; uppercase symbols for matching.
#' @return A matrix with rownames as (unique) gene symbols.
collapse_to_symbol_mean <- function(mat, gene_id, gene_symbol, uppercase = TRUE) {
  stopifnot(!is.null(rownames(mat)))
  stopifnot(is.character(gene_id), is.character(gene_symbol))
  sym <- gene_symbol[match(rownames(mat), gene_id)]
  if (uppercase) {
    sym <- toupper(sym)
  }
  keep <- !is.na(sym) & nzchar(sym)
  mat <- mat[keep, , drop = FALSE]
  sym <- sym[keep]
  if (!nrow(mat)) {
    stop("No rows remain after mapping gene IDs to symbols.")
  }
  if (any(duplicated(sym))) {
    summed <- rowsum(mat, group = sym, reorder = FALSE)
    denom <- as.numeric(table(sym)[rownames(summed)])
    mat <- summed / denom
    rownames(mat) <- rownames(summed)
  } else {
    rownames(mat) <- sym
  }
  mat
}


#' Cached lookup of valid human Ensembl gene identifiers.
#'
#' The lookup uses `org.Hs.eg.db` to determine if an Ensembl identifier is
#' present in the annotation database. Identifiers are cleaned by removing
#' version suffixes and `_PAR_` extensions before matching.
#'
#' @param gene_ids Character vector of gene identifiers.
#' @return Logical vector indicating membership in `org.Hs.eg.db`.
#' @examples
#' is_valid_human_gene(c("ENSG00000121410", "foo"))
is_valid_human_gene <- local({
  cache <- NULL
  function(gene_ids) {
    if (is.null(cache)) {
      cache <<- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, keytype = "ENSEMBL")
    }
    cleaned <- gene_ids %>%
      stringr::str_remove("\\..*$") %>%
      stringr::str_remove("_PAR_.*$")
    cleaned %in% cache
  }
})


#' Extract the numeric portion of Ensembl gene identifiers.
#'
#' @param gene_id Character vector of Ensembl IDs (optionally with version).
#' @return Numeric vector with the integer portion of each identifier.
#' @examples
#' get_gene_number(c("ENSG00000121410.12", "ENSG00000000003"))
get_gene_number <- function(gene_id) {
  gene_id %>%
    sub("\\..*$", "", .) %>%
    sub("^ENSG", "", .) %>%
    sub("^0+", "", .) %>%
    suppressWarnings(as.numeric(.))
}


#' Choose a single gene identifier per symbol using heuristic rules.
#'
#' The function implements a deterministic cascade of rules:
#' 1. Drop `_PAR_` entries when mixed with non-PAR IDs.
#' 2. Prefer identifiers present in `org.Hs.eg.db`.
#' 3. Prefer the smallest numeric Ensembl identifier.
#' 4. Use total counts across all samples as a final tie-breaker.
#'
#' @param df Tibble with at least `gene_symbol` and `gene_id`.
#' @param count_matrix Numeric matrix of counts used for tie-breaking. Rows
#'   must be named by Ensembl identifiers.
#' @return Tibble with a single column `gene_id`.
#' @examples
#' library(tibble)
#' counts <- matrix(runif(6), nrow = 3, dimnames = list(c("ENSG1", "ENSG2", "ENSG3"), letters[1:2]))
#' select_one_gene_id_anyN(
#'   tibble(gene_symbol = "X", gene_id = c("ENSG1", "ENSG2")),
#'   count_matrix = counts
#' )
select_one_gene_id_anyN <- function(df, count_matrix) {
  stopifnot("gene_id" %in% names(df))
  ids <- df$gene_id
  if (length(ids) == 1) {
    return(tibble::tibble(gene_id = ids))
  }
  is_par <- stringr::str_detect(ids, "_PAR_")
  valid <- is_valid_human_gene(ids)
  nums <- get_gene_number(ids)
  candidates <- seq_along(ids)

  if (any(!is_par) && any(is_par)) {
    candidates <- candidates[!is_par[candidates]]
  }
  if (any(valid[candidates])) {
    candidates <- candidates[valid[candidates]]
  }
  available_nums <- nums[candidates]
  if (any(!is.na(available_nums))) {
    candidates <- candidates[which.min(available_nums)]
    return(tibble::tibble(gene_id = ids[candidates]))
  }
  if (!is.null(count_matrix)) {
    idx <- match(ids[candidates], rownames(count_matrix))
    totals <- purrr::map_dbl(idx, ~ if (is.na(.x)) NA_real_ else sum(count_matrix[.x, ], na.rm = TRUE))
    best <- candidates[which.max(totals)]
    return(tibble::tibble(gene_id = ids[best]))
  }
  tibble::tibble(gene_id = ids[1])
}


#' Add numeric age columns derived from AGE or AGE_BIN.
#'
#' @param df Data frame containing at least AGE and/or AGE_BIN columns.
#' @return The input data frame with `age_years` appended when possible.
#' @examples
#' derive_age_numeric(tibble::tibble(AGE_BIN = c("40-49", "80+")))
derive_age_numeric <- function(df) {
  if ("AGE" %in% names(df) && is.numeric(df$AGE)) {
    df$age_years <- df$AGE
    return(df)
  }
  if ("AGE_BIN" %in% names(df)) {
    to_mid <- function(bin) {
      if (is.na(bin)) {
        return(NA_real_)
      }
      match <- stringr::str_match(bin, "([0-9]+)\\D+([0-9]+)")
      if (!any(is.na(match[1, 2:3]))) {
        return(mean(as.numeric(match[1, 2:3])))
      }
      if (grepl("\\+", bin)) {
        return(as.numeric(gsub("\\D", "", bin)) + 5)
      }
      NA_real_
    }
    df$age_years <- vapply(df$AGE_BIN, to_mid, numeric(1))
  } else {
    df$age_years <- NA_real_
  }
  df
}


#' Parse textual age bins into minimum, maximum, and midpoint values.
#'
#' Handles common GTEx encodings such as `"50-59"`, `"60 – 69"`, and `"80+"`.
#'
#' @param x Character vector of age bin labels.
#' @return Tibble with columns `age_min`, `age_max`, and `age_mid`.
#' @examples
#' parse_age_bin(c("50-59", "80+"))
parse_age_bin <- function(x) {
  x <- as.character(x)
  age_min <- rep(NA_real_, length(x))
  age_max <- rep(NA_real_, length(x))
  range_match <- stringr::str_match(x, "\\s*([0-9]{1,3})\\D+([0-9]{1,3})\\s*")
  has_range <- !is.na(range_match[, 2]) & !is.na(range_match[, 3])
  age_min[has_range] <- as.numeric(range_match[has_range, 2])
  age_max[has_range] <- as.numeric(range_match[has_range, 3])
  plus_idx <- grepl("^[^0-9]*([0-9]{1,3})\\s*\\+$", x)
  if (any(plus_idx)) {
    base <- as.numeric(gsub("\\D", "", x[plus_idx]))
    age_min[plus_idx] <- base
    age_max[plus_idx] <- base + 10
  }
  tibble::tibble(
    age_min = age_min,
    age_max = age_max,
    age_mid = ifelse(is.na(age_min) | is.na(age_max), NA_real_, (age_min + age_max) / 2)
  )
}


#' Compute PCA on variance-stabilized counts.
#'
#' @param vst_mat Numeric matrix (genes x samples).
#' @param top_n Number of highest variance genes to retain.
#' @return A `prcomp` object.
#' @examples
#' set.seed(1)
#' mat <- matrix(rnorm(1000), nrow = 100)
#' pca_compute(mat, top_n = 50)
pca_compute <- function(vst_mat, top_n = 5000) {
  stopifnot(is.matrix(vst_mat))
  selected <- head(order(matrixStats::rowVars(vst_mat), decreasing = TRUE), top_n)
  stats::prcomp(t(vst_mat[selected, , drop = FALSE]), center = TRUE, scale. = FALSE)
}


#' Join PCA scores with metadata for plotting.
#'
#' @param pca_obj Result from [pca_compute()].
#' @param meta_df Data frame containing sample-level metadata.
#' @param columns Character vector of metadata columns to include.
#' @return Tibble with PCA scores and requested metadata columns.
#' @examples
#' \dontrun{
#' scores <- pca_scores_df(pca, meta, c("DEPOT", "sex"))
#' }
pca_scores_df <- function(pca_obj, meta_df, columns = c("DEPOT", "sex", "age_mid", "age_years")) {
  stopifnot("x" %in% names(pca_obj))
  scores <- tibble::as_tibble(pca_obj$x, rownames = "SAMPID")
  columns <- intersect(columns, names(meta_df))
  if (length(columns)) {
    scores <- dplyr::left_join(scores, meta_df %>% dplyr::select(dplyr::any_of(c("SAMPID", columns))), by = "SAMPID")
  }
  scores
}


#' Save a PCA scatter plot colored by a metadata column.
#'
#' @param pca_df Tibble from [pca_scores_df()].
#' @param varpct Numeric vector of variance explained percentages.
#' @param color_var Column name used for color mapping.
#' @param title_suffix Character suffix for the plot title.
#' @param outfile_png Path to the PNG file.
#' @param outfile_tsv Optional TSV path to persist the scores.
#' @return Invisibly returns the ggplot object.
#' @examples
#' \dontrun{
#' pca_save_plot(df, varpct, "sex", "(all samples)", "pca.png")
#' }
pca_save_plot <- function(pca_df,
                          varpct,
                          color_var,
                          title_suffix,
                          outfile_png,
                          outfile_tsv = NULL) {
  stopifnot(color_var %in% names(pca_df))
  plt <- ggplot2::ggplot(pca_df, ggplot2::aes(.data$PC1, .data$PC2, shape = .data$DEPOT)) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data[[color_var]]),
      size = 2,
      alpha = 0.9,
      na.rm = TRUE
    ) +
    ggplot2::labs(
      title = glue::glue("GTEx adipose PCA {title_suffix} | color = {color_var}"),
      x = glue::glue("PC1 ({round(varpct[1], 1)}%)"),
      y = glue::glue("PC2 ({round(varpct[2], 1)}%)"),
      color = if (identical(color_var, "age_mid")) "Age (midpoint)" else stringr::str_to_title(color_var)
    ) +
    ggplot2::theme_bw()
  ggplot2::ggsave(outfile_png, plt, width = 7, height = 5, dpi = 300)
  if (!is.null(outfile_tsv)) {
    readr::write_tsv(pca_df, outfile_tsv)
  }
  invisible(plt)
}


#' Summarise PCA principal components against age and sex.
#'
#' @param pca_df Tibble with PCA scores and metadata columns `sex`, `age_mid`,
#'   and optionally `DEPOT`.
#' @return Tibble with correlation and regression statistics per PC.
#' @examples
#' \dontrun{
#' compute_pc_associations(pca_df)
#' }
compute_pc_associations <- function(pca_df) {
  stopifnot(all(c("SAMPID", "sex") %in% names(pca_df)))
  if (!"age_mid" %in% names(pca_df)) {
    pca_df$age_mid <- NA_real_
  }
  pca_df <- pca_df %>%
    dplyr::mutate(
      sex_num = dplyr::case_when(
        sex %in% c("Male", "M", "1") ~ 1,
        sex %in% c("Female", "F", "2") ~ 0,
        TRUE ~ NA_real_
      )
    )
  pc_cols <- grep("^PC[0-9]+$", names(pca_df), value = TRUE)
  res <- purrr::map(pc_cols, function(pc) {
    vec <- pca_df[[pc]]
    age_ok <- !is.na(vec) & !is.na(pca_df$age_mid)
    sex_ok <- !is.na(vec) & !is.na(pca_df$sex_num)
    r_age <- if (sum(age_ok) >= 3) stats::cor(vec[age_ok], pca_df$age_mid[age_ok]) else NA_real_
    p_age <- tryCatch(
      if (sum(age_ok) >= 3) stats::cor.test(vec[age_ok], pca_df$age_mid[age_ok])$p.value else NA_real_,
      error = function(e) NA_real_
    )
    r_sex <- if (sum(sex_ok) >= 3) stats::cor(vec[sex_ok], pca_df$sex_num[sex_ok]) else NA_real_
    p_sex <- tryCatch(
      if (sum(sex_ok) >= 3) stats::cor.test(vec[sex_ok], pca_df$sex_num[sex_ok])$p.value else NA_real_,
      error = function(e) NA_real_
    )
    lm_df <- pca_df %>%
      dplyr::select(dplyr::any_of(c(pc, "age_mid", "sex_num", "DEPOT"))) %>%
      stats::na.omit()
    lm_age_beta <- lm_age_p <- lm_sex_beta <- lm_sex_p <- lm_n <- NA_real_
    if (nrow(lm_df) >= 10) {
      if ("DEPOT" %in% names(lm_df) && dplyr::n_distinct(lm_df$DEPOT) >= 2) {
        fit <- stats::lm(lm_df[[pc]] ~ age_mid + sex_num + DEPOT, data = lm_df)
      } else {
        fit <- stats::lm(lm_df[[pc]] ~ age_mid + sex_num, data = lm_df)
      }
      coefs <- summary(fit)$coefficients
      if ("age_mid" %in% rownames(coefs)) {
        lm_age_beta <- coefs["age_mid", "Estimate"]
        lm_age_p <- coefs["age_mid", "Pr(>|t|)"]
      }
      if ("sex_num" %in% rownames(coefs)) {
        lm_sex_beta <- coefs["sex_num", "Estimate"]
        lm_sex_p <- coefs["sex_num", "Pr(>|t|)"]
      }
      lm_n <- nrow(lm_df)
    }
    tibble::tibble(
      PC = pc,
      n_age = sum(age_ok),
      r_age = r_age,
      p_age = p_age,
      n_sex = sum(sex_ok),
      r_sex = r_sex,
      p_sex = p_sex,
      lm_n = lm_n,
      lm_beta_age_mid = lm_age_beta,
      lm_p_age_mid = lm_age_p,
      lm_beta_sex_num = lm_sex_beta,
      lm_p_sex_num = lm_sex_p
    )
  })
  dplyr::bind_rows(res)
}


#' Save top and bottom PCA loadings for a principal component.
#'
#' @param pca_obj Result from [pca_compute()].
#' @param pc Column name of the principal component (e.g., `"PC1"`).
#' @param n Number of loadings to export from each tail.
#' @param out_path Path to the TSV file to write.
#' @return Invisibly returns the tidy tibble of loadings.
#' @examples
#' \dontrun{
#' save_pc_loadings(pca, "PC1", 100, "pc1.tsv")
#' }
save_pc_loadings <- function(pca_obj, pc, n = 200, out_path) {
  stopifnot(pc %in% colnames(pca_obj$rotation))
  loadings <- pca_obj$rotation[, pc]
  ord <- sort(loadings, decreasing = TRUE)
  top <- head(ord, n)
  bot <- head(rev(ord), n)
  df <- dplyr::bind_rows(
    tibble::tibble(gene_id = names(top), loading = as.numeric(top), tail = "top"),
    tibble::tibble(gene_id = names(bot), loading = as.numeric(bot), tail = "bottom")
  )
  readr::write_tsv(df, out_path)
  invisible(df)
}


#' Plot absolute PCA loadings as a simple bar chart.
#'
#' @param pca_obj Result from [pca_compute()].
#' @param pc Principal component column.
#' @param k Number of top genes to display.
#' @param outfile_png Output PNG path.
#' @return Invisibly returns the ggplot object.
#' @examples
#' \dontrun{
#' plot_top_loadings(pca, "PC1", 20, "pc1_top.png")
#' }
plot_top_loadings <- function(pca_obj, pc, k = 20, outfile_png) {
  loadings <- pca_obj$rotation[, pc]
  sel <- sort(abs(loadings), decreasing = TRUE)
  sel <- sel[seq_len(min(k, length(sel)))]
  df <- tibble::tibble(gene_id = names(sel), abs_loading = as.numeric(sel))
  plt <- ggplot2::ggplot(df, ggplot2::aes(x = stats::reorder(gene_id, abs_loading), y = abs_loading)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = glue::glue("Top |loading| for {pc}"),
      x = "Gene",
      y = "|loading|"
    )
  ggplot2::ggsave(outfile_png, plt, width = 6, height = 7, dpi = 300)
  invisible(plt)
}


#' Scatter PCA scores against age for quick inspection.
#'
#' @param pca_df Tibble with PCA scores, `age_mid`, and `sex`.
#' @param tag Character label used in file naming.
#' @param fig_dir Directory for PNG outputs.
#' @return Invisibly returns a list of ggplot objects.
#' @examples
#' \dontrun{
#' quick_pc_scatter(pca_df, "both", "figs")
#' }
quick_pc_scatter <- function(pca_df, tag, fig_dir) {
  stopifnot("age_mid" %in% names(pca_df))
  plots <- list()
  df <- pca_df %>% dplyr::filter(!is.na(age_mid))
  if (nrow(df) < 10) {
    return(invisible(plots))
  }
  for (pc in c("PC1", "PC2")) {
    if (!pc %in% names(df)) {
      next
    }
    plt <- ggplot2::ggplot(df, ggplot2::aes(age_mid, .data[[pc]], color = sex)) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_smooth(method = "lm", se = FALSE) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = glue::glue("{pc} vs age_mid - {tag}"),
        x = "Age midpoint",
        y = pc
      )
    ggplot2::ggsave(
      filename = file.path(fig_dir, glue::glue("Scatter_{pc}_vs_age_mid_{gsub(' ', '_', tag)}.png")),
      plot = plt,
      width = 6.5,
      height = 4.5,
      dpi = 300
    )
    plots[[pc]] <- plt
  }
  invisible(plots)
}


#' Generate dot plots for selected genes stratified by tissue and sex.
#'
#' @param gene_symbol Gene symbol or Ensembl identifier to plot.
#' @param tissue Single tissue/DEPOT label.
#' @param vst_mat Variance-stabilised expression matrix.
#' @param meta_df Sample metadata with `DEPOT`, `sex`, and age columns.
#' @param fig_dir Output directory for PNG files.
#' @param tab_dir Output directory for CSV files.
#' @return Invisibly returns `NULL`.
#' @examples
#' \dontrun{
#' plot_gene_tissue_dots_from_mats("UCP1", "Subcutaneous", vst, meta, "figs", "tables")
#' }
plot_gene_tissue_dots_from_mats <- function(gene_symbol,
                                            tissue,
                                            vst_mat,
                                            meta_df,
                                            fig_dir,
                                            tab_dir) {
  strip_ver <- function(x) sub("\\.\\d+$", "", x)
  row_ix <- which(rownames(vst_mat) == gene_symbol)
  if (!length(row_ix)) {
    row_ix <- which(strip_ver(rownames(vst_mat)) == strip_ver(gene_symbol))
    if (!length(row_ix)) {
      return(invisible(NULL))
    }
  }
  vec <- as.numeric(vst_mat[row_ix, ])
  df <- tibble::tibble(SAMPID = colnames(vst_mat), vst = vec) %>%
    dplyr::left_join(
      meta_df %>% dplyr::select(dplyr::any_of(c(
        "SAMPID", "DEPOT", "sex", "AGE", "age_bin_label", "age_mid", "age_bin_lo"
      ))),
      by = "SAMPID"
    ) %>%
    dplyr::filter(DEPOT == tissue, !is.na(sex)) %>%
    dplyr::mutate(
      age_bin_lo = dplyr::if_else(is.na(age_bin_lo), Inf, age_bin_lo),
      gene = gene_symbol
    ) %>%
    dplyr::arrange(age_bin_lo, age_mid, vst)

  if (!nrow(df)) {
    return(invisible(NULL))
  }

  df <- df %>%
    dplyr::mutate(
      samp_order = factor(SAMPID, levels = unique(SAMPID[order(age_bin_lo, age_mid, vst)])),
      age_bin_f  = factor(age_bin_label, levels = sort(unique(age_bin_label[order(age_bin_lo)])))
    )

  plt <- ggplot2::ggplot(df, ggplot2::aes(x = samp_order, y = vst, color = age_bin_f)) +
    ggplot2::geom_point(size = 1.8, alpha = 0.85) +
    ggplot2::facet_wrap(~sex, ncol = 1, scales = "free_x") +
    ggplot2::labs(
      title = glue::glue("{gene_symbol} in {tissue}"),
      subtitle = "Samples ordered by age bin, midpoint, and expression.",
      x = "Samples",
      y = "VST expression",
      color = "AGE bin"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92")
    )

  out_png <- file.path(fig_dir, glue::glue("Dots_{gene_symbol}_{gsub(' ', '_', tissue)}.png"))
  out_csv <- file.path(tab_dir, glue::glue("Dots_{gene_symbol}_{gsub(' ', '_', tissue)}.csv"))
  ggplot2::ggsave(out_png, plt, width = 12, height = 5.5, dpi = 300)
  readr::write_csv(
    df %>% dplyr::select(dplyr::any_of(c(
      "SAMPID", "DEPOT", "sex", "AGE", "age_bin_label", "age_mid", "gene", "vst"
    ))),
    out_csv
  )
  invisible(NULL)
}
