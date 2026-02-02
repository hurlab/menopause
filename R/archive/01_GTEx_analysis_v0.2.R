################################################################################
## GTEx v10 Adipose analysis: PCA first, then Female DE (<45 vs >55)
## Outputs: GTEx_v10_AT_analysis_out/{figs,tables,rds}
################################################################################

## -----------------------------------------------------------------------------
## Setup
## -----------------------------------------------------------------------------
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  pth <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) NULL)
  if (!is.null(pth) && nzchar(pth)) setwd(dirname(pth))
}
# Ensure a default CRAN mirror is set for non-interactive runs
.repos <- getOption("repos")
.cran <- unname(.repos["CRAN"])
if (is.null(.cran) || is.na(.cran) || .cran == "@CRAN@" || .cran == "") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
rm(.repos, .cran)


if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  data.table, readr, dplyr, tidyr, stringr, tibble, purrr, glue, magrittr,
  ggplot2, patchwork, DESeq2, vsn, Matrix, matrixStats, AnnotationDbi, org.Hs.eg.db
)



## -----------------------------------------------------------------------------
## Paths
## -----------------------------------------------------------------------------
in_dir   <- "GTExDatav10"
gct_sc   <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
gct_vis  <- file.path(in_dir, "gene_reads_v10_adipose_visceral_omentum.gct.gz")
sattr    <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph   <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")  # optional

out_dir  <- "GTEx_v10_AT_analysis_out"
fig_dir  <- file.path(out_dir, "figs")
tab_dir  <- file.path(out_dir, "tables")
rds_dir  <- file.path(out_dir, "rds")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------
read_gct_v12 <- function(path_gct) {
  stopifnot(file.exists(path_gct))
  hdr <- readr::read_lines(path_gct, n_max = 2)
  if (length(hdr) < 2) stop("GCT header missing or unreadable: ", path_gct)
  # Read table skipping first two lines
  dt <- readr::read_tsv(
    path_gct, skip = 2,
    col_types = readr::cols(.default = readr::col_double(),
                            Name = readr::col_character(),
                            Description = readr::col_character())
  )
  stopifnot(all(c("Name","Description") %in% names(dt)))
  gene_id  <- dt$Name
  gene_sym <- dt$Description
  mat      <- as.matrix(dt[, -(1:2)])
  storage.mode(mat) <- "double"
  rownames(mat) <- gene_id
  list(counts = mat, gene_id = gene_id, gene_symbol = gene_sym, sample_ids = colnames(mat))
}

collapse_dupe_rowsum <- function(mat) {
  if (!any(duplicated(rownames(mat)))) return(mat)
  ag <- rowsum(mat, group = rownames(mat))
  ag[order(rownames(ag)), , drop = FALSE]
}

derive_age_numeric <- function(df) {
  # Prefer numeric AGE if present
  if ("AGE" %in% names(df) && is.numeric(df$AGE)) {
    df$age_years <- df$AGE
    return(df)
  }
  # Fallback to AGE_BIN like "20-29" or "80+"
  to_mid <- function(bin) {
    if (is.na(bin)) return(NA_real_)
    m <- stringr::str_match(bin, "([0-9]+)\\D+([0-9]+)")
    if (!any(is.na(m[1, 2:3]))) return(mean(as.numeric(m[1, 2:3])))
    if (grepl("\\+", bin)) return(as.numeric(gsub("\\D", "", bin)) + 5)
    return(NA_real_)
  }
  if ("AGE_BIN" %in% names(df)) {
    df$age_years <- vapply(df$AGE_BIN, to_mid, numeric(1))
  } else {
    df$age_years <- NA_real_
  }
  df
}


## -----------------------------------------------------------------------------
## Load counts
## -----------------------------------------------------------------------------
x_sc  <- read_gct_v12(gct_sc)
x_vis <- read_gct_v12(gct_vis)

x_sc$counts  <- collapse_dupe_rowsum(x_sc$counts)
x_vis$counts <- collapse_dupe_rowsum(x_vis$counts)

common_genes <- intersect(rownames(x_sc$counts), rownames(x_vis$counts))
x_sc$counts  <- x_sc$counts[common_genes, , drop = FALSE]
x_vis$counts <- x_vis$counts[common_genes, , drop = FALSE]

# Check the dimensions and common genes
message("Subcutaneous: ", nrow(x_sc$counts), " genes x ", ncol(x_sc$counts), " samples")
message("Visceral:    ", nrow(x_vis$counts), " genes x ", ncol(x_vis$counts), " samples")
message("Common genes: ", length(common_genes))

# Combined counts for deduplication logic
counts_all <- cbind(x_sc$counts, x_vis$counts)

## Check the gene annotation
gene_info_df <- data.frame(
  gene_id = x_sc$gene_id,
  gene_symbol = x_sc$gene_symbol
)

valid_ensg_keys <- keys(org.Hs.eg.db, keytype = "ENSEMBL")

is_valid_human_gene <- function(gene_ids) {
  cleaned <- gene_ids %>%
    str_remove("\\..*$") %>%       # drop version
    str_remove("_PAR_.*$") %>%     # drop PAR suffix
    unique()
  cleaned %in% valid_ensg_keys
}

get_gene_number <- function(gene_id) {
  gene_id %>%
    sub("\\..*$", "", .) %>%       # drop version
    sub("^ENSG", "", .) %>%        # drop prefix
    sub("^0+", "", .) %>%          # drop leading zeros
    suppressWarnings(as.numeric(.))
}

select_gene_id <- function(df) {
  # df has columns gene_id, gene_symbol (and maybe others)
  ids <- df$gene_id
  is_par <- str_detect(ids, "_PAR_")
  valid  <- is_valid_human_gene(ids)
  
  # Rule 1: prefer the one without _PAR_
  if (length(unique(is_par)) == 2) {
    return(tibble(gene_id = ids[!is_par][1]))
  }
  
  # Rule 2: if both valid, pick the smaller numeric ENSG
  if (all(valid)) {
    nums <- get_gene_number(ids)
    if (all(!is.na(nums))) {
      return(tibble(gene_id = ids[which.min(nums)]))
    }
  }
  
  # Rule 3: if only one is valid, keep it
  if (length(unique(valid)) == 2) {
    return(tibble(gene_id = ids[valid][1]))
  }
  
  # Rule 4: otherwise return empty
  tibble(gene_id = character())
}



## -----------------------------------------------------------------------------
## Deduplicate gene symbols: keep exactly one gene_id per symbol
## Strategy:
##  1) Prefer IDs without "_PAR_"
##  2) Prefer valid human Ensembl IDs (org.Hs.eg.db)
##  3) Prefer smaller numeric ENSG (older canonical)
##  4) Tie-breaker: higher total raw counts across all samples
## -----------------------------------------------------------------------------

# Helper: total raw counts for a given gene_id using the already built counts_all
total_counts_for_ids <- function(ids) {
  idx <- match(ids, rownames(counts_all))
  sapply(idx, function(i) if (is.na(i)) NA_real_ else sum(counts_all[i, ], na.rm = TRUE))
}

# Generalized selector that can handle 1, 2, or many IDs for the same symbol
select_one_gene_id_anyN <- function(df) {
  # df has at least: gene_id, gene_symbol
  ids <- df$gene_id
  if (length(ids) == 1) return(tibble(gene_id = ids))  # nothing to dedup
  
  is_par <- stringr::str_detect(ids, "_PAR_")
  valid  <- is_valid_human_gene(ids)
  nums   <- get_gene_number(ids)
  
  # Start with all candidates
  cand <- seq_along(ids)
  
  # Rule 1: drop PAR if mixture of PAR and non-PAR
  if (any(!is_par) && any(is_par)) cand <- cand[!is_par[cand]]
  
  # Rule 2: if any valid remain, drop invalid
  if (any(valid[cand])) cand <- cand[valid[cand]]
  
  # Rule 3: if ENSG numbers available for remaining, keep smallest
  if (sum(!is.na(nums[cand])) >= 1) {
    cand_num <- cand[!is.na(nums[cand])]
    if (length(cand_num) >= 1) {
      cand <- cand_num[which.min(nums[cand_num])]
      return(tibble(gene_id = ids[cand]))
    }
  }
  
  # Rule 4: tie-break on total counts
  tc <- total_counts_for_ids(ids[cand])
  cand <- cand[which.max(tc)]
  tibble(gene_id = ids[cand])
}

# Build a full dedup map for ALL symbols (not just the 2-ID ones)
dedup_map <- gene_info_df %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::group_modify(~ select_one_gene_id_anyN(.x)) %>%  # returns tibble(gene_id)
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(gene_id)) %>%              # in case a group returned empty
  dplyr::distinct(gene_symbol, gene_id)           # safety: one row per symbol

# What is kept vs. dropped
kept_ids <- dedup_map$gene_id
all_ids  <- gene_info_df$gene_id
dropped_tbl <- gene_info_df %>%
  dplyr::filter(!gene_id %in% kept_ids) %>%
  dplyr::arrange(gene_symbol, gene_id)

readr::write_csv(dedup_map,   file.path(tab_dir, "dedup_kept_gene_symbol_to_gene_id.csv"))
readr::write_csv(dropped_tbl, file.path(tab_dir, "dedup_dropped_gene_ids.csv"))

# Apply the filter to BOTH matrices and their annotations
keep_sc <- rownames(x_sc$counts)  %in% kept_ids
keep_vi <- rownames(x_vis$counts) %in% kept_ids

x_sc$counts      <- x_sc$counts[keep_sc, , drop = FALSE]
x_sc$gene_id     <- x_sc$gene_id[keep_sc] %>% stringr::str_remove("\\..*$")
x_sc$gene_symbol <- x_sc$gene_symbol[keep_sc]

x_vis$counts      <- x_vis$counts[keep_vi, , drop = FALSE]
x_vis$gene_id     <- x_vis$gene_id[keep_vi] %>% stringr::str_remove("\\..*$")
x_vis$gene_symbol <- x_vis$gene_symbol[keep_vi]

# sanity check: after filtering, ensure the two matrices still share the same rows
common_genes <- intersect(rownames(x_sc$counts), rownames(x_vis$counts))
x_sc$counts  <- x_sc$counts[common_genes, , drop = FALSE]
x_vis$counts <- x_vis$counts[common_genes, , drop = FALSE]

# Rebuild counts_all to reflect deduplication
counts_all <- cbind(x_sc$counts, x_vis$counts)

message("Dedup complete: kept ", length(kept_ids), " gene_ids across ",
        dplyr::n_distinct(dedup_map$gene_symbol), " unique symbols. ",
        "Dropped ", sum(!all_ids %in% kept_ids), " gene_ids.")


## x_vis$gene_id   Remove '.number' and keep the before '.' part only
## get the unique list. 
## see if any gene ID without the .number (version) are duplicated.

## x_vis$gene_id   Remove '.number' and keep the before '.' part only
## get the unique list.
## see if any gene ID without the .number (version) are duplicated.

# 1. Remove '.number' part from gene IDs in x_vis
x_vis$gene_id_clean <- stringr::str_remove(x_vis$gene_id, "\\..*$")

# 2. Update the gene_symbol and gene_id fields in x_vis (optional, but good practice if needed later)
# Note: Rownames typically remain the original full gene_id for tracking,
# but the gene_id column can be updated to the cleaned version if the original is stored elsewhere.
# For simplicity, we'll just check the clean IDs and create a new column.

# 3. Get the unique list of cleaned gene IDs
unique_clean_ids <- unique(x_vis$gene_id_clean)

# 4. Check if any gene ID without the version are duplicated (i.e., multiple *original* IDs map to the same *cleaned* ID)
duplicated_clean_ids <- x_vis$gene_id_clean[duplicated(x_vis$gene_id_clean)]

if (length(duplicated_clean_ids) > 0) {
  message("Warning: ", length(duplicated_clean_ids), " duplicated cleaned gene IDs found in x_vis.")
  # Print examples of duplicates for inspection
  message("Examples of duplicated clean IDs: ", paste(head(unique(duplicated_clean_ids)), collapse = ", "))
  
  # Optional: Inspect the original IDs that map to one of the duplicated clean IDs
  # first_dup_clean_id <- unique(duplicated_clean_ids)[1]
  # message("Original IDs for the first duplicate (", first_dup_clean_id, "): ",
  #         paste(x_vis$gene_id[x_vis$gene_id_clean == first_dup_clean_id], collapse = ", "))
  
} else {
  message("No duplicated cleaned gene IDs found in x_vis.")
}

















## -----------------------------------------------------------------------------
## Build base sample frame with DEPOT
## -----------------------------------------------------------------------------
sample_df <- tibble(
  SAMPID = colnames(counts_all),
  DEPOT  = if_else(SAMPID %in% x_sc$sample_ids, "Subcutaneous", "Visceral"),
  SUBJID = sub("^([A-Z0-9]+\\-[A-Z0-9]+).*", "\\\\1", SAMPID)
)

## -----------------------------------------------------------------------------
## Load sample attributes and join
## -----------------------------------------------------------------------------
stopifnot(file.exists(sattr))
sa <- data.table::fread(sattr, sep = "\t", header = TRUE, data.table = FALSE) %>% as_tibble() %>%
  distinct(SAMPID, .keep_all = TRUE)

meta <- sample_df %>% left_join(sa, by = "SAMPID")

# Optional subject phenotypes
if (file.exists(subjph)) {
  sp <- data.table::fread(subjph, sep = "\t", header = TRUE, data.table = FALSE) %>% as_tibble()
  names(sp) <- toupper(names(sp))
  meta <- meta %>% mutate(SUBJID = sub("^([A-Z0-9]+\\-[A-Z0-9]+).*", "\\\\1", SAMPID))
  meta <- meta %>%
    left_join(sp %>% dplyr::select(SUBJID, SEX, AGE, AGE_BIN = dplyr::any_of("AGE_BIN")), by = "SUBJID")
}

# Harmonize names
names(meta) <- gsub("\\.", "_", names(meta))

# Sex
if (!"SEX" %in% names(meta) && "SMGENDER" %in% names(meta)) meta$SEX <- meta$SMGENDER
meta$sex <- dplyr::case_when(
  as.character(meta$SEX) %in% c("1","M","Male") ~ "Male",
  as.character(meta$SEX) %in% c("2","F","Female") ~ "Female",
  TRUE ~ NA_character_
)

# Age numeric
meta <- derive_age_numeric(meta)

# Keep only samples with counts
meta <- meta %>% filter(SAMPID %in% colnames(counts_all))
counts_all <- counts_all[, meta$SAMPID, drop = FALSE]


readr::write_tsv(meta, file.path(tab_dir, "sample_metadata_full.tsv"))




## -----------------------------------------------------------------------------
## 1) Compute age_min, age_max, age_mid for ALL samples
## -----------------------------------------------------------------------------
parse_age_bin <- function(x) {
  x <- as.character(x)
  age_min <- rep(NA_real_, length(x))
  age_max <- rep(NA_real_, length(x))
  # "A-B" formats like "50-59", "60 – 69"
  m <- stringr::str_match(x, "\\s*([0-9]{1,3})\\D+([0-9]{1,3})\\s*")
  has_range <- !is.na(m[,2]) & !is.na(m[,3])
  age_min[has_range] <- as.numeric(m[has_range, 2])
  age_max[has_range] <- as.numeric(m[has_range, 3])
  # "A+" formats like "80+"
  plus_idx <- grepl("^[^0-9]*([0-9]{1,3})\\s*\\+$", x)
  if (any(plus_idx)) {
    base <- as.numeric(gsub("\\D", "", x[plus_idx]))
    age_min[plus_idx] <- base
    age_max[plus_idx] <- base + 10   # assume a decade width to get a midpoint
  }
  tibble(
    age_min = age_min,
    age_max = age_max,
    age_mid = ifelse(is.na(age_min) | is.na(age_max), NA_real_, (age_min + age_max) / 2)
  )
}

# Always compute midpoints from AGE or AGE_BIN, then coalesce into age_years for downstream
src_age <- if ("AGE" %in% names(meta)) meta$AGE else if ("AGE_BIN" %in% names(meta)) meta$AGE_BIN else NA
ab <- parse_age_bin(src_age)
meta$age_min  <- ab$age_min
meta$age_max  <- ab$age_max
meta$age_mid  <- ab$age_mid
meta$age_years <- dplyr::coalesce(meta$age_years, meta$age_mid)

readr::write_tsv(meta, file.path(tab_dir, "sample_metadata_with_age_mid.tsv"))

## -----------------------------------------------------------------------------
## 2) Build VST on ALL samples (no exclusions) for PCA
## -----------------------------------------------------------------------------
dds_all <- DESeqDataSetFromMatrix(
  countData = round(counts_all),
  colData   = as.data.frame(meta),
  design    = ~ 1
)
keep_all <- rowSums(counts(dds_all) >= 10) >= 20
dds_all  <- dds_all[keep_all, ]
dds_all  <- estimateSizeFactors(dds_all)
vst_all  <- assay(vst(dds_all, blind = TRUE))

## -----------------------------------------------------------------------------
## 3) PCA helpers and plots: color by sex and by age_mid
## -----------------------------------------------------------------------------
.pca_compute <- function(vst_mat) {
  top_var <- head(order(matrixStats::rowVars(vst_mat), decreasing = TRUE), 5000)
  prcomp(t(vst_mat[top_var, , drop = FALSE]), center = TRUE, scale. = FALSE)
}
.pca_scores_df <- function(pca_obj, meta_df) {
  as_tibble(pca_obj$x, rownames = "SAMPID") %>%
    left_join(meta_df %>% dplyr::select(SAMPID, DEPOT, sex, age_years, age_mid), by = "SAMPID")
}
.pca_save_plot <- function(pca_df, varpct, color_var, title_suffix, outfile_png, outfile_tsv = NULL) {
  # color_var should be "sex" or "age_mid"
  if (!color_var %in% names(pca_df)) {
    stop("Column '", color_var, "' not found in pca_df")
  }
  # Make sure age_mid is numeric if used
  if (color_var == "age_mid") {
    pca_df <- pca_df %>% mutate(age_mid = as.numeric(age_mid))
  }
  
  p <- ggplot(pca_df, aes(PC1, PC2, shape = DEPOT)) +
    geom_point(aes(color = .data[[color_var]]), size = 2, alpha = 0.9, na.rm = TRUE) +
    labs(
      title = paste("GTEx adipose PCA", title_suffix, "| color =", color_var),
      x = paste0("PC1 (", round(varpct[1], 1), "%)"),
      y = paste0("PC2 (", round(varpct[2], 1), "%)"),
      color = if (color_var == "age_mid") "Age (midpoint)" else "Sex"
    ) +
    theme_bw()
  
  ggsave(outfile_png, p, width = 7, height = 5, dpi = 300)
  
  if (!is.null(outfile_tsv)) {
    readr::write_tsv(pca_df, outfile_tsv)
  }
  invisible(p)
}


# Compute PCA once on ALL samples
pca_all <- .pca_compute(vst_all)
varpct_all <- (pca_all$sdev^2) / sum(pca_all$sdev^2) * 100
pca_both_df <- .pca_scores_df(pca_all, meta)

# Save objects for later correlation work
pca_both <- list(pca = pca_all, pca_df = pca_both_df, var_explained = varpct_all)

# 3a) ALL samples, color by SEX
.pca_save_plot(
  pca_both_df, varpct_all, color_var = "sex",
  title_suffix = "(all samples, both depots)",
  outfile_png  = file.path(fig_dir, "PCA_PC1_PC2_ALL_both_DEPOT_COLOR_BY_SEX.png"),
  outfile_tsv  = file.path(tab_dir, "PCA_scores_ALL_both_DEPOT_COLOR_BY_SEX.tsv")
)

# 3b) ALL samples, color by AGE midpoint
.pca_save_plot(
  pca_both_df, varpct_all, color_var = "age_mid",
  title_suffix = "(all samples, both depots)",
  outfile_png  = file.path(fig_dir, "PCA_PC1_PC2_ALL_both_DEPOT_COLOR_BY_AGE.png"),
  outfile_tsv  = file.path(tab_dir, "PCA_scores_ALL_both_DEPOT_COLOR_BY_AGE.tsv")
)

# 3c) Depot-specific plots reusing the same PCA rotation (subsetting scores keeps PCs comparable)
for (dep in c("Subcutaneous", "Visceral")) {
  df_dep <- pca_both_df %>% dplyr::filter(DEPOT == dep)
  if (nrow(df_dep) < 10) {
    message("Skipping PCA plots for ", dep, " due to small sample size.")
    next
  }
  .pca_save_plot(
    df_dep, varpct_all, color_var = "sex",
    title_suffix = paste0("(", dep, " only)"),
    outfile_png  = file.path(fig_dir, paste0("PCA_PC1_PC2_ALL_", gsub(" ", "_", dep), "_COLOR_BY_SEX.png")),
    outfile_tsv  = file.path(tab_dir,  paste0("PCA_scores_ALL_", gsub(" ", "_", dep), "_COLOR_BY_SEX.tsv"))
  )
  .pca_save_plot(
    df_dep, varpct_all, color_var = "age_mid",
    title_suffix = paste0("(", dep, " only)"),
    outfile_png  = file.path(fig_dir, paste0("PCA_PC1_PC2_ALL_", gsub(" ", "_", dep), "_COLOR_BY_AGE.png")),
    outfile_tsv  = file.path(tab_dir,  paste0("PCA_scores_ALL_", gsub(" ", "_", dep), "_COLOR_BY_AGE.tsv"))
  )
}

## -----------------------------------------------------------------------------
## 4) PC associations with age and sex (all samples + depot subsets)
## -----------------------------------------------------------------------------
compute_pc_associations <- function(pca_df, out_prefix) {
  stopifnot(all(c("SAMPID","DEPOT","sex") %in% names(pca_df)))
  
  # Ensure age_mid exists
  if (!"age_mid" %in% names(pca_df)) {
    if ("age_years" %in% names(pca_df)) pca_df$age_mid <- pca_df$age_years else pca_df$age_mid <- NA_real_
  }
  
  # Clean coding and drop unused levels
  pca_df <- pca_df %>%
    mutate(
      sex_num = case_when(
        sex %in% c("Male","M","1") ~ 1,
        sex %in% c("Female","F","2") ~ 0,
        TRUE ~ NA_real_
      ),
      DEPOT = droplevels(factor(DEPOT))
    )
  
  pc_cols <- grep("^PC[0-9]+$", names(pca_df), value = TRUE)
  if (length(pc_cols) == 0) stop("No PC columns found in pca_df")
  
  rows <- lapply(pc_cols, function(pc) {
    vec <- pca_df[[pc]]
    
    # Age correlation
    age_ok <- !is.na(vec) & !is.na(pca_df$age_mid)
    r_age <- if (sum(age_ok) >= 3) cor(vec[age_ok], pca_df$age_mid[age_ok], method = "pearson") else NA_real_
    p_age <- tryCatch(if (sum(age_ok) >= 3) cor.test(vec[age_ok], pca_df$age_mid[age_ok])$p.value else NA_real_, error = function(e) NA_real_)
    
    # Sex correlation (point biserial via 0-1 coding)
    sex_ok <- !is.na(vec) & !is.na(pca_df$sex_num)
    r_sex <- if (sum(sex_ok) >= 3) cor(vec[sex_ok], pca_df$sex_num[sex_ok], method = "pearson") else NA_real_
    p_sex <- tryCatch(if (sum(sex_ok) >= 3) cor.test(vec[sex_ok], pca_df$sex_num[sex_ok])$p.value else NA_real_, error = function(e) NA_real_)
    
    # Linear model with optional DEPOT
    fit_df <- pca_df %>% dplyr::select(!!pc, age_mid, sex_num, DEPOT) %>% na.omit()
    fit_df$DEPOT <- droplevels(fit_df$DEPOT)
    
    lm_age_beta <- lm_age_p <- lm_sex_beta <- lm_sex_p <- lm_n <- NA_real_
    
    if (nrow(fit_df) >= 10) {
      if (nlevels(fit_df$DEPOT) >= 2) {
        fit <- lm(fit_df[[pc]] ~ age_mid + sex_num + DEPOT, data = fit_df)
      } else {
        fit <- lm(fit_df[[pc]] ~ age_mid + sex_num, data = fit_df)
      }
      sm <- summary(fit)$coefficients
      if ("age_mid" %in% rownames(sm)) { lm_age_beta <- sm["age_mid","Estimate"]; lm_age_p <- sm["age_mid","Pr(>|t|)"] }
      if ("sex_num" %in% rownames(sm)) { lm_sex_beta <- sm["sex_num","Estimate"]; lm_sex_p <- sm["sex_num","Pr(>|t|)"] }
      lm_n <- nrow(fit_df)
    }
    
    tibble(
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
  
  out <- bind_rows(rows)
  readr::write_tsv(out, file.path(tab_dir, paste0(out_prefix, "_PC_age_sex_associations.tsv")))
  invisible(out)
}


# All samples
assoc_all <- compute_pc_associations(pca_both$pca_df, "ALL_both_DEPOT")

# Subsets in the same PCA space
assoc_sc  <- compute_pc_associations(pca_both$pca_df %>% dplyr::filter(DEPOT == "Subcutaneous"), "ALL_Subcutaneous")
assoc_vis <- compute_pc_associations(pca_both$pca_df %>% dplyr::filter(DEPOT == "Visceral"),      "ALL_Visceral")

assoc_combined <- bind_rows(
  assoc_all %>% mutate(scope = "both_depots"),
  assoc_sc  %>% mutate(scope = "subcutaneous"),
  assoc_vis %>% mutate(scope = "visceral")
)
readr::write_csv(assoc_combined, file.path(tab_dir, "PC_age_sex_associations_combined.csv"))


## -----------------------------------------------------------------------------
## 5) Optional quick visuals: PC1, PC2 vs age_mid, colored by sex
## -----------------------------------------------------------------------------
quick_pc_scatter <- function(pca_df, tag) {
  df <- pca_df %>% dplyr::filter(!is.na(age_mid))
  if (nrow(df) < 10) return(invisible(NULL))
  for (pc in c("PC1","PC2")) {
    if (!pc %in% names(df)) next
    p <- ggplot(df, aes(age_mid, .data[[pc]], color = sex)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE) +
      theme_bw() +
      labs(title = paste0(pc, " vs age_mid - ", tag), x = "Age midpoint", y = pc)
    ggsave(file.path(fig_dir, paste0("Scatter_", pc, "_vs_age_mid_", gsub(" ", "_", tag), ".png")),
           p, width = 6.5, height = 4.5, dpi = 300)
  }
}
quick_pc_scatter(pca_both$pca_df, "both_depots")
quick_pc_scatter(pca_both$pca_df %>% dplyr::filter(DEPOT == "Subcutaneous"), "Subcutaneous")
quick_pc_scatter(pca_both$pca_df %>% dplyr::filter(DEPOT == "Visceral"), "Visceral")



## -----------------------------------------------------------------------------
## 6) PC loadings and quick pathway hints
## -----------------------------------------------------------------------------

# Save top and bottom loadings for PC1 and PC2
save_pc_loadings <- function(pca_obj, pc, n = 200, out_prefix = "PC") {
  stopifnot(pc %in% colnames(pca_obj$rotation))
  rot <- pca_obj$rotation[, pc]
  ord <- sort(rot, decreasing = TRUE)
  top  <- head(ord, n)
  bot  <- head(rev(ord), n)
  top_tbl <- tibble(gene_id = names(top), loading = as.numeric(top), tail = "top")
  bot_tbl <- tibble(gene_id = names(bot), loading = as.numeric(bot), tail = "bottom")
  out <- bind_rows(top_tbl, bot_tbl)
  readr::write_tsv(out, file.path(tab_dir, paste0(out_prefix, pc, "_loadings_top_bottom.tsv")))
  invisible(out)
}

pc1_load <- save_pc_loadings(pca_both$pca, "PC1", n = 200, out_prefix = "ALL_")
pc2_load <- save_pc_loadings(pca_both$pca, "PC2", n = 200, out_prefix = "ALL_")

# Simple barplots of the top 20 absolute loadings
plot_top_loadings <- function(pca_obj, pc, k = 20, outfile_png) {
  rot <- pca_obj$rotation[, pc]
  sel <- sort(abs(rot), decreasing = TRUE)
  sel <- sel[1:min(k, length(sel))]
  df <- tibble(gene_id = names(sel), abs_loading = as.numeric(sel))
  p <- ggplot(df, aes(x = reorder(gene_id, abs_loading), y = abs_loading)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(title = paste0("Top |loading| for ", pc), x = "Gene", y = "|loading|")
  ggsave(outfile_png, p, width = 6, height = 7, dpi = 300)
}

plot_top_loadings(pca_both$pca, "PC1",
                  outfile_png=file.path(fig_dir, "PC1_top_abs_loadings_bar.png"))
plot_top_loadings(pca_both$pca, "PC2",
                  outfile_png=file.path(fig_dir, "PC2_top_abs_loadings_bar.png"))

# Optional: export a ranked vector for external enrichment tools
pc2_rank <- tibble(gene_id = rownames(pca_both$pca$rotation),
                   loading = as.numeric(pca_both$pca$rotation[, "PC2"])) %>%
  arrange(desc(loading))
readr::write_tsv(pc2_rank, file.path(tab_dir, "PC2_ranked_loadings.tsv"))



## -----------------------------------------------------------------------------
## 7) Check technical covariates vs PCs (RIN, ischemic time, center)
## -----------------------------------------------------------------------------

covar_df <- pca_both$pca_df %>%
  left_join(meta %>% dplyr::select(SAMPID, SMRIN, SMTSISCH, SMCENTER, SMNABTCH, SMGEBTCH),
            by = "SAMPID")

# Pairwise Pearson correlations for numeric covariates
num_covs <- c("SMRIN", "SMTSISCH")
pc_cols  <- grep("^PC[0-9]+$", names(covar_df), value = TRUE)
cor_rows <- list()
for (pc in pc_cols) {
  for (cv in num_covs) {
    ok <- !is.na(covar_df[[pc]]) & !is.na(covar_df[[cv]])
    if (sum(ok) >= 10) {
      cor_rows[[length(cor_rows) + 1]] <- tibble(
        PC = pc,
        covar = cv,
        n = sum(ok),
        r = cor(covar_df[[pc]][ok], covar_df[[cv]][ok]),
        p = tryCatch(cor.test(covar_df[[pc]][ok], covar_df[[cv]][ok])$p.value, error = function(e) NA_real_)
      )
    }
  }
}
cov_cor <- bind_rows(cor_rows)
readr::write_tsv(cov_cor, file.path(tab_dir, "PC_numeric_covariate_correlations.tsv"))






# Simple ANOVA for batch-like covariates
library(future.apply)

# choose workers (cap at 16)
n_workers <- min(16, future::availableCores())
future::plan(future::multisession, workers = n_workers)

fac_covs <- c("SMCENTER", "SMNABTCH", "SMGEBTCH")
pc_cols  <- grep("^PC[0-9]+$", names(covar_df), value = TRUE)

# precompute cleaned per-covariate frames to avoid repeating filter work
cov_frames <- lapply(fac_covs, function(cv) {
  df <- covar_df %>%
    dplyr::select(dplyr::all_of(c(pc_cols, cv))) %>%
    dplyr::filter(!is.na(.data[[cv]]))
  list(cv = cv, df = df)
})
names(cov_frames) <- fac_covs

safe_anova <- function(pc, cv) {
  # fetch pre-filtered frame for this covariate
  base <- cov_frames[[cv]]$df
  df <- base %>% dplyr::select(PC = dplyr::all_of(pc), FAC = dplyr::all_of(cv))
  if (nrow(df) < 20 || dplyr::n_distinct(df$FAC) < 2) {
    return(tibble::tibble(PC = pc, covar = cv, n = nrow(df), p = NA_real_))
  }
  df$FAC <- droplevels(factor(df$FAC))
  fit <- tryCatch(lm(PC ~ FAC, data = df), error = function(e) NULL)
  if (is.null(fit)) {
    return(tibble::tibble(PC = pc, covar = cv, n = nrow(df), p = NA_real_))
  }
  p <- tryCatch({
    sm <- anova(fit)
    as.numeric(sm$`Pr(>F)`[1])
  }, error = function(e) NA_real_)
  tibble::tibble(PC = pc, covar = cv, n = nrow(df), p = p)
}

# all combinations
combos <- expand.grid(pc = pc_cols, cv = fac_covs, stringsAsFactors = FALSE)

# run in parallel
anova_rows <- future_lapply(seq_len(nrow(combos)), function(i) {
  safe_anova(combos$pc[i], combos$cv[i])
}, future.seed = TRUE)

cov_anova <- dplyr::bind_rows(anova_rows)
readr::write_csv(cov_anova, file.path(tab_dir, "PC_factor_covariate_ANOVA.csv"))

# return to sequential to avoid affecting later code
future::plan(future::sequential)

message("Parallel ANOVA done with ", n_workers, " workers.")














## -----------------------------------------------------------------------------
## 8) Dot plots of selected genes by tissue, sex panels, ordered by AGE bins
## -----------------------------------------------------------------------------

# Selected genes
sel_genes <- c(
  "TH","DBH","PNMT",
  "ADRB1","ADRB2","ADRB3",
  "NPY","NPY1R",
  "MAOA","MAOB",
  "LIPE","PNPLA2","PLIN1",
  "PRKACA","GNAS",
  "PPARGC1A","UCP1","DIO2",
  "TIE1","ANGPT2","RNASE1","PLVAP","CA2","MPZL2"
)

# Gene map from GCT Description
gene_map <- bind_rows(
  tibble(gene_id = x_sc$gene_id,  gene_symbol = x_sc$gene_symbol),
  tibble(gene_id = x_vis$gene_id, gene_symbol = x_vis$gene_symbol)
) %>%
  filter(!is.na(gene_id) & nzchar(gene_id)) %>%
  mutate(gene_id_clean = sub("\\..*$", "", gene_id)) %>%  # strip version
  distinct(gene_id_clean, .keep_all = TRUE) %>%
  dplyr::select(gene_id_clean, gene_symbol)

# Clean AGE bins and derive numeric helpers
# AGE is a string like "50-59" or "80+" in GTEx
clean_age_bin <- function(age_chr) {
  lab <- as.character(age_chr)
  lab <- gsub("\\s*–\\s*|\\s*-\\s*", "-", lab)  # unify dash variants
  lab <- gsub("^\\s+|\\s+$", "", lab)
  lab
}
age_bounds_from_AGE <- function(age_chr) {
  x <- clean_age_bin(age_chr)
  age_min <- age_max <- rep(NA_real_, length(x))
  # "A-B"
  m <- stringr::str_match(x, "^\\s*([0-9]{1,3})\\s*-\\s*([0-9]{1,3})\\s*$")
  has_range <- !is.na(m[,2]) & !is.na(m[,3])
  age_min[has_range] <- as.numeric(m[has_range,2])
  age_max[has_range] <- as.numeric(m[has_range,3])
  # "A+"
  plus_idx <- grepl("^\\s*([0-9]{1,3})\\s*\\+\\s*$", x)
  if (any(plus_idx)) {
    base <- as.numeric(gsub("\\D", "", x[plus_idx]))
    age_min[plus_idx] <- base
    age_max[plus_idx] <- base + 10
  }
  tibble(
    age_bin_label = ifelse(is.na(x) | x == "", NA_character_, x),
    age_min = age_min,
    age_max = age_max,
    age_mid = ifelse(is.na(age_min) | is.na(age_max), NA_real_, (age_min + age_max) / 2)
  )
}

# Attach AGE helpers to meta if not already present
if (!all(c("age_min","age_max","age_mid") %in% names(meta))) {
  ab <- age_bounds_from_AGE(meta$AGE)
  meta$age_bin_label <- ab$age_bin_label
  meta$age_min  <- ab$age_min
  meta$age_max  <- ab$age_max
  meta$age_mid  <- dplyr::coalesce(meta$age_mid, ab$age_mid)
} else {
  meta$age_bin_label <- clean_age_bin(meta$AGE)
}

# Build tidy VST expression for all samples
expr_all <- as_tibble(t(vst_all), rownames = "SAMPID") %>%
  left_join(
    meta %>% dplyr::select(SAMPID, DEPOT, sex, AGE, age_bin_label, age_min, age_max, age_mid),
    by = "SAMPID"
  ) %>%
  pivot_longer(
    cols = -c(SAMPID, DEPOT, sex, AGE, age_bin_label, age_min, age_max, age_mid),
    names_to = "gene_id", values_to = "vst"
  ) %>%
  mutate(gene_id_clean = sub("\\..*$", "", gene_id)) %>%  # strip version in long table
  left_join(gene_map, by = "gene_id_clean") %>%
  mutate(gene = dplyr::coalesce(gene_symbol, gene_id_clean)) %>%
  mutate(age_bin_lo = suppressWarnings(as.numeric(sub("-.*$", "", age_bin_label))))

# Save which selected genes are available
avail_tbl <- expr_all %>%
  filter(gene %in% sel_genes) %>%
  distinct(gene) %>%
  arrange(gene)
readr::write_csv(avail_tbl, file.path(tab_dir, "selected_genes_available.csv"))









































# 
# 
# # 
# # 
# # # Helper: one plot per gene per tissue, stacked by sex
# # plot_gene_tissue_dots <- function(df_all, gene_symbol, tissue, out_png, out_csv = NULL) {
# #   df <- df_all %>%
# #     filter(DEPOT == tissue, gene == gene_symbol) %>%
# #     filter(!is.na(sex)) %>%
# #     mutate(
# #       age_bin_lo = ifelse(is.na(age_bin_lo), Inf, age_bin_lo)
# #     ) %>%
# #     arrange(age_bin_lo, age_mid, vst)
# #   
# #   if (nrow(df) == 0) {
# #     message("No data for gene ", gene_symbol, " in ", tissue)
# #     return(invisible(NULL))
# #   }
# #   
# #   # Order samples: age bin -> age_mid -> expression
# #   df <- df %>%
# #     mutate(
# #       samp_order = factor(SAMPID, levels = unique(SAMPID[order(age_bin_lo, age_mid, vst)])),
# #       age_bin_f  = factor(age_bin_label, levels = sort(unique(age_bin_label[order(age_bin_lo)])))
# #     )
# #   
# #   p <- ggplot(df, aes(x = samp_order, y = vst, color = age_bin_f)) +
# #     geom_point(size = 1.8, alpha = 0.85) +
# #     facet_wrap(~ sex, ncol = 1, scales = "free_x") +
# #     labs(
# #       title = paste0(gene_symbol, " in ", tissue),
# #       subtitle = "Dots are samples. Panels split by sex. Color is AGE bin. X ordered by bin, then age_mid, then expression.",
# #       x = "Samples", y = "VST expression", color = "AGE bin"
# #     ) +
# #     theme_bw(base_size = 10) +
# #     theme(
# #       axis.text.x = element_blank(),
# #       axis.ticks.x = element_blank(),
# #       strip.background = element_rect(fill = "grey92")
# #     )
# #   
# #   ggsave(out_png, p, width = 12, height = 5.5, dpi = 300)
# #   
# #   if (!is.null(out_csv)) {
# #     # Save the panel data behind the plot for reproducibility
# #     df_out <- df %>%
# #       select(SAMPID, DEPOT, sex, AGE, age_bin_label, age_mid, gene, vst) %>%
# #       arrange(sex, age_bin_label, age_mid, vst)
# #     readr::write_csv(df_out, out_csv)
# #   }
# #   
# #   invisible(p)
# # }
# # 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # Helper that avoids large globals
# plot_gene_tissue_dots_from_mats <- function(gsym, tissue, vst_mat, meta_df, fig_dir, tab_dir) {
#   strip_ver <- function(x) sub("\\.\\d+$", "", x)
#   # locate row by symbol or Ensembl id; adjust if you maintain a symbol map elsewhere
#   row_ix <- which(rownames(vst_mat) == gsym)
#   if (!length(row_ix)) {
#     row_ix <- which(strip_ver(rownames(vst_mat)) == strip_ver(gsym))
#     if (!length(row_ix)) return(invisible(NULL))
#   }
#   vec <- as.numeric(vst_mat[row_ix, ])
#   df  <- tibble(SAMPID = colnames(vst_mat), vst = vec) %>%
#     left_join(meta_df %>% select(SAMPID, DEPOT, sex, AGE, age_bin_label, age_mid, age_bin_lo),
#               by = "SAMPID") %>%
#     filter(DEPOT == tissue, !is.na(sex)) %>%
#     mutate(age_bin_lo = ifelse(is.na(age_bin_lo), Inf, age_bin_lo)) %>%
#     arrange(age_bin_lo, age_mid, vst)
#   
#   if (!nrow(df)) return(invisible(NULL))
#   
#   df <- df %>%
#     mutate(
#       samp_order = factor(SAMPID, levels = unique(SAMPID[order(age_bin_lo, age_mid, vst)])),
#       age_bin_f  = factor(age_bin_label, levels = sort(unique(age_bin_label[order(age_bin_lo)])))
#     )
#   
#   p <- ggplot(df, aes(x = samp_order, y = vst, color = age_bin_f)) +
#     geom_point(size = 1.8, alpha = 0.85) +
#     facet_wrap(~ sex, ncol = 1, scales = "free_x") +
#     labs(
#       title = paste0(gsym, " in ", tissue),
#       subtitle = "Dots are samples. Panels split by sex. Color is AGE bin. X ordered by bin, then age_mid, then expression.",
#       x = "Samples", y = "VST expression", color = "AGE bin"
#     ) +
#     theme_bw(base_size = 10) +
#     theme(
#       axis.text.x = element_blank(),
#       axis.ticks.x = element_blank(),
#       strip.background = element_rect(fill = "grey92")
#     )
#   
#   out_png <- file.path(fig_dir, paste0("Dots_", gsym, "_", gsub(" ", "_", tissue), ".png"))
#   out_csv <- file.path(tab_dir, paste0("Dots_", gsym, "_", gsub(" ", "_", tissue), ".csv"))
#   ggsave(out_png, p, width = 12, height = 5.5, dpi = 300)
#   
#   readr::write_csv(
#     df %>% select(SAMPID, DEPOT, sex, AGE, age_bin_label, age_mid, gene = !!gsym, vst) %>% mutate(gene = gsym),
#     out_csv
#   )
#   invisible(NULL)
# }
# 
# # Cache big objects to avoid exporting multi-GB globals to each worker
# cache_path <- file.path(rds_dir, "dotplot_cache.rds")
# saveRDS(list(vst_all = vst_all, meta = meta), cache_path)
# 
# # Remove any giant intermediates that could be captured by futures
# if (exists("expr_all")) { rm(expr_all); gc() }
# 
# # Parallel plan
# n_workers <- min(16, future::availableCores())
# future::plan(future::multisession, workers = n_workers)
# 
# # Split genes so each worker loads the cache once
# chunks <- split(sel_genes, cut(seq_along(sel_genes), breaks = n_workers, labels = FALSE))
# 
# # Run in parallel
# invisible(
#   future_lapply(chunks, function(gene_chunk, cache_path, fig_dir, tab_dir) {
#     env <- readRDS(cache_path)
#     vst_mat <- env$vst_all
#     meta_df <- env$meta
#     for (g in gene_chunk) {
#       try(plot_gene_tissue_dots_from_mats(g, "Subcutaneous", vst_mat, meta_df, fig_dir, tab_dir), silent = TRUE)
#       try(plot_gene_tissue_dots_from_mats(g, "Visceral",      vst_mat, meta_df, fig_dir, tab_dir), silent = TRUE)
#     }
#     NULL
#   }, cache_path = cache_path, fig_dir = fig_dir, tab_dir = tab_dir, future.seed = TRUE)
# )
# 
# # Return to sequential
# future::plan(future::sequential)
# 
# message("Parallel dot plots complete with ", n_workers, " workers. Files written to: ", normalizePath(fig_dir))
# 
# 
# 
# # ## -----------------------------------------------------------------------------
# # ## Parallel dot plots for selected genes (Section 8 - parallel)
# # ## -----------------------------------------------------------------------------
# # 
# # if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply")
# # library(future.apply)
# # 
# # # Pick worker count
# # n_workers <- min(16, future::availableCores())
# # future::plan(future::multisession, workers = n_workers)
# # 
# # plot_one_gene_both_tissues <- function(g) {
# #   # Subcutaneous
# #   try(
# #     plot_gene_tissue_dots(
# #       df_all   = expr_all,
# #       gene_symbol = g,
# #       tissue   = "Subcutaneous",
# #       out_png  = file.path(fig_dir, paste0("Dots_", g, "_Subcutaneous.png")),
# #       out_csv  = file.path(tab_dir, paste0("Dots_", g, "_Subcutaneous.csv"))
# #     ),
# #     silent = TRUE
# #   )
# #   # Visceral
# #   try(
# #     plot_gene_tissue_dots(
# #       df_all   = expr_all,
# #       gene_symbol = g,
# #       tissue   = "Visceral",
# #       out_png  = file.path(fig_dir, paste0("Dots_", g, "_Visceral.png")),
# #       out_csv  = file.path(tab_dir, paste0("Dots_", g, "_Visceral.csv"))
# #     ),
# #     silent = TRUE
# #   )
# #   invisible(NULL)
# # }
# # 
# # # Run in parallel across genes
# # invisible(
# #   future_lapply(sel_genes, plot_one_gene_both_tissues, future.seed = TRUE)
# # )
# # 
# # # Return to sequential so later code is predictable
# # future::plan(future::sequential)
# # 
# # message("Section 8 complete: gene dot plots written in parallel with ", n_workers, " workers. Tables are CSV.")
# # 
# # 
# # 
# # 
# # 
# # 
# 
# 
# 
# 
# 
# 
# 
# 
# 

# 
# map_sym2id <- setNames(gene_map$gene_id_clean, gene_map$gene_symbol)
# 
# plot_gene_tissue_dots_from_mats <- function(gsym,tissue,vst_mat,meta_df,map_sym2id,fig_dir,tab_dir){
#   strip_ver <- function(x) sub("\\.\\d+$","",x)
#   target <- if(gsym %in% names(map_sym2id)) map_sym2id[[gsym]] else gsym
#   rn <- rownames(vst_mat); row_ix <- which(strip_ver(rn)==strip_ver(target))
#   if(!length(row_ix)){message("No row for ",gsym," (",target,")"); return(invisible(NULL))}
#   df <- tibble(SAMPID=colnames(vst_mat),vst=as.numeric(vst_mat[row_ix,])) %>%
#     left_join(meta_df %>% select(SAMPID,DEPOT,sex,AGE,age_bin_label,age_mid,age_bin_lo),by="SAMPID") %>%
#     filter(DEPOT==tissue,!is.na(sex),!is.na(vst)) %>%
#     mutate(age_bin_lo=ifelse(is.na(age_bin_lo),Inf,age_bin_lo)) %>%
#     arrange(age_bin_lo,age_mid,vst)
#   if(!nrow(df)) return(invisible(NULL))
#   df <- df %>% mutate(samp_order=factor(SAMPID,levels=unique(SAMPID[order(age_bin_lo,age_mid,vst)])),
#                       age_bin_f=factor(age_bin_label,levels=unique(age_bin_label[order(age_bin_lo)])),
#                       gene=gsym)
#   p <- ggplot(df,aes(x=samp_order,y=vst,color=age_bin_f))+
#     geom_point(size=1.8,alpha=0.85)+facet_wrap(~sex,ncol=1,scales="free_x")+
#     labs(title=paste0(gsym," in ",tissue),subtitle="Dots are samples; color=age bin.",x="Samples",y="VST expression",color="Age bin")+
#     theme_bw(base_size=10)+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),strip.background=element_rect(fill="grey92"))
#   if(!dir.exists(fig_dir)) dir.create(fig_dir,recursive=TRUE); if(!dir.exists(tab_dir)) dir.create(tab_dir,recursive=TRUE)
#   ggsave(file.path(fig_dir,paste0("Dots_",gsym,"_",gsub(" ","_",tissue),".png")),p,width=8,height=4,dpi=300)
#   readr::write_csv(df %>% select(SAMPID,DEPOT,sex,AGE,age_bin_label,age_mid,age_bin_lo,gene,vst),
#                    file.path(tab_dir,paste0("Dots_",gsym,"_",gsub(" ","_",tissue),".csv")))
#   invisible(p)
# }
# 
# run_dotplots <- function(sel_genes,vst_mat,meta_df,map_sym2id,fig_dir,tab_dir,tissues=c("Subcutaneous","Visceral")){
#   if(!dir.exists(fig_dir)) dir.create(fig_dir,recursive=TRUE); if(!dir.exists(tab_dir)) dir.create(tab_dir,recursive=TRUE)
#   for(g in sel_genes) for(t in tissues) try(plot_gene_tissue_dots_from_mats(g,t,vst_mat,meta_df,map_sym2id,fig_dir,tab_dir),silent=TRUE)
#   message("Dot plots done for ",length(sel_genes)," genes in ",paste(tissues,collapse=", "))
# }
# 
# # Run
# run_dotplots(sel_genes, vst_all, meta, map_sym2id, "figs", "tables")
# 



## -----------------------------------------------------------------------------
## 8) Dot plots of selected genes (Subcutaneous & Visceral)
## -----------------------------------------------------------------------------

map_sym2id <- setNames(gene_map$gene_id_clean, gene_map$gene_symbol)
rn_clean <- sub("\\..*$", "", rownames(vst_all))
sel_genes_avail <- sel_genes[sel_genes %in% names(map_sym2id) &
                               map_sym2id[sel_genes] %in% rn_clean]
message("Selected genes available in vst_all: ",
        if (length(sel_genes_avail)) paste(sel_genes_avail, collapse = ", ") else "none")

plot_gene_tissue_dots_from_mats <- function(gsym, tissue, vst_all, meta,
                                            fig_dir, tab_dir, map_sym2id) {
  strip_ver <- function(x) sub("\\.\\d+$", "", x)
  target <- if (gsym %in% names(map_sym2id)) map_sym2id[[gsym]] else gsym
  row_ix <- which(strip_ver(rownames(vst_all)) == strip_ver(target))
  if (!length(row_ix)) { message("No row for ", gsym, " (", target, ")"); return(invisible(NULL)) }
  
  df <- tibble(SAMPID = colnames(vst_all), vst = as.numeric(vst_all[row_ix,])) %>%
    left_join(meta %>% dplyr::select(SAMPID, DEPOT, sex, AGE, age_bin_label, age_mid), by = "SAMPID") %>%
    filter(DEPOT == tissue, !is.na(sex), !is.na(vst))
  if (!nrow(df)) { message("No data for ", gsym, " in ", tissue); return(invisible(NULL)) }
  
  df <- df %>%
    mutate(age_bin_lo = suppressWarnings(as.numeric(sub("-.*$", "", age_bin_label))),
           age_bin_lo = ifelse(is.na(age_bin_lo), Inf, age_bin_lo)) %>%
    arrange(age_bin_lo, age_mid, vst) %>%
    mutate(samp_order = factor(SAMPID, levels = unique(SAMPID)),
           age_bin_f = factor(age_bin_label, levels = unique(age_bin_label[order(age_bin_lo)])))
  
  #print(head(df))
  
  p <- ggplot(df, aes(x = samp_order, y = vst, color = age_bin_f)) +
    geom_point(size = 1.5, alpha = 0.85) +
    facet_wrap(~sex, ncol = 1, scales = "free_x") +
    labs(title = paste0(gsym, " in ", tissue),
         subtitle = "Samples ordered by age bin, age_mid, expression",
         x = "Samples", y = "VST expression", color = "Age bin") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "grey92"))
  
  out_png <- file.path(fig_dir, paste0("Dots_", gsym, "_", gsub(" ", "_", tissue), ".png"))
  out_csv <- file.path(tab_dir, paste0("Dots_", gsym, "_", gsub(" ", "_", tissue), ".csv"))
  
  ggsave(out_png, p, width = 8, height = 4, dpi = 300)
  readr::write_csv(df %>% dplyr::select(SAMPID, DEPOT, sex, AGE, age_bin_label, age_mid, vst), out_csv)
  invisible(p)
}

for (g in sel_genes_avail) {
  for (t in c("Subcutaneous", "Visceral")) {
    try(
      plot_gene_tissue_dots_from_mats(
        g, t, vst_all, meta, fig_dir, tab_dir, map_sym2id
      ),
      silent = FALSE
    )
    exit_status <- getOption("try.error")
  }
}


message("Dot plots written to: ", normalizePath(fig_dir))








plot_gene_tissue_violin_from_mats <- function(gsym, tissue, vst_all, meta,
                                              fig_dir, tab_dir, map_sym2id) {
  strip_ver <- function(x) sub("\\.\\d+$", "", x)
  target <- if (gsym %in% names(map_sym2id)) map_sym2id[[gsym]] else gsym
  row_ix <- which(strip_ver(rownames(vst_all)) == strip_ver(target))
  if (!length(row_ix)) {
    message("No row for ", gsym, " (", target, ") for violin")
    return(invisible(NULL))
  }
  
  df <- tibble(
    SAMPID = colnames(vst_all),
    vst    = as.numeric(vst_all[row_ix, ])
  ) %>%
    dplyr::left_join(
      meta %>%
        dplyr::select(SAMPID, DEPOT, sex, AGE, age_bin_label, age_mid),
      by = "SAMPID"
    ) %>%
    dplyr::filter(DEPOT == tissue, !is.na(sex), !is.na(vst))
  if (!nrow(df)) {
    message("No data for ", gsym, " in ", tissue, " for violin")
    return(invisible(NULL))
  }
  
  df <- df %>%
    dplyr::mutate(
      age_bin_label = ifelse(is.na(age_bin_label), "NA", age_bin_label),
      age_bin_f     = factor(age_bin_label,
                             levels = sort(unique(age_bin_label)))
    )
  
  ## stats: only 40-49, 50-59, 60-69
  cmp_bins <- c("40-49", "50-59", "60-69")
  df_stats <- df %>% dplyr::filter(age_bin_label %in% cmp_bins)
  
  stat_rows <- list()
  if (nrow(df_stats) > 0) {
    pairs <- list(
      c("40-49", "50-59"),
      c("40-49", "60-69"),
      c("50-59", "60-69")
    )
    for (sx in sort(unique(df_stats$sex))) {
      df_s <- df_stats %>% dplyr::filter(sex == sx)
      if (nrow(df_s) < 3) next
      for (pr in pairs) {
        g1 <- pr[1]; g2 <- pr[2]
        sub <- df_s %>% dplyr::filter(age_bin_label %in% c(g1, g2))
        n1 <- sum(sub$age_bin_label == g1)
        n2 <- sum(sub$age_bin_label == g2)
        if (n1 >= 3 && n2 >= 3) {
          wt <- suppressWarnings(
            wilcox.test(vst ~ age_bin_label, data = sub, exact = FALSE)
          )
          stat_rows[[length(stat_rows) + 1]] <- tibble::tibble(
            gene        = gsym,
            tissue      = tissue,
            sex         = sx,
            age_bin_1   = g1,
            age_bin_2   = g2,
            n_1         = n1,
            n_2         = n2,
            statistic   = unname(wt$statistic),
            p_value     = wt$p.value
          )
        }
      }
    }
  }
  
  stats_tbl <- if (length(stat_rows)) {
    out <- dplyr::bind_rows(stat_rows)
    out <- out %>%
      dplyr::group_by(gene, tissue, sex) %>%
      dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
      dplyr::ungroup()
    out
  } else {
    tibble::tibble()
  }
  
  ## violin plot
  p <- ggplot(df, aes(x = age_bin_f, y = vst, fill = age_bin_f)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.5) +
    facet_wrap(~sex, ncol = 1) +
    labs(
      title = paste0(gsym, " in ", tissue, " (violin)"),
      x     = "Age bin",
      y     = "VST expression",
      fill  = "Age bin"
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey92"),
      legend.position  = "none"
    )
  
  # optional on-plot p labels if ggpubr is available
  if (requireNamespace("ggpubr", quietly = TRUE) && nrow(df_stats) > 0) {
    cmp_list <- list(
      c("40-49", "50-59"),
      c("40-49", "60-69"),
      c("50-59", "60-69")
    )
    p <- p + ggpubr::stat_compare_means(
      data        = df_stats,
      comparisons = cmp_list,
      method      = "wilcox.test",
      label       = "p.signif",
      hide.ns     = TRUE
    )
  }
  
  out_png <- file.path(
    fig_dir,
    paste0("Violin_", gsym, "_", gsub(" ", "_", tissue), ".png")
  )
  out_csv <- file.path(
    tab_dir,
    paste0("Violin_", gsym, "_", gsub(" ", "_", tissue), ".csv")
  )
  out_stat <- file.path(
    tab_dir,
    paste0("Violin_stats_", gsym, "_", gsub(" ", "_", tissue), ".csv")
  )
  
  ggsave(out_png, p, width = 6, height = 4.5, dpi = 300)
  
  readr::write_csv(
    df %>% dplyr::select(SAMPID, DEPOT, sex, AGE, age_bin_label, age_mid, vst),
    out_csv
  )
  if (nrow(stats_tbl)) readr::write_csv(stats_tbl, out_stat)
  
  invisible(p)
}




plot_gene_tissue_violin_from_mats <- function(gsym, tissue, vst_all, meta,
                                              fig_dir, tab_dir, map_sym2id) {
  strip_ver <- function(x) sub("\\.\\d+$", "", x)
  target <- if (gsym %in% names(map_sym2id)) map_sym2id[[gsym]] else gsym
  row_ix <- which(strip_ver(rownames(vst_all)) == strip_ver(target))
  if (!length(row_ix)) {
    message("No row for ", gsym, " (", target, ") for violin")
    return(invisible(NULL))
  }
  
  df <- tibble(SAMPID = colnames(vst_all),
               vst    = as.numeric(vst_all[row_ix, ])) %>%
    dplyr::left_join(
      meta %>% dplyr::select(SAMPID, DEPOT, sex, AGE, age_bin_label, age_mid),
      by = "SAMPID"
    ) %>%
    dplyr::filter(DEPOT == tissue, !is.na(sex), !is.na(vst))
  if (!nrow(df)) {
    message("No data for ", gsym, " in ", tissue, " for violin")
    return(invisible(NULL))
  }
  
  df <- df %>%
    dplyr::mutate(
      age_bin_label = ifelse(is.na(age_bin_label), "NA", age_bin_label),
      age_bin_f     = factor(age_bin_label,
                             levels = sort(unique(age_bin_label)))
    )
  
  ## stats only for 40-49, 50-59, 60-69
  cmp_bins <- c("40-49", "50-59", "60-69")
  df_stats <- df %>% dplyr::filter(age_bin_label %in% cmp_bins)
  
  stat_rows <- list()
  if (nrow(df_stats) > 0) {
    pairs <- list(
      c("40-49", "50-59"),
      c("40-49", "60-69"),
      c("50-59", "60-69")
    )
    for (sx in sort(unique(df_stats$sex))) {
      df_s <- df_stats %>% dplyr::filter(sex == sx)
      if (nrow(df_s) < 3) next
      for (pr in pairs) {
        g1 <- pr[1]; g2 <- pr[2]
        sub <- df_s %>% dplyr::filter(age_bin_label %in% c(g1, g2))
        n1 <- sum(sub$age_bin_label == g1)
        n2 <- sum(sub$age_bin_label == g2)
        if (n1 >= 3 && n2 >= 3) {
          wt <- suppressWarnings(
            wilcox.test(vst ~ age_bin_label, data = sub, exact = FALSE)
          )
          stat_rows[[length(stat_rows) + 1]] <- tibble::tibble(
            gene      = gsym,
            tissue    = tissue,
            sex       = sx,
            age_bin_1 = g1,
            age_bin_2 = g2,
            n_1       = n1,
            n_2       = n2,
            statistic = unname(wt$statistic),
            p_value   = wt$p.value
          )
        }
      }
    }
  }
  
  stats_tbl <- if (length(stat_rows)) {
    dplyr::bind_rows(stat_rows) %>%
      dplyr::group_by(gene, tissue, sex) %>%
      dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
      dplyr::ungroup()
  } else {
    tibble::tibble(
      gene      = character(),
      tissue    = character(),
      sex       = character(),
      age_bin_1 = character(),
      age_bin_2 = character(),
      n_1       = integer(),
      n_2       = integer(),
      statistic = numeric(),
      p_value   = numeric(),
      p_adj     = numeric()
    )
  }
  
  p <- ggplot(df, aes(x = age_bin_f, y = vst, fill = age_bin_f)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.15, outlier.size = 0.5) +
    facet_wrap(~sex, ncol = 1) +
    labs(
      title = paste0(gsym, " in ", tissue, " (violin)"),
      x     = "Age bin",
      y     = "VST expression",
      fill  = "Age bin"
    ) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey92"),
      legend.position  = "none"
    )
  
  # optional labels on the plot if ggpubr is installed
  if (requireNamespace("ggpubr", quietly = TRUE) && nrow(df_stats) > 0) {
    cmp_list <- list(
      c("40-49", "50-59"),
      c("40-49", "60-69"),
      c("50-59", "60-69")
    )
    p <- p + ggpubr::stat_compare_means(
      data        = df_stats,
      comparisons = cmp_list,
      method      = "wilcox.test",
      label       = "p.signif",
      hide.ns     = TRUE
    )
  }
  
  out_png  <- file.path(fig_dir,  paste0("Violin_", gsym, "_", gsub(" ", "_", tissue), ".png"))
  out_csv  <- file.path(tab_dir,  paste0("Violin_", gsym, "_", gsub(" ", "_", tissue), ".csv"))
  out_stat <- file.path(tab_dir,  paste0("Violin_stats_", gsym, "_", gsub(" ", "_", tissue), ".csv"))
  
  ggsave(out_png, p, width = 6, height = 4.5, dpi = 300)
  readr::write_csv(
    df %>% dplyr::select(SAMPID, DEPOT, sex, AGE, age_bin_label, age_mid, vst),
    out_csv
  )
  readr::write_csv(stats_tbl, out_stat)
  
  invisible(p)
}


for (g in sel_genes_avail) {
  for (t in c("Subcutaneous", "Visceral")) {
    try(
      plot_gene_tissue_violin_from_mats(
        g, t, vst_all, meta, fig_dir, tab_dir, map_sym2id
      ),
      silent = FALSE
    )
  }
}





# adjust if your stats filenames differ
stat_files <- list.files(
  tab_dir,
  pattern = "^Violin_stats_.*\\.csv$",
  full.names = TRUE
)

if (!length(stat_files)) {
  message("No Violin_stats_*.csv files found in: ", tab_dir)
} else {
  violin_stats_all <- map_dfr(
    stat_files,
    ~ read_csv(.x, show_col_types = FALSE)
  )
  
  # keep only your target genes (sel_genes or sel_genes_avail)
  target_genes <- sel_genes_avail
  violin_stats_all <- violin_stats_all %>%
    filter(gene %in% target_genes)
  
  # significant only (BH adjusted p < 0.05)
  violin_stats_sig <- violin_stats_all %>%
    filter(!is.na(p_adj), p_adj < 0.05) %>%
    arrange(gene, tissue, sex, age_bin_1, age_bin_2)
  
  out_all <- file.path(tab_dir, "Violin_stats_MASTER_all_targets.csv")
  out_sig <- file.path(tab_dir, "Violin_stats_MASTER_significant_targets.csv")
  
  write_csv(violin_stats_all, out_all)
  write_csv(violin_stats_sig, out_sig)
  
  message("Master stats written: ",
          "\n  ", out_all,
          "\n  ", out_sig)
}







































# 
# 
# ## -----------------------------------------------------------------------------
# ## 8) Define strict female groups for DE using binned ages
# ## -----------------------------------------------------------------------------
# 
# # Use age_min and age_max to avoid overlap with the 45 and 55 cut points
# meta <- meta %>%
#   mutate(
#     female_group_strict = case_when(
#       sex == "Female" & !is.na(age_max) & age_max < 45 ~ "Female_pre",
#       sex == "Female" & !is.na(age_min) & age_min > 55 ~ "Female_post",
#       TRUE ~ NA_character_
#     ),
#     DEPOT = factor(DEPOT)
#   )
# 
# # Summaries
# group_counts <- meta %>%
#   count(DEPOT, female_group_strict, name = "n") %>%
#   arrange(DEPOT, female_group_strict)
# readr::write_tsv(group_counts, file.path(tab_dir, "female_strict_group_counts_by_DEPOT.tsv"))
# 
# 
# 
# 
# ## -----------------------------------------------------------------------------
# ## 9) DESeq2: Female post vs pre, one model adjusted for DEPOT
# ## -----------------------------------------------------------------------------
# 
# meta_f <- meta %>%
#   filter(!is.na(female_group_strict)) %>%
#   mutate(group_female = factor(female_group_strict, levels = c("Female_pre", "Female_post")))
# 
# counts_f <- counts_all[, meta_f$SAMPID, drop = FALSE]
# 
# dds_f <- DESeqDataSetFromMatrix(
#   countData = round(counts_f),
#   colData   = as.data.frame(meta_f),
#   design    = ~ DEPOT + group_female
# )
# keep_f <- rowSums(counts(dds_f) >= 10) >= 10
# dds_f   <- dds_f[keep_f, ]
# dds_f   <- DESeq(dds_f)
# 
# res_f <- results(dds_f, contrast = c("group_female", "Female_post", "Female_pre"))
# # Shrink LFC for stability
# res_f_shr <- lfcShrink(dds_f, coef = "group_female_Female_post_vs_Female_pre", type = "apeglm")
# 
# res_tbl <- as.data.frame(res_f_shr) %>%
#   rownames_to_column("gene_id") %>%
#   arrange(padj)
# readr::write_tsv(res_tbl, file.path(tab_dir, "DE_Female_post_vs_pre_adjust_DEPOT_strict.tsv"))
# 
# # Volcano
# volc_df <- res_tbl %>% mutate(sig = ifelse(!is.na(padj) & padj < 0.05, "FDR<0.05", "NS"))
# p_volc <- ggplot(volc_df, aes(log2FoldChange, -log10(pvalue), color = sig)) +
#   geom_point(alpha = 0.7, size = 1) +
#   theme_bw() +
#   labs(title = "Female post vs pre (strict, adjusted for DEPOT)", x = "log2FC", y = "-log10 p")
# ggsave(file.path(fig_dir, "Volcano_Female_post_vs_pre_adjust_DEPOT_strict.png"), p_volc,
#        width = 7, height = 5, dpi = 300)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## -----------------------------------------------------------------------------
# ## Age-bin parsing: add age_min, age_max, and age_mid to `meta`
# ## -----------------------------------------------------------------------------
# parse_age_bin <- function(x) {
#   # Handles "50-59", "60 - 69", "80+", "30–39" etc.
#   x <- as.character(x)
#   age_min <- rep(NA_real_, length(x))
#   age_max <- rep(NA_real_, length(x))
#   
#   # "A-B" formats
#   m <- stringr::str_match(x, "\\s*([0-9]{1,3})\\D+([0-9]{1,3})\\s*")
#   has_range <- !is.na(m[, 2]) & !is.na(m[, 3])
#   age_min[has_range] <- as.numeric(m[has_range, 2])
#   age_max[has_range] <- as.numeric(m[has_range, 3])
#   
#   # "A+" formats
#   plus_idx <- grepl("^[^0-9]*([0-9]{1,3})\\s*\\+$", x)
#   if (any(plus_idx)) {
#     base <- as.numeric(gsub("\\D", "", x[plus_idx]))
#     age_min[plus_idx] <- base
#     age_max[plus_idx] <- base + 10  # treat as decade width for midpoint
#   }
#   
#   tibble(age_min = age_min,
#          age_max = age_max,
#          age_mid = ifelse(is.na(age_min) | is.na(age_max), NA_real_, (age_min + age_max) / 2))
# }
# 
# # If meta already has age_years from earlier steps, keep it.
# # Otherwise derive from AGE or AGE_BIN using parse_age_bin.
# if (!("age_years" %in% names(meta)) || all(is.na(meta$age_years))) {
#   src <- if ("AGE" %in% names(meta)) meta$AGE else if ("AGE_BIN" %in% names(meta)) meta$AGE_BIN else NA
#   ab  <- parse_age_bin(src)
#   meta$age_min  <- ab$age_min
#   meta$age_max  <- ab$age_max
#   meta$age_mid  <- ab$age_mid
#   meta$age_years <- meta$age_mid
# } else {
#   # Also compute midpoints alongside existing numeric ages for plotting
#   src <- if ("AGE" %in% names(meta)) meta$AGE else if ("AGE_BIN" %in% names(meta)) meta$AGE_BIN else NA
#   ab  <- parse_age_bin(src)
#   if (!("age_mid" %in% names(meta))) {
#     meta$age_min <- ab$age_min
#     meta$age_max <- ab$age_max
#     meta$age_mid <- ab$age_mid
#   }
# }
# 
# # Rebuild the PCA scores table you saved earlier, but now with age_mid attached.
# # If you saved PCA scores already, you can skip re-writing; otherwise it is harmless.
# 
# ## -----------------------------------------------------------------------------
# ## Extra PCA plots: color by AGE (age_mid)
# ## -----------------------------------------------------------------------------
# pca_from_vst_color_age <- function(vst_mat, meta_df, title_suffix, out_png, out_tsv = NULL) {
#   top_var <- head(order(matrixStats::rowVars(vst_mat), decreasing = TRUE), 5000)
#   pca     <- prcomp(t(vst_mat[top_var, , drop = FALSE]), center = TRUE, scale. = FALSE)
#   varpct  <- (pca$sdev^2) / sum(pca$sdev^2) * 100
#   
#   pca_df <- as_tibble(pca$x, rownames = "SAMPID") %>%
#     left_join(meta_df %>% select(SAMPID, DEPOT, sex, age_mid), by = "SAMPID")
#   
#   p <- ggplot(pca_df, aes(PC1, PC2, color = age_mid, shape = DEPOT)) +
#     geom_point(size = 2, alpha = 0.9, na.rm = TRUE) +
#     labs(
#       title = paste("GTEx adipose PCA", title_suffix, "| color = AGE midpoint"),
#       x = paste0("PC1 (", round(varpct[1], 1), "%)"),
#       y = paste0("PC2 (", round(varpct[2], 1), "%)"),
#       color = "Age (midpoint)"
#     ) +
#     theme_bw()
#   ggsave(out_png, p, width = 7, height = 5, dpi = 300)
#   
#   if (!is.null(out_tsv)) readr::write_tsv(pca_df, out_tsv)
#   invisible(list(pca = pca, pca_df = pca_df, var_explained = varpct))
# }
# 
# # 1) All samples, both depots, color by age
# pca_both_age <- pca_from_vst_color_age(
#   vst_mat  = vst_all,
#   meta_df  = meta,
#   title_suffix = "(all samples, both depots)",
#   out_png  = file.path(fig_dir, "PCA_PC1_PC2_ALL_both_DEPOT_COLOR_BY_AGE.png"),
#   out_tsv  = file.path(tab_dir, "PCA_scores_ALL_both_DEPOT_COLOR_BY_AGE.tsv")
# )
# 
# # 2) Depot-specific, color by age
# for (dep in c("Subcutaneous", "Visceral")) {
#   idx <- meta$DEPOT == dep
#   if (sum(idx) >= 10) {
#     pca_from_vst_color_age(
#       vst_mat  = vst_all[, meta$SAMPID[idx], drop = FALSE],
#       meta_df  = meta[idx, , drop = FALSE],
#       title_suffix = paste0("(all samples, ", dep, " only)"),
#       out_png  = file.path(fig_dir, paste0("PCA_PC1_PC2_ALL_", gsub(" ", "_", dep), "_COLOR_BY_AGE.png")),
#       out_tsv  = file.path(tab_dir, paste0("PCA_scores_ALL_", gsub(" ", "_", dep), "_COLOR_BY_AGE.tsv"))
#     )
#   } else {
#     message("Skipping age-colored PCA for ", dep, " due to small sample size.")
#   }
# }
# # 
# # ## -----------------------------------------------------------------------------
# # ## Optional: derive strict groups for DE with binned ages
# # ## -----------------------------------------------------------------------------
# # # Safe grouping from bins:
# # #   Female_pre  if age_max < 45
# # #   Female_post if age_min > 55
# # #   otherwise NA (ambiguous)
# # meta <- meta %>%
# #   mutate(
# #     female_band_strict = dplyr::case_when(
# #       sex == "Female" & !is.na(age_max) & age_max < 45 ~ "Female_pre",
# #       sex == "Female" & !is.na(age_min) & age_min > 55 ~ "Female_post",
# #       TRUE ~ NA_character_
# #     )
# #   )
# # readr::write_tsv(meta, file.path(tab_dir, "sample_metadata_with_age_mid.tsv"))
# # 
# # 
# 
# 
# 
# ## -----------------------------------------------------------------------------
# ## Correlate PCs with age and sex
# ##   - Pearson r of PCk with age_mid
# ##   - Point-biserial r of PCk with sex (Male=1, Female=0)
# ##   - Linear model PCk ~ age_mid + sex + DEPOT  -> partial associations
# ##   Saves one TSV per analysis and a combined CSV for convenience
# ## -----------------------------------------------------------------------------
# 
# compute_pc_associations <- function(pca_df, out_prefix) {
#   stopifnot(all(c("SAMPID","DEPOT","sex") %in% names(pca_df)))
#   # If age_mid is missing, try age_years as fallback
#   if (!"age_mid" %in% names(pca_df)) {
#     if ("age_years" %in% names(pca_df)) pca_df$age_mid <- pca_df$age_years else pca_df$age_mid <- NA_real_
#   }
#   
#   # Numeric sex coding for simple correlation
#   pca_df <- pca_df %>%
#     mutate(sex_num = case_when(
#       sex %in% c("Male","M","1") ~ 1,
#       sex %in% c("Female","F","2") ~ 0,
#       TRUE ~ NA_real_
#     ),
#     DEPOT = factor(DEPOT))
#   
#   pc_cols <- grep("^PC[0-9]+$", names(pca_df), value = TRUE)
#   if (length(pc_cols) == 0) stop("No PC columns found in pca_df")
#   
#   rows <- lapply(pc_cols, function(pc) {
#     vec <- pca_df[[pc]]
#     # Age correlation
#     age_ok <- !is.na(vec) & !is.na(pca_df$age_mid)
#     r_age <- if (sum(age_ok) >= 3) cor(vec[age_ok], pca_df$age_mid[age_ok], method = "pearson") else NA_real_
#     p_age <- tryCatch(if (sum(age_ok) >= 3) cor.test(vec[age_ok], pca_df$age_mid[age_ok])$p.value else NA_real_, error = function(e) NA_real_)
#     
#     # Sex correlation (point-biserial via Pearson with 0-1 coding)
#     sex_ok <- !is.na(vec) & !is.na(pca_df$sex_num)
#     r_sex <- if (sum(sex_ok) >= 3) cor(vec[sex_ok], pca_df$sex_num[sex_ok], method = "pearson") else NA_real_
#     p_sex <- tryCatch(if (sum(sex_ok) >= 3) cor.test(vec[sex_ok], pca_df$sex_num[sex_ok])$p.value else NA_real_, error = function(e) NA_real_)
#     
#     # Linear model with DEPOT and age_mid adjustment
#     fit_df <- pca_df %>% select(!!pc, age_mid, sex_num, DEPOT) %>% na.omit()
#     lm_age_beta <- lm_sex_beta <- lm_age_p <- lm_sex_p <- lm_n <- NA_real_
#     if (nrow(fit_df) >= 10 && length(unique(fit_df$DEPOT)) >= 1) {
#       fit <- lm(fit_df[[pc]] ~ age_mid + sex_num + DEPOT, data = fit_df)
#       sm <- summary(fit)$coefficients
#       lm_age_beta <- if ("age_mid" %in% rownames(sm)) sm["age_mid","Estimate"] else NA_real_
#       lm_age_p    <- if ("age_mid" %in% rownames(sm)) sm["age_mid","Pr(>|t|)"] else NA_real_
#       lm_sex_beta <- if ("sex_num" %in% rownames(sm)) sm["sex_num","Estimate"] else NA_real_
#       lm_sex_p    <- if ("sex_num" %in% rownames(sm)) sm["sex_num","Pr(>|t|)"] else NA_real_
#       lm_n <- nrow(fit_df)
#     }
#     
#     tibble(
#       PC = pc,
#       n_age = sum(age_ok),
#       r_age = r_age,
#       p_age = p_age,
#       n_sex = sum(sex_ok),
#       r_sex = r_sex,
#       p_sex = p_sex,
#       lm_n = lm_n,
#       lm_beta_age_mid = lm_age_beta,
#       lm_p_age_mid = lm_age_p,
#       lm_beta_sex_num = lm_sex_beta,
#       lm_p_sex_num = lm_sex_p
#     )
#   })
#   
#   out <- bind_rows(rows)
#   readr::write_tsv(out, file.path(tab_dir, paste0(out_prefix, "_PC_age_sex_associations.tsv")))
#   invisible(out)
# }
# 
# ## Run associations for:
# ##  - All samples, both depots (from pca_both)
# ##  - Subcutaneous only
# ##  - Visceral only
# 
# assoc_all <- compute_pc_associations(pca_both$pca_df, "ALL_both_DEPOT")
# 
# # If you created depot-specific PCA scores earlier and saved them, you can recompute quickly.
# # Otherwise, subset pca_both$pca_df by DEPOT to reuse the same rotation space.
# 
# assoc_sc <- compute_pc_associations(
#   pca_both$pca_df %>% filter(DEPOT == "Subcutaneous"),
#   "ALL_Subcutaneous"
# )
# 
# assoc_vis <- compute_pc_associations(
#   pca_both$pca_df %>% filter(DEPOT == "Visceral"),
#   "ALL_Visceral"
# )
# 
# # Combine for convenience
# assoc_combined <- bind_rows(
#   assoc_all %>% mutate(scope = "both_depots"),
#   assoc_sc  %>% mutate(scope = "subcutaneous"),
#   assoc_vis %>% mutate(scope = "visceral")
# )
# readr::write_csv(assoc_combined, file.path(tab_dir, "PC_age_sex_associations_combined.csv"))
# 
# ## Optional quick visuals: PC1 and PC2 vs age, colored by sex
# quick_pc_scatter <- function(pca_df, tag) {
#   for (pc in c("PC1","PC2")) {
#     if (!pc %in% names(pca_df)) next
#     df <- pca_df %>% filter(!is.na(age_mid))
#     if (nrow(df) < 10) next
#     p <- ggplot(df, aes(age_mid, .data[[pc]], color = sex)) +
#       geom_point(alpha = 0.6) +
#       geom_smooth(method = "lm", se = FALSE) +
#       theme_bw() +
#       labs(title = paste0(pc, " vs age_mid - ", tag), x = "Age midpoint", y = pc)
#     ggsave(file.path(fig_dir, paste0("Scatter_", pc, "_vs_age_mid_", gsub(" ", "_", tag), ".png")),
#            p, width = 6.5, height = 4.5, dpi = 300)
#   }
# }
# 
# quick_pc_scatter(pca_both$pca_df, "both_depots")
# quick_pc_scatter(pca_both$pca_df %>% filter(DEPOT == "Subcutaneous"), "Subcutaneous")
# quick_pc_scatter(pca_both$pca_df %>% filter(DEPOT == "Visceral"), "Visceral")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ## -----------------------------------------------------------------------------
# ## PCA on ALL samples (no exclusion by age band)
# ## -----------------------------------------------------------------------------
# dds_all <- DESeqDataSetFromMatrix(
#   countData = round(counts_all),
#   colData   = as.data.frame(meta),
#   design    = ~ 1
# )
# keep_all <- rowSums(counts(dds_all) >= 10) >= 20
# dds_all  <- dds_all[keep_all, ]
# dds_all  <- estimateSizeFactors(dds_all)
# vst_all  <- assay(vst(dds_all, blind = TRUE))
# 
# # PCA helper
# pca_from_vst <- function(vst_mat, meta_df, title_suffix, out_png, out_tsv = NULL) {
#   top_var <- head(order(matrixStats::rowVars(vst_mat), decreasing = TRUE), 5000)
#   pca     <- prcomp(t(vst_mat[top_var, , drop = FALSE]), center = TRUE, scale. = FALSE)
#   varpct  <- (pca$sdev^2) / sum(pca$sdev^2) * 100
#   
#   pca_df <- as_tibble(pca$x, rownames = "SAMPID") %>%
#     left_join(meta_df %>% select(SAMPID, DEPOT, sex, age_years, group4), by = "SAMPID")
#   
#   p <- ggplot(pca_df, aes(PC1, PC2, color = sex, shape = DEPOT)) +
#     geom_point(size = 2, alpha = 0.85) +
#     labs(
#       title = paste("GTEx adipose PCA", title_suffix),
#       x = paste0("PC1 (", round(varpct[1], 1), "%)"),
#       y = paste0("PC2 (", round(varpct[2], 1), "%)")
#     ) +
#     theme_bw()
#   ggsave(out_png, p, width = 7, height = 5, dpi = 300)
#   
#   if (!is.null(out_tsv)) readr::write_tsv(pca_df, out_tsv)
#   invisible(list(pca = pca, pca_df = pca_df, var_explained = varpct))
# }
# 
# # 1) PCA with both depots
# pca_both <- pca_from_vst(
#   vst_mat  = vst_all,
#   meta_df  = meta,
#   title_suffix = "(all samples, both depots)",
#   out_png  = file.path(fig_dir, "PCA_PC1_PC2_ALL_both_DEPOT.png"),
#   out_tsv  = file.path(tab_dir, "PCA_scores_ALL_both_DEPOT.tsv")
# )
# 
# # 2) Depot-specific PCAs
# for (dep in c("Subcutaneous", "Visceral")) {
#   idx <- meta$DEPOT == dep
#   if (sum(idx) >= 10) {
#     pca_from_vst(
#       vst_mat  = vst_all[, meta$SAMPID[idx], drop = FALSE],
#       meta_df  = meta[idx, , drop = FALSE],
#       title_suffix = paste0("(all samples, ", dep, " only)"),
#       out_png  = file.path(fig_dir, paste0("PCA_PC1_PC2_ALL_", gsub(" ", "_", dep), ".png")),
#       out_tsv  = file.path(tab_dir, paste0("PCA_scores_ALL_", gsub(" ", "_", dep), ".tsv"))
#     )
#   } else {
#     message("Skipping PCA for ", dep, " due to small sample size.")
#   }
# }
# 
# 
# 
# 
# 
# 
# 
# # 
# # ## -----------------------------------------------------------------------------
# # ## Define age bands and groups
# # ## -----------------------------------------------------------------------------
# # meta <- meta %>%
# #   mutate(
# #     age_band = case_when(
# #       !is.na(age_years) & age_years < 45 ~ "<45",
# #       !is.na(age_years) & age_years > 55 ~ ">55",
# #       TRUE ~ "50-55_or_unknown"
# #     ),
# #     group4 = case_when(
# #       sex == "Female" & age_band == "<45" ~ "Female_pre",
# #       sex == "Female" & age_band == ">55" ~ "Female_post",
# #       sex == "Male"   & age_band == "<45" ~ "Male_lt45",
# #       sex == "Male"   & age_band == ">55" ~ "Male_gt55",
# #       TRUE ~ NA_character_
# #     )
# #   )
# 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # run_pca_and_plot <- function(count_mat, meta_df, title_suffix, outfile_png, outfile_tsv = NULL) {
# #   dds_tmp <- DESeqDataSetFromMatrix(
# #     countData = round(count_mat),
# #     colData   = as.data.frame(meta_df),
# #     design    = ~ 1
# #   )
# #   keep <- rowSums(counts(dds_tmp) >= 10) >= 20
# #   dds_tmp <- dds_tmp[keep, ]
# #   dds_tmp <- estimateSizeFactors(dds_tmp)
# #   vst_mat <- assay(vst(dds_tmp, blind = TRUE))
# #   
# #   top_var <- head(order(matrixStats::rowVars(vst_mat), decreasing = TRUE), 5000)
# #   pca <- prcomp(t(vst_mat[top_var, ]), center = TRUE, scale. = FALSE)
# #   pc_var <- (pca$sdev^2) / sum(pca$sdev^2) * 100
# #   
# #   pca_df <- as_tibble(pca$x, rownames = "SAMPID") %>%
# #     left_join(meta_df %>% select(SAMPID, group4, sex, DEPOT, age_years), by = "SAMPID")
# #   
# #   p <- ggplot(pca_df, aes(PC1, PC2, color = group4, shape = DEPOT)) +
# #     geom_point(size = 2, alpha = 0.85) +
# #     labs(
# #       title = paste("GTEx adipose PCA", title_suffix),
# #       x = paste0("PC1 (", round(pc_var[1],1), "%)"),
# #       y = paste0("PC2 (", round(pc_var[2],1), "%)")
# #     ) +
# #     theme_bw()
# #   ggsave(outfile_png, p, width = 7, height = 5, dpi = 300)
# #   
# #   if (!is.null(outfile_tsv)) {
# #     readr::write_tsv(pca_df, outfile_tsv)
# #   }
# #   invisible(list(pca = pca, pca_df = pca_df, var_explained = pc_var))
# # }
# 
# 
# # 
# # ## -----------------------------------------------------------------------------
# # ## PCA: both depots together, then per depot
# # ##   Exclude 50–55 or unknown for clearer structure
# # ## -----------------------------------------------------------------------------
# # keep_idx <- meta$age_band %in% c("<45", ">55")
# # meta_pca    <- meta[keep_idx, , drop = FALSE]
# # counts_pca  <- counts_all[, meta_pca$SAMPID, drop = FALSE]
# # 
# # # Both depots
# # run_pca_and_plot(
# #   count_mat   = counts_pca,
# #   meta_df     = meta_pca,
# #   title_suffix = "(both depots)",
# #   outfile_png = file.path(fig_dir, "PCA_PC1_PC2_both_DEPOT.png"),
# #   outfile_tsv = file.path(tab_dir, "PCA_scores_both_DEPOT.tsv")
# # )
# # 
# # # Subcutaneous only
# # idx_sc <- meta_pca$DEPOT == "Subcutaneous"
# # if (sum(idx_sc) >= 10) {
# #   run_pca_and_plot(
# #     count_mat   = counts_pca[, meta_pca$SAMPID[idx_sc], drop = FALSE],
# #     meta_df     = meta_pca[idx_sc, , drop = FALSE],
# #     title_suffix = "(Subcutaneous only)",
# #     outfile_png = file.path(fig_dir, "PCA_PC1_PC2_Subcutaneous.png"),
# #     outfile_tsv = file.path(tab_dir, "PCA_scores_Subcutaneous.tsv")
# #   )
# # } else {
# #   message("Skipping Subcutaneous PCA due to small sample size.")
# # }
# # 
# # # Visceral only
# # idx_vis <- meta_pca$DEPOT == "Visceral"
# # if (sum(idx_vis) >= 10) {
# #   run_pca_and_plot(
# #     count_mat   = counts_pca[, meta_pca$SAMPID[idx_vis], drop = FALSE],
# #     meta_df     = meta_pca[idx_vis, , drop = FALSE],
# #     title_suffix = "(Visceral only)",
# #     outfile_png = file.path(fig_dir, "PCA_PC1_PC2_Visceral.png"),
# #     outfile_tsv = file.path(tab_dir, "PCA_scores_Visceral.tsv")
# #   )
# # } else {
# #   message("Skipping Visceral PCA due to small sample size.")
# # }
# # 
# # ## -----------------------------------------------------------------------------
# # ## DE analysis: Female <45 vs >55, adjusted for DEPOT
# # ##   Excludes 50–55 by design. Runs one model with DEPOT as covariate.
# # ## -----------------------------------------------------------------------------
# # meta_f <- meta %>%
# #   filter(sex == "Female", age_band %in% c("<45", ">55")) %>%
# #   mutate(
# #     group_female = factor(ifelse(age_band == "<45", "Female_pre", "Female_post"),
# #                           levels = c("Female_pre","Female_post")),
# #     DEPOT = factor(DEPOT)
# #   )
# # counts_f <- counts_all[, meta_f$SAMPID, drop = FALSE]
# # 
# # dds_f <- DESeqDataSetFromMatrix(
# #   countData = round(counts_f),
# #   colData   = as.data.frame(meta_f),
# #   design    = ~ DEPOT + group_female
# # )
# # keep_f <- rowSums(counts(dds_f) >= 10) >= 10
# # dds_f <- dds_f[keep_f, ]
# # 
# # dds_f <- DESeq(dds_f)
# # 
# # res_f <- results(dds_f, contrast = c("group_female", "Female_post", "Female_pre"))
# # res_shr <- lfcShrink(dds_f, coef = "group_female_Female_post_vs_Female_pre", type = "apeglm")
# # res_tbl <- as.data.frame(res_shr) %>%
# #   rownames_to_column("gene_id") %>%
# #   arrange(padj)
# # 
# # readr::write_tsv(res_tbl, file.path(tab_dir, "DE_Female_post_vs_pre_adjust_DEPOT.tsv"))
# # 
# # # Quick volcano
# # res_plot <- res_tbl %>%
# #   mutate(sig = ifelse(!is.na(padj) & padj < 0.05, "FDR<0.05", "NS"))
# # p_volc <- ggplot(res_plot, aes(log2FoldChange, -log10(pvalue), color = sig)) +
# #   geom_point(alpha = 0.7, size = 1) +
# #   theme_bw() +
# #   labs(title = "Female >55 vs <45 (adjusted for DEPOT)", x = "log2FC", y = "-log10 p")
# # ggsave(file.path(fig_dir, "Volcano_Female_post_vs_pre.png"), p_volc, width = 7, height = 5, dpi = 300)
# # 
# # ## -----------------------------------------------------------------------------
# # ## Optional: per depot DE runs for females (toggle ON if desired)
# # ## -----------------------------------------------------------------------------
# # run_per_depot_DE <- FALSE
# # if (run_per_depot_DE) {
# #   for (dep in c("Subcutaneous", "Visceral")) {
# #     mf <- meta_f %>% filter(DEPOT == dep)
# #     if (nrow(mf) < 10) {
# #       message("Skip DE in ", dep, " due to small n.")
# #       next
# #     }
# #     cf <- counts_all[, mf$SAMPID, drop = FALSE]
# #     dds_d <- DESeqDataSetFromMatrix(countData = round(cf), colData = as.data.frame(mf), design = ~ group_female)
# #     keep_d <- rowSums(counts(dds_d) >= 10) >= 10
# #     dds_d <- dds_d[keep_d, ]
# #     dds_d <- DESeq(dds_d)
# #     res_d <- lfcShrink(dds_d, coef = "group_female_Female_post_vs_Female_pre", type = "apeglm") %>%
# #       as.data.frame() %>%
# #       rownames_to_column("gene_id") %>%
# #       arrange(padj)
# #     readr::write_tsv(res_d, file.path(tab_dir, paste0("DE_Female_post_vs_pre_", gsub(" ", "_", dep), ".tsv")))
# #     
# #     p_d <- ggplot(res_d %>% mutate(sig = ifelse(!is.na(padj) & padj < 0.05, "FDR<0.05", "NS")),
# #                   aes(log2FoldChange, -log10(pvalue), color = sig)) +
# #       geom_point(alpha = 0.7, size = 1) + theme_bw() +
# #       labs(title = paste0("Female >55 vs <45 in ", dep), x = "log2FC", y = "-log10 p")
# #     ggsave(file.path(fig_dir, paste0("Volcano_Female_post_vs_pre_", gsub(" ", "_", dep), ".png")),
# #            p_d, width = 7, height = 5, dpi = 300)
# #   }
# # }
# # 
# # ## -----------------------------------------------------------------------------
# # ## Save session info
# # ## -----------------------------------------------------------------------------
# # saveRDS(list(meta = meta, meta_f = meta_f), file = file.path(rds_dir, "meta_objects.rds"))
# # writeLines(c(capture.output(sessionInfo())), con = file.path(out_dir, "sessionInfo.txt"))
# # 
# # message("Done. Figures in: ", normalizePath(fig_dir))
# # message("Tables in:  ", normalizePath(tab_dir))




# 
# ## -----------------------------------------------------------------------------
# ## 9) DESeq2: 40-49 vs 50-59 by sex and tissue (parallel, cached)
# ## -----------------------------------------------------------------------------
# ## -----------------------------------------------------------------------------
# ## 9) DESeq2: 40-49 vs 50-59 by sex and tissue (parallel, cached)
# ## -----------------------------------------------------------------------------
# 
# if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply")
# library(DESeq2)
# library(dplyr)
# library(readr)
# library(future.apply)
# 
# # clean age_bin_label once
# meta$age_bin_label <- gsub("\\s*–\\s*|\\s*-\\s*", "-", meta$age_bin_label)
# 
# age_bins <- c("40-49", "50-59")
# 
# meta_de <- meta %>%
#   filter(
#     DEPOT %in% c("Subcutaneous", "Visceral"),
#     sex   %in% c("Female", "Male"),
#     age_bin_label %in% age_bins
#   ) %>%
#   mutate(age_group = factor(age_bin_label, levels = age_bins))
# 
# # gene_map already has gene_id_clean, gene_symbol
# gene_map_de <- gene_map %>%
#   mutate(gene_id_clean = sub("\\..*$", "", gene_id_clean)) %>%
#   distinct(gene_id_clean, gene_symbol, .keep_all = TRUE)
# 
# # cache heavy objects so futures only get a path, not 500MB globals
# deseq_cache_path <- file.path(rds_dir, "deseq2_cache_40vs50.rds")
# saveRDS(
#   list(
#     counts_all = counts_all,
#     meta_de    = meta_de,
#     gene_map_de = gene_map_de,
#     tab_dir    = tab_dir
#   ),
#   deseq_cache_path
# )
# 
# combos <- tidyr::expand_grid(
#   sex    = c("Female", "Male"),
#   tissue = c("Subcutaneous", "Visceral")
# )
# 
# run_one_deseq <- function(sx, dep, cache_path) {
#   env <- readRDS(cache_path)
#   counts_all <- env$counts_all
#   meta_de    <- env$meta_de
#   gene_map_de <- env$gene_map_de
#   tab_dir    <- env$tab_dir
#   
#   meta_s <- meta_de %>% filter(sex == sx, DEPOT == dep)
#   if (nrow(meta_s) < 10) {
#     message("Skip DESeq2 for ", sx, " / ", dep, " (n < 10)")
#     return(NULL)
#   }
#   
#   cts_s <- counts_all[, meta_s$SAMPID, drop = FALSE]
#   keep  <- rowSums(cts_s >= 10) >= 20
#   cts_s <- cts_s[keep, , drop = FALSE]
#   if (nrow(cts_s) == 0) {
#     message("No genes pass filter for ", sx, " / ", dep)
#     return(NULL)
#   }
#   
#   dds <- DESeqDataSetFromMatrix(
#     countData = round(cts_s),
#     colData   = as.data.frame(meta_s),
#     design    = ~ age_group
#   )
#   dds <- DESeq(dds, quiet = TRUE)
#   
#   res <- results(dds, contrast = c("age_group", "50-59", "40-49"))
#   res <- as.data.frame(res)
#   
#   res$gene_id       <- rownames(res)
#   res$gene_id_clean <- sub("\\..*$", "", res$gene_id)
#   res$sex           <- sx
#   res$tissue        <- dep
#   
#   res <- res %>%
#     left_join(
#       gene_map_de %>% dplyr::select(gene_id_clean, gene_symbol),
#       by = "gene_id_clean"
#     ) %>%
#     relocate(gene_id, gene_id_clean, gene_symbol,
#              sex, tissue, .before = baseMean)
#   
#   key <- paste(sx, gsub(" ", "_", dep), sep = "_")
#   out_cmp <- file.path(tab_dir, paste0("DESeq2_40vs50_", key, ".csv"))
#   write_csv(res, out_cmp)
#   message("DESeq2 written: ", out_cmp, " (", nrow(res), " genes)")
#   
#   list(key = key, res = res)
# }
# 
# n_workers <- min(nrow(combos), future::availableCores(), 8)
# future::plan(future::multisession, workers = n_workers)
# 
# res_list_raw <- future_lapply(
#   seq_len(nrow(combos)),
#   function(i, cache_path) {
#     sx  <- combos$sex[i]
#     dep <- combos$tissue[i]
#     run_one_deseq(sx, dep, cache_path)
#   },
#   cache_path = deseq_cache_path,
#   future.seed = TRUE
# )
# 
# future::plan(future::sequential)
# 
# valid <- Filter(Negate(is.null), res_list_raw)
# 
# if (length(valid)) {
#   res_list <- setNames(
#     lapply(valid, `[[`, "res"),
#     vapply(valid, `[[`, character(1), "key")
#   )
#   
#   res_all <- bind_rows(res_list, .id = "comparison_id")
#   
#   out_all <- file.path(tab_dir, "DESeq2_40vs50_ALL_bySexTissue.csv")
#   write_csv(res_all, out_all)
#   
#   res_sig <- res_all %>%
#     filter(!is.na(padj), padj < 0.05) %>%
#     arrange(sex, tissue, padj)
#   
#   out_sig <- file.path(tab_dir, "DESeq2_40vs50_ALL_bySexTissue_significant.csv")
#   write_csv(res_sig, out_sig)
#   
#   message(
#     "DESeq2 master files:",
#     "\n  ", out_all,
#     "\n  ", out_sig
#   )
# } else {
#   message("No DESeq2 results to combine.")
# }
# 


## -----------------------------------------------------------------------------
## 9) DESeq2: 40-49 vs 50-59 by sex and tissue (parallel, cached)
## -----------------------------------------------------------------------------

if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply")
library(DESeq2)
library(dplyr)
library(readr)
library(future.apply)

meta$age_bin_label <- gsub("\\s*–\\s*|\\s*-\\s*", "-", meta$age_bin_label)

age_bins <- c("40-49", "50-59")

meta_de <- meta %>%
  filter(
    DEPOT %in% c("Subcutaneous", "Visceral"),
    sex   %in% c("Female", "Male"),
    age_bin_label %in% age_bins
  ) %>%
  mutate(age_group = factor(age_bin_label, levels = age_bins))

gene_map_de <- gene_map %>%
  mutate(gene_id_clean = sub("\\..*$", "", gene_id_clean)) %>%
  distinct(gene_id_clean, gene_symbol, .keep_all = TRUE)

deseq_cache_path <- file.path(rds_dir, "deseq2_cache_40vs50.rds")
saveRDS(
  list(
    counts_all  = counts_all,
    meta_de     = meta_de,
    gene_map_de = gene_map_de,
    tab_dir     = tab_dir
  ),
  deseq_cache_path
)

combos <- tidyr::expand_grid(
  sex    = c("Female", "Male"),
  tissue = c("Subcutaneous", "Visceral")
)

run_one_deseq <- function(sx, dep, cache_path) {
  env         <- readRDS(cache_path)
  counts_all  <- env$counts_all
  meta_de     <- env$meta_de
  gene_map_de <- env$gene_map_de
  tab_dir     <- env$tab_dir
  
  meta_s <- meta_de %>% filter(sex == sx, DEPOT == dep)
  if (nrow(meta_s) < 10) {
    message("Skip DESeq2 for ", sx, " / ", dep, " (n < 10)")
    return(NULL)
  }
  
  cts_s <- counts_all[, meta_s$SAMPID, drop = FALSE]
  keep  <- rowSums(cts_s >= 10) >= 20
  cts_s <- cts_s[keep, , drop = FALSE]
  if (nrow(cts_s) == 0) {
    message("No genes pass filter for ", sx, " / ", dep)
    return(NULL)
  }
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(cts_s),
    colData   = as.data.frame(meta_s),
    design    = ~ age_group
  )
  dds <- DESeq(dds, quiet = TRUE)
  
  res <- results(dds, contrast = c("age_group", "50-59", "40-49"))
  res <- as.data.frame(res)
  
  res$gene_id       <- rownames(res)
  res$gene_id_clean <- sub("\\..*$", "", res$gene_id)
  res$sex           <- sx
  res$tissue        <- dep
  
  res <- res %>%
    left_join(
      gene_map_de %>% dplyr::select(gene_id_clean, gene_symbol),
      by = "gene_id_clean"
    ) %>%
    relocate(gene_id, gene_id_clean, gene_symbol,
             sex, tissue, .before = baseMean)
  
  key <- paste(sx, gsub(" ", "_", dep), sep = "_")
  
  out_cmp     <- file.path(tab_dir, paste0("DESeq2_40vs50_", key, ".csv"))
  out_cmp_sig <- file.path(tab_dir, paste0("DESeq2_40vs50_", key, "_sig.csv"))
  
  write_csv(res, out_cmp)
  
  res_sig <- res %>%
    filter(!is.na(padj), padj < 0.05)
  write_csv(res_sig, out_cmp_sig)
  
  message("DESeq2 written: ", out_cmp,
          " and ", out_cmp_sig,
          " (", nrow(res), " total, ",
          nrow(res_sig), " sig)")
  
  list(key = key, res = res)
}

n_workers <- min(nrow(combos), future::availableCores(), 8)
future::plan(future::multisession, workers = n_workers)

res_list_raw <- future_lapply(
  seq_len(nrow(combos)),
  function(i, cache_path) {
    sx  <- combos$sex[i]
    dep <- combos$tissue[i]
    run_one_deseq(sx, dep, cache_path)
  },
  cache_path = deseq_cache_path,
  future.seed = TRUE
)

future::plan(future::sequential)

valid <- Filter(Negate(is.null), res_list_raw)

if (length(valid)) {
  res_list <- setNames(
    lapply(valid, `[[`, "res"),
    vapply(valid, `[[`, character(1), "key")
  )
  
  res_all <- bind_rows(res_list, .id = "comparison_id")
  
  out_all <- file.path(tab_dir, "DESeq2_40vs50_ALL_bySexTissue.csv")
  write_csv(res_all, out_all)
  
  res_sig_all <- res_all %>%
    filter(!is.na(padj), padj < 0.05) %>%
    arrange(sex, tissue, padj)
  
  out_sig_all <- file.path(tab_dir, "DESeq2_40vs50_ALL_bySexTissue_significant.csv")
  write_csv(res_sig_all, out_sig_all)
  
  message(
    "DESeq2 master files:",
    "\n  ", out_all,
    "\n  ", out_sig_all
  )
} else {
  message("No DESeq2 results to combine.")
}







## -----------------------------------------------------------------------------
## 10) richR annotations (GO BP and KEGG)
## -----------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db", ask = FALSE)
if (!requireNamespace("reactome.db", quietly = TRUE)) BiocManager::install("reactome.db", ask = FALSE)
if (!requireNamespace("richR", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("guokai8/richR")
}

library(richR)

hsa_go <- buildAnnot(
  species = "human",
  keytype = "SYMBOL",
  anntype = "GO"
)

hsa_ko <- buildAnnot(
  species = "human",
  keytype = "SYMBOL",
  anntype = "KEGG",
  builtin = FALSE
)

saveRDS(hsa_go, file.path(rds_dir, "richR_hsa_go_BP.rds"))
saveRDS(hsa_ko, file.path(rds_dir, "richR_hsa_kegg.rds"))



## -----------------------------------------------------------------------------
## 11) richR enrichment per comparison (GO BP and KEGG)
## -----------------------------------------------------------------------------

library(dplyr)
library(readr)
library(ggplot2)

hsa_go <- readRDS(file.path(rds_dir, "richR_hsa_go_BP.rds"))
hsa_ko <- readRDS(file.path(rds_dir, "richR_hsa_kegg.rds"))

deg_sig_files <- list.files(
  tab_dir,
  pattern = "^DESeq2_40vs50_.*_sig\\.csv$",
  full.names = TRUE
)

go_list  <- list()
kegg_list <- list()

for (f in deg_sig_files) {
  fn <- basename(f)
  cmp <- gsub("^DESeq2_40vs50_|_sig\\.csv$", "", fn)  # e.g. Female_Subcutaneous
  
  dat_sig <- read_csv(f, show_col_types = FALSE)
  if (!"gene_symbol" %in% names(dat_sig)) {
    message("Skip ", cmp, " (no gene_symbol column)")
    next
  }
  
  genes <- dat_sig %>%
    filter(!is.na(padj), padj < 0.05, !is.na(gene_symbol), gene_symbol != "") %>%
    pull(gene_symbol) %>%
    unique()
  
  if (length(genes) < 5) {
    message("Skip enrichment for ", cmp, " (too few DEGs: ", length(genes), ")")
    next
  }
  
  ## GO BP
  GO.res <- try(
    richGO(
      genes,
      godata   = hsa_go,
      ontology = "BP"
    ),
    silent = TRUE
  )
  if (!inherits(GO.res, "try-error") && nrow(GO.res)) {
    GO.sig <- GO.res %>%
      filter(Padj < 0.05) %>%
      arrange(Padj)
    
    if (nrow(GO.sig)) {
      go_list[[cmp]] <- GO.sig
      
      out_go <- file.path(tab_dir, paste0("richR_GO_BP_", cmp, ".csv"))
      write_csv(GO.sig@result, out_go)
      
      # dot plot top 20
      go_top <- GO.sig[seq_len(min(20, nrow(GO.sig))), ]
      ggdot(
        go_top,
        top      = nrow(go_top),
        usePadj  = TRUE,
        filename = file.path(fig_dir, paste0("richR_GO_BP_dot_top20_", cmp, ".pdf")),
        width    = 8,
        height   = 6
      )
    }
  }
  
  ## KEGG
  KEGG.res <- try(
    richKEGG(
      genes,
      kodata  = hsa_ko,
      pvalue  = 0.05,
      builtin = FALSE
    ),
    silent = TRUE
  )
  if (!inherits(KEGG.res, "try-error") && nrow(KEGG.res)) {
    KEGG.sig <- KEGG.res %>%
      filter(Padj < 0.05) %>%
      arrange(Padj)
    
    if (nrow(KEGG.sig)) {
      kegg_list[[cmp]] <- KEGG.sig
      
      out_kegg <- file.path(tab_dir, paste0("richR_KEGG_", cmp, ".csv"))
      write_csv(KEGG.sig@result, out_kegg)
      
      # dot plot top 20
      kegg_top <- KEGG.sig[seq_len(min(20, nrow(KEGG.sig))), ]
      ggdot(
        kegg_top,
        top      = nrow(kegg_top),
        usePadj  = TRUE,
        filename = file.path(fig_dir, paste0("richR_KEGG_dot_top20_", cmp, ".pdf")),
        width    = 8,
        height   = 6
      )
    }
  }
}

message("Per comparison enrichment complete. GO sets: ", length(go_list),
        " KEGG sets: ", length(kegg_list))

# 
# 
# ## -----------------------------------------------------------------------------
# ## 12) Cross comparison enrichment: compareResult + comparedot
# ## -----------------------------------------------------------------------------
# 
# # GO BP
# if (length(go_list) >= 2) {
#   go_comp <- compareResult(go_list)
#   
#   # top 20 terms combined
#   comparedot(
#     go_comp,
#     top      = 20,
#     usePadj  = TRUE,
#     filename = file.path(fig_dir, "richR_GO_BP_comparedot_top20.pdf"),
#     width    = 9,
#     height   = 6
#   )
#   
#   out_go_comp <- file.path(tab_dir, "richR_GO_BP_compareResult.csv")
#   write_csv(go_comp, out_go_comp)
#   
#   message("GO BP compareResult written: ", out_go_comp)
# }
# 
# # KEGG
# if (length(kegg_list) >= 2) {
#   kegg_comp <- compareResult(kegg_list)
#   
#   comparedot(
#     kegg_comp,
#     top      = 20,
#     usePadj  = TRUE,
#     filename = file.path(fig_dir, "richR_KEGG_comparedot_top20.pdf"),
#     width    = 9,
#     height   = 6
#   )
#   
#   out_kegg_comp <- file.path(tab_dir, "richR_KEGG_compareResult.csv")
#   write_csv(kegg_comp, out_kegg_comp)
#   
#   message("KEGG compareResult written: ", out_kegg_comp)
# }
