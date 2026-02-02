################################################################################
## Sex Interaction DE for Age-Bin Comparisons (Menopause-Enriched Signal)
##
## Goal
## - Separate generic aging signal from female-specific (potential menopause-related)
##   signal by using males as a comparator.
##
## Model (Subcutaneous adipose only)
##   design = ~ SMCENTER + sex + age_group + sex:age_group
##
## With baseline levels:
## - sex = Male (baseline), Female
## - age_group = young (baseline), old
##
## Then:
## - "age_group_old_vs_young" is the male aging effect (old vs young).
## - "sexFemale.age_groupold" is the *interaction* (female differential vs male).
## - Female effect (old vs young) = age_group + interaction (via contrast=list()).
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
  readr, dplyr, tidyr, stringr, tibble, purrr, glue, magrittr,
  ggplot2, DESeq2, AnnotationDbi, org.Hs.eg.db
)

source("R/utils.R")

## -----------------------------------------------------------------------------
## Paths
## -----------------------------------------------------------------------------
in_dir <- "GTExDatav10"
gct_sc <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
sattr <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

out_dir <- "GTEx_v10_AT_analysis_out"
fig_dir <- file.path(out_dir, "figs", "sex_interaction_age_bins")
tab_dir <- file.path(out_dir, "tables", "sex_interaction_age_bins")
ensure_dirs(c(fig_dir, tab_dir))

## -----------------------------------------------------------------------------
## Config
## -----------------------------------------------------------------------------
min_n_per_sex_group <- 5

comparisons <- list(
  list(name = "40-49_vs_50-59", young = "40-49", old = "50-59"),
  list(name = "30-39_vs_60-69", young = "30-39", old = "60-69"),
  list(name = "30-39_vs_50-59", young = "30-39", old = "50-59"),
  list(name = "40-49_vs_60-69", young = "40-49", old = "60-69")
)

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------
as_res_df <- function(res_obj, gene_symbol_lookup, comparison_id, effect) {
  as.data.frame(res_obj) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::mutate(
      gene_symbol = gene_symbol_lookup[match(.data$gene_id, names(gene_symbol_lookup))],
      comparison = comparison_id,
      effect = effect
    ) %>%
    dplyr::relocate(comparison, effect, gene_id, gene_symbol)
}

pick_one_coef <- function(rn, candidates, label) {
  if (length(candidates) != 1) {
    stop(
      "Expected exactly 1 coefficient for ", label, ", got ", length(candidates),
      "\nresultsNames(dds) = ", paste(rn, collapse = ", "),
      "\nCandidates = ", paste(candidates, collapse = ", ")
    )
  }
  candidates[[1]]
}

write_sig_list <- function(df, out_path, padj_cutoff = 0.05, lfc_min = 0) {
  sig <- df %>%
    dplyr::filter(!is.na(.data$padj), .data$padj < padj_cutoff, abs(.data$log2FoldChange) >= lfc_min) %>%
    dplyr::arrange(.data$padj) %>%
    dplyr::mutate(gene_symbol = dplyr::if_else(is.na(.data$gene_symbol) | .data$gene_symbol == "", .data$gene_id, .data$gene_symbol))
  readr::write_lines(unique(sig$gene_symbol), out_path)
  nrow(sig)
}

## -----------------------------------------------------------------------------
## Load and prepare data
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("SEX INTERACTION DE (SUBCUTANEOUS)\n")
cat("===========================================\n\n")

cat("Loading sample metadata...\n")
sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything()) %>%
  dplyr::mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])

subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  dplyr::rename_with(~ gsub("^X", "", .), dplyr::everything())

meta_full <- sample_attr %>%
  dplyr::left_join(subject_pheno, by = "SUBJID")

cat("Loading subcutaneous counts...\n")
gct_list_sc <- read_gct_v12(gct_sc)

common_samples <- intersect(colnames(gct_list_sc$counts), meta_full$SAMPID)
meta <- meta_full %>%
  dplyr::filter(SAMPID %in% common_samples, SMTS == "Adipose Tissue", SMTSD == "Adipose - Subcutaneous") %>%
  dplyr::mutate(
    sex = dplyr::case_when(
      SEX == 1 ~ "Male",
      SEX == 2 ~ "Female",
      TRUE ~ NA_character_
    ),
    age_bin_label = AGE
  ) %>%
  dplyr::filter(!is.na(sex), !is.na(age_bin_label), !is.na(SMCENTER)) %>%
  dplyr::bind_cols(parse_age_bin(.$AGE)) %>%
  derive_age_numeric()

counts_sc <- gct_list_sc$counts[, intersect(meta$SAMPID, colnames(gct_list_sc$counts)), drop = FALSE]
counts_sc <- collapse_dupe_rowsum(counts_sc)

gene_symbol_lookup <- gct_list_sc$gene_symbol
names(gene_symbol_lookup) <- gct_list_sc$gene_id

## Gene -> chromosome lookup (used to exclude chrY from "menopause-enriched" sets)
cat("Building Ensembl chromosome lookup (org.Hs.eg.db)...\n")
gene_ids_clean <- unique(stringr::str_remove(rownames(counts_sc), "\\..*$"))
annot <- AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = gene_ids_clean,
  columns = c("CHR"),
  keytype = "ENSEMBL"
)
chr_lookup <- annot %>%
  dplyr::filter(!is.na(.data$CHR) & .data$CHR != "") %>%
  dplyr::group_by(.data$ENSEMBL) %>%
  dplyr::summarise(CHR = .data$CHR[[1]], .groups = "drop") %>%
  tibble::deframe()

cat(sprintf("Loaded %d subcutaneous samples (%d Female, %d Male)\n",
            nrow(meta), sum(meta$sex == "Female"), sum(meta$sex == "Male")))
cat(sprintf("Counts matrix: %d genes x %d samples\n\n", nrow(counts_sc), ncol(counts_sc)))

## -----------------------------------------------------------------------------
## Run comparisons
## -----------------------------------------------------------------------------
summary_rows <- list()

for (cmp in comparisons) {
  comparison_id <- cmp$name
  young_age <- cmp$young
  old_age <- cmp$old

  cat("\n===========================================\n")
  cat(glue::glue("COMPARISON: {young_age} vs {old_age} (both sexes)\n"))
  cat("===========================================\n\n")

  meta_de <- meta %>%
    dplyr::filter(age_bin_label %in% c(young_age, old_age)) %>%
    dplyr::mutate(
      sex = factor(sex, levels = c("Male", "Female")),
      age_group = factor(age_bin_label, levels = c(young_age, old_age)),
      SMCENTER = factor(SMCENTER)
    ) %>%
    dplyr::filter(!is.na(sex), !is.na(age_group), !is.na(SMCENTER))

  n_tbl <- meta_de %>%
    dplyr::count(sex, age_group, name = "n") %>%
    tidyr::complete(sex, age_group, fill = list(n = 0L))
  print(n_tbl)

  if (any(n_tbl$n < min_n_per_sex_group)) {
    warning(glue::glue(
      "Skipping {comparison_id}: insufficient samples (min {min_n_per_sex_group} per sex x age_group)."
    ))
    next
  }

  ## Subset counts to cohort and keep ordering aligned with colData
  cohort_samples <- meta_de$SAMPID
  counts_de <- counts_sc[, cohort_samples, drop = FALSE]

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts_de,
    colData = meta_de,
    design = ~ SMCENTER + sex + age_group + sex:age_group
  )

  keep <- rowSums(DESeq2::counts(dds)) >= 10
  dds <- dds[keep, ]
  cat(sprintf("After filtering: %d genes\n", nrow(dds)))

  cat("Running DESeq2...\n")
  dds <- DESeq2::DESeq(dds)

  rn <- DESeq2::resultsNames(dds)
  cat("\nresultsNames(dds):\n")
  print(rn)

  ## Coefficient discovery (assumes baseline sex = Male, baseline age_group = young)
  age_candidates <- grep("^age_group", rn, value = TRUE)
  age_coef <- pick_one_coef(rn, age_candidates, "age_group (male effect)")

  sex_main <- grep("^sex_", rn, value = TRUE)
  interaction_candidates <- rn[stringr::str_detect(rn, "sex") & stringr::str_detect(rn, "age_group")]
  interaction_candidates <- setdiff(interaction_candidates, c(age_coef, sex_main))
  interaction_coef <- pick_one_coef(rn, interaction_candidates, "sex:age_group interaction")

  ## Effects
  res_male <- DESeq2::results(dds, name = age_coef)
  res_interaction <- DESeq2::results(dds, name = interaction_coef)
  res_female <- DESeq2::results(dds, contrast = list(c(age_coef, interaction_coef)))

  ## Persist (full)
  male_df <- as_res_df(res_male, gene_symbol_lookup, comparison_id, "male_old_vs_young")
  inter_df <- as_res_df(res_interaction, gene_symbol_lookup, comparison_id, "female_minus_male_interaction")
  female_df <- as_res_df(res_female, gene_symbol_lookup, comparison_id, "female_old_vs_young")

  readr::write_tsv(male_df %>% dplyr::arrange(.data$padj),
                   file.path(tab_dir, glue::glue("{comparison_id}__male_effect.tsv")))
  readr::write_tsv(inter_df %>% dplyr::arrange(.data$padj),
                   file.path(tab_dir, glue::glue("{comparison_id}__interaction.tsv")))
  readr::write_tsv(female_df %>% dplyr::arrange(.data$padj),
                   file.path(tab_dir, glue::glue("{comparison_id}__female_effect.tsv")))

  ## Menopause-enriched heuristic set:
  ## - female DE significant
  ## - interaction significant (female differs from male)
  ## - male DE not significant (to reduce generic aging hits)
  combined <- female_df %>%
    dplyr::select(gene_id, gene_symbol, log2FoldChange_female = log2FoldChange, padj_female = padj) %>%
    dplyr::left_join(
      male_df %>% dplyr::select(gene_id, log2FoldChange_male = log2FoldChange, padj_male = padj),
      by = "gene_id"
    ) %>%
    dplyr::left_join(
      inter_df %>% dplyr::select(gene_id, log2FoldChange_interaction = log2FoldChange, padj_interaction = padj),
      by = "gene_id"
    ) %>%
    dplyr::mutate(
      chr = chr_lookup[stringr::str_remove(.data$gene_id, "\\..*$")]
    )

  menopause_enriched <- combined %>%
    dplyr::filter(
      !is.na(.data$padj_female), .data$padj_female < 0.05,
      !is.na(.data$padj_interaction), .data$padj_interaction < 0.05,
      is.na(.data$padj_male) | .data$padj_male >= 0.2,
      ## Exclude chrY artifacts (e.g., Y-linked genes in female effects due to cross-mapping / QC).
      ## org.Hs.eg.db chromosome mapping can be incomplete for some Ensembl IDs, so also fall back
      ## to a conservative symbol suffix heuristic.
      !(
        (!is.na(.data$chr) & .data$chr == "Y") |
          (!is.na(.data$gene_symbol) & stringr::str_detect(.data$gene_symbol, "Y$"))
      )
    ) %>%
    dplyr::arrange(.data$padj_interaction, .data$padj_female)

  men_out <- file.path(tab_dir, glue::glue("{comparison_id}__menopause_enriched.tsv"))
  readr::write_tsv(menopause_enriched, men_out)

  ## Signature lists (symbols)
  n_male_sig <- write_sig_list(
    male_df,
    file.path(tab_dir, glue::glue("{comparison_id}__male_sig_symbols_padj0.05.txt")),
    padj_cutoff = 0.05
  )
  n_female_sig <- write_sig_list(
    female_df,
    file.path(tab_dir, glue::glue("{comparison_id}__female_sig_symbols_padj0.05.txt")),
    padj_cutoff = 0.05
  )
  n_inter_sig <- write_sig_list(
    inter_df,
    file.path(tab_dir, glue::glue("{comparison_id}__interaction_sig_symbols_padj0.05.txt")),
    padj_cutoff = 0.05
  )

  cat(sprintf("\nSignificant (padj<0.05): male=%d, female=%d, interaction=%d\n",
              n_male_sig, n_female_sig, n_inter_sig))
  cat(sprintf("Menopause-enriched (heuristic): %d genes\n", nrow(menopause_enriched)))

  ## Plot: female vs male LFC, highlight interaction hits
  plot_df <- combined %>%
    dplyr::mutate(
      is_interaction_sig = !is.na(.data$padj_interaction) & .data$padj_interaction < 0.05
    )

  p_scatter <- ggplot2::ggplot(plot_df, ggplot2::aes(.data$log2FoldChange_male, .data$log2FoldChange_female)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, alpha = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.2, alpha = 0.4) +
    ggplot2::geom_point(ggplot2::aes(color = .data$is_interaction_sig), alpha = 0.45, size = 1.2, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "#c51b7d")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = glue::glue("Female vs Male log2FC ({young_age} -> {old_age})"),
      subtitle = "Colored by interaction padj<0.05",
      x = "Male log2 fold change (old vs young)",
      y = "Female log2 fold change (old vs young)",
      color = "Interaction sig"
    )

  ggplot2::ggsave(
    filename = file.path(fig_dir, glue::glue("{comparison_id}__female_vs_male_lfc_scatter.png")),
    plot = p_scatter,
    width = 6.5,
    height = 5,
    dpi = 300
  )

  summary_rows[[comparison_id]] <- tibble::tibble(
    comparison = comparison_id,
    young = young_age,
    old = old_age,
    n_total = nrow(meta_de),
    n_female = sum(meta_de$sex == "Female"),
    n_male = sum(meta_de$sex == "Male"),
    n_sig_male = n_male_sig,
    n_sig_female = n_female_sig,
    n_sig_interaction = n_inter_sig,
    n_menopause_enriched = nrow(menopause_enriched)
  )
}

if (length(summary_rows)) {
  summary_tbl <- dplyr::bind_rows(summary_rows) %>%
    dplyr::arrange(.data$comparison)
  readr::write_tsv(summary_tbl, file.path(tab_dir, "sex_interaction_summary.tsv"))
  cat("\nWrote summary table: ", file.path(tab_dir, "sex_interaction_summary.tsv"), "\n", sep = "")
} else {
  cat("\nNo comparisons produced results (insufficient sample sizes or earlier errors).\n")
}
