################################################################################
## Validate Menopause Signature with GSE86244 Dataset
##
## GSE86244 (Shan et al.) contains actual pre/post-menopausal status
## This validates our GTEx-derived 81-gene signature against an independent dataset
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
  data.table, readr, dplyr, tidyr, stringr, tibble, purrr, glue, magrittr,
  ggplot2, patchwork, limma, GEOquery, AnnotationDbi, org.Hs.eg.db
)

## -----------------------------------------------------------------------------
## Source utility functions
## -----------------------------------------------------------------------------
source("R/utils.R")

cat("\nNOTE:\n")
cat("This script previously assumed GSE86244 could be analyzed via a GEOquery ExpressionSet.\n")
cat("In practice, GSE86244 provides RNA-seq count matrices as GEO supplementary files.\n")
cat("Use the updated pipeline instead:\n")
cat("  Rscript R/validate_with_gse86244_counts_age_proxy.R\n\n")

source("R/validate_with_gse86244_counts_age_proxy.R")
quit(save = "no", status = 0)

## -----------------------------------------------------------------------------
## Setup
## -----------------------------------------------------------------------------
out_dir  <- "GTEx_v10_AT_analysis_out"
fig_dir  <- file.path(out_dir, "figs", "gse86244_validation")
tab_dir  <- file.path(out_dir, "tables", "gse86244_validation")
ensure_dirs(c(out_dir, fig_dir, tab_dir))

cat("\n===========================================\n")
cat("GSE86244 VALIDATION ANALYSIS\n")
cat("===========================================\n\n")

## -----------------------------------------------------------------------------
## Load our menopause signature (from Approach 1)
## -----------------------------------------------------------------------------
cat("Loading menopause signature from Approach 1...\n")

sig_file <- file.path("GTEx_v10_AT_analysis_out/tables/approach1_improved_age_proxy/menopause_signature_gene_list.txt")
if (!file.exists(sig_file)) {
  stop("Signature file not found. Run Approach 1 first.")
}

signature_genes <- read_lines(sig_file) %>% str_trim()
cat(sprintf("Loaded %d signature genes from GTEx analysis\n", length(signature_genes)))

## -----------------------------------------------------------------------------
## GSE86244 Dataset Information
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("GSE86244 DATASET INFORMATION\n")
cat("===========================================\n\n")

cat("Title: Transcriptional profiling of adipose tissue from pre- and post-menopausal women\n")
cat("Authors: Shan et al.\n")
cat("Samples: 15 women (12 premenopausal, 3 postmenopausal)\n")
cat("Tissue: Subcutaneous adipose tissue\n")
cat("Platform: GPL570 [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array\n\n")

## -----------------------------------------------------------------------------
## Try to download GSE86244 data from GEO
## -----------------------------------------------------------------------------
cat("Attempting to download GSE86244 from GEO...\n")

gse_dir <- file.path(out_dir, "GSE86244_data")
ensure_dirs(gse_dir)

gse_data_file <- file.path(gse_dir, "GSE86244_series_matrix.txt.gz")

tryCatch({
  cat("Downloading GSE86244 (this may take a few minutes)...\n")

  # Try to get series matrix
  gse <- GEOquery::getGEO("GSE86244", GSEMatrix = TRUE, getGPL = FALSE)

  # Extract expression data and phenotype
  gse_expr <- exprs(gse[[1]])
  gse_pheno <- pData(gse[[1]])

  # Save for future use
  saveRDS(list(expr = gse_expr, pheno = gse_pheno), file.path(gse_dir, "GSE86244_data.rds"))

  cat(sprintf("Successfully downloaded GSE86244: %d samples, %d probes\n",
              ncol(gse_expr), nrow(gse_expr)))

}, error = function(e) {
  cat("GEO download failed:\n")
  cat(sprintf("  Error: %s\n", e$message))
  cat("\nFALLBACK: Using simulated comparison for demonstration\n")
  cat("In production, manually download from:\n")
  cat("  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86244\n\n")
})

## -----------------------------------------------------------------------------
## Check if we have real GSE86244 data or need fallback
## -----------------------------------------------------------------------------
if (file.exists(file.path(gse_dir, "GSE86244_data.rds"))) {
  gse_data <- readRDS(file.path(gse_dir, "GSE86244_data.rds"))
  gse_expr <- gse_data$expr
  gse_pheno <- gse_data$pheno

  cat("\nUsing downloaded GSE86244 data...\n")
  cat(sprintf("Expression matrix: %d probes x %d samples\n",
              nrow(gse_expr), ncol(gse_expr)))

  # Check if we have valid data
  if (nrow(gse_expr) == 0 || ncol(gse_expr) == 0) {
    cat("Warning: Invalid expression data from GEO download\n")
    cat("The GEO query returned 0 probes. This is a known issue with GEOquery.\n")
    cat("Falling back to template mode.\n")
    gse_data_valid <- FALSE
  } else {
    gse_data_valid <- TRUE

    # Show sample characteristics
    cat("Sample characteristics:\n")
    tryCatch({
      print(gse_pheno %>% dplyr::select(dplyr::contains(c("title", "source", "menopause", "age"))))
    }, error = function(e) {
      print(head(gse_pheno))
    })

    ## -----------------------------------------------------------------------------
    ## Analyze GSE86244: Differential expression by menopause status
    ## -----------------------------------------------------------------------------
    cat("\nPerforming differential expression analysis on GSE86244...\n")

    # Extract menopause status from title/source_name
    # Format: "premenopausal woman" or "postmenopausal woman"
    gse_pheno <- gse_pheno %>%
      mutate(
        menopause_status = case_when(
          str_detect(title, "postmenopausal|post-menopausal") ~ "Post",
          str_detect(title, "premenopausal|pre-menopausal") ~ "Pre",
          TRUE ~ NA_character_
        )
      )

  # Check if we have valid menopause status
  if (sum(!is.na(gse_pheno$menopause_status)) < nrow(gse_pheno) * 0.8) {
    cat("Warning: Could not reliably extract menopause status from sample titles\n")
    cat("Sample titles:\n")
    print(gse_pheno$title)
  }

  # Count samples per group
  status_counts <- gse_pheno %>%
    count(menopause_status) %>%
    filter(!is.na(menopause_status))

  cat("\nMenopause status distribution:\n")
  print(status_counts)

  # Map probes to gene symbols
  cat("\nMapping probes to gene symbols...\n")

  # Get platform annotation
  tryCatch({
    gpl <- GEOquery::getGEO("GPL570")

    # Extract gene symbol mapping
    gpl_tbl <- gpl@dataTable@table
    # Typical GPL570 columns include "ID" (probe) and "Gene Symbol".
    sym_col <- intersect(c("Gene Symbol", "GENE_SYMBOL", "Gene symbol"), colnames(gpl_tbl))
    if (!length(sym_col)) {
      stop("Could not find a gene symbol column in GPL570 annotation table.")
    }
    sym_col <- sym_col[[1]]
    probe_map <- gpl_tbl %>%
      dplyr::transmute(
        Probe = .data$ID,
        gene_symbol_raw = .data[[sym_col]]
      ) %>%
      dplyr::filter(!is.na(Probe), Probe != "", !is.na(gene_symbol_raw), gene_symbol_raw != "") %>%
      dplyr::mutate(
        # GPL annotations sometimes contain multiple symbols separated by delimiters.
        gene_symbol = stringr::str_split_fixed(gene_symbol_raw, "\\s*///\\s*|\\s*;\\s*|\\s*,\\s*", 2)[, 1]
      ) %>%
      dplyr::filter(!is.na(gene_symbol), gene_symbol != "")

    # Collapse probes by gene symbol (median of numeric expression columns)
    expr_gene <- gse_expr %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Probe") %>%
      dplyr::left_join(probe_map, by = "Probe") %>%
      dplyr::filter(!is.na(gene_symbol)) %>%
      dplyr::select(gene_symbol, dplyr::where(is.numeric)) %>%
      dplyr::group_by(gene_symbol) %>%
      dplyr::summarise(dplyr::across(dplyr::where(is.numeric), median, na.rm = TRUE), .groups = "drop") %>%
      tibble::column_to_rownames("gene_symbol") %>%
      as.matrix()

    cat(sprintf("Collapsed %d probes to %d genes\n", nrow(gse_expr), nrow(expr_gene)))

    # Perform DE analysis using limma
    if (nrow(status_counts) == 2 && all(status_counts$n >= 2)) {
      keep <- gse_pheno$menopause_status %in% c("Pre", "Post")
      gse_pheno2 <- gse_pheno[keep, , drop = FALSE]
      expr_gene2 <- expr_gene[, rownames(gse_pheno2), drop = FALSE]

      # Create design matrix
      design <- stats::model.matrix(~ 0 + menopause_status, data = gse_pheno2)
      colnames(design) <- sub("^menopause_status", "", colnames(design))

      # Fit linear model
      fit <- limma::lmFit(expr_gene2, design)

      # Make contrast: Post vs Pre
      contrast.matrix <- limma::makeContrasts(Post - Pre, levels = design)
      fit2 <- limma::contrasts.fit(fit, contrast.matrix)
      fit2 <- limma::eBayes(fit2)

      # Extract results
      de_results <- limma::topTable(fit2, number = Inf, adjust.method = "BH")

      # Save results
      readr::write_tsv(
        de_results %>%
          rownames_to_column("gene_symbol") %>%
          mutate(sig = ifelse(adj.P.Val < 0.05, "Yes", "No")),
        file.path(tab_dir, "GSE86244_DE_results_all.tsv")
      )

      # Significant genes
      de_sig <- de_results %>%
        rownames_to_column("gene_symbol") %>%
        filter(adj.P.Val < 0.05)

      cat(sprintf("\nGSE86244 DE analysis: %d significant genes (padj < 0.05)\n", nrow(de_sig)))
      cat("Top 20 genes:\n")
      print(head(de_sig %>% arrange(adj.P.Val), 20))

      # Save significant genes
      readr::write_tsv(
        de_sig,
        file.path(tab_dir, "GSE86244_significant_genes.tsv")
      )

      # Save up and down separately
      de_up <- de_sig %>% filter(logFC > 0)
      de_down <- de_sig %>% filter(logFC < 0)

      readr::write_tsv(
        de_up %>% arrange(adj.P.Val),
        file.path(tab_dir, "GSE86244_upregulated_post_vs_pre.tsv")
      )

      readr::write_tsv(
        de_down %>% arrange(adj.P.Val),
        file.path(tab_dir, "GSE86244_downregulated_post_vs_pre.tsv")
      )

      cat(sprintf("  Upregulated in postmenopausal: %d genes\n", nrow(de_up)))
      cat(sprintf("  Downregulated in postmenopausal: %d genes\n", nrow(de_down)))

      ## -----------------------------------------------------------------------------
      ## Compare with GTEx signature
      ## -----------------------------------------------------------------------------
      cat("\n===========================================\n")
      cat("VALIDATION: GTEx Signature vs GSE86244\n")
      cat("===========================================\n\n")

      # Calculate overlap
      overlap_up <- intersect(signature_genes, de_up$gene_symbol)
      overlap_down <- intersect(signature_genes, de_down$gene_symbol)
      overlap_all <- intersect(signature_genes, de_sig$gene_symbol)

      cat(sprintf("\nGTEx Signature size: %d genes\n", length(signature_genes)))
      cat(sprintf("GSE86244 Significant: %d genes\n", nrow(de_sig)))
      cat(sprintf("Overlap: %d genes (%.1f%%)\n",
                  length(overlap_all),
                  100 * length(overlap_all) / length(signature_genes)))

      if (length(overlap_up) > 0) {
        cat(sprintf("  Upregulated in post (GSE86244): %d genes\n", length(overlap_up)))
        cat("    Genes: ", paste(overlap_up, collapse = ", "), "\n")
      }

      if (length(overlap_down) > 0) {
        cat(sprintf("  Downregulated in post (GSE86244): %d genes\n", length(overlap_down)))
        cat("    Genes: ", paste(overlap_down, collapse = ", "), "\n")
      }

      ## Fisher's exact test for enrichment
      # Create contingency table
      #                   InGSE86244Sig  NotInGSE86244Sig
      # InGTExSig              a                b
      # NotInGTExSig           c                d

      a <- length(overlap_all)
      b <- length(signature_genes) - a
      c <- nrow(de_sig) - a
      d <- nrow(rownames(expr_gene)) - nrow(de_sig) - b

      contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
      rownames(contingency_table) <- c("InGTExSig", "NotInGTExSig")
      colnames(contingency_table) <- c("InGSE86244Sig", "NotInGSE86244Sig")

      cat("\nContingency table:\n")
      print(contingency_table)

      fisher_test <- fisher.test(contingency_table)

      cat(sprintf("\nFisher's exact test:\n"))
      cat(sprintf("  Odds ratio: %.2f\n", fisher_test$estimate))
      cat(sprintf("  P-value: %g\n", fisher_test$p.value))

      if (fisher_test$p.value < 0.05) {
        cat("  *** SIGNIFICANT ENRICHMENT ***\n")
      } else {
        cat("  (Not significant)\n")
      }

      ## -----------------------------------------------------------------------------
      ## Create comparison visualizations
      ## -----------------------------------------------------------------------------
      cat("\nCreating validation visualizations...\n")

      # Venn diagram data
      venn_df <- tibble(
        gene_set = c(rep("GTEx_Sig", length(signature_genes)),
                     rep("GSE86244_Sig", nrow(de_sig))),
        gene = c(signature_genes, de_sig$gene_symbol)
      ) %>%
        group_by(gene) %>%
        mutate(n_sets = n()) %>%
        ungroup() %>%
        mutate(category = case_when(
          n_sets == 2 ~ "Overlap",
          gene_set == "GTEx_Sig" ~ "GTEx only",
          gene_set == "GSE86244_Sig" ~ "GSE86244 only"
        ))

      # Bar plot
      p_bar <- ggplot(venn_df %>% distinct(gene, category), aes(x = category, fill = category)) +
        geom_bar() +
        geom_text(aes(label = after_stat(count)), stat = "count", vjust = -0.5) +
        labs(
          title = "Overlap Between GTEx and GSE86244 Signatures",
          subtitle = sprintf("GTEx: %d genes | GSE86244: %d genes | Overlap: %d genes",
                            length(signature_genes), nrow(de_sig), length(overlap_all)),
          x = "",
          y = "Number of Genes"
        ) +
        theme_bw() +
        theme(legend.position = "none")

      ggsave(file.path(fig_dir, "signature_overlap_barplot.png"), p_bar, width = 8, height = 6, dpi = 300)

      # Overlap genes detailed table
      if (length(overlap_all) > 0) {
        overlap_details <- de_sig %>%
          filter(gene_symbol %in% overlap_all) %>%
          select(gene_symbol, logFC, AveExpr, t, P.Value, adj.P.Val) %>%
          arrange(adj.P.Val)

        readr::write_tsv(
          overlap_details,
          file.path(tab_dir, "GTEx_GSE86244_overlap_genes.tsv")
        )

        cat(sprintf("\nSaved overlap details: %d genes\n", nrow(overlap_details)))
      }

      ## -----------------------------------------------------------------------------
      ## Create validation summary document
      ## -----------------------------------------------------------------------------
      cat("\nCreating validation summary document...\n")

      # Pre-compute conditional strings
      up_section <- if (length(overlap_up) > 0) {
        paste0("### Upregulated in Postmenopausal (GSE86244)\n\n",
               "Genes from GTEx signature also up in GSE86244:\n",
               paste0("- ", overlap_up, collapse = "\n"),
               "\n\n")
      } else {
        ""
      }

      down_section <- if (length(overlap_down) > 0) {
        paste0("### Downregulated in Postmenopausal (GSE86244)\n\n",
               "Genes from GTEx signature also down in GSE86244:\n",
               paste0("- ", overlap_down, collapse = "\n"),
               "\n\n")
      } else {
        ""
      }

      sig_result <- if (fisher_test$p.value < 0.05) {
        "- Result: **SIGNIFICANT** enrichment\n\n"
      } else {
        "- Result: Not significant\n\n"
      }

      interpretation <- if (fisher_test$p.value < 0.05 && length(overlap_all) > 5) {
        paste0("The GTEx signature shows significant overlap with genes differentially\n",
               "expressed in GSE86244, which has actual menopause status. This validates\n",
               "that our age-proxy approach successfully identified menopause-related genes.\n\n")
      } else if (length(overlap_all) > 0) {
        paste0(sprintf("While %d genes overlap (%.1f%%), the enrichment is not statistically\n",
                       length(overlap_all), 100 * length(overlap_all) / length(signature_genes)),
               "significant. This could be due to:\n",
               sprintf("- Small sample size in GSE86244 (n=%d)\n", ncol(gse_expr)),
               "- Platform differences (RNA-seq vs microarray)\n",
               "- Tissue-specific differences\n",
               "- Age proxy imperfection in GTEx\n\n")
      } else {
        paste0("No overlap detected. Possible reasons:\n",
               "- Small sample size in GSE86244\n",
               "- Different platforms (RNA-seq vs microarray)\n",
               "- GTEx age proxy may not capture true menopause signal\n",
               "- Batch effects or other confounders\n\n")
      }

      recommendation <- if (fisher_test$p.value < 0.05) {
        "- GTEx signature is validated by GSE86244\n"
      } else {
        "- Consider expanding signature criteria\n"
      }

      groups_str <- paste(status_counts$menopause_status, status_counts$n, sep = " = ", collapse = ", ")

      validation_summary <- paste0(
        "# GSE86244 Validation Summary\n\n",
        "## Purpose\n\n",
        "Validate the GTEx-derived menopause signature against GSE86244,\n",
        "which contains actual pre/post-menopausal status (not age proxy).\n\n",

        "## Datasets Compared\n\n",
        "### GTEx Signature (Approach 1)\n",
        "- Source: GTEx v10 subcutaneous adipose tissue\n",
        "- Method: Differential expression by age proxy (30-45 vs 55-70)\n",
        "- Genes: ", length(signature_genes), "\n",
        "- Criteria: |log2FC| > 1, padj < 0.05\n\n",

        "### GSE86244 (Shan et al.)\n",
        "- Source: Subcutaneous adipose tissue\n",
        "- Samples: ", ncol(gse_expr), " women\n",
        "- Groups: ", groups_str, "\n",
        "- Platform: Affymetrix HG-U133 Plus 2.0\n\n",

        "## Validation Results\n\n",
        "### Overlap Analysis\n\n",
        "- GTEx signature genes: ", length(signature_genes), "\n",
        "- GSE86244 significant genes: ", nrow(de_sig), "\n",
        "- Overlapping genes: ", length(overlap_all), "\n",
        sprintf("- Overlap percentage: %.1f%%\n\n", 100 * length(overlap_all) / length(signature_genes)),

        up_section,
        down_section,

        "### Statistical Test\n\n",
        "Fisher's exact test for enrichment:\n",
        sprintf("- Odds ratio: %.2f\n", fisher_test$estimate),
        sprintf("- P-value: %g\n", fisher_test$p.value),
        sig_result,

        "## Interpretation\n\n",
        interpretation,

        "## Recommendations\n\n",
        recommendation,
        "- Validate with additional menopause datasets\n",
        "- Test signature on other tissues (blood, skin)\n",
        "- Consider estrogen-responsive gene filter\n\n",

        "---\n",
        "Generated: ", Sys.time(), "\n"
      )

      writeLines(validation_summary, file.path(tab_dir, "GSE86244_validation_summary.md"))

      cat(sprintf("\nSaved validation summary: %s\n",
                  file.path(tab_dir, "GSE86244_validation_summary.md")))

    } else {
      cat("\nCannot perform DE analysis: insufficient samples per group\n")
      cat("Need at least 2 samples per group for limma\n")
    }

  }, error = function(e) {
    cat("Platform analysis failed:\n")
    cat(sprintf("  Error: %s\n", e$message))
  })

  }  # Close the `else { gse_data_valid <- TRUE` block

} else {
  ## -----------------------------------------------------------------------------
  ## Fallback: Simulated validation (for demonstration)
  ## -----------------------------------------------------------------------------
  cat("\n===========================================\n")
  cat("FALLBACK MODE: Simulated Validation\n")
  cat("===========================================\n\n")

  cat("GSE86244 data not available. Using simulated validation.\n\n")

  # In production, this would be replaced with actual GSE86244 analysis
  # For now, create template output showing structure

  cat("Expected validation steps:\n")
  cat("1. Download GSE86244 from GEO\n")
  cat("2. Extract menopause status from sample annotations\n")
  cat("3. Perform differential expression (limma)\n")
  cat("4. Calculate overlap with GTEx signature\n")
  cat("5. Test enrichment with Fisher's exact test\n\n")

  # Create template files
  template_summary <- paste0(
    "# GSE86244 Validation - TEMPLATE\n\n",
    "## Status\n\n",
    "**DATA REQUIRED**: Manual download of GSE86244 from GEO\n\n",
    "## Download Instructions\n\n",
    "1. Download from: https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE86244&format=file\n",
    "2. Or visit: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86244\n",
    "3. Save the series matrix file to: ", gse_dir, "/GSE86244_series_matrix.txt.gz\n",
    "4. Re-run this script\n\n",

    "## Expected Output Structure\n\n",
    "After download, this script will generate:\n",
    "- GSE86244_DE_results_all.tsv\n",
    "- GSE86244_significant_genes.tsv\n",
    "- GTEx_GSE86244_overlap_genes.tsv\n",
    "- GSE86244_validation_summary.md\n\n",

    "## GTEx Signature for Reference\n\n",
    "Genes to validate: ", length(signature_genes), "\n\n",
    paste(signature_genes, collapse = "\n"), "\n\n",

    "---\n",
    "Generated: ", Sys.time(), "\n"
  )

  writeLines(template_summary, file.path(tab_dir, "GSE86244_validation_template.md"))

  cat("Saved template: ", file.path(tab_dir, "GSE86244_validation_template.md"), "\n")
}

## -----------------------------------------------------------------------------
## Summary
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("GSE86244 VALIDATION SUMMARY\n")
cat("===========================================\n\n")

cat(sprintf("GTEx signature size: %d genes\n", length(signature_genes)))

if (file.exists(file.path(tab_dir, "GSE86244_DE_results_all.tsv"))) {
  de_results <- readr::read_tsv(file.path(tab_dir, "GSE86244_DE_results_all.tsv"), show_col_types = FALSE)
  cat(sprintf("GSE86244 DE results: %d genes tested\n", nrow(de_results)))
}

if (file.exists(file.path(tab_dir, "GTEx_GSE86244_overlap_genes.tsv"))) {
  overlap <- readr::read_tsv(file.path(tab_dir, "GTEx_GSE86244_overlap_genes.tsv"), show_col_types = FALSE)
  cat(sprintf("Overlap genes: %d\n", nrow(overlap)))
}

cat("\nOutput files:\n")
cat(sprintf("  Tables: %s\n", tab_dir))
cat(sprintf("  Figures: %s\n", fig_dir))

cat("\n===========================================\n\n")
