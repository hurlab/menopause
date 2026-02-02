################################################################################
## Pathway Enrichment Analysis for Menopause Signature
##
## Analyzes GO, KEGG, and Reactome pathways enriched in the 81-gene signature
## Uses clusterProfiler for enrichment analysis
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
  ggplot2, patchwork, clusterProfiler, org.Hs.eg.db, enrichplot, ReactomePA
)

## -----------------------------------------------------------------------------
## Source utility functions
## -----------------------------------------------------------------------------
source("R/utils.R")

## -----------------------------------------------------------------------------
## Setup
## -----------------------------------------------------------------------------
out_dir  <- "GTEx_v10_AT_analysis_out"
fig_dir  <- file.path(out_dir, "figs", "enrichment_analysis")
tab_dir  <- file.path(out_dir, "tables", "enrichment_analysis")
ensure_dirs(c(out_dir, fig_dir, tab_dir))

cat("\n===========================================\n")
cat("PATHWAY ENRICHMENT ANALYSIS\n")
cat("===========================================\n\n")

## -----------------------------------------------------------------------------
## Load menopause signature genes
## -----------------------------------------------------------------------------
cat("Loading menopause signature genes...\n")

sig_file <- file.path("GTEx_v10_AT_analysis_out/tables/approach1_improved_age_proxy/menopause_signature_gene_list.txt")
if (!file.exists(sig_file)) {
  stop("Signature file not found. Run Approach 1 first.")
}

signature_genes <- read_lines(sig_file)
cat(sprintf("Loaded %d signature genes\n", length(signature_genes)))

## Also load DE results for context
de_file <- file.path(out_dir, "tables", "approach1_improved_age_proxy", "DE_menopause_proxy_improved.tsv")
if (file.exists(de_file)) {
  de_results <- readr::read_tsv(de_file, show_col_types = FALSE)
  cat(sprintf("Loaded %d DE results for context\n", nrow(de_results)))
}

## -----------------------------------------------------------------------------
## Prepare gene list for clusterProfiler (needs Ensembll IDs)
## -----------------------------------------------------------------------------
cat("\nConverting gene symbols to Entrez IDs...\n")

## Get mapping from org.Hs.eg.db
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = signature_genes,
  columns = c("ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "SYMBOL"
)

## Filter to genes found in database
gene_map_filtered <- gene_map %>%
  filter(!is.na(ENTREZID)) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

entrez_ids <- gene_map_filtered$ENTREZID
cat(sprintf("Mapped %d genes to Entrez IDs\n", length(entrez_ids)))

## -----------------------------------------------------------------------------
## GO Biological Process enrichment
## -----------------------------------------------------------------------------
cat("\nRunning GO Biological Process enrichment...\n")

tryCatch({
  go_bp <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  # Convert to tibble
  go_bp_df <- as.data.frame(go_bp) %>%
    mutate(ontology = "GO-BP")

  # Save results
  readr::write_tsv(
    go_bp_df %>%
      dplyr::select(ID, Description, GeneRatio, pvalue, p.adjust, qvalue, Count, ontology),
    file.path(tab_dir, "GO_BP_enrichment.tsv")
  )

  cat(sprintf("GO BP enrichment: %d enriched terms\n", nrow(go_bp_df)))

  # Top 10 terms
  cat("\nTop 10 GO Biological Processes:\n")
  print(head(go_bp_df %>% arrange(p.adjust), 10))

}, error = function(e) {
  cat("GO BP enrichment failed:\n")
  cat(sprintf("  Error: %s\n", e$message))
})

## -----------------------------------------------------------------------------
## GO Cellular Component enrichment
## -----------------------------------------------------------------------------
cat("\nRunning GO Cellular Component enrichment...\n")

tryCatch({
  go_cc <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  go_cc_df <- as.data.frame(go_cc) %>%
    mutate(ontology = "GO-CC")

  readr::write_tsv(
    go_cc_df %>%
      dplyr::select(ID, Description, GeneRatio, pvalue, p.adjust, qvalue, Count, ontology),
    file.path(tab_dir, "GO_CC_enrichment.tsv")
  )

  cat(sprintf("GO CC enrichment: %d enriched terms\n", nrow(go_cc_df)))

}, error = function(e) {
  cat("GO CC enrichment failed:\n")
  cat(sprintf("  Error: %s\n", e$message))
})

## -----------------------------------------------------------------------------
## KEGG pathway enrichment
## -----------------------------------------------------------------------------
cat("\nRunning KEGG pathway enrichment...\n")

tryCatch({
  kegg <- enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2
  )

  kegg_df <- as.data.frame(kegg)

  readr::write_tsv(
    kegg_df %>%
      dplyr::select(ID, Description, GeneRatio, pvalue, p.adjust, qvalue, Count),
    file.path(tab_dir, "KEGG_enrichment.tsv")
  )

  cat(sprintf("KEGG enrichment: %d enriched pathways\n", nrow(kegg_df)))

  cat("\nTop 10 KEGG pathways:\n")
  print(head(kegg_df %>% arrange(p.adjust), 10))

}, error = function(e) {
  cat("KEGG enrichment failed:\n")
  cat(sprintf("  Error: %s\n", e$message))
})

## -----------------------------------------------------------------------------
## Reactome pathway enrichment
## -----------------------------------------------------------------------------
cat("\nRunning Reactome pathway enrichment...\n")

tryCatch({
  reactome <- enrichPathway(
    gene = entrez_ids,
    organism = "human",
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.2,
    readable = TRUE
  )

  reactome_df <- as.data.frame(reactome)

  readr::write_tsv(
    reactome_df %>%
      dplyr::select(dplyr::any_of(c(
        "ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "Count", "geneID"
      ))),
    file.path(tab_dir, "Reactome_enrichment.tsv")
  )

  cat(sprintf("Reactome enrichment: %d enriched pathways\n", nrow(reactome_df)))

}, error = function(e) {
  cat("Reactome enrichment failed:\n")
  cat(sprintf("  Error: %s\n", e$message))
})

## -----------------------------------------------------------------------------
## Create summary document
## -----------------------------------------------------------------------------
cat("\nCreating enrichment summary document...\n")

# Get top terms from each category
go_bp_top <- tryCatch({
  head(go_bp_df %>% arrange(p.adjust), 10)
}, error = function(e) NULL)

go_cc_top <- tryCatch({
  head(go_cc_df %>% arrange(p.adjust), 10)
}, error = function(e) NULL)

kegg_top <- tryCatch({
  head(kegg_df %>% arrange(p.adjust), 10)
}, error = function(e) NULL)

reactome_top <- tryCatch({
  head(reactome_df %>% arrange(p.adjust), 10)
}, error = function(e) NULL)

summary_doc <- paste0(
  "# Pathway Enrichment Analysis - Menopause Signature\n\n",
  "## Signature Overview\n\n",
  "- Total genes: ", length(signature_genes), "\n",
  "- Mapped to Entrez: ", length(entrez_ids), "\n\n",
  "## Enriched Pathways\n\n",
  "Note: Full enrichment results available in separate TSV files.\n\n",
  "## Key Findings\n\n",
  "The enrichment analysis was performed using clusterProfiler.\n",
  "Results include GO Biological Process, GO Cellular Component,\n",
  "KEGG pathways, and Reactome pathways.\n\n",
  "## Interpretation\n\n",
  "The enriched pathways provide insight into biological processes\n",
  "associated with menopause-related transcriptional changes in adipose tissue.\n\n",
  "---\n",
  "Generated: ", Sys.time(), "\n"
)

writeLines(summary_doc, file.path(tab_dir, "enrichment_analysis_summary.md"))

cat(sprintf("Saved summary to: %s\n", file.path(tab_dir, "enrichment_analysis_summary.md")))

## -----------------------------------------------------------------------------
## Create visualizations
## -----------------------------------------------------------------------------
cat("\nCreating visualizations...\n")

if (file.exists(file.path(tab_dir, "GO_BP_enrichment.tsv"))) {
  go_data <- readr::read_tsv(file.path(tab_dir, "GO_BP_enrichment.tsv"), show_col_types = FALSE)

  if (nrow(go_data) > 0) {
    # Get top 20 terms
    go_top_20 <- go_data %>%
      arrange(p.adjust) %>%
      head(20)

    # Create dot plot
    p_dot <- enrichplot::dotplot(go_bp, showCategory = TRUE, title = "GO Biological Processes") +
      ggplot2::theme(legend.position = "right")

    ggsave(file.path(fig_dir, "GO_BP_dotplot.png"), p_dot, width = 12, height = 10, dpi = 300)

    cat("Saved GO BP dotplot\n")
  }
}

cat(sprintf("Saved figures to: %s\n", fig_dir))

## -----------------------------------------------------------------------------
## Summary
## -----------------------------------------------------------------------------
cat("\n===========================================\n")
cat("ENRICHMENT ANALYSIS SUMMARY\n")
cat("===========================================\n\n")

cat(sprintf("Signature size: %d genes\n", length(signature_genes)))

if (file.exists(file.path(tab_dir, "GO_BP_enrichment.tsv"))) {
  n_go_bp <- nrow(readr::read_tsv(file.path(tab_dir, "GO_BP_enrichment.tsv"), show_col_types = FALSE))
  cat(sprintf("GO BP enriched terms: %d\n", n_go_bp))
}

if (file.exists(file.path(tab_dir, "KEGG_enrichment.tsv"))) {
  n_kegg <- nrow(readr::read_tsv(file.path(tab_dir, "KEGG_enrichment.tsv"), show_col_types = FALSE))
  cat(sprintf("KEGG enriched pathways: %d\n", n_kegg))
}

cat("\nOutput files:\n")
cat(sprintf("  GO BP: %s\n", file.path(tab_dir, "GO_BP_enrichment.tsv")))
cat(sprintf("  GO CC: %s\n", file.path(tab_dir, "GO_CC_enrichment.tsv")))
cat(sprintf("  KEGG: %s\n", file.path(tab_dir, "KEGG_enrichment.tsv")))
cat(sprintf("  Reactome: %s\n", file.path(tab_dir, "Reactome_enrichment.tsv")))
cat(sprintf("  Summary: %s\n", file.path(tab_dir, "enrichment_analysis_summary.md")))

cat("===========================================\n\n")
