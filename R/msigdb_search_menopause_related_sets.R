################################################################################
## MSigDB discovery: menopause / estrogen / related gene sets
##
## Uses the `msigdbr` package to list MSigDB gene sets whose names match a set
## of menopause-related keywords (and estrogen-response proxies).
##
## Output:
## - GTEx_v10_AT_analysis_out/tables/msigdb_search_menopause_related/
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
  dplyr, stringr, readr, tibble,
  msigdbr
)

source("R/utils.R")

out_dir <- "GTEx_v10_AT_analysis_out"
tab_dir <- file.path(out_dir, "tables", "msigdb_search_menopause_related")
ensure_dirs(tab_dir)

keywords <- c(
  "MENOPAUS", "POSTMENOPAUS", "PREMENOPAUS", "PERIMENOPAUS",
  "OVARIECTOM", "OOPHORECTOM", "OVX",
  "ESTROGEN", "ESTRADIOL", "ESR1", "ESR2", "ER_ALPHA", "ER_BETA", "PROGESTERONE"
)
pattern <- paste0("(", paste(keywords, collapse = "|"), ")")

msig <- msigdbr::msigdbr(species = "Homo sapiens")

hit_sets <- msig %>%
  dplyr::distinct(gs_name, gs_collection, gs_subcollection, gs_collection_name) %>%
  dplyr::filter(stringr::str_detect(.data$gs_name, stringr::regex(pattern, ignore_case = TRUE)))

hit_summary <- msig %>%
  dplyr::semi_join(hit_sets, by = c("gs_name", "gs_collection", "gs_subcollection", "gs_collection_name")) %>%
  dplyr::count(gs_name, gs_collection, gs_subcollection, gs_collection_name, name = "n_genes") %>%
  dplyr::arrange(desc(.data$n_genes), .data$gs_name)

readr::write_tsv(hit_summary, file.path(tab_dir, "msigdb_keyword_hit_sets.tsv"))
writeLines(hit_summary$gs_name, file.path(tab_dir, "msigdb_keyword_hit_set_names.txt"))

cat("\n===========================================\n")
cat("MSIGDB KEYWORD SEARCH SUMMARY\n")
cat("===========================================\n\n")
cat("Keywords: ", paste(keywords, collapse = ", "), "\n", sep = "")
cat(sprintf("Total gene sets in msigdbr(Hs): %d\n", length(unique(msig$gs_name))))
cat(sprintf("Keyword-hit gene sets: %d\n\n", nrow(hit_summary)))
cat("Wrote:\n")
cat("  - ", file.path(tab_dir, "msigdb_keyword_hit_sets.tsv"), "\n", sep = "")
cat("  - ", file.path(tab_dir, "msigdb_keyword_hit_set_names.txt"), "\n", sep = "")
