################################################################################
## Next steps (2026-02-01): External dataset search + initial triage
##
## Goals:
## - Programmatically query GEO DataSets (NCBI E-utilities) for menopause-related
##   high-throughput expression studies, prioritizing adipose tissue.
## - Query ArrayExpress API for menopause-related experiments.
## - Write a single, curated candidate table with enough metadata to decide what
##   to download/analyze next.
##
## NOTE:
## - Many human menopause studies do not deposit full matrices; this step
##   captures what can be discovered programmatically + a small manual seed list.
##
## Outputs:
## - GTEx_v10_AT_analysis_out/tables/external_dataset_search_2026-02-01/
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
  xml2, jsonlite
)

source("R/utils.R")

out_dir <- "GTEx_v10_AT_analysis_out"
tab_dir <- file.path(out_dir, "tables", "external_dataset_search_2026-02-01")
ensure_dirs(tab_dir)

## -----------------------------------------------------------------------------
## GEO E-utilities helpers (db=gds)
## -----------------------------------------------------------------------------
base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

esearch_gds <- function(term, retmax = 200) {
  url <- paste0(
    base,
    "/esearch.fcgi?db=gds&term=",
    URLencode(term, reserved = TRUE),
    "&retmax=",
    retmax
  )
  doc <- NULL
  for (attempt in 1:6) {
    doc <- tryCatch(xml2::read_xml(url), error = function(e) e)
    if (!inherits(doc, "error")) break
    wait <- min(10, 2^(attempt - 1))
    message("esearch retry (attempt ", attempt, ") after ", wait, "s: ", conditionMessage(doc))
    Sys.sleep(wait)
  }
  if (inherits(doc, "error")) {
    stop("Failed to fetch esearch after retries: ", url)
  }
  ids <- xml2::xml_text(xml2::xml_find_all(doc, ".//IdList/Id"))
  ids
}

esummary_gds <- function(ids) {
  if (!length(ids)) return(tibble::tibble())
  # NCBI recommends batching IDs
  chunks <- split(ids, ceiling(seq_along(ids) / 100))
  out <- list()
  for (ch in chunks) {
    url <- paste0(base, "/esummary.fcgi?db=gds&id=", paste(ch, collapse = ","))
    doc <- NULL
    for (attempt in 1:6) {
      doc <- tryCatch(xml2::read_xml(url), error = function(e) e)
      if (!inherits(doc, "error")) break
      wait <- min(10, 2^(attempt - 1))
      message("esummary retry (attempt ", attempt, ") after ", wait, "s: ", conditionMessage(doc))
      Sys.sleep(wait)
    }
    if (inherits(doc, "error")) {
      stop("Failed to fetch esummary after retries: ", url)
    }
    docs <- xml2::xml_find_all(doc, ".//DocSum")
    for (d in docs) {
      items <- xml2::xml_find_all(d, ".//Item")
      kv <- setNames(xml2::xml_text(items), xml2::xml_attr(items, "Name"))
      out[[length(out) + 1]] <- kv
    }
    Sys.sleep(0.34)
  }
  # normalize to a tibble with common columns
  keys <- unique(unlist(lapply(out, names)))
  mat <- do.call(rbind, lapply(out, function(kv) {
    v <- rep(NA_character_, length(keys))
    names(v) <- keys
    v[names(kv)] <- kv
    v
  }))
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  tibble::as_tibble(df)
}

## -----------------------------------------------------------------------------
## ArrayExpress API helper (JSON v3; experiments endpoint)
## -----------------------------------------------------------------------------
arrayexpress_search <- function(keywords, page_size = 100) {
  # API docs: https://www.ebi.ac.uk/arrayexpress/help/programmatic_access.html
  # v3 JSON endpoint:
  url <- paste0(
    "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?keywords=",
    URLencode(keywords, reserved = TRUE),
    "&pagesize=",
    page_size
  )
  js <- tryCatch(jsonlite::fromJSON(url), error = function(e) NULL)
  empty <- tibble::tibble(
    repository = character(0),
    accession = character(0),
    title = character(0),
    organism = character(0),
    experimental_factor = character(0),
    last_updated = character(0)
  )
  if (is.null(js)) return(empty)
  exps <- js$experiments$experiment
  if (is.null(exps) || !length(exps)) return(empty)
  tibble::tibble(
    repository = "ArrayExpress",
    accession = exps$accession,
    title = exps$name,
    organism = exps$organism,
    experimental_factor = exps$experimentalvariable,
    last_updated = exps$lastupdatedate
  )
}

## -----------------------------------------------------------------------------
## Queries
## -----------------------------------------------------------------------------
queries <- c(
  "\"postmenopausal\" adipose expression profiling",
  "\"premenopausal\" \"postmenopausal\" adipose",
  "premenopausal postmenopausal adipose microarray",
  "menopause adipose tissue expression profiling",
  "menopause adipose microarray",
  "menopause adipose RNA-seq",
  "ovariectomy adipose gene expression",
  "\"postmenopausal\" \"subcutaneous adipose\"",
  "\"postmenopausal\" \"visceral adipose\""
)

cat("Querying GEO DataSets (gds) with E-utilities...\n")
geo_rows <- list()
for (q in queries) {
  ids <- esearch_gds(q, retmax = 200)
  if (!length(ids)) next
  sm <- esummary_gds(ids)
  if (!nrow(sm)) next
  sm <- sm %>%
    dplyr::mutate(
      repository = "GEO (GDS)",
      query = q,
      accession = .data$Accession,
      title = .data$title,
      summary = .data$summary,
      entryType = .data$entryType,
      gse = dplyr::na_if(.data$GSE, ""),
      gpl = dplyr::na_if(.data$GPL, ""),
      organism = .data$taxon,
      n_samples = suppressWarnings(as.integer(.data$n_samples))
    ) %>%
    dplyr::mutate(
      gse = dplyr::if_else(!is.na(gse) & entryType == "GSE", paste0("GSE", gse), NA_character_),
      gpl = dplyr::if_else(!is.na(gpl) & nzchar(gpl), paste0("GPL", gpl), NA_character_)
    ) %>%
    dplyr::filter(entryType == "GSE") %>%
    dplyr::select(repository, query, entryType, accession, gse, gpl, organism, n_samples, title, summary)
  geo_rows[[length(geo_rows) + 1]] <- sm
}
geo_tbl <- dplyr::bind_rows(geo_rows) %>%
  dplyr::distinct(repository, accession, gse, .keep_all = TRUE)

cat("Querying ArrayExpress...\n")
ae_tbl <- dplyr::bind_rows(
  arrayexpress_search("menopause adipose"),
  arrayexpress_search("postmenopausal adipose"),
  arrayexpress_search("menopause"),
  arrayexpress_search("ovariectomy")
) %>%
  dplyr::distinct(accession, .keep_all = TRUE)

## Manual seed list (known from prior project notes + common menopause adipose papers)
manual <- tibble::tibble(
  repository = c("GEO", "GEO", "GEO"),
  accession = c("GSE86244", "GSE44000", NA_character_),
  title = c(
    "Human adipose stromal/stem cells RNA-seq (ASC); age proxy pre<45 vs post>55 (used in this repo)",
    "Adipose tissue dataset referenced in postmenopausal obesity/NETs paper (needs verification of cohort/tissue)",
    "Menopause SAT+VAT microarray in morbidly obese women (PMID:21358552; check if deposited)"
  ),
  notes = c(
    "Already analyzed (supplementary counts); ASC context mismatch vs bulk adipose.",
    "Candidate for mining; menopause status unclear without full metadata check.",
    "May not have public deposition; could rely on supplemental tables or author contact."
  )
)

## Combine into a single candidate table (triage-friendly)
candidates <- geo_tbl %>%
  dplyr::mutate(
    accession_primary = dplyr::coalesce(gse, accession),
    hits_adipose = stringr::str_detect(tolower(paste(title, summary)), "adipose|fat|subcutaneous|omental|visceral"),
    hits_menopause = stringr::str_detect(tolower(paste(title, summary)), "menopaus|postmenopaus|premenopaus|ovariectom"),
    priority = dplyr::case_when(
      hits_adipose & hits_menopause ~ "high",
      hits_menopause ~ "medium",
      TRUE ~ "low"
    )
  ) %>%
  dplyr::select(repository, accession_primary, gse, gpl, organism, n_samples, priority, title, summary, query) %>%
  dplyr::arrange(dplyr::desc(priority), dplyr::desc(n_samples))

ae_candidates <- ae_tbl %>%
  dplyr::mutate(
    accession_primary = accession,
    priority = dplyr::case_when(
      stringr::str_detect(tolower(title), "adipose|fat") & stringr::str_detect(tolower(title), "menopaus|postmenopaus|ovariectom") ~ "high",
      stringr::str_detect(tolower(title), "menopaus|postmenopaus|ovariectom") ~ "medium",
      TRUE ~ "low"
    )
  ) %>%
  dplyr::select(repository, accession_primary, accession, organism, priority, title, experimental_factor, last_updated)

readr::write_tsv(candidates, file.path(tab_dir, "geo_gds_menopause_search_candidates.tsv"))
readr::write_tsv(ae_candidates, file.path(tab_dir, "arrayexpress_menopause_search_candidates.tsv"))
readr::write_tsv(manual, file.path(tab_dir, "manual_seed_candidates.tsv"))

## One merged view for discussion
merged <- dplyr::bind_rows(
  candidates %>%
    dplyr::transmute(
      repository,
      accession = accession_primary,
      organism,
      n_samples,
      priority,
      title,
      source = "GEO(gds)",
      notes = query
    ),
  ae_candidates %>%
    dplyr::transmute(
      repository,
      accession = accession_primary,
      organism,
      n_samples = NA_integer_,
      priority,
      title,
      source = "ArrayExpress",
      notes = experimental_factor
    ),
  manual %>%
    dplyr::transmute(
      repository,
      accession = accession,
      organism = NA_character_,
      n_samples = NA_integer_,
      priority = "manual",
      title,
      source = "Manual",
      notes
    )
) %>%
  dplyr::distinct(repository, accession, title, .keep_all = TRUE) %>%
  dplyr::arrange(dplyr::desc(priority), accession)

readr::write_tsv(merged, file.path(tab_dir, "external_dataset_candidates_merged.tsv"))

md <- c(
  "# External dataset search (2026-02-01)",
  "",
  "This folder contains a programmatic search over GEO DataSets (via NCBI E-utilities; db=gds) and ArrayExpress experiments (JSON API), plus a small manual seed list.",
  "",
  "## Outputs",
  "- `geo_gds_menopause_search_candidates.tsv`: raw GEO(gds) hits with extracted GSE/GPL when possible",
  "- `arrayexpress_menopause_search_candidates.tsv`: raw ArrayExpress hits from keyword searches",
  "- `manual_seed_candidates.tsv`: manually tracked candidates (including datasets already used)",
  "- `external_dataset_candidates_merged.tsv`: merged triage table for team discussion",
  "",
  "## Next triage steps (human-in-the-loop)",
  "1) Filter merged candidates to human + adipose + menopause terms (priority=high/manual).",
  "2) For top candidates, open the GEO/ArrayExpress record and verify: tissue, menopause status labeling, and availability of expression matrix/raw data.",
  "3) Select 1–3 feasible datasets and implement a harmonized analysis in separate folders under `GTEx_v10_AT_analysis_out/tables/external_validation_<ACCESSION>/`.",
  "",
  "Generated:",
  as.character(Sys.time())
)
writeLines(md, file.path(tab_dir, "README_external_dataset_search.md"))

cat("Done. Wrote candidate tables to: ", tab_dir, "\n", sep = "")
