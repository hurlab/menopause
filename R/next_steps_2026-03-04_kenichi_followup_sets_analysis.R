# Ensure a default CRAN mirror is set for non-interactive runs
.repos <- getOption("repos")
.cran <- unname(.repos["CRAN"])
if (is.null(.cran) || is.na(.cran) || .cran == "@CRAN@" || .cran == "") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
rm(.repos, .cran)

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  readr, dplyr, tidyr, stringr, purrr, tibble, glue,
  ggplot2, DESeq2
)
source("R/utils.R")

in_dir <- "GTExDatav10"
gct_sc <- file.path(in_dir, "gene_reads_v10_adipose_subcutaneous.gct.gz")
gct_vis <- file.path(in_dir, "gene_reads_v10_adipose_visceral_omentum.gct.gz")
sattr <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt")
subjph <- file.path(in_dir, "GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt")

sets_tsv <- "references/kenichi_followup_gene_sets_2026-03-04.tsv"
frac_sc <- "GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/cibersort_subcutaneous_perm0_topN50_minCells200_noMALAT1__fractions_with_meta.tsv.gz"
frac_vis <- "GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/cibersort_visceral_perm0_topN50_minCells200_noMALAT1__fractions_with_meta.tsv.gz"

out_dir <- "GTEx_v10_AT_analysis_out"
tab_dir <- file.path(out_dir, "tables", "next_steps_2026-03-04_kenichi_followup")
fig_dir <- file.path(out_dir, "figs", "next_steps_2026-03-04_kenichi_followup")
log_dir <- file.path(out_dir, "logs")
ensure_dirs(c(tab_dir, fig_dir, log_dir))

log_path <- file.path(log_dir, "next_steps_2026-03-04_kenichi_followup.log")
sink(log_path, split = TRUE)
on.exit(sink(), add = TRUE)

score_set_z <- function(expr_sym, genes) {
  genes <- intersect(genes, rownames(expr_sym))
  if (length(genes) < 2) return(rep(NA_real_, ncol(expr_sym)))
  m <- expr_sym[genes, , drop = FALSE]
  mu <- rowMeans(m, na.rm = TRUE)
  sdv <- apply(m, 1, sd, na.rm = TRUE)
  sdv[is.na(sdv) | sdv == 0] <- 1
  z <- sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
  colMeans(z, na.rm = TRUE)
}

lm_contrast <- function(fit, term_a, term_b) {
  cf <- stats::coef(fit)
  vc <- stats::vcov(fit)
  v <- rep(0, length(cf)); names(v) <- names(cf)
  if (!is.null(term_a) && nzchar(term_a) && term_a %in% names(v)) v[term_a] <- v[term_a] + 1
  if (!is.null(term_b) && nzchar(term_b) && term_b %in% names(v)) v[term_b] <- v[term_b] - 1
  est <- sum(v * cf)
  se <- sqrt(as.numeric(t(v) %*% vc %*% v))
  tval <- est / se
  p <- 2 * stats::pt(abs(tval), df = fit$df.residual, lower.tail = FALSE)
  tibble::tibble(estimate = est, se = se, t = tval, p = p)
}

cat("Loading gene sets...\n")
set_df <- readr::read_tsv(sets_tsv, show_col_types = FALSE)
gene_sets <- set_df %>%
  mutate(genes = str_split(toupper(genes), "\\s*,\\s*")) %>%
  select(set_name, category, genes)

cat("Loading metadata...\n")
sample_attr <- readr::read_tsv(sattr, show_col_types = FALSE) %>%
  rename_with(~ gsub("^X", "", .), everything()) %>%
  mutate(SUBJID = stringr::str_match(SAMPID, "^([^-]+-[^-]+)")[, 2])
subject_pheno <- readr::read_tsv(subjph, show_col_types = FALSE) %>%
  rename_with(~ gsub("^X", "", .), everything())
meta_full <- left_join(sample_attr, subject_pheno, by = "SUBJID")

frac_sc_df <- readr::read_tsv(frac_sc, show_col_types = FALSE) %>% mutate(depot = "subcutaneous")
frac_vis_df <- readr::read_tsv(frac_vis, show_col_types = FALSE) %>% mutate(depot = "visceral")
frac_df <- bind_rows(frac_sc_df, frac_vis_df)

depots <- list(
  list(key = "subcutaneous", label = "Adipose - Subcutaneous", gct_path = gct_sc),
  list(key = "visceral", label = "Adipose - Visceral (Omentum)", gct_path = gct_vis)
)

score_rows <- list()
model_rows <- list()
corr_rows <- list()

for (dpt in depots) {
  cat("\n=== Depot:", dpt$key, "===\n")
  gct <- read_gct_v12(dpt$gct_path)
  common_samples <- intersect(colnames(gct$counts), meta_full$SAMPID)

  meta <- meta_full %>%
    filter(SAMPID %in% common_samples, SMTS == "Adipose Tissue", SMTSD == dpt$label) %>%
    mutate(
      sex = case_when(SEX == 1 ~ "Male", SEX == 2 ~ "Female", TRUE ~ NA_character_),
      age_bin_label = AGE
    ) %>%
    filter(!is.na(sex), !is.na(age_bin_label), !is.na(SMCENTER)) %>%
    bind_cols(parse_age_bin(.$AGE)) %>%
    derive_age_numeric() %>%
    mutate(
      depot = dpt$key,
      sex = factor(sex, levels = c("Male", "Female")),
      SMCENTER = factor(SMCENTER),
      menopause_group_3 = case_when(
        is.na(age_mid) ~ NA_character_,
        age_mid < 45 ~ "pre",
        age_mid <= 55 ~ "peri",
        age_mid > 55 ~ "post",
        TRUE ~ NA_character_
      ),
      menopause_group_3 = factor(menopause_group_3, levels = c("pre", "peri", "post"))
    ) %>%
    filter(!is.na(menopause_group_3))

  counts <- gct$counts[, meta$SAMPID, drop = FALSE]
  counts <- collapse_dupe_rowsum(counts)
  counts <- round(counts)
  storage.mode(counts) <- "integer"

  cat("Computing VST...\n")
  dds0 <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = data.frame(row.names = colnames(counts)),
    design = ~ 1
  )
  vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds0, blind = TRUE))
  expr_sym <- collapse_to_symbol_mean(vst_mat, gct$gene_id, gct$gene_symbol, uppercase = TRUE)

  sc_tbl <- tibble(SAMPID = colnames(expr_sym))
  for (i in seq_len(nrow(gene_sets))) {
    nm <- gene_sets$set_name[i]
    gs <- gene_sets$genes[[i]]
    sc_tbl[[paste0("score__", nm)]] <- score_set_z(expr_sym, gs)
  }

  df <- meta %>% left_join(sc_tbl, by = "SAMPID") %>%
    left_join(frac_df %>% select(SAMPID, adipocyte, aspc, macrophage, endothelial, t_cell), by = "SAMPID")

  score_rows[[length(score_rows) + 1]] <- df %>%
    select(SAMPID, depot, sex, AGE, age_mid, menopause_group_3, SMCENTER, starts_with("score__"),
           adipocyte, aspc, macrophage, endothelial, t_cell)

  # Female-focused models for menopause question
  dff <- df %>% filter(sex == "Female")
  for (i in seq_len(nrow(gene_sets))) {
    nm <- gene_sets$set_name[i]
    catg <- gene_sets$category[i]
    col <- paste0("score__", nm)
    d0 <- dff %>% filter(!is.na(.data[[col]]))
    if (nrow(d0) < 30) next

    fit_base <- lm(as.formula(paste0(col, " ~ SMCENTER + menopause_group_3")), data = d0)
    peri_pre_b <- lm_contrast(fit_base, "menopause_group_3peri", "")
    post_pre_b <- lm_contrast(fit_base, "menopause_group_3post", "")

    d1 <- d0 %>% filter(!if_any(c(adipocyte, aspc, macrophage, endothelial, t_cell), is.na))
    if (nrow(d1) >= 30) {
      fit_adj <- lm(as.formula(paste0(col, " ~ SMCENTER + menopause_group_3 + adipocyte + aspc + macrophage + endothelial + t_cell")), data = d1)
      peri_pre_a <- lm_contrast(fit_adj, "menopause_group_3peri", "")
      post_pre_a <- lm_contrast(fit_adj, "menopause_group_3post", "")
    } else {
      peri_pre_a <- tibble(estimate = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_)
      post_pre_a <- tibble(estimate = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_)
    }

    model_rows[[length(model_rows) + 1]] <- tibble(
      depot = dpt$key,
      set_name = nm,
      category = catg,
      n_base = nrow(d0),
      n_adj = nrow(d1),
      peri_vs_pre_est_base = peri_pre_b$estimate,
      peri_vs_pre_p_base = peri_pre_b$p,
      post_vs_pre_est_base = post_pre_b$estimate,
      post_vs_pre_p_base = post_pre_b$p,
      peri_vs_pre_est_adj = peri_pre_a$estimate,
      peri_vs_pre_p_adj = peri_pre_a$p,
      post_vs_pre_est_adj = post_pre_a$estimate,
      post_vs_pre_p_adj = post_pre_a$p
    )

    for (fc in c("adipocyte", "aspc", "macrophage", "endothelial", "t_cell")) {
      dcor <- d0 %>% filter(!is.na(.data[[fc]]))
      if (nrow(dcor) < 20) next
      ct <- suppressWarnings(cor.test(dcor[[col]], dcor[[fc]], method = "spearman"))
      corr_rows[[length(corr_rows) + 1]] <- tibble(
        depot = dpt$key,
        sex = "Female",
        set_name = nm,
        category = catg,
        fraction = fc,
        rho = unname(ct$estimate),
        p = ct$p.value,
        n = nrow(dcor)
      )
    }
  }
}

scores_all <- bind_rows(score_rows)
models <- bind_rows(model_rows)
corrs <- bind_rows(corr_rows)

readr::write_tsv(scores_all, file.path(tab_dir, "kenichi_followup_scores_per_sample.tsv"))
readr::write_tsv(models, file.path(tab_dir, "kenichi_followup_female_model_effects.tsv"))
readr::write_tsv(corrs, file.path(tab_dir, "kenichi_followup_female_fraction_correlations.tsv"))

# Figure 1: Female post-vs-pre effects (base vs +fractions) heatmap
plot1 <- models %>%
  select(depot, set_name, category, post_vs_pre_est_base, post_vs_pre_p_base, post_vs_pre_est_adj, post_vs_pre_p_adj) %>%
  pivot_longer(cols = c(post_vs_pre_est_base, post_vs_pre_est_adj), names_to = "model", values_to = "estimate") %>%
  mutate(model = recode(model,
                        post_vs_pre_est_base = "Base: +SMCENTER",
                        post_vs_pre_est_adj = "Adj: +CIBERSORT"),
         set_label = paste0(set_name, " [", category, "]"))

p1 <- ggplot(plot1, aes(x = model, y = set_label, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", estimate)), size = 2.7) +
  facet_wrap(~ depot, nrow = 1, scales = "free_y") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       name = "Effect\n(post-pre)") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(title = "Female menopause effect by module (post vs pre)",
       subtitle = "Kenichi follow-up sets; values are LM coefficients on module z-scores",
       x = NULL, y = NULL)

ggsave(file.path(fig_dir, "kenichi_followup_female_post_pre_effects_heatmap.png"), p1,
       width = 12.2, height = 6.8, dpi = 240)

# Figure 2: Female correlations with fractions heatmap
plot2 <- corrs %>%
  mutate(set_label = paste0(set_name, " [", category, "]"),
         fraction = factor(fraction, levels = c("adipocyte", "aspc", "macrophage", "endothelial", "t_cell")))

p2 <- ggplot(plot2, aes(x = fraction, y = set_label, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.6) +
  facet_wrap(~ depot, nrow = 1, scales = "free_y") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
                       limits = c(-1, 1), name = "Spearman\nrho") +
  theme_bw(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1)) +
  labs(title = "Female module score correlations with CIBERSORT fractions",
       subtitle = "Composition-linked signals are expected in baseline bulk transcriptomes",
       x = NULL, y = NULL)

ggsave(file.path(fig_dir, "kenichi_followup_female_fraction_correlations_heatmap.png"), p2,
       width = 12.2, height = 6.8, dpi = 240)

# Figure 3: Subcutaneous female selected modules by pre/peri/post
sel <- c("inflammation_cytokines", "fibrosis_ecm", "tgfb_alk7_axis", "adrenergic_beta", "adrenergic_alpha2", "dnl_lipogenesis")
plot3 <- scores_all %>%
  filter(sex == "Female", depot == "subcutaneous") %>%
  select(menopause_group_3, starts_with("score__")) %>%
  pivot_longer(cols = starts_with("score__"), names_to = "set_name", values_to = "score") %>%
  mutate(set_name = sub("^score__", "", set_name)) %>%
  filter(set_name %in% sel)

p3 <- ggplot(plot3, aes(x = menopause_group_3, y = score, fill = menopause_group_3)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.35) +
  geom_jitter(width = 0.18, alpha = 0.28, size = 0.8) +
  facet_wrap(~ set_name, scales = "free_y", ncol = 3) +
  theme_bw(base_size = 9.5) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(title = "Subcutaneous females: selected follow-up module scores", x = NULL, y = "Module z-score")

ggsave(file.path(fig_dir, "kenichi_followup_subcutaneous_female_selected_boxplots.png"), p3,
       width = 12.2, height = 7.0, dpi = 240)

cat("Done.\n")
