################################################################################
## Institutional Batch Effect Analysis
## Investigates batch differences between sequencing centers (B1, C1, D1)
################################################################################
# Ensure a default CRAN mirror is set for non-interactive runs
.repos <- getOption("repos")
.cran <- unname(.repos["CRAN"])
if (is.null(.cran) || is.na(.cran) || .cran == "@CRAN@" || .cran == "") {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}
rm(.repos, .cran)


if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(data.table, dplyr, tidyr, ggplot2, readr, patchwork, stringr)

## Load metadata
meta <- readr::read_tsv("GTEx_v10_AT_analysis_out/tables/sample_metadata_full.tsv",
                       show_col_types = FALSE)

## Clean age bin labels
meta$age_bin_label <- gsub("\\s*–\\s*|\\s*-\\s*", "-", meta$AGE)

## Derive age_mid from AGE column
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
  (age_min + age_max) / 2
}

meta$age_mid <- parse_age_bin(meta$AGE)

## Filter to relevant samples
meta_filt <- meta %>%
  filter(DEPOT %in% c("Subcutaneous", "Visceral"))

cat("\n===========================================\n")
cat("INSTITUTIONAL BATCH EFFECT ANALYSIS\n")
cat("===========================================\n\n")

## 1. Sample distribution by sequencing center
cat("1. Sample Distribution by Sequencing Center (SMCENTER):\n")
cat("   ---------------------------------------------------\n")

center_dist <- meta_filt %>%
  group_by(SMCENTER) %>%
  summarise(
    n = n(),
    subcutaneous = sum(DEPOT == "Subcutaneous"),
    visceral = sum(DEPOT == "Visceral"),
    female = sum(sex == "Female"),
    male = sum(sex == "Male"),
    age_mean = mean(age_mid, na.rm = TRUE),
    age_sd = sd(age_mid, na.rm = TRUE),
    .groups = "drop"
  )

print(center_dist)
cat("\n")

## 2. Age distribution within each center
cat("2. Age Group Distribution by Center:\n")
cat("   ----------------------------------\n")

age_center <- table(meta_filt$age_bin_label, meta_filt$SMCENTER)
print(age_center)
cat("\n")

# Chi-square test
chi_age_center <- chisq.test(age_center)
cat(sprintf("Chi-square p-value: %g\n", chi_age_center$p.value))
if (chi_age_center$p.value < 0.05) {
  cat("*** WARNING: Age groups are unevenly distributed across centers! ***\n")
}
cat("\n")

## 3. Sex distribution within each center
cat("3. Sex Distribution by Center:\n")
cat("   ----------------------------\n")

sex_center <- table(meta_filt$sex, meta_filt$SMCENTER)
print(sex_center)
cat("\n")

chi_sex_center <- tryCatch(chisq.test(sex_center), error = function(e) NULL)
if (!is.null(chi_sex_center)) {
  cat(sprintf("Chi-square p-value: %g\n", chi_sex_center$p.value))
  if (chi_sex_center$p.value < 0.05) {
    cat("*** WARNING: Sex is unevenly distributed across centers! ***\n")
  }
}
cat("\n")

## 4. Tissue depot distribution by center
cat("4. Tissue Depot Distribution by Center:\n")
cat("   ------------------------------------\n")

depot_center <- table(meta_filt$DEPOT, meta_filt$SMCENTER)
print(depot_center)
cat("\n")

chi_depot_center <- chisq.test(depot_center)
cat(sprintf("Chi-square p-value: %g\n", chi_depot_center$p.value))
if (chi_depot_center$p.value < 0.05) {
  cat("*** WARNING: Tissue depots are unevenly distributed across centers! ***\n")
}
cat("\n")

## 5. Create visualizations
cat("5. Creating visualizations...\n")

fig_dir <- "GTEx_v10_AT_analysis_out/figs"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Figure 1: Sample counts by center and depot
p1 <- ggplot(meta_filt, aes(x = SMCENTER, fill = DEPOT)) +
  geom_bar(position = "dodge") +
  labs(
    title = "Sample Distribution by Sequencing Center and Tissue Depot",
    x = "Sequencing Center (SMCENTER)",
    y = "Sample Count",
    fill = "Tissue Depot"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(fig_dir, "institution_sample_distribution.png"),
  p1, width = 8, height = 5, dpi = 300
)

# Figure 2: Age group distribution by center
p2 <- ggplot(meta_filt, aes(x = SMCENTER, fill = age_bin_label)) +
  geom_bar(position = "dodge") +
  labs(
    title = "Age Group Distribution by Sequencing Center",
    x = "Sequencing Center (SMCENTER)",
    y = "Sample Count",
    fill = "Age Group"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(fig_dir, "institution_age_distribution.png"),
  p2, width = 8, height = 5, dpi = 300
)

# Figure 3: Sex distribution by center
p3 <- ggplot(meta_filt, aes(x = SMCENTER, fill = sex)) +
  geom_bar(position = "dodge") +
  labs(
    title = "Sex Distribution by Sequencing Center",
    x = "Sequencing Center (SMCENTER)",
    y = "Sample Count",
    fill = "Sex"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(fig_dir, "institution_sex_distribution.png"),
  p3, width = 8, height = 5, dpi = 300
)

# Figure 4: Age midpoint by center and depot (boxplot)
p4 <- ggplot(meta_filt, aes(x = SMCENTER, y = age_mid, fill = DEPOT)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(
    title = "Age Distribution by Sequencing Center and Tissue Depot",
    x = "Sequencing Center (SMCENTER)",
    y = "Age Midpoint",
    fill = "Tissue Depot"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  file.path(fig_dir, "institution_age_boxplot.png"),
  p4, width = 8, height = 5, dpi = 300
)

# Figure 5: Combined summary panel
combined <- (p1 / p2) | (p3 / p4)
ggsave(
  file.path(fig_dir, "institution_batch_panel.png"),
  combined, width = 12, height = 10, dpi = 300
)

cat("   Saved visualizations to", fig_dir, "\n\n")

## 6. Summary statistics for batch effect consideration
cat("===========================================\n")
cat("SUMMARY: INSTITUTIONAL BATCH EFFECTS\n")
cat("===========================================\n\n")

cat("Sample Distribution:\n")
for (center in unique(meta_filt$SMCENTER)) {
  n <- sum(meta_filt$SMCENTER == center, na.rm = TRUE)
  pct <- 100 * n / nrow(meta_filt)
  cat(sprintf("  %s: %d samples (%.1f%%)\n", center, n, pct))
}
cat("\n")

cat("Potential Batch Effects to Consider:\n")
cat("  1. Center-specific technical biases (library prep, sequencing)\n")
cat("  2. Center-specific protocols or reagent lots\n")
cat("  3. Temporal effects (samples processed at different times)\n\n")

cat("Recommendations:\n")
if (chi_age_center$p.value < 0.05) {
  cat("  *** Age distribution differs by center - INCLUDE SMCENTER as covariate ***\n")
}
if (!is.null(chi_sex_center) && chi_sex_center$p.value < 0.05) {
  cat("  *** Sex distribution differs by center - CHECK for sex-specific batch effects ***\n")
}
if (chi_depot_center$p.value < 0.05) {
  cat("  *** Depot distribution differs by center - CHECK for depot-specific batch effects ***\n")
}

cat("  General: Consider ComBat or SVA for batch correction in DE analysis\n")
cat("  General: Always include SMCENTER as a covariate in linear models\n")
cat("===========================================\n\n")
