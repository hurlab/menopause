# Outputs For TODOs #2, #3, #4 (Run Date: 2026-02-28)

This note records what was generated for the meeting TODOs (from `docs/todos_from_me_2026-02-24.md`):

## Item #2: Adrenergic receptor expression (alpha vs beta; depot/sex/age)

Script:
- `R/next_steps_2026-02-28_adrenergic_receptor_expression_heatmaps.R`

Figures:
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_receptor_heatmaps/heatmap_ADRB1_ADRB2_ADRB3_meanVST_by_depot_sex_agebin.png`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_receptor_heatmaps/heatmap_ADRA2A_ADRA2B_ADRA2C_meanVST_by_depot_sex_agebin.png`

Table:
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_receptor_heatmaps/adrenergic_receptor_expression_group_means_vst.tsv`

## Item #3: Slides 7-8 style (Kenichi panel) with Male vs Female

Script:
- `R/next_steps_2026-02-28_kenichi_panel_by_sex_slides7_8_style.R`

Figures (sex-faceted, Slides 7-8 style):
- Subcutaneous:
  - `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_kenichi_panel_by_sex/figure_candidate__lipolysis_thermogenesis__by_sex.png`
- Visceral:
  - `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_kenichi_panel_by_sex/figure_candidate__lipolysis_thermogenesis__by_sex_visceral.png`

Table (LM p-values per sex):
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_kenichi_panel_by_sex/kenichi_panel_lm_pvalues_by_sex.tsv`

## Item #4: Cell fractions (CIBERSORT; first pass)

Reference used:
- GSE176171 (because it ships a ready-to-use Hs10X count matrix + a cell-level metadata table with broad cell types).

Script:
- `R/next_steps_2026-02-28_cibersort_cell_fractions_adipose_reference_GSE176171.R`

Reference signature (Top-N markers per cell type; CPM; Ensembl IDs without version):
- `external_data/GSE176171/derived/signature_topN50_cpm_noversion.tsv.gz`

CIBERSORT outputs:
- Subcutaneous:
  - `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/cibersort_subcutaneous_perm0_topN50_minCells200__cibersort_fractions.csv`
  - `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/cibersort_subcutaneous_perm0_topN50_minCells200__fractions_with_meta.tsv.gz`
  - `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_cibersort_fractions/cibersort_subcutaneous_perm0_topN50_minCells200__fractions_mean_by_sex_agebin.png`
- Visceral:
  - `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/cibersort_visceral_perm0_topN50_minCells200__cibersort_fractions.csv`
  - `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/cibersort_visceral_perm0_topN50_minCells200__fractions_with_meta.tsv.gz`
  - `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_cibersort_fractions/cibersort_visceral_perm0_topN50_minCells200__fractions_mean_by_sex_agebin.png`

