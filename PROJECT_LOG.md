# Project Log (Append-Only)

Session 2026-02-28 08:56
- Coding CLI used: OpenAI Codex CLI
- Phase(s) worked on: TODO #2 (receptor heatmaps); TODO #3 (Slides 7-8 sex display); TODO #4 (CIBERSORT fractions).
- Concrete changes implemented: Added receptor expression heatmaps (group-mean VST) for ADRB1/2/3 and ADRA2A/B/C (faceted by `sex x depot`, x-axis = age bins).
- Concrete changes implemented: Added Slide 7-8 style reproduction of the Kenichi panel plots (lipolysis_core, thermogenesis_program) with sex facets and per-sex LM p-values.
- Concrete changes implemented: Implemented a minimal CIBERSORT workflow using GSE176171 as reference (downloaded MatrixMarket counts + metadata; built CPM pseudo-bulk signature; ran CIBERSORT on GTEx subcutaneous and visceral mixtures).
- Concrete changes implemented: Downloaded GSE281356 raw 10X `.h5` matrices for potential future reference use (not used yet).
- Concrete changes implemented: Added an output index markdown to point to the newly generated figures/tables.
- Files/modules/functions touched: Added `R/next_steps_2026-02-28_adrenergic_receptor_expression_heatmaps.R`.
- Files/modules/functions touched: Added `R/next_steps_2026-02-28_kenichi_panel_by_sex_slides7_8_style.R`.
- Files/modules/functions touched: Added `R/next_steps_2026-02-28_cibersort_cell_fractions_adipose_reference_GSE176171.R`.
- Files/modules/functions touched: Added `docs/next_steps_2026-02-28_items_2_3_4_outputs.md`.
- Files/modules/functions touched: Added `PROJECT_HANDOFF.md` and `PROJECT_LOG.md`.
- External data added/updated: `external_data/GSE176171/` (cell metadata + Hs10X MatrixMarket inputs).
- External data added/updated: `external_data/GSE281356/GSE281356_RAW.tar`.
- Key technical decisions and rationale: Reused existing GTEx loading + VST workflows to keep figures consistent with prior decks.
- Key technical decisions and rationale: For Slide 7-8 male/female display, computed z-scoring within `depot+sex` to keep scales within sex comparable.
- Key technical decisions and rationale: Chose GSE176171 as the smallest practical reference because it includes both a combined matrix and cell-type labels; GSE281356 does not ship curated labels in GEO.
- Key technical decisions and rationale: Stripped Ensembl version suffixes in reference and GTEx to fix near-zero overlap.
- Key technical decisions and rationale: Used a fast/stable first pass for CIBERSORT (`TopN=50`, `perm=0`, `mc.cores=1`) due to runtime and parallelism fragility at larger settings.
- Problems encountered and resolutions: `pacman` install failed due to a stale lock and missing CRAN mirror; fixed by explicitly setting CRAN and reinstalling.
- Problems encountered and resolutions: Matrix coercion deprecation (`dgTMatrix` -> `dgCMatrix`); switched to `as(., "CsparseMatrix")`.
- Problems encountered and resolutions: `IOBR::CIBERSORT` returns a matrix; wrapped output with `as.data.frame()` + rownames -> `Mixture` column before writing CSV.
- Items explicitly completed, resolved, or superseded in this session: Completed TODO #2 and TODO #3; partially completed TODO #4 (pipeline + outputs generated; validation/tuning remain).
- Verification performed: Ran `Rscript --vanilla R/next_steps_2026-02-28_adrenergic_receptor_expression_heatmaps.R`.
- Verification performed: Ran `Rscript --vanilla R/next_steps_2026-02-28_kenichi_panel_by_sex_slides7_8_style.R`.
- Verification performed: Ran `CIBERSORT_SIGNATURE_TOP_N=50 CIBERSORT_PERM=0 CIBERSORT_CORES=1 Rscript --vanilla R/next_steps_2026-02-28_cibersort_cell_fractions_adipose_reference_GSE176171.R`.
- Verification performed: Confirmed outputs created under `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_*` and `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_*`.

Session 2026-02-28 10:04
- Coding CLI used: OpenAI Codex CLI
- Phase(s) worked on: Glycolysis hypothesis synthesis + PPTX deliverables; project deck update.
- Concrete changes implemented: Created a new glycolysis-focused hypothesis deck embedding the latest GTEx hallmark plots and a model-summary table (subq/visc; base vs +CIBERSORT).
- Concrete changes implemented: Created an updated project summary deck by copying `docs/menopause_sns_gtex_enhanced_2026-02-02_v4.1.2.pptx`, updating the title-slide date, and appending new “2026-02-28 updates” slides for receptor heatmaps, sex-faceted Kenichi panels, CIBERSORT fraction plots, and glycolysis/hallmarks plots.
- Concrete changes implemented: Fixed a pre-existing out-of-bounds textbox on Slide 21 during the deck update; ensured newly added shapes stay within slide boundaries.
- Files/modules/functions touched: Added `scripts/make_pptx_2026-02-28_glycolysis_and_update_project.py`.
- Files/modules/functions touched: Wrote `docs/menopause_glycolysis_hypothesis_2026-02-28.pptx`.
- Files/modules/functions touched: Wrote `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.1.3.pptx`.
- Key technical decisions and rationale: Used the existing v4.1.2 PPTX as a theme/template to keep formatting consistent; kept edits to the main deck append-only to preserve established slide content.
- Problems encountered and resolutions: `python-pptx` was not available; installed it, then fixed slide-deletion logic (relationship IDs) and adjusted schematic/table slide geometry to avoid out-of-bounds shapes.
- Items explicitly completed, resolved, or superseded in this session: Completed PPTX deliverables requested in this session (glycolysis hypothesis deck; updated project summary deck).
- Verification performed: Opened both PPTX via `python-pptx` and ran a shape bounding-box check to confirm no shapes extend past slide boundaries (2026-02-28 10:04).

Session 2026-02-28 15:02
- Coding CLI used: OpenAI Codex CLI
- Phase(s) worked on: Deck revisions requested from review of v4.1.3 (receptors, Kenichi sex display, CIBERSORT plots/reference); rebuilt main deck as v4.2.
- Concrete changes implemented: Extended receptor heatmap outputs to include (a) group-mean heatmaps with mean VST values printed in cells (2 decimals) and (b) per-sample VST heatmaps ordered by age-bin within each (sex × depot) facet.
- Concrete changes implemented: Added a Kenichi panel re-plot with female/male interleaved groups on a single x-axis (F-pre, M-pre, F-peri, M-peri, F-post, M-post) for easier cross-sex comparison.
- Concrete changes implemented: Fixed the CIBERSORT plotting bug where the output column header is `P-value` (dash) not `P.value` (dot), which caused the `P-value` column to be treated as a “cell type” in stacked bar plots; regenerated “clean” fraction plots using signature cell-type columns only.
- Concrete changes implemented: Added a CIBERSORT reference summary plot (GSE176171 cell counts by cell type) to support a new reference-background slide.
- Concrete changes implemented: Generated updated project deck `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.pptx` by copying v4.1.2 and appending an expanded update section: receptor (means + values + per-sample), Kenichi (sex-faceted + interleaved), CIBERSORT reference slide + cleaned fraction plots, glycolysis/hallmarks, updated outputs list.
- Files/modules/functions touched: Updated `R/next_steps_2026-02-28_adrenergic_receptor_expression_heatmaps.R`.
- Files/modules/functions touched: Added `R/next_steps_2026-02-28_kenichi_panel_FM_interleaved_single_plot.R`.
- Files/modules/functions touched: Added `R/next_steps_2026-02-28_cibersort_reference_slide_and_clean_fraction_plots.R`.
- Files/modules/functions touched: Updated `scripts/make_pptx_2026-02-28_glycolysis_and_update_project.py` (added v4.2 deck builder).
- Files/modules/functions touched: Wrote `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.pptx`.
- Problems encountered and resolutions: R plotting helper used a non-constant `slice_head(n = n()-1)`; replaced with a group-wise `row_number()==n()` filter.
- Verification performed: Regenerated figures and deck; ran a shape bounding-box check to confirm no shapes extend past slide boundaries in `v4.2` (2026-02-28 15:02).

Session 2026-02-28 15:16
- Coding CLI used: OpenAI Codex CLI
- Phase(s) worked on: Deck polish from review (VST comparability questions + slide improvements requested).
- Concrete changes implemented: Updated per-sample receptor heatmaps to include age-bin labels (annotated above the heatmap) while keeping individual sample IDs hidden.
- Concrete changes implemented: Generated interleaved-sex Kenichi plots with median values annotated (2 decimals) and added as an extra slide.
- Concrete changes implemented: Split CIBERSORT fraction plots into two separate, larger slides (subcutaneous and visceral) for readability.
- Concrete changes implemented: Added the GSE176171 publication citation to the reference slide.
- Concrete changes implemented: Generated updated deck `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx`.
- Files/modules/functions touched: Updated `R/next_steps_2026-02-28_adrenergic_receptor_expression_heatmaps.R` (age-bin labels on per-sample heatmap).
- Files/modules/functions touched: Updated `R/next_steps_2026-02-28_kenichi_panel_FM_interleaved_single_plot.R` (median annotations; new output filenames).
- Files/modules/functions touched: Updated `scripts/make_pptx_2026-02-28_glycolysis_and_update_project.py` (v4.2.1 deck updates: extra Kenichi medians slide; split CIBERSORT slides; citation; output v4.2.1).
- Verification performed: Computed GTEx adipose sample counts by sex to confirm male>female (subcutaneous: 874 M vs 459 F; visceral: 758 M vs 370 F) and ran a shape bounding-box check to confirm no shapes extend past slide boundaries in `v4.2.1` (2026-02-28 15:16).

Session 2026-02-28 15:25
- Coding CLI used: OpenAI Codex CLI
- Phase(s) worked on: Deck tweak requested (add cell-type label explainer slide after CIBERSORT reference).
- Concrete changes implemented: Inserted a new slide after Slide 30 in `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx` explaining the GSE176171 broad cell-type labels used in the CIBERSORT signature (`aspc`, `lec`, `smc`, etc.).
- Files/modules/functions touched: Updated `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx` in place (no new file).
- Verification performed: Ran a shape bounding-box check to confirm no shapes extend past slide boundaries after insertion (2026-02-28 15:25).

Session 2026-02-28 18:15
- Coding CLI used: OpenAI Codex CLI
- Phase(s) worked on: Deck tweak requested (add clarification slide(s) on lipolysis vs glycolysis confusion).
- Concrete changes implemented: Added two interpretive slides after the glycolysis/hallmarks slide in `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx`: (1) definitions + interpretability cautions, (2) summary of what our results show and why the signals can be mixed up.
- Concrete changes implemented: Updated the final “new outputs” slide with a short deck note indicating the clarification slides were added.
- Files/modules/functions touched: Updated `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx` in place (no new file).
- Verification performed: Ran a shape bounding-box check to confirm no shapes extend past slide boundaries after insertion (2026-02-28 18:15).

Session 2026-02-28 18:19
- Coding CLI used: OpenAI Codex CLI
- Phase(s) worked on: Deck formatting fix (last slides were mis-positioned).
- Concrete changes implemented: Recreated the last three text-heavy slides using the same Title+Body layout geometry as Slide 2 (title at top, body below) and enabled text autofit with reduced font size to prevent overflow.
- Concrete changes implemented: Re-ordered the last slides so the “new outputs” slide remains last.
- Files/modules/functions touched: Updated `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx` in place (no new file).
- Verification performed: Ran a shape bounding-box check to confirm no shapes extend past slide boundaries after the formatting fix (2026-02-28 18:19).

Session 2026-02-28 19:22
- Coding CLI used: OpenAI Codex CLI
- Phase(s) worked on: CIBERSORT debugging (adipocyte fraction missing/near-zero in adipose tissues).
- Concrete findings: In the original CIBERSORT run (`cibersort_*_perm0_topN50_minCells200`), adipocyte fraction was essentially absent: subcutaneous mean ~4.5e-05 with ~99.4% exact zeros; visceral mean ~2e-04 with ~98.1% zeros.
- Concrete findings: Despite near-zero adipocyte fractions, adipocyte marker genes (e.g., ADIPOQ, PLIN1, LPL) are present and high in the CIBERSORT mixture matrix, indicating an estimation artifact rather than truly missing adipocyte signal.
- Key technical root cause: The signature matrix included MALAT1 (`ENSG00000251562`), an extreme high-expression gene that dominates scaling inside `IOBR::CIBERSORT` and collapses adipocyte weights to ~0.
- Concrete changes implemented: Added configurable exclusion of Ensembl IDs from signature+mixture in `R/next_steps_2026-02-28_cibersort_cell_fractions_adipose_reference_GSE176171.R` via `CIBERSORT_EXCLUDE_ENSG` (default excludes MALAT1). The script now produces a tagged run: `*_noMALAT1__*`.
- Concrete changes implemented: Re-ran CIBERSORT (perm=0) with MALAT1 excluded; adipocyte fraction recovered (subq mean ~0.083, median ~0.070; visc mean ~0.067, median ~0.054).
- Concrete changes implemented: Added a plotting helper `R/next_steps_2026-02-28_cibersort_plot_from_fractions_with_meta.R` and generated cleaned noMALAT1 plots:
  - `.../cibersort_subcutaneous_*_noMALAT1__fractions_mean_by_sex_agebin_clean.png`
  - `.../cibersort_visceral_*_noMALAT1__fractions_mean_by_sex_agebin_clean.png`
- Concrete QC: Verified adipocyte fraction correlates strongly with adipocyte markers in the mixture (e.g., ADIPOQ/PLIN1/LPL r~0.81–0.88) and ASPC fraction correlates with DCN/COL1A2.
- Files/modules/functions touched: Updated `R/next_steps_2026-02-28_cibersort_cell_fractions_adipose_reference_GSE176171.R`.
- Files/modules/functions touched: Added `R/next_steps_2026-02-28_cibersort_debug_adipocyte_missing.R`.
- Files/modules/functions touched: Added `R/next_steps_2026-02-28_cibersort_plot_from_fractions_with_meta.R`.

Session 2026-02-28 19:59
- Coding CLI used: OpenAI Codex CLI
- Phase(s) worked on: Add CIBERSORT debugging conclusions into the main deck (v4.2.1).
- Concrete changes implemented: Generated a reference-vs-bulk comparison scatter plot (reference cell-count % vs CIBERSORT mean fractions; noMALAT1) and saved it under `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_cibersort_debug/ref_pct_vs_cibersort_mean_noMALAT1.png`.
- Concrete changes implemented: Inserted 3 slides after the visceral fraction slide in `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx`:
  - Slide 34: issue summary + original adipocyte near-zero stats + root cause (MALAT1 scaling artifact)
  - Slide 35: mitigation summary + corrected adipocyte stats + marker coherence correlations
  - Slide 36: corrected noMALAT1 fraction plots (subq+visc) + reference-vs-mean scatter + interpretation note
- Concrete changes implemented: Appended a note to the final “new outputs” slide that CIBERSORT QC/corrections were added.
- Files/modules/functions touched: Added `R/next_steps_2026-02-28_cibersort_reference_vs_bulk_comparison_plot.R`.
- Files/modules/functions touched: Updated `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx` in place (no new file).
- Verification performed: Ran a shape bounding-box check to confirm no shapes extend past slide boundaries after insertion (2026-02-28 19:59).

Session 2026-03-04 18:53
- Coding CLI used: OpenAI Codex CLI
- Phase(s) worked on: Kenichi follow-up (inflammation/fibrosis/TGF-β axes; catecholamine resistance mechanisms).
- Concrete changes implemented:
  - Added a curated gene-sets TSV derived from `references/GeneForJunguk_030226.docx` plus review-derived mechanism axes from `references/nrendo.2017.90.pdf` (Fig 3): `references/kenichi_gene_sets_immune_fibrosis_2026-03-04.tsv`.
  - Implemented gene-set scoring + key-gene extraction + correlation QC vs CIBERSORT fractions (noMALAT1): `R/next_steps_2026-03-04_kenichi_inflammation_fibrosis_CA_resistance.R`.
  - Generated new figures and tables under `GTEx_v10_AT_analysis_out/{figs,tables}/next_steps_2026-03-04_kenichi_followup/`:
    - `kenichi_followup_gene_set_trends_by_age_sex__both_depots.png`
    - `kenichi_followup_correlations__subcutaneous.png`, `kenichi_followup_correlations__visceral.png`
    - score tables: `kenichi_followup_scores_and_key_genes_vst__*.tsv`
  - Created a new dated deck by copying `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx` and appending a “Kenichi follow-up” section at the end: `docs/menopause_sns_gtex_enhanced_2026-03-04_v4.3.0.pptx` (6 new slides appended, including literature mechanism summary + the new plots + next-steps slide).
- Files/modules/functions touched:
  - Added: `references/kenichi_gene_sets_immune_fibrosis_2026-03-04.tsv`
  - Added: `R/next_steps_2026-03-04_kenichi_inflammation_fibrosis_CA_resistance.R`
  - Added: `docs/menopause_sns_gtex_enhanced_2026-03-04_v4.3.0.pptx`
- Key technical decisions and rationale:
  - Used z-score module scoring on VST (per depot) to make heterogeneous gene sets comparable and to focus on relative shifts across age/sex bins.
  - Focused the mechanism axes on the review’s Fig 3 nodes (TGF-β/ALK7 and TNF→NF-κB→IKKε/TBK1→PDE3B) because these are directly testable in bulk RNA and integrate naturally with the project’s existing adrenergic/lipolysis framing.
  - Kept the existing v4.2.1 deck immutable; created a new dated deck per user preference for new additions.
- Problems encountered and resolutions:
  - Section Header slide title in the new deck was initially blank due to placeholder identity checks; fixed by targeting placeholder `idx==0` (title) and `idx==1` (subtitle) and re-saving the deck.
- Items explicitly completed, resolved, or superseded in this session:
  - Completed initial Kenichi follow-up analysis pass and deck append (no changes to prior completed items).
- Verification performed (if any):
  - Verified R script runs end-to-end and writes expected outputs (2026-03-04 18:48).
  - Verified new deck exists and includes appended slides (2026-03-04 18:53).
