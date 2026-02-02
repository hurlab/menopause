# Run Journal (2026-02-01)
## Menopause / SNS / Lipolysis project (GTEx v10)

Purpose: track long-running analysis steps so the work can be resumed cleanly if interrupted.

---

### 2026-02-01

- Started implementing the next-steps plan in `docs/next_steps_plan_menopause_sns_gtex_2026-02-01.md`.
- Planned outputs are grouped into dated “next_steps_2026-02-01_*” subfolders under `GTEx_v10_AT_analysis_out/{tables,figs,logs}/`.
- Completed covariate sensitivity ladder (both depots): `R/next_steps_2026-02-01_covariate_sensitivity.R`
  - Outputs: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_covariate_sensitivity/` and `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_covariate_sensitivity/`
  - Note: adding gene-expression batch (SMGEBTCH) can strongly attenuate/flip effects in some contrasts; interpret as sensitivity/collinearity check rather than “best” model.
- Completed marker-based composition proxy analysis (both depots): `R/next_steps_2026-02-01_composition_markers.R`
  - Outputs: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_composition/` and `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/`
  - Key observation: subcutaneous adipocyte marker score (ADIPOQ/PLIN1/FABP4) decreases in post vs pre; adjusting Kenichi gene-set models for composition proxies largely removes the subcutaneous “lipolysis_core down” signal (suggesting bulk composition/adipocyte content is a major driver).
- Completed expanded adrenergic/cAMP/desensitization module scoring (both depots): `R/next_steps_2026-02-01_adrenergic_modules.R`
  - Outputs: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_adrenergic_modules/` and `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_adrenergic_modules/`
  - Key observation (subcutaneous): cAMP core module trends down in post vs pre (even after composition adjustment), and a desensitization module shows a modest downshift after composition adjustment (consistent with altered responsiveness rather than acute activation at baseline).
- Implemented first-pass menopause signature refinement + scoring (v1 data-driven + v2 curated): `R/next_steps_2026-02-01_signature_refinement.R`
  - Outputs: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_signature_refinement/` and `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_signature_refinement/`
  - Signature files: `references/menopause_signature_senMayoLike_v1.tsv`, `references/menopause_signature_senMayoLike_v2.tsv`
  - Evidence skeletons (to fill with PMIDs/DOIs): `references/menopause_signature_senMayoLike_v1_evidence.tsv`, `references/menopause_signature_senMayoLike_v2_evidence.tsv`
  - Note: v2 is the intended “biology-facing” SenMayo-like candidate (estrogen/Wnt/ECM/adipose genes); next step is to populate gene-by-gene evidence from the literature and update the evidence TSV.
- Populated gene-by-gene evidence (PMID-based) for the curated signature:
  - `references/menopause_signature_senMayoLike_v2_evidence.tsv`
  - Dropped GALNT5 from v2 due to insufficient direct evidence during this pass.
- External dataset processed:
  - GSE44000 (postmenopausal human subcutaneous adipocytes; obese vs lean): `R/external_validation_GSE44000_postmenopausal_adipocytes_obese_vs_lean.R`
  - Outputs: `GTEx_v10_AT_analysis_out/tables/external_validation_GSE44000/` and `GTEx_v10_AT_analysis_out/figs/external_validation_GSE44000/`
  - Note: not a pre/post menopause validation; used to contextualize adipocyte program behavior without bulk composition confounding.
- Wrote a biology-facing master summary: `docs/master_summary_menopause_sns_gtex_2026-02-01.md`
- Created an updated slide outline and generated a pptx via pandoc:
  - Outline: `docs/ppt_outline_menopause_sns_gtex_master_2026-02-01.md`
  - PPTX: `docs/menopause_sns_gtex_master_2026-02-01.pptx`
- Created an enhanced, figure-forward PPTX (images + table PNGs):
  - Table-image generator: `R/make_ppt_tables_2026-02-01.R`
  - Deck markdown: `docs/ppt_deck_enhanced_menopause_sns_gtex_2026-02-01.md`
  - PPTX: `docs/menopause_sns_gtex_enhanced_2026-02-01.pptx`
- Updated enhanced deck with additional method explanations + a summary slide and exported senMayoLike v2 details:
  - Updated deck markdown: `docs/ppt_deck_enhanced_menopause_sns_gtex_2026-02-02.md`
  - Updated PPTX: `docs/menopause_sns_gtex_enhanced_2026-02-02.pptx`
  - senMayoLike v2 details: `docs/senMayoLike_v2_signature_details_2026-02-02.csv`, `docs/senMayoLike_v2_signature_details_2026-02-02.xlsx`
- Added male vs female context analyses (score-level) and integrated into a refreshed figure-forward slide deck:
  - Sex×age-bin interaction table + heatmap: `R/next_steps_2026-02-02_sex_interaction_gene_set_scores.R`
  - Sex-by-age-bin trajectories plot (key scores): `R/next_steps_2026-02-02_sex_age_trends_plots.R`
  - Updated deck markdown: `docs/ppt_deck_enhanced_menopause_sns_gtex_2026-02-02_v2.md`
  - Updated PPTX: `docs/menopause_sns_gtex_enhanced_2026-02-02_v2.pptx`
- Added explicit “batch/center confounding” slides near the front of the deck:
  - Figures: `GTEx_v10_AT_analysis_out/figs/institution_batch_panel.png`, `GTEx_v10_AT_analysis_out/figs/batch_correction_summary_panel.png`, `GTEx_v10_AT_analysis_out/figs/PCA_combat_correction_comparison.png`
  - Updated deck markdown: `docs/ppt_deck_enhanced_menopause_sns_gtex_2026-02-02_v3.md`
  - Updated PPTX: `docs/menopause_sns_gtex_enhanced_2026-02-02_v3.pptx`
