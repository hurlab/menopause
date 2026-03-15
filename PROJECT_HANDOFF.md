# Project Handoff (Living State)

## 1. Project Overview
- Purpose and scope: Analyze GTEx v10 adipose (subcutaneous + visceral omentum) to characterize menopause/SNS-related transcriptional patterns; add sex-stratified views, adrenergic receptor readouts, and robust cell composition estimation (CIBERSORT) to reduce bulk-mixture confounding.
- Last updated: 2026-03-04 18:53
- Last coding CLI used (informational): OpenAI Codex CLI

## 2. Current State
- Adrenergic receptor expression heatmaps (beta vs alpha2; depot/sex/age): Status = Completed; Completed in Session 2026-02-28 08:56; Outputs indexed in `docs/next_steps_2026-02-28_items_2_3_4_outputs.md`.
- Adrenergic receptor heatmaps with mean VST values in cells + per-sample heatmaps ordered by age bins: Status = Completed; Completed in Session 2026-02-28 15:02; Outputs under `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_receptor_heatmaps/`.
- Slides 7-8 style Kenichi panel with Male vs Female visualization: Status = Completed; Completed in Session 2026-02-28 08:56; Outputs indexed in `docs/next_steps_2026-02-28_items_2_3_4_outputs.md`.
- Kenichi panel with Female/Male interleaved groups in a single plot (F-pre, M-pre, F-peri, M-peri, F-post, M-post): Status = Completed; Completed in Session 2026-02-28 15:02; Outputs under `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_kenichi_panel_by_sex/`.
- Cell-fraction estimation (CIBERSORT; first pass using adipose sc reference): Status = In progress; Last updated = 2026-02-28 08:56; First-pass uses TopN=50 signature and perm=0 to unblock downstream modeling quickly.
- CIBERSORT plots cleaned to exclude statistics columns + reference summary plot: Status = Completed; Completed in Session 2026-02-28 15:02; Outputs under `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_cibersort_fractions/`.
- CIBERSORT adipocyte-missing issue: Status = Mitigated; Last updated = 2026-02-28 19:22; Finding = adipocyte fraction was ~0 due to MALAT1 dominating IOBR::CIBERSORT scaling; Fix = exclude `ENSG00000251562` (MALAT1); New outputs under `.../cibersort_*_noMALAT1__fractions_with_meta.tsv.gz` + `*_noMALAT1__fractions_mean_by_sex_agebin_clean.png`.
- CIBERSORT QC slides added to main deck: Status = Completed; Completed in Session 2026-02-28 19:59; Added after the visceral fraction slide in `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx` (problem, mitigation, corrected plots + reference comparison).
- Glycolysis hypothesis deck (PPTX): Status = Completed; Completed in Session 2026-02-28 10:04; Output = `docs/menopause_glycolysis_hypothesis_2026-02-28.pptx`.
- Updated project summary deck (PPTX): Status = Completed; Completed in Session 2026-02-28 10:04; Output = `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.1.3.pptx` (copied from v4.1.2; new results appended). Note = Superseded by v4.2.
- Updated project summary deck (PPTX; v4.2): Status = Completed; Completed in Session 2026-02-28 15:02; Output = `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.pptx` (expanded update section; fixed CIBERSORT plotting issue; added reference slide).
- Updated project summary deck (PPTX; v4.2.1): Status = Completed; Completed in Session 2026-02-28 15:16; Output = `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx` (added Kenichi median-annotated slide; split CIBERSORT subq/visc for readability; added GSE176171 publication citation; improved per-sample heatmap age-bin labeling).
- CIBERSORT cell type label explainer slide (aspc, lec, smc, etc.): Status = Completed; Completed in Session 2026-02-28 15:25; Added after Slide 30 in `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx`.
- Lipolysis vs glycolysis clarification slides: Status = Completed; Completed in Session 2026-02-28 18:15; Added after the glycolysis/hallmarks slide in `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx`.
- Fixed formatting of the last slides (title at top, body below, autofit): Status = Completed; Completed in Session 2026-02-28 18:19; Applied to the last slides in `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx`.
- Kenichi follow-up (inflammation/fibrosis/TGF-β axes + catecholamine resistance mechanism readouts): Status = Completed; Completed in Session 2026-03-04 18:53; Outputs under `GTEx_v10_AT_analysis_out/{figs,tables}/next_steps_2026-03-04_kenichi_followup/` and deck update in `docs/menopause_sns_gtex_enhanced_2026-03-04_v4.3.0.pptx`.

## 3. Execution Plan Status
- Phase A (Item #2 receptor expression visualization): Status = Completed; Last updated = 2026-02-28 08:56; Note = group-mean VST heatmaps faceted by `sex x depot` and age-bin.
- Phase B (Item #3 sex-stratified display for Slide 7-8 results): Status = Completed; Last updated = 2026-02-28 08:56; Note = Slide 7-8 style reproduced with sex facets; per-sex LM p-values exported.
- Phase C (Item #4 CIBERSORT fractions): Status = In progress; Last updated = 2026-02-28 08:56; Note = minimal CIBERSORT pipeline using GSE176171; next is validate fractions and optionally re-run with perm>0 and/or alternate reference.
- Phase D (Glycolysis discrepancy hypothesis + updated decks): Status = Completed; Last updated = 2026-02-28 19:59; Note = deck matured through v4.2.1 including glycolysis context + receptor/sex/CIBERSORT/QC materials.
- Phase E (Kenichi follow-up inflammation/fibrosis/TGF-β axes): Status = Completed (initial pass); Last updated = 2026-03-04 18:53; Note = gene-set scoring + correlation QC plots + appended deck section in v4.3.0.

## 4. Outstanding Work
- Validate and tune CIBERSORT fractions (TopN, perm, cell types included): Status = In progress; Last updated = 2026-02-28 08:56; Most recent session ref = PROJECT_LOG.md (Session 2026-02-28 08:56).
- Integrate CIBERSORT fractions into the existing GTEx models (replace marker-proxy covariates): Status = Not started; Last updated = 2026-02-28 08:56; Most recent session ref = PROJECT_LOG.md (Session 2026-02-28 08:56).
- Optionally extend receptor readouts beyond ADRA2/ADRB (other alpha/beta families; downstream signaling genes): Status = Not started; Last updated = 2026-02-28 08:56; Most recent session ref = PROJECT_LOG.md (Session 2026-02-28 08:56).
- Glycolysis follow-up: gene-node modeling and fraction correlations (beyond Hallmark score): Status = Not started; Last updated = 2026-02-28 10:04; Most recent session ref = PROJECT_LOG.md (Session 2026-02-28 10:04).
- Kenichi follow-up modeling: test menopause/sex effects on inflammation/fibrosis/TGF-β axes and whether they explain (or track with) adrenergic receptor/lipolysis readouts after composition adjustment: Status = Not started; Last updated = 2026-03-04 18:53; Most recent session ref = PROJECT_LOG.md (Session 2026-03-04 18:53).

## 5. Risks, Open Questions, and Assumptions
- CIBERSORT reference choice and comparability to GTEx bulk adipose: Status = Open; Date opened = 2026-02-28; Default assumption = GSE176171 broad cell types are acceptable for first pass; validate and run sensitivity checks (TopN/perm/QN/reference).
- CIBERSORT adipocyte fraction near-zero (artifact): Status = Resolved; Date opened = 2026-02-28; Resolution = exclude MALAT1 (`ENSG00000251562`) from signature+mixture; adipocyte fraction becomes non-zero and correlates with ADIPOQ/PLIN1/LPL expression.
- Permutations disabled in first pass (perm=0): Status = Open; Date opened = 2026-02-28; Default assumption = fractions usable for exploratory covariate sensitivity but no permutation p-values; rerun with perm>0 when runtime is acceptable/stable.
- Heatmaps are group means on VST scale, not per-sample heatmaps: Status = Open; Date opened = 2026-02-28; Default assumption = group-mean heatmaps sufficient for initial receptor-level readout; add per-sample heatmaps if needed.
- Glycolysis “mRNA vs protein contradiction” claim: Status = Resolved; Date opened = 2026-02-28; Resolution = meeting summary explicitly discussed mRNA vs protein paradox for lipolysis genes/proteins (not glycolysis); avoid claiming glycolysis protein discordance without separate evidence.
- Kenichi DOCX gene list includes mouse-style marker synonyms (for example “f4/80”, “CD11c”): Status = Open; Date opened = 2026-03-04; Default assumption = map to human gene symbols where unambiguous (ADGRE1, ITGAX) and treat as coarse markers, not definitive cell identity.

## 6. Verification Status
- Verified: Item #2 scripts executed and produced expected output files; Method = `Rscript --vanilla R/next_steps_2026-02-28_adrenergic_receptor_expression_heatmaps.R`; Result = heatmaps + group-mean table written; Verified = 2026-02-28 08:56.
- Verified: Item #3 scripts executed and produced expected output files; Method = `Rscript --vanilla R/next_steps_2026-02-28_kenichi_panel_by_sex_slides7_8_style.R`; Result = sex-faceted Slide 7-8 style figures + p-value table written; Verified = 2026-02-28 08:56.
- Verified: Item #4 first-pass CIBERSORT pipeline executed and produced expected output files; Method = `CIBERSORT_SIGNATURE_TOP_N=50 CIBERSORT_PERM=0 CIBERSORT_CORES=1 Rscript --vanilla R/next_steps_2026-02-28_cibersort_cell_fractions_adipose_reference_GSE176171.R`; Result = signature + fractions tables + mean-by-sex-age plots for both depots written; Verified = 2026-02-28 08:56.
- Verified: Glycolysis hypothesis PPTX created; Method = `python scripts/make_pptx_2026-02-28_glycolysis_and_update_project.py` + shape-in-bounds check; Result = `docs/menopause_glycolysis_hypothesis_2026-02-28.pptx` written; Verified = 2026-02-28 10:04.
- Verified: Updated project summary PPTX created; Method = `python scripts/make_pptx_2026-02-28_glycolysis_and_update_project.py` + shape-in-bounds check; Result = `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.1.3.pptx` written; Verified = 2026-02-28 10:04.
- Verified: Updated project summary PPTX v4.2 created; Method = `python scripts/make_pptx_2026-02-28_glycolysis_and_update_project.py` + shape-in-bounds check; Result = `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.pptx` written; Verified = 2026-02-28 15:02.
- Verified: Updated project summary PPTX v4.2 created; Method = `python scripts/make_pptx_2026-02-28_glycolysis_and_update_project.py` + shape-in-bounds check; Result = `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.pptx` written; Verified = 2026-02-28 15:02.
- Verified: Updated project summary PPTX v4.2.1 created; Method = `python scripts/make_pptx_2026-02-28_glycolysis_and_update_project.py` + shape-in-bounds check; Result = `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx` written; Verified = 2026-02-28 15:16.
- Verified: Inserted CIBERSORT cell-type explainer slide into v4.2.1; Method = `python-pptx` slide insertion + shape-in-bounds check; Result = `docs/menopause_sns_gtex_enhanced_2026-02-28_v4.2.1.pptx` updated; Verified = 2026-02-28 15:25.
- Not yet verified: Biological plausibility of fractions and stability across settings; Reason = needs comparative checks (marker proxies, expected adipocyte dominance, sensitivity to TopN/perm/reference).
- Verified: CIBERSORT adipocyte-missing diagnosis + mitigation; Method = summarize original fractions (adipocyte ~0) and re-run CIBERSORT excluding MALAT1; Result = adipocyte fraction median ~0.07 (subq) and ~0.054 (visc), and correlates with ADIPOQ/PLIN1/LPL; Verified = 2026-02-28 19:22.
- Verified: Kenichi follow-up gene-set scoring + correlation plots; Method = `Rscript R/next_steps_2026-03-04_kenichi_inflammation_fibrosis_CA_resistance.R`; Result = trends plot + correlation heatmaps + joined score tables under `GTEx_v10_AT_analysis_out/{figs,tables}/next_steps_2026-03-04_kenichi_followup/`; Verified = 2026-03-04 18:48.
- Verified: Deck updated via copy-and-append (no edits to v4.2.1); Method = `python` + `python-pptx`; Result = `docs/menopause_sns_gtex_enhanced_2026-03-04_v4.3.0.pptx` written with 6 new slides appended; Verified = 2026-03-04 18:53.

## 7. Restart Instructions
- Starting point: Read `docs/next_steps_2026-02-28_items_2_3_4_outputs.md` for the output index.
- Starting point: Prefer the latest main deck `docs/menopause_sns_gtex_enhanced_2026-03-04_v4.3.0.pptx` (appends Kenichi follow-up section at end; v4.2.1 remains unchanged).
- Starting point: If extending the Kenichi follow-up, start from `R/next_steps_2026-03-04_kenichi_inflammation_fibrosis_CA_resistance.R` and the gene sets file `references/kenichi_gene_sets_immune_fibrosis_2026-03-04.tsv`.
- Starting point: Review `R/next_steps_2026-02-28_adrenergic_receptor_expression_heatmaps.R`.
- Starting point: Review `R/next_steps_2026-02-28_kenichi_panel_by_sex_slides7_8_style.R`.
- Starting point: Review `R/next_steps_2026-02-28_cibersort_cell_fractions_adipose_reference_GSE176171.R`.
- Recommended next action (1): Validate fraction outputs quickly (inspect distributions; sanity-check adipocyte/ASPC dominance; compare to marker-proxy trends).
- Recommended next action (2): Rerun CIBERSORT with `CIBERSORT_PERM=100` (or higher) if runtime is acceptable; consider `CIBERSORT_SIGNATURE_TOP_N=100`.
- Recommended next action (3): Use the `__fractions_with_meta.tsv.gz` outputs as covariates in existing score/DE models (replace marker-based composition proxies).
- Last updated: 2026-03-04 18:53
