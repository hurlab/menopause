**Findings (Current Risks / Notes)**

- GSE86244 is RNA-seq with processed count matrices in GEO supplementary files; GEOquery `ExpressionSet` downloads can return 0 probes. Use `R/validate_with_gse86244_counts_age_proxy.R` (the old `R/validate_with_gse86244.R` now delegates to it).
- GTEx menopause signature vs GSE86244 age-proxy DE (ASC-focused) currently shows **0 overlapping significant genes (padj < 0.05)**; interpret as “different cohort/cell context” rather than validation.
- MSigDB contains many estrogen/progesterone response sets but no obvious direct “menopause” hallmark; proxy scoring in GTEx shows little shift for hallmark estrogen-response sets across menopause proxy groups.
- The batch-correction DE script `R/run_DE_with_batch_correction.R` may still be logically broken if it intersects subcutaneous and visceral sample IDs (likely empty intersection).
- Prefer the newer VST+GSVA scoring workflow in `R/menopause_signature_scoring_gsva.R` over the older `R/create_menopause_scoring_function.R` for SenMayo-style scoring.

**How The R Codebase Works (Data Flow + Entry Points)**

- Inputs (assumed in-repo): `GTExDatav10/` with GCT count matrices and annotations:
  - `gene_reads_v10_adipose_subcutaneous.gct.gz`
  - `gene_reads_v10_adipose_visceral_omentum.gct.gz`
  - `GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt`
  - `GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt`
- Outputs: everything is written under `GTEx_v10_AT_analysis_out/` into `figs/`, `tables/`, and sometimes `rds/`.

**Core utilities**

- `R/utils.R` provides shared helpers:
  - I/O + prep: `ensure_dirs()`, `read_gct_v12()`, `collapse_dupe_rowsum()`
  - Expression: `collapse_to_symbol_mean()` (maps gene IDs → symbols, collapses duplicates by mean)
  - Annotation + dedup heuristics: `is_valid_human_gene()`, `select_one_gene_id_anyN()`
  - Age parsing: `derive_age_numeric()`, `parse_age_bin()`
  - PCA helpers + plots + PC associations, and a “gene dot plot” helper.

**Big monolithic analysis**

- `R/01_GTEx_analysis_v0.2.R` is the “kitchen sink” script:
  - Reads both depots (subcutaneous + visceral), intersects gene sets, deduplicates gene symbols (keeps a single Ensembl ID per symbol with heuristic rules).
  - Builds metadata (sample attributes + subject phenotypes), derives age midpoints, writes `GTEx_v10_AT_analysis_out/tables/sample_metadata_full.tsv`.
  - Runs VST + PCA, saves PCA plots/scores/loadings and PC~covariate associations.
  - Produces lots of per-gene dotplots.
  - Runs DESeq2 for 40-49 vs 50-59 across sex × depot combinations (notably: design is `~ age_group`, i.e., no batch covariate here), then runs GO/KEGG enrichment via `richR`.

**Batch/bias diagnostics**

- `R/check_batch_confound.R` and `R/check_institution_batch.R` read `GTEx_v10_AT_analysis_out/tables/sample_metadata_full.tsv` (i.e., they depend on `R/01_GTEx_analysis_v0.2.R` having been run) and do chi-square cross-tabs + barplots for batch confounding.

**Focused DE + comparisons**

- `R/DE_wider_age_groups.R`: runs multiple female subcutaneous age-bin contrasts using DESeq2 with batch in the model (`design = ~ SMCENTER + age_group`), saves combined results + summaries/figures.
- `R/compare_batch_correction.R`: runs the same contrast two ways (naive vs with batch covariate) and classifies genes as “lost/gained/stable”; produces volcano/scatter/category plots (MA plot saving is currently wrong).
- `R/run_DE_with_batch_correction.R`: intended to be “correct” DE (batch covariate) + ComBat-corrected PCA for visualization, but currently likely produces an empty cohort due to the sample-ID intersection bug.

**Menopause signature “approaches”**

- `R/menopause_signature_v1_improved_age_proxy.R`: main signature discovery.
  - Defines menopause proxy by age midpoints: Pre (<45), Post (>55), excludes 45–55 “transition”.
  - Runs DESeq2 `design = ~ SMCENTER + menopause_proxy` on female subcutaneous.
  - Writes DE table + a “signature gene list” filtered by `padj < 0.05` and `|log2FC| > 1`.
- `R/menopause_signature_v2_estrogen_filter.R`: runs DE (female 30-39 vs 60-69, batch-corrected model) then filters results to a curated estrogen-responsive gene list and writes an “estrogen-focused signature”.
- `R/menopause_signature_v3_gse86244.R`: intended GEO-based approach (template), but currently not runnable (missing `ensure_dirs` and doesn’t create `GSE86244_data/` before saving).
- `R/menopause_signature_v4_literature.R`: intended “literature signature overlap” script, but currently not runnable (missing `ensure_dirs`) and points at the wrong DE results location.

**Orchestration + downstream**

- `R/run_all_menopause_approaches.R`: spawns the four approach scripts via `system2("Rscript", ...)`, collects status, and writes a markdown summary. Currently brittle because Approach 3/4 aren’t runnable and the signature-file `list.files()` regex is invalid.
- `R/create_menopause_scoring_function.R`: loads the Approach 1 signature list and computes a per-sample score (default “zscore” across signature genes) for female subcutaneous samples; writes score tables + plots. Optional GSVA/ssGSEA branches exist but GSVA branch is currently incorrect.
- `R/perform_enrichment_analysis.R`: runs clusterProfiler GO/KEGG/Reactome enrichment on the Approach 1 signature; writes TSVs + a dotplot. (The optional “load DE results for context” path is currently pointed at the wrong directory.)
- `R/validate_with_gse86244.R`: intended to download GSE86244 and validate overlap vs the GTEx signature using limma + Fisher’s exact test + plots, but probe mapping currently can’t run due to the syntax error noted above.
- `R/batch_correction_guide.R`: just prints best-practice notes; no pipeline effects.

**Natural next steps**

1) If you want, I can turn this into a runnable “happy path” by fixing the hard errors (Approach 3/4, orchestrator regex, GSE86244 validation, batch-corrected DE script).
2) Or I can map a minimal recommended run order (from scratch vs “using existing outputs”) based on what you want to reproduce (PCA/QC vs signature discovery vs validation).

---

## Progress Updates (Jan 26, 2026)

- Fixed the hard runtime errors listed above across the `R/` scripts (regex in `R/run_all_menopause_approaches.R`, missing `ensure_dirs` via `source("R/utils.R")`, broken MA plot saving, broken SUBJID extraction, DESeq2 shrink robustness, enrichment outputs, and the GEO syntax issue).
- Archived the legacy monolithic script: `R/archive/01_GTEx_analysis_v0.2.R` (kept for reference; no longer required for the main pipeline).
- Removed the unused `.claude/` directory from the repo.
- Reran the main analysis scripts; logs are under `GTEx_v10_AT_analysis_out/logs/` and results under `GTEx_v10_AT_analysis_out/tables/` + `GTEx_v10_AT_analysis_out/figs/`.
- Implemented a male-comparator interaction analysis:
  - New script: `R/DE_sex_interaction_age_bins.R`
  - Outputs: `GTEx_v10_AT_analysis_out/tables/sex_interaction_age_bins/` and `GTEx_v10_AT_analysis_out/figs/sex_interaction_age_bins/`
  - Summary table: `GTEx_v10_AT_analysis_out/tables/sex_interaction_age_bins/sex_interaction_summary.tsv`

## Collaborator Inputs (Kenichi, Jan 2026 emails)

- Interpretation: prior results suggested stronger age-associated gene expression changes in males than females; this was considered unexpected and worth framing carefully (male-vs-female comparisons may miss complexity).
- Caveats: men also undergo age-related hormonal changes; older age bins (50s-70s), especially in men, may reflect social/psychological stress effects that could confound adipose transcriptional changes.
- Recommendation: even if less "striking", focus remains valuable on **younger vs older women** alone and demonstrate menopause-related changes, particularly in **SNS-related gene expression**, while acknowledging limitations of baseline (unstimulated) human adipose.
- SNS marker nuance:
  - Acute beta-adrenergic stimulation markers (rapid/transient): **NR4A1/NR4A2/NR4A3**, and immediate-early genes **FOS/JUN/EGR1** (per Myers et al., 2009; PMID: 19465082).
  - Many downstream metabolic genes reported after isoproterenol (e.g., lipid/glucose metabolism genes) were noted as modest/non-specific in that context and not significant after BH FDR in that paper; treat as contextual/supportive rather than primary SNS readouts.
  - Chronic SNS stimulation can reduce **ADRB3** (catecholamine resistance), so higher expression is not necessarily higher SNS tone; interpretation must be cautious.
- Manuscript framing goal (timeline ~1 month): strengthen an estrogen-deficiency insulin-resistance mechanism model (central estrogen -> dopamine/D2 -> SNS outflow -> adipose lipolysis) by demonstrating increased adrenergic signaling / lipolysis-related expression in **postmenopausal human adipose**, consistent with OVX mouse mechanistic data.

## Implications For Next Analysis Steps

- Add targeted, hypothesis-driven readouts alongside genome-wide DE:
  - Curate a small **SNS/adrenergic/lipolysis/thermogenic** gene panel (acute markers vs chronic program separated).
  - Quantify effects in women (pre vs post proxy; and 40-49 vs 50-59) with batch covariates, and report effect sizes + uncertainty (not just padj counts).
  - Use the existing sex-interaction framework to identify signals enriched in women vs men (while explicitly stating that men are not a perfect negative control).

## Targeted Kenichi Panel Results (Subcutaneous, Female-only)

- New script: `R/kenichi_sns_panel_analysis.R`
  - Outputs: `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/` and `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/`
  - Report: `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/kenichi_sns_panel_report.md`
  - Draft message to Kenichi: `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/kenichi_draft_message.md`
  - Log: `GTEx_v10_AT_analysis_out/logs/kenichi_sns_panel_analysis.log`
  - Optional input for collaborator-provided gene lists: `references/kenichi_extra_gene_sets.tsv` (set_name + gene_symbol; comments allowed)
- Menopause proxy definitions used:
  - 3-group model (includes transition): pre (<45), peri (45-55), post (>55) based on age-bin midpoints.
  - Strict 2-group model: pre vs post only (excludes 45-55 / peri).
- Initial readouts (VST z-score gene-set activity; adjusted for SMCENTER):
  - `lipolysis_core`: lower in peri and post vs pre (see `gene_set_activity_stats.tsv`).
  - `acute_beta_adrenergic` (NR4A1/2/3, FOS/JUN/EGR1): small/non-significant shift at the gene-set level; EGR1 shows a detectable increase in strict pre vs post in the panel DE table.
  - `adrenergic_receptors` (includes ADRB3): modest/non-significant shift overall; ADRB3 is lowly expressed and not strongly DE in this baseline GTEx setting.
- Menopause-adjacent bin check (Female 40-49 vs 50-59; adjusted for SMCENTER):
  - Gene-set stats written to `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/gene_set_activity_stats_40-49_vs_50-59.tsv` show no strong shift across the curated sets at this tight boundary in this dataset.
  - DE table for this contrast: `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/DE_bins_40-49_vs_50-59_all_genes.tsv`

- Candidate manuscript figures produced:
  - `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__lipolysis_thermogenesis.png`
  - `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__acute_beta_adrenergic.png`

## Kenichi To-Dos Completed (Jan 26, 2026)

1) Review the panel outputs and decide claim language:
   - Snapshot + key p-values/top hits are embedded in `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/kenichi_sns_panel_report.md`.
   - Draft collaborator message (editable) with suggested conservative wording is in `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/kenichi_draft_message.md`.
   - Current best-supported wording from GTEx baseline adipose: **no robust acute SNS marker upregulation at baseline**; instead **program-level decreases** in lipolysis/thermogenesis across menopause proxy groups, consistent with metabolic remodeling and compatible with chronic SNS/catecholamine-resistance caveats.

2) Pick 1–2 readouts detectable in baseline GTEx adipose and generate figure candidates:
   - Primary: lipolysis + thermogenesis gene-set activity (3-group pre/peri/post, SMCENTER-adjusted):
     - `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__lipolysis_thermogenesis.png`
   - Supporting: acute beta-adrenergic immediate-early gene-set activity:
     - `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__acute_beta_adrenergic.png`

3) Incorporate additional gene lists (Christoph/Kenichi) and regenerate outputs:
   - A drop-in TSV was added: `references/kenichi_extra_gene_sets.tsv` (columns: set_name, gene_symbol; `#` comments allowed).
   - `R/kenichi_sns_panel_analysis.R` auto-loads this file (if present) and regenerates the same tables/figs/report in-place.

## Kenichi-Aligned Deliverable (Manuscript-Support Package)

- Rationale (from Kenichi emails): produce a human adipose result package that:
  - Focuses on women-only comparisons (primary evidence).
  - Treats SNS markers correctly (acute/transient vs chronic/downstream; ADRB3 caveat).
  - Explicitly includes the 45-55 transition range as a sensitivity/initial analysis.
- Implemented:
  - Script: `R/kenichi_sns_panel_analysis.R` (Subcutaneous only; Female only)
  - Models/contrasts:
    - 3-group model including transition: pre (<45), peri (45-55), post (>55) by age-bin midpoint
    - Strict 2-group model (excludes peri): pre (<45) vs post (>55)
    - Menopause-adjacent bins: 40-49 vs 50-59
  - Outputs:
    - Tables: `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/`
    - Figures: `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/`
    - Collaborator-facing report: `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/kenichi_sns_panel_report.md`
    - Log: `GTEx_v10_AT_analysis_out/logs/kenichi_sns_panel_analysis.log`

## Next Steps (Kenichi / ~1 Month Timeline)

1) Review the key outputs and agree on claim language:
   - `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/gene_set_activity_stats.tsv`
   - `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/gene_set_activity_stats_40-49_vs_50-59.tsv`
   - `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/DE_strict2_post_vs_pre_PANEL.tsv`
   - Decide: "increased SNS" vs "altered adrenergic/lipolysis program consistent with chronic stimulation/catecholamine resistance".
2) If a stronger "human supports OVX mechanism" figure is needed:
   - Pick 1-2 readouts that are detectable in baseline GTEx adipose (acute NR4A1/2/3 do not move much here; lipolysis/thermogenesis sets trend lower).
   - Frame explicitly with Kenichi's caveats (baseline tissue; chronic SNS may reduce ADRB3).
3) Incorporate additional materials:
   - Add Christoph's exact gene list(s) and any Kenichi-approved adrenergic/lipolysis gene lists.
   - Rerun `R/kenichi_sns_panel_analysis.R` to regenerate the report/figs in-place.

---

## Progress Updates (Jan 27, 2026)

- Added a versioned/curated menopause signature list: `references/menopause_signature_curated.tsv` (includes `keep` + simple flags).
- Added a SenMayo-style scoring workflow using VST + GSVA methods:
  - Script: `R/menopause_signature_scoring_gsva.R`
  - Outputs: `GTEx_v10_AT_analysis_out/tables/menopause_signature_scoring_gsva/` and `GTEx_v10_AT_analysis_out/figs/menopause_signature_scoring_gsva/`
- Added MSigDB discovery + proxy scoring:
  - Search script: `R/msigdb_search_menopause_related_sets.R`
  - Outputs: `GTEx_v10_AT_analysis_out/tables/msigdb_search_menopause_related/`
  - Proxy set names: `references/msigdb_menopause_proxy_set_names.txt`
  - Scoring script: `R/msigdb_proxy_scoring_gtex.R`
  - Outputs: `GTEx_v10_AT_analysis_out/tables/msigdb_proxy_scoring/` and `GTEx_v10_AT_analysis_out/figs/msigdb_proxy_scoring/`
- Implemented a runnable GSE86244 age-proxy validation from GEO supplementary count matrices:
  - Script: `R/validate_with_gse86244_counts_age_proxy.R`
  - Outputs: `GTEx_v10_AT_analysis_out/tables/gse86244_validation_counts_age_proxy/` and `GTEx_v10_AT_analysis_out/figs/gse86244_validation_counts_age_proxy/`
  - Note: overlap with the current GTEx signature is 0 (padj < 0.05) in this run; see `GTEx_v10_AT_analysis_out/tables/gse86244_validation_counts_age_proxy/GTEx_vs_GSE86244_overlap_summary.tsv`.
- Updated `R/validate_with_gse86244.R` to delegate to the new validation script (avoids the unreliable `ExpressionSet` path).
- Updated Approach 4 to use the new GSE86244 age-proxy gene lists and removed the placeholder “SASP genes” entry:
  - Script: `R/menopause_signature_v4_literature.R`
  - Outputs: `GTEx_v10_AT_analysis_out/tables/approach4_literature/`

---

## Progress Updates (Jan 29, 2026)

- Expanded Kenichi panel analysis to include **both depots** (still female-only; SMCENTER-adjusted):
  - Script: `R/kenichi_sns_panel_analysis.R`
  - Subcutaneous outputs remain unchanged (no suffix).
  - Visceral outputs are written with `_visceral` suffix under the same folders:
    - Tables: `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/*_visceral.tsv`
    - Figures: `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/*_visceral.png`
    - Report: `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/kenichi_sns_panel_report_visceral.md`
    - Draft message: `GTEx_v10_AT_analysis_out/tables/kenichi_sns_panel/kenichi_draft_message_visceral.md`
- Expanded menopause signature scoring (VST + GSVA methods) to include visceral:
  - Script: `R/menopause_signature_scoring_gsva.R`
  - Visceral outputs: `GTEx_v10_AT_analysis_out/tables/menopause_signature_scoring_gsva/*_visceral.tsv` and `GTEx_v10_AT_analysis_out/figs/menopause_signature_scoring_gsva/*_visceral.png`
- Expanded MSigDB estrogen/menopause-proxy scoring to include visceral:
  - Script: `R/msigdb_proxy_scoring_gtex.R`
  - Visceral outputs: `GTEx_v10_AT_analysis_out/tables/msigdb_proxy_scoring/*_visceral.tsv` and `GTEx_v10_AT_analysis_out/figs/msigdb_proxy_scoring/*_visceral.png`
- Saved collaborator-facing deliverables:
  - PPT outline: `docs/ppt_outline_menopause_sns_gtex_v10_15_slides.md`
  - Email draft: `docs/email_to_kenichi_update_gtex_adipose_2026-01-29.md`

---

## Progress Updates (Feb 1, 2026)

- Reviewed key project documents (docx + md), with special attention to collaborator constraints and caveats in `RecentEmailFromKenichi.txt` (baseline/unstimulated tissue; catecholamine resistance; men as an imperfect “negative control”).
- Identified primary remaining gaps to address before making any stronger “human supports ↑SNS/↑lipolysis” claim:
  - **Technical confounding** (ischemic time, RIN, batch) is substantial in GTEx; Kenichi panel models currently adjust for SMCENTER but should be sensitivity-tested with additional covariates (see `A quick summary.docx`).
  - **Cell-type composition** (adipocyte vs SVF) may drive menopause/age signals in bulk adipose; composition checks and/or deconvolution are needed for interpretability.
  - **Interpretation risk:** current GTEx baseline expression is not showing acute SNS markers and shows lipolysis/thermogenesis program activity trending down; this may reflect chronic adaptation/catecholamine resistance and/or composition rather than “no SNS effect.”
- Added a prioritized next-steps analysis + deliverables plan: `docs/next_steps_plan_menopause_sns_gtex_2026-02-01.md`
- Implemented covariate sensitivity ladder (RIN/ischemic time/Hardy/batch) for key Kenichi gene sets and panel genes, both depots:
  - Script: `R/next_steps_2026-02-01_covariate_sensitivity.R`
  - Tables: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_covariate_sensitivity/`
  - Figures: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_covariate_sensitivity/`
  - Quick take: lipolysis/thermogenesis program effects remain directionally negative in subcutaneous under RIN+ischemic adjustment, but are attenuated when adding Hardy+expression batch (likely collinearity / over-adjustment sensitivity); visceral effects are less stable across the most aggressive models.
- Implemented marker-based composition proxy analysis and re-tested Kenichi gene-set effects with composition adjustment (both depots):
  - Script: `R/next_steps_2026-02-01_composition_markers.R`
  - Tables: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_composition/`
  - Figures: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/`
  - Key finding: subcutaneous adipocyte marker score (ADIPOQ/PLIN1/FABP4) is lower in peri/post vs pre; after adjusting for composition proxies, the subcutaneous “lipolysis_core down” gene-set effect becomes ~0 (suggesting bulk composition/adipocyte content is a major driver). Visceral thermogenesis_program downshift is more robust to composition adjustment.
- Added expanded adrenergic/cAMP/desensitization module scoring to better align with the “baseline tissue + catecholamine resistance” caveat (both depots):
  - Script: `R/next_steps_2026-02-01_adrenergic_modules.R`
  - Tables: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_adrenergic_modules/`
  - Figures: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_adrenergic_modules/`
  - Note (subcutaneous): cAMP core module is lower in post vs pre (remains negative after composition adjustment); receptor desensitization module shows a modest negative shift after composition adjustment. These baseline patterns are more compatible with altered responsiveness/chronic adaptation than acute SNS activation.
- Implemented an initial “SenMayo-like” menopause signature refinement workflow and produced two candidate signatures (both scored in subcutaneous + visceral):
  - Script: `R/next_steps_2026-02-01_signature_refinement.R`
  - Outputs: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_signature_refinement/` and `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_signature_refinement/`
  - v1 (data-driven robustness screen): `references/menopause_signature_senMayoLike_v1.tsv`
  - v2 (curated, biology-facing): `references/menopause_signature_senMayoLike_v2.tsv`
  - Evidence skeletons (to be populated with PMIDs/DOIs): `references/menopause_signature_senMayoLike_v1_evidence.tsv`, `references/menopause_signature_senMayoLike_v2_evidence.tsv`
  - Populated gene-by-gene evidence (PMID-based) for v2: `references/menopause_signature_senMayoLike_v2_evidence.tsv` (and removed GALNT5 from v2 due to insufficient direct evidence in this pass).
- Implemented an external adipocyte-enriched context dataset (not menopause status, but all postmenopausal) to help interpret adipocyte programs without bulk composition:
  - Script: `R/external_validation_GSE44000_postmenopausal_adipocytes_obese_vs_lean.R`
  - Outputs: `GTEx_v10_AT_analysis_out/tables/external_validation_GSE44000/` and `GTEx_v10_AT_analysis_out/figs/external_validation_GSE44000/`
  - Note: GSE44000 compares obese vs lean purified adipocytes in postmenopausal women; thermogenesis/lipolysis scores trend lower in obesity (n=14; underpowered), providing context for interpreting adipocyte program directionality in human adipose.
- Wrote a master biology-facing summary and an updated presentation package:
  - Master summary: `docs/master_summary_menopause_sns_gtex_2026-02-01.md`
  - Updated slide outline: `docs/ppt_outline_menopause_sns_gtex_master_2026-02-01.md`
  - Generated pptx (pandoc): `docs/menopause_sns_gtex_master_2026-02-01.pptx`
- Created an enhanced, figure-forward PPTX that embeds the generated plots and table images:
  - Table images: `GTEx_v10_AT_analysis_out/figs/ppt_tables_2026-02-01/table_subq_kenichi_sets_composition_adjustment.png`, `GTEx_v10_AT_analysis_out/figs/ppt_tables_2026-02-01/table_senMayoLike_v2_gene_list.png`
  - Deck markdown: `docs/ppt_deck_enhanced_menopause_sns_gtex_2026-02-01.md`
  - PPTX: `docs/menopause_sns_gtex_enhanced_2026-02-01.pptx`
- Updated the enhanced PPTX to add explicit method explanations (pre/peri/post definitions; z-score calculation; score definitions) and inserted a quick summary slide after the confounding section:
  - Updated deck markdown: `docs/ppt_deck_enhanced_menopause_sns_gtex_2026-02-02.md`
  - Updated PPTX: `docs/menopause_sns_gtex_enhanced_2026-02-02.pptx`
- Exported senMayoLike v2 signature full details (evidence-backed) to shareable formats:
  - CSV: `docs/senMayoLike_v2_signature_details_2026-02-02.csv`
  - Excel: `docs/senMayoLike_v2_signature_details_2026-02-02.xlsx`

---

## Progress Updates (Feb 2, 2026)

- Added male vs female context analyses (Kenichi caveat: “men are not a perfect negative control”):
  - Sex×age-bin interaction (score-level; both depots): `R/next_steps_2026-02-02_sex_interaction_gene_set_scores.R`
    - Table: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-02_sex_interaction_scores/sex_interaction_score_effects.tsv`
    - Heatmap: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-02_sex_interaction_scores/sex_interaction_scores_heatmap_base.png`
  - Sex-by-age-bin trajectories (key scores; both depots): `R/next_steps_2026-02-02_sex_age_trends_plots.R`
    - Plot: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-02_sex_interaction_scores/sex_agebin_score_trajectories_key_scores.png`
- Updated biology-facing summary to include the male vs female context section:
  - `docs/master_summary_menopause_sns_gtex_2026-02-01.md`
- Created a refreshed figure-forward slide deck that:
  - Re-states the age-bin mapping for “pre/peri/post” directly on the Kenichi panel slides,
  - Moves the confounding summary to immediately follow the confounding figures,
  - Adds 2 male vs female context slides (trajectory plot + interaction heatmap),
  - Links the evidence-backed senMayoLike v2 Excel/CSV.
  - Deck markdown: `docs/ppt_deck_enhanced_menopause_sns_gtex_2026-02-02_v2.md`
  - PPTX: `docs/menopause_sns_gtex_enhanced_2026-02-02_v2.pptx`
- Updated the enhanced deck again to add “batch/center” findings near the beginning (after the cohort slide):
  - Rationale: “batch correction” is not plug-and-play because sequencing center/batch is correlated with major covariates (age-bin distribution, depot, sex).
  - Key figures: `GTEx_v10_AT_analysis_out/figs/institution_batch_panel.png`, `GTEx_v10_AT_analysis_out/figs/batch_correction_summary_panel.png`, `GTEx_v10_AT_analysis_out/figs/PCA_combat_correction_comparison.png`
  - Deck markdown: `docs/ppt_deck_enhanced_menopause_sns_gtex_2026-02-02_v3.md`
  - PPTX: `docs/menopause_sns_gtex_enhanced_2026-02-02_v3.pptx`
- Revised collaborator email draft to reflect batch/center + composition sensitivity framing:
  - `docs/email_to_kenichi_update_gtex_adipose_2026-02-02.md`
- Created updated slide decks using your manually revised v4 as the base, preserving its style while adding method/interpretation clarifications:
  - Kenichi panel gene-set definitions file: `references/kenichi_panel_gene_sets.tsv`
  - v4.1.1 (adds requested clarifications + adds a confounding intro slide): `docs/menopause_sns_gtex_enhanced_2026-02-02_v4.1.1.pptx`
    - Notes: `docs/slide_revision_notes_v4.1.1.md`
  - v4.1.2 (same as v4.1.1, but with a slightly more forward-leaning framing slide): `docs/menopause_sns_gtex_enhanced_2026-02-02_v4.1.2.pptx`
    - Notes: `docs/slide_revision_notes_v4.1.2.md`
