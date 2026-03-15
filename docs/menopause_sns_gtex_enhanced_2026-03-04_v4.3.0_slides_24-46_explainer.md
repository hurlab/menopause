# Explainer Note for `menopause_sns_gtex_enhanced_2026-03-04_v4.3.0.pptx` (Slides 24-46)

Last updated: 2026-03-14

This note explains the back half of the deck
`docs/menopause_sns_gtex_enhanced_2026-03-04_v4.3.0.pptx`, starting at Slide 24 and continuing through the end.

Scope of this note:
- Slides 24-40 are the appended 2026-02-28 update block added after the original v4.1.2 deck.
- Slides 41-46 are the appended 2026-03-04 Kenichi follow-up block.
- The goal here is not only to restate slide titles, but to explain what each slide is trying to communicate, how it should be presented, and what caveats must be stated out loud.

Important framing:
- This deck is based on baseline GTEx bulk RNA-seq, not stimulated tissue and not direct functional lipolysis assays.
- "Menopause" in the GTEx analyses is a proxy based on age bins, not true observed menopause status.
- Several results in the update section are intentionally shown in sequence to document how interpretation evolved. In particular, the original CIBERSORT fraction plots on Slides 32-33 are not the final corrected interpretation; Slides 34-36 explain the artifact and the corrected noMALAT1 rerun.

## Overall Story Arc for Slides 24-46

Slides 24-29 add detail around receptor expression and sex-stratified versions of the core Kenichi panel. The point of this block is to move from broad pathway scores to more granular transcript-level views and to show where the female-centric menopause-proxy framing sits relative to male aging trajectories.

Slides 30-36 introduce the first-pass CIBERSORT deconvolution, explain the adipocyte-near-zero artifact, and show the corrected rerun after excluding MALAT1. The main message is that composition-aware analysis is necessary, but the first-pass fractions needed QC before they could be trusted biologically.

Slides 37-40 address the glycolysis question directly. The conclusion is conservative: bulk Hallmark Glycolysis is not robustly shifted, while other hallmarks do move, and lipolysis should not be conflated with glycolysis.

Slides 41-46 add the Kenichi follow-up mechanism block. This section shifts from "what changed with menopause proxy" toward "what remodeling and catecholamine-resistance mechanisms might explain the transcript patterns."

## Slide-by-Slide Explainer

## Slide 24
### Title
`Update (2026-02-28): adrenergic receptor expression (depot/sex/age bins)`

### What is on the slide
- Two heatmaps of mean VST expression by group.
- Left: beta-adrenergic receptors (`ADRB1`, `ADRB2`, `ADRB3`).
- Right: alpha2-adrenergic receptors (`ADRA2A`, `ADRA2B`, `ADRA2C`).
- Grouping is by depot x sex x GTEx age bin.

### Why this slide exists
- Earlier parts of the deck discuss adrenergic signaling largely at the module level.
- This slide asks a more direct transcript-level question: do the receptor genes themselves show depot-, sex-, or age-patterning consistent with the broader catecholamine-resistance/remodeling story?

### How to explain it
- Start by saying these are group means on the VST scale, not model coefficients.
- Emphasize that the goal is pattern recognition, not formal inference from the heatmap alone.
- Point out that separating beta and alpha2 receptors is biologically useful because they play different regulatory roles in adipose adrenergic responsiveness.
- Use the slide as a visual bridge to the later Kenichi follow-up where reduced adrenergic receptor expression becomes mechanistically relevant.

### Take-home message
- Receptor expression is not uniform across depot, sex, and age context.
- These genes are plausible transcript-level readouts to carry into more formal composition-aware follow-up models.

### Caveats to say out loud
- Mean VST values are descriptive.
- VST values are safest to compare within a depot, not as absolute cross-depot measurements.
- Receptor expression in bulk tissue can still reflect cell-composition shifts, not only per-cell regulation.

### Backing files
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_receptor_heatmaps/heatmap_ADRB1_ADRB2_ADRB3_meanVST_by_depot_sex_agebin.png`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_receptor_heatmaps/heatmap_ADRA2A_ADRA2B_ADRA2C_meanVST_by_depot_sex_agebin.png`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_receptor_heatmaps/adrenergic_receptor_expression_group_means_vst.tsv`
- Script: `R/next_steps_2026-02-28_adrenergic_receptor_expression_heatmaps.R`

## Slide 25
### Title
`Update (2026-02-28): adrenergic receptor heatmaps with mean VST values (2 decimals)`

### What is on the slide
- The same receptor heatmaps as Slide 24, but with the mean VST printed inside each cell.

### Why this slide exists
- Slide 24 is better for visual pattern recognition.
- Slide 25 is better for precise cross-cell comparison when a viewer wants to know whether an apparent color difference is actually large or modest.

### How to explain it
- Use this slide when a collaborator wants concrete numeric values rather than only relative coloring.
- It is especially useful when discussing whether age-bin changes are gradual or step-like.

### Take-home message
- The receptor heatmaps are not just color impressions; the underlying group means are available directly on the slide.

### Caveats to say out loud
- These are still descriptive group means, not adjusted effect estimates.
- Small numeric shifts can look important if overinterpreted without model context.

### Backing files
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_receptor_heatmaps/heatmap_ADRB1_ADRB2_ADRB3_meanVST_by_depot_sex_agebin_with_values.png`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_receptor_heatmaps/heatmap_ADRA2A_ADRA2B_ADRA2C_meanVST_by_depot_sex_agebin_with_values.png`
- Script: `R/next_steps_2026-02-28_adrenergic_receptor_expression_heatmaps.R`

## Slide 26
### Title
`Update (2026-02-28): adrenergic receptor expression per sample (ordered by age bin)`

### What is on the slide
- Per-sample receptor-expression heatmaps rather than group means.
- Samples are ordered by GTEx age bin within each sex x depot facet.

### Why this slide exists
- Group means can hide heterogeneity.
- This slide shows whether the apparent age/sex/depot patterns are broad shifts across many samples or driven by a small subset.

### How to explain it
- Say that this is the "heterogeneity check" for the receptor story.
- Draw attention to whether signal appears smoothly distributed or patchy.
- Use it to remind the audience that GTEx sample sizes are uneven across sex and age bins, especially in older female bins.

### Take-home message
- There is substantial sample-level variability, so any receptor interpretation needs modeling and not just visual pattern reading.
- Still, the overall structure seen in the group means is not obviously an artifact of a tiny number of outliers.

### Caveats to say out loud
- These per-sample panels are dense and mainly for QC/context.
- Ordering by age bin is for readability, not for suggesting a continuous time course.

### Backing files
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_receptor_heatmaps/heatmap_ADRB1_ADRB2_ADRB3_VST_per_sample_ordered_by_agebin.png`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_receptor_heatmaps/heatmap_ADRA2A_ADRA2B_ADRA2C_VST_per_sample_ordered_by_agebin.png`
- Script: `R/next_steps_2026-02-28_adrenergic_receptor_expression_heatmaps.R`

## Slide 27
### Title
`Update (2026-02-28): Kenichi panel (Slides 7-8 style), faceted by sex`

### What is on the slide
- The original Kenichi-style lipolysis and thermogenesis panels, now plotted separately for female and male within each depot.

### Why this slide exists
- The earlier deck emphasized female menopause-proxy comparisons.
- This slide adds male aging context so that female trends can be interpreted against a non-menopause control axis.

### How to explain it
- Present this as a context slide, not the final causal answer.
- Explain that the same scoring logic used in Slides 7-8 is retained; only the plotting layout changes.
- Use it to ask whether the female trend is unique, stronger, or merely one part of a broader age-associated pattern shared by males.

### Take-home message
- Sex stratification matters.
- Some bulk transcript trends are not purely "female menopause only" and need to be interpreted against male aging and sample-composition differences.

### Caveats to say out loud
- Male and female GTEx adipose sample sizes are not balanced.
- This slide is descriptive and score-based; formal sex interaction testing belongs to modeling slides/tables, not this panel alone.

### Backing files
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_kenichi_panel_by_sex/figure_candidate__lipolysis_thermogenesis__by_sex.png`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_kenichi_panel_by_sex/figure_candidate__lipolysis_thermogenesis__by_sex_visceral.png`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_kenichi_panel_by_sex/kenichi_panel_lm_pvalues_by_sex.tsv`
- Script: `R/next_steps_2026-02-28_kenichi_panel_by_sex_slides7_8_style.R`

## Slide 28
### Title
`Update (2026-02-28): Kenichi panel with female/male interleaved groups (single plot)`

### What is on the slide
- Female and male groups are interleaved in one plot with ordering:
  `F-pre, M-pre, F-peri, M-peri, F-post, M-post`.

### Why this slide exists
- Separate sex facets are useful but make direct side-by-side comparison cumbersome.
- Interleaving makes it easier to see whether female group changes exceed, mirror, or differ from male aging shifts.

### How to explain it
- Use this slide when the discussion is explicitly about sex contrast at roughly matched age/proxy stages.
- Point out that the display is optimized for visual comparison, not for statistical testing.

### Take-home message
- The viewer can compare female and male groups at each stage directly without flipping between panels.

### Caveats to say out loud
- Interleaving improves readability but can overemphasize between-sex comparison if sample-size imbalance is ignored.
- The group labels still represent age-proxy bins, not observed menopause status.

### Backing files
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_kenichi_panel_by_sex/figure_candidate__lipolysis_thermogenesis__F_M_interleaved.png`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_kenichi_panel_by_sex/figure_candidate__lipolysis_thermogenesis__F_M_interleaved_visceral.png`
- Script: `R/next_steps_2026-02-28_kenichi_panel_FM_interleaved_single_plot.R`

## Slide 29
### Title
`Update (2026-02-28): Kenichi panel (interleaved) with medians labeled`

### What is on the slide
- The same interleaved female/male display as Slide 28, but with each group median printed above the distribution.

### Why this slide exists
- This is a readability slide for discussion with collaborators who want exact group summaries without looking up tables.

### How to explain it
- Use it as the most presentation-friendly version of the interleaved panel.
- When discussing a perceived female-vs-male difference, refer to the median labels so the claim is anchored in numbers.

### Take-home message
- The interleaved plots are not just stylistic; the group medians give a quick sense of effect direction and magnitude.

### Caveats to say out loud
- Median labels are summaries, not adjusted model outputs.
- Distribution shape and sample size still matter.

### Backing files
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_kenichi_panel_by_sex/figure_candidate__lipolysis_thermogenesis__F_M_interleaved_with_medians.png`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_kenichi_panel_by_sex/figure_candidate__lipolysis_thermogenesis__F_M_interleaved_with_medians_visceral.png`
- Script: `R/next_steps_2026-02-28_kenichi_panel_FM_interleaved_single_plot.R`

## Slide 30
### Title
`Update (2026-02-28): CIBERSORT reference (GSE176171) used for first-pass fractions`

### What is on the slide
- A reference cell-type count plot for GSE176171.
- Text describing the reference composition, publication, and first-pass CIBERSORT settings.

### Why this slide exists
- Deconvolution results are only interpretable if the audience knows what the reference is.
- This slide makes the reference explicit: broad adipose cell classes from GSE176171, human cells only, pseudo-bulk CPM signature, TopN=50, minCells/type=200, perm=0, QN=false.

### How to explain it
- Say that this is a practical first-pass adipose reference chosen because it ships both a ready-to-use matrix and broad cell labels.
- Emphasize that this is a broad-class deconvolution, not a fine-grained adipocyte substate atlas.
- Mention that the reference supports the main goal here: getting composition covariates good enough to improve interpretation of bulk GTEx adipose.

### Take-home message
- The project moved from crude marker proxies to an explicit reference-based composition estimate.

### Caveats to say out loud
- Reference mismatch is still possible.
- Cell-count percentages in the reference are not expected to equal bulk CIBERSORT output fractions.
- The first-pass run intentionally used `perm=0` to unblock downstream work quickly.

### Backing files
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_cibersort_fractions/gse176171_reference_celltype_counts.png`
- `external_data/GSE176171/derived/reference_cell_type_counts.tsv`
- Script: `R/next_steps_2026-02-28_cibersort_cell_fractions_adipose_reference_GSE176171.R`

## Slide 31
### Title
`Update (2026-02-28): CIBERSORT cell type labels (GSE176171 broad classes)`

### What is on the slide
- Definitions for the broad cell-type labels used in the CIBERSORT reference:
  `adipocyte`, `aspc`, `endothelial`, `lec`, `pericyte`, `smc`, `mesothelium`, `macrophage`, `monocyte`, `dendritic_cell`, `mast_cell`, `t_cell`, `nk_cell`, `b_cell`.

### Why this slide exists
- Broad label names such as `aspc` or `lec` are not universally familiar.
- Without a label glossary, the fraction plots on the next slides are harder to read correctly.

### How to explain it
- Present this as a legend slide.
- For adipose interpretation, highlight `adipocyte`, `aspc`, endothelial classes, and immune classes because those are the most biologically relevant for the composition-remodeling story.

### Take-home message
- The deconvolution is framed at broad adipose cell-class resolution, which is appropriate for composition adjustment but not for strong claims about fine subpopulations.

### Caveats to say out loud
- These are reference labels, not direct annotations of GTEx bulk samples.
- Broad classes hide heterogeneity within immune and stromal compartments.

### Backing context
- Reference metadata from GSE176171.
- Publication cited on Slide 30: Emont et al., Nature 2022.

## Slide 32
### Title
`Update (2026-02-28): CIBERSORT fractions (subcutaneous; cleaned to cell types only)`

### What is on the slide
- The first-pass subcutaneous CIBERSORT fraction plot, cleaned so only cell-type columns remain and the statistics columns (`P-value`, `Correlation`, `RMSE`) are removed from the display.

### Why this slide exists
- It was the first readable presentation of the deconvolution outputs.
- The cleaning step removed distracting non-fraction columns from the plot.

### How to explain it
- Present it historically: this was the initial first-pass fraction view before QC uncovered the adipocyte artifact.
- Use it to show the workflow progression, not as the final trusted result.

### Take-home message
- The project had moved to a deconvolution-based composition view, but this specific first-pass output later required correction.

### Caveats to say out loud
- This slide is superseded for biological interpretation by Slides 34-36.
- Do not present the apparent adipocyte scarcity here as a real biological result.

### Backing files
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_cibersort_fractions/cibersort_subcutaneous_perm0_topN50_minCells200__fractions_mean_by_sex_agebin_clean.png`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/cibersort_subcutaneous_perm0_topN50_minCells200__fractions_with_meta.tsv.gz`
- Script: `R/next_steps_2026-02-28_cibersort_plot_from_fractions_with_meta.R`

## Slide 33
### Title
`Update (2026-02-28): CIBERSORT fractions (visceral; cleaned to cell types only)`

### What is on the slide
- The analogous first-pass cleaned fraction plot for visceral adipose.

### Why this slide exists
- Same reason as Slide 32, but for the visceral depot.

### How to explain it
- Keep it paired with Slide 32 and use both as the "before QC" reference point.

### Take-home message
- The initial deconvolution outputs looked tidy enough to plot, but plausibility checking was still required.

### Caveats to say out loud
- Same as Slide 32: this is not the final corrected fraction interpretation.

### Backing files
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_cibersort_fractions/cibersort_visceral_perm0_topN50_minCells200__fractions_mean_by_sex_agebin_clean.png`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_fractions/cibersort_visceral_perm0_topN50_minCells200__fractions_with_meta.tsv.gz`
- Script: `R/next_steps_2026-02-28_cibersort_plot_from_fractions_with_meta.R`

## Slide 34
### Title
`CIBERSORT QC: why adipocytes were ~0 in the first-pass results`

### What is on the slide
- A text-only QC summary documenting the first-pass problem:
  adipocyte fractions were essentially absent in most samples, which is biologically implausible for adipose tissue.

### Why this slide exists
- This is the critical honesty slide.
- Without it, Slides 32-33 could mislead a reader into thinking adipocytes are truly near-zero in GTEx adipose, which would be nonsense biologically.

### How to explain it
- State the problem plainly:
  the original run produced almost all-zero adipocyte fractions in both depots.
- Then state why this cannot be real:
  canonical adipocyte marker genes such as `ADIPOQ`, `PLIN1`, and `LPL` are clearly abundant in the bulk mixture matrix.
- Finally state the technical diagnosis:
  MALAT1 (`ENSG00000251562`) was dominating the scaling in `IOBR::CIBERSORT`.

### Take-home message
- Composition estimation is necessary, but it must be QC'd against biological plausibility.
- The adipocyte-near-zero result was a technical artifact, not a biological finding.

### Caveats to say out loud
- This is a method-specific issue, not a claim that CIBERSORT is unusable in general.
- The diagnosis is specific to the current signature/mixture setup and the presence of MALAT1.

### Backing files
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_debug/fraction_summary_adipocyte.tsv`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_debug/canonical_presence_in_used_matrices.tsv`
- Script: `R/next_steps_2026-02-28_cibersort_debug_adipocyte_missing.R`

## Slide 35
### Title
`Mitigation + sanity checks (exclude MALAT1 from signature + mixture)`

### What is on the slide
- A text summary of the rerun after excluding MALAT1 from both signature and mixture.
- It reports recovered adipocyte fractions and marker-coherence checks.

### Why this slide exists
- After diagnosing the artifact, the deck needs to show that the mitigation produced biologically more plausible outputs.

### How to explain it
- Say that the rerun used the same overall pipeline but removed MALAT1 from both sides.
- Point out the key improvement:
  adipocyte fractions became non-zero and aligned with known adipocyte markers.
- Mention the coherence checks:
  adipocyte fraction tracks `ADIPOQ`, `PLIN1`, `LPL`;
  ASPC tracks `DCN`, `COL1A2`;
  endothelial tracks `VWF`, `PECAM1`.

### Take-home message
- The corrected rerun restores biological plausibility and makes the fractions usable for exploratory composition adjustment.

### Caveats to say out loud
- This is still an exploratory first-pass deconvolution.
- The handoff still treats sensitivity checks and tuning (`TopN`, `perm`, possibly other exclusions) as unfinished work.

### Backing files
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_debug/rerun_no_malat1__summary_subcutaneous.tsv`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_debug/rerun_no_malat1__summary_visceral.tsv`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_debug/rerun_no_malat1__fractions_subcutaneous.tsv.gz`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_debug/rerun_no_malat1__fractions_visceral.tsv.gz`
- Script: `R/next_steps_2026-02-28_cibersort_debug_adipocyte_missing.R`

## Slide 36
### Title
`Corrected CIBERSORT outputs (noMALAT1) + reference comparison`

### What is on the slide
- Corrected noMALAT1 subcutaneous fraction plot.
- Corrected noMALAT1 visceral fraction plot.
- A reference-vs-bulk comparison plot contrasting GSE176171 cell-count proportions with mean CIBERSORT fractions.

### Why this slide exists
- This is the "after correction" visual that replaces the interpretive role of Slides 32-33.
- It also addresses a predictable audience question:
  should the bulk fractions match the reference cell-count percentages?

### How to explain it
- Say explicitly that these are the corrected plots that should be used for interpretation going forward.
- Explain the reference comparison carefully:
  reference cell-count percentages are not priors or targets for CIBERSORT, so mismatch is not automatically a failure.
- The useful criterion is whether the corrected output is biologically coherent and usable as a covariate, not whether it exactly reproduces atlas cell-count abundances.

### Take-home message
- The noMALAT1 rerun yields a much more plausible adipose composition picture.
- The corrected fractions are the working version for downstream composition-aware analyses.

### Caveats to say out loud
- These fractions are still estimated mixture weights, not literal histology percentages.
- RNA content per cell type and reference mismatch can distort absolute fraction interpretation.

### Backing files
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_cibersort_fractions/cibersort_subcutaneous_perm0_topN50_minCells200_noMALAT1__fractions_mean_by_sex_agebin_clean.png`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_cibersort_fractions/cibersort_visceral_perm0_topN50_minCells200_noMALAT1__fractions_mean_by_sex_agebin_clean.png`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_cibersort_debug/ref_pct_vs_cibersort_mean_noMALAT1.png`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_cibersort_debug/ref_pct_vs_cibersort_mean_noMALAT1.tsv`
- Script: `R/next_steps_2026-02-28_cibersort_reference_vs_bulk_comparison_plot.R`

## Slide 37
### Title
`Update (2026-02-28): glycolysis and related hallmarks (baseline GTEx)`

### What is on the slide
- Two hallmark-summary figures:
  female menopause-proxy group boxplots and sex-by-age-bin trajectory plots.

### Why this slide exists
- The project discussion had drifted toward "glycolysis" as a possible interpretation.
- This slide tests that directly in the same general scoring/modeling framework used elsewhere in the repo.

### How to explain it
- Lead with the negative result:
  Hallmark Glycolysis does not show a robust female pre/peri/post shift in baseline GTEx adipose.
- Then immediately add the more nuanced result:
  OXPHOS, inflammatory response, and adipogenesis do move, especially in subcutaneous female adipose.
- Frame this as evidence for remodeling rather than a simple bulk glycolysis increase or decrease.

### Take-home message
- The data do not support a strong bulk transcriptomic glycolysis shift.
- Other pathway changes are more convincing and align better with tissue remodeling/composition change.

### Caveats to say out loud
- Hallmark scores are broad summaries and may miss gene-node-specific behavior.
- Bulk RNA does not directly measure metabolic flux.

### Backing files
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_glycolysis/female_menopause_groups_hallmarks_boxplots.png`
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-28_glycolysis/sex_agebin_trajectories_hallmarks.png`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_glycolysis/female_pre_peri_post_lm_stats_hallmarks.tsv`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-28_glycolysis/female_pre_peri_post_lm_stats_hallmarks_visceral.tsv`
- Memo: `docs/glycolysis_hypothesis_memo_2026-02-28.md`
- Script: `R/next_steps_2026-02-28_glycolysis_hypoxia_oxphos_analysis.R`

## Slide 38
### Title
`Clarifying Lipolysis vs Glycolysis (GTEx adipose, baseline bulk RNA-seq)`

### What is on the slide
- A text explanation distinguishing lipolysis from glycolysis and clarifying the meaning of the meeting-summary "mRNA vs protein paradox."

### Why this slide exists
- Earlier discussion risked conflating two different biological processes:
  lipolysis and glycolysis.
- The meeting-summary caution was about lipolysis proteins in a high-fat-diet context, not direct GTEx glycolysis evidence.

### How to explain it
- State clearly:
  lipolysis is triglyceride breakdown and is heavily regulated post-transcriptionally through adrenergic signaling and phosphorylation.
- Glycolysis is glucose breakdown and is influenced by different biological drivers, including cell type and inflammatory state.
- This slide resets the vocabulary so downstream claims do not become sloppy.

### Take-home message
- Bulk RNA evidence for lipolysis and glycolysis must be interpreted separately.
- The deck should not claim a glycolysis "mRNA vs protein contradiction" without separate evidence.

### Caveats to say out loud
- GTEx baseline bulk RNA cannot directly resolve pathway flux for either process.
- Composition shifts can distort both pathways in bulk.

### Backing context
- `Menopause_meeting_summary_2026-02-24.txt`
- `docs/glycolysis_hypothesis_memo_2026-02-28.md`

## Slide 39
### Title
`What our results actually show (and why it is easy to mix them up)`

### What is on the slide
- A plain-language interpretation slide summarizing what the project has and has not shown regarding lipolysis, glycolysis, and remodeling.

### Why this slide exists
- By this point in the update block, the audience has seen several kinds of evidence that could easily be overinterpreted.
- This slide prevents overstatement and puts the results back into a conservative narrative.

### How to explain it
- Stress the distinction between:
  "bulk lipolysis-core transcript signal is lower" and
  "actual adipocyte lipolysis flux is lower."
- Then say that Hallmark Glycolysis is largely flat, while inflammatory response, adipogenesis, and OXPHOS show stronger bulk shifts.
- End by saying that if glycolysis remains biologically interesting, the right next step is node-level analysis plus fraction-aware models.

### Take-home message
- The strongest bulk story is adipose remodeling and composition sensitivity, not a clean glycolysis shift.

### Caveats to say out loud
- This slide is a synthesis, not a new statistical result.
- It relies on the correction and interpretation of earlier slides.

### Backing context
- `docs/glycolysis_hypothesis_memo_2026-02-28.md`
- `docs/menopause_sns_gtex_standalone_summary_2026-03-04.md`

## Slide 40
### Title
`Update (2026-02-28): new outputs added since v4.1.2`

### What is on the slide
- A file map of the new receptor, Kenichi, CIBERSORT, and glycolysis outputs added in the February 28 update cycle.

### Why this slide exists
- This is the reproducibility and collaborator-handoff slide.
- It tells other project participants where to find the exact artifacts behind the update section.

### How to explain it
- Use it as a logistical slide:
  "if you want any figure or supporting table from this update block, these are the directories."
- Mention that the CIBERSORT QC/correction work was added after the initial fraction plots.

### Take-home message
- The appended update block is backed by concrete figures, tables, scripts, and a separate glycolysis memo.

### Caveats to say out loud
- The file list includes both first-pass and corrected CIBERSORT materials; for interpretation, the corrected noMALAT1 outputs are the relevant ones.

### Backing context
- `docs/next_steps_2026-02-28_items_2_3_4_outputs.md`
- `docs/glycolysis_hypothesis_memo_2026-02-28.md`

## Slide 41
### Title
`Update (2026-03-04): Kenichi follow-up`

### What is on the slide
- A section header introducing the March 4 follow-up block:
  inflammation, fibrosis, TGF-beta-related signaling, and obesity-linked catecholamine-resistance mechanisms.

### Why this slide exists
- The earlier deck showed descriptive and composition-aware context.
- This new section asks a deeper mechanistic question:
  could an inflammatory/fibrotic remodeled adipose state help explain reduced adrenergic responsiveness?

### How to explain it
- Use this as a transition slide.
- Say that the focus is shifting from broad pathway patterns to candidate mechanism axes motivated by Kenichi's follow-up note and the adipose inflammation review literature.

### Take-home message
- The deck is moving from "what changed?" to "what mechanism could connect remodeling with adrenergic/lipolysis phenotypes?"

### Backing context
- `references/kenichi_gene_sets_immune_fibrosis_2026-03-04.tsv`
- `references/nrendo.2017.90.pdf`

## Slide 42
### Title
`Hypothesis and analysis strategy (Kenichi follow-up)`

### What is on the slide
- A text summary of the working hypothesis and the analysis design.
- It ties inflammation/fibrosis/TGF-beta/catecholamine-resistance literature to new GTEx tests.

### Why this slide exists
- The audience needs the logic before seeing the trend and correlation panels.
- It explains why `ACVR1C`, `PDE3B`, `TBK1`, `IKBKE`, and related axes were brought into the project.

### How to explain it
- Start with the biological question:
  obesity-associated catecholamine resistance is known; whether estrogen deficiency induces an analogous state is unclear.
- Then describe the literature-linked candidate axes:
  `TGF-beta -> ALK7 (ACVR1C) -> reduced ADRB3`
  and
  `TNF -> NF-kB -> IKKepsilon/TBK1 -> PDE3B -> dampened cAMP signaling`.
- Then say what was actually done in GTEx:
  score curated gene sets, extract key mechanism genes, and correlate them with receptor expression and corrected CIBERSORT fractions.

### Take-home message
- The follow-up block is a mechanistic screening analysis, not a definitive causal test.

### Caveats to say out loud
- The preliminary associations on this slide are not yet composition-adjusted formal models.
- The literature mechanism comes from obesity/inflammation context, not directly from menopause studies.

### Backing files
- `references/kenichi_gene_sets_immune_fibrosis_2026-03-04.tsv`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-03-04_kenichi_followup/kenichi_and_review_gene_sets.tsv`
- Script: `R/next_steps_2026-03-04_kenichi_inflammation_fibrosis_CA_resistance.R`

## Slide 43
### Title
`GTEx adipose: Kenichi gene-set trends (age bins x sex x depot)`

### What is on the slide
- A faceted line-and-error-bar panel showing mean gene-set scores by age bin, sex, and depot.
- Included sets:
  inflammation cytokines,
  macrophage markers,
  fibrosis/ECM,
  senescence,
  TGFb/ALK7 axis,
  TNF/NFkB/IKKepsilon/TBK1/PDE3B axis.

### Why this slide exists
- It shows whether these candidate remodeling/mechanism axes have visible age- and sex-structure before formal regression modeling.

### How to explain it
- Explain the encoding first:
  mean z-score by group, error bars are about 95 percent CI, faceted by set and depot, with male and female trajectories overlaid.
- Then give the intended interpretation:
  this is a screening view to identify which axes appear most age- or sex-sensitive and whether subcutaneous looks more dynamic than visceral.
- Connect it back to the rest of the deck:
  if senescence, fibrosis, or inflammatory axes rise where adrenergic/lipolysis programs fall, that supports a remodeling-based story.

### Take-home message
- The follow-up block is not limited to one pathway; it compares several plausible remodeling and catecholamine-resistance axes across both depots and sexes.

### Caveats to say out loud
- Female counts are smaller than male counts in GTEx adipose, especially in older bins.
- These are still descriptive means, not adjusted causal estimates.

### Backing files
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-03-04_kenichi_followup/kenichi_followup_gene_set_trends_by_age_sex__both_depots.png`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-03-04_kenichi_followup/sample_counts_by_sex_agebin__subcutaneous.tsv`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-03-04_kenichi_followup/sample_counts_by_sex_agebin__visceral.tsv`
- Script: `R/next_steps_2026-03-04_kenichi_inflammation_fibrosis_CA_resistance.R`

## Slide 44
### Title
`GTEx subcutaneous: correlations (gene sets, mechanism genes, CIBERSORT fractions)`

### What is on the slide
- A Pearson correlation heatmap for the subcutaneous depot.
- Variables include gene-set scores, key mechanism genes, adrenergic receptors, and selected noMALAT1 CIBERSORT fractions.

### Why this slide exists
- This is the key integrative slide for subcutaneous adipose, the depot that has generally shown the strongest menopause-proxy/remodeling signals in the project.

### How to explain it
- Tell the audience how to read the panel:
  red means positive correlation, blue means negative correlation, and the rows/columns include both pathways and composition variables.
- Then call out the main qualitative patterns:
  inflammation, fibrosis, and related remodeling programs tend to oppose adipocyte-rich signatures and can track inversely with adrenergic receptor expression.
- A useful concrete example from the saved analyses is that in subcutaneous tissue, fibrosis and inflammation-related scores are negatively correlated with `ADRB2`, and inflammatory/fibrotic scores are also negatively associated with adipocyte fraction.

### Take-home message
- In subcutaneous adipose, the candidate mechanism axes are embedded in the same correlation structure as composition shifts and receptor-expression changes.

### Caveats to say out loud
- Correlation is not adjustment.
- Some relationships may be driven mainly by composition rather than per-cell transcriptional regulation.
- This slide motivates the adjusted models; it does not replace them.

### Backing files
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-03-04_kenichi_followup/kenichi_followup_correlations__subcutaneous.png`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-03-04_kenichi_followup/kenichi_followup_correlations__subcutaneous.tsv`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-03-04_kenichi_followup/kenichi_followup_female_fraction_correlations.tsv`
- Script: `R/next_steps_2026-03-04_kenichi_inflammation_fibrosis_CA_resistance.R`

## Slide 45
### Title
`GTEx visceral: correlations (gene sets, mechanism genes, CIBERSORT fractions)`

### What is on the slide
- The same correlation analysis as Slide 44, but for the visceral depot.

### Why this slide exists
- The project repeatedly sees depot differences.
- This slide tests whether the subcutaneous correlation structure generalizes to visceral adipose or weakens there.

### How to explain it
- Present it as the depot comparison slide.
- Ask whether the same remodeling-vs-adrenergic opposition is present, weaker, or differently organized in visceral tissue.
- This helps avoid overgeneralizing a subcutaneous-driven story to all adipose depots.

### Take-home message
- Depot matters. Any catecholamine-resistance/remodeling model coming out of this project should be checked separately in subcutaneous and visceral adipose.

### Caveats to say out loud
- Sample sizes and correlation structures differ by depot.
- Lack of a strong pattern in visceral does not invalidate subcutaneous findings; it may reflect genuine depot biology.

### Backing files
- `GTEx_v10_AT_analysis_out/figs/next_steps_2026-03-04_kenichi_followup/kenichi_followup_correlations__visceral.png`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-03-04_kenichi_followup/kenichi_followup_correlations__visceral.tsv`
- Script: `R/next_steps_2026-03-04_kenichi_inflammation_fibrosis_CA_resistance.R`

## Slide 46
### Title
`Concrete next analyses / validations`

### What is on the slide
- A forward-looking action slide listing the next computational and experimental tests.

### Why this slide exists
- The prior five slides raise a mechanistic hypothesis but do not finish the modeling.
- This slide converts the follow-up block into an explicit next-step plan.

### How to explain it
- Present the top computational priority first:
  within each depot, fit models such as
  `(gene set score or key gene VST) ~ menopause_bin * sex + SMCENTER + adipocyte_fraction + optional macrophage_fraction`.
- Then explain the interpretive goal:
  determine whether inflammation/fibrosis/TGF-beta-linked axes still track receptor or lipolysis readouts after composition adjustment.
- End with the experimental bridge:
  qPCR and ex vivo lipolysis / beta-agonist response assays are the most direct ways to test whether the transcript patterns reflect true catecholamine resistance.

### Take-home message
- The March 4 section is an initial mechanistic screen.
- The real decision point is whether the suggested axes remain associated after explicit composition adjustment and whether that can be validated functionally.

### Caveats to say out loud
- This slide intentionally lists hypotheses and validations, not results.
- The models described here were still marked as not yet started in the handoff at the time this deck was created.

### Backing files
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-03-04_kenichi_followup/kenichi_followup_female_model_effects.tsv`
- `GTEx_v10_AT_analysis_out/tables/next_steps_2026-03-04_kenichi_followup/kenichi_followup_female_fraction_correlations.tsv`
- Script: `R/next_steps_2026-03-04_kenichi_inflammation_fibrosis_CA_resistance.R`

## Practical Presentation Notes

If presenting Slides 24-46 to a collaborator or PI, the safest sequence is:

1. Slides 24-29:
   receptor detail plus sex-stratified context.
2. Slides 30-36:
   composition estimation, artifact, correction, and why the corrected fractions matter.
3. Slides 37-39:
   glycolysis clarification so the deck does not overclaim.
4. Slide 40:
   reproducibility/file map.
5. Slides 41-46:
   mechanism-oriented follow-up and next-step plan.

Three points that should always be stated explicitly:

1. GTEx here is baseline bulk RNA-seq, so transcript changes do not equal lipolysis flux or acute SNS activation.
2. Composition confounding is a central part of the story, not a side note.
3. The noMALAT1 CIBERSORT rerun is the biologically interpretable deconvolution result; the earlier first-pass plots are historical context only.

## Short Summary in One Paragraph

Slides 24-46 document the project's transition from descriptive bulk transcript patterns toward a more composition-aware and mechanism-aware interpretation. The receptor and sex-stratified Kenichi slides show why simple female-only bulk summaries are incomplete. The CIBERSORT block demonstrates that deconvolution is necessary but technically fragile, and that the MALAT1-corrected noMALAT1 rerun is the usable version. The glycolysis block narrows the claim to remodeling and pathway-specific caution rather than a strong glycolysis conclusion. The final March 4 block proposes an inflammation/fibrosis/TGF-beta/catecholamine-resistance framework that is suggestive in baseline GTEx and sets up the next composition-adjusted and experimental validation steps.
