# Q&A: Composition, Gene Sets, and “Adipose Scores” (GTEx v10 menopause-proxy deck v4.1.2)

This document answers your questions about the analyses shown in:
- `docs/slide_revision_notes_v4.1.2.md`
- `docs/menopause_sns_gtex_enhanced_2026-02-02_v4.1.2.pptx` (Slides 10–15; plus Slide 9 for the methods preface)

## Quick context (important for interpretation)

- Data: GTEx v10 bulk RNA-seq adipose, analyzed separately by depot:
  - Subcutaneous adipose
  - Visceral adipose (omentum)
- “Menopause status” is a **proxy** based on GTEx age bins (GTEx does not provide true menopause status):
  - `pre`: age-bin midpoint < 45
  - `peri`: 45–55
  - `post`: > 55
- Many figures use **per-sample gene-set scores** computed from VST expression (not classical GSEA output).

## (1) How were the composition analyses performed? (Slides 10–13)

No bulk deconvolution method (e.g., **CIBERSORT**, MuSiC, Bisque, etc.) was used in the Slides 10–13 workflow.

Instead, the deck uses **marker-based composition proxy scores** and then asks a targeted question:
“Do our key gene-set results survive after we include composition proxies as covariates?”

### 1A. What were the “composition proxies”?

Defined as intentionally small, interpretable marker panels in:
- `R/next_steps_2026-02-01_composition_markers.R`

Marker sets:
- `adipocyte`: ADIPOQ, PLIN1, FABP4
- `macrophage_mono`: LST1, C1QC, TYROBP
- `endothelial`: PECAM1, VWF, KDR
- `fibroblast_ecm`: COL1A1, COL3A1, DCN
- `t_cell`: TRAC, CD3D, CD3E

### 1B. How was each marker score computed from bulk RNA-seq?

Per depot, the script:
1. Runs DESeq2 VST on raw counts (`DESeq2::vst(..., blind=TRUE)`).
2. Collapses Ensembl IDs to gene symbols (duplicate symbols averaged).
3. For each marker gene, computes a z-score across samples (within that depot’s cohort).
4. For each sample, the marker score is the **mean** of the gene-wise z-scores across the marker genes in that set.

This is the same “mean gene-wise z across set genes” scoring framework described on Slide 9.

### 1C. Was it “fraction estimation” or “deconvolution”?

It is **not** fraction estimation. It is a **proxy** that tends to move with:
- adipocyte enrichment vs stromal/vascular/immune enrichment,
- immune infiltration,
- fibrosis/ECM remodeling signals.

The main value is interpretability + a low-assumption sensitivity test, not absolute cell fractions.

### 1D. How were composition proxies used statistically?

The key sensitivity test shown on Slides 12–13 is:
- Base model: `score ~ SMCENTER + group`
- Composition-adjusted: `score ~ SMCENTER + group + adipocyte + macrophage_mono + fibroblast_ecm`
- Optional addition of technical covariates: `+ SMRIN + SMTSISCH`

This is implemented directly in:
- `R/next_steps_2026-02-01_composition_markers.R`

## (2) Conclusions of the cell composition analysis (and what it implies for “downregulation”)

### Main empirical findings from Slides 11–13

- Subcutaneous bulk adipose shows a **decrease in the adipocyte proxy score** in peri/post vs pre (Slide 11).
- In subcutaneous adipose, the apparent “lipolysis_core down” signal in bulk is **strongly composition-sensitive**:
  - In the base model it is negative (peri/post < pre).
  - After adding composition proxies, the residual post–pre effect becomes ~0 (Slide 12), leading to the Slide 13 interpretation:
    bulk data alone cannot support a strong “per-adipocyte downregulation” claim.

Concrete example (from `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_composition/kenichi_gene_set_effects_with_composition_adjustment.tsv`):
- Subcutaneous `lipolysis_core`:
  - base post_vs_pre: estimate ≈ −0.46 (p ≈ 4e−4)
  - adj_comp post_vs_pre: estimate ≈ +0.006 (p ≈ 0.84)

### What you can and cannot conclude

What the deck supports (subcutaneous):
- The bulk “lipolysis program down” signal is **consistent with** a shift toward lower adipocyte fraction / higher SVF contribution in post vs pre.
- Therefore, the conservative interpretation is: **composition likely dominates that bulk signal**.

What it does NOT prove:
- It does **not** prove “no per-adipocyte transcriptional change”.
- Bulk + marker proxies cannot uniquely separate:
  - “fewer adipocytes” vs
  - “true per-adipocyte downregulation” vs
  - “both at once”.

Slide 13 explicitly frames the mixture issue as:
`bulk_expr ≈ p_adip·expr_adip + p_SVF·expr_SVF`.

Visceral note:
- Visceral effects are generally weaker/noisier in this deck; thermogenesis_program looks more consistent than lipolysis_core (see also Slide 22 summary and `..._visceral.tsv` outputs).

## (3) How many gene sets are used, and where did the “pathways” on Slide 14 come from?

### 3A. “Kenichi panel” gene sets (hypothesis-driven; not GSEA)

The “Kenichi panel” is **not one** gene set; it is a collection of curated sets saved in:
- `references/kenichi_panel_gene_sets.tsv`

Core sets used in the main Kenichi/SNS figures and sensitivity analyses:
- `acute_beta_adrenergic`
- `adrenergic_receptors`
- `lipolysis_core`
- `thermogenesis_program`

An additional optional context set exists in code (`myers2009_context`) but is not emphasized in the main slides.

These sets are **curated hypothesis modules** (Kenichi email + adipose biology readouts), implemented in:
- `R/kenichi_sns_panel_analysis.R`

### 3B. Composition proxy “gene sets”

Slides 11–13 also use small marker sets (composition proxies). These are also “gene sets” in the scoring sense:
- adipocyte, macrophage/mono, endothelial, fibroblast/ECM, T-cell
Defined in `R/next_steps_2026-02-01_composition_markers.R`.

### 3C. Adrenergic/cAMP/desensitization “modules” (Slide 17)

Slide 17’s “modules” are additional curated gene sets (5 total in v1), defined in:
- `R/next_steps_2026-02-01_adrenergic_modules.R`

Modules:
- `camp_pka_feedback` (PDE/RGS/DUSP/ATF3)
- `receptor_desensitization` (GRKs/ARRBs/RGS2)
- `alpha2_antilipolytic_axis` (ADRA2* + PDE3B)
- `camp_core` (ADCY/GNAS/PKA subunits)
- `immediate_early` (NR4A/FOS/JUN/EGR1 plus feedback genes)

### 3D. Menopause signature gene set (Slides 19–20)

The senMayoLike v2 signature is a **13-gene** curated menopause-proxy signature:
- `references/menopause_signature_senMayoLike_v2.tsv`

### 3E. Slide 14 specifically: where did those “pathways” come from?

Slide 14 (“Subcutaneous: key table (post vs pre) — base vs composition-adjusted”) is a **model summary table** for the curated Kenichi gene sets under:
- base vs composition-adjusted linear models,
not an enrichment result.

It is sourced from the same computations as Slides 11–13:
- `R/next_steps_2026-02-01_composition_markers.R`
- Table data: `GTEx_v10_AT_analysis_out/tables/next_steps_2026-02-01_composition/kenichi_gene_set_effects_with_composition_adjustment.tsv`
- Rendered PPT table image: `GTEx_v10_AT_analysis_out/figs/ppt_tables_2026-02-01/table_subq_kenichi_sets_composition_adjustment.png`

So the “pathways” on Slide 14 are best understood as **curated hypothesis-driven modules**, not GSEA/MSigDB-discovered pathways.

## (4) Slide 15: how was the “adipose score” calculated?

Slide 15 plots “key adipose scores” across GTEx age bins by sex, separately by depot.

Implementation:
- `R/next_steps_2026-02-02_sex_age_trends_plots.R`

What “score” means on Slide 15:
- The same VST-based **gene-set mean z-score** approach is used (`score_set_z`).
- For example, the “Adipocyte proxy” curve is the score of the gene set:
  - ADIPOQ, PLIN1, FABP4

This is **not** GSVA/ssGSEA for Slide 15.

Related note:
- Slide 20 computes the senMayoLike v2 score in two ways:
  - a signed z-score (`z_signed`)
  - and an ssGSEA score (`GSVA::ssgsea`)
But Slide 15’s “trajectory” figure uses the z-score framework (see script).

## (5) Overall summary (in-depth; connecting the pieces)

### What the deck tries to do (and why it’s careful)

Goal: provide human baseline adipose transcriptomic evidence relevant to a “SNS/adrenergic → lipolysis/thermogenesis” hypothesis, while being explicit about three major constraints:
1. GTEx adipose is **baseline/unstimulated** tissue (acute SNS markers may be transient/weak).
2. GTEx adipose is **bulk tissue** (composition shifts can dominate).
3. Technical covariates (center, RIN, ischemic time, etc.) can be non-trivially confounded with age bins.

Therefore, most “findings” are framed as **stress-tested** patterns rather than single-model claims.

### The core result chain (Slides 7–22)

1. Hypothesis-driven Kenichi panel scoring in bulk adipose suggests:
   - Acute SNS immediate-early signature is weak at baseline.
   - Lipolysis/thermogenesis programs show “down” patterns with age/menopause proxy in some comparisons.
2. Technical sensitivity (Slide 10):
   - Adding RIN/ischemic/Hardy/batch covariates can attenuate effects, so results are reported as sensitivity-aware.
3. Composition proxies (Slides 11–13):
   - Subcutaneous shows a clear adipocyte proxy decrease in peri/post vs pre.
   - The subcutaneous “lipolysis_core down” signal essentially vanishes after composition adjustment.
   - Net: bulk subcutaneous “lipolysis down” is likely driven by composition (or at least heavily confounded by it).
4. Adrenergic/cAMP/desensitization modules (Slides 17–18):
   - Even if acute IEG-like SNS activation isn’t seen, there are module shifts compatible with altered adrenergic/cAMP responsiveness and chronic adaptation (a “catecholamine resistance / remodeling” framing).
5. A separate, positive menopause-proxy signal (Slides 19–20):
   - The evidence-backed senMayoLike v2 signature (13 genes) separates menopause-proxy groups, especially in subcutaneous adipose.
6. External adipocyte context (Slide 21):
   - GSE44000 is used only as adipocyte-enriched context (obese vs lean within postmenopausal women), not menopause-status validation.

### Bottom-line interpretation (Slide 22 wording)

Most defensible claims from these analyses:
- Baseline bulk GTEx adipose does not show robust acute SNS immediate-early activation across menopause-proxy groups.
- Subcutaneous lipolysis program differences are strongly composition-sensitive; avoid per-adipocyte downregulation claims from bulk alone.
- Module patterns are compatible with altered adrenergic/cAMP responsiveness and remodeling (possible chronic SNS adaptation/catecholamine resistance framing, with appropriate caution).
- senMayoLike v2 is a clear menopause-proxy-associated signature in GTEx adipose (stronger in subcutaneous).

## Pointers to the exact code and outputs behind Slides 10–15

- Slide 10 (technical ladder): `R/next_steps_2026-02-01_covariate_sensitivity.R`
- Slides 11–14 (composition proxies + adjustment): `R/next_steps_2026-02-01_composition_markers.R`
- Slide 15 (trajectories by sex/age bin): `R/next_steps_2026-02-02_sex_age_trends_plots.R`

