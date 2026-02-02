# Menopause proxy + SNS/lipolysis in human adipose (GTEx v10)
## With technical + composition confounding addressed

Junguk Hur — 2026-02-02

---

## Collaborator goal (Kenichi) + constraints

- Goal: add *human adipose* transcriptomic evidence consistent with adrenergic/lipolysis differences in postmenopause.
- Constraints:
  - GTEx is baseline/unstimulated tissue → acute SNS transcripts may be weak.
  - Chronic SNS can yield catecholamine resistance → directionality not straightforward.
  - Bulk adipose is mixed cell types → composition can dominate signals.

---

## Data + cohorts (what “pre / peri / post” mean here)

- Dataset: **GTEx v10** bulk adipose RNA-seq (counts), analyzed separately by depot.
- Women-only samples:
  - Subcutaneous: n=233
  - Visceral (omentum): n=180
- Menopause proxy by **age-bin midpoint** (GTEx provides age bins, not true menopause status):
  - **pre**: age_mid < 45
  - **peri**: 45 ≤ age_mid ≤ 55
  - **post**: age_mid > 55
- Note (important for interpretation): GTEx age is **10-year bins**; with this midpoint-based rule, groups map to **pre=20–49**, **peri=50–59**, **post=60+**.
- “Strict pre vs post” analyses exclude peri (50–59) to focus on likely extremes.

---

## How gene-set “z-scores” were computed (for Slides 4–6)

For each depot (women-only), using VST-normalized expression:
- Step 1: DESeq2 VST on counts (blind).
- Step 2: map Ensembl → gene symbols; collapse duplicates by mean.
- Step 3: for each gene in a gene set, compute per-sample z-score across samples:
  - z(gene, sample) = (expr − mean_over_samples) / sd_over_samples
- Step 4: gene-set score per sample = mean of z(gene, sample) across genes in the set.

Interpretation: a higher score means **higher relative expression of the genes in that set** *within that depot’s cohort* (not an absolute activity measurement).

---

## Kenichi panel — subcutaneous baseline gene-set activity

**Conclusion (SMCENTER-adjusted):**
- Groups shown: pre=20–49, peri=50–59, post=60+ (age-bin proxy; not true menopause status).
- Score shown: mean gene-wise z-score across the genes in each set (VST expression; z computed across women in this depot).
- Acute SNS immediate-early set: not robust at baseline
- Lipolysis and thermogenesis programs: lower in peri/post vs pre (bulk signal; see composition adjustment later)

![](GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__lipolysis_thermogenesis.png){width=90%}

---

## Kenichi panel — visceral baseline gene-set activity

**Conclusion (SMCENTER-adjusted):**
- Groups shown: pre=20–49, peri=50–59, post=60+ (age-bin proxy; not true menopause status).
- Score shown: mean gene-wise z-score across the genes in each set (VST expression; z computed across women in this depot).
- Acute SNS immediate-early set: not robust at baseline
- Thermogenesis program: lower in peri/post vs pre
- Lipolysis core: not robust in this depot

![](GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__lipolysis_thermogenesis_visceral.png){width=90%}

---

## Technical confounding sensitivity (RIN / ischemic time / Hardy / batch)

**Purpose:** test whether menopause-proxy effects are stable when adding GTEx technical covariates.

![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_covariate_sensitivity/gene_set_covariate_sensitivity.png){width=49%}
![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_covariate_sensitivity/gene_set_covariate_sensitivity_visceral.png){width=49%}

---

## Bulk adipose composition proxies shift with menopause proxy

**Key point:** adipocyte-content proxy decreases in peri/post vs pre in subcutaneous bulk tissue.

![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/marker_scores_by_group.png){width=49%}
![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/marker_scores_by_group_visceral.png){width=49%}

---

## Subcutaneous: composition adjustment explains “lipolysis_core down”

**Key point:** the subcutaneous lipolysis_core “down” signal is strongly composition-sensitive.

![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/kenichi_gene_set_effects_with_composition_adjustment.png){width=90%}

---

## Quick summary (Slides 6–9): what changes after stress-tests?

- Technical covariates (RIN/ischemic/Hardy/batch) can attenuate effects → not all baseline signals are stable.
- Subcutaneous bulk tissue shows a clear **adipocyte-content proxy decrease** in peri/post vs pre.
- Once we adjust for composition proxies, the subcutaneous **lipolysis_core** menopause-proxy difference becomes ~0.
- Therefore: bulk subcutaneous “lipolysis program down” is likely dominated by **composition/adipocyte fraction**, not necessarily per-adipocyte regulation.

---

## Subcutaneous: key table (post vs pre) — base vs composition-adjusted

![](GTEx_v10_AT_analysis_out/figs/ppt_tables_2026-02-01/table_subq_kenichi_sets_composition_adjustment.png){width=95%}

---

## Male vs female context (Kenichi caveat) — trajectories across age bins

- Kenichi note: male–female comparisons are not a perfect control (men also have hormonal aging + psychosocial confounds), but can provide context.
- Here we plot key adipose scores across GTEx age bins in **both sexes** (separate depots).

![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-02_sex_interaction_scores/sex_agebin_score_trajectories_key_scores.png){width=95%}

---

## Male vs female context — sex×age interaction (score-level)

- Model: score ~ SMCENTER + sex + age_group + sex:age_group (base); comparisons: 40–49 vs 60–69 and 30–39 vs 60–69.
- Color = −log10(p) for interaction term (female aging effect minus male aging effect).
- Summary: several aging effects appear stronger in males (consistent with the “male aging may be larger” caveat); female-specific interaction is modest for most scores after composition-aware adjustment.

![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-02_sex_interaction_scores/sex_interaction_scores_heatmap_base.png){width=90%}

---

## Adrenergic/cAMP/desensitization modules (subcutaneous)

**Interpretation:** more compatible with altered responsiveness/chronic adaptation than acute SNS activation at baseline.

![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_adrenergic_modules/adrenergic_module_effect_sensitivity_post_vs_pre.png){width=90%}

---

## Adrenergic/cAMP/desensitization modules (visceral)

![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_adrenergic_modules/adrenergic_module_effect_sensitivity_post_vs_pre_visceral.png){width=90%}

---

## Menopause signature: senMayoLike v2 (curated, evidence-backed)

**Signature size:** 13 total genes (12 up in post, 1 down in post)
- Full per-gene evidence + sources: `docs/senMayoLike_v2_signature_details_2026-02-02.xlsx` (or CSV).

![](GTEx_v10_AT_analysis_out/figs/ppt_tables_2026-02-01/table_senMayoLike_v2_gene_list.png){width=95%}

---

## senMayoLike v2 score in GTEx (subcutaneous + visceral)

**How the score is computed (Slide 13):**
- `z_signed` = mean(z of “up” genes) − mean(z of “down” genes), using the same per-gene z-scoring as earlier.
- `ssgsea_all` = ssGSEA enrichment score (GSVA::ssgsea) using VST expression.
- Conclusion: menopause-proxy groups separate most clearly in subcutaneous (post > pre); visceral is weaker/noisier.

![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_signature_refinement/senMayoLike_v2_scores_by_group.png){width=49%}
![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_signature_refinement/senMayoLike_v2_scores_by_group_visceral.png){width=49%}

---

## External context (not menopause-status): GSE44000 adipocytes (postmeno only)

![](GTEx_v10_AT_analysis_out/figs/external_validation_GSE44000/scores_by_group.png){width=90%}

---

## Recommended manuscript framing (what we can claim)

- Baseline GTEx adipose does not show strong acute SNS immediate-early activation across menopause proxy groups.
- Bulk subcutaneous “lipolysis program down” is largely explained by adipocyte-content/composition proxies.
- Patterns are compatible with chronic remodeling / altered adrenergic-cAMP responsiveness (possible catecholamine resistance).
- senMayoLike v2 menopause signature separates menopause proxy groups (stronger in subcutaneous).

---

## Files to share with collaborators

- Master summary: `docs/master_summary_menopause_sns_gtex_2026-02-01.md`
- senMayoLike v2 (full details + evidence):
  - TSV: `references/menopause_signature_senMayoLike_v2_evidence.tsv`
  - CSV: `docs/senMayoLike_v2_signature_details_2026-02-02.csv`
  - Excel: `docs/senMayoLike_v2_signature_details_2026-02-02.xlsx`
- Key figures:
  - `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__lipolysis_thermogenesis.png`
  - `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/kenichi_gene_set_effects_with_composition_adjustment.png`
  - `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_signature_refinement/senMayoLike_v2_scores_by_group.png`
