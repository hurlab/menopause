# Menopause proxy + SNS/lipolysis in human adipose (GTEx v10)
## With technical + composition confounding addressed

Junguk Hur — 2026-02-01

---

## Collaborator goal (Kenichi) + constraints

- Goal: add *human adipose* transcriptomic evidence consistent with adrenergic/lipolysis differences in postmenopause.
- Constraints:
  - GTEx is baseline/unstimulated tissue → acute SNS transcripts may be weak.
  - Chronic SNS can yield catecholamine resistance → directionality not straightforward.
  - Bulk adipose is mixed cell types → composition can dominate signals.

---

## Data + cohorts

- GTEx v10: Adipose Subcutaneous + Visceral (Omentum), women-only.
- Menopause proxy (age-bin midpoint): pre <45, peri 45–55, post >55.
- Key outputs summarized in: `docs/master_summary_menopause_sns_gtex_2026-02-01.md`

---

## Kenichi panel — subcutaneous baseline gene-set activity

![](GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__lipolysis_thermogenesis.png){width=90%}

---

## Kenichi panel — visceral baseline gene-set activity

![](GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__lipolysis_thermogenesis_visceral.png){width=90%}

---

## Technical confounding sensitivity (RIN / ischemic time / Hardy / batch)

![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_covariate_sensitivity/gene_set_covariate_sensitivity.png){width=49%}
![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_covariate_sensitivity/gene_set_covariate_sensitivity_visceral.png){width=49%}

---

## Bulk adipose composition proxies shift with menopause proxy

![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/marker_scores_by_group.png){width=49%}
![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/marker_scores_by_group_visceral.png){width=49%}

---

## Subcutaneous: composition adjustment explains “lipolysis_core down”

![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/kenichi_gene_set_effects_with_composition_adjustment.png){width=90%}

---

## Subcutaneous: key table (post vs pre) — base vs composition-adjusted

![](GTEx_v10_AT_analysis_out/figs/ppt_tables_2026-02-01/table_subq_kenichi_sets_composition_adjustment.png){width=95%}

---

## Adrenergic/cAMP/desensitization modules (subcutaneous)

![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_adrenergic_modules/adrenergic_module_effect_sensitivity_post_vs_pre.png){width=90%}

---

## Adrenergic/cAMP/desensitization modules (visceral)

![](GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_adrenergic_modules/adrenergic_module_effect_sensitivity_post_vs_pre_visceral.png){width=90%}

---

## Menopause signature: senMayoLike v2 (curated, evidence-backed)

![](GTEx_v10_AT_analysis_out/figs/ppt_tables_2026-02-01/table_senMayoLike_v2_gene_list.png){width=95%}

---

## senMayoLike v2 score in GTEx (subcutaneous + visceral)

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
- Key figures:
  - `GTEx_v10_AT_analysis_out/figs/kenichi_sns_panel/figure_candidate__lipolysis_thermogenesis.png`
  - `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_composition/kenichi_gene_set_effects_with_composition_adjustment.png`
  - `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_signature_refinement/senMayoLike_v2_scores_by_group.png`
- senMayoLike v2 signature + evidence:
  - `references/menopause_signature_senMayoLike_v2.tsv`
  - `references/menopause_signature_senMayoLike_v2_evidence.tsv`

