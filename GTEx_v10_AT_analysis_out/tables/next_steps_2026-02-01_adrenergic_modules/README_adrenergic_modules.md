# Adrenergic module scoring summary (Adipose - Subcutaneous; female-only)

This folder scores small adrenergic/cAMP/desensitization modules intended to better match baseline tissue and chronic SNS adaptation caveats.

## Modules (v1)
- camp_pka_feedback: PDE4D, PDE4B, PDE3B, RGS2, DUSP1, DUSP2, ATF3
- receptor_desensitization: ADRBK1, ADRBK2, ARRB1, ARRB2, RGS2
- alpha2_antilipolytic_axis: ADRA2A, ADRA2B, ADRA2C, PDE3B
- camp_core: ADCY3, ADCY6, GNAS, PRKACA, PRKACB, PRKAR1A, PRKAR2B
- immediate_early: NR4A1, NR4A2, NR4A3, FOS, JUN, EGR1, DUSP1, ATF3

## Outputs
- Per-sample module scores: `adrenergic_module_scores_per_sample.tsv`
- Module effect stats (base vs adjusted models): `adrenergic_module_effect_stats.tsv`
- Module boxplots: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_adrenergic_modules/adrenergic_module_scores_by_group.png`
- Effect sensitivity plot: `GTEx_v10_AT_analysis_out/figs/next_steps_2026-02-01_adrenergic_modules/adrenergic_module_effect_sensitivity_post_vs_pre.png`

Generated:
2026-02-01 00:56:16.221177
